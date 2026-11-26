/*!
 * @file OpenFIRE_Scope_Video.cpp
 * @brief Video of the scope: newest frame, roll rotation and display
 * @n CPP Video of the scope: newest frame, roll rotation and display
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#include <Arduino.h>
#include "OpenFIRE_Scope_Video.h"

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32) && defined(SCOPE) && defined(USES_DISPLAY) && defined(USE_LOVYAN_GFX)

#include <string.h>
#include "freertos/semphr.h"
#include "esp_heap_caps.h"
#include "OpenFIRE_Wireless_Scope.h"
#include "OpenFIRE_Scope_Decoder.h"

// EN: From the main: state of the display (Command of the S3) and roll of the lightgun.
// IT: Dal main: stato del display (Comando dell'S3) e roll della lightgun.
extern volatile bool   scope_display_on;
extern volatile int8_t scope_roll;
extern volatile bool   scope_logo_spin_on;   // sprite of the turning logo in memory / sprite del logo che gira in memoria

#define SLICE_BYTES   (SCOPE_DEC_W * SCOPE_DEC_SLICE_H * 2)   // 7680
#define ALL_SLICES    ((1u << SCOPE_DEC_SLICES) - 1)

static LGFX             *g_tft = NULL;
static SemaphoreHandle_t g_panelMutex = NULL;
static bool              g_panelOn = true;        // power applied to the panel, under the lock / alimentazione data al pannello, sotto il blocco
static volatile bool     g_noSignal = false;

// EN: Buffers: decoded frame (15 blocks of 240x16, one per slice, and the table of its rows), work slice of the
//     decoder, output slices of the rotation (DMA).
// IT: Buffer: frame decodificato (15 blocchi da 240x16, uno per fascia, e la tabella delle sue righe), fascia di lavoro
//     del decoder, fasce di uscita della rotazione (DMA).
static bool      g_buffers = false;
static uint16_t *g_slice[SCOPE_DEC_SLICES];
static const uint16_t *g_rows[SCOPE_DEC_H];
static uint16_t *g_work = NULL;
static uint16_t *g_out[SCOPE_VIDEO_OUT_SLICES];
static uint8_t   g_nOut = 0;

#ifdef SCOPE_VIDEO_DEBUG
static uint32_t d_frames = 0, d_slices = 0, d_corrupt = 0, d_old = 0, d_dup = 0, d_decUs = 0, d_rotUs = 0;
#endif // SCOPE_VIDEO_DEBUG

// ===================================================================================
// PANEL LOCK / BLOCCO DEL PANNELLO
// ===================================================================================
static bool g_owned = false;   // the video has the panel / il video ha il pannello

bool scopePanelTryLock() { return g_panelMutex != NULL && xSemaphoreTake(g_panelMutex, 0) == pdTRUE; }
void scopePanelUnlock()  { xSemaphoreGive(g_panelMutex); }

void scopePanelPowerApply() {
  const bool on = scope_display_on;
  if (on == g_panelOn) return;
  g_panelOn = on;
  if (on) g_tft->wakeup(); else g_tft->sleep();
}

bool scopeVideoNoSignal() { return g_noSignal; }

// ===================================================================================
// OUTPUT TO THE DISPLAY / USCITA VERSO IL DISPLAY
// ===================================================================================
// EN: The rotated slices wait in a queue; outPump() starts the next one as soon as the DMA is free, without waiting.
//     It is called after every MCU of the decoder and between the packets, so the DMA keeps working while the
//     next frame is decoded. A buffer is free again when its DMA has finished.
// IT: Le fasce ruotate aspettano in una coda; outPump() fa partire la successiva appena il DMA e' libero, senza aspettare.
//     Si chiama dopo ogni MCU del decoder e tra un pacchetto e l'altro, cosi' il DMA continua a lavorare mentre si
//     decodifica il frame successivo. Un buffer torna libero quando il suo DMA ha finito.
static bool    g_bufBusy[SCOPE_VIDEO_OUT_SLICES];
static uint8_t g_qSy[SCOPE_VIDEO_OUT_SLICES], g_qBuf[SCOPE_VIDEO_OUT_SLICES];
static uint8_t g_qHead = 0, g_qCount = 0;
static int     g_inflight = -1;   // buffer on the DMA / buffer sul DMA

static void outPump() {
  if (g_inflight >= 0 && !g_tft->dmaBusy()) { g_bufBusy[g_inflight] = false; g_inflight = -1; }
  if (g_qCount == 0 || g_inflight >= 0) return;
  const uint8_t sy = g_qSy[g_qHead], b = g_qBuf[g_qHead];
  g_qHead = (uint8_t)((g_qHead + 1) % g_nOut);
  g_qCount--;
  // EN: Data already in the byte order of the display (swap565): straight to the DMA, no conversion.
  // IT: Dati gia' nell'ordine di byte del display (swap565): dritti al DMA, nessuna conversione.
  g_tft->pushImageDMA(scopeRotWindowX(sy), sy * SCOPE_DEC_SLICE_H, scopeRotWindowW(sy), SCOPE_DEC_SLICE_H,
                      (const lgfx::swap565_t *)g_out[b]);
  g_inflight = b;
}

static bool outIdle() { return g_qCount == 0 && g_inflight < 0; }

static void outFlush() { while (!outIdle()) outPump(); }

static uint8_t outFreeBuffer() {
  for (;;) {
    for (uint8_t b = 0; b < g_nOut; b++) if (!g_bufBusy[b]) return b;
    outPump();   // all busy: waits for the DMA / tutti occupati: aspetta il DMA
  }
}

// EN: Roll of the image from the roll received, once per frame (see SCOPE_ROLL_ZERO_ZONE in the .h).
//     A change of FAST degrees or more starts a fast turn: for FOLLOW_MS after each change applied at once, every
//     further degree in the same direction is applied at once too, so the image keeps up with the turn.
// IT: Roll dell'immagine dal roll ricevuto, una volta per frame (vedi SCOPE_ROLL_ZERO_ZONE nel .h).
//     Un cambio di FAST gradi o piu' avvia una rotazione veloce: per FOLLOW_MS dopo ogni cambio applicato subito, ogni
//     grado in piu' nella stessa direzione e' applicato subito anche lui, cosi' l'immagine tiene il passo.
static int rollShown(int roll) {
  static const int FAST = 5;              // degrees applied at once / gradi applicati subito
  static const uint32_t LEVEL_MS = 250;   // within +-1 degree for this long -> straight / entro +-1 grado per questo tempo -> dritta
  static const uint32_t FOLLOW_MS = 150;  // fast turn followed degree by degree / rotazione veloce seguita grado per grado
  static const uint32_t GAP_MS = 200;     // longer pause of the video: start again / pausa del video piu' lunga: si riparte da capo
  static int shown = 0;                   // roll of the image / roll dell'immagine
  static int dir = 0;                     // direction of the fast turn / direzione della rotazione veloce
  static int waitJump = 0;                // change waited / cambio atteso
  static bool level = false, wait = false, follow = false, started = false;
  static uint32_t levelMs = 0, waitMs = 0, followMs = 0, lastMs = 0;
  const uint32_t now = millis();
  if (started && now - lastMs > GAP_MS) { wait = false; level = false; follow = false; }   // old readings / letture vecchie
  started = true;
  lastMs = now;
  if (roll >= -1 && roll <= 1) { if (!level) { level = true; levelMs = now; } }
  else level = false;
  if (follow && now - followMs > FOLLOW_MS) follow = false;
  int want = shown;
  bool back = false;
  if (SCOPE_ROLL_ZERO_ZONE > 0 && shown == 0) {
    if (roll > SCOPE_ROLL_ZERO_ZONE || roll < -SCOPE_ROLL_ZERO_ZONE) want = roll;   // leaves the zone / esce dalla zona
  } else if (SCOPE_ROLL_ZERO_ZONE > 0 && level && now - levelMs >= LEVEL_MS) {
    want = 0;                                                                       // back straight / di nuovo dritta
    back = true;                                                                    // already stable for 250 ms / gia' stabile da 250 ms
  } else {
    const int diff = roll - shown;
    if (diff >= SCOPE_ROLL_HYSTERESIS || diff <= -SCOPE_ROLL_HYSTERESIS || (follow && diff * dir >= 1)) want = roll;
  }
  if (want == shown) { wait = false; return shown; }
  const int jump = want - shown;
  const bool big = jump >= FAST || jump <= -FAST;
  if (SCOPE_ROLL_STABLE_MS == 0 || big || back || (follow && jump * dir > 0)) {   // at once / subito
    if (big || follow) { follow = true; followMs = now; }
    dir = jump > 0 ? 1 : -1;
    shown = want;
    wait = false;
    return shown;
  }
  if (wait && (jump > 0) != (waitJump > 0)) wait = false;   // the direction changed: start again / la direzione e' cambiata: si riparte
  if (!wait) { wait = true; waitMs = now; waitJump = jump; }   // the change must last / il cambio deve durare
  if (now - waitMs >= SCOPE_ROLL_STABLE_MS) { shown = want; wait = false; }
  return shown;
}

// EN: Rotation of the whole frame with the last angle received: when it returns the frame can be overwritten,
//     the last slices are still going to the display.
// IT: Rotazione dell'intero frame con l'ultimo angolo ricevuto: quando torna il frame si puo' sovrascrivere,
//     le ultime fasce stanno ancora andando al display.
static void rotateFrame() {
  #ifdef SCOPE_VIDEO_DEBUG
  const uint32_t t0 = micros();
  #endif // SCOPE_VIDEO_DEBUG
  int roll = scope_roll;
  if (roll > SCOPE_ROLL_MAX) roll = SCOPE_ROLL_MAX;
  if (roll < -SCOPE_ROLL_MAX) roll = -SCOPE_ROLL_MAX;
  roll = rollShown(roll);
  #ifdef SCOPE_ROLL_INVERT
    const int angle = roll;
  #else
    const int angle = -roll;   // the image turns opposite to the lightgun / l'immagine gira al contrario della lightgun
  #endif // SCOPE_ROLL_INVERT
  int32_t c, s;
  scopeRotAngle(angle, &c, &s);
  for (uint8_t sy = 0; sy < SCOPE_DEC_SLICES; sy++) {
    const uint8_t b = outFreeBuffer();
    g_bufBusy[b] = true;
    scopeRotSlice(g_rows, sy, c, s, g_out[b]);
    const uint8_t tail = (uint8_t)((g_qHead + g_qCount) % g_nOut);
    g_qSy[tail] = sy;
    g_qBuf[tail] = b;
    g_qCount++;
    outPump();
  }
  #ifdef SCOPE_VIDEO_DEBUG
  d_rotUs += micros() - t0;
  #endif // SCOPE_VIDEO_DEBUG
}

static void panelTake() {
  if (g_owned) return;
  xSemaphoreTake(g_panelMutex, portMAX_DELAY);
  g_tft->startWrite();   // one transaction for the whole video: the DMA is not waited at every slice / una transazione per tutto il video: il DMA non si aspetta a ogni fascia
  g_owned = true;
  scopePanelPowerApply();
}

static void panelGive() {
  if (!g_owned) return;
  outFlush();
  g_tft->endWrite();
  g_owned = false;
  xSemaphoreGive(g_panelMutex);
}

// ===================================================================================
// NEWEST FRAME / FRAME PIU' RECENTE
// ===================================================================================
static bool     g_active = false;     // video on screen (phase 3) / video sullo schermo (fase 3)
static uint8_t  g_session = 0;        // 0 = none since the pairing / 0 = nessuna dall'abbinamento
static bool     g_haveLast = false;   // a frame of this session already rotated / un frame di questa sessione gia' ruotato
static uint16_t g_lastFrame = 0;
static bool     g_building = false;   // frame being received / frame in ricezione
static uint16_t g_buildFrame = 0;
static uint16_t g_got = 0;            // slices received of the frame being received / fasce ricevute del frame in ricezione
static uint32_t g_firstMs = 0;        // first slice of the frame being received / prima fascia del frame in ricezione
static uint32_t g_sliceOkMs[SCOPE_DEC_SLICES];   // last valid decoding of every slice / ultima decodifica valida di ogni fascia
static uint16_t g_black = 0;          // slices already black in the frame / fasce gia' nere nel frame
static bool     g_refresh = false;    // display back on: the last frame again / display riacceso: di nuovo l'ultimo frame

static inline uint16_t pktFrame(const uint8_t *p) { return (uint16_t)(p[SCOPE_VIDEO_B_FRAME] | (p[SCOPE_VIDEO_B_FRAME + 1] << 8)); }
static inline bool frameNewer(uint16_t a, uint16_t b) { return (int16_t)(a - b) > 0; }   // 16-bit with sign / a 16 bit con segno
// EN: Sessions 1-255: d = (new - current) mod 255; 0 same, 1-127 newer, 128-254 older.
// IT: Sessioni 1-255: d = (nuova - attuale) mod 255; 0 stessa, 1-127 piu' recente, 128-254 vecchia.
static inline bool sessionNewer(uint8_t n, uint8_t cur) { const int d = ((int)n - (int)cur + 255) % 255; return d >= 1 && d <= 127; }

static void stopVideo() {
  if (!g_active) return;
  g_active = false;
  WirelessScope.videoActive(false);
  panelGive();
}

static void showFrame() {
  g_building = false;
  g_lastFrame = g_buildFrame;
  g_haveLast = true;
  // EN: A slice missing from every frame of the last 100 ms goes black instead of an image piece standing still.
  // IT: Una fascia mancante da tutti i frame degli ultimi 100 ms diventa nera invece di un pezzo di immagine fermo.
  const uint32_t now = millis();
  for (uint8_t sy = 0; sy < SCOPE_DEC_SLICES; sy++) {
    if (!(g_black & (1u << sy)) && now - g_sliceOkMs[sy] > SCOPE_VIDEO_BLACK_MS) {
      scopeDecBlackSlice(sy, g_slice[sy]);
      g_black |= (uint16_t)(1u << sy);
    }
  }
  if (!g_active) {
    g_active = true;
    panelTake();
  }
  WirelessScope.videoActive(true);
  #ifdef SCOPE_VIDEO_DEBUG
  d_frames++;
  #endif // SCOPE_VIDEO_DEBUG
  // EN: Display off: the decoding goes on, only the panel stays off. / IT: display spento: la decodifica continua, solo il pannello resta spento.
  if (g_panelOn) rotateFrame();
}

// EN: Slice of the session and frame already chosen. The decoding goes into the work slice and is copied into the
//     frame only if valid: a corrupted slice leaves the previous valid one in its place.
// IT: Fascia della sessione e del frame gia' scelti. La decodifica va nella fascia di lavoro e si copia nel frame
//     solo se valida: una fascia corrotta lascia al suo posto quella valida precedente.
static void handleSlice(const uint8_t *p, uint16_t len) {
  const uint8_t  session = p[SCOPE_VIDEO_B_SESSION];
  const uint16_t frame = pktFrame(p);
  const uint8_t  sy = p[SCOPE_VIDEO_B_SLICE];

  if (session != g_session) {   // newer session: the first frame is accepted even with a lower number / sessione nuova: il primo frame si accetta anche con un numero piu' basso
    if (g_building) showFrame();
    g_session = session;
    g_haveLast = false;
  }
  if (g_building && frame != g_buildFrame) {
    if (!frameNewer(frame, g_buildFrame)) return;
    showFrame();   // a packet of a newer frame: the current one is shown / un pacchetto di un frame piu' nuovo: si mostra quello attuale
  }
  if (!g_building) {
    // EN: From the start of the rotation of a frame, the same or older numbers are discarded.
    // IT: Dall'inizio della rotazione di un frame si scartano i numeri uguali o piu' vecchi.
    if (g_haveLast && !frameNewer(frame, g_lastFrame)) {
      #ifdef SCOPE_VIDEO_DEBUG
      d_old++;
      #endif // SCOPE_VIDEO_DEBUG
      return;
    }
    g_building = true;
    g_buildFrame = frame;
    g_got = 0;
    g_firstMs = millis();
  }
  if (g_got & (1u << sy)) {   // duplicate / duplicato
    #ifdef SCOPE_VIDEO_DEBUG
    d_dup++;
    #endif // SCOPE_VIDEO_DEBUG
    return;
  }
  g_got |= (uint16_t)(1u << sy);

  #ifdef SCOPE_VIDEO_DEBUG
  const uint32_t t0 = micros();
  #endif // SCOPE_VIDEO_DEBUG
  const bool ok = scopeDecSlice(p + OPENFIRE_SCOPE_VIDEO_HEADER, (uint16_t)(len - OPENFIRE_SCOPE_VIDEO_HEADER), sy,
                                p[SCOPE_VIDEO_B_Q], g_work, outPump);
  if (ok) {
    scopeDecCopySlice(g_work, sy, g_slice[sy]);
    g_sliceOkMs[sy] = millis();
    g_black &= (uint16_t)~(1u << sy);
  }
  #ifdef SCOPE_VIDEO_DEBUG
  d_decUs += micros() - t0;
  d_slices++;
  if (!ok) d_corrupt++;
  #endif // SCOPE_VIDEO_DEBUG

  if (g_got == ALL_SLICES) showFrame();
}

// EN: The decoder always takes the packets of the newest frame in the ring: first the newest session and frame among
//     the waiting packets, then only their slices, in order. The others are discarded.
// IT: Il decoder prende sempre i pacchetti del frame piu' recente nell'anello: prima la sessione e il frame piu' recenti
//     tra i pacchetti in attesa, poi solo le loro fasce, in ordine. Gli altri si scartano.
static void processPending(const int *pend, int n) {
  uint8_t S = g_session;
  for (int i = 0; i < n; i++) {
    uint16_t len;
    const uint8_t s = WirelessScope.videoSlot(pend[i], &len)[SCOPE_VIDEO_B_SESSION];
    if (S == 0 || sessionNewer(s, S)) S = s;
  }
  bool haveF = false;
  uint16_t F = 0;
  for (int i = 0; i < n; i++) {
    uint16_t len;
    const uint8_t *p = WirelessScope.videoSlot(pend[i], &len);
    if (p[SCOPE_VIDEO_B_SESSION] != S) continue;
    const uint16_t f = pktFrame(p);
    if (!haveF || frameNewer(f, F)) { F = f; haveF = true; }
  }
  for (int i = 0; i < n; i++) {
    uint16_t len;
    const uint8_t *p = WirelessScope.videoSlot(pend[i], &len);
    if (p[SCOPE_VIDEO_B_SESSION] == S && pktFrame(p) == F) handleSlice(p, len);
    #ifdef SCOPE_VIDEO_DEBUG
    else d_old++;
    #endif // SCOPE_VIDEO_DEBUG
    WirelessScope.videoRelease(pend[i]);
  }
}

// ===================================================================================
// MEMORY / MEMORIA
// ===================================================================================
// EN: Taken at the first pairing, after the sprite of the turning logo has gone: frame in 15 blocks of 7.5 KB (no
//     need of 115 KB in one piece), work slice 7.5 KB, output slices 7.5 KB each (DMA memory). With little memory
//     fewer output slices, at least 2.
// IT: Presa al primo abbinamento, dopo che lo sprite del logo che gira e' sparito: frame in 15 blocchi da 7,5 KB (non
//     servono 115 KB in un pezzo solo), fascia di lavoro 7,5 KB, fasce di uscita da 7,5 KB l'una (memoria DMA).
//     Con poca memoria meno fasce di uscita, almeno 2.
static void freeBuffers(uint8_t nSlice, uint8_t nOut) {
  for (uint8_t i = 0; i < nSlice; i++) heap_caps_free(g_slice[i]);
  for (uint8_t i = 0; i < nOut; i++) heap_caps_free(g_out[i]);
  if (g_work) heap_caps_free(g_work);
  g_work = NULL;
}

static bool allocBuffers() {
  if (g_buffers) return true;
  if (scope_logo_spin_on) return false;   // its 51 KB are not free yet: it tries again / i suoi 51 KB non sono ancora liberi: riprova
  uint8_t ns = 0, n = 0;
  while (ns < SCOPE_DEC_SLICES) {
    g_slice[ns] = (uint16_t *)heap_caps_calloc(SCOPE_DEC_W * SCOPE_DEC_SLICE_H, 2, MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT);   // black / nero
    if (g_slice[ns] == NULL) break;
    ns++;
  }
  g_work = (ns == SCOPE_DEC_SLICES) ? (uint16_t *)heap_caps_malloc(SLICE_BYTES, MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT) : NULL;
  if (g_work != NULL) {
    while (n < SCOPE_VIDEO_OUT_SLICES) {
      g_out[n] = (uint16_t *)heap_caps_malloc(SLICE_BYTES, MALLOC_CAP_DMA | MALLOC_CAP_INTERNAL);
      if (g_out[n] == NULL) break;
      n++;
    }
  }
  if (g_work == NULL || n < 2) {
    freeBuffers(ns, n);
    #ifdef SCOPE_VIDEO_DEBUG
    Serial.printf("VIDEO no memory: free %u, largest block %u\n", (unsigned)heap_caps_get_free_size(MALLOC_CAP_INTERNAL),
                  (unsigned)heap_caps_get_largest_free_block(MALLOC_CAP_INTERNAL));
    #endif // SCOPE_VIDEO_DEBUG
    return false;
  }
  for (int y = 0; y < SCOPE_DEC_H; y++) g_rows[y] = g_slice[y / SCOPE_DEC_SLICE_H] + (y % SCOPE_DEC_SLICE_H) * SCOPE_DEC_W;
  g_nOut = n;
  for (uint8_t i = 0; i < n; i++) g_bufBusy[i] = false;
  g_buffers = true;
  #ifdef SCOPE_VIDEO_DEBUG
  Serial.printf("VIDEO %u output slices, free heap %u\n", n, (unsigned)heap_caps_get_free_size(MALLOC_CAP_INTERNAL));
  #endif // SCOPE_VIDEO_DEBUG
  return true;
}

// ===================================================================================
// VIDEO TASK / TASK DEL VIDEO
// ===================================================================================
static void scopeVideoTask(void *param) {
  static int pend[OPENFIRE_SCOPE_RX_RING];
  bool paired = false;
  uint32_t pairMs = 0;
  #ifdef SCOPE_VIDEO_DEBUG
  uint32_t statMs = millis();
  #endif // SCOPE_VIDEO_DEBUG

  scopeDecInit();

  for (;;) {
    // 1. EN: Link: lost -> the video stops and the state starts from zero at the next pairing.
    //    IT: Collegamento: perso -> il video si ferma e lo stato riparte da zero al prossimo abbinamento.
    if (WirelessScope.link_state < SCOPE_LINK_STATE::SCOPE_LINK_PAIRED) {
      if (paired) {
        paired = false;
        stopVideo();
      }
      g_noSignal = false;
      int s;
      while ((s = WirelessScope.videoReceive(0)) >= 0) WirelessScope.videoRelease(s);
      vTaskDelay(pdMS_TO_TICKS(10));
      continue;
    }
    if (!paired) {
      paired = true;
      pairMs = millis();
      g_session = 0;
      g_haveLast = false;
      g_building = false;
      for (uint8_t sy = 0; sy < SCOPE_DEC_SLICES; sy++) g_sliceOkMs[sy] = pairMs;
      // EN: Nothing of the previous connection on screen: a slice missing from the first frame stays black.
      // IT: Niente della connessione precedente sullo schermo: una fascia mancante dal primo frame resta nera.
      if (g_buffers) {
        for (uint8_t sy = 0; sy < SCOPE_DEC_SLICES; sy++) scopeDecBlackSlice(sy, g_slice[sy]);
        g_black = ALL_SLICES;
      }
    }
    if (!allocBuffers()) {   // memory not free yet (sprite of the logo): it tries again / memoria non ancora libera (sprite del logo): riprova
      int s;
      while ((s = WirelessScope.videoReceive(0)) >= 0) WirelessScope.videoRelease(s);
      vTaskDelay(pdMS_TO_TICKS(100));
      continue;
    }

    // 2. EN: Packets: waits for the first only if the DMA has nothing to do, then takes all those waiting.
    //    IT: Pacchetti: aspetta il primo solo se il DMA non ha niente da fare, poi prende tutti quelli in attesa.
    uint32_t wait = 20;
    if (!outIdle()) wait = 0;
    else if (g_building) {
      const uint32_t el = millis() - g_firstMs;
      wait = (el >= SCOPE_VIDEO_SHOW_MS) ? 0 : SCOPE_VIDEO_SHOW_MS - el;
    }
    int n = 0;
    int s = WirelessScope.videoReceive(wait);
    while (s >= 0) {
      pend[n++] = s;
      if (n == OPENFIRE_SCOPE_RX_RING) break;
      s = WirelessScope.videoReceive(0);
    }
    if (n) processPending(pend, n);
    if (g_owned) outPump();

    // 3. EN: Frame not complete 40 ms after its first slice: shown as it is. / IT: frame non completo 40 ms dopo la sua prima fascia: si mostra cosi' com'e'.
    if (g_building && millis() - g_firstMs >= SCOPE_VIDEO_SHOW_MS) showFrame();

    // 4. EN: Display on/off while the video has the panel; back on, the last frame again.
    //    IT: Display acceso/spento mentre il video ha il pannello; riacceso, di nuovo l'ultimo frame.
    if (g_owned && scope_display_on != g_panelOn) {
      outFlush();
      scopePanelPowerApply();
      if (g_panelOn && g_haveLast) g_refresh = true;
    }
    if (g_refresh && g_owned && g_panelOn && !g_building) {
      g_refresh = false;
      rotateFrame();
    }

    // 5. EN: No signal: 500 ms without slices and no valid pause (500 ms from the last presence with q = 1,
    //        closed by a newer slice). Right after the pairing it waits 500 ms as well.
    //    IT: Nessun segnale: 500 ms senza fasce e nessuna pausa valida (500 ms dall'ultima presenza con q = 1,
    //        chiusa da una fascia piu' nuova). Subito dopo l'abbinamento aspetta anche lui 500 ms.
    //    EN: The times written by the reception are read before now: read after, a newer time would make
    //        now - sliceMs negative (a huge unsigned value) and stop the video in the middle of the game.
    //    IT: I tempi scritti dalla ricezione si leggono prima di now: letti dopo, un tempo piu' nuovo renderebbe
    //        now - sliceMs negativo (un enorme valore senza segno) e fermerebbe il video in pieno gioco.
    const bool sliceSeen = WirelessScope.video_slice_seen;
    const uint32_t sliceMs = WirelessScope.video_slice_ms;
    const bool pauseSeen = WirelessScope.video_pause_seen;
    const uint32_t pauseMs = WirelessScope.video_pause_ms;
    const uint32_t now = millis();
    const bool recent = sliceSeen && (now - sliceMs < SCOPE_VIDEO_NOSIGNAL_MS);
    const bool pause = pauseSeen && (now - pauseMs < OPENFIRE_SCOPE_PAUSE_VALID_MS) &&
                       (!sliceSeen || (int32_t)(pauseMs - sliceMs) >= 0);
    const bool noSignal = !recent && !pause && (now - pairMs >= SCOPE_VIDEO_NOSIGNAL_MS);
    if (noSignal && g_active && !g_building) stopVideo();
    g_noSignal = noSignal && !g_active;
    // EN: No video: the next session is accepted whatever its number. Every START of the PC moves the session of
    //     the dongle forward also without video, and after 128 of them a new session would look older.
    // IT: Nessun video: la sessione successiva si accetta qualunque sia il suo numero. Ogni START del PC fa avanzare
    //     la sessione del dongle anche senza video, e dopo 128 una sessione nuova sembrerebbe piu' vecchia.
    if (noSignal && !g_active && !g_building) g_session = 0;

    #ifdef SCOPE_VIDEO_DEBUG
    if (now - statMs >= 1000) {
      statMs = now;
      Serial.printf("VIDEO fps=%lu slices=%lu corrupt=%lu old=%lu dup=%lu dec=%luus rot=%luus rx=%lu full=%lu bad=%lu out=%u heap=%u\n",
                    (unsigned long)d_frames, (unsigned long)d_slices, (unsigned long)d_corrupt, (unsigned long)d_old, (unsigned long)d_dup,
                    (unsigned long)(d_frames ? d_decUs / d_frames : 0), (unsigned long)(d_frames ? d_rotUs / d_frames : 0),
                    (unsigned long)WirelessScope.video_rx_slices, (unsigned long)WirelessScope.video_rx_full,
                    (unsigned long)WirelessScope.video_rx_bad, g_nOut, (unsigned)heap_caps_get_free_size(MALLOC_CAP_INTERNAL));
      d_frames = d_slices = d_corrupt = d_old = d_dup = d_decUs = d_rotUs = 0;
    }
    #endif // SCOPE_VIDEO_DEBUG
  }
}

void scopeVideoBegin(LGFX *display) {
  g_tft = display;
  g_panelMutex = xSemaphoreCreateMutex();
  xTaskCreate(scopeVideoTask, "scopeVideo", 6144, NULL, SCOPE_VIDEO_TASK_PRIO, NULL);
}

#endif // OPENFIRE_WIRELESS_ENABLE && ARDUINO_ARCH_ESP32 && SCOPE && USES_DISPLAY && USE_LOVYAN_GFX
