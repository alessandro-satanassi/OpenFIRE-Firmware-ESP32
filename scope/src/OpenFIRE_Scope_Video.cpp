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

#define BAND_BYTES    (SCOPE_DEC_W * SCOPE_ROT_BAND_H * 2)    // 3840: output band / fascia di uscita
#define ALL_SLICES    ((1u << SCOPE_DEC_SLICES) - 1)

static LGFX             *g_tft = NULL;
static SemaphoreHandle_t g_panelMutex = NULL;
static bool              g_panelOn = true;        // power applied to the panel, under the lock / alimentazione data al pannello, sotto il blocco
static volatile bool     g_noSignal = false;

// EN: Buffers: decoded frame in YCbCr (fixed array with the coded MCU, and the table of its MCU), output bands of
//     8 rows of the rotation (DMA).
//     The array is shared with the sprite of the turning logo: the logo turns only during the search, the frame is
//     used only while paired. Whoever wants it first sets its own flag and then looks at the other's: if the other
//     has it, it takes its own flag back and tries again later (one core: the two flags are never both kept).
// IT: Buffer: frame decodificato in YCbCr (array fisso con gli MCU codificati, e la tabella dei suoi MCU), fasce di
//     uscita da 8 righe della rotazione (DMA).
//     L'array e' condiviso con lo sprite del logo che gira: il logo gira solo durante la ricerca, il frame si usa solo
//     da abbinati. Chi lo vuole prima alza il suo flag e poi guarda quello dell'altro: se l'altro lo ha, riabbassa il
//     suo e riprova dopo (un core solo: i due flag non restano mai alzati insieme).
static_assert(SCOPE_VIDEO_SHARED_BYTES == SCOPE_DEC_MCU_CODED * SCOPE_DEC_MCU_BYTES, "frame size / misura del frame");
alignas(4) static uint8_t g_shared[SCOPE_VIDEO_SHARED_BYTES];
static volatile bool g_frameInUse = false;   // the video has the array / il video ha l'array
static bool      g_buffers = false;
static uint8_t  *g_mcu[SCOPE_DEC_MCU_TABLE];
static uint16_t *g_out[SCOPE_VIDEO_OUT_SLICES];
static uint8_t   g_nOut = 0;

#ifdef SCOPE_VIDEO_DEBUG
// EN: Counters of the last video, from its start to its end (shown when the video is over).
// IT: Contatori dell'ultimo video, dal suo inizio alla sua fine (mostrati quando il video e' finito).
static uint32_t d_frames = 0, d_slices = 0, d_corrupt = 0, d_old = 0, d_dup = 0, d_black = 0, d_decUs = 0, d_rotUs = 0, d_rotMax = 0;
static uint32_t d_startMs = 0, d_endMs = 0, d_rx0 = 0, d_full0 = 0, d_bad0 = 0, d_rx1 = 0, d_full1 = 0, d_bad1 = 0;
static char     d_mem[32] = "MEM ?";   // state of the memory of the video / stato della memoria del video
#endif // SCOPE_VIDEO_DEBUG

// ===================================================================================
// PANEL LOCK / BLOCCO DEL PANNELLO
// ===================================================================================
static bool g_owned = false;   // the video has the panel / il video ha il pannello

bool scopePanelTryLock() {
  #ifdef SCOPE_VIDEO_DEBUG
  // EN: The debug lines take the panel for a moment: the state screens wait for them, otherwise they would take them for
  //     the video and draw everything again (the lines would flash).
  // IT: Le righe di debug prendono il pannello per un attimo: le schermate di stato le aspettano, altrimenti le
  //     scambierebbero per il video e ridisegnerebbero tutto (le righe lampeggerebbero).
  if (g_panelMutex != NULL && !g_owned) return xSemaphoreTake(g_panelMutex, pdMS_TO_TICKS(50)) == pdTRUE;
  #endif // SCOPE_VIDEO_DEBUG
  return g_panelMutex != NULL && xSemaphoreTake(g_panelMutex, 0) == pdTRUE;
}
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
// EN: The rotated bands wait in a queue; outPump() starts the next one as soon as the DMA is free, without waiting.
//     It is called after every MCU of the decoder and between the packets, so the DMA keeps working while the
//     next frame is decoded. A buffer is free again when its DMA has finished.
// IT: Le fasce ruotate aspettano in una coda; outPump() fa partire la successiva appena il DMA e' libero, senza aspettare.
//     Si chiama dopo ogni MCU del decoder e tra un pacchetto e l'altro, cosi' il DMA continua a lavorare mentre si
//     decodifica il frame successivo. Un buffer torna libero quando il suo DMA ha finito.
static bool    g_bufBusy[SCOPE_VIDEO_OUT_SLICES];
static uint8_t g_qBand[SCOPE_VIDEO_OUT_SLICES], g_qBuf[SCOPE_VIDEO_OUT_SLICES];
static uint8_t g_qHead = 0, g_qCount = 0;
static int     g_inflight = -1;   // buffer on the DMA / buffer sul DMA

static void outPump() {
  if (g_inflight >= 0 && !g_tft->dmaBusy()) { g_bufBusy[g_inflight] = false; g_inflight = -1; }
  if (g_qCount == 0 || g_inflight >= 0) return;
  const uint8_t band = g_qBand[g_qHead], b = g_qBuf[g_qHead];
  g_qHead = (uint8_t)((g_qHead + 1) % g_nOut);
  g_qCount--;
  // EN: Data already in the byte order of the display (swap565): straight to the DMA, no conversion.
  // IT: Dati gia' nell'ordine di byte del display (swap565): dritti al DMA, nessuna conversione.
  g_tft->pushImageDMA(scopeRotBandX(band), band * SCOPE_ROT_BAND_H, scopeRotBandW(band), SCOPE_ROT_BAND_H,
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
  for (uint8_t band = 0; band < SCOPE_ROT_BANDS; band++) {
    const uint8_t b = outFreeBuffer();
    g_bufBusy[b] = true;
    scopeRotBandYcc(g_mcu, band, c, s, g_out[b]);
    const uint8_t tail = (uint8_t)((g_qHead + g_qCount) % g_nOut);
    g_qBand[tail] = band;
    g_qBuf[tail] = b;
    g_qCount++;
    outPump();
  }
  #ifdef SCOPE_VIDEO_DEBUG
  const uint32_t rt = micros() - t0;
  d_rotUs += rt;
  if (rt > d_rotMax) d_rotMax = rt;
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
  #ifdef SCOPE_VIDEO_DEBUG
  d_endMs = millis();
  d_rx1 = WirelessScope.video_rx_slices;
  d_full1 = WirelessScope.video_rx_full;
  d_bad1 = WirelessScope.video_rx_bad;
  #endif // SCOPE_VIDEO_DEBUG
  WirelessScope.videoActive(false);
  panelGive();
}

static void showFrame() {
  g_building = false;
  g_lastFrame = g_buildFrame;
  g_haveLast = true;
  // EN: A slice missing from every frame of the last SCOPE_VIDEO_BLACK_MS goes black instead of an image piece standing still.
  // IT: Una fascia mancante da tutti i frame degli ultimi SCOPE_VIDEO_BLACK_MS diventa nera invece di un pezzo di immagine fermo.
  const uint32_t now = millis();
  for (uint8_t sy = 0; sy < SCOPE_DEC_SLICES; sy++) {
    if (!(g_black & (1u << sy)) && now - g_sliceOkMs[sy] > SCOPE_VIDEO_BLACK_MS) {
      scopeDecBlackSliceYcc(sy, g_mcu);
      g_black |= (uint16_t)(1u << sy);
      #ifdef SCOPE_VIDEO_DEBUG
      if (g_active) d_black++;
      #endif // SCOPE_VIDEO_DEBUG
    }
  }
  if (!g_active) {
    #ifdef SCOPE_VIDEO_DEBUG
    d_frames = d_slices = d_corrupt = d_old = d_dup = d_black = d_decUs = d_rotUs = d_rotMax = 0;
    d_startMs = millis();
    d_rx0 = WirelessScope.video_rx_slices;
    d_full0 = WirelessScope.video_rx_full;
    d_bad0 = WirelessScope.video_rx_bad;
    #endif // SCOPE_VIDEO_DEBUG
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

// EN: Slice of the session and frame already chosen. The decoding goes straight into the frame: a corrupted slice
//     (practically never: the Wi-Fi discards the damaged packets) becomes black.
// IT: Fascia della sessione e del frame gia' scelti. La decodifica va direttamente nel frame: una fascia corrotta
//     (praticamente mai: il Wi-Fi scarta i pacchetti danneggiati) diventa nera.
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
  const bool ok = scopeDecSliceYcc(p + OPENFIRE_SCOPE_VIDEO_HEADER, (uint16_t)(len - OPENFIRE_SCOPE_VIDEO_HEADER), sy,
                                   p[SCOPE_VIDEO_B_Q], g_mcu, outPump);
  if (ok) {
    g_sliceOkMs[sy] = millis();
    g_black &= (uint16_t)~(1u << sy);
  } else {
    scopeDecBlackSliceYcc(sy, g_mcu);
    g_black |= (uint16_t)(1u << sy);
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
// EN: Taken at the first pairing: table of the frame in the fixed array, then the output bands of 3.75 KB (DMA
//     memory): 2 always, the others (up to SCOPE_VIDEO_OUT_SLICES) only while SCOPE_VIDEO_HEAP_MARGIN stays free for
//     the Wi-Fi. They are kept until the reset.
// IT: Presa al primo abbinamento: tabella del frame nell'array fisso, poi le fasce di uscita da 3,75 KB (memoria DMA):
//     2 sempre, le altre (fino a SCOPE_VIDEO_OUT_SLICES) solo se restano liberi SCOPE_VIDEO_HEAP_MARGIN per il Wi-Fi.
//     Restano fino al reset.
static bool allocBuffers() {
  if (g_buffers) return true;
  if (scopeDecMcuCount() > SCOPE_DEC_MCU_CODED) {   // never with the 240x240 circle / mai col cerchio 240x240
    #ifdef SCOPE_VIDEO_DEBUG
    snprintf(d_mem, sizeof(d_mem), "NO MEM MCU %u", (unsigned)scopeDecMcuCount());
    #endif // SCOPE_VIDEO_DEBUG
    return false;
  }
  uint8_t n = 0;
  while (n < SCOPE_VIDEO_OUT_SLICES) {
    if (n >= 2 && heap_caps_get_free_size(MALLOC_CAP_INTERNAL) < (size_t)SCOPE_VIDEO_HEAP_MARGIN + BAND_BYTES) break;
    g_out[n] = (uint16_t *)heap_caps_malloc(BAND_BYTES, MALLOC_CAP_DMA | MALLOC_CAP_INTERNAL);
    if (g_out[n] == NULL) break;
    n++;
  }
  if (n < 2) {
    for (uint8_t i = 0; i < n; i++) heap_caps_free(g_out[i]);
    #ifdef SCOPE_VIDEO_DEBUG
    snprintf(d_mem, sizeof(d_mem), "NO MEM F%u B%u", (unsigned)(heap_caps_get_free_size(MALLOC_CAP_INTERNAL) / 1024),
             (unsigned)(heap_caps_get_largest_free_block(MALLOC_CAP_INTERNAL) / 1024));
    #endif // SCOPE_VIDEO_DEBUG
    return false;
  }
  uint8_t *const piece = g_shared;
  const uint8_t count = (uint8_t)scopeDecMcuCount();
  scopeDecMcuTable(&piece, &count, 1, g_mcu);
  g_nOut = n;
  for (uint8_t i = 0; i < n; i++) g_bufBusy[i] = false;
  g_buffers = true;
  return true;
}

// EN: Array of the frame: the video takes it at the pairing and gives it back when the link is lost.
// IT: Array del frame: il video lo prende all'abbinamento e lo restituisce quando il collegamento si perde.
static bool frameTake() {
  g_frameInUse = true;
  __sync_synchronize();
  if (scope_logo_spin_on) {   // the logo is turning in it / il logo ci sta girando
    g_frameInUse = false;
    return false;
  }
  return true;
}

static void frameGive() {
  __sync_synchronize();
  g_frameInUse = false;
}

uint8_t *scopeLogoMemoryTake() {
  scope_logo_spin_on = true;
  __sync_synchronize();
  if (g_frameInUse) {   // the video has it / ce l'ha il video
    scope_logo_spin_on = false;
    return NULL;
  }
  return g_shared;
}

void scopeLogoMemoryGive() {
  __sync_synchronize();
  scope_logo_spin_on = false;
}

#ifdef SCOPE_VIDEO_DEBUG
// ===================================================================================
// DEBUG ON THE DISPLAY / DEBUG SUL DISPLAY
// ===================================================================================
// EN: Every second: 7 lines in large text over the logo, drawn only while the video does not have the panel and the
//     logo is not turning, with the lock of the panel and without waiting for it. The counters start again every second
//     (rx, full and bad count from the power on).
// IT: Ogni secondo: 7 righe a caratteri grandi sopra il logo, disegnate solo mentre il video non ha il pannello e il
//     logo non gira, con il blocco del pannello e senza aspettarlo. I contatori ripartono ogni secondo
//     (rx, full e bad contano dall'accensione).
static void debugStats() {
  static uint32_t statMs = 0;
  const uint32_t now = millis();
  if (now - statMs < 1000) return;
  statMs = now;
  char l[6][24];
  // EN: last video: average FPS and frames, then per frame the average decoding and rotation (us) and the longest
  //     rotation (RM); slices received, lost with the ring full and discarded during the video.
  // IT: ultimo video: FPS medi e frame, poi per frame la decodifica e la rotazione medie (us) e la rotazione piu'
  //     lunga (RM); fasce ricevute, perse con l'anello pieno e scartate durante il video.
  const uint32_t endMs = g_active ? now : d_endMs;
  const uint32_t ms = endMs - d_startMs;
  const uint32_t rx = (g_active ? WirelessScope.video_rx_slices : d_rx1) - d_rx0;
  snprintf(l[0], sizeof(l[0]), "FPS %lu FR %lu", (unsigned long)(ms ? (uint64_t)d_frames * 1000 / ms : 0), (unsigned long)d_frames);
  snprintf(l[1], sizeof(l[1]), "COR %lu OLD %lu", (unsigned long)d_corrupt, (unsigned long)d_old);
  snprintf(l[2], sizeof(l[2]), "DUP %lu BLK %lu", (unsigned long)d_dup, (unsigned long)d_black);
  snprintf(l[3], sizeof(l[3]), "DEC %lu ROT %lu", (unsigned long)(d_frames ? d_decUs / d_frames : 0), (unsigned long)(d_frames ? d_rotUs / d_frames : 0));
  snprintf(l[4], sizeof(l[4]), "RM %lu RX %lu", (unsigned long)d_rotMax, (unsigned long)rx);
  snprintf(l[5], sizeof(l[5]), "FULL %lu BAD %lu", (unsigned long)((g_active ? WirelessScope.video_rx_full : d_full1) - d_full0),
           (unsigned long)((g_active ? WirelessScope.video_rx_bad : d_bad1) - d_bad0));
  // EN: buffers taken: output bands and free memory now (KB). / IT: buffer presi: fasce di uscita e memoria libera ora (KB).
  if (g_buffers) snprintf(d_mem, sizeof(d_mem), "MEM OK O%u H%u", g_nOut, (unsigned)(heap_caps_get_free_size(MALLOC_CAP_INTERNAL) / 1024));

  if (g_owned || scope_logo_spin_on || !scopePanelTryLock()) return;
  if (g_panelOn) {
    g_tft->fillRect(0, 43, 240, 155, TFT_BLACK);   // between "NO SIGNAL" and the state line / tra "NO SIGNAL" e la riga di stato
    g_tft->setTextSize(2);
    g_tft->setTextDatum(textdatum_t::middle_center);
    g_tft->setTextColor(TFT_YELLOW, TFT_BLACK);
    g_tft->drawString(d_mem, 120, 54);
    g_tft->setTextColor(TFT_WHITE, TFT_BLACK);
    for (uint8_t i = 0; i < 6; i++) g_tft->drawString(l[i], 120, 76 + 22 * i);
  }
  scopePanelUnlock();
}
#endif // SCOPE_VIDEO_DEBUG

// ===================================================================================
// VIDEO TASK / TASK DEL VIDEO
// ===================================================================================
static void scopeVideoTask(void *param) {
  static int pend[OPENFIRE_SCOPE_RX_PACKETS];
  bool paired = false;
  uint32_t pairMs = 0;

  scopeDecInit();

  for (;;) {
    // 1. EN: Link: lost -> the video stops and the state starts from zero at the next pairing.
    //    IT: Collegamento: perso -> il video si ferma e lo stato riparte da zero al prossimo abbinamento.
    if (WirelessScope.link_state < SCOPE_LINK_STATE::SCOPE_LINK_PAIRED) {
      if (paired) {
        paired = false;
        frameGive();   // first: the logo can take it as soon as it has the panel / prima: il logo puo' prenderlo appena ha il pannello
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
    }
    if (!g_frameInUse) {
      if (!allocBuffers() || !frameTake()) {   // no bands or the logo still turning in the array: it tries again / niente fasce o il logo gira ancora nell'array: riprova
        #ifdef SCOPE_VIDEO_DEBUG
        debugStats();
        #endif // SCOPE_VIDEO_DEBUG
        int s;
        while ((s = WirelessScope.videoReceive(0)) >= 0) WirelessScope.videoRelease(s);
        vTaskDelay(pdMS_TO_TICKS(100));
        continue;
      }
      // EN: Array just taken back from the logo, nothing of the previous connection: a slice missing from the first
      //     frame stays black.
      // IT: Array appena ripreso dal logo, niente della connessione precedente: una fascia mancante dal primo frame
      //     resta nera.
      for (uint8_t sy = 0; sy < SCOPE_DEC_SLICES; sy++) scopeDecBlackSliceYcc(sy, g_mcu);
      g_black = ALL_SLICES;
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
      if (n == OPENFIRE_SCOPE_RX_PACKETS) break;
      s = WirelessScope.videoReceive(0);
    }
    if (n) processPending(pend, n);
    if (g_owned) outPump();

    // 3. EN: Frame not complete SCOPE_VIDEO_SHOW_MS after its first slice: shown as it is. / IT: frame non completo SCOPE_VIDEO_SHOW_MS dopo la sua prima fascia: si mostra cosi' com'e'.
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
    debugStats();
    #endif // SCOPE_VIDEO_DEBUG
  }
}

void scopeVideoBegin(LGFX *display) {
  g_tft = display;
  g_panelMutex = xSemaphoreCreateMutex();
  xTaskCreate(scopeVideoTask, "scopeVideo", 6144, NULL, SCOPE_VIDEO_TASK_PRIO, NULL);
}

#endif // OPENFIRE_WIRELESS_ENABLE && ARDUINO_ARCH_ESP32 && SCOPE && USES_DISPLAY && USE_LOVYAN_GFX
