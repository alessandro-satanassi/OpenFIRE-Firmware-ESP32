/*!
 * @file OpenFIRE_Scope_Dongle_USB.cpp
 * @brief USB link between the PC daemon and the scope dongle
 * @n CPP USB link between the PC daemon and the scope dongle
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#include <Arduino.h>
#include "OpenFIRE_Scope_Dongle_USB.h"

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32) && defined(SCOPE_DONGLE)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "freertos/semphr.h"
#include "OpenFIRE_Wireless_Scope.h"
#include "OpenFIRE-SCOPE-DONGLE-version.h"

#define USB_SYNC_0            0xA5
#define USB_SYNC_1            0x5A
#define USB_TYPE_VIDEO        OPENFIRE_SCOPE_VIDEO_MARKER   // 0x4D 'M'
#define USB_TYPE_CONTROL      0x43                          // 'C'
#define USB_TYPE_TEST         0x54                          // 'T'
#define USB_MAX_CONTENT       OPENFIRE_SCOPE_VIDEO_MAX_PACKET   // 1408
#define USB_RX_BUFFER         8192   // ~8 ms of USB at 1 MB/s / circa 8 ms di USB a 1 MB/s
#define USB_STAT_EVERY_MS     500
#define USB_TEXT_MAX          96

static SemaphoreHandle_t g_usbTxMutex = NULL;
static volatile uint32_t g_usbRxBytes = 0;   // OF_STAT RX: every byte received, bandwidth test included / ogni byte ricevuto, prova di banda compresa

// ===================================================================================
// TRANSMISSION TO THE PC / TRASMISSIONE AL PC
// ===================================================================================
static void usbSendContent(const uint8_t *content, uint16_t len) {
  const uint8_t head[4] = { USB_SYNC_0, USB_SYNC_1, (uint8_t)(len & 0xFF), (uint8_t)(len >> 8) };
  xSemaphoreTake(g_usbTxMutex, portMAX_DELAY);   // whole frames, never mixed / trame intere, mai mescolate
  // EN: PC connected but without room for the whole frame (it is not reading): the frame is skipped, never cut and
  //     never waited for. PC not connected: written as before, the core keeps the newest bytes.
  // IT: PC collegato ma senza posto per tutta la trama (non sta leggendo): la trama si salta, mai tagliata e mai
  //     aspettata. PC non collegato: scritta come prima, il core tiene i byte piu' nuovi.
  if (!Serial || Serial.availableForWrite() >= (int)(sizeof(head) + len)) {
    Serial.write(head, sizeof(head));
    Serial.write(content, len);
  }
  xSemaphoreGive(g_usbTxMutex);
}

static void usbSendControl(const char *text) {
  uint8_t buf[1 + USB_TEXT_MAX];
  size_t n = strlen(text);
  if (n > USB_TEXT_MAX) n = USB_TEXT_MAX;
  buf[0] = USB_TYPE_CONTROL;
  memcpy(&buf[1], text, n);
  usbSendContent(buf, (uint16_t)(n + 1));
}

// EN: State for WHOAMI and for the messages sent by itself. / IT: stato per WHOAMI e per i messaggi mandati da solo.
static void usbStateText(char *out, size_t size) {
  const uint8_t state = WirelessScope.link_state;
  if (state == SCOPE_LINK_STATE::SCOPE_LINK_WAIT) {   // channel search / ricerca del canale
    snprintf(out, size, "OF_SCOPE:BOOT");
  } else if (state >= SCOPE_LINK_STATE::SCOPE_LINK_PAIRED) {
    snprintf(out, size, "OF_SCOPE:P%u;FW=%u.%u.%u;CH=%u", WirelessScope.player,
             OPENFIRE_SCOPE_DONGLE_VERSION_MAJOR, OPENFIRE_SCOPE_DONGLE_VERSION_MINOR, OPENFIRE_SCOPE_DONGLE_VERSION_PATCH,
             WirelessScope.channel);
  } else {
    snprintf(out, size, "OF_SCOPE:UNPAIRED;FW=%u.%u.%u;CH=%u",
             OPENFIRE_SCOPE_DONGLE_VERSION_MAJOR, OPENFIRE_SCOPE_DONGLE_VERSION_MINOR, OPENFIRE_SCOPE_DONGLE_VERSION_PATCH,
             WirelessScope.channel);
  }
}

// ===================================================================================
// CONTROL MESSAGES FROM THE PC / MESSAGGI DI CONTROLLO DAL PC
// ===================================================================================
// EN: text starts with cmd, followed by nothing or by ';' and optional fields. / IT: il testo inizia con cmd, seguito da niente o da ';' e campi facoltativi.
static bool usbTextIs(const char *text, uint16_t len, const char *cmd) {
  const size_t l = strlen(cmd);
  return len >= l && memcmp(text, cmd, l) == 0 && (len == l || text[l] == ';');
}

// EN: SKIP=36,40 -> bit i = 1 for the channel i of OPENFIRE_SCOPE_CHANNELS. / IT: SKIP=36,40 -> bit i = 1 per il canale i di OPENFIRE_SCOPE_CHANNELS.
static uint8_t usbSkipMask(const char *text) {
  static const uint8_t channels[OPENFIRE_SCOPE_N_CHANNELS] = OPENFIRE_SCOPE_CHANNELS;
  uint8_t mask = 0;
  const char *p = strstr(text, "SKIP=");
  if (p == NULL) return 0;
  p += 5;
  while (*p) {
    char *end;
    const long ch = strtol(p, &end, 10);
    if (end == p) break;
    for (uint8_t i = 0; i < OPENFIRE_SCOPE_N_CHANNELS; i++) if (channels[i] == ch) mask |= (1 << i);
    p = end;
    if (*p != ',') break;
    p++;
  }
  return mask;
}

static void usbHandleControl(const char *data, uint16_t len) {
  char text[USB_TEXT_MAX + 1];
  if (len > USB_TEXT_MAX) len = USB_TEXT_MAX;
  memcpy(text, data, len);
  text[len] = '\0';

  if (usbTextIs(text, len, "WHOAMI")) {
    char answer[64];
    usbStateText(answer, sizeof(answer));
    usbSendControl(answer);
  } else if (usbTextIs(text, len, "START")) {   // video of the slot starts: new session / parte il video dello slot: sessione nuova
    WirelessScope.videoStart();
  } else if (usbTextIs(text, len, "PAUSE")) {   // game paused, every 200 ms / gioco in pausa, ogni 200 ms
    WirelessScope.videoPause();
  } else if (usbTextIs(text, len, "RESCAN")) {  // ignored with a scope paired / ignorato con un mirino abbinato
    WirelessScope.requestRescan(usbSkipMask(text));
  }
  // EN: unknown messages are ignored. / IT: i messaggi sconosciuti si ignorano.
}

static void usbHandleContent(const uint8_t *content, uint16_t len) {
  switch (content[0]) {
    case USB_TYPE_VIDEO:   WirelessScope.videoPush(content, len); break;   // straight to the radio queue / dritto nella coda della radio
    case USB_TYPE_CONTROL: usbHandleControl((const char *)&content[1], len - 1); break;
    case USB_TYPE_TEST:    break;   // counted in RX, never on the air / contato in RX, mai in aria
    default:               break;
  }
}

static inline bool usbTypeValid(uint8_t type) {
  return type == USB_TYPE_VIDEO || type == USB_TYPE_CONTROL || type == USB_TYPE_TEST;
}

// ===================================================================================
// USB TASK / TASK DELLA USB
// ===================================================================================
static void scopeUsbTask(void *param) {
  static uint8_t content[USB_MAX_CONTENT];
  uint8_t chunk[512];
  enum { WAIT_SYNC_0, WAIT_SYNC_1, WAIT_LEN_0, WAIT_LEN_1, WAIT_DATA } st = WAIT_SYNC_0;
  uint16_t len = 0;
  uint16_t pos = 0;

  char lastState[64] = "";
  uint8_t lastRescan = WirelessScope.rescan_count;
  uint32_t lastStatMs = millis();

  for (;;) {
    // 1. EN: Reception: every packet is passed on as soon as it is complete, without waiting for the end of the frame.
    //    IT: Ricezione: ogni pacchetto si passa avanti appena e' completo, senza aspettare la fine del frame.
    int n = Serial.available();
    if (n > 0) {
      if (n > (int)sizeof(chunk)) n = sizeof(chunk);
      n = Serial.read(chunk, n);
      if (n > 0) g_usbRxBytes = g_usbRxBytes + n;
      for (int i = 0; i < n; i++) {
        const uint8_t b = chunk[i];
        switch (st) {
          case WAIT_SYNC_0: if (b == USB_SYNC_0) st = WAIT_SYNC_1; break;
          case WAIT_SYNC_1: st = (b == USB_SYNC_1) ? WAIT_LEN_0 : (b == USB_SYNC_0 ? WAIT_SYNC_1 : WAIT_SYNC_0); break;
          case WAIT_LEN_0:  len = b; st = WAIT_LEN_1; break;
          case WAIT_LEN_1:
            len |= (uint16_t)b << 8;
            pos = 0;
            // EN: length 0 or over 1408: not a frame, search the next sync. / IT: lunghezza 0 o oltre 1408: non e' una trama, si cerca il sincronismo successivo.
            st = (len == 0 || len > USB_MAX_CONTENT) ? WAIT_SYNC_0 : WAIT_DATA;
            break;
          case WAIT_DATA: {
            // EN: Unknown type: discarded at once, and this byte may be the start of the next sync.
            // IT: Tipo sconosciuto: scartato subito, e questo byte puo' essere l'inizio del sincronismo successivo.
            if (pos == 0 && !usbTypeValid(b)) {
              st = (b == USB_SYNC_0) ? WAIT_SYNC_1 : WAIT_SYNC_0;
              break;
            }
            uint16_t take = len - pos;
            if (take > (uint16_t)(n - i)) take = (uint16_t)(n - i);
            memcpy(&content[pos], &chunk[i], take);
            pos += take;
            i += take - 1;
            if (pos == len) {
              usbHandleContent(content, len);
              st = WAIT_SYNC_0;
            }
            break;
          }
        }
      }
    }

    // 2. EN: State to the PC by itself at every change (UNPAIRED with the channel, Pn at the pairing), and after a RESCAN.
    //    IT: Stato al PC da solo a ogni cambio (UNPAIRED con il canale, Pn all'abbinamento), e dopo un RESCAN.
    if (WirelessScope.link_state != SCOPE_LINK_STATE::SCOPE_LINK_WAIT) {
      char state[64];
      usbStateText(state, sizeof(state));
      const uint8_t rescan = WirelessScope.rescan_count;
      if (strcmp(state, lastState) != 0 || rescan != lastRescan) {
        strcpy(lastState, state);
        lastRescan = rescan;
        usbSendControl(state);
      }
    }

    // 3. EN: Statistics every 500 ms (counters from power on). / IT: statistiche ogni 500 ms (contatori dall'accensione).
    const uint32_t now = millis();
    if (now - lastStatMs >= USB_STAT_EVERY_MS) {
      lastStatMs = now;
      char stat[USB_TEXT_MAX];
      snprintf(stat, sizeof(stat), "OF_STAT;RX=%lu;TX=%lu;DROP=%lu;QFULL=%lu",
               (unsigned long)g_usbRxBytes, (unsigned long)WirelessScope.video_tx_bytes,
               (unsigned long)WirelessScope.video_drop, (unsigned long)WirelessScope.video_qfull);
      usbSendControl(stat);
    }

    if (n <= 0) vTaskDelay(pdMS_TO_TICKS(1));
  }
}

void scopeUsbBegin() {
  Serial.setRxBufferSize(USB_RX_BUFFER);   // before begin() / prima di begin()
  Serial.setTxTimeoutMs(0);                // never waits for a PC that does not read / non aspetta mai un PC che non legge
  Serial.begin();
  g_usbTxMutex = xSemaphoreCreateMutex();
  xTaskCreate(scopeUsbTask, "scopeUsb", 4096, NULL, 5, NULL);   // below the video task (6) / sotto il task del video (6)
}

#endif // OPENFIRE_WIRELESS_ENABLE && ARDUINO_ARCH_ESP32 && SCOPE_DONGLE
