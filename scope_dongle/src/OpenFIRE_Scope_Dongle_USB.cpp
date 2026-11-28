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
#include "driver/usb_serial_jtag.h"
#include "hal/usb_serial_jtag_ll.h"
#include "esp_rom_sys.h"
#include "OpenFIRE_Wireless_Scope.h"
#include "OpenFIRE-SCOPE-DONGLE-version.h"

#define USB_SYNC_0            0xA5
#define USB_SYNC_1            0x5A
#define USB_TYPE_VIDEO        OPENFIRE_SCOPE_VIDEO_MARKER   // 0x4D 'M'
#define USB_TYPE_CONTROL      0x43                          // 'C'
#define USB_TYPE_TEST         0x54                          // 'T'
#define USB_MAX_CONTENT       OPENFIRE_SCOPE_VIDEO_MAX_PACKET   // 1408
#define USB_RX_BUFFER         8192   // ~8 ms of USB at 1 MB/s / circa 8 ms di USB a 1 MB/s
#define USB_TX_BUFFER         256    // control messages only / solo messaggi di controllo
#define USB_READ_CHUNK        1024
#define USB_READ_WAIT_MS      10     // the task wakes at least this often for the state and the statistics / il task si sveglia almeno cosi' spesso per lo stato e le statistiche
#define USB_STAT_EVERY_MS     500
#define USB_TEXT_MAX          96
#define USB_RESYNC_MS         100    // frame left half way for so long: back to the sync / trama rimasta a meta' per tanto: si torna al sincronismo
#define USB_STUCK_FLUSH_MS    50     // IN buffer busy so long while the PC asks for data: the flush again / buffer IN occupato da tanto mentre il PC chiede dati: di nuovo il flush
#define USB_STUCK_RESET_MS    100    // still busy so long after that flush: the IN FIFO is reset / ancora occupato da tanto dopo quel flush: si azzera la FIFO IN
#define USB_STUCK_RETRY_MS    1000   // still busy after the reset: everything again after so long / ancora occupato dopo l'azzeramento: tutto di nuovo dopo tanto

static SemaphoreHandle_t g_usbTxMutex = NULL;
static volatile uint32_t g_usbRxBytes = 0;   // OF_STAT RX: every byte received, bandwidth test included / ogni byte ricevuto, prova di banda compresa
static portMUX_TYPE g_usbKickMux = portMUX_INITIALIZER_UNLOCKED;
#ifdef SCOPE_DONGLE_USB_DEBUG
static volatile uint16_t d_whoami = 0, d_answerLost = 0;
static volatile uint8_t d_flush = 0, d_reset = 0, d_busReset = 0;
#endif // SCOPE_DONGLE_USB_DEBUG

// ===================================================================================
// TRANSMISSION TO THE PC / TRASMISSIONE AL PC
// ===================================================================================
// EN: Returns false when the frame did not go out whole (no room: the PC is not reading). Only control messages go
//     to the PC: header and content in one buffer, one write into the ring buffer of the ESP-IDF driver, which takes
//     the frame whole or not at all (never cut, never waited for).
// IT: Restituisce false se la trama non e' uscita intera (niente posto: il PC non sta leggendo). Al PC vanno solo
//     messaggi di controllo: intestazione e contenuto in un buffer, una scrittura nel buffer circolare del driver
//     ESP-IDF, che prende la trama intera o niente (mai tagliata, mai aspettata).
static bool usbSendContent(const uint8_t *content, uint16_t len) {
  uint8_t frame[4 + 1 + USB_TEXT_MAX];
  if (len > sizeof(frame) - 4) return false;
  frame[0] = USB_SYNC_0;
  frame[1] = USB_SYNC_1;
  frame[2] = (uint8_t)(len & 0xFF);
  frame[3] = (uint8_t)(len >> 8);
  memcpy(&frame[4], content, len);
  const size_t total = 4 + len;

  xSemaphoreTake(g_usbTxMutex, portMAX_DELAY);   // whole frames, never mixed / trame intere, mai mescolate
  const bool sent = usb_serial_jtag_write_bytes(frame, total, 0) == (int)total;
  // EN: Kick: the driver sends only when the host has taken the previous packet (IN_EMPTY). After a reset without
  //     a USB bus reset (e.g. after the upload) that signal may never come, with the IN buffer empty: nothing would
  //     ever go out and the PC would get no answer to WHOAMI. A free IN buffer here, with the interrupt just
  //     enabled, means exactly that case: an empty packet makes the host pick it up and the signal comes back.
  //     Otherwise the buffer is busy and nothing is done. With interrupts off, the ISR cannot fall in between.
  // IT: Spinta: il driver invia solo quando l'host ha preso il pacchetto precedente (IN_EMPTY). Dopo un reset senza
  //     reset del bus USB (per esempio dopo il caricamento) quel segnale puo' non arrivare mai, con il buffer IN
  //     vuoto: non uscirebbe piu' niente e il PC non avrebbe risposta a WHOAMI. Un buffer IN libero qui, con
  //     l'interrupt appena abilitato, e' proprio quel caso: un pacchetto vuoto lo fa prendere all'host e il segnale
  //     torna. Altrimenti il buffer e' occupato e non si fa niente. Con gli interrupt spenti l'ISR non puo' cadere in mezzo.
  portENTER_CRITICAL(&g_usbKickMux);
  if (usb_serial_jtag_ll_txfifo_writable()) usb_serial_jtag_ll_txfifo_flush();
  portEXIT_CRITICAL(&g_usbKickMux);
  xSemaphoreGive(g_usbTxMutex);
  return sent;
}

static bool usbSendControl(const char *text) {
  uint8_t buf[1 + USB_TEXT_MAX];
  size_t n = strlen(text);
  if (n > USB_TEXT_MAX) n = USB_TEXT_MAX;
  buf[0] = USB_TYPE_CONTROL;
  memcpy(&buf[1], text, n);
  return usbSendContent(buf, (uint16_t)(n + 1));
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
    #ifdef SCOPE_DONGLE_USB_DEBUG
    d_whoami = d_whoami + 1;
    if (!usbSendControl(answer)) d_answerLost = d_answerLost + 1;
    #else
    usbSendControl(answer);
    #endif // SCOPE_DONGLE_USB_DEBUG
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
// STUCK USB OUTPUT / USCITA USB BLOCCATA
// ===================================================================================
// EN: IN buffer (64 bytes towards the PC) busy while the PC keeps asking for data: a ready packet would leave within
//     1 ms, so this one is stuck in the peripheral and the answers to WHOAMI stay in the dongle. Serial of Arduino
//     sent the flush (WR_DONE) again at every write while the PC was not reading; the ESP-IDF driver sends it only
//     once per packet, so a lost "packet ready" is never sent again. Here the flush again first (nothing is lost);
//     if that is not enough, the IN FIFO is reset (at most the old packet is lost) and an empty packet restarts the
//     sending. With the port closed the PC asks for nothing and nothing is done; in normal use the buffer is never
//     busy for so long.
// IT: Buffer IN (64 byte verso il PC) occupato mentre il PC continua a chiedere dati: un pacchetto pronto uscirebbe
//     entro 1 ms, quindi questo e' bloccato nella periferica e le risposte a WHOAMI restano nel dongle. Il Serial di
//     Arduino rimandava il flush (WR_DONE) a ogni scrittura mentre il PC non leggeva; il driver ESP-IDF lo manda una
//     volta sola per pacchetto, quindi un "pacchetto pronto" perso non viene piu' rimandato. Qui prima di nuovo il
//     flush (non si perde niente); se non basta, si azzera la FIFO IN (si perde al massimo il vecchio pacchetto) e un
//     pacchetto vuoto fa ripartire l'invio. Con la porta chiusa il PC non chiede niente e non si fa niente; nell'uso
//     normale il buffer non resta mai occupato cosi' a lungo.
// EN: Returns true after a reset of the IN FIFO: the lost packet may have held the state, which is then announced again.
// IT: Restituisce true dopo un azzeramento della FIFO IN: il pacchetto perso poteva contenere lo stato, che si riannuncia.
static bool usbUnstick() {
  static bool busy = false;
  static bool asked = false;
  static uint8_t step = 0;
  static uint32_t stepMs = 0;

  // EN: IN token: the PC asks for data, only with the port open. Bit not used by the ESP-IDF driver.
  // IT: token IN: il PC chiede dati, solo con la porta aperta. Bit non usato dal driver ESP-IDF.
  const uint32_t raw = usb_serial_jtag_ll_get_intraw_mask();
  const bool token = (raw & USB_SERIAL_JTAG_INTR_TOKEN_REC_IN_EP1) != 0;
  if (token) usb_serial_jtag_ll_clr_intsts_mask(USB_SERIAL_JTAG_INTR_TOKEN_REC_IN_EP1);
  #ifdef SCOPE_DONGLE_USB_DEBUG
  if (raw & USB_SERIAL_JTAG_INTR_BUS_RESET) {
    usb_serial_jtag_ll_clr_intsts_mask(USB_SERIAL_JTAG_INTR_BUS_RESET);
    if (d_busReset < 9) d_busReset = d_busReset + 1;
  }
  #endif // SCOPE_DONGLE_USB_DEBUG

  if (usb_serial_jtag_ll_txfifo_writable()) {
    busy = false;
    return false;
  }
  const uint32_t now = millis();
  if (!busy) {   // just busy: the tokens from before do not count / appena occupato: i token di prima non contano
    busy = true;
    asked = false;
    step = 0;
    stepMs = now;
    return false;
  }
  if (token) asked = true;
  if (!asked) return false;   // port closed: the PC asks for nothing / porta chiusa: il PC non chiede niente

  if (step == 0 && now - stepMs >= USB_STUCK_FLUSH_MS) {
    portENTER_CRITICAL(&g_usbKickMux);
    usb_serial_jtag_ll_txfifo_flush();
    portEXIT_CRITICAL(&g_usbKickMux);
    step = 1;
    stepMs = now;
    #ifdef SCOPE_DONGLE_USB_DEBUG
    if (d_flush < 9) d_flush = d_flush + 1;
    #endif // SCOPE_DONGLE_USB_DEBUG
  } else if (step == 1 && now - stepMs >= USB_STUCK_RESET_MS) {
    // EN: both clock domains of the IN FIFO (CPU and USB), then the empty packet: when the PC takes it, the driver
    //     sends the rest. With interrupts off, the ISR cannot fall in between. Buffer checked again in here: freed just
    //     before, it holds nothing to lose.
    // IT: tutti e due i domini di clock della FIFO IN (CPU e USB), poi il pacchetto vuoto: quando il PC lo prende, il
    //     driver manda il resto. Con gli interrupt spenti l'ISR non puo' cadere in mezzo. Buffer ricontrollato qui
    //     dentro: liberato appena prima, non c'e' niente da perdere.
    portENTER_CRITICAL(&g_usbKickMux);
    const bool reset = !usb_serial_jtag_ll_txfifo_writable();
    if (reset) {
      USB_SERIAL_JTAG.ser_afifo_config.serial_in_afifo_reset_wr = 1;
      USB_SERIAL_JTAG.ser_afifo_config.serial_in_afifo_reset_rd = 1;
      esp_rom_delay_us(1);
      USB_SERIAL_JTAG.ser_afifo_config.serial_in_afifo_reset_wr = 0;
      USB_SERIAL_JTAG.ser_afifo_config.serial_in_afifo_reset_rd = 0;
      usb_serial_jtag_ll_txfifo_flush();
    }
    portEXIT_CRITICAL(&g_usbKickMux);
    step = 2;
    stepMs = now;
    #ifdef SCOPE_DONGLE_USB_DEBUG
    if (reset && d_reset < 9) d_reset = d_reset + 1;
    #endif // SCOPE_DONGLE_USB_DEBUG
    return reset;
  } else if (step == 2 && now - stepMs >= USB_STUCK_RETRY_MS) {
    step = 0;   // still stuck: from the flush again, when the PC asks / ancora bloccato: di nuovo dal flush, quando il PC chiede
    stepMs = now;
    asked = false;
  }
  return false;
}

// ===================================================================================
// USB TASK / TASK DELLA USB
// ===================================================================================
// EN: The bytes come from the ring buffer of the ESP-IDF USB Serial/JTAG driver: its interrupt copies every USB packet
//     (64 bytes) in one go, the task takes up to 1 KB per read (Serial of Arduino moved them one byte at a time
//     through a queue, about 155 KB/s at most).
// IT: I byte arrivano dal buffer circolare del driver USB Serial/JTAG di ESP-IDF: il suo interrupt copia ogni pacchetto
//     USB (64 byte) in un colpo, il task prende fino a 1 KB per lettura (il Serial di Arduino li spostava un byte alla
//     volta in una coda, circa 155 KB/s al massimo).
static void scopeUsbTask(void *param) {
  static uint8_t content[USB_MAX_CONTENT];
  static uint8_t chunk[USB_READ_CHUNK];
  enum { WAIT_SYNC_0, WAIT_SYNC_1, WAIT_LEN_0, WAIT_LEN_1, WAIT_DATA } st = WAIT_SYNC_0;
  uint16_t len = 0;
  uint16_t pos = 0;
  uint32_t lastRxMs = millis();

  char lastState[64] = "";
  uint8_t lastRescan = WirelessScope.rescan_count;
  uint32_t lastStatMs = millis();

  for (;;) {
    // 1. EN: Reception: every packet is passed on as soon as it is complete, without waiting for the end of the frame.
    //    IT: Ricezione: ogni pacchetto si passa avanti appena e' completo, senza aspettare la fine del frame.
    //    EN: It waits for the data at most USB_READ_WAIT_MS, without polling.
    //    IT: Aspetta i dati al massimo USB_READ_WAIT_MS, senza interrogare a vuoto.
    int n = usb_serial_jtag_read_bytes(chunk, sizeof(chunk), pdMS_TO_TICKS(USB_READ_WAIT_MS));
    // EN: A frame cut by the PC (e.g. its program closed while writing) must not swallow the next ones.
    // IT: Una trama troncata dal PC (per esempio il suo programma chiuso mentre scriveva) non deve mangiare le successive.
    if (n <= 0 && st != WAIT_SYNC_0 && millis() - lastRxMs >= USB_RESYNC_MS) st = WAIT_SYNC_0;
    if (n > 0) {
      lastRxMs = millis();
      g_usbRxBytes = g_usbRxBytes + n;
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
      // EN: Marked as sent only when the whole frame went out: a skipped one (e.g. Pn at the pairing) is tried again.
      // IT: Segnato come inviato solo se la trama e' uscita intera: una saltata (per esempio Pn all'abbinamento) si ritenta.
      if ((strcmp(state, lastState) != 0 || rescan != lastRescan) && usbSendControl(state)) {
        strcpy(lastState, state);
        lastRescan = rescan;
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

    // 4. EN: Output towards the PC stuck: see usbUnstick(). After a reset of the IN FIFO the state is sent again: after
    //       the first answer the daemon asks WHOAMI only during BOOT, and the same state again changes nothing for it.
    //    IT: uscita verso il PC bloccata: vedi usbUnstick(). Dopo un azzeramento della FIFO IN lo stato si rimanda: dopo
    //       la prima risposta il demone chiede WHOAMI solo durante BOOT, e lo stesso stato di nuovo per lui non cambia niente.
    if (usbUnstick()) lastState[0] = '\0';

    // EN: Less than a whole chunk: one tick for the next packets to gather, so the task wakes about once per
    //     millisecond and not at every USB packet of 64 bytes (1 ms is about 1 KB, the ring holds 8).
    // IT: Meno di un pezzo intero: un tick perche' i pacchetti successivi si accumulino, cosi' il task si sveglia
    //     circa una volta al millisecondo e non a ogni pacchetto USB da 64 byte (1 ms e' circa 1 KB, l'anello ne tiene 8).
    if (n > 0 && n < (int)sizeof(chunk)) vTaskDelay(1);
  }
}

#ifdef SCOPE_DONGLE_USB_DEBUG
void scopeUsbDebugText(char *out, size_t size) {
  snprintf(out, size, "W%u F%u K%u Q%u B%u E%u S%04X", (unsigned)(d_whoami < 99 ? d_whoami : 99),
           (unsigned)(d_answerLost < 99 ? d_answerLost : 99), (unsigned)d_flush, (unsigned)d_reset, (unsigned)d_busReset,
           (unsigned)usb_serial_jtag_ll_txfifo_writable(), (unsigned)(USB_SERIAL_JTAG.in_ep1_st.val & 0xFFFF));
}
#endif // SCOPE_DONGLE_USB_DEBUG

// EN: Driver of ESP-IDF in place of Serial (Serial must not be started: both would use the same peripheral).
// IT: Driver di ESP-IDF al posto di Serial (Serial non va avviato: userebbero tutti e due la stessa periferica).
void scopeUsbBegin() {
  if (!usb_serial_jtag_is_driver_installed()) {
    usb_serial_jtag_driver_config_t cfg = USB_SERIAL_JTAG_DRIVER_CONFIG_DEFAULT();
    cfg.rx_buffer_size = USB_RX_BUFFER;
    cfg.tx_buffer_size = USB_TX_BUFFER;
    if (usb_serial_jtag_driver_install(&cfg) != ESP_OK) return;   // no USB: the radio works alone / niente USB: la radio lavora da sola
  }
  g_usbTxMutex = xSemaphoreCreateMutex();
  xTaskCreate(scopeUsbTask, "scopeUsb", 4096, NULL, 5, NULL);   // below the video task (6) / sotto il task del video (6)
}

#endif // OPENFIRE_WIRELESS_ENABLE && ARDUINO_ARCH_ESP32 && SCOPE_DONGLE
