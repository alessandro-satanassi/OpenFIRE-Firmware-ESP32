/*!
 * @file OpenFIRE_Wireless_Scope.cpp
 * @brief Library for OpenFIRE wireless scope with Wifi6 5Ghz UDP no ACK
 * @n CPP Library for OpenFIRE wireless scope with Wifi6 5Ghz UDP no ACK
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)

#include "OpenFIRE_Wireless_Scope.h"
#include "freertos/semphr.h"

WirelessScope_ WirelessScope;

static const uint8_t g_scopeChannels[OPENFIRE_SCOPE_N_CHANNELS] = OPENFIRE_SCOPE_CHANNELS;

// ===================================================================================
// COMMON RADIO SETTINGS / IMPOSTAZIONI RADIO COMUNI
// ===================================================================================
// EN: Called after WiFi.mode(), on the interface in use (STA on the scope, AP on the dongle).
// IT: Chiamata dopo WiFi.mode(), sull'interfaccia in uso (STA sul mirino, AP sul dongle).
void WirelessScope_::applyRadioSettings(wifi_interface_t ifx) {
  // EN: Country: with "IT" the channels 36-48 are allowed. Without ieee80211d the setting never changes.
  // IT: Paese: con "IT" i canali 36-48 sono permessi. Senza ieee80211d l'impostazione non cambia mai.
  esp_wifi_set_country_code(OPENFIRE_WIFI_COUNTRY, false);

  esp_wifi_set_band_mode(WIFI_BAND_MODE_5G_ONLY);   // 5 GHz only / solo 5 GHz

  // EN: Wi-Fi 6 (11ax) and the older 5 GHz protocols. / IT: Wi-Fi 6 (11ax) e i protocolli a 5 GHz precedenti.
  wifi_protocols_t protocols = {};
  protocols.ghz_5g = WIFI_PROTOCOL_11A | WIFI_PROTOCOL_11N | WIFI_PROTOCOL_11AC | WIFI_PROTOCOL_11AX;
  esp_wifi_set_protocols(ifx, &protocols);

  wifi_bandwidths_t bandwidths = {};
  bandwidths.ghz_5g = WIFI_BW20;                     // 20 MHz
  esp_wifi_set_bandwidths(ifx, &bandwidths);

  esp_wifi_set_max_tx_power(OPENFIRE_WIFI_POWER);

  // EN: No power save: otherwise the access point holds the packets until the scope wakes up (100 ms and more).
  // IT: Niente risparmio energetico: altrimenti l'access point trattiene i pacchetti finche' il mirino non si sveglia (100 ms e oltre).
  esp_wifi_set_ps(WIFI_PS_NONE);
}

#ifdef SCOPE_DONGLE
static void scopeVideoInit();   // video queue and task, below / coda e task del video, piu' sotto
#endif // SCOPE_DONGLE

#ifdef SCOPE
// EN: The scope leaves the dongle's network after the pairing (link lost): back to the search, the main
//     repeats the same pairing. During the search and the pairing the disconnections are handled by connection_scope().
// IT: Il mirino esce dalla rete del dongle dopo l'abbinamento (collegamento perso): torna alla ricerca, il main
//     ripete lo stesso abbinamento. Durante ricerca e abbinamento gli scollegamenti li gestisce connection_scope().
static void scopeWiFiEvent(arduino_event_id_t event, arduino_event_info_t info) {
  if (event == ARDUINO_EVENT_WIFI_STA_DISCONNECTED && WirelessScope.link_state >= SCOPE_LINK_STATE::SCOPE_LINK_PAIRED) {
    WirelessScope.link_state = SCOPE_LINK_STATE::SCOPE_LINK_SEARCH;
  }
}
#endif // SCOPE

// ===================================================================================
// BEGIN: COMMON INITIALIZATION / INIZIALIZZAZIONE COMUNE
// ===================================================================================
void WirelessScope_::begin() {
  WiFi.persistent(false);   // Wi-Fi settings only in RAM / impostazioni Wi-Fi solo in RAM
  WiFi.setSleep(false);     // also Arduino's saved preference / anche la preferenza memorizzata da Arduino
  link_state = SCOPE_LINK_STATE::SCOPE_LINK_WAIT;

  #ifdef SCOPE_DONGLE
    #ifdef OPENFIRE_AUTO_CHANNEL_WIFI
      channel = findBestChannel();
    #else
      channel = OPENFIRE_WIFI_CHANNEL;
    #endif // OPENFIRE_AUTO_CHANNEL_WIFI
    startAccessPoint(channel);
    scopeVideoInit();
  #endif // SCOPE_DONGLE

  #ifdef SCOPE
    WiFi.mode(WIFI_STA);
    applyRadioSettings(WIFI_IF_STA);
    // EN: No automatic reconnection by Arduino: a reconnection without Confirm would be dropped by the dongle after 2 s.
    //     Every pairing, the first as the next ones, goes through connection_scope().
    // IT: Niente ricollegamento automatico di Arduino: un ricollegamento senza Conferma verrebbe scollegato dal dongle dopo 2 s.
    //     Ogni abbinamento, il primo come i successivi, passa da connection_scope().
    WiFi.setAutoReconnect(false);
    WiFi.disconnect();
    esp_wifi_get_mac(WIFI_IF_STA, mac_local);
    WiFi.onEvent(scopeWiFiEvent, ARDUINO_EVENT_WIFI_STA_DISCONNECTED);
  #endif // SCOPE
}

// ===================================================================================
// SCOPE DONGLE: CHANNEL SEARCH, ACCESS POINT, PAIRING / DONGLE: RICERCA CANALE, ACCESS POINT, ABBINAMENTO
// ===================================================================================
#ifdef SCOPE_DONGLE

// EN: Vendor element in the dongle's beacons: it stays also when the network is hidden, so a dongle that
//     is searching its channel recognizes the channels already taken by other OpenFIRE dongles.
// IT: Elemento vendor nei beacon del dongle: resta anche quando la rete e' nascosta, cosi' un dongle che
//     cerca il suo canale riconosce i canali gia' presi da altri dongle OpenFIRE.
#define SCOPE_VND_OUI_0   0x4F   // 'O'
#define SCOPE_VND_OUI_1   0x46   // 'F'
#define SCOPE_VND_OUI_2   0x53   // 'S'
#define SCOPE_VND_TYPE    0x01
static const uint8_t g_scopeVendorIe[] = { 0xDD, 4 + 5, SCOPE_VND_OUI_0, SCOPE_VND_OUI_1, SCOPE_VND_OUI_2, SCOPE_VND_TYPE,
                                           'S', 'C', 'O', 'P', 'E' };

// ----------------- channel search / ricerca del canale -----------------
#define BEST_CH_RANDOM_WAIT_MS   500   // two dongles plugged in together do not listen at the same time / due dongle inseriti insieme non ascoltano nello stesso momento
#define BEST_CH_SETTLE_MS        10    // after a channel change / dopo un cambio canale
#define BEST_CH_LISTEN_MS        200   // listening per channel / ascolto per canale
#define BEST_CH_FRAME_US         40    // fixed cost of a frame on air (preamble, gap) / costo fisso di un frame in aria (preambolo, pausa)
#define BEST_CH_RATE_MBPS        24    // rate assumed for the frames heard / velocita' ipotizzata per i frame sentiti
// EN: Weights of the cost, to be tuned at the bench. / IT: pesi del costo, da regolare al banco.
#define BEST_CH_W_BUSY           1.0f  // per % of busy air / per % di etere occupato
#define BEST_CH_W_NOISE          1.0f  // per unit of noise floor above the quietest channel / per unita' di rumore sopra il canale piu' silenzioso
#define BEST_CH_W_RSSI           0.2f  // per dBm of the strongest network above -95 dBm / per dBm della rete piu' forte sopra -95 dBm

typedef struct {
  uint32_t frames;
  uint32_t airUs;      // estimated busy air time / tempo di etere occupato stimato
  int32_t  nfSum;      // noise floor sum / somma del rumore di fondo
  uint32_t nfCount;
  int8_t   rssiMax;    // strongest beacon / beacon piu' forte
  bool     openfire;   // beacon of another OpenFIRE dongle heard / sentito il beacon di un altro dongle OpenFIRE
} BestChMeasure;

static volatile BestChMeasure g_meas;
static volatile bool g_listening = false;

// EN: True if the beacon carries the OpenFIRE scope vendor element. Beacon: 24 bytes of header,
//     12 of fixed fields, then the elements; the last 4 bytes are the FCS.
// IT: Vero se il beacon contiene l'elemento vendor del mirino OpenFIRE. Beacon: 24 byte di intestazione,
//     12 di campi fissi, poi gli elementi; gli ultimi 4 byte sono l'FCS.
static bool beaconHasScopeIe(const uint8_t *p, int len) {
  int i = 24 + 12;
  while (i + 2 <= len) {
    const uint8_t id = p[i];
    const uint8_t l  = p[i + 1];
    if (i + 2 + l > len) break;
    if (id == 0xDD && l >= 4 && p[i + 2] == SCOPE_VND_OUI_0 && p[i + 3] == SCOPE_VND_OUI_1 &&
        p[i + 4] == SCOPE_VND_OUI_2 && p[i + 5] == SCOPE_VND_TYPE) return true;
    i += 2 + l;
  }
  return false;
}

static void bestChPromiscuousCb(void *buf, wifi_promiscuous_pkt_type_t type) {
  if (!g_listening) return;
  const wifi_promiscuous_pkt_t *pkt = (const wifi_promiscuous_pkt_t *)buf;
  const uint32_t len = pkt->rx_ctrl.sig_len;

  g_meas.frames = g_meas.frames + 1;
  g_meas.airUs  = g_meas.airUs + BEST_CH_FRAME_US + (len * 8) / BEST_CH_RATE_MBPS;
  g_meas.nfSum  = g_meas.nfSum + pkt->rx_ctrl.noise_floor;
  g_meas.nfCount = g_meas.nfCount + 1;

  if (type == WIFI_PKT_MGMT && len > 4 && pkt->payload[0] == 0x80) {   // beacon
    if (pkt->rx_ctrl.rssi > g_meas.rssiMax) g_meas.rssiMax = pkt->rx_ctrl.rssi;
    if (!g_meas.openfire && beaconHasScopeIe(pkt->payload, (int)len - 4)) g_meas.openfire = true;
  }
}

uint8_t WirelessScope_::findBestChannel(uint8_t skipMask) {
  link_state = SCOPE_LINK_STATE::SCOPE_LINK_WAIT;
  vTaskDelay(pdMS_TO_TICKS(esp_random() % (BEST_CH_RANDOM_WAIT_MS + 1)));

  WiFi.mode(WIFI_STA);
  applyRadioSettings(WIFI_IF_STA);
  WiFi.disconnect();

  BestChMeasure meas[OPENFIRE_SCOPE_N_CHANNELS] = {};

  wifi_promiscuous_filter_t filter = {};
  filter.filter_mask = WIFI_PROMIS_FILTER_MASK_MGMT | WIFI_PROMIS_FILTER_MASK_CTRL | WIFI_PROMIS_FILTER_MASK_DATA |
                       WIFI_PROMIS_FILTER_MASK_MISC | WIFI_PROMIS_FILTER_MASK_FCSFAIL;
  wifi_promiscuous_filter_t ctrlFilter = {};
  ctrlFilter.filter_mask = WIFI_PROMIS_CTRL_FILTER_MASK_ALL;
  esp_wifi_set_promiscuous_filter(&filter);
  esp_wifi_set_promiscuous_ctrl_filter(&ctrlFilter);
  esp_wifi_set_promiscuous_rx_cb(&bestChPromiscuousCb);
  esp_wifi_set_promiscuous(true);

  for (uint8_t i = 0; i < OPENFIRE_SCOPE_N_CHANNELS; i++) {
    if (skipMask & (1 << i)) continue;
    esp_wifi_set_channel(g_scopeChannels[i], WIFI_SECOND_CHAN_NONE);
    vTaskDelay(pdMS_TO_TICKS(BEST_CH_SETTLE_MS));
    g_meas.frames = 0; g_meas.airUs = 0; g_meas.nfSum = 0; g_meas.nfCount = 0;
    g_meas.rssiMax = -128; g_meas.openfire = false;
    g_listening = true;
    vTaskDelay(pdMS_TO_TICKS(BEST_CH_LISTEN_MS));
    g_listening = false;
    memcpy(&meas[i], (const void *)&g_meas, sizeof(BestChMeasure));
  }

  esp_wifi_set_promiscuous(false);
  esp_wifi_set_promiscuous_rx_cb(NULL);
  vTaskDelay(pdMS_TO_TICKS(20));   // a callback still running ends / una callback ancora in corso termina

  // EN: Quietest noise floor among the channels heard: the noise of a channel is counted above it.
  // IT: Rumore di fondo piu' basso tra i canali ascoltati: il rumore di un canale si conta sopra questo.
  float nfMin = 0.0f;
  bool nfMinSet = false;
  for (uint8_t i = 0; i < OPENFIRE_SCOPE_N_CHANNELS; i++) {
    if ((skipMask & (1 << i)) || meas[i].nfCount == 0) continue;
    const float nf = (float)meas[i].nfSum / meas[i].nfCount;
    if (!nfMinSet || nf < nfMin) { nfMin = nf; nfMinSet = true; }
  }

  // EN: Two rounds: first without the channels of other OpenFIRE dongles, then (if all are taken) with them.
  // IT: Due giri: prima senza i canali di altri dongle OpenFIRE, poi (se sono tutti presi) anche con quelli.
  for (uint8_t round = 0; round < 2; round++) {
    int8_t best = -1;
    float bestCost = 0.0f;
    for (uint8_t i = 0; i < OPENFIRE_SCOPE_N_CHANNELS; i++) {
      if (skipMask & (1 << i)) continue;
      if (round == 0 && meas[i].openfire) continue;
      float cost = BEST_CH_W_BUSY * (100.0f * meas[i].airUs / (BEST_CH_LISTEN_MS * 1000.0f));
      if (meas[i].nfCount) cost += BEST_CH_W_NOISE * ((float)meas[i].nfSum / meas[i].nfCount - nfMin);
      if (meas[i].rssiMax > -95) cost += BEST_CH_W_RSSI * (meas[i].rssiMax + 95);
      if (best < 0 || cost < bestCost) { best = i; bestCost = cost; }
    }
    if (best >= 0) return g_scopeChannels[best];
  }
  return OPENFIRE_WIFI_CHANNEL;   // every channel skipped / tutti i canali esclusi
}

// ----------------- access point / access point -----------------

// EN: Station connected to the access point, not yet confirmed. / IT: stazione collegata all'access point, non ancora confermata.
static volatile bool     g_staConnected = false;
// EN: Station present on the access point, confirmed or not. / IT: stazione presente sull'access point, confermata o no.
static volatile bool     g_staPresent = false;
static volatile uint32_t g_staConnMs = 0;
static volatile uint16_t g_staAid = 0;
static bool g_eventRegistered = false;

// EN: Hidden network = busy dongle: it does not answer the scopes' searches.
//     TO BE CHECKED AT THE BENCH: changing the access point configuration must not disconnect the scope.
// IT: Rete nascosta = dongle occupato: non risponde alle ricerche dei mirini.
//     DA VERIFICARE AL BANCO: cambiare la configurazione dell'access point non deve scollegare il mirino.
static void setApHidden(bool hidden) {
  wifi_config_t cfg = {};
  if (esp_wifi_get_config(WIFI_IF_AP, &cfg) != ESP_OK) return;
  if (cfg.ap.ssid_hidden == (hidden ? 1 : 0)) return;
  cfg.ap.ssid_hidden = hidden ? 1 : 0;
  esp_wifi_set_config(WIFI_IF_AP, &cfg);
}

static void scopeDongleWiFiEvent(arduino_event_id_t event, arduino_event_info_t info) {
  switch (event) {
    case ARDUINO_EVENT_WIFI_AP_STACONNECTED:
      memcpy(WirelessScope.mac_peer, info.wifi_ap_staconnected.mac, 6);
      g_staAid = info.wifi_ap_staconnected.aid;
      g_staConnMs = millis();
      g_staConnected = true;
      g_staPresent = true;
      break;
    case ARDUINO_EVENT_WIFI_AP_STADISCONNECTED:
      g_staConnected = false;
      g_staPresent = false;
      // EN: The scope left: the dongle is free and visible again. / IT: il mirino se n'e' andato: il dongle torna libero e visibile.
      if (WirelessScope.link_state >= SCOPE_LINK_STATE::SCOPE_LINK_PAIRED) {
        WirelessScope.link_state = SCOPE_LINK_STATE::SCOPE_LINK_SEARCH;
        setApHidden(false);
      }
      break;
    default:
      break;
  }
}

bool WirelessScope_::startAccessPoint(uint8_t ch) {
  WiFi.mode(WIFI_AP);
  applyRadioSettings(WIFI_IF_AP);
  WiFi.softAPConfig(IPAddress(OPENFIRE_SCOPE_DONGLE_IP), IPAddress(OPENFIRE_SCOPE_DONGLE_IP), IPAddress(OPENFIRE_SCOPE_NETMASK));
  if (!WiFi.softAP(OPENFIRE_SCOPE_SSID, NULL, ch, 0, 1)) return false;

  // EN: Values forced on the driver too: open network, visible, one station only, channel of the pair.
  // IT: Valori forzati anche sul driver: rete aperta, visibile, una sola stazione, canale della coppia.
  wifi_config_t cfg = {};
  esp_wifi_get_config(WIFI_IF_AP, &cfg);
  cfg.ap.channel = ch;
  cfg.ap.authmode = WIFI_AUTH_OPEN;
  cfg.ap.ssid_hidden = 0;
  cfg.ap.max_connection = 1;
  esp_wifi_set_config(WIFI_IF_AP, &cfg);

  // EN: During the video the scope does not transmit: the access point must not disconnect it for inactivity
  //     (65535 s is the maximum of the API).
  // IT: Durante il video il mirino non trasmette: l'access point non deve scollegarlo per inattivita'
  //     (65535 s e' il massimo dell'API).
  esp_wifi_set_inactive_time(WIFI_IF_AP, 65535);

  esp_wifi_set_vendor_ie(false, WIFI_VND_IE_TYPE_BEACON, WIFI_VND_IE_ID_0, NULL);   // also after a RESCAN / anche dopo un RESCAN
  esp_wifi_set_vendor_ie(true, WIFI_VND_IE_TYPE_BEACON, WIFI_VND_IE_ID_0, g_scopeVendorIe);

  esp_wifi_get_mac(WIFI_IF_AP, mac_local);
  if (!g_eventRegistered) {
    WiFi.onEvent(scopeDongleWiFiEvent);
    g_eventRegistered = true;
  }
  channel = ch;
  g_staConnected = false;   // new access point: no station yet / access point nuovo: ancora nessuna stazione
  g_staPresent = false;
  link_state = SCOPE_LINK_STATE::SCOPE_LINK_SEARCH;   // free / libero
  return true;
}

// ----------------- pairing / abbinamento -----------------
bool WirelessScope_::connection_scope_dongle() {
  link_state = SCOPE_LINK_STATE::SCOPE_LINK_SEARCH;
  setApHidden(false);
  udp_pair.begin(OPENFIRE_SCOPE_PAIR_PORT);

  while (true) {
    // EN: RESCAN from the PC: the main does the new channel search. / IT: RESCAN dal PC: la nuova ricerca del canale la fa il main.
    if (rescan_request) {
      udp_pair.stop();
      return false;
    }

    // EN: Connected without a Confirm for 2 s: the station is disconnected and the place is free again.
    // IT: Collegato senza Conferma da 2 s: la stazione viene scollegata e il posto torna libero.
    if (g_staConnected && (millis() - g_staConnMs > OPENFIRE_SCOPE_DONGLE_CONFIRM_MS)) {
      g_staConnected = false;
      esp_wifi_deauth_sta(g_staAid);
    }

    if (udp_pair.parsePacket() == sizeof(Scope_Pair_Msg)) {
      Scope_Pair_Msg msg;
      udp_pair.read((uint8_t *)&msg, sizeof(msg));
      // EN: Only from the station waiting for its Confirm: a Confirm still in the buffer when the 2 s check
      //     has just removed the station is ignored (otherwise the dongle would stay busy without a scope).
      // IT: Solo dalla stazione che aspetta la sua Conferma: una Conferma ancora nel buffer quando il controllo
      //     dei 2 s ha appena tolto la stazione viene ignorata (altrimenti il dongle resterebbe occupato senza mirino).
      if (g_staConnected && msg.marker == OPENFIRE_SCOPE_PAIR_MARKER && msg.version == OPENFIRE_SCOPE_PAIR_VERSION &&
          msg.type == SCOPE_PAIR_MSG::SCOPE_PAIR_CONFIRM && msg.player >= 1 && msg.player <= 4) {
        // EN: The first valid Confirm wins: from now on the dongle is busy. / IT: vince la prima Conferma valida: da qui il dongle e' occupato.
        player = msg.player;
        ip_peer = udp_pair.remoteIP();
        const uint16_t port = udp_pair.remotePort();
        link_state = SCOPE_LINK_STATE::SCOPE_LINK_PAIRED;
        g_staConnected = false;   // no more 2 s check / niente piu' controllo dei 2 s

        msg.type = SCOPE_PAIR_MSG::SCOPE_PAIR_OK;   // same nonce and player / stesso numero casuale e giocatore
        for (uint8_t i = 0; i < OPENFIRE_SCOPE_OK_REPEAT; i++) {
          if (i > 0) vTaskDelay(pdMS_TO_TICKS(OPENFIRE_SCOPE_OK_INTERVAL_MS));
          udp_pair.beginPacket(ip_peer, port);
          udp_pair.write((const uint8_t *)&msg, sizeof(msg));
          udp_pair.endPacket();
        }
        setApHidden(true);
        udp_pair.stop();
        // EN: The scope left while the OKs were going out: the event may have found the dongle not yet busy.
        //     Back to free: the main calls connection_scope_dongle() again, which makes the network visible.
        // IT: Il mirino e' uscito mentre partivano gli OK: l'evento puo' aver trovato il dongle non ancora occupato.
        //     Si torna liberi: il main richiama connection_scope_dongle(), che rende di nuovo visibile la rete.
        if (!g_staPresent) link_state = SCOPE_LINK_STATE::SCOPE_LINK_SEARCH;
        return (link_state == SCOPE_LINK_STATE::SCOPE_LINK_PAIRED);
      }
    }
    vTaskDelay(pdMS_TO_TICKS(1));
  }
}

// ----------------- RESCAN from the PC / RESCAN dal PC -----------------
// EN: Never with a scope paired: an active pair is never interrupted. / IT: mai con un mirino abbinato: una coppia attiva non si interrompe mai.
bool WirelessScope_::requestRescan(uint8_t skipMask) {
  if (link_state != SCOPE_LINK_STATE::SCOPE_LINK_SEARCH) return false;
  rescan_skip_mask = skipMask;
  rescan_request = true;
  return true;
}

bool WirelessScope_::rescanChannel() {
  rescan_request = false;
  if (link_state != SCOPE_LINK_STATE::SCOPE_LINK_SEARCH) return false;   // a scope paired in the meantime / un mirino si e' abbinato nel frattempo
  channel = findBestChannel(rescan_skip_mask);   // WAIT meanwhile: WHOAMI answers BOOT / intanto WAIT: WHOAMI risponde BOOT
  const bool ok = startAccessPoint(channel);
  rescan_count = rescan_count + 1;
  return ok;
}

// ----------------- video: queue and transmission / video: coda e trasmissione -----------------
// EN: The USB task puts the packets of the PC in the queue, the video task passes them to the Wi-Fi one at a time.
//     When the first slice of a newer frame arrives, the slices of the older frames still waiting are discarded:
//     the most recent frame always wins. The session and bit 3 are written when the packet enters the queue,
//     so after a START a slice of the old session never gets the new number.
// IT: Il task della USB mette in coda i pacchetti del PC, il task del video li passa al Wi-Fi uno alla volta.
//     Quando arriva la prima fascia di un frame piu' nuovo, le fasce dei frame vecchi ancora in attesa si scartano:
//     vince sempre il frame piu' recente. Sessione e bit 3 si scrivono quando il pacchetto entra in coda,
//     cosi' dopo uno START una fascia della sessione vecchia non riceve mai il numero nuovo.
typedef struct {
  uint16_t len;
  uint8_t  data[OPENFIRE_SCOPE_VIDEO_MAX_PACKET];
} ScopeVideoSlot;

static ScopeVideoSlot     g_vq[OPENFIRE_SCOPE_DONGLE_QUEUE];
static uint8_t            g_vqHead = 0;
static uint8_t            g_vqCount = 0;
static uint16_t           g_vqLastFrame = 0;          // frame of the last slice queued / frame dell'ultima fascia in coda
static SemaphoreHandle_t  g_vqMutex = NULL;
static TaskHandle_t       g_videoTask = NULL;

static volatile uint32_t  g_lastVideoMs = 0;          // last video packet from the PC / ultimo pacchetto video dal PC
static volatile uint32_t  g_lastPauseMs = 0;          // last PAUSE from the PC / ultima PAUSE dal PC
static volatile bool      g_pauseSeen = false;
static volatile bool      g_sessionFlagPending = false;   // START received, first video not yet / START ricevuto, primo video non ancora
static volatile uint32_t  g_sessionVideoMs = 0;       // first video of the session / primo video della sessione
static volatile uint16_t  g_lastFrameSent = 0;

static inline uint16_t videoFrame(const uint8_t *p) {
  return (uint16_t)p[SCOPE_VIDEO_B_FRAME] | ((uint16_t)p[SCOPE_VIDEO_B_FRAME + 1] << 8);
}
// EN: a more recent than b: 16 bit difference read with sign. / IT: a piu' recente di b: differenza a 16 bit letta con segno.
static inline bool videoFrameNewer(uint16_t a, uint16_t b) { return (int16_t)(a - b) > 0; }

static inline uint8_t videoFlags(uint8_t flags) {
  // EN: Version and the flags of the PC, bit 3 only in the first second of video of the session.
  // IT: Versione e flag del PC, bit 3 solo nel primo secondo di video della sessione.
  flags &= (uint8_t)~OPENFIRE_SCOPE_VIDEO_FLAG_SESSION;
  if (!g_sessionFlagPending && (millis() - g_sessionVideoMs < OPENFIRE_SCOPE_SESSION_FLAG_MS)) flags |= OPENFIRE_SCOPE_VIDEO_FLAG_SESSION;
  return flags;
}

bool WirelessScope_::videoPush(const uint8_t *pkt, uint16_t len) {
  if (g_vqMutex == NULL || len < OPENFIRE_SCOPE_VIDEO_HEADER || len > OPENFIRE_SCOPE_VIDEO_MAX_PACKET) return false;
  g_lastVideoMs = millis();
  if (link_state < SCOPE_LINK_STATE::SCOPE_LINK_PAIRED) return false;   // no scope: discarded / nessun mirino: scartato

  const uint16_t frame = videoFrame(pkt);
  xSemaphoreTake(g_vqMutex, portMAX_DELAY);
  if (g_vqCount > 0 && videoFrameNewer(frame, g_vqLastFrame)) {   // older frames still waiting / frame vecchi ancora in attesa
    video_drop = video_drop + g_vqCount;
    g_vqCount = 0;
  }
  if (g_vqCount == OPENFIRE_SCOPE_DONGLE_QUEUE) {   // full: the oldest slice goes out / piena: esce la fascia piu' vecchia
    g_vqHead = (g_vqHead + 1) % OPENFIRE_SCOPE_DONGLE_QUEUE;
    g_vqCount--;
    video_drop = video_drop + 1;
  }
  if (g_sessionFlagPending) {   // first video after START / primo video dopo START
    g_sessionFlagPending = false;
    g_sessionVideoMs = millis();
  }
  ScopeVideoSlot &s = g_vq[(g_vqHead + g_vqCount) % OPENFIRE_SCOPE_DONGLE_QUEUE];
  memcpy(s.data, pkt, len);
  s.len = len;
  s.data[SCOPE_VIDEO_B_SESSION] = video_session;
  s.data[SCOPE_VIDEO_B_FLAGS] = videoFlags(s.data[SCOPE_VIDEO_B_FLAGS]);
  g_vqCount++;
  g_vqLastFrame = frame;
  xSemaphoreGive(g_vqMutex);

  if (g_videoTask != NULL) xTaskNotifyGive(g_videoTask);
  return true;
}

void WirelessScope_::videoStart() {
  if (g_vqMutex == NULL) return;
  xSemaphoreTake(g_vqMutex, portMAX_DELAY);
  video_drop = video_drop + g_vqCount;   // the queue is emptied / la coda si svuota
  g_vqCount = 0;
  video_session = (video_session >= 255) ? 1 : (uint8_t)(video_session + 1);
  g_sessionFlagPending = true;
  xSemaphoreGive(g_vqMutex);
}

void WirelessScope_::videoPause() {
  g_lastPauseMs = millis();
  g_pauseSeen = true;
}

static bool videoPop(uint8_t *out, uint16_t *len) {
  bool ok = false;
  xSemaphoreTake(g_vqMutex, portMAX_DELAY);
  if (g_vqCount > 0) {
    const ScopeVideoSlot &s = g_vq[g_vqHead];
    memcpy(out, s.data, s.len);
    *len = s.len;
    g_vqHead = (g_vqHead + 1) % OPENFIRE_SCOPE_DONGLE_QUEUE;
    g_vqCount--;
    ok = true;
  }
  xSemaphoreGive(g_vqMutex);
  return ok;
}

static bool videoNewerWaiting(uint16_t frame) {
  xSemaphoreTake(g_vqMutex, portMAX_DELAY);
  const bool newer = (g_vqCount > 0 && videoFrameNewer(g_vqLastFrame, frame));
  xSemaphoreGive(g_vqMutex);
  return newer;
}

static void scopeVideoTask(void *param) {
  WiFiUDP udp_video;
  udp_video.begin(OPENFIRE_SCOPE_VIDEO_PORT);
  static uint8_t pkt[OPENFIRE_SCOPE_VIDEO_MAX_PACKET];
  uint16_t len = 0;
  uint32_t lastPresenceMs = 0;

  for (;;) {
    ulTaskNotifyTake(pdTRUE, pdMS_TO_TICKS(10));   // a new packet or 10 ms / un pacchetto nuovo oppure 10 ms

    // 1. EN: Video: one packet at a time to the Wi-Fi stack. / IT: video: un pacchetto alla volta allo stack Wi-Fi.
    while (videoPop(pkt, &len)) {
      if (WirelessScope.link_state < SCOPE_LINK_STATE::SCOPE_LINK_PAIRED) continue;   // the scope left / il mirino se n'e' andato
      const uint16_t frame = videoFrame(pkt);
      g_lastFrameSent = frame;
      for (uint8_t attempt = 0; attempt < OPENFIRE_SCOPE_TX_RETRY; attempt++) {
        udp_video.beginPacket(WirelessScope.ip_peer, OPENFIRE_SCOPE_VIDEO_PORT);
        udp_video.write(pkt, len);
        if (udp_video.endPacket()) {
          WirelessScope.video_tx_bytes = WirelessScope.video_tx_bytes + len;
          break;
        }
        // EN: Wi-Fi transmit queue full: a new attempt in 1 ms, unless a newer frame is already waiting.
        // IT: Coda di trasmissione Wi-Fi piena: nuovo tentativo fra 1 ms, a meno che un frame piu' nuovo sia gia' in attesa.
        if (attempt == 0) WirelessScope.video_qfull = WirelessScope.video_qfull + 1;
        if (videoNewerWaiting(frame) || attempt == OPENFIRE_SCOPE_TX_RETRY - 1) {
          xSemaphoreTake(g_vqMutex, portMAX_DELAY);
          WirelessScope.video_drop = WirelessScope.video_drop + 1;
          xSemaphoreGive(g_vqMutex);
          break;
        }
        vTaskDelay(pdMS_TO_TICKS(1));
      }
    }

    // 2. EN: Presence: no video from the PC for 200 ms -> header only, slice 0xFF, every 200 ms.
    //        q = 1 (game paused) only with a PAUSE in the last 500 ms, otherwise q = 0 (no source).
    //    IT: Presenza: niente video dal PC per 200 ms -> solo intestazione, fascia 0xFF, ogni 200 ms.
    //        q = 1 (gioco in pausa) solo con una PAUSE negli ultimi 500 ms, altrimenti q = 0 (nessuna sorgente).
    const uint32_t now = millis();
    if (WirelessScope.link_state >= SCOPE_LINK_STATE::SCOPE_LINK_PAIRED &&
        now - g_lastVideoMs >= OPENFIRE_SCOPE_PRESENCE_IDLE_MS && now - lastPresenceMs >= OPENFIRE_SCOPE_PRESENCE_EVERY_MS) {
      const uint16_t frame = g_lastFrameSent;
      uint8_t presence[OPENFIRE_SCOPE_VIDEO_HEADER];
      presence[SCOPE_VIDEO_B_MARKER]    = OPENFIRE_SCOPE_VIDEO_MARKER;
      presence[SCOPE_VIDEO_B_FLAGS]     = videoFlags(OPENFIRE_SCOPE_VIDEO_VERSION << 4);
      presence[SCOPE_VIDEO_B_FRAME]     = (uint8_t)(frame & 0xFF);
      presence[SCOPE_VIDEO_B_FRAME + 1] = (uint8_t)(frame >> 8);
      presence[SCOPE_VIDEO_B_SLICE]     = OPENFIRE_SCOPE_VIDEO_PRESENCE;
      presence[SCOPE_VIDEO_B_Q]         = (g_pauseSeen && now - g_lastPauseMs < OPENFIRE_SCOPE_PAUSE_VALID_MS) ? 1 : 0;
      presence[SCOPE_VIDEO_B_PLAYER]    = WirelessScope.player;
      presence[SCOPE_VIDEO_B_SESSION]   = WirelessScope.video_session;
      udp_video.beginPacket(WirelessScope.ip_peer, OPENFIRE_SCOPE_VIDEO_PORT);
      udp_video.write(presence, sizeof(presence));
      if (udp_video.endPacket()) WirelessScope.video_tx_bytes = WirelessScope.video_tx_bytes + sizeof(presence);
      lastPresenceMs = now;
    }
  }
}

static void scopeVideoInit() {
  if (g_vqMutex == NULL) g_vqMutex = xSemaphoreCreateMutex();
  if (g_videoTask == NULL) xTaskCreate(scopeVideoTask, "scopeVideo", 4096, NULL, 6, &g_videoTask);   // above the USB task (5): every packet goes out at once / sopra il task della USB (5): ogni pacchetto parte subito
}

#endif // SCOPE_DONGLE

// ===================================================================================
// SCOPE: DONGLE SEARCH AND PAIRING / MIRINO: RICERCA DEL DONGLE E ABBINAMENTO
// ===================================================================================
#ifdef SCOPE

bool WirelessScope_::connection_scope(uint8_t player_number) {
  player = player_number;
  WiFi.disconnect();   // from a lost link: no connection still in progress / da un collegamento perso: nessun collegamento ancora in corso

  while (true) {
    link_state = SCOPE_LINK_STATE::SCOPE_LINK_SEARCH;

    // 1. EN: Active search on 36, 40, 44, 48. The probe request has no network name (wildcard): a hidden
    //        network, that is a busy dongle, does not answer it. The first free dongle wins, not the strongest.
    //    IT: Ricerca attiva su 36, 40, 44, 48. La probe request non ha il nome di rete (generica): una rete
    //        nascosta, cioe' un dongle occupato, non risponde. Vince il primo dongle libero, non il piu' forte.
    uint8_t bssid[6];
    uint8_t ch = 0;
    for (uint8_t i = 0; i < OPENFIRE_SCOPE_N_CHANNELS && ch == 0; i++) {
      const int16_t n = WiFi.scanNetworks(false, false, false, OPENFIRE_SCOPE_SCAN_MS_PER_CHANNEL, g_scopeChannels[i]);
      for (int16_t j = 0; j < n; j++) {
        if (WiFi.SSID(j) == OPENFIRE_SCOPE_SSID && WiFi.channel(j) == g_scopeChannels[i]) {
          memcpy(bssid, WiFi.BSSID(j), 6);
          ch = g_scopeChannels[i];
          break;
        }
      }
      WiFi.scanDelete();
    }
    if (ch == 0) continue;   // nothing found: next round / niente trovato: giro successivo

    // 2. EN: Connection to that precise dongle (BSSID), fixed address, power save off.
    //    IT: Collegamento proprio a quel dongle (BSSID), indirizzo fisso, risparmio energetico spento.
    WiFi.config(IPAddress(OPENFIRE_SCOPE_IP), IPAddress(OPENFIRE_SCOPE_DONGLE_IP), IPAddress(OPENFIRE_SCOPE_NETMASK));
    WiFi.setSleep(false);
    WiFi.begin(OPENFIRE_SCOPE_SSID, NULL, ch, bssid);
    uint32_t t0 = millis();
    while (WiFi.status() != WL_CONNECTED && (millis() - t0 < OPENFIRE_SCOPE_PAIR_TIMEOUT_MS)) vTaskDelay(pdMS_TO_TICKS(5));
    if (WiFi.status() != WL_CONNECTED) {
      WiFi.disconnect();
      continue;
    }
    // EN: Check that power save is really off once connected. / IT: verifica che il risparmio energetico sia davvero spento a collegamento avvenuto.
    wifi_ps_type_t ps;
    if (esp_wifi_get_ps(&ps) == ESP_OK && ps != WIFI_PS_NONE) esp_wifi_set_ps(WIFI_PS_NONE);

    // 3. EN: Confirm with the player, once. / IT: Conferma con il giocatore, una volta sola.
    Scope_Pair_Msg msg = { OPENFIRE_SCOPE_PAIR_MARKER, OPENFIRE_SCOPE_PAIR_VERSION, SCOPE_PAIR_MSG::SCOPE_PAIR_CONFIRM, player, esp_random() };
    udp_pair.begin(OPENFIRE_SCOPE_PAIR_PORT);
    udp_pair.beginPacket(IPAddress(OPENFIRE_SCOPE_DONGLE_IP), OPENFIRE_SCOPE_PAIR_PORT);
    udp_pair.write((const uint8_t *)&msg, sizeof(msg));
    udp_pair.endPacket();

    // 4. EN: Paired at the first OK with the same random number; the copies are ignored. No OK within 3 s: disconnect and search again.
    //    IT: Abbinato al primo OK con lo stesso numero casuale; le copie si ignorano. Niente OK entro 3 s: si scollega e ricerca.
    bool paired = false;
    t0 = millis();
    // EN: Also stops if the dongle has already removed the scope (Confirm lost: it does so after 2 s).
    // IT: Si ferma anche se il dongle ha gia' tolto il mirino (Conferma persa: lo fa dopo 2 s).
    while (!paired && WiFi.status() == WL_CONNECTED && (millis() - t0 < OPENFIRE_SCOPE_PAIR_TIMEOUT_MS)) {
      if (udp_pair.parsePacket() == sizeof(Scope_Pair_Msg)) {
        Scope_Pair_Msg rx;
        udp_pair.read((uint8_t *)&rx, sizeof(rx));
        paired = (rx.marker == OPENFIRE_SCOPE_PAIR_MARKER && rx.version == OPENFIRE_SCOPE_PAIR_VERSION &&
                  rx.type == SCOPE_PAIR_MSG::SCOPE_PAIR_OK && rx.nonce == msg.nonce);
      }
      if (!paired) vTaskDelay(pdMS_TO_TICKS(1));
    }
    udp_pair.stop();

    if (paired) {
      channel = ch;
      memcpy(mac_peer, bssid, 6);
      ip_peer = IPAddress(OPENFIRE_SCOPE_DONGLE_IP);
      link_state = SCOPE_LINK_STATE::SCOPE_LINK_PAIRED;
      // EN: The link may have dropped while the OK was arriving: the event found the scope still searching.
      // IT: Il collegamento puo' essere caduto mentre arrivava l'OK: l'evento ha trovato il mirino ancora in ricerca.
      if (WiFi.status() == WL_CONNECTED) return true;
      link_state = SCOPE_LINK_STATE::SCOPE_LINK_SEARCH;
    }
    WiFi.disconnect();   // disconnecting also frees the dongle / scollegandosi libera anche il dongle
  }
}

#endif // SCOPE

#endif // OPENFIRE_WIRELESS_ENABLE
