/*!
 * @file OpenFIRE_Wireless_Scope.h
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

#pragma once

// ===================================================================================
// SHARED LIBRARY SCOPE / SCOPE_DONGLE / LIBRERIA CONDIVISA SCOPE / SCOPE_DONGLE
// ===================================================================================
// EN: The same library is built for both devices, the role is chosen by the build flag:
//     -D SCOPE        -> the scope (station): searches a free dongle, connects, sends the player
//     -D SCOPE_DONGLE -> the scope dongle (access point): chooses the channel, waits for a scope
//     Link: Wi-Fi 6 at 5 GHz, 20 MHz, open network, UDP. No ESP-NOW.
// IT: La stessa libreria e' compilata per entrambi i dispositivi, il ruolo lo sceglie il flag di compilazione:
//     -D SCOPE        -> il mirino (stazione): cerca un dongle libero, si collega, comunica il giocatore
//     -D SCOPE_DONGLE -> il dongle del mirino (access point): sceglie il canale, aspetta un mirino
//     Collegamento: Wi-Fi 6 a 5 GHz, 20 MHz, rete aperta, UDP. Niente ESP-NOW.

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)

#include <Arduino.h>
#include <WiFi.h>
#include <WiFiUdp.h>
#include "esp_wifi.h"

// ===================================================================================
// COMMON SETTINGS / IMPOSTAZIONI COMUNI
// ===================================================================================

// EN: Network name, the same for every dongle: the scope connects to a precise dongle by its BSSID.
// IT: Nome di rete, uguale per tutti i dongle: il mirino si collega a un dongle preciso tramite il suo BSSID.
#define OPENFIRE_SCOPE_SSID          "OpenFIRE_SCOPE"

// EN: UDP ports: pairing and video use different ports. / IT: porte UDP: abbinamento e video usano porte diverse.
#define OPENFIRE_SCOPE_PAIR_PORT     4210
#define OPENFIRE_SCOPE_VIDEO_PORT    4211

// EN: Fixed addresses: one access point, one station. A fixed address on the scope avoids waiting for DHCP.
// IT: Indirizzi fissi: un access point, una stazione. L'indirizzo fisso sul mirino evita l'attesa del DHCP.
#define OPENFIRE_SCOPE_DONGLE_IP     192, 168, 4, 1
#define OPENFIRE_SCOPE_IP            192, 168, 4, 2
#define OPENFIRE_SCOPE_NETMASK       255, 255, 255, 0

// EN: 5 GHz channels without DFS (Europe). / IT: canali a 5 GHz senza DFS (Europa).
#define OPENFIRE_SCOPE_CHANNELS      { 36, 40, 44, 48 }
#define OPENFIRE_SCOPE_N_CHANNELS    4

#ifndef OPENFIRE_WIFI_CHANNEL
  #define OPENFIRE_WIFI_CHANNEL 36   // default channel of the dongle / canale predefinito del dongle
#endif // OPENFIRE_WIFI_CHANNEL

#ifndef OPENFIRE_WIFI_POWER
  #define OPENFIRE_WIFI_POWER 68     // 0.25 dBm units: 68 = 17 dBm / unita' da 0,25 dBm: 68 = 17 dBm
#endif // OPENFIRE_WIFI_POWER

#ifndef OPENFIRE_WIFI_COUNTRY
  #define OPENFIRE_WIFI_COUNTRY "IT"
#endif // OPENFIRE_WIFI_COUNTRY

// EN: Pairing times (same logic as lightgun and dongle). / IT: tempi dell'abbinamento (stessa logica di lightgun e dongle).
#define OPENFIRE_SCOPE_SCAN_MS_PER_CHANNEL   30    // scope: active scan per channel / mirino: ricerca attiva per canale
#define OPENFIRE_SCOPE_PAIR_TIMEOUT_MS       3000  // scope: no OK within 3 s -> disconnect and search again / mirino: niente OK in 3 s -> si scollega e ricerca
#define OPENFIRE_SCOPE_DONGLE_CONFIRM_MS     2000  // dongle: no Confirm within 2 s from the connection -> frees the place / dongle: niente Conferma in 2 s dal collegamento -> libera il posto
#define OPENFIRE_SCOPE_OK_REPEAT             3     // dongle: OK sent 3 times... / dongle: OK inviato 3 volte...
#define OPENFIRE_SCOPE_OK_INTERVAL_MS        69    // ...69 ms apart / ...a 69 ms l'una dall'altra

// ===================================================================================
// PAIRING MESSAGES (UDP) / MESSAGGI DI ABBINAMENTO (UDP)
// ===================================================================================
// EN: Scope -> dongle: CONFIRM, once, with a random number and the player.
//     Dongle -> scope: OK, 3 times, with the same random number.
// IT: Mirino -> dongle: CONFERMA, una volta, con un numero casuale e il giocatore.
//     Dongle -> mirino: OK, 3 volte, con lo stesso numero casuale.

#define OPENFIRE_SCOPE_PAIR_MARKER   0x50   // 'P'
#define OPENFIRE_SCOPE_PAIR_VERSION  1

enum SCOPE_PAIR_MSG {
  SCOPE_PAIR_CONFIRM = 1,   // scope -> dongle / mirino -> dongle
  SCOPE_PAIR_OK      = 2    // dongle -> scope / dongle -> mirino
};

typedef struct __attribute__ ((packed)) {
  uint8_t  marker;    // OPENFIRE_SCOPE_PAIR_MARKER
  uint8_t  version;   // OPENFIRE_SCOPE_PAIR_VERSION
  uint8_t  type;      // SCOPE_PAIR_MSG
  uint8_t  player;    // 1-4 (CONFIRM); copied back in OK / ricopiato nell'OK
  uint32_t nonce;     // random number of the Confirm / numero casuale della Conferma
} Scope_Pair_Msg;

// ===================================================================================
// VIDEO PACKETS (UDP) / PACCHETTI VIDEO (UDP)
// ===================================================================================
// EN: One slice per UDP packet: 8 bytes of header, then at most 1400 bytes of compressed data.
//     The dongle forwards the packets of the PC without touching the data: it writes only the session (byte 7)
//     and the start of session flag (bit 3 of byte 1).
// IT: Una fascia per pacchetto UDP: 8 byte di intestazione, poi al massimo 1400 byte di dati compressi.
//     Il dongle inoltra i pacchetti del PC senza toccare i dati: scrive solo la sessione (byte 7)
//     e il flag di inizio sessione (bit 3 del byte 1).
#define OPENFIRE_SCOPE_VIDEO_MARKER        0x4D   // 'M'
#define OPENFIRE_SCOPE_VIDEO_VERSION       1      // 4 high bits of byte 1 / 4 bit alti del byte 1
#define OPENFIRE_SCOPE_VIDEO_FLAG_SESSION  0x08   // bit 3 of byte 1: start of session / bit 3 del byte 1: inizio sessione
#define OPENFIRE_SCOPE_VIDEO_HEADER        8
#define OPENFIRE_SCOPE_VIDEO_MAX_DATA      1400
#define OPENFIRE_SCOPE_VIDEO_MAX_PACKET    (OPENFIRE_SCOPE_VIDEO_HEADER + OPENFIRE_SCOPE_VIDEO_MAX_DATA)
#define OPENFIRE_SCOPE_VIDEO_SLICES        15
#define OPENFIRE_SCOPE_VIDEO_PRESENCE      0xFF   // slice index of the presence packet / indice di fascia del pacchetto di presenza

// EN: Header byte positions. / IT: posizione dei byte dell'intestazione.
#define SCOPE_VIDEO_B_MARKER    0
#define SCOPE_VIDEO_B_FLAGS     1   // version (high) and flags (low) / versione (alti) e flag (bassi)
#define SCOPE_VIDEO_B_FRAME     2   // 2-3, little-endian
#define SCOPE_VIDEO_B_SLICE     4
#define SCOPE_VIDEO_B_Q         5
#define SCOPE_VIDEO_B_PLAYER    6
#define SCOPE_VIDEO_B_SESSION   7

// EN: Dongle times and sizes. / IT: tempi e dimensioni del dongle.
#define OPENFIRE_SCOPE_DONGLE_QUEUE          32    // slices waiting for the radio (2 frames) / fasce in attesa della radio (2 frame)
#define OPENFIRE_SCOPE_PRESENCE_IDLE_MS      200   // no video from the PC for 200 ms -> presence / niente video dal PC per 200 ms -> presenza
#define OPENFIRE_SCOPE_PRESENCE_EVERY_MS     200   // one presence every 200 ms / una presenza ogni 200 ms
#define OPENFIRE_SCOPE_PAUSE_VALID_MS        500   // q = 1 only with a PAUSE in the last 500 ms / q = 1 solo con una PAUSE negli ultimi 500 ms
#define OPENFIRE_SCOPE_SESSION_FLAG_MS       1000  // bit 3 during the first second of video / bit 3 nel primo secondo di video
#define OPENFIRE_SCOPE_TX_RETRY              5     // attempts, 1 ms apart, with the Wi-Fi queue full / tentativi, a 1 ms, con la coda Wi-Fi piena

// EN: Scope sizes: ring of the received slices in bytes, each slice followed by 4 bytes 0xFF for the bit reader of
//     the decoder and taking only its length: 12 slices of the maximum size, or many more of the usual ones.
// IT: Dimensioni del mirino: anello delle fasce ricevute in byte, ogni fascia seguita da 4 byte 0xFF per il lettore di
//     bit del decoder e grande solo quanto la sua lunghezza: 12 fasce della misura massima, o molte di piu' di quelle normali.
#define OPENFIRE_SCOPE_RX_RING               12    // slices of the maximum size, 17 KB: decoded as they arrive / fasce della misura massima, 17 KB: decodificate appena arrivano
#define OPENFIRE_SCOPE_RX_MARGIN             4
#define OPENFIRE_SCOPE_RX_SLOT               (OPENFIRE_SCOPE_VIDEO_MAX_PACKET + OPENFIRE_SCOPE_RX_MARGIN)   // 1412
#define OPENFIRE_SCOPE_RX_BYTES              (OPENFIRE_SCOPE_RX_RING * OPENFIRE_SCOPE_RX_SLOT)              // 16944
#define OPENFIRE_SCOPE_RX_PACKETS            64    // slices in the ring at most / fasce nell'anello al massimo
#define OPENFIRE_SCOPE_RX_TASK_PRIO          19    // above lwIP (18), below the Wi-Fi (23) / sopra lwIP (18), sotto il Wi-Fi (23)

// ===================================================================================
// LINK STATE / STATO DEL COLLEGAMENTO
// ===================================================================================
// EN: Same numbers as the "Scope state" phase sent to the S3. / IT: stessi numeri della fase di "Stato del mirino" mandata all'S3.
enum SCOPE_LINK_STATE {
  SCOPE_LINK_WAIT   = 0,   // waiting (scope: no player yet; dongle: channel search) / in attesa (mirino: niente giocatore; dongle: ricerca canale)
  SCOPE_LINK_SEARCH = 1,   // scope: searching a dongle; dongle: free, access point visible / mirino: cerca il dongle; dongle: libero, access point visibile
  SCOPE_LINK_PAIRED = 2,   // paired / abbinato
  SCOPE_LINK_VIDEO  = 3    // scope: paired and video on screen / mirino: abbinato e video sullo schermo
};

// ===================================================================================
// CLASS / CLASSE
// ===================================================================================
class WirelessScope_ {
  public:

  volatile uint8_t link_state = SCOPE_LINK_STATE::SCOPE_LINK_WAIT;
  uint8_t channel = OPENFIRE_WIFI_CHANNEL;     // channel of the pair / canale della coppia
  uint8_t player = 0;                          // 1-4, fixed at pairing / fissato all'abbinamento
  uint8_t mac_local[6];                        // MAC of this board / MAC di questa scheda
  uint8_t mac_peer[6];                         // MAC of the other board / MAC dell'altra scheda
  IPAddress ip_peer;                           // IP of the other board / IP dell'altra scheda

  // EN: Common radio start: Wi-Fi 6, 5 GHz only, 20 MHz, power, no power save.
  //     The dongle also chooses the channel and starts its access point.
  // IT: Avvio comune della radio: Wi-Fi 6, solo 5 GHz, 20 MHz, potenza, niente risparmio energetico.
  //     Il dongle sceglie anche il canale e avvia il suo access point.
  void begin();

  #ifdef SCOPE_DONGLE
  // EN: Listens to the 4 channels and returns the freest one. skipMask: bit i = 1 excludes channel i of
  //     OPENFIRE_SCOPE_CHANNELS (for RESCAN;SKIP from the PC daemon).
  // IT: Ascolta i 4 canali e restituisce il piu' libero. skipMask: bit i = 1 esclude il canale i di
  //     OPENFIRE_SCOPE_CHANNELS (per RESCAN;SKIP dal demone del PC).
  uint8_t findBestChannel(uint8_t skipMask = 0);
  bool startAccessPoint(uint8_t ch);
  // EN: Blocking: waits for a scope's Confirm, answers OK and becomes busy. Returns false at a RESCAN request.
  // IT: Bloccante: aspetta la Conferma di un mirino, risponde OK e diventa occupato. Torna false a una richiesta di RESCAN.
  bool connection_scope_dongle();

  // EN: RESCAN;SKIP from the PC daemon, only while free: requestRescan() from the USB, rescanChannel() from the main
  //     (blocking: new channel search without the skipped channels, then the access point again).
  // IT: RESCAN;SKIP dal demone del PC, solo da libero: requestRescan() dalla USB, rescanChannel() dal main
  //     (bloccante: nuova ricerca del canale senza i canali esclusi, poi di nuovo l'access point).
  volatile bool    rescan_request = false;
  volatile uint8_t rescan_skip_mask = 0;
  volatile uint8_t rescan_count = 0;   // +1 at every rescan done, the USB tells the channel again / +1 a ogni rescan fatto, la USB ricomunica il canale
  bool requestRescan(uint8_t skipMask);
  bool rescanChannel();

  // EN: Video: the USB passes every packet of the PC to videoPush(); a task of the library sends it to the scope.
  // IT: Video: la USB passa ogni pacchetto del PC a videoPush(); un task della libreria lo manda al mirino.
  volatile uint8_t  video_session = 1;    // 1-255, +1 at every START / 1-255, +1 a ogni START
  volatile uint32_t video_tx_bytes = 0;   // bytes passed to the Wi-Fi stack (OF_STAT TX) / byte passati allo stack Wi-Fi (OF_STAT TX)
  volatile uint32_t video_drop = 0;       // slices of old frames discarded (OF_STAT DROP) / fasce di frame vecchi scartate (OF_STAT DROP)
  volatile uint32_t video_qfull = 0;      // Wi-Fi transmit queue full (OF_STAT QFULL) / coda di trasmissione Wi-Fi piena (OF_STAT QFULL)
  bool videoPush(const uint8_t *pkt, uint16_t len);
  void videoStart();    // START from the PC / START dal PC
  void videoPause();    // PAUSE from the PC / PAUSE dal PC
  #endif // SCOPE_DONGLE

  #ifdef SCOPE
  // EN: Blocking: searches the first free dongle, connects and sends the player; returns when paired.
  // IT: Bloccante: cerca il primo dongle libero, si collega e comunica il giocatore; torna quando e' abbinato.
  bool connection_scope(uint8_t player_number);

  // EN: Video reception: a task of the library receives the UDP packets of the dongle and puts in the ring only the valid
  //     slices (marker, version, flags 0-2 at 0, player, slice 0-14, q 1-100, session 1-255). The presence packets do
  //     not go in the ring: they only update video_pause_ms (q = 1). Session and frame are checked by the decoder.
  // IT: Ricezione del video: un task della libreria riceve i pacchetti UDP del dongle e mette nell'anello solo le fasce
  //     valide (marcatore, versione, flag 0-2 a 0, giocatore, fascia 0-14, q 1-100, sessione 1-255). I pacchetti di
  //     presenza non vanno nell'anello: aggiornano solo video_pause_ms (q = 1). Sessione e frame li controlla il decoder.
  volatile uint32_t video_rx_slices = 0;      // valid slices put in the ring / fasce valide messe nell'anello
  volatile uint32_t video_rx_full = 0;        // valid slices lost with the ring full / fasce valide perse con l'anello pieno
  volatile uint32_t video_rx_bad = 0;         // packets discarded / pacchetti scartati
  volatile bool     video_slice_seen = false; // a valid slice since the pairing / una fascia valida dall'abbinamento
  volatile uint32_t video_slice_ms = 0;       // millis() of the last valid slice / millis() dell'ultima fascia valida
  volatile bool     video_pause_seen = false; // a presence with q = 1 since the pairing / una presenza con q = 1 dall'abbinamento
  volatile uint32_t video_pause_ms = 0;       // millis() of the last presence with q = 1 / millis() dell'ultima presenza con q = 1

  // EN: Number of a slice in the ring (-1 none within wait_ms), in order of arrival. After the decoding its space goes
  //     back with videoRelease().
  // IT: Numero di una fascia nell'anello (-1 nessuna entro wait_ms), in ordine di arrivo. Dopo la decodifica il suo spazio
  //     torna libero con videoRelease().
  int videoReceive(uint32_t wait_ms);
  const uint8_t *videoSlot(int slot, uint16_t *len);   // whole packet, header included / pacchetto intero, intestazione compresa
  void videoRelease(int slot);

  // EN: The video module tells when the video is on screen (phase 3 for the S3). Only PAIRED <-> VIDEO: never over a lost link.
  // IT: Il modulo video comunica quando il video e' sullo schermo (fase 3 per l'S3). Solo PAIRED <-> VIDEO: mai sopra un collegamento perso.
  void videoActive(bool on);
  #endif // SCOPE

  private:
  WiFiUDP udp_pair;
  void applyRadioSettings(wifi_interface_t ifx);
};

extern WirelessScope_ WirelessScope;

#endif // OPENFIRE_WIRELESS_ENABLE
