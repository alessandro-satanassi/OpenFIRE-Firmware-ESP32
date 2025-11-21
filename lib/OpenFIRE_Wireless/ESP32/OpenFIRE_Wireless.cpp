#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)

#include "OpenFIRE_Wireless.h"
#include "esp_idf_version.h"

#ifdef DONGLE
  #ifdef USES_DISPLAY
    #ifdef USE_LOVYAN_GFX
      #include <LovyanGFX.hpp>
      #include "../../src/LGFX_096_ST7735S_80x160.hpp"
      extern LGFX tft;
    #else
      #include <Adafruit_ST7735.h>
      extern Adafruit_ST7735 tft;
    #endif // USE_LOVYAN_GFX
  #endif // USES_DISPLAY
#endif // DONGLE

#ifdef GUN
  #ifdef USES_DISPLAY
    //#include "../../src/OpenFIREdisplay.h"
    //#include "../../src/OpenFIREcommon.h"
    #ifdef USE_LOVYAN_GFX
      #include <LovyanGFX.hpp>
      #include "../../src/LGFX_096_SSD1306_64x128.hpp"
      extern LGFX_SSD1306 *display_OLED;
    #else
      #include <Adafruit_SSD1306.h>
      extern Adafruit_SSD1306 *display_OLED;
      //extern ExtDisplay OLED;
    #endif // USE_LOVYAN_GFX
  #endif // USES_DISPLAY
#endif // GUN



#ifndef OPENFIRE_ESPNOW_WIFI_CHANNEL
  #define OPENFIRE_ESPNOW_WIFI_CHANNEL 12 // canale sul quale si sintonizza la lightgun di default
#endif // OPENFIRE_ESPNOW_WIFI_CHANNEL

#ifndef OPENFIRE_ESPNOW_WIFI_POWER
  // la potenza di trasmissione può andare da 8 a 84, dove 84 è il valore massimo che corrisponde a 20 db
  #define OPENFIRE_ESPNOW_WIFI_POWER 84 
#endif //OPENFIRE_ESPNOW_WIFI_POWER

uint8_t espnow_wifi_channel = OPENFIRE_ESPNOW_WIFI_CHANNEL;  // FATTA VARIABILE PER FUTURA CONFIGURAZIONE TRAMITE APP O OLED
uint8_t espnow_wifi_power = OPENFIRE_ESPNOW_WIFI_POWER;      // FATTA VARIABILE PER FUTURA CONFIGURAZIONE TRAMITE APP O OLED


USB_Data_GUN_Wireless usb_data_wireless = {
  "OpenFIRE_DONGLE",  // MANIFACTURES
  "FIRECon",          // NAME
  0xF143,             // VID
  0x1998, // 0x0001   // PID
  1,                  // PLAYER
  espnow_wifi_channel
  //OPENFIRE_ESPNOW_WIFI_CHANNEL // CHANNEL
  //,""               // ????
};

SerialWireless_ SerialWireless;
extern Adafruit_USBD_HID usbHid;
SemaphoreHandle_t tx_sem = NULL; // usato per callbalck dio fine trasmissione espnow


SemaphoreHandle_t mutex_tx_serial = NULL; // usato per trasmettere buffer seriale
SemaphoreHandle_t mutex_writer_bin = NULL; // usato per trasmettere buffer seriale

uint8_t buffer_espnow[ESP_NOW_MAX_DATA_LEN];
esp_now_peer_info_t peerInfo; // deve stare fuori funzioni da funzioni -- globale --variabile di utilità per configurazione

static void _esp_now_rx_cb(const esp_now_recv_info_t *info, const uint8_t *data, int len); // callback esp_now


#if ESP_IDF_VERSION_MAJOR == 5 && ESP_IDF_VERSION_MINOR <= 4
  // Codice per versioni fino alla 5.4
  static void _esp_now_tx_cb(const uint8_t *mac_addr, esp_now_send_status_t status); // callback esp_now
#elif ESP_IDF_VERSION_MAJOR > 5 || (ESP_IDF_VERSION_MAJOR == 5 && ESP_IDF_VERSION_MINOR >= 5)
  // Codice per versioni 5.5 o superiori
  static void _esp_now_tx_cb(const esp_now_send_info_t *tx_info, esp_now_send_status_t status); // callback esp_now // poer nuova versione IDF 55
#endif


void packet_callback_read_dongle(); // callback packet 
void packet_callback_read_gun(); // callback packet

///////////////////////////////////////////////////////////////////
//#define SIZE_BASE_AUX 13
uint8_t aux_buffer[13 + sizeof(usb_data_wireless)];
// 1 BYTE INDICATA L'AVANZAMENTO DELLO STABILIMENTO DELLA CONNESSIONE,A CHE PUNTO SIAMO
// 6 BYTE IL MAC ADDRESS DEL DISPOSITIVO CHE TRASMETTE
// 6 BYTE IL MAC ADDRESS DEL DISPOSITIVO A CUI E' INDIRIZZATO IL PACCHETTO
// ALL'ULTIMO INVIO LA GUN INVIA ANCHE I DATI DELLA STRUTTURA USB_DATA_WIRELESS CHE CONTIENE VID, PID, ECC. DELLA GUN IN MODO CHE IL DONGLE INZIALIZZI LA CONSSESIONE USB CON QUEI DATI


//uint8_t mac_esp_another_card[6];
//uint8_t stato_conneddione_wireless = 0;
uint8_t tipo_connessione = 0;

// stato della connessione 

// in caso di DONGLE
// 0 = non stabilita
// 1 = inviata richiesta blroad cast
// 2 = qualcuno ha risposto alla connessione ed ha chiesto di essere conesso
// 3 = inviato conferma connessione;

// in caso di GUN
// 0 = non stabilita
// 1 = arrivata richiesta broadcast di connessione
// 2 = accettata ed inviata richiesta di connessione
// 3 = arrivara accettazione di connessione

const uint8_t BROADCAST_ADDR[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
#if defined(DONGLE)
  ///////////////////////uint8_t peerAddress[6] = {0x24, 0x58, 0x7C, 0xDA, 0x38, 0xA0}; // quello montato con piedini su breackboard
  uint8_t peerAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}; // broadcast

  const functionPtr callbackArr[] = { packet_callback_read_dongle };
  uint8_t lastDongleAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  uint8_t lastDongleChannel = OPENFIRE_ESPNOW_WIFI_CHANNEL;
  bool lastDongleSave = false; // se true significa che abbiamo un indirizzo dell'ultimo dongle altrimenti false
#elif defined(GUN)
  //const uint8_t peerAddress[6] = {0xA0, 0x85, 0xE3, 0xE8, 0x0F, 0xB8}; // espe32s3 con piedini (riceve ma non trasmette)
  ///////////////////////////uint8_t peerAddress[6] = {0xA0, 0x85, 0xE3, 0xE7, 0x65, 0xD4}; // espe32s3 senz apiedini
  uint8_t peerAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}; // broadcast
  //const uint8_t peerAddress[6] = {0xF0, 0x9E, 0x9E, 0x28, 0x9E, 0xB8}; // DONGLE LILYGO

  const functionPtr callbackArr[] = { packet_callback_read_gun };

  // nel caso si spenga la pistola ed abbiamo memorizzato l'ultimo dongle connesso, prova a riconnettersi subito
  uint8_t lastDongleAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  uint8_t lastDongleChannel = OPENFIRE_ESPNOW_WIFI_CHANNEL;
  bool lastDongleSave = false; // se true significa che abbiamo un indirizzo dell'ultimo dongle altrimenti false
  // ==========================================================================================================

#endif
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
#if defined(GUN) && defined(USES_DISPLAY)
void animTaskLink(void *pvParameters) {
  unsigned long lastChange = 0;
  int8_t currentIndex = 0;
  int8_t direzione = 1;
  const char* word = "Scanning WiFi (((o)))"; // Waiting link 
  const char* rotazione = "-\\|/-\\|/";
  int8_t rotazioneIndex = 0;
  const uint8_t baseX = 1;
  const uint8_t baseY = 2;
  const uint8_t charWidth = 6;
  uint8_t len_word = strlen(word);

  // Setup iniziale display
  display_OLED->setCursor(baseX, baseY);
  display_OLED->setTextSize(1);
  display_OLED->setTextColor(WHITE, BLACK);
  display_OLED->fillRect(0, 0, 128, 16, BLACK);
  display_OLED->drawFastHLine(0, 15, 128, WHITE);
  char buffer[30];
  sprintf(buffer, "Ch:%2d Waiting link", espnow_wifi_channel);
  display_OLED->print(buffer);
  display_OLED->display();

  // Loop del task
  for (;;) {    
    display_OLED->setCursor(baseX + (19 * charWidth), baseY);
    display_OLED->write(rotazione[rotazioneIndex]);
    rotazioneIndex = (rotazioneIndex + 1) % 8;
    display_OLED->display();
    vTaskDelay(pdMS_TO_TICKS(150)); // piccolo delay per non saturare la CPU
  }
}
#endif // GUN


//////////////////////////////////////////////////////////////////

#if defined(GUN) && defined(USES_DISPLAY)
void animTask(void *pvParameters) {
  int8_t currentIndex = 0;
  const uint8_t baseX = 1;
  const uint8_t baseY = 2;
  const uint8_t charWidth = 6;

  // Setup iniziale display
  display_OLED->setCursor(baseX, baseY);
  display_OLED->setTextSize(1);
  display_OLED->setTextColor(WHITE, BLACK);
  display_OLED->fillRect(0, 0, 128, 16, BLACK);
  display_OLED->drawFastHLine(0, 15, 128, WHITE);
  display_OLED->print("Scanning WiFi        ");
  display_OLED->display();

  // Loop del task
  for (;;) {  
    display_OLED->setCursor(baseX + (14 * charWidth), baseY);
    switch (currentIndex)
    {
    case 0:
      /* code */
      display_OLED->print("   o   ");
      break;
    case 1:
      /* code */
      display_OLED->print("  (o)  ");
      break;
    case 2:
      /* code */
      display_OLED->print(" ((o)) ");
      break;
    case 3:
      /* code */
      display_OLED->print("(((o)))");
      break;
    
    default:
      
    break;
    }
    currentIndex = (currentIndex + 1) % 4;
    display_OLED->display();
    vTaskDelay(pdMS_TO_TICKS(150)); // piccolo delay per non saturare la CPU
  }
}
#endif //GUN

// ===============================================================
// ESP-NOW OPTIMAL CHANNEL FINDER - VERSIONE PERFETTA
// Tempo totale: ~14 secondi | Accuratezza: massima
// ===============================================================

// ================= VARIABILI GLOBALI PER CALLBACK =================
static volatile uint32_t g_packetCounter = 0;
static volatile bool g_sniffing = false;

// ================= CALLBACK PROMISCUOUS MODE =================
void IRAM_ATTR promiscuousCallback(void *buf, wifi_promiscuous_pkt_type_t type) {
  if (!g_sniffing) return;
  g_packetCounter = g_packetCounter + 1;
}

// ================= FUNZIONE PRINCIPALE =================
uint8_t findBestChannel() {
  
  #if defined(GUN) && defined(USES_DISPLAY)
  TaskHandle_t animTaskHandle = NULL;  
  // Avvio animazione
    if (animTaskHandle == NULL) {
      xTaskCreatePinnedToCore(
        animTask,          // funzione del task
        "AnimTask",        // nome
        4096,              // stack size
        NULL,              // parametri
        1,                 // priorità
        &animTaskHandle,   // handle
        APP_CPU_NUM        // core (puoi usare 0 o 1)
      );
    }
  #endif // USES_DISPLAY
  
  // ================= STRUTTURA PER STATISTICHE CANALE =================
  typedef struct {
    uint16_t networks;   // Numero di AP rilevati
    float avgRSSI;       // RSSI medio degli AP
    int8_t maxRSSI;      // RSSI massimo (interferenza più forte)
    uint32_t packets;    // Pacchetti catturati in 600ms
    float noise;         // Stima noise floor
    float score;         // Score base del canale
  } ChannelStats;
  
  ChannelStats channelStats[14];
  memset(channelStats, 0, sizeof(channelStats));
  g_packetCounter = 0;
  g_sniffing = false;
  for (int i = 0; i < 14; i++) channelStats[i].maxRSSI = -128;

  // ================= FASE 1: SCAN RETI WiFi =================
  // Tempo: ~2 secondi (150ms × 13 canali)
  WiFi.mode(WIFI_STA);
  WiFi.disconnect();
  vTaskDelay(pdMS_TO_TICKS(100));

  // Scan ottimizzato: 150ms per canale cattura tutti i beacon
  // (beacon interval tipico: 100ms)
  int n = WiFi.scanNetworks(false, true, false, 150UL, 0);
  if (n < 0) n = 0;

  if (n > 0) {
    for (int i = 0; i < n; i++) {
      int ch = WiFi.channel(i);
      int rssi = WiFi.RSSI(i);
      if (ch >= 1 && ch <= 13) {
        if (channelStats[ch].networks < 65535) channelStats[ch].networks++;
        channelStats[ch].avgRSSI += rssi;
        if (rssi > channelStats[ch].maxRSSI) {
          channelStats[ch].maxRSSI = (int8_t)rssi;
        }
      }
    }
    
    // Calcola RSSI medio
    for (int ch = 1; ch <= 13; ch++) {
      if (channelStats[ch].networks > 0) {
        channelStats[ch].avgRSSI /= channelStats[ch].networks;
      } else {
        channelStats[ch].maxRSSI = -100; // Canale vuoto
      }
    }
  }
  
  WiFi.scanDelete();
  vTaskDelay(pdMS_TO_TICKS(50));

  // ================= FASE 2: SNIFF TRAFFICO E NOISE =================
  // Tempo: ~12 secondi (920ms × 13 canali)
  
  vTaskDelay(pdMS_TO_TICKS(50));
  esp_wifi_set_promiscuous(true);
  esp_wifi_set_promiscuous_rx_cb(&promiscuousCallback);

  for (int ch = 1; ch <= 13; ch++) {

    esp_wifi_set_channel(ch, WIFI_SECOND_CHAN_NONE);
    vTaskDelay(pdMS_TO_TICKS(70)); // Stabilizzazione radio
    
    // --- Misurazione traffico: 600ms cattura pattern completi ---
    // Include: beacon (100ms), data burst, retry packets
    g_packetCounter = 0;
    g_sniffing = true;
    uint32_t startCount = g_packetCounter;
    vTaskDelay(pdMS_TO_TICKS(600));
    uint32_t endCount = g_packetCounter;
    channelStats[ch].packets = endCount - startCount;
    
    // --- Misurazione Noise Floor: 15 campioni per media stabile ---
    int32_t noiseSum = 0;
    const int samples = 15;
    for (int i = 0; i < samples; i++) {
      uint32_t start = g_packetCounter;
      vTaskDelay(pdMS_TO_TICKS(17)); // Campionamento 17ms (totale 255ms)
      uint32_t activity = g_packetCounter - start;
      
      // Formula noise migliorata: scala logaritmica
      if (activity == 0) {
        noiseSum += -100; // Silenzio completo
      } else if (activity <= 3) {
        noiseSum += -95 + (activity * 3);
      } else if (activity <= 10) {
        noiseSum += -86 + (activity - 3) * 2;
      } else if (activity <= 30) {
        noiseSum += -72 + (activity - 10);
      } else {
        noiseSum += -52 + (activity - 30) / 2;
      }
    }
    
    channelStats[ch].noise = noiseSum / (float)samples;
    channelStats[ch].noise = constrain(channelStats[ch].noise, -100, -50);
  }

  g_sniffing = false;
  esp_wifi_set_promiscuous(false);
  esp_wifi_set_promiscuous_rx_cb(NULL);

  // ================= FASE 3: CALCOLO SCORE OTTIMIZZATO =================
  // Pesi basati su analisi empirica e documentazione Espressif
  // 
  // PRIORITÀ PER ESP-NOW:
  // 1. RSSI forte = interferenza diretta che abbassa SNR (KILLER #1)
  // 2. Traffico reale = collisioni pacchetti (KILLER #2)
  // 3. Numero reti = congestione potenziale (moderato)
  // 4. Noise floor = basso impatto (ESP-NOW usa MCS0/802.11b robusto)
  
  const float wRSSI  = 0.40f;  // Interferenza diretta = priorità massima
  const float wPkt   = 0.35f;  // Collisioni reali = molto critico
  const float wNet   = 0.20f;  // Congestione = moderato (scala log)
  const float wNoise = 0.05f;  // Rumore = minimo (modulazione robusta)

  for (int ch = 1; ch <= 13; ch++) {
    
    // --- Score RSSI: interferenza diretta più pericolosa ---
    // AP con segnale forte (-40dBm) vicino = disastro garantito
    // -100dBm (lontano) = 0 punti | -40dBm (vicino) = 40 punti
    float sRSSI = 0;
    if (channelStats[ch].maxRSSI > -100) {
      int rssiMapped = constrain(channelStats[ch].maxRSSI, -100, -30);
      sRSSI = map(rssiMapped, -100, -30, 0, 100) * wRSSI;
    }
    
    // --- Score Pacchetti: traffico reale misurato ---
    // Normalizzato su 600ms: 0 pkt = 0 | 100 pkt = 35 punti
    float sPkt = (channelStats[ch].packets / 10.0) * wPkt;
    if (sPkt > 100 * wPkt) sPkt = 100 * wPkt; // clamp max
    
    // --- Score Reti: scala logaritmica per evitare dominanza ---
    // 0 reti = 0 | 1 rete = 4.2 | 5 reti = 9.2 | 10 reti = 13.8
    float sNet = (channelStats[ch].networks > 0) 
                 ? (log(channelStats[ch].networks + 1) * 20.0) * wNet 
                 : 0;
    
    // --- Score Noise: peso minimo (ESP-NOW resistente) ---
    // -100dBm (silenzio) = 0 | -50dBm (rumoroso) = 5 punti
    int noiseMapped = constrain((int)channelStats[ch].noise, -100, -50);
    float sNoise = map(noiseMapped, -100, -50, 0, 100) * wNoise;
    
    channelStats[ch].score = sRSSI + sPkt + sNet + sNoise;
  }

  // ================= PENALITÀ OVERLAP CANALI ADIACENTI =================
  // WiFi 20MHz: ogni canale si sovrappone con ±4 canali
  // Riferimento: 802.11b usa 22MHz, 802.11g/n usa 20MHz + guardband
  // Canali non-overlapping: 1, 6, 11 (separati di 5 canali)
  //
  // Penalità progressive per simulare overlap reale:
  // ±1: 50% (overlap massimo, circa 75% spettro condiviso)
  // ±2: 35% (overlap significativo, circa 50% spettro)
  // ±3: 20% (overlap moderato, circa 25% spettro)
  // ±4: 10% (overlap minimo, bordi spettro)
  
  float finalScores[14];
  for (int ch = 1; ch <= 13; ch++) {
    finalScores[ch] = channelStats[ch].score;
    
    // Penalità overlap completo modello 20MHz
    if (ch > 1) finalScores[ch] += channelStats[ch-1].score * 0.50f;
    if (ch > 2) finalScores[ch] += channelStats[ch-2].score * 0.35f;
    if (ch > 3) finalScores[ch] += channelStats[ch-3].score * 0.20f;
    if (ch > 4) finalScores[ch] += channelStats[ch-4].score * 0.10f;
    
    if (ch < 13) finalScores[ch] += channelStats[ch+1].score * 0.50f;
    if (ch < 12) finalScores[ch] += channelStats[ch+2].score * 0.35f;
    if (ch < 11) finalScores[ch] += channelStats[ch+3].score * 0.20f;
    if (ch < 10) finalScores[ch] += channelStats[ch+4].score * 0.10f;
  
    // Dopo il calcolo finalScores, dare bonus ai canali ideali
    //if (ch == 1 || ch == 6 || ch == 11) {
    //  finalScores[ch] *= 0.92f; // -8% bonus (preferenza)
    //}  
  }

  // ================= FASE 4: SELEZIONE MIGLIOR CANALE =================
  uint8_t bestCh = 1;
  float minScore = 999999.0;
  
  for (uint8_t ch = 1; ch <= 13; ch++) {
    if (finalScores[ch] < minScore) {
      minScore = finalScores[ch];
      bestCh = ch;
    }
  }
  
  // ================= PULIZIA FINALE =================
  vTaskDelay(pdMS_TO_TICKS(50));
  g_packetCounter = 0;
  g_sniffing = false;
  
  #if defined(GUN) && defined(USES_DISPLAY)
    if (animTaskHandle != NULL) {
      vTaskDelete(animTaskHandle);
      animTaskHandle = NULL;
    }
  #endif // USES_DISPLAY
  
  return bestCh;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*****************************
 *   SERIAL WIRELESS SECTION
 *****************************/

SerialWireless_::operator bool() {
  return true;
}

int SerialWireless_::peek() {
  if (lenBufferSerialRead) {
    return bufferSerialRead[_readerSerialRead];
  }
  return -1;
}

int SerialWireless_::peekBin() {
  if (_readLen) {
    return _queue[_reader];
  }
  return -1;
}

int SerialWireless_::read() {
  if (lenBufferSerialRead) {
    uint8_t ret = bufferSerialRead[_readerSerialRead];
    //lenBufferSerialRead --;
    lenBufferSerialRead = lenBufferSerialRead - 1;
    //_readerSerialRead ++;
    _readerSerialRead = _readerSerialRead + 1;
    if (_readerSerialRead >= FIFO_SIZE_READ_SERIAL) {
      _readerSerialRead -= FIFO_SIZE_READ_SERIAL;
    }
    return ret;
  }
  return -1;  
}

bool SerialWireless_::checkForRxPacket() {
  uint16_t numAvailableBin = _readLen;
  uint8_t dato;
  for (uint16_t i = 0; i<numAvailableBin; i++) {
    dato = (uint8_t) readBin();
    if (dato == START_BYTE) packet.reset(); //resetta inizio pacchetto - // controllo dato .. se è uguale a packet::start_byte .. azzera tutto e fai partire da capo altrimenti
    packet.parse(dato, true);
  }
  return true;
}

int SerialWireless_::readBin() {
  if (_readLen) {
    uint8_t ret = _queue[_reader];
    //_readLen --;
    _readLen = _readLen - 1;
    //_reader ++;
    _reader = _reader + 1;
    if (_reader >= FIFO_SIZE_READ) {
      _reader -= FIFO_SIZE_READ;
    }
    return ret;
  }
  return -1;  
}

int SerialWireless_::available() {
  return lenBufferSerialRead;
}

int SerialWireless_::availableBin() {
  return _readLen;
}

int SerialWireless_::availableForWrite() {
  return FIFO_SIZE_WRITE_SERIAL - lenBufferSerialWrite;
}

int SerialWireless_::availableForWriteBin() {
  return BUFFER_SIZE - _writeLen;
}

// bloccante fino a quando il buffer seriale non è vuoto (non bloccante fino a quando sono stati inviati effettivamente i dati)
void SerialWireless_::flush() {
  esp_timer_stop(timer_handle_serial); // spegni timer
  xSemaphoreTake(mutex_tx_serial, portMAX_DELAY); // prende semaforo / attende finchè è libero
  while (lenBufferSerialWrite) {
    flush_sem();
    taskYIELD(); //yield();
  }
  
  xSemaphoreGive(mutex_tx_serial);  // Rilascio il semaforo dopo la callback
  
}

// decidere se controllare se è vuoto anche il buffer_bin
bool SerialWireless_::flush_sem() { // è bloccante e non esce fino a quando il buffer di uscita è completamente vuoto // in virtual com con TinyUISB non è bloccante .. invia il pacchetto più grande che può e ritorna
  if (lenBufferSerialWrite) {
    if ((BUFFER_SIZE - _writeLen) > (lenBufferSerialWrite + PREAMBLE_SIZE + POSTAMBLE_SIZE)) {
      memcpy(&packet.txBuff[PREAMBLE_SIZE], bufferSerialWrite, lenBufferSerialWrite);
      packet.constructPacket(lenBufferSerialWrite, PACKET_TX::SERIAL_TX);
      writeBin(packet.txBuff, lenBufferSerialWrite + PREAMBLE_SIZE+POSTAMBLE_SIZE);
      lenBufferSerialWrite = 0;
      SendData(); // completata la transizione lasciare solo questo
      return true;
    }
    else SendData();
  }
  return false;
}

void SerialWireless_::flushBin() { // mai usata
  SendData();
}

void SerialWireless_::SendData() {
  if (xSemaphoreTake(tx_sem, 0) == pdTRUE) {
    xSemaphoreTake(mutex_writer_bin, portMAX_DELAY);
    if (_writeLen > 0) {
       SendData_sem();
    } else xSemaphoreGive(tx_sem);
    xSemaphoreGive(mutex_writer_bin);
  }
}

void SerialWireless_::SendData_sem() {
  uint16_t len_tx = _writeLen > ESP_NOW_MAX_DATA_LEN ? ESP_NOW_MAX_DATA_LEN : _writeLen;
  uint16_t bytesToSendEnd = len_tx > (BUFFER_SIZE - readIndex) ? BUFFER_SIZE - readIndex : len_tx;
  memcpy(buffer_espnow, &buffer[readIndex], bytesToSendEnd);
  if (len_tx > bytesToSendEnd) {
    memcpy(&buffer_espnow[bytesToSendEnd], &buffer[0], len_tx - bytesToSendEnd);
  }
  esp_err_t result = esp_now_send(peerAddress, buffer_espnow, len_tx);
  if (result == ESP_OK) { // verificare se esistono casi in cui seppur non risporta ESP_OK vengono trasmessi lo steso
    readIndex += len_tx; 
    if (readIndex >= BUFFER_SIZE) { readIndex -= BUFFER_SIZE; }
    _writeLen -= len_tx;
  }
}

size_t SerialWireless_::write(uint8_t c) {
  return write(&c, 1);
}

size_t SerialWireless_::writeBin(uint8_t c) {
  return writeBin(&c, 1);
}

size_t SerialWireless_::write(const uint8_t *data, size_t len) {
  size_t len_remainer = len;
  size_t pos_remainer = 0;

  xSemaphoreTake(mutex_tx_serial, portMAX_DELAY);
  while (len_remainer > 0) {
    size_t available_space = FIFO_SIZE_WRITE_SERIAL - lenBufferSerialWrite;
    if (available_space >= len_remainer) {
      memcpy(&bufferSerialWrite[lenBufferSerialWrite], &data[pos_remainer], len_remainer);
      if (lenBufferSerialWrite == 0) {
        esp_timer_start_once(timer_handle_serial, TIMER_HANDLE_SERIAL_DURATION_MICROS);
      }
      lenBufferSerialWrite += len_remainer;
      len_remainer = 0;
    } else {
      if (len_remainer == len) esp_timer_stop(timer_handle_serial);
      if (available_space) {
        memcpy(&bufferSerialWrite[lenBufferSerialWrite], &data[pos_remainer], available_space);
        lenBufferSerialWrite = FIFO_SIZE_WRITE_SERIAL;
        pos_remainer += available_space;
        len_remainer -= available_space;
      }
      flush_sem();
    }
    taskYIELD();
  }
  xSemaphoreGive(mutex_tx_serial);
  return len;
}

size_t SerialWireless_::writeBin(const uint8_t *data, size_t len) {
  
  xSemaphoreTake(mutex_writer_bin, portMAX_DELAY);
  
  if ((BUFFER_SIZE - _writeLen) >= len) {
    size_t firstChunk = BUFFER_SIZE - writeIndex;
    if (firstChunk < len) {
      memcpy(&buffer[writeIndex], data, firstChunk);
      writeIndex = len - firstChunk;  
      memcpy(&buffer[0], data + firstChunk, writeIndex);
    }
    else {
      memcpy(&buffer[writeIndex], data, len);
      writeIndex += len;
      if (writeIndex == BUFFER_SIZE) writeIndex = 0;
    }
    _writeLen += len;
    xSemaphoreGive(mutex_writer_bin);
    return len;
  }
  else {
    _overflow_write = true;
    // TODO: gestire overflow
    //Serial.println("Overflow nello scrivere nel BUFFER scrittura");   
    xSemaphoreGive(mutex_writer_bin);
    return 0;
  }
}

void SerialWireless_::SendPacket(const uint8_t *data, const uint8_t &len,const uint8_t &packetID) { 
  /*
  TODO: gestire fifo separate per controllo ??? valutare
  switch (packetID)
  {
  case PACKET_TX::MOUSE_TX:
    // gestire una fifo separata ?
    break;
  case PACKET_TX::KEYBOARD_TX:
    //
    break;
  case PACKET_TX::GAMEPADE_TX:
    //
    break;
  default:
    break;
  }
  */
  if ((BUFFER_SIZE - _writeLen) > (len + PREAMBLE_SIZE + POSTAMBLE_SIZE)) {
    memcpy(&packet.txBuff[PREAMBLE_SIZE], data, len);
    packet.constructPacket(len, packetID);
    writeBin(packet.txBuff, len + PREAMBLE_SIZE + POSTAMBLE_SIZE);
    SendData(); // try_Send
  }
}

void SerialWireless_::write_on_rx_serialBuffer(const uint8_t *data, int len) {
  if ((FIFO_SIZE_READ_SERIAL - lenBufferSerialRead) >= len) {
    size_t firstChunk = FIFO_SIZE_READ_SERIAL - _writerSerialRead;
    if (firstChunk < len) {
      memcpy(&bufferSerialRead[_writerSerialRead], data, firstChunk);
      _writerSerialRead = len - firstChunk;  
      memcpy(&bufferSerialRead[0], data + firstChunk, _writerSerialRead);
    }
    else {
      memcpy(&bufferSerialRead[_writerSerialRead], data, len);
      _writerSerialRead += len;
      if (_writerSerialRead == FIFO_SIZE_READ_SERIAL) _writerSerialRead = 0;
    }
    lenBufferSerialRead += len;
  }
  else {
    _overflow_bufferSerialRead = true;
    // TODO: gestire overflow
    //Serial.println("Overflow nello scrivere nel BUFFER lettura del SERIAL BUFFER");
  }
}

int SerialWireless_::availablePacket() {
  return numAvailablePacket;
}


void SerialWireless_::init_wireless() {
  configST myConfig; // variabile di utilità per configurazione
  // ============ inizializzazione semafori =============
  tx_sem = xSemaphoreCreateBinary();
  mutex_tx_serial = xSemaphoreCreateMutex();
  mutex_writer_bin = xSemaphoreCreateMutex();

  xSemaphoreGive(tx_sem);
  xSemaphoreGive(mutex_tx_serial);
  xSemaphoreGive(mutex_writer_bin);
  // ==== fine inizializzazione semafori

  myConfig.port         = &Serial; // questo andrà tolta - rimasta solo per contabilità =========================================
  myConfig.debug        = false; //true; //false; //true;
  myConfig.debugPort    = &Serial;
  myConfig.timeout      = DEFAULT_TIMEOUT; // 50ms
  myConfig.callbacks    = callbackArr;
  //myConfig.callbacks    = NULL;
  myConfig.callbacksLen = sizeof(callbackArr) / sizeof(functionPtr);
  //myConfig.callbacksLen = 0;
  packet.begin(myConfig);

  setupTimer(); // crea i timer .. timer per invio dati seriali

}


void SerialWireless_::begin() {
  #ifdef COMMENTO
  configST myConfig; // variabile di utilità per configurazione
  // ============ inizializzazione semafori =============
  tx_sem = xSemaphoreCreateBinary();
  mutex_tx_serial = xSemaphoreCreateMutex();
  mutex_writer_bin = xSemaphoreCreateMutex();

  xSemaphoreGive(tx_sem);
  xSemaphoreGive(mutex_tx_serial);
  xSemaphoreGive(mutex_writer_bin);
  // ==== fine inizializzazione semafori
  #endif //COMMENTO

  #ifdef GUN
  if (lastDongleSave) espnow_wifi_channel=lastDongleChannel;
    else {  
      #ifdef OPENFIRE_AUTO_CHANNEL_ESPNOW_WIFI
        espnow_wifi_channel = findBestChannel(); //12;
               
        #ifdef USES_DISPLAY
          display_OLED->setCursor(10, 2);
          display_OLED->setTextSize(1);
          display_OLED->setTextColor(WHITE, BLACK);
          display_OLED->fillRect(0, 0, 128, 16, BLACK);
          display_OLED->drawFastHLine(0, 15, 128, WHITE);
          display_OLED->print("Best Channel Ready");
          display_OLED->display();
          vTaskDelay(pdMS_TO_TICKS(500));
        #endif // USES_DISPLAY
    
      #endif // OPENFIRE_AUTO_CHANNEL_ESPNOW_WIFI
      
      #ifdef COMMENTO
      #ifdef USES_DISPLAY
        //FW_Common::OLED.TopPanelUpdate(" ... CONNECTION ...");
        //FW_Common::OLED.display->setTextSize(1);

        display_OLED->fillRect(0, 0, 128, 16, BLACK);
        display_OLED->drawFastHLine(0, 15, 128, WHITE);
        display_OLED->setCursor(2, 2);
        display_OLED->setTextSize(1);
        display_OLED->setTextColor(WHITE, BLACK);
        display_OLED->print("Ricerca canale");
        display_OLED->display();
      
      
        //FW_Common::OLED.ScreenModeChange(ExtDisplay::Screen_Init);
        //FW_Common::OLED.TopPanelUpdate(" ... CONNECTION ...");
      #endif // USES_DISPLAY
      //espnow_wifi_channel = findBestChannel(); //12;
      #ifdef USES_DISPLAY
        char buffer[50];
        sprintf(buffer, "Canale: %2d ", espnow_wifi_channel);
        // FW_Common::OLED.TopPanelUpdate(buffer);
        display_OLED->fillRect(0, 0, 128, 16, BLACK);
        display_OLED->drawFastHLine(0, 15, 128, WHITE);
        display_OLED->setCursor(2, 2);
        display_OLED->setTextSize(1);
        display_OLED->setTextColor(WHITE, BLACK);
        display_OLED->print(buffer);
        display_OLED->display();
      #endif // USES_DISPLAY
      #endif // COMMENTO

    }
  //findBestChannel();
  #endif //GUN

  WiFi.mode(WIFI_STA); 
  
  /*
  esp_err_t err = esp_wifi_start();
  if (err != ESP_OK) {
    log_e("WiFi not started! 0x%x)", err);
    return false;
  }
  */

  //WiFi.macAddress(mac_esp_inteface); // registra il mac addres della scheda

  esp_err_t err = esp_wifi_get_mac(WIFI_IF_STA, mac_esp_inteface);
  if (err != ESP_OK) {
    //Serial.println("Failed to read MAC address");
  }
  
  esp_wifi_set_promiscuous(true);
  err = esp_wifi_set_channel(espnow_wifi_channel, WIFI_SECOND_CHAN_NONE);
  if (err != ESP_OK) {
    //Serial.printf("esp_wifi_set_channel failed! 0x%x", err);
  }
  esp_wifi_set_promiscuous(false);

  err = esp_wifi_set_max_tx_power(espnow_wifi_power); // tra 8 e 84 corrispondenti a 2dbm a 20 dbm);
  if (err != ESP_OK) {
    //Serial.printf("esp_wifi_set_max_tx_power failed! 0x%x", err);
  }

  WiFi.disconnect();  // ??? si va messo

  vTaskDelay(pdMS_TO_TICKS(1000)); // delay(1000);
    
  err = esp_now_init();
  if (err != ESP_OK) {
    //Serial.printf("esp_now_init failed! 0x%x", err);
  }

  if (lastDongleSave) {
    memcpy(peerAddress, lastDongleAddress,6);
    stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION_DONGLE;
  } else {
    memcpy(peerAddress, BROADCAST_ADDR, 6);
    stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
  }
  
  memcpy(peerInfo.peer_addr, peerAddress, 6);
  peerInfo.channel = espnow_wifi_channel;
  peerInfo.encrypt = false;
  if (esp_now_add_peer(&peerInfo) != ESP_OK) {
    //Serial.println("Errore nell'aggiunta del peer");
  }

  err = esp_now_register_recv_cb(_esp_now_rx_cb);
  if (err != ESP_OK) {
    //Serial.printf("esp_now_register_recv_cb failed! 0x%x", err);
  }

  err = esp_now_register_send_cb(_esp_now_tx_cb);
  if (err != ESP_OK) {
    //Serial.printf("esp_now_register_send_cb failed! 0x%x", err);
  }

  #ifdef COMMENTO
  myConfig.port         = &Serial; // questo andrà tolta - rimasta solo per contabilità =========================================
  myConfig.debug        = false; //true; //false; //true;
  myConfig.debugPort    = &Serial;
  myConfig.timeout      = DEFAULT_TIMEOUT; // 50ms
  myConfig.callbacks    = callbackArr;
  //myConfig.callbacks    = NULL;
  myConfig.callbacksLen = sizeof(callbackArr) / sizeof(functionPtr);
  //myConfig.callbacksLen = 0;
  
  packet.begin(myConfig);
  #endif // COMMENTO

  TinyUSBDevices.wireless_mode = WIRELESS_MODE::ENABLE_ESP_NOW_TO_DONGLE;
  #ifdef COMMENTO
  setupTimer(); // crea i timer .. timer per invio dati seriali
  #endif // COMMENTO
}

bool SerialWireless_::end() {

  //esp_now_del_peer(peerAddress); // non serve lo fa direttamente esp_now_deinit()
  esp_err_t err = esp_now_deinit();
  if (err != ESP_OK) {
    //Serial.printf("esp_now_deinit failed! 0x%x", err);
    return false;
  }
  WiFi.disconnect(true);
  WiFi.mode(WIFI_OFF); // Disattiva il WiFi
  esp_timer_delete(timer_handle_serial);
  return true;
}

bool SerialWireless_::connection_dongle() {

  uint8_t channel = espnow_wifi_channel;   // tra e 1 e 13 (il 14 dovrebbe essere riservato)
  #define TIMEOUT_TX_PACKET 500 // in millisecondi
  #define TIMEOUT_CHANGE_CHANNEL 2000 // in millisecondi - cambia canale ogni
  #define TIMEOUT_DIALOGUE 6000 // in millisecondi - tempo massimo per completare operazione accoppiamento
  unsigned long lastMillis_tx_packet = millis ();
  unsigned long lastMillis_change_channel = millis ();
  unsigned long lastMillis_start_dialogue = millis ();
  uint8_t aux_buffer_tx[14]; // aggiunto un byte per trasmettere anche il canale di trasmissione
                             // durante tx pacchetto pubblicazione presenza, mette anche il canale
  
  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
  aux_buffer_tx[0] = CONNECTION_STATE::TX_DONGLE_SEARCH_GUN_BROADCAST;
  memcpy(&aux_buffer_tx[1], SerialWireless.mac_esp_inteface, 6);
  memcpy(&aux_buffer_tx[7], peerAddress, 6);
  //aux_buffer_tx[13] = espnow_wifi_channel;
  
  
  #ifdef DONGLE
    #ifdef USES_DISPLAY
      tft.fillRect(95,60,50,20,0/*BLACK*/);
      tft.setCursor(95, 60);  
      tft.printf("%2d", channel);
    #endif //USES_DISPLAY
  #endif //DONGLE
  
  
  while (stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
    if (stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
      if (((millis() - lastMillis_change_channel) > TIMEOUT_CHANGE_CHANNEL) && 
         ((millis() - (lastMillis_tx_packet-50))) > TIMEOUT_TX_PACKET) {  // aggiunta impostato 50 ms come margine, per evitare che quando invia pacchetto cambi subito casnale senza dare possibilità risposta
        channel++;
        if (channel >13) channel = 1;
        aux_buffer_tx[13] = channel;
        #ifdef DONGLE
          #ifdef USES_DISPLAY
            tft.fillRect(95,60,50,20,0/*BLACK*/);  
            tft.setCursor(95, 60);
            tft.printf("%2d", channel);   
          #endif //USES_DISPLAY
        #endif //DONGLE
        esp_wifi_set_promiscuous(true);
        if (esp_wifi_set_channel(channel, WIFI_SECOND_CHAN_NONE) != ESP_OK) {
          //Serial.printf("DONGLE - esp_wifi_set_channel failed!");
        }
        esp_wifi_set_promiscuous(false);
        peerInfo.channel = channel;
        if (esp_now_mod_peer(&peerInfo) != ESP_OK) {  // modifica il canale del peer
          //Serial.println("DONGLE - Errore nella modifica del canale");
        }             
        lastMillis_change_channel = millis ();
        lastMillis_tx_packet = 0; // per fare inviare subito un pacchetto sul nuovo canale
      }
      if ((millis() - lastMillis_tx_packet) > TIMEOUT_TX_PACKET) {
        SerialWireless.SendPacket((const uint8_t *)aux_buffer_tx, 14, PACKET_TX::CONNECTION); // aggiunto un byte per trasmettere anche il canale di trasmissione
        //Serial.print("DONGLE - inviato pacchetto broadcast sul canale: ");
        //Serial.println(channel);
        lastMillis_tx_packet = millis (); 
      }
      lastMillis_start_dialogue = millis();
    }
    else {
      if (((millis() - lastMillis_start_dialogue) > TIMEOUT_DIALOGUE) && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
        stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
        //Serial.println("DONGLE - Non si è conclusa la negoziazione tra DONGLE/GUN e si riparte da capo");
        lastMillis_change_channel = millis ();
      }  
    }
    yield(); // in attesa dello stabilimento di una connessione
  }
  
  //Serial.println("DONGLE - Negosazione completata - associazione dei dispositivi GUN/DONGLE");
  if (esp_now_del_peer(peerAddress) != ESP_OK) {  // cancella il broadcast dai peer
    //Serial.println("DONGLE - Errore nella cancellazione del peer broadcast");
  }
  memcpy(peerAddress, mac_esp_another_card, 6);
  memcpy(peerInfo.peer_addr, peerAddress, 6);
  espnow_wifi_channel=usb_data_wireless.channel;
  esp_wifi_set_promiscuous(true);
  esp_wifi_set_channel(espnow_wifi_channel, WIFI_SECOND_CHAN_NONE);
  esp_wifi_set_promiscuous(false);
  peerInfo.channel = espnow_wifi_channel;
  //peerInfo.channel = ESPNOW_WIFI_CHANNEL;
  //peerInfo.encrypt = false;              
  if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserisce il dongle nei peer
    //Serial.println("DONGLE - Errore nell'aggiunta del nuovo peer della GUN");
  }
  TinyUSBDevices.onBattery = true;
  return true;
}

bool SerialWireless_::connection_gun_at_last_dongle() {
  #define TIMEOUT_TX_PACKET_LAST_DONGLE 300 // in millisecondi - tempo di invio pacchetti ogni millisecondi quindi 4-5 pacchetti
  #define TIMEOUT_DIALOGUE_LAST_DONGLE 2000 // in millisecondi - tempo massimo per ricerca ultimo dongle
  unsigned long lastMillis_tx_packet_last_dongle = 0;
  unsigned long lastMillis_start_dialogue_last_dongle = millis ();

  uint8_t aux_buffer_tx[13];
  
  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION_DONGLE;
  aux_buffer_tx[0] = CONNECTION_STATE::TX_CHECK_CONNECTION_LAST_DONGLE;
  memcpy(&aux_buffer_tx[1], SerialWireless.mac_esp_inteface, 6);
  memcpy(&aux_buffer_tx[7], peerAddress, 6);
  lastMillis_tx_packet_last_dongle = 0;
  lastMillis_start_dialogue_last_dongle = millis();
  while (!TinyUSBDevice.mounted() && 
         (stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED_WITH_LAST_DONGLE) &&
         ((millis() - lastMillis_start_dialogue_last_dongle) < TIMEOUT_DIALOGUE_LAST_DONGLE)) { 
      if ((millis() - lastMillis_tx_packet_last_dongle) > TIMEOUT_TX_PACKET_LAST_DONGLE)
      {
        SerialWireless.SendPacket((const uint8_t *)aux_buffer_tx, 13, PACKET_TX::CHECK_CONNECTION_LAST_DONGLE); 
        lastMillis_tx_packet_last_dongle = millis();
      }
    yield();
  }
  if (stato_connessione_wireless == CONNECTION_STATE::DEVICES_CONNECTED_WITH_LAST_DONGLE) {    
    stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED;
    TinyUSBDevices.onBattery = true;
    return true;
  } else {
    stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
    lastDongleSave=false;
    esp_now_deinit();
    begin();
    return false;
  }
}


bool SerialWireless_::connection_gun() {
  
  #if defined(GUN) && defined(USES_DISPLAY)
  TaskHandle_t animTaskHandleLink = NULL;  
  // Avvio animazione
    if (animTaskHandleLink == NULL) {
      xTaskCreatePinnedToCore(
        animTaskLink,          // funzione del task
        "AnimTaskLink",        // nome
        4096,              // stack size
        NULL,              // parametri
        1,                 // priorità
        &animTaskHandleLink,   // handle
        APP_CPU_NUM        // core (puoi usare 0 o 1)
      );
    }
  #endif // USES_DISPLAY



  #define TIMEOUT_GUN_DIALOGUE 1000 // in millisecondi
  unsigned long lastMillis_start_dialogue = millis ();
  
  lastMillis_start_dialogue = millis ();
  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
  //IMPOSTARE IL PEER BOARCAST QUI ======================
  if (esp_now_del_peer(peerAddress) != ESP_OK) {  // cancella il broadcast dai peer
    //Serial.println("Errore nella cancellazione del peer");
  }
  memcpy(peerAddress, BROADCAST_ADDR, 6);
  memcpy(peerInfo.peer_addr, peerAddress, 6);
  //peerInfo.channel = ESPNOW_WIFI_CHANNEL;
  //peerInfo.encrypt = false;              
  if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserisce il dongle nei peer
    //Serial.println("Errore nell'aggiunta del nuovo peer");
  }                       
  // ====================================================

  while(!TinyUSBDevice.mounted() && stato_connessione_wireless != CONNECTION_STATE::TX_GUN_TO_DONGLE_CONFERM) { 
    if (stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
      lastMillis_start_dialogue = millis ();
    }
    if (((millis() - lastMillis_start_dialogue) > TIMEOUT_GUN_DIALOGUE) && stato_connessione_wireless != CONNECTION_STATE::TX_GUN_TO_DONGLE_CONFERM) {
      stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
    }
  }

  if (stato_connessione_wireless == CONNECTION_STATE::TX_GUN_TO_DONGLE_CONFERM) {
    if (esp_now_del_peer(peerAddress) != ESP_OK) {  // cancella il broadcast dai peer
      //Serial.println("Errore nella cancellazione del peer");
    }
    
    memcpy(peerAddress, mac_esp_another_card, 6);
    memcpy(peerInfo.peer_addr, peerAddress, 6);
    //peerInfo.channel = ESPNOW_WIFI_CHANNEL;
    //peerInfo.encrypt = false;              
    if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserisce il dongle nei peer
      //Serial.println("Errore nell'aggiunta del nuovo peer");
    }                       
    
    TinyUSBDevices.onBattery = true;

    #if defined(GUN) && defined(USES_DISPLAY)
      if (animTaskHandleLink != NULL) {
        vTaskDelete(animTaskHandleLink);
        animTaskHandleLink = NULL;
      }
    #endif // USES_DISPLAY

    return true;
  }
  
  #if defined(GUN) && defined(USES_DISPLAY)
    if (animTaskHandleLink != NULL) {
      vTaskDelete(animTaskHandleLink);
      animTaskHandleLink = NULL;
    }
  #endif // USES_DISPLAY
   
  return false;
}


void packet_callback_read_dongle() {
  switch (SerialWireless.packet.currentPacketID()) {
    case PACKET_TX::SERIAL_TX:
      Serial.write(&SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
      Serial.flush(); // ????
      break;
    case PACKET_TX::MOUSE_TX :
      usbHid.sendReport(HID_RID_e::HID_RID_MOUSE, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);  
      break;
    case PACKET_TX::GAMEPADE_TX:
      usbHid.sendReport(HID_RID_e::HID_RID_GAMEPAD, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
      break;
    case PACKET_TX::KEYBOARD_TX:
      usbHid.sendReport(HID_RID_e::HID_RID_KEYBOARD, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
      break; 
    case PACKET_TX::CHECK_CONNECTION_LAST_DONGLE:
      // CODICE
      memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
      switch (aux_buffer[0])
      {
      case CONNECTION_STATE::TX_CHECK_CONNECTION_LAST_DONGLE:
        if ((memcmp(&aux_buffer[1],peerAddress,6) == 0) && 
           (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
           SerialWireless.stato_connessione_wireless == CONNECTION_STATE::DEVICES_CONNECTED) { 
              aux_buffer[0] = CONNECTION_STATE::TX_CONFERM_CONNECTION_LAST_DONGLE; 
              memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
              memcpy(&aux_buffer[7], peerAddress, 6);
              // valutare se inviarlo un paio di volte il pacchetto o solo una volta
              // lo invia 3 volte - una volta ogni 70ms
              for (uint8_t i = 0; i<3; i++) {
                if (i>0) {
                  //vTaskDelay(pdMS_TO_TICKS(1000)); // equivalente a delay(1000) ma non bloccante su esp32
                  unsigned long lastMillis_tx_packet_connection_last_dongle = millis();
                  while ((millis() - lastMillis_tx_packet_connection_last_dongle) < 70) yield(); 
                }
                SerialWireless.SendPacket((const uint8_t *)aux_buffer, 13, PACKET_TX::CHECK_CONNECTION_LAST_DONGLE);
              }
        }
        break;
      default:
        break;
      }     
      break;
    case PACKET_TX::CONNECTION:
      memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead); //sizeof(aux_buffer));
      //Serial.println("DONGLE - arrivato richiesta di connessione");
      switch (aux_buffer[0]) {
        case CONNECTION_STATE::TX_GUN_TO_DONGLE_PRESENCE:
          if ((memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) && SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
            memcpy(SerialWireless.mac_esp_another_card, &aux_buffer[1], 6);
            aux_buffer[0] = CONNECTION_STATE::TX_DONGLE_TO_GUN_ACCEPT;
            memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
            memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
            SerialWireless.SendPacket((const uint8_t *)aux_buffer, 13, PACKET_TX::CONNECTION);
            //Serial.println("DONGLE - inviato pacchetto di conferma di connessione");

            // assicurati che i dati siano stati spediti
 
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::TX_DONGLE_TO_GUN_ACCEPT;
          }
          break;
        case CONNECTION_STATE::TX_GUN_TO_DONGLE_CONFERM:
          if ((memcmp(&aux_buffer[1],SerialWireless.mac_esp_another_card,6) == 0) && 
             (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
             SerialWireless.stato_connessione_wireless == CONNECTION_STATE::TX_DONGLE_TO_GUN_ACCEPT) {           
            // SALVA I DATI RELATIVI ALLA GUN VID, PID, PLAYER , ECC.ECC.
            memcpy(&usb_data_wireless, &aux_buffer[13], sizeof(usb_data_wireless));
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED;
            //Serial.println("DONGLE - ricevuto pacchetto di avvenuta connessione con dati della GUN");
        }
          break;
        default:
          break;
      }
      break;  
    default:
      break;
  }
}

void packet_callback_read_gun() {
  switch (SerialWireless.packet.currentPacketID()) {
    case PACKET_TX::SERIAL_TX:
      SerialWireless.write_on_rx_serialBuffer(&SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
      break;
    case PACKET_TX::MOUSE_TX :
      /* code */
      break;
    case PACKET_TX::GAMEPADE_TX:
      /* code */
      break;
    case PACKET_TX::KEYBOARD_TX:
      /* code */  
      break; 
    case PACKET_TX::CHECK_CONNECTION_LAST_DONGLE:
      // CODICE VALUTARE SE FARE CONTROLLO SE PISTOLA ANCORA CONNESSA .. SE NON CONNESSA RIAVVIA LA SCHEDA ?
      memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
      switch (aux_buffer[0])
      {
      case CONNECTION_STATE::TX_CONFERM_CONNECTION_LAST_DONGLE:
        if ((memcmp(&aux_buffer[1],peerAddress,6) == 0) && 
           (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
           SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION_DONGLE) { 
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED_WITH_LAST_DONGLE;
        }
        break;
      default:
        break;
      }     
      break;
    case PACKET_TX::CONNECTION:
      memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead); //13); //sizeof(aux_buffer)); // qui va bene anche 13 come dati da copiare
      switch (aux_buffer[0]) {
        case CONNECTION_STATE::TX_DONGLE_SEARCH_GUN_BROADCAST:
          if ((SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) &&
              (aux_buffer[13] == espnow_wifi_channel))
          { // prende la prima dongle disposnibile
            memcpy(SerialWireless.mac_esp_another_card, &aux_buffer[1], 6);
            // invia richiesta connessione
            aux_buffer[0] = CONNECTION_STATE::TX_GUN_TO_DONGLE_PRESENCE;
            memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
            memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
            SerialWireless.SendPacket((const uint8_t *)aux_buffer, 13, PACKET_TX::CONNECTION);
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::TX_GUN_TO_DONGLE_PRESENCE;
            // assicurati che i dati siano stati spediti

          }
          break;
        case 2:
          /* code */
          break;
        case CONNECTION_STATE::TX_DONGLE_TO_GUN_ACCEPT:
          if ((memcmp(&aux_buffer[1],SerialWireless.mac_esp_another_card,6) == 0) && 
             (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
              SerialWireless.stato_connessione_wireless == CONNECTION_STATE::TX_GUN_TO_DONGLE_PRESENCE) {
              aux_buffer[0] = CONNECTION_STATE::TX_GUN_TO_DONGLE_CONFERM;
              memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
              memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
              // INVIA ANCHE DATI RELATIVI A VID, PID, ECC,ECC, DELLA GUN
              usb_data_wireless.channel = espnow_wifi_channel;
              memcpy(&aux_buffer[13], &usb_data_wireless, sizeof(usb_data_wireless));

              // =========================================================
              // INVIARE IL PACCHETTO FINALE PIU' VOLTE  
              // =========================================================
              for (uint8_t i = 0; i<3; i++) {
                if (i>0) {
                  //vTaskDelay(pdMS_TO_TICKS(1000)); // equivalente a delay(1000) ma non bloccante su esp32
                  unsigned long lastMillis_tx_packet_gun_to_dongle_conferm = millis();
                  while ((millis() - lastMillis_tx_packet_gun_to_dongle_conferm) < 70) yield(); 
                }
                SerialWireless.SendPacket((const uint8_t *)aux_buffer, sizeof(aux_buffer), PACKET_TX::CONNECTION);
              }
              SerialWireless.stato_connessione_wireless = CONNECTION_STATE::TX_GUN_TO_DONGLE_CONFERM;
          }
          break;        
        default:
          break;
      }
      break;
    default:
      break;
  }
}

static void _esp_now_rx_cb(const esp_now_recv_info_t *info, const uint8_t *data, int len) {
  if ((FIFO_SIZE_READ - SerialWireless._readLen) >= len) {
    size_t firstChunk = FIFO_SIZE_READ - SerialWireless._writer;
    if (firstChunk < len) {
      memcpy(&SerialWireless._queue[SerialWireless._writer], data, firstChunk);
      SerialWireless._writer = len - firstChunk;  
      memcpy(&SerialWireless._queue[0], data + firstChunk, SerialWireless._writer);
    }
    else {
      memcpy(&SerialWireless._queue[SerialWireless._writer], data, len);
      SerialWireless._writer += len;
      if (SerialWireless._writer == FIFO_SIZE_READ) SerialWireless._writer = 0;
    }
    SerialWireless._readLen += len;
  }
  else {
    SerialWireless._overflow_read = true;
    //Serial.println("Overflow nello scrivere nel BUFFER lettura");
  }
  SerialWireless.checkForRxPacket();
}

#if ESP_IDF_VERSION_MAJOR == 5 && ESP_IDF_VERSION_MINOR <= 4
  // Codice per versioni fino alla 5.4
  static void _esp_now_tx_cb(const uint8_t *mac_addr, esp_now_send_status_t status) {
#elif ESP_IDF_VERSION_MAJOR > 5 || (ESP_IDF_VERSION_MAJOR == 5 && ESP_IDF_VERSION_MINOR >= 5)
  // Codice per versioni 5.5 o superiori
  static void _esp_now_tx_cb(const esp_now_send_info_t *tx_info, esp_now_send_status_t status) { // callback esp_now // poer nuova versione IDF 55
#endif

  xSemaphoreTake(mutex_writer_bin, portMAX_DELAY);

  if (SerialWireless._writeLen > 0 ) {
    SerialWireless.SendData_sem();
  }
  else xSemaphoreGive(tx_sem);
  
  xSemaphoreGive(mutex_writer_bin);
}

// ================================== TIMER ===================================
// CALLBACK
void timer_callback_serial(void* arg) {
  //Serial.println("Timer scaduto!");
    
  xSemaphoreTake(mutex_tx_serial, portMAX_DELAY);
  
  while (SerialWireless.lenBufferSerialWrite) {
    SerialWireless.flush_sem();
    taskYIELD(); //yield();
  }
  
  xSemaphoreGive(mutex_tx_serial);  // Rilascio il semaforo dopo la callback
}
// ===============================================================================

void SerialWireless_::setupTimer() {
    esp_timer_create_args_t timer_args = {
        .callback = &timer_callback_serial,
        .arg = NULL,
        .dispatch_method = ESP_TIMER_TASK,
        .name = "timer_serial"
    };
    
    esp_timer_create(&timer_args, &timer_handle_serial);
}

void SerialWireless_::stopTimer_serial() {
    esp_timer_stop(timer_handle_serial);
}

void SerialWireless_::resetTimer_serial(uint64_t duration_us) {
    esp_timer_stop(timer_handle_serial);
    esp_timer_start_once(timer_handle_serial, duration_us);
}


// ==========================   FINE TIMER ====================================

#endif //OPENFIRE_WIRELESS_ENABLE