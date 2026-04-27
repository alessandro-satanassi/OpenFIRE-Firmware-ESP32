

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)

#include "OpenFIRE_Wireless.h"
#include "esp_idf_version.h"

#ifdef DONGLE
  #ifdef USES_DISPLAY
    #ifdef USE_LOVYAN_GFX
      #define LGFX_USE_V1
      #include <LovyanGFX.hpp>
      #include "LGFX_096_ST7735S_80x160.hpp"
      extern LGFX tft;
      #define WHITE            TFT_WHITE
      #define BLACK            TFT_BLACK
      #define BLUE             TFT_BLUE
      #define RED              TFT_RED
      #define MAGENTA          TFT_MAGENTA
      #define GREEN            TFT_GREEN
      #define CYAN             TFT_CYAN
      #define YELLOW           TFT_YELLOW
      #define BROWN            TFT_BROWN
      #define GRAY             TFT_LIGHTGRAY
      #define ORANGE           TFT_ORANGE
    #else
      #include <Adafruit_ST7735.h>
      extern Adafruit_ST7735 tft;
      #define WHITE            ST7735_WHITE
      #define BLACK            ST7735_BLACK
      #define BLUE             ST7735_BLUE
      #define RED              ST7735_RED
      #define MAGENTA          ST7735_MAGENTA
      #define GREEN            ST7735_GREEN
      #define CYAN             ST7735_CYAN
      #define YELLOW           ST7735_YELLOW
      #define BROWN            0XBC40
      #define GRAY             0X8430
      #define ORANGE           ST7735_ORANGE
    #endif // USE_LOVYAN_GFX
  #endif // USES_DISPLAY
#endif // DONGLE

#ifdef GUN
  #ifdef USES_DISPLAY
    //#include "../../src/OpenFIREdisplay.h"
    //#include "../../src/OpenFIREcommon.h"
    #ifdef USE_LOVYAN_GFX
      #include <LovyanGFX.hpp>
      #include "LGFX_096_SSD1306_64x128.hpp"
      extern LGFX_SSD1306 *&display_OLED;
    #else
      #include <Adafruit_SSD1306.h>
      extern Adafruit_SSD1306 *&display_OLED;
      //extern ExtDisplay OLED;
    #endif // USE_LOVYAN_GFX
  #endif // USES_DISPLAY
#endif // GUN

extern bool display_init;


#ifndef OPENFIRE_ESPNOW_WIFI_CHANNEL
  #define OPENFIRE_ESPNOW_WIFI_CHANNEL 11 // canale sul quale si sintonizza la lightgun di default
#endif // OPENFIRE_ESPNOW_WIFI_CHANNEL

#ifndef OPENFIRE_ESPNOW_WIFI_POWER
  // la potenza di trasmissione può andare da 8 a 84, dove 84 è il valore massimo che corrisponde a 21 db (x4)
  //#define OPENFIRE_ESPNOW_WIFI_POWER WIFI_POWER_15dBm  // corrisponde a 60 .. potenza più conservativa ma con segnale eccellente
  #define OPENFIRE_ESPNOW_WIFI_POWER WIFI_POWER_17dBm  // corrisponde a 68 .. massima potenza senza distursioni con segnale eccellette
  
#endif //OPENFIRE_ESPNOW_WIFI_POWER

uint8_t espnow_wifi_channel = OPENFIRE_ESPNOW_WIFI_CHANNEL;  // FATTA VARIABILE PER FUTURA CONFIGURAZIONE TRAMITE APP O OLED
uint8_t espnow_wifi_power = OPENFIRE_ESPNOW_WIFI_POWER;      // FATTA VARIABILE PER FUTURA CONFIGURAZIONE TRAMITE APP O OLED


hid_abs_mouse_report_t absmouse5Report_last_wifi = {0,0,0,0,0};
hid_keyboard_report_t  keyReport_last_wifi = {0,0,{0,0,0,0,0,0}};
hid_gamepad16_report_t gamepad16Report_last_wifi = {0,0,0,0,0,0,0,0};


volatile bool absmouse5Report_pending = false; 
volatile bool keyReport_pending = false;
volatile bool gamepad16Report_pending = false;


#ifdef OPENFIRE_ESPNOW_WIFI_POWER_AUTO
  bool espnow_wifi_power_auto = true;
#endif // OPENFIRE_ESPNOW_WIFI_POWER_AUTO


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

portMUX_TYPE mux_write_bin = portMUX_INITIALIZER_UNLOCKED;

// Creazione del Mutex globale (da chiamare possibilmente in setup() se non lo crei globalmente)
//SemaphoreHandle_t txMutex = xSemaphoreCreateMutex();

static portMUX_TYPE mux_serial_tx = portMUX_INITIALIZER_UNLOCKED;
static portMUX_TYPE mux_radio_tx = portMUX_INITIALIZER_UNLOCKED;


uint8_t buffer_espnow[ESP_NOW_MAX_DATA_LEN];
esp_now_peer_info_t peerInfo; // deve stare fuori funzioni da funzioni -- globale --variabile di utilità per configurazione

static void _esp_now_rx_cb(const esp_now_recv_info_t *info, const uint8_t *data, int len); // callback esp_now
static void _esp_now_tx_cb(const esp_now_send_info_t *tx_info, esp_now_send_status_t status);

#ifdef COMMENTO
#if ESP_IDF_VERSION_MAJOR == 5 && ESP_IDF_VERSION_MINOR <= 4
  // Codice per versioni fino alla 5.4
  static void _esp_now_tx_cb(const uint8_t *mac_addr, esp_now_send_status_t status); // callback esp_now
#elif ESP_IDF_VERSION_MAJOR > 5 || (ESP_IDF_VERSION_MAJOR == 5 && ESP_IDF_VERSION_MINOR >= 5)
  // Codice per versioni 5.5 o superiori
  static void _esp_now_tx_cb(const esp_now_send_info_t *tx_info, esp_now_send_status_t status); // callback esp_now // poer nuova versione IDF 55
#endif
#endif // COMMENTO

void packet_callback_read_dongle(); // callback packet 
void packet_callback_read_gun(); // callback packet
void packet_callback_read_pedal(); // callback packet

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
  #elif defined(PEDAL)
  ///////////////////////uint8_t peerAddress[6] = {0x24, 0x58, 0x7C, 0xDA, 0x38, 0xA0}; // quello montato con piedini su breackboard
  uint8_t peerAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}; // broadcast

  const functionPtr callbackArr[] = { packet_callback_read_pedal };
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

// ========= per gestione PEDAL ===========
uint8_t peerAddress_pedal[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}; // broadcast
//unsigned long lastMillis_packet_pedal = 0;
#define MAX_TIMEOUT_LAST_PACKET (uint64_t)100000   // = 100  ms  .. il valore va specificato in microsecondi
void setupTimerPedal();
//void stopTimer_pedal();
//void resetTimer_pedal(uint64_t duration_us);
esp_timer_handle_t timer_handle_pedal;
//bool last_pedal;
//bool last_pedal2;
// ===================================================


///////////////////////////////////////////////////////////////////
// Definizione globale della configurazione del rate
esp_now_rate_config_t rate_config = {
    .phymode = WIFI_PHY_MODE_11G, // Forza lo standard G
    .rate = WIFI_PHY_RATE_12M,    // Forza i 12 Mbps
    .ersu = false,                // Non necessario per 11g
    .dcm = false                  // Non necessario per 11g
};
////////////////////////////////////////////////////////////////////

#ifdef OPENFIRE_ESPNOW_WIFI_POWER_AUTO
// ---------------- INIZIO LOGICA TPC ----------------
// Variabile globale per salvare l'RSSI istantaneo del pacchetto appena arrivato
volatile int8_t last_rssi_percepito = -100;
int8_t ultimo_rssi_trasmesso = TARGET_RSSI;
int8_t espnow_rssi_ricevuto = 0;
//volatile int8_t last_rssi_percepito = -100;
// ultimo_rssi_trasmesso = last_rssi_percepito;

uint8_t calcolaPotenzaOttimale(int8_t rssi_remoto) {
    // Usiamo direttamente int16_t per evitare cast impliciti continui nei calcoli successivi
    int16_t potenza_attuale = espnow_wifi_power;

    // Calcolo in un'unica riga: 
    // Usiamo lo shift bit a bit (<< 1) al posto della moltiplicazione (* 2) per massima velocità hardware
    potenza_attuale += (TARGET_RSSI - rssi_remoto) * 2; // * 4; 

    // AGGIUNGIAMO UN MARGINE DI SICUREZZA DI +2 dBm (8 unità) 
    // visto che la potenza non verrà più corretta in volo!
    //potenza_attuale += 8;
    
    #if defined(GUN)
        // La GUN scende per risparmiare batteria
        potenza_attuale = constrain(potenza_attuale, WIFI_POWER_MIN, WIFI_POWER_MAX);
    #elif defined(DONGLE)
        // Il DONGLE resta mediamente alto (minimo 40) perché ha la corrente via USB
        potenza_attuale = constrain(potenza_attuale, 40, WIFI_POWER_MAX);
    #elif defined(PEDAL)
        // Il DONGLE resta mediamente alto (minimo 40) perché ha la corrente via USB
        potenza_attuale = constrain(potenza_attuale, 40, WIFI_POWER_MAX);   
    #endif

    return (uint8_t)potenza_attuale;
}


// ---------------- FINE LOGICA TPC ----------------
#endif //OPENFIRE_ESPNOW_WIFI_POWER_AUTO

///////////////////////////////////////////////////////////////////

volatile uint8_t channel_display = espnow_wifi_channel;
volatile uint8_t seconds_display = 30;
volatile bool broadcast_receiver = false;  // ??????? non serve

//////////////////////////////////////////////////////////////////////////////
TaskHandle_t xUSBTaskHandle = NULL;

void usbTask(void *pvParameters) {
  // Teniamo traccia di chi ha il turno. 
  // 0 = Mouse, 1 = Tastiera, 2 = Gamepad
  static int currentRotation = 0; 

  for (;;) {
    // 1. METTI IN PAUSA IL TASK
    // Si sveglierà istantaneamente se:
    // - La radio riceve un dato nuovo (il tuo xTaskNotifyGive)
    // - La USB ha finito di trasmettere (la callback tud_hid_report_complete_cb)
    ulTaskNotifyTake(pdTRUE, portMAX_DELAY);

    // 2. CONTROLLO HARDWARE
    // Procediamo solo se l'Endpoint HID è fisicamente libero
    if (usbHid.ready()) {

      // 3. ROTAZIONE ROUND-ROBIN
      // Facciamo un ciclo di massimo 3 tentativi per trovare il primo dispositivo con dati pendenti
      for (int i = 0; i < 3; i++) {
        int target = (currentRotation + i) % 3;

        // --- TURNO 0: MOUSE ---
        if (target == 0 && absmouse5Report_pending) {
          if (usbHid.sendReport(HID_RID_e::HID_RID_MOUSE, &absmouse5Report_last_wifi, sizeof(absmouse5Report_last_wifi))) {
            absmouse5Report_pending = false; // Dato smaltito
            currentRotation = 1;            // Prossimo giro tocca alla Tastiera
            break; // IMPORTANTE: Usciamo dal for. Il buffer USB ora è pieno.
          }
        }
        
        // --- TURNO 1: TASTIERA ---
        else if (target == 1 && keyReport_pending) {
          if (usbHid.sendReport(HID_RID_e::HID_RID_KEYBOARD, &keyReport_last_wifi, sizeof(keyReport_last_wifi))) {
            keyReport_pending = false; 
            currentRotation = 2;            // Prossimo giro tocca al Gamepad
            break; 
          }
        }
        
        // --- TURNO 2: GAMEPAD ---
        else if (target == 2 && gamepad16Report_pending) {
          if (usbHid.sendReport(HID_RID_e::HID_RID_GAMEPAD, &gamepad16Report_last_wifi, sizeof(gamepad16Report_last_wifi))) {
            gamepad16Report_pending = false;
            currentRotation = 0;            // Prossimo giro si riparte dal Mouse
            break; 
          }
        }
        
      } // Fine del ciclo for
    }
    // Se arriviamo qui (che la USB fosse occupata o che abbiamo inviato un report),
    // il task finisce il giro e torna a dormire all'inizio del for(;;)
  }
}

/////////////////////////////////////////////////////////////////////////////
TaskHandle_t xRadioTaskHandle = NULL;
volatile bool radioFree = true;


void radioTask(void *pvParameters) {
  while (1) {
    ulTaskNotifyTake(pdTRUE, portMAX_DELAY);
    
    bool progress;
    
    do {
      progress = false;

      // 1. GESTIONE RICEZIONE (Invariata)
      if (SerialWireless._reader != SerialWireless._writer) {
        SerialWireless.checkForRxPacket();
        progress = true;
      }

      // 2. GESTIONE SERIALE -> RECOVERY MODE!
      // Interveniamo SOLO se il timer ha fallito in precedenza
      if (SerialWireless._serial_needs_recovery && (SerialWireless.availableBufferSerialWrite() > 0)) {
        if (SerialWireless.flush_sem()) {
          // Salvataggio riuscito! Abbassiamo la bandierina.
          // SerialWireless._serial_needs_recovery = false; // flush_sem() lo fa già da sola prima di restituire true.
          progress = true; 
        }
      }

      // 3. GESTIONE TRASMISSIONE RADIO (Invariata)
      if (radioFree && (SerialWireless.readIndex != SerialWireless.writeIndex)) {
        SerialWireless.SendData_sem();
        progress = true;
      }

    } while (progress); 
  } 
}
// =====================================


/////////////////////////////////////////////////////////////////////////////
// VERSIONE GRAFICA RICERCA CANALE PER GUN //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#if defined(GUN) && defined(USES_DISPLAY)
void animTaskLink_pedal(void *pvParameters) {
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
  uint8_t seconds_tread = seconds_display;
  
  // Setup iniziale display
  display_OLED->setCursor(baseX, baseY);
  display_OLED->setTextSize(1);
  display_OLED->setTextColor(WHITE, BLACK);
  display_OLED->fillRect(0, 0, 128, 16, BLACK);
  display_OLED->drawFastHLine(0, 15, 128, WHITE);
  char buffer[30];
  sprintf(buffer, "Search Pedal [%2ds]", seconds_tread);                  
  display_OLED->print(buffer);
  display_OLED->display();

  // Loop del task
  for (;;) {
    if (seconds_display != seconds_tread)
    {
      display_OLED->setCursor((14*6)+1, 2);
              //display_OLED->setTextSize(1);
              //display_OLED->setTextColor(WHITE, BLACK);
              ///////////display_OLED->fillRect(0, 0, 128, 16, BLACK);
              //display_OLED->drawFastHLine(0, 15, 128, WHITE);
      seconds_tread = seconds_display;
      sprintf(buffer, "%2d", seconds_tread);
      display_OLED->print(buffer);
              //display_OLED->display();          
      

      //seconds_tread = seconds_display;
    }
    display_OLED->setCursor(baseX + (19 * charWidth), baseY);
    display_OLED->write(rotazione[rotazioneIndex]);
    rotazioneIndex = (rotazioneIndex + 1) % 8;
    display_OLED->display();
    vTaskDelay(pdMS_TO_TICKS(150)); // piccolo delay per non saturare la CPU
  } 
}
#endif // GUN



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
  uint8_t channel_tread = espnow_wifi_channel;
  
  // Setup iniziale display
  display_OLED->setCursor(baseX, baseY);
  display_OLED->setTextSize(1);
  display_OLED->setTextColor(WHITE, BLACK);
  display_OLED->fillRect(0, 0, 128, 16, BLACK);
  display_OLED->drawFastHLine(0, 15, 128, WHITE);
  char buffer[30];
  sprintf(buffer, "Ch:%2d Waiting link", channel_tread);
  display_OLED->print(buffer);
  display_OLED->display();

  // Loop del task
  for (;;) {
    if (channel_display != channel_tread)
    {
      display_OLED->setCursor(19, 2);
              //display_OLED->setTextSize(1);
              //display_OLED->setTextColor(WHITE, BLACK);
              ///////////display_OLED->fillRect(0, 0, 128, 16, BLACK);
              //display_OLED->drawFastHLine(0, 15, 128, WHITE);
      channel_tread = channel_display;
      sprintf(buffer, "%2d", channel_tread);
      display_OLED->print(buffer);
              //display_OLED->display();          
      

      //channel_tread = channel_display;
    }
    display_OLED->setCursor(baseX + (19 * charWidth), baseY);
    display_OLED->write(rotazione[rotazioneIndex]);
    rotazioneIndex = (rotazioneIndex + 1) % 8;
    display_OLED->display();
    vTaskDelay(pdMS_TO_TICKS(150)); // piccolo delay per non saturare la CPU
  } 
}
#endif // GUN



/////////////////////////////////////////////////////////////////////////////
// FINE VERSIONE GRAFICA RICERCA CANALE PER GUN /////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// VERSIONE GRAFICA RICERCA CANALE PER DONGLE ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#if defined(DONGLE) && defined(USES_DISPLAY)
void animTaskLink(void *pvParameters) {
  unsigned long lastChange = 0;
  int8_t currentIndex = 0;
  int8_t direzione = 1;
  const char* word = "Scanning WiFi (((o)))"; // Waiting link 
  const char* rotazione = "-\\|/-\\|/";
  int8_t rotazioneIndex = 0;
  const uint8_t baseX = 75;
  const uint8_t baseY = 58;
  const uint8_t charWidth = 6;
  uint8_t len_word = strlen(word);

  
  // Setup iniziale display
  //tft.setCursor(baseX, baseY);
  tft.setTextSize(1);
  tft.setTextColor(WHITE, BLACK);
  tft.fillRect(70,20,100,60,BLACK);
  ////////tft.fillRect(0, 0, 128, 16, BLACK);
  ////////tft.drawFastHLine(0, 15, 128, WHITE);
  tft.setCursor(baseX, baseY - 30);
  tft.print("Selected");
  tft.setCursor(baseX, baseY - 19);
  tft.print("channel:");
  tft.setCursor(baseX + (9 * charWidth)+3, baseY - 26);
  tft.setTextSize(2);
  tft.setTextColor(RED, BLACK);
  char buffer[10];
  sprintf(buffer, "%2d", espnow_wifi_channel);
  tft.print(buffer);
  tft.setTextSize(1);
  tft.setTextColor(WHITE, BLACK);
  tft.setCursor(baseX, baseY);
  tft.print("Connecting");
  tft.display();
  tft.setTextSize(2);
  tft.setTextColor(BLUE, BLACK);

  // Loop del task
  for (;;) {    
    tft.setCursor(baseX + (11 * charWidth) + 3, baseY);
    tft.write(rotazione[rotazioneIndex]);
    rotazioneIndex = (rotazioneIndex + 1) % 8;
    tft.display();
    vTaskDelay(pdMS_TO_TICKS(150)); // piccolo delay per non saturare la CPU
  } 
}
#endif // DONGLE


//////////////////////////////////////////////////////////////////

#if defined(DONGLE) && defined(USES_DISPLAY)
void animTask(void *pvParameters) {
  int8_t currentIndex = 0;
  const uint8_t baseX = 75;
  const uint8_t baseY = 45;
  const uint8_t charWidth = 6;

  // Setup iniziale display
  tft.setCursor(baseX, baseY);
  tft.setTextSize(1);
  tft.setTextColor(WHITE, BLACK);
  //tft.fillRect(0, 0, 128, 16, BLACK);
  //tft.drawFastHLine(0, 15, 128, WHITE);
  //tft.setTextColor(RED, BLACK);
  tft.setCursor(baseX + 12, baseY - 20);
  tft.print("Searching");
  tft.setCursor(baseX+3, baseY - 10);
  tft.print("best channel");
  //tft.setTextColor(BLUE, BLACK);
  tft.setCursor(baseX, baseY);
  tft.print("Scanning WiFi");
  tft.display();
  tft.setTextColor(RED, BLACK);

  //tft.setTextColor(WHITE, BLACK);
  // Loop del task
  for (;;) {  
    tft.setCursor(baseX + (3 * charWidth), baseY + 18);
    switch (currentIndex)
    {
    case 0:
      /* code */
      tft.print("   o   ");
      break;
    case 1:
      /* code */
      tft.print("  (o)  ");
      break;
    case 2:
      /* code */
      tft.print(" ((o)) ");
      break;
    case 3:
      /* code */
      tft.print("(((o)))");
      break;
    
    default:
      
    break;
    }
    currentIndex = (currentIndex + 1) % 4;
    tft.display();
    vTaskDelay(pdMS_TO_TICKS(150)); // piccolo delay per non saturare la CPU
  }
}
#endif //DONGLE


/////////////////////////////////////////////////////////////////////////////
// FINE VERSIONE GRAFICA RICERCA CANALE PER DONGLE //////////////////////////
/////////////////////////////////////////////////////////////////////////////


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
  
  #if defined(DONGLE) && defined(USES_DISPLAY)
  TaskHandle_t animTaskHandle = NULL;  
  if(display_init) {
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
  }
  #endif // USES_DISPLAY


  // ================= STRUTTURA PER STATISTICHE CANALE PER RICERCA BEST CHANNEL=================
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
    if (ch == 1 || ch == 6 || ch == 11) {
      finalScores[ch] *= 0.92f; // -8% bonus (preferenza) 
    }  
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
  
  #if defined(DONGLE) && defined(USES_DISPLAY)
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
  if (_writerSerialRead == _readerSerialRead) return -1;
  return bufferSerialRead[_readerSerialRead];
}

int SerialWireless_::peekBin() {
  if (writeIndex == readIndex) return -1;
  return buffer[readIndex];
}

int SerialWireless_::read() {
  uint16_t r = _readerSerialRead;
  if (r == _writerSerialRead) return -1;

  uint8_t ret = bufferSerialRead[r];
  _readerSerialRead = (r + 1) & MASK_READ_SERIAL;
  return (int)ret;
}

int SerialWireless_::readBin() {
  uint16_t r = _reader;
  uint16_t w = _writer;

  if (r == w) return -1;
  uint8_t ret = _queue[r];

  const uint16_t MASK = FIFO_SIZE_READ - 1;
  _reader = (r + 1) & MASK;
  return (int)ret;
}

bool SerialWireless_::checkForRxPacket() {
    // 1. Snapshot locale degli indici in questo preciso istante
    uint16_t r = _reader;
    uint16_t w = _writer; // Guardiamo dove è arrivata la callback ESP-NOW
    const uint16_t MASK = FIFO_SIZE_READ - 1;

    // 2. Calcolo Lock-Free dei dati disponibili
    // Questa formula gestisce il wrap-around da sola in un ciclo di clock
    uint16_t available = (w - r) & MASK;
    
    // Uscita rapida se non c'è nulla da leggere
    if (available == 0) return false;

    // 3. Cache locale per massime prestazioni
    uint8_t* queuePtr = _queue;

    // 4. Ciclo di "svuotamento" rapido e Parsing
    for (uint16_t i = 0; i < available; i++) {
        uint8_t dato = queuePtr[r];
        
        // Il parser lavora in tempo reale.
        packet.parse(dato, true); 

        // Avanzamento con wrap-around bitwise
        r = (r + 1) & MASK;
    }

    // --- BARRIERA DEL COMPILATORE ---
    // Assicura che la CPU finisca fisicamente di leggere tutti i byte dalla RAM 
    // PRIMA di dire all'ISR che lo spazio è di nuovo libero. 
    asm volatile ("memw" : : : "memory");

    // 5. Sincronizzazione finale ATOMICA
    // Aggiorniamo _reader in RAM una sola volta. Appena viene eseguita questa riga,
    // la callback RX capisce matematicamente che c'è nuovo spazio libero.
    _reader = r;

    return true;
}

int SerialWireless_::available() {
  // Distanza tra chi scrive e chi legge
  return (_writerSerialRead - _readerSerialRead) & MASK_READ_SERIAL;
}

int SerialWireless_::availableBin() {
  return (writeIndex - readIndex) & (BUFFER_SIZE - 1);
}

// Quanti byte LIBERI rimangono nel buffer (Standard Arduino API)
int SerialWireless_::availableForWrite() {
  return (_readerSerialWrite - _writerSerialWrite - 1) & MASK_WRITE_SERIAL;
}

// funzione di comodo
// Quanti byte sono attualmente nel buffer in attesa di essere inviati
int SerialWireless_::availableBufferSerialWrite() {
  return (_writerSerialWrite - _readerSerialWrite) & MASK_WRITE_SERIAL;
}

int SerialWireless_::availableForWriteBin() {
  return (readIndex - writeIndex - 1) & (BUFFER_SIZE - 1);
}

void SerialWireless_::flush() {
  // 1. Il flush deve attendere finché ci sono dati nel buffer di scrittura seriale
  while (availableBufferSerialWrite() > 0) {
    
    // 2. Chiamiamo flush_sem(). 
    // Non serve prendere il Mutex qui perché flush_sem() al suo interno 
    // usa già lo Spinlock (portENTER_CRITICAL) per proteggere i puntatori.
    if (!flush_sem()) {
      // Se flush_sem() restituisce false, significa che il buffer RADIO è pieno.
      // Dobbiamo dare tempo al radioTask di trasmettere per liberare spazio.
      taskYIELD(); 
    }
    
    // Un piccolo respiro per non "sequestrare" il Core se il buffer è molto grande
    taskYIELD(); 
  }
}

bool SerialWireless_::flush_sem() {
  bool packed_data = false;
  
  // 1. Lucchetto Seriale
  portENTER_CRITICAL(&mux_serial_tx);

  uint16_t w = _writerSerialWrite;
  uint16_t r = _readerSerialWrite;
  uint16_t available_tx = (w - r) & MASK_WRITE_SERIAL;

  if (available_tx > 0) {
    // 2. LUCCHETTO RADIO PRIMA DI CONTARE LO SPAZIO!
    // Ora nessuno può rubarci lo spazio mentre calcoliamo il pacchetto
    portENTER_CRITICAL(&mux_radio_tx);

    uint16_t space_in_radio = availableForWriteBin();
    const uint16_t overhead = PREAMBLE_SIZE + POSTAMBLE_SIZE;

    if (space_in_radio > overhead) {
      uint16_t max_payload = space_in_radio - overhead;
      uint16_t len_to_pack = std::min((uint16_t)available_tx, std::min((uint16_t)200, max_payload));

      uint16_t firstChunk = FIFO_SIZE_WRITE_SERIAL - r;
      if (firstChunk < len_to_pack) {
        memcpy(&packet.txBuff[PREAMBLE_SIZE], &bufferSerialWrite[r], firstChunk);
        memcpy(&packet.txBuff[PREAMBLE_SIZE + firstChunk], &bufferSerialWrite[0], len_to_pack - firstChunk);
      } else {
        memcpy(&packet.txBuff[PREAMBLE_SIZE], &bufferSerialWrite[r], len_to_pack);
      }

      packet.constructPacket(len_to_pack, PACKET_TX::SERIAL_TX);
      
      writeBin(packet.txBuff, len_to_pack + overhead);
      
      _readerSerialWrite = (r + len_to_pack) & MASK_WRITE_SERIAL;
      packed_data = true;
    }
    
    // Rilasciamo la radio
    portEXIT_CRITICAL(&mux_radio_tx);
  }

  // Rilasciamo la seriale
  portEXIT_CRITICAL(&mux_serial_tx); 

  // API di sistema
  if (packed_data) {
    esp_timer_stop(timer_handle_serial);
    xTaskNotifyGive(xRadioTaskHandle);
    _serial_needs_recovery = false;
    return true;
  }

  return false; 
}

void SerialWireless_::SendData_sem() {
  // 1. Snapshot locale degli indici
  const uint16_t r = readIndex;
  const uint16_t w = writeIndex; // Guardiamo dove è arrivato il produttore in questo istante
  const uint16_t MASK = BUFFER_SIZE - 1;

  // 2. Calcolo Lock-Free dei dati disponibili
  // Questa formula gestisce il wrap-around automaticamente
  const uint16_t available_data = (w - r) & MASK;

  // 3. Uscita rapida se non c'è nulla da spedire
  if (available_data == 0) return;

  uint16_t len_tx = available_data;

  // 4. ALLINEAMENTO PACCHETTI (Prevenzione "taglio")
  if (available_data > ESP_NOW_MAX_DATA_LEN) {
    uint16_t aligned_len = 0;
    
    // Scansione a ritroso dal limite fisico.
    for (uint16_t i = ESP_NOW_MAX_DATA_LEN; i > 0; i--) {
      if (buffer[(r + i) & MASK] == START_BYTE) {
        aligned_len = i; 
        break;
      }
    }

    len_tx = (aligned_len > 0) ? aligned_len : ESP_NOW_MAX_DATA_LEN;
  }

  // 5. COPIA BATCH
  const uint16_t firstChunk = BUFFER_SIZE - r;
  
  if (firstChunk < len_tx) {
    const uint16_t secondChunk = len_tx - firstChunk;
    memcpy(buffer_espnow, &buffer[r], firstChunk);
    memcpy(&buffer_espnow[firstChunk], &buffer[0], secondChunk);
  } else {
    memcpy(buffer_espnow, &buffer[r], len_tx);
  }

// 6. INVIO E AGGIORNAMENTO
  esp_err_t result = esp_now_send(peerAddress, buffer_espnow, len_tx);
  
  if (result == ESP_OK) {
    radioFree = false; 
    
    // Aggiornamento ATOMICO
    readIndex = (r + len_tx) & MASK; 
  } else {
    // IL FREEZE SANO:
    // L'invio è fallito (es. peer spento/interferenza).
    // Non avanziamo readIndex per non perdere i dati seriali vitali.
    // Ma diamo 1 ms di pausa al sistema per non far intervenire il Watchdog
    // e prevenire un riavvio forzato della Lightgun.
    vTaskDelay(pdMS_TO_TICKS(1));
  }
}

size_t SerialWireless_::write(uint8_t c) {
  return write(&c, 1);
}

size_t SerialWireless_::writeBin(uint8_t c) {
  return writeBin(&c, 1);
}

size_t SerialWireless_::write(const uint8_t *data, size_t len) {
  size_t written = 0;
  bool start_timer = false; // Flag di sicurezza

  while (written < len) {
    start_timer = false; // Resettiamo a ogni ciclo
    
    portENTER_CRITICAL(&mux_serial_tx); // --- INIZIO LOCK ---

    uint16_t w = _writerSerialWrite;
    uint16_t r = _readerSerialWrite;
    uint16_t free_space = (r - w - 1) & MASK_WRITE_SERIAL;

    if (free_space > 0) {
      size_t toWrite = std::min((size_t)free_space, len - written);
      uint16_t firstChunk = FIFO_SIZE_WRITE_SERIAL - w;
      
      if (firstChunk < toWrite) {
        memcpy(&bufferSerialWrite[w], &data[written], firstChunk);
        memcpy(&bufferSerialWrite[0], &data[written + firstChunk], toWrite - firstChunk);
      } else {
        memcpy(&bufferSerialWrite[w], &data[written], toWrite);
      }

      if (w == r) {
        start_timer = true; // Alziamo la flag invece di chiamare l'API
      }

      asm volatile ("memw" : : : "memory");
      _writerSerialWrite = (w + toWrite) & MASK_WRITE_SERIAL;
      
      portEXIT_CRITICAL(&mux_serial_tx); // --- FINE LOCK ---

      // Ora siamo al sicuro: avviamo il timer se necessario!
      if (start_timer) {
        esp_timer_start_once(timer_handle_serial, TIMER_HANDLE_SERIAL_DURATION_MICROS);
      }

      written += toWrite;
      xTaskNotifyGive(xRadioTaskHandle);
    } 
    else {
      portEXIT_CRITICAL(&mux_serial_tx);
      
      if (!flush_sem()) {
        taskYIELD(); 
      }
    }
  }
  return written;
}

size_t SerialWireless_::writeBin(const uint8_t *data, size_t len) {
  // 1. Usa lo stesso lucchetto di SendPacket per coerenza e velocità
  portENTER_CRITICAL(&mux_radio_tx);

  uint16_t w = writeIndex;
  uint16_t r = readIndex;
  const uint16_t MASK = BUFFER_SIZE - 1;

  // 2. Calcolo dello spazio (Regola del Meno Uno) - PERFETTO
  uint16_t free_space = (r - w - 1) & MASK;

  if (free_space < len) {
    _overflow_write = true;
    portEXIT_CRITICAL(&mux_radio_tx); 
    return 0; 
  }

  // 3. COPIA BATCH (Wrap-around logic) - PERFETTA
  uint16_t firstChunk = BUFFER_SIZE - w;

  if (firstChunk < len) {
    memcpy(&buffer[w], data, firstChunk);
    memcpy(&buffer[0], data + firstChunk, len - firstChunk);
  } else {
    memcpy(&buffer[w], data, len);
  }

  // 4. BARRIERA DI MEMORIA - Fondamentale per ESP32/Xtensa
  asm volatile ("memw" : : : "memory");

  // 5. Aggiornamento dell'indice - Atomico
  writeIndex = (w + len) & MASK; 

  portEXIT_CRITICAL(&mux_radio_tx);

  return len;
}

void SerialWireless_::SendPacket(const uint8_t *data, uint8_t len, uint8_t packetID) { 
  const uint16_t totalLen = len + PREAMBLE_SIZE + POSTAMBLE_SIZE;

  // 1. Entriamo in sezione critica (Spinlock hardware)
  // Velocissimo, protegge il "tavolo da lavoro" packet.txBuff da altri Core o Task
  portENTER_CRITICAL(&mux_radio_tx);

  // 2. Snapshot Lock-Free dello spazio disponibile
  uint16_t w = writeIndex;
  uint16_t r = readIndex;
  const uint16_t MASK = BUFFER_SIZE - 1;
  uint16_t free_space = (r - w - 1) & MASK; 

  // 3. Se c'è spazio, costruiamo e scriviamo subito
  if (free_space >= totalLen) {
    uint8_t* pDest = packet.txBuff;

    // Copiamo i dati e costruiamo il pacchetto nel buffer temporaneo
    memcpy(&pDest[PREAMBLE_SIZE], data, len);
    packet.constructPacket(len, packetID);
    
    // Inseriamo nel Ring Buffer principale
    // (Dato che siamo già in sezione critica, writeBin sarà un fulmine)
    writeBin(pDest, totalLen);
    
    // Usciamo dalla sezione critica prima di notificare
    portEXIT_CRITICAL(&mux_radio_tx);

    // 4. Svegliamo il radioTask
    xTaskNotifyGive(xRadioTaskHandle);
  } 
  else {
    portEXIT_CRITICAL(&mux_radio_tx);
    _overflow_write = true;
  }
}

void SerialWireless_::write_on_rx_serialBuffer(const uint8_t *data, int len) {
  uint16_t w = _writerSerialRead;
  uint16_t r = _readerSerialRead;
  
  // Calcolo spazio libero (Regola del Meno Uno)
  uint16_t free_space = (r - w - 1) & MASK_READ_SERIAL;

  if (free_space >= len) {
    uint16_t firstChunk = FIFO_SIZE_READ_SERIAL - w;
    if (firstChunk < len) {
      memcpy(&bufferSerialRead[w], data, firstChunk);
      memcpy(&bufferSerialRead[0], data + firstChunk, len - firstChunk);
    } else {
      memcpy(&bufferSerialRead[w], data, len);
    }
    
    // Barriera di memoria prima di aggiornare il puntatore
    asm volatile ("memw" : : : "memory");
    _writerSerialRead = (w + len) & MASK_READ_SERIAL;
  } else {
    _overflow_bufferSerialRead = true;
  }
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

  setupTimerSerial(); // crea i timer .. timer per invio dati seriali

  // ========= AGGIUNTA PER NUOVA GESTIONE TRASMISSIONE RICEZIONE RADIO X TUTTI===============
  xTaskCreatePinnedToCore(
        radioTask,          // funzione del task
        "radioTask",        // nome
        4096,              // stack size
        NULL,              // parametri
        15,                 // priorità ALTA DA 1 A 24 ?
        &xRadioTaskHandle,   // handle
        APP_CPU_NUM        // core (puoi usare 0 o 1)
      );
  // ========== FINE AGGIUNTA PER NUOVA GESTIONE =================

  // ========= AGGIUNTA PER NUOVA GESTIONE DONGLE USB (SERVE SOLO AL DONGLE) ===============
  xTaskCreatePinnedToCore(
        usbTask,          // funzione del task
        "usbTask",        // nome
        4096,              // stack size
        NULL,              // parametri
        10,                 // priorità ALTA DA 1 A 24 ?
        &xUSBTaskHandle,   // handle
        APP_CPU_NUM        // core (puoi usare 0 o 1) APP_CPU_NUM = 1 (dove gira il loop) PRO_CPU_NUM = 0 (dove gira freertos, wifi, ecc.)
      );
  // ========== FINE AGGIUNTA PER NUOVA GESTIONE =================

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
      else espnow_wifi_channel=OPENFIRE_ESPNOW_WIFI_CHANNEL;
  #endif //GUN

  #ifdef DONGLE
    
      #ifdef OPENFIRE_AUTO_CHANNEL_ESPNOW_WIFI
        espnow_wifi_channel = findBestChannel(); //12;             
      #else
        espnow_wifi_channel=OPENFIRE_ESPNOW_WIFI_CHANNEL;
      #endif // OPENFIRE_AUTO_CHANNEL_ESPNOW_WIFI
     
  #endif // DONGLE

  #ifdef PEDAL
      espnow_wifi_channel=OPENFIRE_ESPNOW_WIFI_CHANNEL;
  #endif //GUN


  WiFi.mode(WIFI_STA); 
  WiFi.disconnect();  // ??? SI va messo
  

  esp_err_t err = esp_wifi_set_protocol(WIFI_IF_STA, WIFI_PROTOCOL_11G); // | WIFI_PROTOCOL_11N);
  if (err != ESP_OK) {
    //Serial.printf("esp_wifi_set_protocol failed! 0x%x", err);
  }
  
  err = esp_wifi_set_bandwidth(WIFI_IF_STA, WIFI_BW_HT20);
  if (err != ESP_OK) {
    //Serial.printf("esp_wifi_set_bandwidth failed! 0x%x", err);
  }
  
  err = esp_wifi_get_mac(WIFI_IF_STA, mac_esp_inteface);
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

  err = esp_wifi_set_ps(WIFI_PS_NONE);  // DISABILITA IL POWER SAVE PER MASSIMA REATTIVITà
  if (err != ESP_OK) {
    //Serial.printf("esp_wifi_set_ps failed! 0x%x", err);
  }
  
  vTaskDelay(pdMS_TO_TICKS(500)); // delay(1000);
    
  err = esp_now_init();
  if (err != ESP_OK) {
    //Serial.printf("esp_now_init failed! 0x%x", err);
  }
    
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
  } else esp_now_set_peer_rate_config(peerAddress, &rate_config);

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
  setupTimerSerial(); // crea i timer .. timer per invio dati seriali
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

// ============ NUOVA IMPLEMNTAZIONE == LA GUN FA IL FARO ================================
bool SerialWireless_::connection_dongle() {
    
  #if defined(DONGLE) && defined(USES_DISPLAY)
  TaskHandle_t animTaskHandleLink = NULL;  
  if(display_init) {
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
  }
  #endif // USES_DISPLAY

  #define TIMEOUT_DONGLE_DIALOGUE 2000 //  1000 // in millisecondi   ?????????????????????????
  unsigned long lastMillis_start_dialogue = millis ();
  
  
  
  lastMillis_start_dialogue = millis ();
  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
    
  // ====================================================

  broadcast_receiver = true;

  while(/*!TinyUSBDevice.mounted() && */stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) { 
    if (stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
      lastMillis_start_dialogue = millis ();
    }
    if (((millis() - lastMillis_start_dialogue) > TIMEOUT_DONGLE_DIALOGUE) && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
      stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
    }
    taskYIELD();
    //yield();
    //vTaskDelay(pdMS_TO_TICKS(10));
  }
  
  broadcast_receiver = false;

  //Serial.println("DONGLE - Negosazione completata - associazione dei dispositivi GUN/DONGLE");
  if (esp_now_del_peer(peerAddress) != ESP_OK) {  // cancella il broadcast dai peer
    //Serial.println("DONGLE - Errore nella cancellazione del peer broadcast");
  }
  memcpy(peerAddress, mac_esp_another_card, 6);
  memcpy(peerInfo.peer_addr, peerAddress, 6);
  ////////////espnow_wifi_channel=usb_data_wireless.channel;
  ////////////esp_wifi_set_promiscuous(true);
  ////////////esp_wifi_set_channel(espnow_wifi_channel, WIFI_SECOND_CHAN_NONE);
  ////////////esp_wifi_set_promiscuous(false);
  //peerInfo.channel = espnow_wifi_channel;
  //peerInfo.channel = ESPNOW_WIFI_CHANNEL;
  //peerInfo.encrypt = false;              
  if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserisce il dongle nei peer
    //Serial.println("DONGLE - Errore nell'aggiunta del nuovo peer della GUN");
  } else esp_now_set_peer_rate_config(peerAddress, &rate_config);

  #ifdef OPENFIRE_ESPNOW_WIFI_POWER_AUTO // ==========================================================
    if (espnow_wifi_power_auto) {
      SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
      SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
      SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
      //espnow_wifi_power = 45;
    }
  #endif //OPENFIRE_ESPNOW_WIFI_POWER_AUTO // ==========================================================


  TinyUSBDevices.onBattery = true;
  
    #if defined(DONGLE) && defined(USES_DISPLAY)
      if (animTaskHandleLink != NULL) {
        vTaskDelete(animTaskHandleLink);
        animTaskHandleLink = NULL;
      }
    #endif // USES_DISPLAY
  
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
      
      taskYIELD();
      //vTaskDelay(pdMS_TO_TICKS(10));
      //yield();
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

// ======================= NUOVA IMPLEMENTAZIONE DOVE LA GUN FA IL FARO ===============
bool SerialWireless_::connection_gun() {
  
  channel_display = espnow_wifi_channel;

  #if defined(GUN) && defined(USES_DISPLAY)
  TaskHandle_t animTaskHandleLink = NULL;  
  if(display_OLED != nullptr) {
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
  }
  #endif // USES_DISPLAY
  
  
  uint8_t channel = espnow_wifi_channel;   // tra e 1 e 13 (il 14 dovrebbe essere riservato)
  #define TIMEOUT_GUN_TX_PACKET 80  // 500 // in millisecondi
  #define TIMEOUT_GUN_CHANGE_CHANNEL 250  // 2000 // in millisecondi - cambia canale ogni
  #define TIMEOUT_GUN_DIALOGUE 3000  //  6000 // in millisecondi - tempo massimo per completare operazione accoppiamento
  unsigned long lastMillis_tx_packet = millis ();
  unsigned long lastMillis_change_channel = millis ();
  unsigned long lastMillis_start_dialogue = millis ();
  uint8_t aux_buffer_tx[14]; // aggiunto un byte per trasmettere anche il canale di trasmissione
                             // durante tx pacchetto pubblicazione presenza, mette anche il canale
  
  
  

  // ===== IMPOSTARE IL PEER BOARCAST QUI ====================== SERVE SE SI RIPARTE DA CERCA ULTIMO DONGLE, o forse no, non ricordo .... ????????? verificare  ????
  if (esp_now_del_peer(peerAddress) != ESP_OK) {  // cancella il broadcast dai peer
    //Serial.println("Errore nella cancellazione del peer");
  }
  memcpy(peerAddress, BROADCAST_ADDR, 6);
  memcpy(peerInfo.peer_addr, peerAddress, 6);
  //peerInfo.channel = ESPNOW_WIFI_CHANNEL;
  //peerInfo.encrypt = false;              
  if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserisce il dongle nei peer
    //Serial.println("Errore nell'aggiunta del nuovo peer");
  } else esp_now_set_peer_rate_config(peerAddress, &rate_config);                       
  // ====================================================
  
  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
  aux_buffer_tx[0] = CONNECTION_STATE::TX_GUN_SEARCH_DONGLE_BROADCAST;
  memcpy(&aux_buffer_tx[1], SerialWireless.mac_esp_inteface, 6);
  memcpy(&aux_buffer_tx[7], peerAddress, 6);
  //aux_buffer_tx[13] = espnow_wifi_channel;
 
  broadcast_receiver = true;

  while (!TinyUSBDevice.mounted() && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
    if (stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
      if (((millis() - lastMillis_change_channel) > TIMEOUT_GUN_CHANGE_CHANNEL) && 
         ((millis() - lastMillis_tx_packet)) >= (TIMEOUT_GUN_TX_PACKET - 50)) {  // aggiunta impostato 50 ms come margine, per evitare che quando invia pacchetto cambi subito casnale senza dare possibilità risposta
        //channel++;
        //if (channel >13) channel = 1;
        aux_buffer_tx[13] = channel;
        channel_display = channel;      
        
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
        lastMillis_tx_packet = 0;
        
        channel++;
        if (channel >13) channel = 1;// per fare inviare subito un pacchetto sul nuovo canale
      }
      if ((millis() - lastMillis_tx_packet) > TIMEOUT_GUN_TX_PACKET) {
        SerialWireless.SendPacket((const uint8_t *)aux_buffer_tx, 14, PACKET_TX::CONNECTION); // aggiunto un byte per trasmettere anche il canale di trasmissione
        //Serial.print("DONGLE - inviato pacchetto broadcast sul canale: ");
        //Serial.println(channel);
        lastMillis_tx_packet = millis (); 
      }
      lastMillis_start_dialogue = millis();
    }
    else {
      if (((millis() - lastMillis_start_dialogue) > TIMEOUT_GUN_DIALOGUE) && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
        stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
        //Serial.println("DONGLE - Non si è conclusa la negoziazione tra DONGLE/GUN e si riparte da capo");
        lastMillis_change_channel = millis ();
      }  
    }
    taskYIELD();
    //yield(); // in attesa dello stabilimento di una connessione
    //vTaskDelay(pdMS_TO_TICKS(10));
  }

  broadcast_receiver = false;  

  if (stato_connessione_wireless == CONNECTION_STATE::DEVICES_CONNECTED) {
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
    } else esp_now_set_peer_rate_config(peerAddress, &rate_config);

    #ifdef OPENFIRE_ESPNOW_WIFI_POWER_AUTO // ==========================================================
      if (espnow_wifi_power_auto) {
        SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
        SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
        SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
        //espnow_wifi_power = 45;
      }
    #endif //OPENFIRE_ESPNOW_WIFI_POWER_AUTO // ==========================================================


    #if defined(GUN) && defined(USES_DISPLAY)
      if (animTaskHandleLink != NULL) {
        vTaskDelete(animTaskHandleLink);
        animTaskHandleLink = NULL;
      }
    #endif // USES_DISPLAY


    TinyUSBDevices.onBattery = true;
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

// ======================= NUOVA IMPLEMENTAZIONE DOVE LA GUN FA IL FARO per connetersi al pedal ===============
bool SerialWireless_::connection_gun_at_pedal() {
  
  uint8_t seconds = 10; 
  seconds_display = seconds;

  channel_display = espnow_wifi_channel;

  #if defined(GUN) && defined(USES_DISPLAY)
  TaskHandle_t animTaskHandleLink_pedal = NULL;  
  if(display_OLED != nullptr) {
  // Avvio animazione
    if (animTaskHandleLink_pedal == NULL) {
      xTaskCreatePinnedToCore(
        animTaskLink_pedal,          // funzione del task
        "AnimTaskLink_pedal",        // nome
        4096,              // stack size
        NULL,              // parametri
        1,                 // priorità
        &animTaskHandleLink_pedal,   // handle
        APP_CPU_NUM        // core (puoi usare 0 o 1)
      );
    }
  }
  #endif // USES_DISPLAY
  
  
  //uint8_t seconds = 10;  
  #define TIMEOUT_GUN_AT_PEDAL_TX_PACKET 50 //500 // in millisecondi
  #define TIMEOUT_GUN_AT_PEDAL_CHANGE_SECONDS 1000 // in millisecondi - cambia canale ogni
  #define TIMEOUT_GUN_AT_PEDAL_DIALOGUE 3000 //6000 // in millisecondi - tempo massimo per completare operazione accoppiamento
  unsigned long lastMillis_tx_packet = millis ();
  unsigned long lastMillis_change_seconds = millis ();
  unsigned long lastMillis_start_dialogue = millis ();
  uint8_t aux_buffer_tx[14]; // aggiunto un byte per trasmettere anche il canale di trasmissione
                             // durante tx pacchetto pubblicazione presenza, mette anche il canale
  uint8_t peerAddress_copy[6];
  
  
  
  ////////// ===== IMPOSTARE IL PEER BOARCAST QUI ====================== SERVE SE SI RIPARTE DA CERCA ULTIMO DONGLE, o forse no, non ricordo .... ????????? verificare  ????
  //////////if (esp_now_del_peer(peerAddress) != ESP_OK) {  // cancella il broadcast dai peer
  //////////  //Serial.println("Errore nella cancellazione del peer");
  //////////}
  memcpy(peerAddress_copy, peerAddress, 6);

  memcpy(peerAddress, BROADCAST_ADDR, 6);        
  memcpy(peerInfo.peer_addr, peerAddress, 6);
  peerInfo.channel = espnow_wifi_channel;
  //peerInfo.encrypt = false;              
  if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserisce il dongle nei peer
    //Serial.println("Errore nell'aggiunta del nuovo peer");
  } else esp_now_set_peer_rate_config(peerAddress, &rate_config);                       
  // ====================================================
  
  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
  aux_buffer_tx[0] = CONNECTION_STATE::TX_GUN_SEARCH_PEDAL_BROADCAST;
  memcpy(&aux_buffer_tx[1], SerialWireless.mac_esp_inteface, 6);
  memcpy(&aux_buffer_tx[7], peerAddress, 6);
  aux_buffer_tx[13] = espnow_wifi_channel;
 

  broadcast_receiver = true;
  
  while ( (seconds > 1) && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
    if (stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
      if (((millis() - lastMillis_change_seconds) > TIMEOUT_GUN_AT_PEDAL_CHANGE_SECONDS) && 
         ((millis() - lastMillis_tx_packet)) >= (TIMEOUT_GUN_AT_PEDAL_TX_PACKET - 50)) {  // aggiunta impostato 50 ms come margine, per evitare che quando invia pacchetto cambi subito casnale senza dare possibilità risposta
        //channel++;
        //if (channel >13) channel = 1;
        aux_buffer_tx[13] = espnow_wifi_channel;
        seconds_display = seconds;
                        
        lastMillis_change_seconds = millis ();
        lastMillis_tx_packet = 0;
        
        seconds--;
      }
      if ((millis() - lastMillis_tx_packet) > TIMEOUT_GUN_AT_PEDAL_TX_PACKET) {
        // ??????????????????????? siatemare per gestione pedal di SendPacket ????????????????? .. ci vorrà un buffer per ogni peer o gestire diversamente le cose a seconda di chi si trasmettono i dati oppure si sposta provvisoriamente il peeraddres solo durante ricerca pedale
        SerialWireless.SendPacket((const uint8_t *)aux_buffer_tx, 14, PACKET_TX::CONNECTION_PEDAL); // aggiunto un byte per trasmettere anche il canale di trasmissione
        //Serial.print("DONGLE - inviato pacchetto broadcast sul canale: ");
        //Serial.println(channel);
        lastMillis_tx_packet = millis (); 
      }
      lastMillis_start_dialogue = millis();
    }
    else {
      if (((millis() - lastMillis_start_dialogue) > TIMEOUT_GUN_AT_PEDAL_DIALOGUE) && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
        stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
        //Serial.println("DONGLE - Non si è conclusa la negoziazione tra DONGLE/GUN e si riparte da capo");
        lastMillis_change_seconds = millis ();
      }  
    }
    taskYIELD();
    //yield(); // in attesa dello stabilimento di una connessione
    //vTaskDelay(pdMS_TO_TICKS(10));
  }
  
  broadcast_receiver = false;

  if (esp_now_del_peer(peerAddress) != ESP_OK) {  // cancella il broadcast dai peer
      //Serial.println("DONGLE - Errore nella cancellazione del peer broadcast");
  }

  memcpy(peerAddress, peerAddress_copy, 6);

  for (uint8_t i=0; i<3; i++) {
    SerialWireless.SendPacket((const uint8_t *)&TinyUSBDevices.is_pedal_wireless, 1, PACKET_TX::PEDAL_TX);
    vTaskDelay(pdMS_TO_TICKS(69));
  }

  if (stato_connessione_wireless == CONNECTION_STATE::DEVICES_CONNECTED) {
    //Serial.println("DONGLE - Negosazione completata - associazione dei dispositivi GUN/DONGLE");
    /*
    if (esp_now_del_peer(peerAddress) != ESP_OK) {  // cancella il broadcast dai peer
      //Serial.println("DONGLE - Errore nella cancellazione del peer broadcast");
    }
    */

    memcpy(peerAddress_pedal, mac_esp_another_card, 6);
    memcpy(peerInfo.peer_addr, peerAddress_pedal, 6);
    /////espnow_wifi_channel=usb_data_wireless.channel;
    /////esp_wifi_set_promiscuous(true);
    /////esp_wifi_set_channel(espnow_wifi_channel, WIFI_SECOND_CHAN_NONE);
    /////esp_wifi_set_promiscuous(false);
    peerInfo.channel = espnow_wifi_channel;
    //peerInfo.channel = ESPNOW_WIFI_CHANNEL;
    //peerInfo.encrypt = false;              
    if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserisce il dongle nei peer
      //Serial.println("DONGLE - Errore nell'aggiunta del nuovo peer della GUN");
    } else esp_now_set_peer_rate_config(peerAddress_pedal, &rate_config);

    setupTimerPedal();  // CREA IL TIMER PER GESTIRE IN SICUREZZA IL PEDALE

    #ifdef OPENFIRE_ESPNOW_WIFI_POWER_AUTO // ==========================================================
      if (espnow_wifi_power_auto) {
        SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
        SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
        SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
        //espnow_wifi_power = 45;
      }
    #endif //OPENFIRE_ESPNOW_WIFI_POWER_AUTO // ==========================================================


    #if defined(GUN) && defined(USES_DISPLAY)
      if (animTaskHandleLink_pedal != NULL) {
        vTaskDelete(animTaskHandleLink_pedal);
        animTaskHandleLink_pedal = NULL;
      }
    #endif // USES_DISPLAY

    TinyUSBDevices.is_pedal_wireless = true;
    //TinyUSBDevices.onBattery = true;
    return true;  
  }
  
  #if defined(GUN) && defined(USES_DISPLAY)
    if (animTaskHandleLink_pedal != NULL) {
      vTaskDelete(animTaskHandleLink_pedal);
      animTaskHandleLink_pedal = NULL;
    }
  #endif // USES_DISPLAY

  TinyUSBDevices.is_pedal_wireless = false;
  return false; 

}


// ======================= IMPLEMENTAZIONE PEDAL - RIMANE IN ASCOLTO SU TUTTI I CANALI ===============
bool SerialWireless_::connection_pedal() {
  
    
  uint8_t channel = espnow_wifi_channel;   // tra e 1 e 13 (il 14 dovrebbe essere riservato)
  
  #define TIMEOUT_PEDAL_CHANGE_CHANNEL 150 // 2000 // in millisecondi - cambia canale ogni
  #define TIMEOUT_PEDAL_DIALOGUE 3000 // 6000 // in millisecondi - tempo massimo per completare operazione accoppiamento
  unsigned long lastMillis_change_channel = millis ();
  unsigned long lastMillis_start_dialogue = millis ();
  
   
  
  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
  
  broadcast_receiver = true;

  while (stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
    if (stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
      if (((millis() - lastMillis_change_channel) > TIMEOUT_PEDAL_CHANGE_CHANNEL)) { 
         
        esp_wifi_set_promiscuous(true);
        if (esp_wifi_set_channel(channel, WIFI_SECOND_CHAN_NONE) != ESP_OK) {
          //Serial.printf("DONGLE - esp_wifi_set_channel failed!");
        }
        esp_wifi_set_promiscuous(false);
        
        espnow_wifi_channel = channel;
        
        peerInfo.channel = channel;  // si potrebbe anche impostare sempre a 0 che usa quello del wifi
        if (esp_now_mod_peer(&peerInfo) != ESP_OK) {  // modifica il canale del peer
          //Serial.println("DONGLE - Errore nella modifica del canale");
        }             
        lastMillis_change_channel = millis ();
        
        
        channel++;
        if (channel >13) channel = 1;
      }
      
      lastMillis_start_dialogue = millis();
    }
    else {
      if (((millis() - lastMillis_start_dialogue) > TIMEOUT_PEDAL_DIALOGUE) && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
        stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
        //Serial.println("DONGLE - Non si è conclusa la negoziazione tra DONGLE/GUN e si riparte da capo");
        lastMillis_change_channel = millis ();
      }  
    }
    taskYIELD();
    //vTaskDelay(pdMS_TO_TICKS(10));
    //yield(); // in attesa dello stabilimento di una connessione
  }
    
  broadcast_receiver = false;
    
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
    } else esp_now_set_peer_rate_config(peerAddress, &rate_config);

    #ifdef OPENFIRE_ESPNOW_WIFI_POWER_AUTO // ==========================================================
      if (espnow_wifi_power_auto) {
        SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
        SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
        SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_FEEDBACK);
        //espnow_wifi_power = 45;
      }
    #endif //OPENFIRE_ESPNOW_WIFI_POWER_AUTO // ==========================================================
  
    TinyUSBDevices.onBattery = true;
    return true;  
      
}

void packet_callback_read_dongle() {
  switch (SerialWireless.packet.currentPacketID()) {
    case PACKET_TX::SERIAL_TX:
      Serial.write(&SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
      Serial.flush();
      break;
    case PACKET_TX::MOUSE_TX : {
      uint8_t* ptr = &SerialWireless.packet.rxBuff[PREAMBLE_SIZE];
      if (memcmp(ptr, &absmouse5Report_last_wifi, sizeof(absmouse5Report_last_wifi))) {
        memcpy(&absmouse5Report_last_wifi, ptr, sizeof(absmouse5Report_last_wifi));
        absmouse5Report_pending = true;    
        xTaskNotifyGive(xUSBTaskHandle);
      }
      //usbHid.sendReport(HID_RID_e::HID_RID_MOUSE, ptr, sizeof(absmouse5Report_last_wifi)); // per inviare a USB
      break;
    }
    case PACKET_TX::GAMEPADE_TX: {
      uint8_t* ptr = &SerialWireless.packet.rxBuff[PREAMBLE_SIZE];
      // --- GAMEPAD ---
      if (memcmp(ptr, &gamepad16Report_last_wifi, sizeof(gamepad16Report_last_wifi))) {
        memcpy(&gamepad16Report_last_wifi, ptr, sizeof(gamepad16Report_last_wifi));
        gamepad16Report_pending = true;
        xTaskNotifyGive(xUSBTaskHandle);
      }
      //usbHid.sendReport(HID_RID_e::HID_RID_GAMEPAD, ptr, sizeof(gamepad16Report_last_wifi)); // per inviare a USB
      break;
    }
    case PACKET_TX::KEYBOARD_TX: {
      uint8_t* ptr = &SerialWireless.packet.rxBuff[PREAMBLE_SIZE];
      if (memcmp(ptr, &keyReport_last_wifi, sizeof(keyReport_last_wifi))) {
        memcpy(&keyReport_last_wifi, ptr, sizeof(keyReport_last_wifi));
        keyReport_pending = true;
        xTaskNotifyGive(xUSBTaskHandle);
      }
      //usbHid.sendReport(HID_RID_e::HID_RID_KEYBOARD, ptr, sizeof(keyReport_last_wifi));  // per inviare a USB
      break;
    }
    // ====== poi va tolto ================
    case PACKET_TX::PEDAL_TX:
      {
        //memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
        //TinyUSBDevices.is_pedal_wireless = (bool)aux_buffer[0];
        SerialWireless.is_pedal_wireless_comunication = true;
      }
      break;
    // ====== fine poi va tolto  
    #ifdef OPENFIRE_USE_ESPNOW_UNIFIED_PACKET
    case PACKET_TX::MOUSE_KEY_PAD_TX: { // ATTENZIONE: Le parentesi graffe qui sono OBBLIGATORIE in C++ per dichiarare variabili dentro un 'case'
      
      // Calcoliamo l'indirizzo base UNA sola volta e lo mettiamo nel puntatore 'ptr'
      uint8_t* ptr = &SerialWireless.packet.rxBuff[PREAMBLE_SIZE];

      // --- MOUSE ---
      // Confrontiamo il buffer con l'ultimo stato noto
      if (memcmp(ptr, &absmouse5Report_last_wifi, sizeof(absmouse5Report_last_wifi))) {
        memcpy(&absmouse5Report_last_wifi, ptr, sizeof(absmouse5Report_last_wifi));
        absmouse5Report_pending = true;    
        //xTaskNotifyGive(xUSBTaskHandle);
      }   

      ptr += sizeof(absmouse5Report_last_wifi); // Fai scorrere il puntatore in avanti alla fine dei dati mouse

      // --- TASTIERA ---
      if (memcmp(ptr, &keyReport_last_wifi, sizeof(keyReport_last_wifi))) {
        memcpy(&keyReport_last_wifi, ptr, sizeof(keyReport_last_wifi));
        keyReport_pending = true;
        //xTaskNotifyGive(xUSBTaskHandle);
      }
      
      ptr += sizeof(keyReport_last_wifi); // Fai scorrere il puntatore in avanti alla fine dei dati tastiera

      // --- GAMEPAD ---
      if (memcmp(ptr, &gamepad16Report_last_wifi, sizeof(gamepad16Report_last_wifi))) {
        memcpy(&gamepad16Report_last_wifi, ptr, sizeof(gamepad16Report_last_wifi));
        gamepad16Report_pending = true;
        //xTaskNotifyGive(xUSBTaskHandle);
      }
      
      if (absmouse5Report_pending || keyReport_pending || gamepad16Report_pending) xTaskNotifyGive(xUSBTaskHandle);

      break;  
    } // Chiusura del blocco del case
    #endif // OPENFIRE_USE_ESPNOW_UNIFIED_PACKET
    #ifdef OPENFIRE_ESPNOW_WIFI_POWER_AUTO
    case PACKET_TX::DUMMY_PACKET:
      /* pacchetto vuoto - non fare nulla - per eventuali test - da implementare se necessario */
      /* code */  
      break; 
    case PACKET_TX::RSSI_FEEDBACK:
      /* richiesta di inviare rssi quindi inviare un pacchetto rssi_report con l'rssi */
      aux_buffer[0] = last_rssi_percepito;
      SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_REPORT);
      break; 
    case PACKET_TX::RSSI_REPORT:
      /* ricevuto ul pacchetto rssi, quindi adattare la potenza di trasmissione in base a tale valore */
      {
        memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
        int8_t rssi_misurato_dall_altro = (int8_t)aux_buffer[0];
        uint8_t mia_nuova_potenza = calcolaPotenzaOttimale(rssi_misurato_dall_altro);
        if (esp_wifi_set_max_tx_power(mia_nuova_potenza) == ESP_OK) {
          espnow_wifi_power = mia_nuova_potenza;
          espnow_rssi_ricevuto = rssi_misurato_dall_altro;
          ultimo_rssi_trasmesso = last_rssi_percepito;
          //espnow_wifi_power = 41;
        }
      }
      break; 
    #endif //OPENFIRE_ESPNOW_WIFI_POWER_AUTO
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
                  while ((millis() - lastMillis_tx_packet_connection_last_dongle) < 70) taskYIELD(); ///////yield(); 
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
        // =========================== NUOVI ==================================
        case CONNECTION_STATE::TX_GUN_SEARCH_DONGLE_BROADCAST:
          if ((SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) &&
              (aux_buffer[13] == espnow_wifi_channel))
          { // prende la prima gun disposnibile
            memcpy(SerialWireless.mac_esp_another_card, &aux_buffer[1], 6);
            // invia richiesta connessione
            aux_buffer[0] = CONNECTION_STATE::TX_DONGLE_TO_GUN_PRESENCE;
            memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
            memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
            aux_buffer[13] = espnow_wifi_channel;
            SerialWireless.SendPacket((const uint8_t *)aux_buffer, 14, PACKET_TX::CONNECTION);
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::TX_DONGLE_TO_GUN_PRESENCE;
            // assicurati che i dati siano stati spediti

          }
          break;
        case CONNECTION_STATE::TX_GUN_TO_DONGLE_ACCEPT:
          if ((memcmp(&aux_buffer[1],SerialWireless.mac_esp_another_card,6) == 0) && 
             (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
              SerialWireless.stato_connessione_wireless == CONNECTION_STATE::TX_DONGLE_TO_GUN_PRESENCE) {
              
                // SALVA I DATI RELATIVI ALLA GUN VID, PID, PLAYER , ECC.ECC.
              memcpy(&usb_data_wireless, &aux_buffer[13], sizeof(usb_data_wireless));
              
              aux_buffer[0] = CONNECTION_STATE::TX_DONGLE_TO_GUN_CONFERM;
              memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
              memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
                            
              // =========================================================
              // INVIARE IL PACCHETTO FINALE PIU' VOLTE  
              // =========================================================
              for (uint8_t i = 0; i<3; i++) {
                if (i>0) {
                  //vTaskDelay(pdMS_TO_TICKS(1000)); // equivalente a delay(1000) ma non bloccante su esp32
                  unsigned long lastMillis_tx_packet_gun_to_dongle_conferm = millis();
                  while ((millis() - lastMillis_tx_packet_gun_to_dongle_conferm) < 70) taskYIELD(); /////yield(); 
                }
                SerialWireless.SendPacket((const uint8_t *)aux_buffer, 13, PACKET_TX::CONNECTION);
              }
              SerialWireless.stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED;
          }
          break;          
        // =========================== FINE NUOVI =============================

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
    case PACKET_TX::PEDAL_TX:
      {
 
  // =====================================================
        
        //uint8_t pedali;
        //  pedal attivo = 00000001 - pedal2 attivo = 00000010 - entrambi pedali attivi = 00000011 - tutti i pedali non premuti = 00000000
        //pedali = aux_buffer[0];
        memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
        TinyUSBDevices.pedals_wireless_state = aux_buffer[0];
        //TinyUSBDevices.pedal_wireless = (pedali & 0b00000001) != 0;  //false = pedale non premuto;
        //TinyUSBDevices.pedal2_wireless = (pedali & 0b00000010) != 0; //false = pedale non premuto;
        esp_timer_stop(timer_handle_pedal); 
        if (TinyUSBDevices.pedals_wireless_state) esp_timer_start_once(timer_handle_pedal, MAX_TIMEOUT_LAST_PACKET);
      }
      break;
    #ifdef OPENFIRE_ESPNOW_WIFI_POWER_AUTO
    case PACKET_TX::DUMMY_PACKET:
      /* pacchetto vuoto - non fare nulla - per eventuali test - da implementare se necessario */
      /* code */  
      break; 
    case PACKET_TX::RSSI_FEEDBACK:
      /* richiesta di inviare rssi quindi inviare un pacchetto rssi_report con l'rssi */
      aux_buffer[0] = last_rssi_percepito;
      SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_REPORT);
      break; 
    case PACKET_TX::RSSI_REPORT:
      /* ricevuto ul pacchetto rssi, quindi adattare la potenza di trasmissione in base a tale valore */
      {
        memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
        int8_t rssi_misurato_dall_altro = (int8_t)aux_buffer[0];
        uint8_t mia_nuova_potenza = calcolaPotenzaOttimale(rssi_misurato_dall_altro);
        if (esp_wifi_set_max_tx_power(mia_nuova_potenza) == ESP_OK) {
          espnow_wifi_power = mia_nuova_potenza;
          espnow_rssi_ricevuto = rssi_misurato_dall_altro;
          ultimo_rssi_trasmesso = last_rssi_percepito;
          //espnow_wifi_power = 41;
        }
      }      
      break;
    #endif //OPENFIRE_ESPNOW_WIFI_POWER_AUTO
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
        // ================== NUOVI ==================== // AL PRIMO PACCHETTO IL DONGLE TRASMETTE ANCHE IL CANALE PER EVITARE CHE NEL FRATTEMPO IL CICLO VADA AVANTI ED IL CANALE 
        case CONNECTION_STATE::TX_DONGLE_TO_GUN_PRESENCE:
          if ((memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) && SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
            memcpy(SerialWireless.mac_esp_another_card, &aux_buffer[1], 6);
            usb_data_wireless.channel= aux_buffer[13];
            aux_buffer[0] = CONNECTION_STATE::TX_GUN_TO_DONGLE_ACCEPT;
            memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
            memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
            // INVIA ANCHE DATI RELATIVI A VID, PID, ECC,ECC, DELLA GUN
            memcpy(&aux_buffer[13], &usb_data_wireless, sizeof(usb_data_wireless));
            SerialWireless.SendPacket((const uint8_t *)aux_buffer, sizeof(aux_buffer), PACKET_TX::CONNECTION);
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::TX_GUN_TO_DONGLE_ACCEPT;
          }
          break;
        case CONNECTION_STATE::TX_DONGLE_TO_GUN_CONFERM:
          if ((memcmp(&aux_buffer[1],SerialWireless.mac_esp_another_card,6) == 0) && 
             (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
             SerialWireless.stato_connessione_wireless == CONNECTION_STATE::TX_GUN_TO_DONGLE_ACCEPT) {           
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED;
        }
          break;
        // ================== FINE NUOVI =================
        default:
          break;
      }
      break;
    case PACKET_TX::CONNECTION_PEDAL:
      memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead); //13); //sizeof(aux_buffer)); // qui va bene anche 13 come dati da copiare
      switch (aux_buffer[0]) {
        // ================== NUOVI ==================== // AL PRIMO PACCHETTO IL DONGLE TRASMETTE ANCHE IL CANALE PER EVITARE CHE NEL FRATTEMPO IL CICLO VADA AVANTI ED IL CANALE 
        case CONNECTION_STATE::TX_PEDAL_TO_GUN_PRESENCE:
          if ((memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) && SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
            memcpy(SerialWireless.mac_esp_another_card, &aux_buffer[1], 6);
            usb_data_wireless.channel= espnow_wifi_channel;  // non dovrebbe servire
            aux_buffer[0] = CONNECTION_STATE::TX_GUN_TO_PEDAL_ACCEPT;
            memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
            memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
            // INVIA ANCHE DATI RELATIVI A VID, PID, ECC,ECC, DELLA GUN
            memcpy(&aux_buffer[13], &usb_data_wireless, sizeof(usb_data_wireless));
            SerialWireless.SendPacket((const uint8_t *)aux_buffer, sizeof(aux_buffer), PACKET_TX::CONNECTION_PEDAL);
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::TX_GUN_TO_PEDAL_ACCEPT;
          }
          break;
        case CONNECTION_STATE::TX_PEDAL_TO_GUN_CONFERM:
          if ((memcmp(&aux_buffer[1],SerialWireless.mac_esp_another_card,6) == 0) && 
             (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
             SerialWireless.stato_connessione_wireless == CONNECTION_STATE::TX_GUN_TO_PEDAL_ACCEPT) {           
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED;
        }
          break;
        // ================== FINE NUOVI =================
        default:
          break;
      }
      break;
    default:
      break;
  }
}

void packet_callback_read_pedal() {
  switch (SerialWireless.packet.currentPacketID()) {
    case PACKET_TX::SERIAL_TX:
      break;
    case PACKET_TX::MOUSE_TX :
      break;
    case PACKET_TX::GAMEPADE_TX:
      break;
    case PACKET_TX::KEYBOARD_TX:
      break;
    #ifdef OPENFIRE_USE_ESPNOW_UNIFIED_PACKET
    case PACKET_TX::MOUSE_KEY_PAD_TX: { // ATTENZIONE: Le parentesi graffe qui sono OBBLIGATORIE in C++ per dichiarare variabili dentro un 'case'
      break;  
    } // Chiusura del blocco del case
    #endif // OPENFIRE_USE_ESPNOW_UNIFIED_PACKET
    #ifdef OPENFIRE_ESPNOW_WIFI_POWER_AUTO
    case PACKET_TX::DUMMY_PACKET:
      /* pacchetto vuoto - non fare nulla - per eventuali test - da implementare se necessario */
      /* code */  
      break; 
    case PACKET_TX::RSSI_FEEDBACK:
      /* richiesta di inviare rssi quindi inviare un pacchetto rssi_report con l'rssi */
      aux_buffer[0] = last_rssi_percepito;
      SerialWireless.SendPacket((const uint8_t *)aux_buffer, 1, PACKET_TX::RSSI_REPORT);
      break; 
    case PACKET_TX::RSSI_REPORT:
      /* ricevuto ul pacchetto rssi, quindi adattare la potenza di trasmissione in base a tale valore */
      {
        memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
        int8_t rssi_misurato_dall_altro = (int8_t)aux_buffer[0];
        uint8_t mia_nuova_potenza = calcolaPotenzaOttimale(rssi_misurato_dall_altro);
        if (esp_wifi_set_max_tx_power(mia_nuova_potenza) == ESP_OK) {
          espnow_wifi_power = mia_nuova_potenza;
          espnow_rssi_ricevuto = rssi_misurato_dall_altro;
          ultimo_rssi_trasmesso = last_rssi_percepito;
          //espnow_wifi_power = 41;
        }
      }
      break; 
    #endif //OPENFIRE_ESPNOW_WIFI_POWER_AUTO
    case PACKET_TX::CHECK_CONNECTION_LAST_PEDAL:
      // CODICE
      memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
      switch (aux_buffer[0])
      {
      case CONNECTION_STATE::TX_CHECK_CONNECTION_LAST_PEDAL:
        if ((memcmp(&aux_buffer[1],peerAddress,6) == 0) && 
           (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
           SerialWireless.stato_connessione_wireless == CONNECTION_STATE::DEVICES_CONNECTED) { 
              aux_buffer[0] = CONNECTION_STATE::TX_CONFERM_CONNECTION_LAST_PEDAL; 
              memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
              memcpy(&aux_buffer[7], peerAddress, 6);
              // valutare se inviarlo un paio di volte il pacchetto o solo una volta
              // lo invia 3 volte - una volta ogni 70ms
              for (uint8_t i = 0; i<3; i++) {
                if (i>0) {
                  //vTaskDelay(pdMS_TO_TICKS(1000)); // equivalente a delay(1000) ma non bloccante su esp32
                  unsigned long lastMillis_tx_packet_connection_last_dongle = millis();
                  while ((millis() - lastMillis_tx_packet_connection_last_dongle) < 70) taskYIELD();  ////yield(); 
                }
                SerialWireless.SendPacket((const uint8_t *)aux_buffer, 13, PACKET_TX::CHECK_CONNECTION_LAST_PEDAL);
              }
        }
        break;
      default:
        break;
      }     
      break;
    case PACKET_TX::CONNECTION_PEDAL:
      memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead); //sizeof(aux_buffer));
      //Serial.println("DONGLE - arrivato richiesta di connessione");
      switch (aux_buffer[0]) {
        // =========================== NUOVI ==================================
        case CONNECTION_STATE::TX_GUN_SEARCH_PEDAL_BROADCAST:
          if ((SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) &&
              (aux_buffer[13] == espnow_wifi_channel))  // sistema espnow_wifi_channel
          { // prende la prima gun disposnibile
            memcpy(SerialWireless.mac_esp_another_card, &aux_buffer[1], 6);
            // invia richiesta connessione
            aux_buffer[0] = CONNECTION_STATE::TX_PEDAL_TO_GUN_PRESENCE;
            memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
            memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
            /////aux_buffer[13] = espnow_wifi_channel;  // no??????????????????????????????????????????????????????????????????????
            SerialWireless.SendPacket((const uint8_t *)aux_buffer, 13, PACKET_TX::CONNECTION_PEDAL);
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::TX_PEDAL_TO_GUN_PRESENCE;
            // assicurati che i dati siano stati spediti

          }
          break;
        case CONNECTION_STATE::TX_GUN_TO_PEDAL_ACCEPT:
          if ((memcmp(&aux_buffer[1],SerialWireless.mac_esp_another_card,6) == 0) && 
             (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
              SerialWireless.stato_connessione_wireless == CONNECTION_STATE::TX_PEDAL_TO_GUN_PRESENCE) {
              
                // SALVA I DATI RELATIVI ALLA GUN VID, PID, PLAYER , ECC.ECC.
              memcpy(&usb_data_wireless, &aux_buffer[13], sizeof(usb_data_wireless));
              
              aux_buffer[0] = CONNECTION_STATE::TX_PEDAL_TO_GUN_CONFERM;
              memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
              memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
                            
              // =========================================================
              // INVIARE IL PACCHETTO FINALE PIU' VOLTE  
              // =========================================================
              for (uint8_t i = 0; i<3; i++) {
                if (i>0) {
                  //vTaskDelay(pdMS_TO_TICKS(1000)); // equivalente a delay(1000) ma non bloccante su esp32
                  unsigned long lastMillis_tx_packet_gun_to_dongle_conferm = millis();
                  while ((millis() - lastMillis_tx_packet_gun_to_dongle_conferm) < 70) taskYIELD();  ////yield(); 
                }
                SerialWireless.SendPacket((const uint8_t *)aux_buffer, 13, PACKET_TX::CONNECTION_PEDAL);
              }
              SerialWireless.stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED;
          }
          break;          
        // =========================== FINE NUOVI =============================

        default:
          break;
      }
      break;  
    default:
      break;
  }
}

static void _esp_now_rx_cb(const esp_now_recv_info_t *info, const uint8_t *data, int len) {
  
  // 1. Filtro broadcast (invariato)
  if ((info->des_addr[0] == 0xFF) && (broadcast_receiver == false)) return;

  // 2. Cache locale degli indici (Il cuore del Lock-Free)
  uint16_t w = SerialWireless._writer;
  uint16_t r = SerialWireless._reader; // Leggiamo dove si trova il task in questo istante
  const uint16_t MASK = FIFO_SIZE_READ - 1; 

  // 3. Calcolo matematico dello spazio disponibile (Regola del Meno Uno)
  // Questa formula magica gestisce il wrap-around da sola e garantisce che il 
  // buffer non diventi mai pieno al 100%, lasciando sempre 1 byte di cuscinetto
  // per distinguere lo stato "vuoto" da quello "pieno".
  uint16_t free_space = (r - w - 1) & MASK;

  // 4. Controllo capienza
  if (free_space >= len) {
    uint16_t firstChunk = FIFO_SIZE_READ - w;

    // COPIA BATCH
    if (firstChunk < len) {
      // Caso a capo: dobbiamo dividere la copia in due
      uint16_t secondChunk = len - firstChunk;
      memcpy(SerialWireless._queue + w, data, firstChunk);
      memcpy(SerialWireless._queue, data + firstChunk, secondChunk);
    } 
    else {
      // Caso lineare: copia singola
      memcpy(SerialWireless._queue + w, data, len);
    }

    // --- BARRIERA DI MEMORIA ---
    // Fondamentale: garantisce che la RAM abbia salvato i dati della memcpy 
    // PRIMA di segnalare al task che ci sono nuovi dati disponibili.
    asm volatile ("memw" : : : "memory");

    // 5. Aggiornamento ATOMICO del solo indice _writer
    // Nessun lock necessario! Appena questa riga viene eseguita, 
    // il task (che sta guardando _writer) vede magicamente apparire i nuovi dati.
    SerialWireless._writer = (w + len) & MASK;
    
    // ADDIO SerialWireless._readLen!
  } 
  else {
    // Spazio insufficiente (o pacchetto scartato per evitare l'incrocio dei puntatori)
    SerialWireless._overflow_read = true;
  }

  // 6. Notifica Task (Campanello invariato)
  BaseType_t xHigherPriorityTaskWoken = pdFALSE;
  vTaskNotifyGiveFromISR(xRadioTaskHandle, &xHigherPriorityTaskWoken);
  if (xHigherPriorityTaskWoken) {
    portYIELD_FROM_ISR();
  }
}

static void _esp_now_tx_cb(const esp_now_send_info_t *tx_info, esp_now_send_status_t status) { 
  // 1. Liberiamo la radio per il prossimo pacchetto
  radioFree = true;

  // 2. Controllo Lock-Free: Il buffer TX ha ancora dati da spedire?
  if (SerialWireless.writeIndex != SerialWireless.readIndex) {   
    
    // 3. Svegliamo il radioTask per fargli svuotare il resto
    BaseType_t xHigherPriorityTaskWoken = pdFALSE;
    vTaskNotifyGiveFromISR(xRadioTaskHandle, &xHigherPriorityTaskWoken);
    
    if (xHigherPriorityTaskWoken) {
      portYIELD_FROM_ISR();
    }
  }
}

// ================================== TIMER PER SERIALE ===================================
// CALLBACK

void timer_callback_serial(void* arg) {
  while (SerialWireless.availableBufferSerialWrite() > 0) {
    if (!SerialWireless.flush_sem()) {
      // LA RADIO È PIENA! Il timer non può fare yield.
      // Chiediamo aiuto al radioTask per quando si libererà spazio.
      SerialWireless._serial_needs_recovery = true; 
      break; 
    }
  }
}

// ===============================================================================

void SerialWireless_::setupTimerSerial() {
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
    //esp_timer_stop(timer_handle_serial);
    //esp_timer_start_once(timer_handle_serial, duration_us);
    esp_timer_restart(timer_handle_serial, duration_us);
}

// ==========================   FINE TIMER PER SERIALE ====================================


// ================================== TIMER PER PEDAL ===================================
// CALLBACK
void timer_callback_pedal(void* arg) {
  
  TinyUSBDevices.pedals_wireless_state = 0;
 
}
// ===============================================================================

void setupTimerPedal() {
    esp_timer_create_args_t timer_args = {
        .callback = &timer_callback_pedal,
        .arg = NULL,
        .dispatch_method = ESP_TIMER_TASK,
        .name = "timer_pedal"
    };
    
    esp_timer_create(&timer_args, &timer_handle_pedal);
}

// ======================== FINE TIMER PER PEDAL =============================

// NON VA DEFINITA IN QUANTO SOSTITUISCE QUELLA DI DEFAULT CHE NON FA NULLA

#ifdef __cplusplus
extern "C" {
#endif

void tud_hid_report_complete_cb(uint8_t instance, uint8_t const* report, uint16_t len) {
    // Il PC ha appena liberato il buffer. Svegliamo usbTask!
    BaseType_t xHigherPriorityTaskWoken = pdFALSE;
    vTaskNotifyGiveFromISR(xUSBTaskHandle, &xHigherPriorityTaskWoken);
    if (xHigherPriorityTaskWoken) {
        portYIELD_FROM_ISR();
    }
}

#ifdef __cplusplus
}
#endif


#endif //OPENFIRE_WIRELESS_ENABLE