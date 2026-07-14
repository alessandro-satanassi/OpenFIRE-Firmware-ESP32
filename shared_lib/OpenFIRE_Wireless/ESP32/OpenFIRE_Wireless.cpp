/*!
 * @file OpenFIRE_Wireless.cpp
 * @brief Library for OpenFIRE wireless with Esp-Now
 * @n CPP Library for OpenFIRE wireless with Esp-Now
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)

#include "OpenFIRE_Wireless.h"
#include "esp_idf_version.h"

// ===================================================================================
// ROLE DEPENDENT HARDWARE INITIALIZATION / INIZIALIZZAZIONE HARDWARE DIPENDENTE DAL RUOLO (GUN vs DONGLE)
// ===================================================================================
// We selectively include display drivers based on the build target
// (#ifdef DONGLE or #ifdef GUN) to save precious space on Flash memory
// and do not load graphics libraries on devices that do not have them.
// /
// Includiamo selettivamente i driver del display in base al target di compilazione
// (#ifdef DONGLE o #ifdef GUN) per risparmiare prezioso spazio sulla memoria Flash
// e non caricare librerie grafiche su dispositivi che ne sono sprovvisti.

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
    // #include "../../src/OpenFIREdisplay.h"
    // #include "../../src/OpenFIREcommon.h"
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


// ===================================================================================
// ESP-NOW Configuration and Global Variables / CONFIGURAZIONE ESP-NOW E VARIABILI GLOBALI
// ===================================================================================

#ifndef OPENFIRE_ESPNOW_WIFI_CHANNEL
  #define OPENFIRE_ESPNOW_WIFI_CHANNEL 11 // default channel for the light gun  / canale sul quale si sintonizza la lightgun di default
#endif // OPENFIRE_ESPNOW_WIFI_CHANNEL

#ifndef OPENFIRE_ESPNOW_WIFI_POWER
  // transmit power can range from 8 to 84, where 84 is the maximum value corresponding to 21 dB (x4) / la potenza di trasmissione può andare da 8 a 84, dove 84 è il valore massimo che corrisponde a 21 db (x4)
  // #define OPENFIRE_ESPNOW_WIFI_POWER WIFI_POWER_15dBm  // corresponds to 60... more conservative power but with an excellent signal / corrisponde a 60 .. potenza più conservativa ma con segnale eccellente
  #define OPENFIRE_ESPNOW_WIFI_POWER WIFI_POWER_17dBm  // corresponds to 68... maximum power without distortion and with an excellent signal / corrisponde a 68 .. massima potenza senza distursioni con segnale eccellette
  
#endif //OPENFIRE_ESPNOW_WIFI_POWER

uint8_t espnow_wifi_channel = OPENFIRE_ESPNOW_WIFI_CHANNEL;  // Configurable for future setup via app or OLED. / FATTA VARIABILE PER FUTURA CONFIGURAZIONE TRAMITE APP O OLED
uint8_t espnow_wifi_power = OPENFIRE_ESPNOW_WIFI_POWER;      // Configurable for future setup via app or OLED. / FATTA VARIABILE PER FUTURA CONFIGURAZIONE TRAMITE APP O OLED

// HID SHADOW STATE: We maintain a local copy of the last state sent to the PC.
// We update the USB host *only* when this data changes, avoiding saturating
// the USB endpoint with unnecessary polling.
// /
// SHADOW STATE HID: Manteniamo una copia locale dell'ultimo stato inviato al PC.
// Aggiorniamo l'host USB *solo* quando questi dati cambiano, evitando di saturare
// l'Endpoint USB con polling inutili.
hid_abs_mouse_report_t absmouse5Report_last_wifi = {0,0,0,0,0};
hid_keyboard_report_t  keyReport_last_wifi = {0,0,{0,0,0,0,0,0}};
hid_gamepad16_report_t gamepad16Report_last_wifi = {0,0,0,0,0,0,0,0};

volatile bool absmouse5Report_pending = false; 
volatile bool keyReport_pending = false;
volatile bool gamepad16Report_pending = false;

// SPOOFING USB: 
// Default dongle data, which will be overwritten by the Lightgun's
// data during the handshake to fool the OS.
// /
//Dati di default del Dongle, che verranno sovrascritti dai dati della
// Lightgun durante la negoziazione (Handshake) per ingannare il SO.
USB_Data_GUN_Wireless usb_data_wireless = {
  "OpenFIRE_DONGLE",  // MANIFACTURES
  "FIRECon",          // NAME
  0xF143,             // VID
  0x1998, // 0x0001   // PID
  1,                  // PLAYER
  espnow_wifi_channel
  //OPENFIRE_ESPNOW_WIFI_CHANNEL // CHANNEL
  //,""               // 
};

SerialWireless_ SerialWireless;
extern Adafruit_USBD_HID usbHid;

// Spinlocks (Hardware Mutexes): Essential for ensuring that reads and writes
// to circular buffers are never interrupted mid-operation by tasks on different cores
// or by WiFi interrupts, preventing fatal RAM corruption.
// /
// Spinlocks (Mutex Hardware): Fondamentali per garantire che letture e scritture
// sui buffer circolari non vengano mai interrotte a metà da task su core diversi 
// o da interrupt WiFi, prevenendo la corruzione letale della RAM.
static portMUX_TYPE mux_serial_tx = portMUX_INITIALIZER_UNLOCKED;
static portMUX_TYPE mux_radio_tx = portMUX_INITIALIZER_UNLOCKED;


uint8_t buffer_espnow[ESP_NOW_MAX_DATA_LEN];
esp_now_peer_info_t peerInfo; // It must be placed outside the functions—a global utility variable for configuration. / deve stare fuori funzioni da funzioni -- globale --variabile di utilità per configurazione

static void _esp_now_rx_cb(const esp_now_recv_info_t *info, const uint8_t *data, int len); // callback esp_now
static void _esp_now_tx_cb(const esp_now_send_info_t *tx_info, esp_now_send_status_t status);

#ifdef COMMENTO
#if ESP_IDF_VERSION_MAJOR == 5 && ESP_IDF_VERSION_MINOR <= 4
  // Code for versions up to 5.4 / Codice per versioni fino alla 5.4
  static void _esp_now_tx_cb(const uint8_t *mac_addr, esp_now_send_status_t status); // callback esp_now
#elif ESP_IDF_VERSION_MAJOR > 5 || (ESP_IDF_VERSION_MAJOR == 5 && ESP_IDF_VERSION_MINOR >= 5)
  // Code for versions 5.5 or higher / Codice per versioni 5.5 o superiori
  static void _esp_now_tx_cb(const esp_now_send_info_t *tx_info, esp_now_send_status_t status); // callback esp_now // for the new IDF 55 version / poer nuova versione IDF 55
#endif
#endif // COMMENTO

void packet_callback_read_dongle(); // callback packet 
void packet_callback_read_gun(); // callback packet
void packet_callback_read_pedal(); // callback packet

// ===================================================================================
// Addressing Management and Network State Machine / GESTIONE INDIRIZZAMENTO E MACCHINA A STATI DI RETE
// ===================================================================================

///////////////////////////////////////////////////////////////////
// #define SIZE_BASE_AUX 13
uint8_t aux_buffer[13 + sizeof(usb_data_wireless)];
// 1 byte indicating connection establishment progress (current status)
// 6 bytes: MAC address of the transmitting device
// 6 bytes: MAC address of the device to which the packet is addressed
// On the final transmission, the gun also sends data from the USB_DATA_WIRELESS structure (containing the gun's VID, PID, etc.) so the dongle can initialize the USB connection using that data
// /
// 1 BYTE INDICATA L'AVANZAMENTO DELLO STABILIMENTO DELLA CONNESSIONE,A CHE PUNTO SIAMO
// 6 BYTE IL MAC ADDRESS DEL DISPOSITIVO CHE TRASMETTE
// 6 BYTE IL MAC ADDRESS DEL DISPOSITIVO A CUI E' INDIRIZZATO IL PACCHETTO
// ALL'ULTIMO INVIO LA GUN INVIA ANCHE I DATI DELLA STRUTTURA USB_DATA_WIRELESS CHE CONTIENE VID, PID, ECC. DELLA GUN IN MODO CHE IL DONGLE INZIALIZZI LA CONSSESIONE USB CON QUEI DATI

// Connection status
// /
// stato della connessione 

// For the GUN
// 0 = Not established
// 1 = Broadcast request sent
// 2 = Connection response received (request to connect)
// 3 = Connection confirmation sent
// /
// in caso di GUN
// 0 = non stabilita
// 1 = inviata richiesta blroad cast
// 2 = qualcuno ha risposto alla connessione ed ha chiesto di essere conesso
// 3 = inviato conferma connessione;

// For the DONGLE
// 0 = Not established
// 1 = Broadcast connection request received
// 2 = Accepted; connection request sent
// 3 = Connection acceptance received
// /
// in caso di DONGLE
// 0 = non stabilita
// 1 = arrivata richiesta broadcast di connessione
// 2 = accettata ed inviata richiesta di connessione
// 3 = arrivara accettazione di connessione

const uint8_t BROADCAST_ADDR[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};

#if defined(DONGLE)
  uint8_t peerAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}; // broadcast
  const functionPtr callbackArr[] = { packet_callback_read_dongle };
  uint8_t lastDongleAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  uint8_t lastPedalAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  uint8_t lastDongleChannel = OPENFIRE_ESPNOW_WIFI_CHANNEL;
  uint8_t lastPedalChannel = OPENFIRE_ESPNOW_WIFI_CHANNEL;
  bool lastDongleSave = false; // If true, it means we have the address of the last dongle; otherwise, it is false. / se true significa che abbiamo un indirizzo dell'ultimo dongle altrimenti false
  bool lastPedalSave = false; 
#elif defined(PEDAL)
  uint8_t peerAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}; // broadcast
  const functionPtr callbackArr[] = { packet_callback_read_pedal };
  uint8_t lastDongleAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  uint8_t lastPedalAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  uint8_t lastDongleChannel = OPENFIRE_ESPNOW_WIFI_CHANNEL;
  uint8_t lastPedalChannel = OPENFIRE_ESPNOW_WIFI_CHANNEL;
  bool lastDongleSave = false; // If true, it means we have the address of the last dongle; otherwise, it is false. / se true significa che abbiamo un indirizzo dell'ultimo dongle altrimenti false
  bool lastPedalSave = false; 
#elif defined(GUN)
  uint8_t peerAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF}; // broadcast
  const functionPtr callbackArr[] = { packet_callback_read_gun };
  uint8_t lastDongleAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  uint8_t lastPedalAddress[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  uint8_t lastDongleChannel = OPENFIRE_ESPNOW_WIFI_CHANNEL;
  uint8_t lastPedalChannel = OPENFIRE_ESPNOW_WIFI_CHANNEL;
  bool lastDongleSave = false; // if true, it means the last dongle configuration was saved / se true significa
  bool lastPedalSave = false; // if true, it means the last pedal configuration was saved / se true significa
  // ==========================================================================================================
#endif

// ========= per gestione PEDAL ===========
uint8_t peerAddress_pedal[6] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00}; // broadcast
#define MAX_TIMEOUT_LAST_PACKET (uint64_t)100000   // = 100  ms  .. the value must be specified in microseconds. / il valore va specificato in microsecondi
void setupTimerPedal();
esp_timer_handle_t timer_handle_pedal;
// ===================================================

// ===================================================================================
// HARDWARE SETTINGS: RADIO FREQUENCY AND PHY / FREQUENZA RADIO E PHY
// ===================================================================================
///////////////////////////////////////////////////////////////////
// Overall definition of the rate configuration / Definizione globale della configurazione del rate
esp_now_rate_config_t rate_config = {
    .phymode = WIFI_PHY_MODE_11G, // Uses the G standard (best compromise between penetration and speed) / Forza lo standard G (Miglior compromesso penetrazione/velocità)
    .rate = WIFI_PHY_RATE_12M,    // Forces 12 Mbps (Ensures link stability at the expense of unnecessary bandwidth) / Forza i 12 Mbps (Garantisce stabilità del link a discapito della banda larga non necessaria)
    .ersu = false,                // Not necessary for 11g / Non necessario per 11g
    .dcm = false                  // Not necessary for 11g / Non necessario per 11g
};
////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////

volatile uint8_t channel_display = espnow_wifi_channel;
volatile uint8_t seconds_display = 30;
volatile bool broadcast_receiver = false;  

// ===================================================================================
// FREERTOS ASYNCHRONOUS TASKS (USB AND RADIO) / TASKS ASINCRONI FREERTOS (USB E RADIO)
// ===================================================================================

//////////////////////////////////////////////////////////////////////////////
TaskHandle_t xUSBTaskHandle = NULL;
#ifdef DONGLE

// The USB task manages the transmission of HID data to the PC, ensuring that the endpoint
// never becomes saturated. It employs a round-robin algorithm to rotate priority
// among the mouse, keyboard, and gamepad, preventing "chatty" peripherals (such as the mouse)
// from overshadowing the other controls.
// /
// Il task USB gestisce l'invio dei dati HID al PC garantendo che l'Endpoint non 
// venga mai saturato. Utilizza un algoritmo Round-Robin per dare priorità a rotazione 
// tra Mouse, Tastiera e Gamepad, evitando che periferiche "loquaci" (come il mouse) 
// oscurino gli altri controlli.

void usbTask(void *pvParameters) {
  // Keep track of whose turn it is. 
  // 0 = Mouse, 1 = Keyboard, 2 = Gamepad
  // /
  // Teniamo traccia di chi ha il turno. 
  // 0 = Mouse, 1 = Tastiera, 2 = Gamepad
  static int currentRotation = 0; 

  for (;;) {
    // PAUSE THE TASK
    // It will wake up immediately if:
    // - The radio receives new data (your xTaskNotifyGive)
    // - USB transmission is complete (the tud_hid_report_complete_cb callback)
    // /
    // METTI IN PAUSA IL TASK
    // Si sveglierà istantaneamente se:
    // - La radio riceve un dato nuovo (il tuo xTaskNotifyGive)
    // - La USB ha finito di trasmettere (la callback tud_hid_report_complete_cb)
    ulTaskNotifyTake(pdTRUE, portMAX_DELAY);

    // HARDWARE CHECK
    // Proceed only if the HID endpoint is physically free
    // /
    // CONTROLLO HARDWARE
    // Procediamo solo se l'Endpoint HID è fisicamente libero
    if (usbHid.ready()) {

      // ROUND-ROBIN ROTATION
      // Loop up to 3 times to find the first device with pending data
      // /
      // ROTAZIONE ROUND-ROBIN
      // Facciamo un ciclo di massimo 3 tentativi per trovare il primo dispositivo con dati pendenti
      for (int i = 0; i < 3; i++) {
        int target = (currentRotation + i) % 3;

        // --- TURN 0: MOUSE / TURNO 0: MOUSE ---
        if (target == 0 && absmouse5Report_pending) {
          if (usbHid.sendReport(HID_RID_e::HID_RID_MOUSE, &absmouse5Report_last_wifi, sizeof(absmouse5Report_last_wifi))) {
            absmouse5Report_pending = false; // Data processed / Dato smaltito
            currentRotation = 1;            // Next turn goes to the keyboard / Prossimo giro tocca alla Tastiera
            break; // IMPORTANT: Exit the for loop. The USB buffer is now full. / IMPORTANTE: Usciamo dal for. Il buffer USB ora è pieno.
          }
        }
        
        // --- TURN 1: KEYBOARD / TURNO 1: TASTIERA ---
        else if (target == 1 && keyReport_pending) {
          if (usbHid.sendReport(HID_RID_e::HID_RID_KEYBOARD, &keyReport_last_wifi, sizeof(keyReport_last_wifi))) {
            keyReport_pending = false; 
            currentRotation = 2;            // Next up: the Gamepad / Prossimo giro tocca al Gamepad
            break; 
          }
        }
        
        // --- TURN 2: GAMEPAD / TURNO 2: GAMEPAD ---
        else if (target == 2 && gamepad16Report_pending) {
          if (usbHid.sendReport(HID_RID_e::HID_RID_GAMEPAD, &gamepad16Report_last_wifi, sizeof(gamepad16Report_last_wifi))) {
            gamepad16Report_pending = false;
            currentRotation = 0;            // Next round, we start again from the Mouse. / Prossimo giro si riparte dal Mouse
            break; 
          }
        }
        
      } // End of the for loop / Fine del ciclo for
    }
    // If we reach this point (whether the USB was busy or we sent a report),
    // the task completes the iteration and goes back to sleep at the start of the for(;;) loop.
    // /
    // Se arriviamo qui (che la USB fosse occupata o che abbiamo inviato un report),
    // il task finisce il giro e torna a dormire all'inizio del for(;;)
  }
}
#endif // DONGLE
/////////////////////////////////////////////////////////////////////////////

TaskHandle_t xRadioTaskHandle = NULL;
volatile bool radioFree = true;

// The Radio task acts as the traffic controller. It remains "asleep" until
// it receives a notification (ESP-NOW ISR or Serial Timer).
// It handles both the parsing of incoming packets (bypassing the FSM)
// and the transmission of queued packets over the air.
// /
// Il task della Radio è il controllore del traffico. Rimane "addormentato" fino
// a quando non riceve una notifica (ISR ESP-NOW o Timer Seriale).
// Gestisce sia il Parsing dei pacchetti in arrivo (scavalcando la FSM), 
// sia l'emissione dei pacchetti in coda verso l'etere.

void radioTask(void *pvParameters) {
  while (1) {
    ulTaskNotifyTake(pdTRUE, portMAX_DELAY);
    
    bool progress;
    
    do {
      progress = false;

      // RECEPTION MANAGEMENT / GESTIONE RICEZIONE
      if (SerialWireless._reader != SerialWireless._writer) {
        SerialWireless.checkForRxPacket();
        progress = true;
      }

      // SERIAL MANAGEMENT -> RECOVERY MODE! / GESTIONE SERIALE -> RECOVERY MODE!
      // intervene ONLY if the timer has previously failed. / Interveniamo SOLO se il timer ha fallito in precedenza
      if (SerialWireless._serial_needs_recovery && (SerialWireless.availableBufferSerialWrite() > 0)) {
        if (SerialWireless.flush_sem()) {
          // Recovery successful! Lower the flag. 
          // SerialWireless._serial_needs_recovery = false; // flush_sem() already does this itself before returning true.
          // /
          // Salvataggio riuscito! Abbassiamo la bandierina.
          // SerialWireless._serial_needs_recovery = false; // flush_sem() lo fa già da sola prima di restituire true.
          progress = true; 
        }
      }

      // GESTIONE TRASMISSIONE RADIO
      if (radioFree && (SerialWireless.readIndex != SerialWireless.writeIndex)) {
        SerialWireless.SendData_sem();
        progress = true;
      }

    } while (progress); 
  } 
}
// =====================================

// ===================================================================================
// ASYNCHRONOUS GRAPHIC ANIMATIONS (TFT/OLED DISPLAY) / ANIMAZIONI GRAFICHE ASINCRONE (DISPLAY TFT/OLED)
// ===================================================================================
// These functions are executed on separate cores during "blocking" network phases
// to ensure the screen remains smooth and responsive, providing clear visual
// feedback to the user regarding the pairing process (UX pattern).
// /
// Queste funzioni sono lanciate su Core separati durante le fasi "bloccanti" di rete 
// per garantire che lo schermo sia fluido e reattivo, fornendo un chiaro feedback
// visuale sull'accoppiamento all'utente (UX Pattern).

/////////////////////////////////////////////////////////////////////////////
// GRAPHICAL CHANNEL SEARCH FOR GUN / VERSIONE GRAFICA RICERCA CANALE PER GUN //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#if defined(GUN) && defined(USES_DISPLAY)
void animTaskLink_pedal(void *pvParameters) {
  const char* rotazione = "-\\|/-\\|/";
  int8_t rotazioneIndex = 0;
  const uint8_t baseX = 1;
  const uint8_t baseY = 2;
  const uint8_t charWidth = 6;
  uint8_t seconds_tread = seconds_display;
  
  // Initial display setup / Setup iniziale display
  display_OLED->setCursor(baseX, baseY);
  display_OLED->setTextSize(1);
  display_OLED->setTextColor(WHITE, BLACK);
  display_OLED->fillRect(0, 0, 128, 16, BLACK);
  display_OLED->drawFastHLine(0, 15, 128, WHITE);
  char buffer[30];
  sprintf(buffer, "Search Pedal [%2ds]", seconds_tread);                  
  display_OLED->print(buffer);
  display_OLED->display();

  // Task loop / Loop del task
  for (;;) {
    if (seconds_display != seconds_tread)
    {
      display_OLED->setCursor((14*6)+1, 2);
      seconds_tread = seconds_display;
      sprintf(buffer, "%2d", seconds_tread);
      display_OLED->print(buffer);
    }
    display_OLED->setCursor(baseX + (19 * charWidth), baseY);
    display_OLED->write(rotazione[rotazioneIndex]);
    rotazioneIndex = (rotazioneIndex + 1) % 8;
    display_OLED->display();
    vTaskDelay(pdMS_TO_TICKS(150)); // a short delay to avoid saturating the CPU / piccolo delay per non saturare la CPU
  } 
}

void animTaskLink(void *pvParameters) {
  const char* rotazione = "-\\|/-\\|/";
  int8_t rotazioneIndex = 0;
  const uint8_t baseX = 1;
  const uint8_t baseY = 2;
  const uint8_t charWidth = 6;
  uint8_t channel_tread = espnow_wifi_channel;
  
  // Initial display setup / Setup iniziale display
  display_OLED->setCursor(baseX, baseY);
  display_OLED->setTextSize(1);
  display_OLED->setTextColor(WHITE, BLACK);
  display_OLED->fillRect(0, 0, 128, 16, BLACK);
  display_OLED->drawFastHLine(0, 15, 128, WHITE);
  char buffer[30];
  sprintf(buffer, "Ch:%2d Waiting link", channel_tread);
  display_OLED->print(buffer);
  display_OLED->display();

  // Task loop / Loop del task
  for (;;) {
    if (channel_display != channel_tread)
    {
      display_OLED->setCursor(19, 2);
      channel_tread = channel_display;
      sprintf(buffer, "%2d", channel_tread);
      display_OLED->print(buffer);
    }
    display_OLED->setCursor(baseX + (19 * charWidth), baseY);
    display_OLED->write(rotazione[rotazioneIndex]);
    rotazioneIndex = (rotazioneIndex + 1) % 8;
    display_OLED->display();
    vTaskDelay(pdMS_TO_TICKS(150)); // a short delay to avoid saturating the CPU / piccolo delay per non saturare la CPU
  } 
}
#endif // GUN



/////////////////////////////////////////////////////////////////////////////
// END OF GRAPHIC VERSION CHANNEL SEARCH FOR GUN / FINE VERSIONE GRAFICA RICERCA CANALE PER GUN /////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Graphical Channel Search Version for Dongle / VERSIONE GRAFICA RICERCA CANALE PER DONGLE ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#if defined(DONGLE) && defined(USES_DISPLAY)
void animTaskLink(void *pvParameters) {
  const char* rotazione = "-\\|/-\\|/";
  int8_t rotazioneIndex = 0;
  const uint8_t baseX = 75;
  const uint8_t baseY = 58;
  const uint8_t charWidth = 6;
  
  // Initial display setup / Setup iniziale display
  tft.setTextSize(1);
  tft.setTextColor(WHITE, BLACK);
  tft.fillRect(70,20,100,60,BLACK);
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
  #ifdef USE_LOVYAN_GFX
  tft.display();
  #endif // USE_LOVYAN_GFX
  tft.setTextSize(2);
  tft.setTextColor(BLUE, BLACK);

  // Task loop / Loop del task
  for (;;) {    
    tft.setCursor(baseX + (11 * charWidth) + 3, baseY);
    tft.write(rotazione[rotazioneIndex]);
    rotazioneIndex = (rotazioneIndex + 1) % 8;
    #ifdef USE_LOVYAN_GFX
    tft.display();
    #endif // USE_LOVYAN_GFX
    vTaskDelay(pdMS_TO_TICKS(150)); // piccolo delay per non saturare la CPU
  } 
}

void animTask(void *pvParameters) {
  int8_t currentIndex = 0;
  const uint8_t baseX = 75;
  const uint8_t baseY = 45;
  const uint8_t charWidth = 6;

  // Initial display setup / Setup iniziale display
  tft.setCursor(baseX, baseY);
  tft.setTextSize(1);
  tft.setTextColor(WHITE, BLACK);
  tft.setCursor(baseX + 12, baseY - 20);
  tft.print("Searching");
  tft.setCursor(baseX+3, baseY - 10);
  tft.print("best channel");
  tft.setCursor(baseX, baseY);
  tft.print("Scanning WiFi");
  #ifdef USE_LOVYAN_GFX
  tft.display();
  #endif // USE_LOVYAN_GFX
  tft.setTextColor(RED, BLACK);

  // Task loop / Loop del task
  for (;;) {  
    tft.setCursor(baseX + (3 * charWidth), baseY + 18);
    switch (currentIndex)
    {
    case 0:
      tft.print("   o   ");
      break;
    case 1:
      tft.print("  (o)  ");
      break;
    case 2:
      tft.print(" ((o)) ");
      break;
    case 3:
      tft.print("(((o)))");
      break;
    default:
    break;
    }
    currentIndex = (currentIndex + 1) % 4;
    #ifdef USE_LOVYAN_GFX
    tft.display();
    #endif // USE_LOVYAN_GFX
    vTaskDelay(pdMS_TO_TICKS(150)); // piccolo delay per non saturare la CPU
  }
}
#endif //DONGLE


/////////////////////////////////////////////////////////////////////////////
// End of graphical channel search version for dongle / FINE VERSIONE GRAFICA RICERCA CANALE PER DONGLE //////////////////////////
/////////////////////////////////////////////////////////////////////////////

// ===================================================================================
// Heuristic Best Channel Search (Auto-Channel Configuration) / RICERCA EURISTICA CANALE MIGLIORE (AUTO-CHANNEL CONFIGURATION)
// ===================================================================================
// This function, executed by the dongle, activates the radio chip's "promiscuous" mode
// to sniff the airwaves. It assigns a score to each channel based not only on the
// number of networks but, more importantly, on the actual noise floor and the
// packet traffic present at that moment. It selects the channel offering the
// best Signal-to-Noise Ratio (SNR) for ESP-NOW.
// /
// Questa funzione eseguita dal Dongle attiva la modalità "Promiscuous" del chip radio
// per annusare l'etere (sniffing). Assegna uno score a ciascun canale basato non
// solo sul numero di reti, ma soprattutto sul rumore reale (Noise Floor) e sul
// traffico effettivo dei pacchetti che transitano in quell'istante. Sceglie il 
// canale che offre il SNR (Signal-to-Noise Ratio) migliore per l'ESP-NOW.

// ================= BEST CHANNEL SEARCH FUNCTIONS - USED ONLY BY THE DONGLE / FUNZIONI PER RICERCA MIGLIORE CANALE - USATO SOLO DAL DONGLE =========================
#ifdef DONGLE

// ================= GLOBAL VARIABLES FOR CALLBACKS / VARIABILI GLOBALI PER CALLBACK =================
static volatile uint32_t g_packetCounter = 0;
static volatile bool g_sniffing = false;

// ================= CALLBACK PROMISCUOUS MODE =================
void IRAM_ATTR promiscuousCallback(void *buf, wifi_promiscuous_pkt_type_t type) {
  if (!g_sniffing) return;
  g_packetCounter = g_packetCounter + 1;
}

// ================= MAIN FUNCTION / FUNZIONE PRINCIPALE =================
uint8_t findBestChannel() {
  
  #if defined(DONGLE) && defined(USES_DISPLAY)
  TaskHandle_t animTaskHandle = NULL;  
  if(display_init) {
  // Start animation / Avvio animazione
    if (animTaskHandle == NULL) {
      xTaskCreatePinnedToCore(
        animTask,          // task function / funzione del task
        "AnimTask",        // name / nome
        4096,              // stack size
        NULL,              // parameters / parametri
        1,                 // priority / priorità
        &animTaskHandle,   // handle
        APP_CPU_NUM        // core (you can use 0 or 1) / (puoi usare 0 o 1)
      );
    }
  }
  #endif // USES_DISPLAY


  // ================= STRUCTURE FOR CHANNEL STATISTICS & CHANNEL SEARCH / STRUTTURA PER STATISTICHE CANALE PER RICERCA BEST CHANNEL=================
  typedef struct {
    uint16_t networks;   // Number of APs detected / Numero di AP rilevati
    float avgRSSI;       // Average AP RSSI / RSSI medio degli AP
    int8_t maxRSSI;      // Maximum RSSI (strongest interference) / RSSI massimo (interferenza più forte)
    uint32_t packets;    // Packets captured in 600ms / Pacchetti catturati in 600ms
    float noise;         // Noise floor estimate / Stima noise floor
    float score;         // Base channel score / Score base del canale
  } ChannelStats;
  
  ChannelStats channelStats[14];
  memset(channelStats, 0, sizeof(channelStats));
  g_packetCounter = 0;
  g_sniffing = false;
  for (int i = 0; i < 14; i++) channelStats[i].maxRSSI = -128;

  // ================= PHASE 1: WiFi Network Scan / FASE 1: SCAN RETI WiFi =================
  // Time: ~2 seconds (150ms × 13 channels)
  // /
  // Tempo: ~2 secondi (150ms × 13 canali)
  WiFi.mode(WIFI_STA);
  WiFi.disconnect();
  vTaskDelay(pdMS_TO_TICKS(100));

  // Optimized scan: 150ms per channel captures all beacons
  // (typical beacon interval: 100ms)
  // /
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
    
    // Calculate average RSSI
    // /
    // Calcola RSSI medio
    for (int ch = 1; ch <= 13; ch++) {
      if (channelStats[ch].networks > 0) {
        channelStats[ch].avgRSSI /= channelStats[ch].networks;
      } else {
        channelStats[ch].maxRSSI = -100; // Empty channel / Canale vuoto
      }
    }
  }
  
  WiFi.scanDelete();
  vTaskDelay(pdMS_TO_TICKS(50));

  // ================= PHASE 2: TRAFFIC SNIFFING AND NOISE / FASE 2: SNIFF TRAFFICO E NOISE =================
  // Time: ~12 seconds (920ms × 13 channels)
  // /
  // Tempo: ~12 secondi (920ms × 13 canali)
  
  vTaskDelay(pdMS_TO_TICKS(50));
  esp_wifi_set_promiscuous(true);
  esp_wifi_set_promiscuous_rx_cb(&promiscuousCallback);

  for (int ch = 1; ch <= 13; ch++) {

    esp_wifi_set_channel(ch, WIFI_SECOND_CHAN_NONE);
    vTaskDelay(pdMS_TO_TICKS(70)); // Radio stabilization / Stabilizzazione radio
    
    // --- Traffic measurement: 600ms full pattern capture ---
    // Includes: beacon (100ms), data burst, retry packets
    // /
    // --- Misurazione traffico: 600ms cattura pattern completi ---
    // Include: beacon (100ms), data burst, retry packets
    g_packetCounter = 0;
    g_sniffing = true;
    uint32_t startCount = g_packetCounter;
    vTaskDelay(pdMS_TO_TICKS(600));
    uint32_t endCount = g_packetCounter;
    channelStats[ch].packets = endCount - startCount;
    
    // --- Noise Floor Measurement: 15 samples for a stable average ---
    // /
    // --- Misurazione Noise Floor: 15 campioni per media stabile ---
    int32_t noiseSum = 0;
    const int samples = 15;
    for (int i = 0; i < samples; i++) {
      uint32_t start = g_packetCounter;
      vTaskDelay(pdMS_TO_TICKS(17)); // // Sampling 17ms (total 255ms) / Campionamento 17ms (totale 255ms)
      uint32_t activity = g_packetCounter - start;
      
      // Improved noise formula: logarithmic scale / Formula noise migliorata: scala logaritmica
      if (activity == 0) {
        noiseSum += -100; // Complete silence / Silenzio completo
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

  // ================= PHASE 3: OPTIMIZED SCORE CALCULATION / FASE 3: CALCOLO SCORE OTTIMIZZATO =================
  // Weights based on empirical analysis and Espressif documentation
  //
  // PRIORITIES FOR ESP-NOW:
  // 1. Strong RSSI = direct interference lowering SNR (KILLER #1)
  // 2. Actual traffic = packet collisions (KILLER #2)
  // 3. Number of networks = potential congestion (moderate)
  // 4. Noise floor = low impact (ESP-NOW uses robust MCS0/802.11b)
  // /
  // Pesi basati su analisi empirica e documentazione Espressif
  // 
  // priority / priorità PER ESP-NOW:
  // 1. RSSI forte = interferenza diretta che abbassa SNR (KILLER #1)
  // 2. Traffico reale = collisioni pacchetti (KILLER #2)
  // 3. Numero reti = congestione potenziale (moderato)
  // 4. Noise floor = basso impatto (ESP-NOW usa MCS0/802.11b robusto)
  
  const float wRSSI  = 0.40f;  // Direct interference = highest priority / Interferenza diretta = priorità massima
  const float wPkt   = 0.35f;  // Actual collisions = highly critical / Collisioni reali = molto critico
  const float wNet   = 0.20f;  // Congestion = moderate (log scale) / Congestione = moderato (scala log)
  const float wNoise = 0.05f;  // Noise = minimal (robust modulation) / Rumore = minimo (modulazione robusta)

  for (int ch = 1; ch <= 13; ch++) {
    
    // --- RSSI Score: most dangerous direct interference ---
    // AP with strong signal (-40dBm) nearby = guaranteed disaster
    // -100dBm (far away) = 0 points | -40dBm (nearby) = 40 points
    // /
    // --- Score RSSI: interferenza diretta più pericolosa ---
    // AP con segnale forte (-40dBm) vicino = disastro garantito
    // -100dBm (lontano) = 0 punti | -40dBm (vicino) = 40 punti
    float sRSSI = 0;
    if (channelStats[ch].maxRSSI > -100) {
      int rssiMapped = constrain(channelStats[ch].maxRSSI, -100, -30);
      sRSSI = map(rssiMapped, -100, -30, 0, 100) * wRSSI;
    }
    
    // --- Packet Score: measured real traffic ---
    // Normalized to 600ms: 0 pkts = 0 | 100 pkts = 35 points
    // /
    // --- Score Pacchetti: traffico reale misurato ---
    // Normalizzato su 600ms: 0 pkt = 0 | 100 pkt = 35 punti
    float sPkt = (channelStats[ch].packets / 10.0) * wPkt;
    if (sPkt > 100 * wPkt) sPkt = 100 * wPkt; // clamp max
    
    // --- Goals Score: logarithmic scale to avoid dominance ---
    // 0 goals = 0 | 1 goal = 4.2 | 5 goals = 9.2 | 10 goals = 13.8
    // /
    // --- Score Reti: scala logaritmica per evitare dominanza ---
    // 0 reti = 0 | 1 rete = 4.2 | 5 reti = 9.2 | 10 reti = 13.8
    float sNet = (channelStats[ch].networks > 0) 
                 ? (log(channelStats[ch].networks + 1) * 20.0) * wNet 
                 : 0;

    // --- RSSI Noise: minimum weight (ESP-NOW resistant) ---
    // -100dBm (silence) = 0 | -50dBm (noisy) = 5 points
    // /
    // --- Score Noise: peso minimo (ESP-NOW resistente) ---
    // -100dBm (silenzio) = 0 | -50dBm (rumoroso) = 5 punti
    int noiseMapped = constrain((int)channelStats[ch].noise, -100, -50);
    float sNoise = map(noiseMapped, -100, -50, 0, 100) * wNoise;
    
    channelStats[ch].score = sRSSI + sPkt + sNet + sNoise;
  }

  // ================= ADJACENT CHANNEL OVERLAP PENALTY / PENALITÀ OVERLAP CANALI ADIACENTI =================
  // WiFi 20MHz: each channel overlaps with ±4 channels
  // Reference: 802.11b uses 22MHz, 802.11g/n uses 20MHz + guardband
  // Non-overlapping channels: 1, 6, 11 (separated by 5 channels)
  //
  // Progressive penalties to simulate actual overlap:
  // ±1: 50% (maximum overlap, approx. 75% shared spectrum)
  // ±2: 35% (significant overlap, approx. 50% spectrum)
  // ±3: 20% (moderate overlap, approx. 25% spectrum)
  // ±4: 10% (minimal overlap, spectrum edges)
  // /
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
    
    // Full overlap penalty, 20MHz model / Penalità overlap completo modello 20MHz
    if (ch > 1) finalScores[ch] += channelStats[ch-1].score * 0.50f;
    if (ch > 2) finalScores[ch] += channelStats[ch-2].score * 0.35f;
    if (ch > 3) finalScores[ch] += channelStats[ch-3].score * 0.20f;
    if (ch > 4) finalScores[ch] += channelStats[ch-4].score * 0.10f;
    
    if (ch < 13) finalScores[ch] += channelStats[ch+1].score * 0.50f;
    if (ch < 12) finalScores[ch] += channelStats[ch+2].score * 0.35f;
    if (ch < 11) finalScores[ch] += channelStats[ch+3].score * 0.20f;
    if (ch < 10) finalScores[ch] += channelStats[ch+4].score * 0.10f;
  
    // After calculating finalScores, award bonuses to the ideal channels / Dopo il calcolo finalScores, dare bonus ai canali ideali
    if (ch == 1 || ch == 6 || ch == 11) {
      finalScores[ch] *= 0.92f; // -8% bonus(preference)/(preferenza) 
    }  
  }

  // ================= PHASE 4: SELECTION OF THE BEST CHANNEL / FASE 4: SELEZIONE MIGLIOR CANALE =================
  uint8_t bestCh = 1;
  float minScore = 999999.0;
  
  for (uint8_t ch = 1; ch <= 13; ch++) {
    if (finalScores[ch] < minScore) {
      minScore = finalScores[ch];
      bestCh = ch;
    }
  }
  
  // ================= Final cleaning / PULIZIA FINALE =================
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
#endif // DONGLE
// ================= End of best channel search functions – used only by the dongle. / FINE FUNZIONI PER RICERCA MIGLIORE CANALE - USATO SOLO DAL DONGLE =========================

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*****************************
 *   SERIAL WIRELESS SECTION
 *****************************/
// ===================================================================================
// STREAM INTERFACE: ABSTRACTION OF THE RADIO AS A VIRTUAL SERIAL PORT / INTERFACCIA STREAM: ASTRAZIONE DELLA RADIO COME PORTA SERIALE VIRTUALE
// ===================================================================================
// These methods override the standard Arduino `Stream` interface.
// The architectural goal is to allow the rest of the code (e.g. parsing logic
// or external libraries) to read and write data from/to the ESP-NOW network
// exactly as if it were a normal serial cable connection.
// Latency is masked by the use of Ring Buffers and Lock-Free Indexes.
// /
// Questi metodi sovrascrivono l'interfaccia standard `Stream` di Arduino. 
// L'obiettivo architettonico è permettere al resto del codice (es. logica di parsing 
// o librerie esterne) di leggere e scrivere dati dalla/verso la rete ESP-NOW 
// esattamente come se fosse una normale connessione seriale via cavo.
// La latenza è mascherata dall'uso di buffer circolari (Ring Buffers) e indici Lock-Free.

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

// Extracts a single byte from the circular receive buffer.
// Using the bitwise operator (& MASK_READ_SERIAL) replaces the costly modulo operator (%),
// ensuring index wrap-around in a single clock cycle.
// /
// Estrae un singolo byte dal buffer circolare di ricezione. 
// L'uso dell'operatore bitwise (& MASK_READ_SERIAL) sostituisce il costoso operatore modulo (%),
// garantendo il wrap-around dell'indice in un singolo ciclo di clock.
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

// ===================================================================================
// ASYNCHRONOUS PARSING ENGINE AND LOCK-FREE MEMORY MANAGEMENT / MOTORE DI PARSING ASINCRONO E GESTIONE MEMORIA LOCK-FREE
// ===================================================================================

// This function is called by the `radioTask`. It inspects the raw data queue
// (populated by the radio ISR) and immediately passes it to the COBS protocol
// state machine (`packet.parse`).
// /
// Questa funzione viene chiamata dal `radioTask`. Ispeziona la coda dati grezzi 
// (riempita dalla ISR della radio) e la passa istantaneamente alla macchina a stati 
// del protocollo COBS (`packet.parse`).
bool SerialWireless_::checkForRxPacket() {
    // 1. Local snapshot of the indices at this precise moment
    // This is the heart of the lock-free pattern: we read the current state without blocking
    // the hardware interrupt that might be writing simultaneously.
    // /
    // 1. Snapshot locale degli indici in questo preciso istante
    // Questo è il cuore del pattern Lock-Free: leggiamo lo stato corrente senza bloccare 
    // l'hardware interrupt che potrebbe star scrivendo contemporaneamente.
    uint16_t r = _reader;
    uint16_t w = _writer; // Let's look at where the ESP-NOW callback has ended up. / Guardiamo dove è arrivata la callback ESP-NOW
    const uint16_t MASK = FIFO_SIZE_READ - 1;

    // 2. Lock-free calculation of available data
    // This formula handles wrap-around automatically within a single clock cycle.
    // /
    // 2. Calcolo Lock-Free dei dati disponibili
    // Questa formula gestisce il wrap-around da sola in un ciclo di clock
    uint16_t available = (w - r) & MASK;
    
    // Quick exit if there is nothing to read / Uscita rapida se non c'è nulla da leggere
    if (available == 0) return false;

    // 3. Local cache for maximum performance / Cache locale per massime prestazioni
    uint8_t* queuePtr = _queue;

    // 4. Rapid "draining" cycle and parsing / Ciclo di "svuotamento" rapido e Parsing
    for (uint16_t i = 0; i < available; i++) {
        uint8_t dato = queuePtr[r];
        
        // The parser operates in real time. / Il parser lavora in tempo reale.
        packet.parse(dato, true); 

        // Advancement with bitwise wrap-around / Avanzamento con wrap-around bitwise
        r = (r + 1) & MASK;
    }

    // --- COMPILER BARRIER / BARRIERA DEL COMPILATORE ---
    // Essential on Xtensa/ARM architectures: Ensures the CPU physically finishes
    // transferring all bytes from cached RAM to main RAM
    // BEFORE telling the ISR (Interrupt Service Routine) that the space is free again.
    // /
    // Essenziale sulle architetture Xtensa/ARM: Assicura che la CPU finisca 
    // fisicamente di trasferire tutti i byte dalla RAM cacheata alla RAM principale
    // PRIMA di dire all'ISR (Interrupt Service Routine) che lo spazio è di nuovo libero. 
    asm volatile ("memw" : : : "memory");

    // 5. Final ATOMIC synchronization / Sincronizzazione finale ATOMICA
    // We update _reader in RAM exactly once. As soon as this line executes,
    // the RX callback mathematically determines that there is new free space for writing.
    // /
    // Aggiorniamo _reader in RAM una sola volta. Appena viene eseguita questa riga,
    // la callback RX capisce matematicamente che c'è nuovo spazio libero per scrivere.
    _reader = r;

    return true;
}

int SerialWireless_::available() {
  // Distance between writer and reader / Distanza tra chi scrive e chi legge
  return (_writerSerialRead - _readerSerialRead) & MASK_READ_SERIAL;
}

int SerialWireless_::availableBin() {
  return (writeIndex - readIndex) & (BUFFER_SIZE - 1);
}

// How many FREE bytes remain in the buffer (Standard Arduino API) / Quanti byte LIBERI rimangono nel buffer (Standard Arduino API)
int SerialWireless_::availableForWrite() {
  return (_readerSerialWrite - _writerSerialWrite - 1) & MASK_WRITE_SERIAL;
}

// helper function / funzione di comodo
// How many bytes are currently in the buffer waiting to be sent
// /
// Quanti byte sono attualmente nel buffer in attesa di essere inviati
int SerialWireless_::availableBufferSerialWrite() {
  return (_writerSerialWrite - _readerSerialWrite) & MASK_WRITE_SERIAL;
}

int SerialWireless_::availableForWriteBin() {
  return (readIndex - writeIndex - 1) & (BUFFER_SIZE - 1);
}

// ===================================================================================
// DATA TRANSMISSION: PRODUCER-CONSUMER PATTERN / TRASMISSIONE DATI: PRODUCER-CONSUMER PATTERN
// ===================================================================================

// In a hardware serial interface, `flush` blocks the thread until the last bit leaves the TX pin.
// Here, however, `flush` packages the data into the ring buffer and "wakes up" the radio task.
// If the network buffer is saturated, we use `taskYIELD` to voluntarily yield control
// to the operating system (FreeRTOS), allowing the radio task to transmit and free up space,
// thereby preventing Watchdog Timeouts and keeping the system responsive.
// /
// In una seriale hardware, `flush` blocca il thread finché l'ultimo bit non lascia il pin TX. 
// Qui invece il flush impacchetta i dati nel ring buffer e "sveglia" il task radio. 
// Se il buffer di rete è saturo, usiamo taskYIELD per cedere volontariamente il controllo 
// al sistema operativo (FreeRTOS) permettendo al task radio di trasmettere e liberare spazio, 
// prevenendo così i Watchdog Timeout e mantenendo il sistema responsivo.
void SerialWireless_::flush() {
  // The flush must wait as long as there is data in the serial write buffer. / Il flush deve attendere finché ci sono dati nel buffer di scrittura seriale
  while (availableBufferSerialWrite() > 0) {
    
    // We call flush_sem(). 
    // There is no need to acquire the mutex here because flush_sem()
    // internally uses the spinlock (portENTER_CRITICAL) to protect the pointers.
    // /
    // Chiamiamo flush_sem(). 
    // Non serve prendere il Mutex qui perché flush_sem() al suo interno 
    // usa già lo Spinlock (portENTER_CRITICAL) per proteggere i puntatori.
    if (!flush_sem()) {
      // If flush_sem() returns false, it means the RADIO buffer is full. 
      // We need to give radioTask time to transmit in order to free up space.
      // /
      // Se flush_sem() restituisce false, significa che il buffer RADIO è pieno.
      // Dobbiamo dare tempo al radioTask di trasmettere per liberare spazio.
      taskYIELD(); 
    }
    
    // A brief pause to avoid "hogging" the core if the buffer is very large
    // /
    // Un piccolo respiro per non "sequestrare" il Core se il buffer è molto grande
    taskYIELD(); 
  }
}

// Preparation phase: Retrieve pending serial data, apply COBS stuffing,
// encapsulate the result in a packet with an ID and CRC, and place it in the buffer,
// ready for transmission over the air.
// /
// Fase di preparazione: Preleva i dati seriali in attesa, applica lo stuffing COBS, 
// incapsula il tutto in un pacchetto con ID e CRC, e lo deposita nel buffer pronto 
// per essere spedito in etere.
bool SerialWireless_::flush_sem() {
  bool packed_data = false;
  
  // 1. Serial Lock / Lucchetto Seriale
  portENTER_CRITICAL(&mux_serial_tx);

  uint16_t w = _writerSerialWrite;
  uint16_t r = _readerSerialWrite;
  uint16_t available_tx = (w - r) & MASK_WRITE_SERIAL;

  if (available_tx > 0) {
    // 2. Lock the radio before calculating the space! / LUCCHETTO RADIO PRIMA DI CONTARE LO SPAZIO!
    // Now no one can steal the space from us while we calculate the packet. / Ora nessuno può rubarci lo spazio mentre calcoliamo il pacchetto
    portENTER_CRITICAL(&mux_radio_tx);

    uint16_t space_in_radio = availableForWriteBin();
    const uint16_t overhead = PREAMBLE_SIZE + POSTAMBLE_SIZE;

    if (space_in_radio > overhead) {
      uint16_t max_payload = space_in_radio - overhead;
      
      // Dynamic fragmentation strategy:
      // If the data exceeds the physical limit of a single ESP-NOW packet or the payload
      // capacity (200 bytes for safety), we split the transmission into smaller chunks.
      // /
      // Strategia di frammentazione dinamica:
      // Se i dati superano il limite fisico di un singolo pacchetto ESP-NOW o la capacità 
      // del payload (200 byte per sicurezza), spezziamo l'invio in chunk più piccoli.
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
    
    // Release the radio / Rilasciamo la radio
    portEXIT_CRITICAL(&mux_radio_tx);
  }

  // Release the serial port / Rilasciamo la seriale
  portEXIT_CRITICAL(&mux_serial_tx); 

  // System API / API di sistema
  // Signal the asynchronous Radio Task that there are new packets ready for transmission
  // /
  // Segnaliamo al Task Radio asincrono che ci sono nuovi pacchetti pronti al decollo.
  if (packed_data) {
    esp_timer_stop(timer_handle_serial);
    xTaskNotifyGive(xRadioTaskHandle);
    _serial_needs_recovery = false;
    return true;
  }

  return false; 
}

// Launch phase: This function physically "pushes" the bits into the ESP-NOW radio hardware.
// /
// Fase di lancio: Questa funzione "spinge" fisicamente i bit nell'hardware radio ESP-NOW.
void SerialWireless_::SendData_sem() {
  // Snapshot locale degli indici
  const uint16_t r = readIndex;
  const uint16_t w = writeIndex; // Let's look at where the producer has reached at this moment. / Guardiamo dove è arrivato il produttore in questo istante
  const uint16_t MASK = BUFFER_SIZE - 1;

  // Lock-free calculation of available data
  // This formula automatically handles wrap-around
  // /
  // Calcolo Lock-Free dei dati disponibili
  // Questa formula gestisce il wrap-around automaticamente
  const uint16_t available_data = (w - r) & MASK;

  // Early exit if there is nothing to send
  // /
  // Uscita rapida se non c'è nulla da spedire
  if (available_data == 0) return;

  uint16_t len_tx = available_data;

  // PACKET ALIGNMENT ("Cut" prevention)
  // ESP-NOW has a hard hardware limit of 250 bytes per payload. 
  // Truncating a frame mid-stream breaks the COBS protocol on the receiving end. 
  // Therefore, we scan the buffer backwards for a terminator byte (0xFF)
  // and align the cut to ensure that logically intact packets are sent.
  // /
  // ALLINEAMENTO PACCHETTI (Prevenzione "taglio")
  // ESP-NOW ha un limite hardware rigido di 250 byte per payload.
  // Se tronchiamo un frame a metà distruggiamo il protocollo COBS in ricezione.
  // Pertanto, analizziamo il buffer a ritroso cercando un byte terminatore (0xFF) 
  // e allineiamo il taglio garantendo l'invio di pacchetti logicamente integri.
  if (available_data > ESP_NOW_MAX_DATA_LEN) {
    uint16_t aligned_len = 0;
    
    // Backward scan from the physical limit. / Scansione a ritroso dal limite fisico.
    for (uint16_t i = ESP_NOW_MAX_DATA_LEN; i > 0; i--) {
      if (buffer[(r + i) & MASK] == START_BYTE) {
        aligned_len = i; 
        break;
      }
    }

    len_tx = (aligned_len > 0) ? aligned_len : ESP_NOW_MAX_DATA_LEN;
  }

  // BATCH COPY / COPIA BATCH
  const uint16_t firstChunk = BUFFER_SIZE - r;
  
  if (firstChunk < len_tx) {
    const uint16_t secondChunk = len_tx - firstChunk;
    memcpy(buffer_espnow, &buffer[r], firstChunk);
    memcpy(&buffer_espnow[firstChunk], &buffer[0], secondChunk);
  } else {
    memcpy(buffer_espnow, &buffer[r], len_tx);
  }

// SUBMISSION AND UPDATING / INVIO E AGGIORNAMENTO
  esp_err_t result = esp_now_send(peerAddress, buffer_espnow, len_tx);
  
  if (result == ESP_OK) {
    radioFree = false; 
    
    // ATOMIC Update / Aggiornamento ATOMICO
    readIndex = (r + len_tx) & MASK; 
  } else {
    // THE HEALTHY FREEZE:
    // Transmission failed (e.g., peer powered off/interference). 
    // We do not advance readIndex to avoid losing vital serial data. 
    // However, we pause the system for 1 ms to prevent the Watchdog from triggering
    // and to avoid a forced reboot of the Lightgun.
    // /
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
  bool start_timer = false; // Security flag / Flag di sicurezza

  while (written < len) {
    start_timer = false; // We reset at every cycle. / Resettiamo a ogni ciclo
    
    portENTER_CRITICAL(&mux_serial_tx); // --- LOCK START / INIZIO LOCK ---

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
        start_timer = true; // Let's raise the flag instead of calling the API. / Alziamo la flag invece di chiamare l'API
      }

      asm volatile ("memw" : : : "memory");
      _writerSerialWrite = (w + toWrite) & MASK_WRITE_SERIAL;
      
      portEXIT_CRITICAL(&mux_serial_tx); // --- END LOCK / FINE LOCK ---

      // We are safe now: start the timer if necessary! 
      // We use a hardware timer to ensure that even small payloads (which do not
      // saturate the buffer) are regularly flushed without getting stuck.
      // /
      // Ora siamo al sicuro: avviamo il timer se necessario!
      // Usiamo un timer hardware per garantire che anche piccoli payload (che non 
      // saturano il buffer) vengano flushati regolarmente senza rimanere incastrati.
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
  // Use the same lock as SendPacket for consistency and speed
  // A hardware mutex momentarily disables interrupts: use it
  // for the absolute minimum time necessary.
  // /
  // Usa lo stesso lucchetto di SendPacket per coerenza e velocità
  // Un mutex hardware disabilita gli interrupt momentaneamente: usiamolo
  // per il minimo tempo indispensabile.
  portENTER_CRITICAL(&mux_radio_tx);

  uint16_t w = writeIndex;
  uint16_t r = readIndex;
  const uint16_t MASK = BUFFER_SIZE - 1;

  // Space calculation (Minus-One Rule) / Calcolo dello spazio (Regola del Meno Uno)
  uint16_t free_space = (r - w - 1) & MASK;

  if (free_space < len) {
    _overflow_write = true;
    portEXIT_CRITICAL(&mux_radio_tx); 
    return 0; 
  }

  // BATCH COPY / COPIA BATCH (Wrap-around logic)
  uint16_t firstChunk = BUFFER_SIZE - w;

  if (firstChunk < len) {
    memcpy(&buffer[w], data, firstChunk);
    memcpy(&buffer[0], data + firstChunk, len - firstChunk);
  } else {
    memcpy(&buffer[w], data, len);
  }

  // MEMORY BARRIER - Essential for ESP32/Xtensa / BARRIERA DI MEMORIA - Fondamentale per ESP32/Xtensa
  asm volatile ("memw" : : : "memory");

  // Index Update - Atomic / Aggiornamento dell'indice - Atomico
  writeIndex = (w + len) & MASK; 

  portEXIT_CRITICAL(&mux_radio_tx);

  return len;
}

void SerialWireless_::SendPacket(const uint8_t *data, uint8_t len, uint8_t packetID) { 
  const uint16_t totalLen = len + PREAMBLE_SIZE + POSTAMBLE_SIZE;

  // Enter critical section (hardware spinlock)
  // Extremely fast; protects the packet.txBuff "workspace" from other cores or tasks
  // /
  // Entriamo in sezione critica (Spinlock hardware)
  // Velocissimo, protegge il "tavolo da lavoro" packet.txBuff da altri Core o Task
  portENTER_CRITICAL(&mux_radio_tx);

  // Lock-free snapshot of available space / Snapshot Lock-Free dello spazio disponibile
  uint16_t w = writeIndex;
  uint16_t r = readIndex;
  const uint16_t MASK = BUFFER_SIZE - 1;
  uint16_t free_space = (r - w - 1) & MASK; 

  // If there is space, we build and write right away. / Se c'è spazio, costruiamo e scriviamo subito
  if (free_space >= totalLen) {
    uint8_t* pDest = packet.txBuff;

    // We copy the data and build the packet in the temporary buffer. / Copiamo i dati e costruiamo il pacchetto nel buffer temporaneo
    memcpy(&pDest[PREAMBLE_SIZE], data, len);
    packet.constructPacket(len, packetID);
    
    // Insert into the main Ring Buffer
    // (Since we are already in a critical section, writeBin will be lightning-fast)
    // /
    // Inseriamo nel Ring Buffer principale
    // (Dato che siamo già in sezione critica, writeBin sarà un fulmine)
    writeBin(pDest, totalLen);
    
    // Exit the critical section before notifying
    // /
    // Usciamo dalla sezione critica prima di notificare
    portEXIT_CRITICAL(&mux_radio_tx);

    // Wake up the radioTask / Svegliamo il radioTask
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
  
  // Calculation of free space (Minus-One Rule) / Calcolo spazio libero (Regola del Meno Uno)
  uint16_t free_space = (r - w - 1) & MASK_READ_SERIAL;

  if (free_space >= len) {
    uint16_t firstChunk = FIFO_SIZE_READ_SERIAL - w;
    if (firstChunk < len) {
      memcpy(&bufferSerialRead[w], data, firstChunk);
      memcpy(&bufferSerialRead[0], data + firstChunk, len - firstChunk);
    } else {
      memcpy(&bufferSerialRead[w], data, len);
    }
    
    // Memory barrier before updating the pointer / Barriera di memoria prima di aggiornare il puntatore
    asm volatile ("memw" : : : "memory");
    _writerSerialRead = (w + len) & MASK_READ_SERIAL;
  } else {
    _overflow_bufferSerialRead = true;
  }
}

bool SerialWireless_::end() {

  vTaskDelete(xRadioTaskHandle);
  // vTaskDelete(xUSBTaskHandle); // This only goes into the dongle. / questo va solo nel dongle
  
  esp_err_t err = esp_now_deinit();
  if (err != ESP_OK) {
    //Serial.printf("esp_now_deinit failed! 0x%x", err);
    return false;
  }
  WiFi.disconnect(true);
  WiFi.mode(WIFI_OFF); // Turn off Wi-Fi / Disattiva il WiFi
  esp_timer_delete(timer_handle_serial);
  return true;
}


void SerialWireless_::init_wireless() {
  configST myConfig; // configuration utility variable / variabile di utilità per configurazione
  
  myConfig.timeout      = DEFAULT_TIMEOUT; // 50ms
  myConfig.callbacks    = callbackArr;
  myConfig.callbacksLen = sizeof(callbackArr) / sizeof(functionPtr);
  packet.begin(myConfig);

  setupTimerSerial(); // Create the timers for sending serial data. / crea i timer .. timer per invio dati seriali

  // ========= ADDITION FOR NEW RADIO TRANSMISSION/RECEPTION MANAGEMENT FOR ALL / AGGIUNTA PER NUOVA GESTIONE TRASMISSIONE RICEZIONE RADIO X TUTTI===============
  xTaskCreatePinnedToCore(
        radioTask,          // task function / funzione del task
        "radioTask",        // name / nome
        4096,              // stack size
        NULL,              // parameters / parametri
        15,                 // HIGH priority 1 to 24? / priorità ALTA DA 1 A 24 ?
        &xRadioTaskHandle,   // handle
        APP_CPU_NUM        // core (you can use 0 or 1) / (puoi usare 0 o 1)
      );
  // ========== END ADDITION FOR NEW RADIO TRANSMISSION/RECEPTION MANAGEMENT FOR ALL / FINE AGGIUNTA PER NUOVA GESTIONE =================

  #ifdef DONGLE
  // ========= Addition for new USB dongle management (applies only to the dongle) / AGGIUNTA PER NUOVA GESTIONE DONGLE USB (SERVE SOLO AL DONGLE) ===============
  xTaskCreatePinnedToCore(
        usbTask,          // task function / funzione del task
        "usbTask",        // name / nome
        4096,              // stack size
        NULL,              // parameters / parametri
        10,                 // HIGH priority 1 to 24? / priorità ALTA DA 1 A 24 ?
        &xUSBTaskHandle,   // handle
        APP_CPU_NUM        // core (you can use 0 or 1) / (puoi usare 0 o 1) APP_CPU_NUM = 1 (dove gira il loop) PRO_CPU_NUM = 0 (dove gira freertos, wifi, ecc.)
      );
  // ========== END Addition for new USB dongle management (applies only to the dongle) / FINE AGGIUNTA PER NUOVA GESTIONE =================
  #endif // DONGLE

}


void SerialWireless_::begin() {

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
  WiFi.disconnect();
  

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

  err = esp_wifi_set_max_tx_power(espnow_wifi_power); // between 8 and 84 (corresponding to 2 dBm to 20 dBm); / tra 8 e 84 corrispondenti a 2dbm a 20 dbm);
  if (err != ESP_OK) {
    //Serial.printf("esp_wifi_set_max_tx_power failed! 0x%x", err);
  }

  err = esp_wifi_set_ps(WIFI_PS_NONE);  // Disable power saving for maximum responsiveness. / DISABILITA IL POWER SAVE PER MASSIMA REATTIVITà
  if (err != ESP_OK) {
    //Serial.printf("esp_wifi_set_ps failed! 0x%x", err);
  }
  
  vTaskDelay(pdMS_TO_TICKS(500)); // delay(1000);
    
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
    //Serial.println("Error adding peer");
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

  TinyUSBDevices.wireless_mode = WIRELESS_MODE::ENABLE_ESP_NOW_TO_DONGLE;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

// ===================================================================================
// HANDSHAKE PROTOCOL: CONNECTION STATE MACHINE / MACCHINA A STATI DI CONNESSIONE
// ===================================================================================
// This implements the complex pairing dance between the devices.
// The logic is asymmetric:
// - The GUN acts as a "Beacon" (Master) and broadcasts packets across various channels.
// - The DONGLE and PEDAL act as "Listeners" (Slaves) in Promiscuous Mode.
// Once discovered, they exchange MAC addresses (to switch from insecure broadcast
// to point-to-point Unicast communication) and USB information (Spoofing).
// /
// Questa è l'implementazione della complessa danza di accoppiamento tra i dispositivi.
// La logica è asimmetrica: 
// - La GUN fa da "Faro" (Master) e spara pacchetti in broadcast sui vari canali.
// - Il DONGLE e il PEDAL fanno da "Ascoltatori" (Slave) in Promiscuous Mode.
// Una volta trovati, si scambiano i MAC Address (per passare dall'insicuro broadcast 
// a una comunicazione Unicast punto-a-punto) e le informazioni USB (Spoofing).

// ============ NUOVA IMPLEMNTAZIONE DONGLE == LA GUN FA IL FARO ================================
#ifdef DONGLE
bool SerialWireless_::connection_dongle() {
    
  #if defined(DONGLE) && defined(USES_DISPLAY)
  TaskHandle_t animTaskHandleLink = NULL;  
  if(display_init) {
  // Avvio animazione
    if (animTaskHandleLink == NULL) {
      xTaskCreatePinnedToCore(
        animTaskLink,          // task function / funzione del task
        "AnimTaskLink",        // name / nome
        4096,              // stack size
        NULL,              // parameters / parametri
        1,                 // priority / priorità
        &animTaskHandleLink,   // handle
        APP_CPU_NUM        // core (you can use 0 or 1) / (puoi usare 0 o 1)
      );
    }
  }
  #endif // USES_DISPLAY

  #define TIMEOUT_DONGLE_DIALOGUE 2000 // in milliseconds / in millisecondi
  unsigned long lastMillis_start_dialogue = millis ();
  
  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
    
  // ====================================================
  
  // WAITING FOR THE BEACON / ATTESA DEL FARO
  // The dongle opens the ports to accept broadcast packets from the gun.
  // /
  // Il Dongle apre le porte per accettare i pacchetti Broadcast della Gun
  broadcast_receiver = true;

  while(stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) { 
    if (stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
      lastMillis_start_dialogue = millis ();
    }
    // Security Timeout: If the Gun disappears mid-handshake, we reset the state. / Timeout di Sicurezza: Se la Gun sparisce a metà handshake, resettiamo lo stato
    if (((millis() - lastMillis_start_dialogue) > TIMEOUT_DONGLE_DIALOGUE) && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
      stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
    }
    taskYIELD();
  }
  
  // UNICAST LOCK-IN
  // Connection established. Closing broadcast ports for security and to
  // reduce CPU load from processing stray packets.
  // /
  // Connessione stabilita. Chiudiamo le porte al broadcast per sicurezza e per 
  // alleggerire il carico CPU sulla ricezione di pacchetti vaganti.
  broadcast_receiver = false;

  //Serial.println("DONGLE - Negotiation complete - pairing GUN/DONGLE devices");
  //Serial.println("DONGLE - Negosazione completata - associazione dei dispositivi GUN/DONGLE");
  if (esp_now_del_peer(peerAddress) != ESP_OK) {  // remove the broadcast from the peers / cancella il broadcast dai peer
    //Serial.println("DONGLE - Error deleting broadcast peer");
    //Serial.println("DONGLE - Errore nella cancellazione del peer broadcast");
  }
  
  // We register the Gun as the only valid peer for communication. / Registriamo la Gun come unico peer valido per la comunicazione
  memcpy(peerAddress, mac_esp_another_card, 6);
  memcpy(peerInfo.peer_addr, peerAddress, 6);
  if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // adds the dongle to the peers / inserisce il dongle nei peer
    //Serial.println("DONGLE - Error adding the new GUN peer");
    //Serial.println("DONGLE - Errore nell'aggiunta del nuovo peer della GUN");
  } else esp_now_set_peer_rate_config(peerAddress, &rate_config);

  
  TinyUSBDevices.onBattery = true;
  
    #if defined(DONGLE) && defined(USES_DISPLAY)
      if (animTaskHandleLink != NULL) {
        vTaskDelete(animTaskHandleLink);
        animTaskHandleLink = NULL;
      }
    #endif // USES_DISPLAY
  
  return true; 
}
#endif // DONGLE
// ============ END OF NEW DONGLE IMPLEMENTATION == LAGANN ACTS AS A BEACON / FINE NUOVA IMPLEMNTAZIONE DONGLE == LA GUN FA IL FARO ================================

// ===================================================================================
// FAST RECOVERY / RICONNESSIONE RAPIDA (FAST RECOVERY)
// ===================================================================================
// Functions used to resume communication without repeating the entire
// slow handshake procedure. If the devices already know each other's
// MAC addresses (e.g., saved in EEPROM/Flash in the future), they "knock"
// directly on the specific port.
// /
// Funzioni utilizzate per riprendere la comunicazione senza dover rifare tutta
// la lenta procedura di Handshake. Se i dispositivi conoscono già i rispettivi 
// MAC Address (es. salvati in EEPROM/Flash in futuro), "bussano" direttamente 
// alla porta specifica.

#ifdef GUN
bool SerialWireless_::connection_gun_at_last_dongle() {
  #define TIMEOUT_TX_PACKET_LAST_DONGLE 300 // in milliseconds - Packet transmission interval in milliseconds, resulting in 4–5 packets./ in millisecondi - tempo di invio pacchetti ogni millisecondi quindi 4-5 pacchetti
  #define TIMEOUT_DIALOGUE_LAST_DONGLE 2000 // in milliseconds - Maximum time to search for the last dongle / in millisecondi - tempo massimo per ricerca ultimo dongle
  unsigned long lastMillis_tx_packet_last_dongle = 0;
  unsigned long lastMillis_start_dialogue_last_dongle = millis ();

  uint8_t aux_buffer_tx[13];
  
  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION_DONGLE;
  aux_buffer_tx[0] = CONNECTION_STATE::TX_CHECK_CONNECTION_LAST_DONGLE;
  memcpy(&aux_buffer_tx[1], SerialWireless.mac_esp_inteface, 6);
  memcpy(&aux_buffer_tx[7], peerAddress, 6);
  while (!TinyUSBDevice.mounted() && 
         (stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED_WITH_LAST_DONGLE) &&
         ((millis() - lastMillis_start_dialogue_last_dongle) < TIMEOUT_DIALOGUE_LAST_DONGLE)) { 
      if ((millis() - lastMillis_tx_packet_last_dongle) > TIMEOUT_TX_PACKET_LAST_DONGLE)
      {
        SerialWireless.SendPacket((const uint8_t *)aux_buffer_tx, 13, PACKET_TX::CHECK_CONNECTION_LAST_DONGLE); 
        lastMillis_tx_packet_last_dongle = millis();
      }
      
      taskYIELD();
  }
  if (stato_connessione_wireless == CONNECTION_STATE::DEVICES_CONNECTED_WITH_LAST_DONGLE) {    
    stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED;
    TinyUSBDevices.onBattery = true;
    return true;
  } else {
    // Fast Recovery failure: preparing for a Full Handshake. / Fallimento Fast Recovery: prepariamo il terreno per il Full Handshake.
    stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
    lastDongleSave=false;
    esp_now_deinit();
    begin();
    return false;
  }
}
#endif // GUN

#ifdef GUN
bool SerialWireless_::connection_gun_at_last_pedal() {
  #define TIMEOUT_TX_PACKET_LAST_PEDAL 300 // in milliseconds - Packet transmission interval in milliseconds, resulting in 4–5 packets./ in millisecondi - tempo di invio pacchetti ogni millisecondi quindi 4-5 pacchetti
  #define TIMEOUT_DIALOGUE_LAST_PEDAL 2000 // in milliseconds - Maximum time to search for the last dongle / in millisecondi - tempo massimo per ricerca ultimo pedal
  unsigned long lastMillis_tx_packet_last_pedal = 0;
  unsigned long lastMillis_start_dialogue_last_pedal = millis ();
  uint8_t aux_buffer_tx[13];
  
  // =============================================
  uint8_t peerAddress_copy[6]; 
  memcpy(peerAddress_copy, peerAddress, 6);

  memcpy(peerAddress, lastPedalAddress, 6);        
  memcpy(peerInfo.peer_addr, peerAddress, 6);
  peerInfo.channel = espnow_wifi_channel;
  if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserts the dongle into the peers / inserisce il dongle nei peer
    //Serial.println("Error adding the new peer");
    //Serial.println("Errore nell'aggiunta del nuovo peer");
  } else esp_now_set_peer_rate_config(peerAddress, &rate_config);         
  // =============================================

  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION_PEDAL;
  aux_buffer_tx[0] = CONNECTION_STATE::TX_CHECK_CONNECTION_LAST_PEDAL;
  memcpy(&aux_buffer_tx[1], SerialWireless.mac_esp_inteface, 6);
  memcpy(&aux_buffer_tx[7], peerAddress, 6);
  while ((stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED_WITH_LAST_PEDAL) &&
         ((millis() - lastMillis_start_dialogue_last_pedal) < TIMEOUT_DIALOGUE_LAST_PEDAL)) { 
      if ((millis() - lastMillis_tx_packet_last_pedal) > TIMEOUT_TX_PACKET_LAST_PEDAL)
      {
        SerialWireless.SendPacket((const uint8_t *)aux_buffer_tx, 13, PACKET_TX::CHECK_CONNECTION_LAST_PEDAL); 
        lastMillis_tx_packet_last_pedal = millis();
      }
      
      taskYIELD();
  }
  
  if (esp_now_del_peer(peerAddress) != ESP_OK) {  // remove the broadcast from the peers / cancella il broadcast dai peer
      //Serial.println("DONGLE - Error removing the peer broadcast");
      //Serial.println("DONGLE - Errore nella cancellazione del peer broadcast");
  }

  memcpy(peerAddress, peerAddress_copy, 6);

  if (stato_connessione_wireless == CONNECTION_STATE::DEVICES_CONNECTED_WITH_LAST_PEDAL) {    
    //stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED;
    
    memcpy(peerAddress_pedal, lastPedalAddress, 6);
    memcpy(peerInfo.peer_addr, peerAddress_pedal, 6);
    peerInfo.channel = espnow_wifi_channel;
    if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserts the dongle into the peers / inserisce il dongle nei peer
      //Serial.println("DONGLE - Error adding the new GUN peer");
      //Serial.println("DONGLE - Errore nell'aggiunta del nuovo peer della GUN");
    } else esp_now_set_peer_rate_config(peerAddress_pedal, &rate_config);

    setupTimerPedal();  // Create the timer to safely manage the pedal / CREA IL TIMER PER GESTIRE IN SICUREZZA IL PEDA

    
    TinyUSBDevices.is_pedal_wireless = true;

    return true;
  } else {
    //stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
    TinyUSBDevices.is_pedal_wireless = false;
    
    lastPedalSave=false;
    return false;
  }
}
#endif // GUN

// ===================================================================================
// MASTER SIDE (GUN): COMPLETE HANDSHAKE PROCEDURE
// /
// LATO MASTER (GUN): PROCEDURA COMPLETA DI HANDSHAKE
// ===================================================================================

// ======================= IMPLEMENTATION WHERE THE GUN IS THE BEACON. / NUOVA IMPLEMENTAZIONE DOVE LA GUN FA IL FARO ===============
#ifdef GUN
bool SerialWireless_::connection_gun() {
  
  channel_display = espnow_wifi_channel;

  #if defined(GUN) && defined(USES_DISPLAY)
  TaskHandle_t animTaskHandleLink = NULL;  
  if(display_OLED != nullptr) {
  // Avvio animazione
    if (animTaskHandleLink == NULL) {
      xTaskCreatePinnedToCore(
        animTaskLink,          // task function / funzione del task
        "AnimTaskLink",        // name / nome
        4096,              // stack size
        NULL,              // parameters / parametri
        1,                 // priority / priorità
        &animTaskHandleLink,   // handle
        APP_CPU_NUM        // core (you can use 0 or 1) / (puoi usare 0 o 1)
      );
    }
  }
  #endif // USES_DISPLAY
  
  
  uint8_t channel = espnow_wifi_channel;   // between 1 and 13 (14 should be reserved) / tra e 1 e 13 (il 14 dovrebbe essere riservato)
  #define TIMEOUT_GUN_TX_PACKET 80  // 500 // in milliseconds / in millisecondi
  #define TIMEOUT_GUN_CHANGE_CHANNEL 250  // 2000 // in milliseconds - changes channel every / in millisecondi - cambia canale ogni
  #define TIMEOUT_GUN_DIALOGUE 3000  //  6000 // in milliseconds - maximum time to complete the pairing operation / in millisecondi - tempo massimo per completare operazione accoppiamento
  unsigned long lastMillis_tx_packet = millis ();
  unsigned long lastMillis_change_channel = millis ();
  unsigned long lastMillis_start_dialogue = millis ();
  uint8_t aux_buffer_tx[14]; // added a byte to also transmit the transmission channel 
                             // when transmitting the presence advertisement packet, the channel is also included
                             // /
                             // aggiunto un byte per trasmettere anche il canale di trasmissione
                             // durante tx pacchetto pubblicazione presenza, mette anche il canale
  
 

  // ===== IMPOSTARE IL PEER BOARCAST ===============
  if (esp_now_del_peer(peerAddress) != ESP_OK) {  // remove the broadcast from the peers / cancella il broadcast dai peer
    //Serial.println("Error deleting peer");
    //Serial.println("Errore nella cancellazione del peer");
  }
  memcpy(peerAddress, BROADCAST_ADDR, 6);
  memcpy(peerInfo.peer_addr, peerAddress, 6);
  if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserts the dongle into the peers / inserisce il dongle nei peer
    //Serial.println("Error adding the new peer");
    //Serial.println("Errore nell'aggiunta del nuovo peer");
  } else esp_now_set_peer_rate_config(peerAddress, &rate_config);                       
  // ====================================================
  
  stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
  aux_buffer_tx[0] = CONNECTION_STATE::TX_GUN_SEARCH_DONGLE_BROADCAST;
  memcpy(&aux_buffer_tx[1], SerialWireless.mac_esp_inteface, 6);
  memcpy(&aux_buffer_tx[7], peerAddress, 6);
 
  broadcast_receiver = true;

  // SEARCH LOOP (CHANNEL HOPPING) / LOOP DI RICERCA (CHANNEL HOPPING)
  // The Gun physically changes the WiFi channel every TIMEOUT_GUN_CHANGE_CHANNEL (250ms)
  // and transmits beacon packets (80ms). If it receives a response, it transitions to the next state.
  // /
  // La Gun cambia fisicamente il canale WiFi ogni TIMEOUT_GUN_CHANGE_CHANNEL (250ms)
  // e spara pacchetti beacon (80ms). Se riceve risposta passa allo stato successivo.
  while (!TinyUSBDevice.mounted() && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
    if (stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
      if (((millis() - lastMillis_change_channel) > TIMEOUT_GUN_CHANGE_CHANNEL) && 
         ((millis() - lastMillis_tx_packet)) >= (TIMEOUT_GUN_TX_PACKET - 50)) {  // Added a 50 ms margin to prevent a random channel switch immediately after sending a packet, which would otherwise leave no time for a response. / aggiunta impostato 50 ms come margine, per evitare che quando invia pacchetto cambi subito casnale senza dare possibilità risposta
        aux_buffer_tx[13] = channel;
        channel_display = channel;      
        
        esp_wifi_set_promiscuous(true);
        if (esp_wifi_set_channel(channel, WIFI_SECOND_CHAN_NONE) != ESP_OK) {
          //Serial.printf("DONGLE - esp_wifi_set_channel failed!");
        }
        esp_wifi_set_promiscuous(false);
        peerInfo.channel = channel;
        if (esp_now_mod_peer(&peerInfo) != ESP_OK) {  // modify the peer's channel / modifica il canale del peer
          //Serial.println("DONGLE - Error changing channel");
          //Serial.println("DONGLE - Errore nella modifica del canale");
        }             
        lastMillis_change_channel = millis ();
        lastMillis_tx_packet = 0;
        
        channel++;
        if (channel >13) channel = 1;// to send a packet immediately on the new channel / per fare inviare subito un pacchetto sul nuovo canale
      }
      if ((millis() - lastMillis_tx_packet) > TIMEOUT_GUN_TX_PACKET) {
        SerialWireless.SendPacket((const uint8_t *)aux_buffer_tx, 14, PACKET_TX::CONNECTION); // added a byte to also transmit the transmission channel / aggiunto un byte per trasmettere anche il canale di trasmissione
        //Serial.print("DONGLE - broadcast packet sent on channel: ");
        //Serial.print("DONGLE - inviato pacchetto broadcast sul canale: ");
        //Serial.println(channel);
        lastMillis_tx_packet = millis (); 
      }
      lastMillis_start_dialogue = millis();
    }
    else {
      // Timeout if the handshake stalls midway / Timeout se l'handshake si blocca a metà
      if (((millis() - lastMillis_start_dialogue) > TIMEOUT_GUN_DIALOGUE) && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
        stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
        //Serial.println("DONGLE - Negotiation between DONGLE/GUN did not complete; restarting from the beginning");
        //Serial.println("DONGLE - Non si è conclusa la negoziazione tra DONGLE/GUN e si riparte da capo");
        lastMillis_change_channel = millis ();
      }  
    }
    taskYIELD();
  }

  broadcast_receiver = false;  

  if (stato_connessione_wireless == CONNECTION_STATE::DEVICES_CONNECTED) {
    //Serial.println("DONGLE - Negotiation complete - pairing GUN/DONGLE devices");
    //Serial.println("DONGLE - Negosazione completata - associazione dei dispositivi GUN/DONGLE");
    if (esp_now_del_peer(peerAddress) != ESP_OK) {  // remove the broadcast from the peers / cancella il broadcast dai peer
      //Serial.println("DONGLE - Error deleting broadcast peer");
      //Serial.println("DONGLE - Errore nella cancellazione del peer broadcast");
    }
    memcpy(peerAddress, mac_esp_another_card, 6);
    memcpy(peerInfo.peer_addr, peerAddress, 6);
    espnow_wifi_channel=usb_data_wireless.channel;
    esp_wifi_set_promiscuous(true);
    esp_wifi_set_channel(espnow_wifi_channel, WIFI_SECOND_CHAN_NONE);
    esp_wifi_set_promiscuous(false);
    peerInfo.channel = espnow_wifi_channel;
    if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserts the dongle into the peers / inserisce il dongle nei peer
      //Serial.println("DONGLE - Error adding the new GUN peer");
      //Serial.println("DONGLE - Errore nell'aggiunta del nuovo peer della GUN");
    } else esp_now_set_peer_rate_config(peerAddress, &rate_config);



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
#endif // GUN
// ======================= END IMPLEMENTATION WHERE THE GUN IS THE BEACON. / FINE NUOVA IMPLEMENTAZIONE DOVE LA GUN FA IL FARO ===============
// ============ Packet transmission to the pedal connection dongle completed. / INVIO PACCHETTO AL DONGLE DI CONNESSIONE AL PEDAL COMPLETATA =============================
#ifdef GUN
void SerialWireless_::tx_gun_at_dongle_pedal_ready() {
  for (uint8_t i=0; i<3; i++) {
    SerialWireless.SendPacket((const uint8_t *)&TinyUSBDevices.is_pedal_wireless, 1, PACKET_TX::PEDAL_TX);
    vTaskDelay(pdMS_TO_TICKS(69));
  }
}
#endif // GUN
// ============ END Packet transmission to the pedal connection dongle completed. / INVIO PACCHETTO AL DONGLE DI CONNESSIONE AL PEDAL COMPLETATA =============================
// ======================= NEW IMPLEMENTATION WHERE THE GUN IS THE BEACON to connect to the pedal / NUOVA IMPLEMENTAZIONE DOVE LA GUN FA IL FARO per connetersi al pedal ===============
#ifdef GUN
bool SerialWireless_::connection_gun_at_pedal() {
  
  uint8_t seconds = MAX_SECONDI_CONNESSIONE_PEDAL;  //10; 
  seconds_display = seconds;

  channel_display = espnow_wifi_channel;

  #if defined(GUN) && defined(USES_DISPLAY)
  TaskHandle_t animTaskHandleLink_pedal = NULL;  
  if(display_OLED != nullptr) {
  // Avvio animazione
    if (animTaskHandleLink_pedal == NULL) {
      xTaskCreatePinnedToCore(
        animTaskLink_pedal,          // task function / funzione del task
        "AnimTaskLink_pedal",        // name / nome
        4096,              // stack size
        NULL,              // parameters / parametri
        1,                 // priority / priorità
        &animTaskHandleLink_pedal,   // handle
        APP_CPU_NUM        // core (you can use 0 or 1) / (puoi usare 0 o 1)
      );
    }
  }
  #endif // USES_DISPLAY
  
  
  #define TIMEOUT_GUN_AT_PEDAL_TX_PACKET 50 //500 // in milliseconds / in millisecondi
  #define TIMEOUT_GUN_AT_PEDAL_CHANGE_SECONDS 1000 // in milliseconds = 1 second / in millisecondi =  1 secondo
  #define TIMEOUT_GUN_AT_PEDAL_DIALOGUE 3000 //6000 // in milliseconds - maximum time to complete pairing operation / in millisecondi - tempo massimo per completare operazione accoppiamento
  unsigned long lastMillis_tx_packet = millis ();
  unsigned long lastMillis_change_seconds = millis ();
  unsigned long lastMillis_start_dialogue = millis ();
  uint8_t aux_buffer_tx[14]; // added a byte to also broadcast the broadcast channel 
                             // during tx presence publication package, also puts the channel
                             // /
                             // aggiunto un byte per trasmettere anche il canale di trasmissione
                             // durante tx pacchetto pubblicazione presenza, mette anche il canale
  uint8_t peerAddress_copy[6];
  
  
  memcpy(peerAddress_copy, peerAddress, 6);

  memcpy(peerAddress, BROADCAST_ADDR, 6);        
  memcpy(peerInfo.peer_addr, peerAddress, 6);
  peerInfo.channel = espnow_wifi_channel;
  if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserts the dongle into the peers / inserisce il dongle nei peer
    //Serial.println("Error adding the new peer");
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
         ((millis() - lastMillis_tx_packet)) >= (TIMEOUT_GUN_AT_PEDAL_TX_PACKET - 50)) {  // addition set 50 ms as a margin, to avoid that when sending a packet it immediately changes the channel without giving the possibility of a response / aggiunta impostato 50 ms come margine, per evitare che quando invia pacchetto cambi subito casnale senza dare possibilità risposta
        aux_buffer_tx[13] = espnow_wifi_channel;
        seconds_display = seconds;
                        
        lastMillis_change_seconds = millis ();
        lastMillis_tx_packet = 0;
        
        seconds--;
      }
      if ((millis() - lastMillis_tx_packet) > TIMEOUT_GUN_AT_PEDAL_TX_PACKET) {
        SerialWireless.SendPacket((const uint8_t *)aux_buffer_tx, 14, PACKET_TX::CONNECTION_PEDAL); // added a byte to also transmit the broadcast channel / aggiunto un byte per trasmettere anche il canale di trasmissione
        //Serial.print("GUN -> PEDAL - broadcast packet sent on channel: ");
        //Serial.print("GUN -> PEDAL - inviato pacchetto broadcast sul canale: ");
        //Serial.println(channel);
        lastMillis_tx_packet = millis (); 
      }
      lastMillis_start_dialogue = millis();
    }
    else {
      if (((millis() - lastMillis_start_dialogue) > TIMEOUT_GUN_AT_PEDAL_DIALOGUE) && stato_connessione_wireless != CONNECTION_STATE::DEVICES_CONNECTED) {
        stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION;
        //Serial.println("GUN - The negotiation between GUN/PEDAL was not concluded and we start again from scratch");
        //Serial.println("GUN - Non si è conclusa la negoziazione tra GUN/PEDAL e si riparte da capo");
        lastMillis_change_seconds = millis ();
      }  
    }
    taskYIELD();
  }
  
  broadcast_receiver = false;

  if (esp_now_del_peer(peerAddress) != ESP_OK) {  // remove the broadcast from the peers / cancella il broadcast dai peer
      //Serial.println("GUN - Error canceling peer broadcast");
      //Serial.println("GUN - Errore nella cancellazione del peer broadcast");
  }

  memcpy(peerAddress, peerAddress_copy, 6);

  // ================================
  // Force the pedal into the "alive" state by sending a fake signal
  // /
  // Forza al pedale lo stato "vivo" inviando un segnale fasullo
  SerialWireless.tx_gun_at_dongle_pedal_ready();
  /*
  for (uint8_t i=0; i<3; i++) {
    SerialWireless.SendPacket((const uint8_t *)&TinyUSBDevices.is_pedal_wireless, 1, PACKET_TX::PEDAL_TX);
    vTaskDelay(pdMS_TO_TICKS(69));
  }
  */
  // ================================

  if (stato_connessione_wireless == CONNECTION_STATE::DEVICES_CONNECTED) {
    //Serial.println("DONGLE - Negotiation completed - pairing of GUN/DONGLE devices");
    //Serial.println("DONGLE - Negosazione completata - associazione dei dispositivi GUN/DONGLE");

    memcpy(peerAddress_pedal, mac_esp_another_card, 6);
    memcpy(peerInfo.peer_addr, peerAddress_pedal, 6);
    peerInfo.channel = espnow_wifi_channel;
    if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserts the dongle into the peers / inserisce il dongle nei peer
      //Serial.println("DONGLE - Error adding new GUN peer");
      //Serial.println("DONGLE - Errore nell'aggiunta del nuovo peer della GUN");
    } else esp_now_set_peer_rate_config(peerAddress_pedal, &rate_config);

    setupTimerPedal();  // CREATE THE TIMER TO SAFELY MANAGE THE PEDAL / CREA IL TIMER PER GESTIRE IN SICUREZZA IL PEDALE



    #if defined(GUN) && defined(USES_DISPLAY)
      if (animTaskHandleLink_pedal != NULL) {
        vTaskDelete(animTaskHandleLink_pedal);
        animTaskHandleLink_pedal = NULL;
      }
    #endif // USES_DISPLAY

    TinyUSBDevices.is_pedal_wireless = true;
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
#endif // GUN
// ======================= END NEW IMPLEMENTATION WHERE THE GUN IS THE BEACON to connect to the pedal / FINE NUOVA IMPLEMENTAZIONE DOVE LA GUN FA IL FARO per connetersi al pedal ===============

// ======================= PEDAL IMPLEMENTATION - STAY LISTENING ON ALL CHANNELS / IMPLEMENTAZIONE PEDAL - RIMANE IN ASCOLTO SU TUTTI I CANALI ===============
#ifdef PEDAL
bool SerialWireless_::connection_pedal() {
  
    
  uint8_t channel = espnow_wifi_channel;   // between and 1 and 13 (14 should be reserved) / tra e 1 e 13 (il 14 dovrebbe essere riservato)
  
  #define TIMEOUT_PEDAL_CHANGE_CHANNEL 150 // 2000 // in milliseconds- change channel every / in millisecondi - cambia canale ogni
  #define TIMEOUT_PEDAL_DIALOGUE 3000 // 6000 // in milliseconds - maximum time to complete pairing operation / in millisecondi - tempo massimo per completare operazione accoppiamento
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
        
        peerInfo.channel = channel;  // you could also always set it to 0 which uses the wifi one / si potrebbe anche impostare sempre a 0 che usa quello del wifi
        if (esp_now_mod_peer(&peerInfo) != ESP_OK) {  // modify the peer's channel / modifica il canale del peer
          //Serial.println("DONGLE - Error changing channel");
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
        //Serial.println("DONGLE - Negotiation between DONGLE/GUN did not complete; restarting from the beginning");
        //Serial.println("DONGLE - Non si è conclusa la negoziazione tra DONGLE/GUN e si riparte da capo");
        lastMillis_change_channel = millis ();
      }  
    }
    taskYIELD();
  }
    
  broadcast_receiver = false;
    
    //Serial.println("DONGLE - Negotiation completed - pairing of GUN/DONGLE devices");
    //Serial.println("DONGLE - Negosazione completata - associazione dei dispositivi GUN/DONGLE");
    if (esp_now_del_peer(peerAddress) != ESP_OK) {  // remove the broadcast from the peers / cancella il broadcast dai peer
      //Serial.println("DONGLE - Error deleting peer broadcast");
      //Serial.println("DONGLE - Errore nella cancellazione del peer broadcast");
    }
    memcpy(peerAddress, mac_esp_another_card, 6);
    memcpy(peerInfo.peer_addr, peerAddress, 6);
    espnow_wifi_channel=usb_data_wireless.channel;
    esp_wifi_set_promiscuous(true);
    esp_wifi_set_channel(espnow_wifi_channel, WIFI_SECOND_CHAN_NONE);
    esp_wifi_set_promiscuous(false);
    peerInfo.channel = espnow_wifi_channel;
    if (esp_now_add_peer(&peerInfo) != ESP_OK) {  // inserts the dongle into the peers / inserisce il dongle nei peer
      //Serial.println("DONGLE - Error adding new GUN peer");
      //Serial.println("DONGLE - Errore nell'aggiunta del nuovo peer della GUN");
    } else esp_now_set_peer_rate_config(peerAddress, &rate_config);

  
    TinyUSBDevices.onBattery = true;
    return true;  
      
}
#endif // PEDAL
// ======================= END PEDAL IMPLEMENTATION - STAY LISTENING ON ALL CHANNELS / FINE IMPLEMENTAZIONE PEDAL - RIMANE IN ASCOLTO SU TUTTI I CANALI ===============

// ===================================================================================
// CALLBACK LEVEL: DECODING OF RECEIVED PACKETS / DECODIFICA DEI PACCHETTI RICEVUTI
// ===================================================================================
// This block contains the parsing brain. By target (#ifdef)
// packets are routed to different functions (e.g. the Dongle will send
// the HID Report to the PC, while the Gun will decode the Force Feedback to move the motors).
// Runs in a task (not interrupt) but still needs to be lightning fast.
// /
// Questo blocco contiene il cervello del parsing. In base al target (#ifdef)
// i pacchetti vengono smistati verso funzioni diverse (es. il Dongle invierà
// l'HID Report al PC, mentre la Gun decodificherà il Force Feedback per muovere i motori).
// Viene eseguito in un task (non in interrupt) ma deve comunque essere velocissimo.

// ======================= CALLBACK PACKET DONGLE =========================================
#ifdef DONGLE
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
      break;
    }
    case PACKET_TX::KEYBOARD_TX: {
      uint8_t* ptr = &SerialWireless.packet.rxBuff[PREAMBLE_SIZE];
      if (memcmp(ptr, &keyReport_last_wifi, sizeof(keyReport_last_wifi))) {
        memcpy(&keyReport_last_wifi, ptr, sizeof(keyReport_last_wifi));
        keyReport_pending = true;
        xTaskNotifyGive(xUSBTaskHandle);
      }
      break;
    }
    // ====== consider whether to change name or change implementation / valutare se cambiare nome o cambiare implementazione ================
    case PACKET_TX::PEDAL_TX:
      {
        SerialWireless.is_pedal_wireless_comunication = true;
      }
      break;
    // ============================================================================= 
    #ifdef OPENFIRE_USE_ESPNOW_UNIFIED_PACKET
    case PACKET_TX::MOUSE_KEY_PAD_TX: {
      
      // We calculate the base address ONCE and put it in the 'ptr' pointer / Calcoliamo l'indirizzo base UNA sola volta e lo mettiamo nel puntatore 'ptr'
      uint8_t* ptr = &SerialWireless.packet.rxBuff[PREAMBLE_SIZE];

      // --- MOUSE ---
      // We compare the buffer to the last known state / Confrontiamo il buffer con l'ultimo stato noto
      if (memcmp(ptr, &absmouse5Report_last_wifi, sizeof(absmouse5Report_last_wifi))) {
        memcpy(&absmouse5Report_last_wifi, ptr, sizeof(absmouse5Report_last_wifi));
        absmouse5Report_pending = true;    
      }   

      ptr += sizeof(absmouse5Report_last_wifi); // Slide the pointer forward to the end of the mouse data / Fai scorrere il puntatore in avanti alla fine dei dati mouse

      // --- KEYBOARD / TASTIERA ---
      if (memcmp(ptr, &keyReport_last_wifi, sizeof(keyReport_last_wifi))) {
        memcpy(&keyReport_last_wifi, ptr, sizeof(keyReport_last_wifi));
        keyReport_pending = true;
      }
      
      ptr += sizeof(keyReport_last_wifi); // Slide the pointer forward to the end of the keyboard data / Fai scorrere il puntatore in avanti alla fine dei dati tastiera

      // --- GAMEPAD ---
      if (memcmp(ptr, &gamepad16Report_last_wifi, sizeof(gamepad16Report_last_wifi))) {
        memcpy(&gamepad16Report_last_wifi, ptr, sizeof(gamepad16Report_last_wifi));
        gamepad16Report_pending = true;
      }
      
      if (absmouse5Report_pending || keyReport_pending || gamepad16Report_pending) xTaskNotifyGive(xUSBTaskHandle);

      break;  
    } // Case panel lock / Chiusura del blocco del case
    #endif // OPENFIRE_USE_ESPNOW_UNIFIED_PACKET
    case PACKET_TX::CHECK_CONNECTION_LAST_DONGLE:
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
              
              // ARCHITECTURAL ISSUE: / CRITICITA' ARCHITETTURALE:
              // The Dongle must send 3 consecutive acknowledgments, but
              // must not perform any other actions in the meantime to ensure clean pairing. 
              // Consider whether to send the packet a couple of times or just once;
              // it sends it 3 times—once every 70ms.
              // /
              // Il Dongle deve inviare 3 conferme consecutive, ma
              // nel frattempo non deve fare nient'altro per assicurare l'accoppiamento pulito.
              // valutare se inviarlo un paio di volte il pacchetto o solo una volta
              // lo invia 3 volte - una volta ogni 70ms
              for (uint8_t i = 0; i<3; i++) {
                if (i>0) {
                  vTaskDelay(pdMS_TO_TICKS(69)); // equivalent to delay(1000) but non-blocking on ESP32 / equivalente a delay(1000) ma non bloccante su esp32
                  //unsigned long lastMillis_tx_packet_connection_last_dongle = millis();
                  //while ((millis() - lastMillis_tx_packet_connection_last_dongle) < 70) taskYIELD(); ///////yield(); 
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
      //Serial.println("DONGLE - connection request received");
      //Serial.println("DONGLE - arrivato richiesta di connessione");
      switch (aux_buffer[0]) {
        case CONNECTION_STATE::TX_GUN_SEARCH_DONGLE_BROADCAST:
          if ((SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) &&
              (aux_buffer[13] == espnow_wifi_channel))
          { // takes the first available gun / prende la prima gun disposnibile
            memcpy(SerialWireless.mac_esp_another_card, &aux_buffer[1], 6);
            // send connection request / invia richiesta connessione
            aux_buffer[0] = CONNECTION_STATE::TX_DONGLE_TO_GUN_PRESENCE;
            memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
            memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
            aux_buffer[13] = espnow_wifi_channel;
            SerialWireless.SendPacket((const uint8_t *)aux_buffer, 14, PACKET_TX::CONNECTION);
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::TX_DONGLE_TO_GUN_PRESENCE;
            // ensure that the data has been sent / assicurati che i dati siano stati spediti

          }
          break;
        case CONNECTION_STATE::TX_GUN_TO_DONGLE_ACCEPT:
          if ((memcmp(&aux_buffer[1],SerialWireless.mac_esp_another_card,6) == 0) && 
             (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
              SerialWireless.stato_connessione_wireless == CONNECTION_STATE::TX_DONGLE_TO_GUN_PRESENCE) {
              
              // Save data related to the gun VID, PID, player, etc. (USB spoofing) / SALVA I DATI RELATIVI ALLA GUN VID, PID, PLAYER , ECC.ECC. (USB Spoofing)
              memcpy(&usb_data_wireless, &aux_buffer[13], sizeof(usb_data_wireless));
              
              aux_buffer[0] = CONNECTION_STATE::TX_DONGLE_TO_GUN_CONFERM;
              memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
              memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
                            
              // =========================================================
              // Send the final package multiple times. / INVIARE IL PACCHETTO FINALE PIU' VOLTE  
              // =========================================================
              for (uint8_t i = 0; i<3; i++) {
                if (i>0) {
                  vTaskDelay(pdMS_TO_TICKS(69)); // equivalent to delay(1000) but non-blocking on ESP32 / equivalente a delay(1000) ma non bloccante su esp32
                  //unsigned long lastMillis_tx_packet_gun_to_dongle_conferm = millis();
                  //while ((millis() - lastMillis_tx_packet_gun_to_dongle_conferm) < 70) taskYIELD(); /////yield(); 
                }
                SerialWireless.SendPacket((const uint8_t *)aux_buffer, 13, PACKET_TX::CONNECTION);
              }
              SerialWireless.stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED;
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
#endif // DONGLE
// ======================= FINE CALLBACK PACKET DONGLE =========================================

// ======================= CALLBACK PACKET GUN =========================================
#ifdef GUN
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
        
        // pedal active = 00000001 - pedal2 active = 00000010 - both pedals active = 00000011 - no pedals pressed = 00000000
        // pedal attivo = 00000001 - pedal2 attivo = 00000010 - entrambi pedali attivi = 00000011 - tutti i pedali non premuti = 00000000
        memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
        TinyUSBDevices.pedals_wireless_state = aux_buffer[0];
        esp_timer_stop(timer_handle_pedal); 
        if (TinyUSBDevices.pedals_wireless_state) esp_timer_start_once(timer_handle_pedal, MAX_TIMEOUT_LAST_PACKET);
      }
      break;
    case PACKET_TX::CHECK_CONNECTION_LAST_DONGLE:
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
    case PACKET_TX::CHECK_CONNECTION_LAST_PEDAL:
      memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead);
      switch (aux_buffer[0])
      {
      case CONNECTION_STATE::TX_CONFERM_CONNECTION_LAST_PEDAL:
        if ((memcmp(&aux_buffer[1],peerAddress,6) == 0) && 
           (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
           SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION_PEDAL) { 
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED_WITH_LAST_PEDAL;
        }
        break;
      default:
        break;
      }     
      break;
    case PACKET_TX::CONNECTION:
      memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead); //13); //sizeof(aux_buffer)); // 13 is also fine here for the data to be copied / qui va bene anche 13 come dati da copiare
      switch (aux_buffer[0]) {
        // With the first packet, the dongle also transmits the channel to prevent the cycle from advancing and the channel... / AL PRIMO PACCHETTO IL DONGLE TRASMETTE ANCHE IL CANALE PER EVITARE CHE NEL FRATTEMPO IL CICLO VADA AVANTI ED IL CANALE 
        case CONNECTION_STATE::TX_DONGLE_TO_GUN_PRESENCE:
          if ((memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) && SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
            memcpy(SerialWireless.mac_esp_another_card, &aux_buffer[1], 6);
            usb_data_wireless.channel= aux_buffer[13];
            aux_buffer[0] = CONNECTION_STATE::TX_GUN_TO_DONGLE_ACCEPT;
            memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
            memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
            // Also sends data related to the gun's VID, PID, etc. / INVIA ANCHE DATI RELATIVI A VID, PID, ECC,ECC, DELLA GUN
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
        default:
          break;
      }
      break;
    case PACKET_TX::CONNECTION_PEDAL:
      memcpy(aux_buffer, &SerialWireless.packet.rxBuff[PREAMBLE_SIZE], SerialWireless.packet.bytesRead); //13); //sizeof(aux_buffer)); // 13 is also fine here for the data to be copied / qui va bene anche 13 come dati da copiare
      switch (aux_buffer[0]) {
        // With the first packet, the dongle also transmits the channel to prevent the cycle from advancing and the channel... / AL PRIMO PACCHETTO IL DONGLE TRASMETTE ANCHE IL CANALE PER EVITARE CHE NEL FRATTEMPO IL CICLO VADA AVANTI ED IL CANALE 
        case CONNECTION_STATE::TX_PEDAL_TO_GUN_PRESENCE:
          if ((memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) && SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) {
            memcpy(SerialWireless.mac_esp_another_card, &aux_buffer[1], 6);
            usb_data_wireless.channel= espnow_wifi_channel;  // shouldn't be necessary / non dovrebbe servire
            aux_buffer[0] = CONNECTION_STATE::TX_GUN_TO_PEDAL_ACCEPT;
            memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
            memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
            // Also sends data related to the gun's VID, PID, etc. / INVIA ANCHE DATI RELATIVI A VID, PID, ECC,ECC, DELLA GUN
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
        default:
          break;
      }
      break;
    default:
      break;
  }
}
#endif // GUN
// ======================= END CALLBACK PACKET GUN / FINE CALLBACK PACKET GUN =========================================

// ======================= CALLBACK PACKET PEDAL =========================================
#ifdef PEDAL
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
    case PACKET_TX::MOUSE_KEY_PAD_TX: { 
      break;  
    } 
    #endif // OPENFIRE_USE_ESPNOW_UNIFIED_PACKET
    case PACKET_TX::CHECK_CONNECTION_LAST_PEDAL:
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
              // consider whether to send the packet a couple of times or just once / valutare se inviarlo un paio di volte il pacchetto o solo una volta
              // sends it 3 times – once every 70ms / lo invia 3 volte - una volta ogni 70ms
              for (uint8_t i = 0; i<3; i++) {
                if (i>0) {
                  vTaskDelay(pdMS_TO_TICKS(69)); // equivalent to delay(1000) but non-blocking on ESP32 / equivalente a delay(1000) ma non bloccante su esp32
                  //unsigned long lastMillis_tx_packet_connection_last_pedal = millis();
                  //while ((millis() - lastMillis_tx_packet_connection_last_pedal) < 70) taskYIELD();  ////yield(); 
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
      //Serial.println("DONGLE - connection request received");
      //Serial.println("DONGLE - arrivato richiesta di connessione");
      switch (aux_buffer[0]) {
        case CONNECTION_STATE::TX_GUN_SEARCH_PEDAL_BROADCAST:
          if ((SerialWireless.stato_connessione_wireless == CONNECTION_STATE::NONE_CONNECTION) &&
              (aux_buffer[13] == espnow_wifi_channel))  // system espnow_wifi_channel / sistema espnow_wifi_channel
          { // takes the first available gun / prende la prima gun disposnibile
            memcpy(SerialWireless.mac_esp_another_card, &aux_buffer[1], 6);
            // send connection request / invia richiesta connessione
            aux_buffer[0] = CONNECTION_STATE::TX_PEDAL_TO_GUN_PRESENCE;
            memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
            memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
            SerialWireless.SendPacket((const uint8_t *)aux_buffer, 13, PACKET_TX::CONNECTION_PEDAL);
            SerialWireless.stato_connessione_wireless = CONNECTION_STATE::TX_PEDAL_TO_GUN_PRESENCE;
            // ensure that the data has been sent / assicurati che i dati siano stati spediti

          }
          break;
        case CONNECTION_STATE::TX_GUN_TO_PEDAL_ACCEPT:
          if ((memcmp(&aux_buffer[1],SerialWireless.mac_esp_another_card,6) == 0) && 
             (memcmp(&aux_buffer[7],SerialWireless.mac_esp_inteface,6) == 0) &&
              SerialWireless.stato_connessione_wireless == CONNECTION_STATE::TX_PEDAL_TO_GUN_PRESENCE) {
              
                // Save data related to the gun (VID, PID, player, etc.) / SALVA I DATI RELATIVI ALLA GUN VID, PID, PLAYER , ECC.ECC.
              memcpy(&usb_data_wireless, &aux_buffer[13], sizeof(usb_data_wireless));
              
              aux_buffer[0] = CONNECTION_STATE::TX_PEDAL_TO_GUN_CONFERM;
              memcpy(&aux_buffer[1], SerialWireless.mac_esp_inteface, 6);
              memcpy(&aux_buffer[7], SerialWireless.mac_esp_another_card, 6);
                            
              // =========================================================
              // Send the final package multiple times. / INVIARE IL PACCHETTO FINALE PIU' VOLTE  
              // =========================================================
              for (uint8_t i = 0; i<3; i++) {
                if (i>0) {
                  vTaskDelay(pdMS_TO_TICKS(69)); // equivalent to delay(1000) but non-blocking on ESP32 / equivalente a delay(1000) ma non bloccante su esp32
                  //unsigned long lastMillis_tx_packet_gun_to_pedal_conferm = millis();
                  //while ((millis() - lastMillis_tx_packet_gun_to_pedal_conferm) < 70) taskYIELD();  ////yield(); 
                }
                SerialWireless.SendPacket((const uint8_t *)aux_buffer, 13, PACKET_TX::CONNECTION_PEDAL);
              }
              SerialWireless.stato_connessione_wireless = CONNECTION_STATE::DEVICES_CONNECTED;
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
#endif // PEDAL
// ======================= FINE CALLBACK PACKET PEDAL =========================================

// ===================================================================================
// LOW-LEVEL ESP-NOW INTERRUPT (ISR)
// ===================================================================================
// These two functions are executed directly by the WiFi subsystem (FreeRTOS ISR).
// It is strictly forbidden to use blocking functions or semaphores within them.
// /
// Queste due funzioni sono eseguite direttamente dal sottosistema WiFi (FreeRTOS ISR).
// È assolutamente vietato usare funzioni bloccanti o semafori al loro interno.

// =================== ESP-NOW Receive and Send Callbacks / CALLBACK DI RICEZIONE ED INVIO DI ESP-NOW ========================

// CALLBACK RICEZIONE ESP-NOW
static void _esp_now_rx_cb(const esp_now_recv_info_t *info, const uint8_t *data, int len) {
  
  // // Broadcast filter / Filtro broadcast
  if ((info->des_addr[0] == 0xFF) && (broadcast_receiver == false)) return;

  // Local index cache (The heart of Lock-Free) / Cache locale degli indici (Il cuore del Lock-Free)
  uint16_t w = SerialWireless._writer;
  uint16_t r = SerialWireless._reader; // Let's see where the task is located right now. / Leggiamo dove si trova il task in questo istante
  const uint16_t MASK = FIFO_SIZE_READ - 1; 

  // Mathematical calculation of available space (Minus-One Rule)
  // This magic formula handles the wrap-around automatically and ensures that the
  // buffer never becomes 100% full, always leaving a 1-byte buffer
  // to distinguish the "empty" state from the "full" state.
  // /
  // Calcolo matematico dello spazio disponibile (Regola del Meno Uno)
  // Questa formula magica gestisce il wrap-around da sola e garantisce che il 
  // buffer non diventi mai pieno al 100%, lasciando sempre 1 byte di cuscinetto
  // per distinguere lo stato "vuoto" da quello "pieno".
  uint16_t free_space = (r - w - 1) & MASK;

  // Capacity check / Controllo capienza
  if (free_space >= len) {
    uint16_t firstChunk = FIFO_SIZE_READ - w;

    // BATCH COPY / COPIA BATCH
    if (firstChunk < len) {
      // Line break case: we need to split the copy in two. / Caso a capo: dobbiamo dividere la copia in due
      uint16_t secondChunk = len - firstChunk;
      memcpy(SerialWireless._queue + w, data, firstChunk);
      memcpy(SerialWireless._queue, data + firstChunk, secondChunk);
    } 
    else {
      // Linear case: single copy / Caso lineare: copia singola
      memcpy(SerialWireless._queue + w, data, len);
    }

    // --- MEMORY BARRIER / BARRIERA DI MEMORIA ---
    // Crucial: ensures that RAM has saved the memcpy data
    // BEFORE signaling to the task that new data is available.
    // /
    // Fondamentale: garantisce che la RAM abbia salvato i dati della memcpy 
    // PRIMA di segnalare al task che ci sono nuovi dati disponibili.
    asm volatile ("memw" : : : "memory");

    // ATOMIC update of the _writer index only
    // No lock required! As soon as this line executes,
    // the task (watching _writer) magically sees the new data appear.
    // /
    // Aggiornamento ATOMICO del solo indice _writer
    // Nessun lock necessario! Appena questa riga viene eseguita, 
    // il task (che sta guardando _writer) vede magicamente apparire i nuovi dati.
    SerialWireless._writer = (w + len) & MASK;   
    
  } 
  else {
    // Insufficient space (or packet discarded to avoid pointer crossing)
    // Spazio insufficiente (o pacchetto scartato per evitare l'incrocio dei puntatori)
    SerialWireless._overflow_read = true;
  }
  
  // Task Notification (Bell icon unchanged)
  // Notifica Task (Campanello invariato)
  BaseType_t xHigherPriorityTaskWoken = pdFALSE;
  vTaskNotifyGiveFromISR(xRadioTaskHandle, &xHigherPriorityTaskWoken);
  if (xHigherPriorityTaskWoken) {
    portYIELD_FROM_ISR();
  }
}

// ESP-NOW SEND CALLBACK / CALLBACK SPEDIZIONE ESP-NOW
static void _esp_now_tx_cb(const esp_now_send_info_t *tx_info, esp_now_send_status_t status) { 
  // Let's free up the radio for the next packet. / Liberiamo la radio per il prossimo pacchetto
  radioFree = true;

  // Lock-free check: Does the TX buffer still have data to send? / Controllo Lock-Free: Il buffer TX ha ancora dati da spedire?
  if (SerialWireless.writeIndex != SerialWireless.readIndex) {   
    
    // Let's wake up the radioTask to have it empty the rest. / Svegliamo il radioTask per fargli svuotare il resto
    BaseType_t xHigherPriorityTaskWoken = pdFALSE;
    vTaskNotifyGiveFromISR(xRadioTaskHandle, &xHigherPriorityTaskWoken);
    
    if (xHigherPriorityTaskWoken) {
      portYIELD_FROM_ISR();
    }
  }
}
// =================== END ESP-NOW Receive and Send Callbacks / FINE CALLBACK DI RICEZIONE ED INVIO DI ESP-NOW ========================


// ===================================================================================
// HARDWARE TIMERS (SERIAL AND PEDAL) / (SERIALE E PEDALE)
// ===================================================================================
// Using esp_timer (high-resolution hardware timer) instead of delay() or millis()
// ensures precise, high-priority execution—essential for unblocking a
// congested serial connection or declaring a pedal "dead" if packets are lost.
// /
// L'uso di esp_timer (alta risoluzione hardware) invece di delay() o millis() 
// garantisce un'esecuzione prioritaria e precisa, essenziale per sbloccare la 
// seriale intasata o per dichiarare un pedale "morto" se si perdono pacchetti.

// ================================== Serial Port Timer / TIMER PER SERIALE ===================================
// CALLBACK per timer Serial

void timer_callback_serial(void* arg) {
  while (SerialWireless.availableBufferSerialWrite() > 0) {
    if (!SerialWireless.flush_sem()) {
      // The radio is full! The timer cannot yield. 
      // We ask the radioTask for help for when space becomes available.
      // /
      // LA RADIO È PIENA! Il timer non può fare yield.
      // Chiediamo aiuto al radioTask per quando si libererà spazio.
      SerialWireless._serial_needs_recovery = true; 
      break; 
    }
  }
}

void SerialWireless_::setupTimerSerial() {
    esp_timer_create_args_t timer_args = {
        .callback = &timer_callback_serial,
        .arg = NULL,
        .dispatch_method = ESP_TIMER_TASK,
        .name = "timer_serial"
    };
    
    esp_timer_create(&timer_args, &timer_handle_serial);
}


// ==========================   END Serial Port Timer / FINE TIMER PER SERIALE ====================================


// ================================== PEDAL TIMER / TIMER PER PEDAL ===================================
// CALLBACK per timer Pedal
void timer_callback_pedal(void* arg) {
  
  TinyUSBDevices.pedals_wireless_state = 0;
 
}

void setupTimerPedal() {
    esp_timer_create_args_t timer_args = {
        .callback = &timer_callback_pedal,
        .arg = NULL,
        .dispatch_method = ESP_TIMER_TASK,
        .name = "timer_pedal"
    };
    
    esp_timer_create(&timer_args, &timer_handle_pedal);
}
// ======================== END PEDAL TIMER / FINE TIMER PER PEDAL =============================

// ===================================================================================
// SYSTEM CALLBACK: TINYUSB / CALLBACK DI SISTEMA: TINYUSB
// ===================================================================================

// ============ TinyUSB Send Callback / CALLBACK PER SPEDIZIONE TINYUSB =============================
#ifdef DONGLE
// It does not need to be defined, as it replaces the default one that does nothing (it is already defined within the tinyUSB library)
// /
// NON VA DEFINITA IN QUANTO SOSTITUISCE QUELLA DI DEFAULT CHE NON FA NULLA (e' gia' definita dentro la libreria tinyUSB)
#ifdef __cplusplus
extern "C" {
#endif

// This callback is invoked by the low-level USB driver whenever
// the PC "consumes" the HID data sent by the Dongle. We intercept it to
// wake up the usbTask and immediately send the next packet.
// /
// Questa callback viene invocata dal driver USB a basso livello ogni volta che 
// il PC "assorbe" i dati HID inviati dal Dongle. Noi la intercettiamo per
// risvegliare l'usbTask e inviare immediatamente il pacchetto successivo.
void tud_hid_report_complete_cb(uint8_t instance, uint8_t const* report, uint16_t len) {
    // The PC has just freed the buffer. Let's wake up usbTask! / Il PC ha appena liberato il buffer. Svegliamo usbTask!
    BaseType_t xHigherPriorityTaskWoken = pdFALSE;
    vTaskNotifyGiveFromISR(xUSBTaskHandle, &xHigherPriorityTaskWoken);
    if (xHigherPriorityTaskWoken) {
        portYIELD_FROM_ISR();
    }
}

#ifdef __cplusplus
}
#endif

#endif // DONGLE
// ============ END TinyUSB Send Callback / FINE CALLBACK PER SPEDIZIONE TINYUSB =============================


#endif //OPENFIRE_WIRELESS_ENABLE