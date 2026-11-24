/*!
 * @file main.cpp
 * @brief main SCOPE
 * @n CPP main SCOPE
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#include <Arduino.h>
#include <SPI.h>

////////////#include "TinyUSB_Devices.h"

// ===================================================================================
// HARDWARE ABSTRACTION: OPTIONAL GRAPHICS MANAGEMENT AND COMPATIBILITY / GESTIONE GRAFICA OPZIONALE E COMPATIBILITA'
// ===================================================================================
// The dongle supports various hardware configurations (with or without a TFT display)
// and is agnostic regarding the graphics library in use (LovyanGFX vs. Adafruit).
// Extensive use of preprocessor directives (#ifdef) ensures that the final binary
// includes only the code strictly necessary for the target board, reducing
// Flash memory usage.
// /
// Il Dongle supporta diverse configurazioni hardware (con o senza display TFT) 
// ed è agnostico rispetto alla libreria grafica in uso (LovyanGFX vs Adafruit). 
// L'uso intensivo di direttive preprocessore (#ifdef) assicura che il binario finale 
// includa solo il codice strettamente necessario alla board target, riducendo le 
// dimensioni della Flash.

#ifdef USES_DISPLAY
  #include "OpenFIRE_logo.h"
  #ifdef USE_LOVYAN_GFX
    #define LGFX_USE_V1
    #include <LovyanGFX.hpp>
    #include "LGFX_128_GC9A01_240x240.hpp"

    LGFX tft;

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
  #endif // USE_LOVYAN_GFX

#endif // USES_DISPLAY

//#include "OpenFIRE-SCOPE-version.h"

/*
CPU frequency: 80 MHz is the lowest the radio (Wi-Fi/ESP-NOW) supports. The dongle only forwards a few hundred small
packets per second between ESP-NOW and USB, so it has plenty of margin, and a lower clock makes the chip run a little cooler.
If start-up, pairing or play ever misbehave, try 160. Commented out = 240 MHz, as before.
/
Frequenza della CPU: 80 MHz è la minima ammessa dalla radio (Wi-Fi/ESP-NOW). Il dongle inoltra solo poche centinaia di
piccoli pacchetti al secondo tra ESP-NOW e USB, quindi ha margine in abbondanza, e un clock più basso fa scaldare un po' meno il chip.
Se avvio, accoppiamento o gioco dovessero dare problemi, provare 160. Commentato = 240 MHz, come prima.
Da ulcune prove sembra che tra 240mhz e 80mhz ci siano 4 gradi di differenza e non risultano differenze di prestazioni
*/
//#define SCOPE_CPU_FREQ_MHZ 80

/*
DEBUG_TEMP: shows the CPU frequency and the chip's internal temperature in the top right corner of the display,
refreshed every DEBUG_TEMP_INTERVAL_MS. Only for tests: the internal sensor is approximate, compare readings on the same dongle.
Not drawn with DEBUG_WIFI_ESPNOW: its counters use the whole display.
/
DEBUG_TEMP: mostra la frequenza della CPU e la temperatura interna del chip nell'angolo in alto a destra del display,
aggiornate ogni DEBUG_TEMP_INTERVAL_MS. Solo per prove: il sensore interno è approssimativo, confrontare letture sullo stesso dongle.
Non disegnato con DEBUG_WIFI_ESPNOW: i suoi contatori usano tutto il display.
*/
//#define DEBUG_TEMP

#if defined(DEBUG_TEMP) && defined(USES_DISPLAY) && !defined(DEBUG_WIFI_ESPNOW)
  #define DEBUG_TEMP_INTERVAL_MS 30000UL
  unsigned long debug_temp_millis = 0;
  bool debug_temp_first = true;
#endif // DEBUG_TEMP && USES_DISPLAY

bool display_init = false;

#if defined(DEBUG_WIFI_ESPNOW) && defined(USES_DISPLAY) && defined(OPENFIRE_WIRELESS_ENABLE) && defined(OPENFIRE_WIRELESS_DEVICE_ESPNOW)
// EN: Radio counters on the whole display, only with -D DEBUG_WIFI_ESPNOW: they replace the final screen of setup().
//     Five lines, font size 2, counted since power on:
//       TX   <datagrams sent to the gun: PC data (app, MAMEHooker) and connection replies>
//       lost <not acknowledged by the gun after the radio retries>                          (red if not 0)
//       RX   <datagrams received from the gun>
//       ERR  <refused by esp_now_send(), sent again later>/<packets dropped, transmit queue full>/
//            <datagrams dropped, receive queue full>                                       (red if one is not 0)
//       <gun signal, dBm: worst/last/best; "--" if nothing for 2 s>  (green, yellow from -68, red below -80)
//     A line that does not fit in 13 characters is drawn with font size 1 (26 characters): no digit is cut.
//     Recomposed once a second, only the characters that changed are drawn. loop() has the lowest priority:
//     the radio and USB tasks interrupt it.
// IT: Contatori radio su tutto il display, solo con -D DEBUG_WIFI_ESPNOW: sostituiscono la schermata finale di setup().
//     Cinque righe, font dimensione 2, contate dall'accensione:
//       TX   <datagrammi inviati alla pistola: dati del PC (app, MAMEHooker) e risposte di connessione>
//       lost <non confermati dalla pistola dopo i tentativi della radio>                     (rosso se non 0)
//       RX   <datagrammi ricevuti dalla pistola>
//       ERR  <rifiutati da esp_now_send(), rispediti dopo>/<pacchetti scartati, coda di trasmissione piena>/
//            <datagrammi scartati, coda di ricezione piena>                                (rosso se uno non e' 0)
//       <segnale della pistola, dBm: peggiore/ultimo/migliore; "--" se niente da 2 s>  (verde, giallo da -68, rosso sotto -80)
//     Una riga che non sta in 13 caratteri e' disegnata con il font dimensione 1 (26 caratteri): nessuna cifra e' tagliata.
//     Ricomposte una volta al secondo, si disegnano solo i caratteri cambiati. loop() ha la priorita' piu' bassa:
//     i task della radio e dell'USB lo interrompono.
#define DBG_WIFI_LINES 5
#define DBG_WIFI_COLS 13                   // size 2: 12x16 pixels, 13 x 12 = 156 of 160 / dimensione 2: 12x16 pixel, 13 x 12 = 156 su 160
static void DebugWifiEspNowDisplay()
{
  static const char *const label[DBG_WIFI_LINES] = { "TX", "lost", "RX", "ERR", "" };
  static const uint16_t labelColor[DBG_WIFI_LINES] = { CYAN, CYAN, GREEN, GRAY, GRAY };
  static char shown[DBG_WIFI_LINES][2 * DBG_WIFI_COLS + 1];   // on the display (empty = to draw) / sul display (vuota = da disegnare)
  static uint16_t shownColor[DBG_WIFI_LINES];
  static uint8_t shownSize[DBG_WIFI_LINES];   // font size on the display (0 = none yet) / dimensione del font sul display (0 = ancora nessuna)
  static uint32_t lastMs = 0;
  static bool first = true;
  const uint32_t now = millis();
  if (!display_init || (!first && now - lastMs < 1000)) return;
  if (first) {                             // replaces the final screen / sostituisce la schermata finale
    first = false;
    tft.fillScreen(BLACK);
  }
  lastMs = now;
  const uint32_t fail = dbgWifiTxFail;
  const uint32_t sendErr = dbgWifiSendErr, txFull = dbgWifiTxQueueFull, rxFull = dbgWifiRxQueueFull;
  const DbgWifiRx_t &rx = dbgWifiRx[0];     // [0]: the paired gun (peerAddress) / la pistola accoppiata (peerAddress)
  for (uint8_t i = 0; i < DBG_WIFI_LINES; ++i) {
    char value[40], text[40];
    uint16_t valueColor = WHITE;
    switch (i) {
    case 0: snprintf(value, sizeof(value), "%lu", (unsigned long)(dbgWifiTxOk + fail)); break;
    case 1: snprintf(value, sizeof(value), "%lu", (unsigned long)fail); if (fail) valueColor = RED; break;
    case 2: snprintf(value, sizeof(value), "%lu", (unsigned long)rx.count); break;
    case 3: snprintf(value, sizeof(value), "%lu/%lu/%lu", (unsigned long)sendErr, (unsigned long)txFull, (unsigned long)rxFull);
            if (sendErr || txFull || rxFull) valueColor = RED;
            break;
    default: {
      valueColor = GRAY;
      if (!rx.count) { snprintf(value, sizeof(value), "--/--/--"); break; }   // count first: the callback sets it last / prima count: la callback lo aggiorna per ultimo
      char last[6] = "--";
      // EN: Signed: a datagram arrived after "now" is recent too. / IT: Con segno: anche un datagramma arrivato dopo "now" e' recente.
      if ((int32_t)(now - rx.lastMs) < 2000) {
        const int rssi = rx.rssi;
        snprintf(last, sizeof(last), "%d", rssi);
        valueColor = (rssi >= -67) ? GREEN : (rssi >= -80) ? YELLOW : RED;
      }
      snprintf(value, sizeof(value), "%d/%s/%d", (int)rx.rssiMin, last, (int)rx.rssiMax);
      } break;
    }
    // EN: Value on the right of the label (rssi on the left), spaces clear the old text. Font size 2 if label, a space
    //     and value fit in 13 characters, otherwise size 1 (26 characters, in the middle of the line's 16 pixels).
    // IT: Valore a destra dell'etichetta (rssi a sinistra), gli spazi cancellano il testo vecchio. Font dimensione 2 se
    //     etichetta, uno spazio e valore stanno in 13 caratteri, altrimenti dimensione 1 (26 caratteri, al centro dei 16 pixel della riga).
    const uint8_t labelLen = strlen(label[i]);
    const uint8_t size = (labelLen + (labelLen ? 1 : 0) + strlen(value) <= DBG_WIFI_COLS) ? 2 : 1;
    const uint8_t cols = (size == 2) ? DBG_WIFI_COLS : 2 * DBG_WIFI_COLS;
    if (labelLen) snprintf(text, sizeof(text), "%s%*s", label[i], cols - labelLen, value);
    else snprintf(text, sizeof(text), "%-*s", cols, value);
    text[cols] = 0;
    if (size != shownSize[i]) {             // size changed: clear the line / dimensione cambiata: cancella la riga
      if (shownSize[i]) tft.fillRect(0, 16 * i, 160, 16, BLACK);
      shown[i][0] = 0;
    }
    if (valueColor != shownColor[i]) shown[i][0] = 0;   // colour changed: whole line / colore cambiato: tutta la riga
    const bool all = !shown[i][0];
    tft.setTextSize(size);
    for (uint8_t k = 0; k < cols; ++k) {
      if (!all && text[k] == shown[i][k]) continue;
      tft.setTextColor((k < labelLen) ? labelColor[i] : valueColor, BLACK);
      tft.setCursor(2 + 6 * size * k, 16 * i + 8 - 4 * size);
      tft.print(text[k]);
    }
    memcpy(shown[i], text, cols + 1);
    shownColor[i] = valueColor;
    shownSize[i] = size;
  }
}
#undef DBG_WIFI_COLS
#undef DBG_WIFI_LINES
#endif // DEBUG_WIFI_ESPNOW
// ===================================================================================
// MAIN SETUP
// ===================================================================================

// The main show!
void setup() { 

  #ifdef SCOPE_CPU_FREQ_MHZ
    setCpuFrequencyMhz(SCOPE_CPU_FREQ_MHZ);
  #endif // SCOPE_CPU_FREQ_MHZ

  // 1. HARDWARE INITIALIZATION / INIZIALIZZAZIONE HARDWARE (Video Feedback)
  #ifdef USES_DISPLAY
    #ifdef USE_LOVYAN_GFX
      display_init = tft.init();
      if(display_init) {
      tft.setSwapBytes(true);
      tft.setColorDepth(16);
      }
    #endif // USE_LOVYAN_GFX
    
    if(display_init) tft.setRotation(3);
    
    if(display_init) {
    tft.fillScreen(BLACK);
    #ifdef USE_LOVYAN_GFX
      tft.pushImage(10,1,LOGO_RGB_ALPHA_OPEN_WIDTH, LOGO_RGB_ALPHA_OPEN_HEIGHT, (uint16_t *)logo_rgb_alpha_open);
    #endif //USE_LOVYAN_GFX
    
    vTaskDelay(pdMS_TO_TICKS(3000));

    tft.fillScreen(BLACK);
    tft.setTextSize(2);
    tft.setCursor(0, 0);
    tft.setTextColor(RED);
    tft.drawBitmap(40, 0, customSplashBanner, CUSTSPLASHBANN_WIDTH, CUSTSPLASHBANN_HEIGHT, BLUE);
    // tft.drawBitmap(10, 25, customSplash, CUSTSPLASH_WIDTH, CUSTSPLASH_HEIGHT, RED);
    #ifdef USE_LOVYAN_GFX
      tft.pushImage(9, 24, LOGO_RGB_ALPHA_WIDTH, LOGO_RGB_ALPHA_HEIGHT, (uint16_t *)logo_rgb_alpha);
    #endif
    tft.setTextSize(2);
    tft.setTextColor(WHITE);
    }
  #endif //USES_DISPLAY

  // ===================================================================================
  // Wireless Connection Management and Device Pairing (Blocking) / GESTIONE CONNESSIONE WIRELESS E ASSOCIAZIONE DISPOSITIVO (Bloccante)
  // ===================================================================================
  // The dongle searches for the best channel to tune into and then waits
  // to establish the ESP-NOW link with the Lightgun. Until it obtains
  // the correct MAC address and identification parameters, it CANNOT proceed.
  // /
  // Il Dongle cerca il canale migliore su cui sintonizzarsi  e poi attende 
  // di stabilire il link ESP-NOW con la Lightgun. Finché non ottiene
  // il MAC Address corretto e i parametri identificativi, NON può procedere.
  
  SerialWireless.init_wireless();
  SerialWireless.begin();
  SerialWireless.connection_dongle();

  #ifdef USES_DISPLAY   
  if(display_init) {
    tft.fillScreen(BLACK);
    tft.drawBitmap(40, 0, customSplashBanner, CUSTSPLASHBANN_WIDTH, CUSTSPLASHBANN_HEIGHT, BLUE);
    tft.setTextSize(2);
    tft.setTextColor(RED);
    tft.setCursor(0, 40);
    tft.printf("AWAIT PEDAL");
  }
  #endif // USES_DISPLAY

  // =========== waiting for the pedal signal from the light gun / attesa segnale del pedale da parte della lightgun ===========
  const unsigned long start_millis = millis();
  const unsigned long timeout_ms = (MAX_SECONDI_CONNESSIONE_PEDAL + 2) * 1000;
  SerialWireless.is_pedal_wireless_comunication = false;
  while (!SerialWireless.is_pedal_wireless_comunication && (millis() - start_millis < timeout_ms)) {
    vTaskDelay(pdMS_TO_TICKS(100));
  }

  // ===================================================================================
  // USB CONNECTION – set the VID and PID to match those provided by the lightgun. / CONNESSIONE USB - imposta VID e PID come quello che gli passa la pistola
  // ===================================================================================
  // To ensure full compatibility with front-ends (MAME, Batocera, Windows),
  // the dongle dynamically adopts the identity (VID, PID, and Manufacturer) previously
  // sent by the lightgun, acting as an invisible hardware proxy.
  // /
  // Per garantire compatibilità totale con i front-end (MAME, Batocera, Windows),
  // il Dongle assume dinamicamente l'identità (VID, PID e Manufacturer) inviata 
  // precedentemente dalla Lightgun, comportandosi come un proxy hardware invisibile.
  
  if (!TinyUSBDevice.isInitialized()) { 
    TinyUSBDevice.begin(0);
  }
      
  TinyUSBDevice.setManufacturerDescriptor(usb_data_wireless.deviceManufacturer);
  TinyUSBDevice.setProductDescriptor(usb_data_wireless.deviceName);
  TinyUSBDevice.setID(usb_data_wireless.deviceVID, usb_data_wireless.devicePID);
  
  // Adds the USB port description serial number so that the same COM port is always mapped for each lightgun.
  // Formats the byte array into a hexadecimal string.
  // /
  // Aggiunge il serial number della descrizione della porta USB in modo che venga mappata sempre la stessa COM per ogni lightgun
  // Formatta l'array di byte in una stringa esadecimale.
  static char serialNumber[13]; 
  snprintf(serialNumber, sizeof(serialNumber), "%02X%02X%02X%02X%02X%02X", MAC2STR(peerAddress)); // SerialWireless.mac_esp_another_card  
  TinyUSBDevice.setSerialDescriptor(serialNumber);
  
  // Initializing the USB devices chunk.
  TinyUSBDevices.begin(1);
  Serial.begin(9600);
  Serial.setTimeout(0);
    
  // ===================================================================================
  // DISPLAYING LIGHT GUN DATA ON THE SCREEN / VISUALIZZAZIONE DATI DELLA LIGHTGUN SUL DISPLAY
  // ===================================================================================
  // Final visual feedback confirming the successful cloning of the identity.
  // /
  // Feedback visivo finale che conferma l'avvenuta clonazione dell'identità.
 
   #ifdef USES_DISPLAY   
    if(display_init) {
    tft.fillScreen(BLACK);
    tft.drawBitmap(40, 0, customSplashBanner, CUSTSPLASHBANN_WIDTH, CUSTSPLASHBANN_HEIGHT, BLUE);
    tft.setTextSize(2);
    tft.setCursor(0, 20);
    tft.setTextColor(RED);
    tft.println(usb_data_wireless.deviceName);
    tft.setTextColor(RED);
    tft.setCursor(0, 40);
    tft.printf("Player: %d", usb_data_wireless.devicePlayer);
    tft.setTextSize(2);
    tft.setCursor(0, 60);
    tft.setTextColor(GRAY); 
    tft.printf("Channel: %d", usb_data_wireless.channel);
    }
  #endif // USES_DISPLAY
}

// ===================================================================================
// SERIAL BUFFERING PARAMETERS -> RADIO / PARAMETRI BUFFERIZZAZIONE SERIALE -> RADIO
// ===================================================================================
// Setting FIFO_SIZE to 200 ensures the handling of rapid
// data bursts from Mamehooker (e.g., multiple force feedback effects like Rumble + Solenoid)
// without causing bottlenecks or COM port desynchronization.
// /
// FIFO_SIZE a 200 garantisce l'assorbimento di rapidi
// treni di dati dal Mamehooker (es. force feedback multipli come Rumble + Solenoide) 
// senza causare colli di bottiglia o desincronizzazione della porta COM

#define FIFO_SIZE_READ_SER 200
#define TIME_OUT_AVALAIBLE 2
#define TIME_OUT_SERIAL_MICRO 1000 // 1000 microseconds = 1 millisecond / 1000 microsecondi = 1 millisecondo

int rx_avalaible = 0;
unsigned long startTime = 0; // = millis();

uint64_t timer_serial_micro = 0;
uint8_t buffer_aux[FIFO_SIZE_READ_SER];

// ===================================================================================
// MAIN LOOP: THE SERIAL-RADIO BRIDGE (MAMEHOOKER/FORCE FEEDBACK) / IL BRIDGE SERIALE-RADIO (MAMEHOOKER/FORCE FEEDBACK)
// ===================================================================================
// IMPORTANT ARCHITECTURAL NOTE:
// This loop handles ONLY the PC -> DONGLE -> LIGHTGUN direction.
// Incoming data from the host (serial force feedback commands) is "fetched"
// in blocks, placed into the radio Ring Buffer, and signaled to the asynchronous tasks
// via the `SerialWireless` library.
// Reverse radio reception (Lightgun -> Dongle -> USB HID PC) is handled
// asynchronously—invisibly and independently of the delegated tasks.
// /
// NOTA ARCHITETTURALE IMPORTANTE:
// Questo loop si occupa SOLTANTO della direzione PC -> DONGLE -> LIGHTGUN
// I dati in ingresso dall'host (comandi force feedback via seriale) vengono "pescati" 
// in blocchi, inseriti nel Ring Buffer radio e notificati ai task asincroni tramite
// las librei `SerialWireless`.
// La ricezione radio inversa (Lightgun -> Dongle -> USB HID PC) è gestita in modo 
// asincrono invisibile e indipendente dai Task delegati.

void loop()
{
  vTaskDelay(pdMS_TO_TICKS(1)); // Let's give the microcontroller some breathing room to avoid resets triggered by the Task Watchdog. / diamo un po' di respiro al micro per non incorrere in reset dal Task Watchdog
  rx_avalaible = Serial.available();
  if (rx_avalaible > FIFO_SIZE_READ_SER) rx_avalaible = FIFO_SIZE_READ_SER;
  if (rx_avalaible)
  {
    Serial.readBytes(buffer_aux, rx_avalaible);
    SerialWireless.write(buffer_aux, rx_avalaible);
    
    // The flush wakes up the asynchronous radio task to physically transmit them over the air without blocking the loop. / Il flush sveglia il task radio asincrono per spedirli fisicamente in aria senza bloccare il loop.
    SerialWireless.flush(); 
  }

  #if defined(DEBUG_TEMP) && defined(USES_DISPLAY) && !defined(DEBUG_WIFI_ESPNOW)
  // Right of the banner (x 120-159, y 0-15), free on every screen after setup. / A destra del banner (x 120-159, y 0-15), libero in ogni schermata dopo il setup.
  if (display_init && (debug_temp_first || (millis() - debug_temp_millis >= DEBUG_TEMP_INTERVAL_MS))) {
    debug_temp_first = false;
    debug_temp_millis = millis();
    char debug_temp_buf[8];
    tft.setTextSize(1);
    tft.setTextColor(YELLOW, BLACK);
    snprintf(debug_temp_buf, sizeof(debug_temp_buf), "%3luMHz", (unsigned long)getCpuFrequencyMhz());
    tft.setCursor(122, 0);
    tft.print(debug_temp_buf);
    snprintf(debug_temp_buf, sizeof(debug_temp_buf), "%5.1fC", temperatureRead());
    tft.setCursor(122, 8);
    tft.print(debug_temp_buf);
  }
  #endif // DEBUG_TEMP && USES_DISPLAY
  #if defined(DEBUG_WIFI_ESPNOW) && defined(USES_DISPLAY) && defined(OPENFIRE_WIRELESS_ENABLE) && defined(OPENFIRE_WIRELESS_DEVICE_ESPNOW)
  DebugWifiEspNowDisplay();
  #endif // DEBUG_WIFI_ESPNOW
}