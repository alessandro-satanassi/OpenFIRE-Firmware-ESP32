/*!
 * @file main.cpp
 * @brief main DONGLE
 * @n CPP main DONGLE
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

#include "TinyUSB_Devices.h"

// ===================================================================================
// HARDWARE ABSTRACTION: GESTIONE GRAFICA OPZIONALE E COMPATIBILITA'
// ===================================================================================
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
    #include "LGFX_096_ST7735S_80x160.hpp"

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
  #else
    #include <Adafruit_ST7735.h>
    Adafruit_ST7735 tft = Adafruit_ST7735(TFT_CS, TFT_DC, TFT_MOSI, TFT_SCLK, TFT_RST);

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

#include "OpenFIRE-DONGLE-version.h"

bool display_init = false;

// ===================================================================================
// MAIN SETUP
// ===================================================================================

// The main show!
void setup() { 

  // 1. INIZIALIZZAZIONE HARDWARE (Video Feedback)
  #ifdef USES_DISPLAY
    #ifdef USE_LOVYAN_GFX
      display_init = tft.init();
      if(display_init) {
      tft.setSwapBytes(true);
      tft.setColorDepth(16);
      }
    #else
      tft.initR(INITR_MINI160x80_PLUGIN); // non riporta alcun valore, quindi display_init va inizializzato da solo
      display_init = true;
      if(display_init) {  // Init ST7735S mini display
      pinMode(TFT_PIN_BL, OUTPUT);
      digitalWrite(TFT_PIN_BL, 0); // accende retroilluminazione del display
      }
    #endif // USE_LOVYAN_GFX
    
    if(display_init) tft.setRotation(3);
    
    if(display_init) {
    tft.fillScreen(BLACK);
    #ifdef USE_LOVYAN_GFX
      tft.pushImage(10,1,LOGO_RGB_ALPHA_OPEN_WIDTH, LOGO_RGB_ALPHA_OPEN_HEIGHT, (uint16_t *)logo_rgb_alpha_open);
    #else
      tft.drawRGBBitmap(10,1,(uint16_t *)logo_rgb_alpha_open,LOGO_RGB_ALPHA_OPEN_WIDTH,LOGO_RGB_ALPHA_OPEN_HEIGHT);
    #endif //USE_LOVYAN_GFX
    
    vTaskDelay(pdMS_TO_TICKS(3000));

    tft.fillScreen(BLACK);
    tft.setTextSize(2);
    tft.setCursor(0, 0);
    tft.setTextColor(RED);
    tft.drawBitmap(40, 0, customSplashBanner, CUSTSPLASHBANN_WIDTH, CUSTSPLASHBANN_HEIGHT, BLUE);
    tft.drawBitmap(10, 25, customSplash, CUSTSPLASH_WIDTH, CUSTSPLASH_HEIGHT, RED);
    tft.setTextSize(2);
    tft.setTextColor(WHITE);
    }
  #endif //USES_DISPLAY

  // ===================================================================================
  // GESTIONE CONNESSIONE WIRELESS E ASSOCIAZIONE DISPOSITIVO (Bloccante)
  // ===================================================================================
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

  // =========== attesa segnale del pedale da parte della lightgun
  const unsigned long start_millis = millis();
  const unsigned long timeout_ms = (MAX_SECONDI_CONNESSIONE_PEDAL + 2) * 1000;
  SerialWireless.is_pedal_wireless_comunication = false;
  while (!SerialWireless.is_pedal_wireless_comunication && (millis() - start_millis < timeout_ms)) {
    vTaskDelay(pdMS_TO_TICKS(100));
  }

  // ===================================================================================
  // CONNESSIONE USB - imposta VID e PID come quello che gli passa la pistola
  // ===================================================================================
  // Per garantire compatibilità totale con i front-end (MAME, Batocera, Windows),
  // il Dongle assume dinamicamente l'identità (VID, PID e Manufacturer) inviata 
  // precedentemente dalla Lightgun, comportandosi come un proxy hardware invisibile.
  
  if (!TinyUSBDevice.isInitialized()) { 
    TinyUSBDevice.begin(0);
  }
      
  TinyUSBDevice.setManufacturerDescriptor(usb_data_wireless.deviceManufacturer);
  TinyUSBDevice.setProductDescriptor(usb_data_wireless.deviceName);
  TinyUSBDevice.setID(usb_data_wireless.deviceVID, usb_data_wireless.devicePID);
  
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
  // VISUALIZZAZIONE DATI DELLA LIGHTGUN SUL DISPLAY
  // ===================================================================================
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
// PARAMETRI BUFFERIZZAZIONE SERIALE -> RADIO
// ===================================================================================
// FIFO_SIZE a 200 garantisce l'assorbimento di rapidi
// treni di dati dal Mamehooker (es. force feedback multipli come Rumble + Solenoide) 
// senza causare colli di bottiglia o desincronizzazione della porta COM

#define FIFO_SIZE_READ_SER 200
#define TIME_OUT_AVALAIBLE 2
#define TIME_OUT_SERIAL_MICRO 1000 // 1000 microsecondi = 1 millisecondo

int rx_avalaible = 0;
unsigned long startTime = 0; // = millis();

uint64_t timer_serial_micro = 0;
uint8_t buffer_aux[FIFO_SIZE_READ_SER];

// ===================================================================================
// MAIN LOOP: IL BRIDGE SERIALE-RADIO (MAMEHOOKER/FORCE FEEDBACK)
// ===================================================================================
// NOTA ARCHITETTURALE IMPORTANTE:
// Questo loop si occupa SOLTANTO della direzione PC -> DONGLE -> LIGHTGUN
// I dati in ingresso dall'host (comandi force feedback via seriale) vengono "pescati" 
// in blocchi, inseriti nel Ring Buffer radio e notificati ai task asincroni tramite
// las librei `SerialWireless`.
// La ricezione radio inversa (Lightgun -> Dongle -> USB HID PC) è gestita in modo 
// asincrono invisibile e indipendente dai Task delegati.

void loop()
{
  vTaskDelay(pdMS_TO_TICKS(1)); // diamo un po' di respiro al micro per non incorrere in reset dal Task Watchdog
  rx_avalaible = Serial.available();
  if (rx_avalaible > FIFO_SIZE_READ_SER) rx_avalaible = FIFO_SIZE_READ_SER;
  if (rx_avalaible)
  {
    Serial.readBytes(buffer_aux, rx_avalaible);
    SerialWireless.write(buffer_aux, rx_avalaible);
    
    // Il flush sveglia il task radio asincrono per spedirli fisicamente in aria senza bloccare il loop.
    SerialWireless.flush(); 
  }
}