// by SATANASSI Alessandro

#include <Arduino.h>
#include <SPI.h>

#include "TinyUSB_Devices.h"


#include "OpenFIRE-PEDAL-version.h"

// ================== GESTIONE DUAL CORE ============================
#if defined(DUAL_CORE) && defined(ESP_PLATFORM) && false
        void setup1();
        void loop1();
        TaskHandle_t task_loop1;
        void esploop1(void* pvParameters) {
            setup1();
            for (;;) loop1();
        }
#endif //DUAL_CORE
// ========================= FINE GESTIONE DUAL CORE ======================================

bool display_init = false;

const int8_t leds[4] = {PIN_LED1, 
                        PIN_LED2, 
                        PIN_LED3, 
                        PIN_LED4}; // I tuoi 4 GPIO


void animTaskLink(void *pvParameters) {
  
  for (uint8_t i = 0; i < 4; i++) {
    pinMode(leds[i], OUTPUT);
    digitalWrite(leds[i], LOW); 
  }
  
  // Variabili di stato dell'animazione (fuori dal loop infinito)
  int currentLed = 0;
  int direction = 1;

  // Loop del task
  for (;;) {
    // ========= inserire codice di animazione led 'kit supercar' qui ========
    // 1. Accende SOLO il led attuale
    digitalWrite(leds[currentLed], HIGH);

    // 2. Mantiene l'animazione in pausa per far vedere il led acceso
    vTaskDelay(pdMS_TO_TICKS(150)); 

    // 3. Il tempo è scaduto: spegne SOLO il led attuale PRIMA di passare al prossimo
    digitalWrite(leds[currentLed], LOW);

    // 4. Prepara l'indice per il prossimo giro
    currentLed += direction;

    // 5. Se tocca i bordi (0 o 3), inverte la marcia
    if (currentLed >= 3) {
      direction = -1;
    } else if (currentLed <= 0) {
      direction = 1;
    }
    // =======================================================================
    //vTaskDelay(pdMS_TO_TICKS(150)); // piccolo delay per non saturare la CPU
  } 
}

// The main show!
void setup() {

  // Configurazione Pedali (Ingressi con Pull-Up)
  pinMode(PIN_PEDAL, INPUT_PULLUP);
  pinMode(PIN_PEDAL2, INPUT_PULLUP);

  TaskHandle_t animTaskHandleLink = NULL;  
  xTaskCreatePinnedToCore(
        animTaskLink,          // funzione del task
        "AnimTaskLink",        // nome
        4096,              // stack size
        NULL,              // parametri
        1,                 // priorità
        &animTaskHandleLink,   // handle
        APP_CPU_NUM        // core (puoi usare 0 o 1)
  );  
  
  // ====== gestione connessione wireless ====================
  SerialWireless.init_wireless();
  SerialWireless.begin();
  SerialWireless.connection_pedal();
  // ====== fine gestione wireless .. va avanti solo dopo che si è accoppiato il dispositivo =======


  ///////////////////////// POI ANDRA' TOLTO ////////////////////////////////////////////////
  // ====== connessione USB ====== imposta VID e PID come quello che gli passa la pistola ===============
  if (!TinyUSBDevice.isInitialized()) { // aggiunto ..funzionava lo stesso, ma così è più sicuro .. sicuramente serve per Esp32 con libreria non integrfata nel core
    TinyUSBDevice.begin(0);
  }
      
  TinyUSBDevice.setManufacturerDescriptor(usb_data_wireless.deviceManufacturer);
  TinyUSBDevice.setProductDescriptor(usb_data_wireless.deviceName);
  TinyUSBDevice.setID(usb_data_wireless.deviceVID, usb_data_wireless.devicePID);

  // Initializing the USB devices chunk.
  TinyUSBDevices.begin(1);
  //////////////////////////// FINE POI ANDRA' TOLDO /////////////////////////////////////////


  Serial.begin(9600);
  Serial.setTimeout(0);
  ////////////////////Serial.setTxTimeoutMs(0);
  // ====== fine connessione USB ==========================================================================
  #ifdef OPENFIRE_ESPNOW_WIFI_POWER_AUTO
  vTaskDelay(pdMS_TO_TICKS(100)); // per aggiornare espnow_wifi_power
  #endif //OPENFIRE_ESPNOW_WIFI_POWER_AUTO
  
  if (animTaskHandleLink != NULL) {
      vTaskDelete(animTaskHandleLink);
      animTaskHandleLink = NULL;
  }

  vTaskDelay(pdMS_TO_TICKS(150));

  for (uint8_t i = 0; i < 4; i++) {
    pinMode(leds[i], OUTPUT);
    digitalWrite(leds[i], LOW); 
  }

  digitalWrite(leds[usb_data_wireless.devicePlayer -1 ], HIGH); // accende fisso il led del player corrispondente 1,2,3,4


}

#define FIFO_SIZE_READ_SER 200  // l'originale era 32
#define TIME_OUT_AVALAIBLE 2
#define TIME_OUT_SERIAL_MICRO 1000 // 1000 microsecondi = 1 millisecondo

int rx_avalaible = 0;
unsigned long startTime = 0; // = millis();

uint64_t timer_serial_micro = 0;

#ifdef DEBUG_OPENFIRE
#define TEMPO_OGNI_VISUALIZZAZIONE 1000000 // microsecondi = 1 secondo
uint64_t timer_ultimo_pacchetto_scompattato = 0;
uint64_t timer_ultimo_pacchetto_ricevuto_da_callback = 0;

uint64_t tempo_ultimo_pacchetto_scompattato = 0;
uint64_t tempo_ultimo_pacchetto_ricevuto_da_callback = 0;


uint64_t timer_ultimo_dato_visualizzato = 0;
uint64_t timer_tempo_minimo_pacchetto_ricevuto_da_callback = __UINT64_MAX__;
uint64_t timer_tempo_minimo_pacchetto_scompattato = __UINT64_MAX__;
uint64_t timer_tempo_massimo_pacchetto_ricevuto_da_callback = 0;
uint64_t timer_tempo_massimo_pacchetto_scompattato = 0;
uint16_t pacchetti_scompattati = 0;
uint16_t pacchetti_ricevuti_da_callback = 0;
#endif


//uint16_t len_aux;
uint8_t buffer_aux[FIFO_SIZE_READ_SER];
void loop()
{
  
#ifdef DEBUG_OPENFIRE
if (esp_timer_get_time() - timer_ultimo_dato_visualizzato > TEMPO_OGNI_VISUALIZZAZIONE) {
  Serial.print("tempo_minimo_pacchetto_ricevuto_da_callback: "), Serial.println(timer_tempo_minimo_pacchetto_ricevuto_da_callback);
  Serial.print("tempo_minimo_pacchetto_scompattato: "), Serial.println(timer_tempo_minimo_pacchetto_scompattato);

  Serial.print("tempo_massimo_pacchetto_ricevuto_da_callback: "), Serial.println(timer_tempo_massimo_pacchetto_ricevuto_da_callback);
  Serial.print("tempo_massimo_pacchetto_scompattato: "), Serial.println(timer_tempo_massimo_pacchetto_scompattato);

  Serial.print("pacchetti_scompattati: "), Serial.println(pacchetti_scompattati);
  Serial.print("pacchetti_ricevuti_da_callback: "), Serial.println(pacchetti_ricevuti_da_callback);

  Serial.println("================================");

  timer_ultimo_dato_visualizzato = esp_timer_get_time();
  
  //timer_ultimo_pacchetto_scompattato = 0;
  //timer_ultimo_pacchetto_ricevuto_da_callback = 0;
  
  timer_tempo_minimo_pacchetto_ricevuto_da_callback = __UINT64_MAX__;
  timer_tempo_minimo_pacchetto_scompattato = __UINT64_MAX__;
  
  timer_tempo_massimo_pacchetto_ricevuto_da_callback = 0;
  timer_tempo_massimo_pacchetto_scompattato = 0;
  
  pacchetti_scompattati = 0;
  pacchetti_ricevuti_da_callback = 0;
}
#endif
vTaskDelay(pdMS_TO_TICKS(1));
rx_avalaible = Serial.available();
if (rx_avalaible > FIFO_SIZE_READ_SER) rx_avalaible = FIFO_SIZE_READ_SER;
if (rx_avalaible)
{
    Serial.readBytes(buffer_aux, rx_avalaible);
    SerialWireless.write(buffer_aux, rx_avalaible);
    SerialWireless.flush(); // provare a togliere
}
//yield();
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef COMMENTO
  
  //micros();
  //esp_timer_get_time();
  if (esp_timer_get_time() - timer_serial_micro > TIME_OUT_SERIAL_MICRO) // controlla ogni millisecondo più o meno a 9600 bps
  {
  if (Serial.available() > rx_avalaible) {
    startTime = millis();
    rx_avalaible = Serial.available();
  }     
  if (rx_avalaible && (millis() > startTime + TIME_OUT_AVALAIBLE)) { 
    if (rx_avalaible > FIFO_SIZE_READ_SER) rx_avalaible = FIFO_SIZE_READ_SER;
    Serial.readBytes(buffer_aux, rx_avalaible);
    SerialWireless.write(buffer_aux, rx_avalaible);
    //SerialWireless.lenBufferSerialWrite = rx_avalaible;
    //SerialWireless.flush();
    rx_avalaible = 0;
  } 
  timer_serial_micro = esp_timer_get_time();
}
#endif


///////////////////////////////////////////////////////////////////////////////////////////////////////////


  /*
  len_aux = Serial.available();
  if (len_aux) {
    if (len_aux > FIFO_SIZE_WRITE_SERIAL) len_aux = FIFO_SIZE_WRITE_SERIAL;
    Serial.readBytes(buffer_aux, len_aux);
    SerialWireless.write(buffer_aux,len_aux);
  }
  */

  /*
  SerialWireless.lenBufferSerialWrite = Serial.available();
  if (SerialWireless.lenBufferSerialWrite) {
    if (SerialWireless.lenBufferSerialWrite > FIFO_SIZE_WRITE_SERIAL) SerialWireless.lenBufferSerialWrite = FIFO_SIZE_WRITE_SERIAL;
    Serial.readBytes(SerialWireless.bufferSerialWrite, SerialWireless.lenBufferSerialWrite);
    SerialWireless.flush();
  }
  */

#ifdef COMMENTO  
  if (Serial.available() > rx_avalaible) {
    startTime = millis();
    rx_avalaible = Serial.available();
  }     
  if (rx_avalaible && (millis() > startTime + TIME_OUT_AVALAIBLE)) { 
    if (rx_avalaible > FIFO_SIZE_READ_SER) rx_avalaible= FIFO_SIZE_READ_SER;
    Serial.readBytes(SerialWireless.bufferSerialWrite, rx_avalaible);
    SerialWireless.lenBufferSerialWrite = rx_avalaible;
    SerialWireless.flush();
    rx_avalaible = 0;
  } 
#endif //COMMENTO
}