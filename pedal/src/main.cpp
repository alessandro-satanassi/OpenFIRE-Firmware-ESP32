// by SATANASSI Alessandro

#include <Arduino.h>

#include "TinyUSB_Devices.h"

#include "OpenFIRE-PEDAL-version.h"

bool display_init = false;

const int8_t leds[4] = {PIN_LED1, 
                        PIN_LED2, 
                        PIN_LED3, 
                        PIN_LED4}; // I tuoi 4 GPIO


//uint8_t buffer_aux[5];                        volatile uint8_t buttons_state = 0;
uint8_t buttons_state = 0;
volatile uint8_t buttons_last_state = 0;
volatile bool send_packet_pedal = false; // true se bisogna spedire un pacchetto
#define DEBOUNCE_DELAY 15 // o meglio 20 ms ??// tempo di deboincing per pulsanti
#define TIME_REPEAT_SEND (uint64_t)69000 // 50 ms // reinvias il pacchetto ogni tot ms solo se qualche buttone è premuto .. il valore va specificato in microsecondi quidi ms x 1000
#define NUM_BUTTONS 2

struct Button {
    int8_t pin;
    bool currentState;
    bool lastState;
    uint8_t debounceTime; // 0 vuole dire che non abbisogna di deboucing
    unsigned long lastDebounceTime;
};

Button buttons[NUM_BUTTONS] = {
    {PIN_PEDAL, LOW, LOW, 0, 0}, // Pedale 1
    {PIN_PEDAL2, LOW, LOW, 0, 0}  // Pedale 2
};

void animTaskLink(void *pvParameters) {
  
  /*
  for (uint8_t i = 0; i < 4; i++) {
    pinMode(leds[i], OUTPUT);
    digitalWrite(leds[i], LOW); 
  }
  */
  
  int currentLed = 0;
  int direction = 1;

  
  for (;;) {
    digitalWrite(leds[currentLed], HIGH);
    vTaskDelay(pdMS_TO_TICKS(150)); 
    digitalWrite(leds[currentLed], LOW);
    currentLed += direction;
    if (currentLed >= 3) {
      direction = -1;
    } else if (currentLed <= 0) {
      direction = 1;
    }
  } 
}

// ================================== TIMER REINVIO AUTOMATICO ===================================
// CALLBACK
void timer_callback_send_repeat(void* arg) {
  
  if (buttons_last_state) send_packet_pedal = true;
 
}

esp_timer_handle_t timer_handle_send_repeat;
// ===============================================================================

// The main show!
void setup() {
  //TinyUSBDevices.begin(1);
  //TinyUSBDevice.begin(0);
  //Serial.begin(115200);

  // Configurazione Pedali (Ingressi con Pull-Up)
  for (uint8_t i = 0; i < NUM_BUTTONS; i++) {
    pinMode(buttons[i].pin, INPUT_PULLUP);
  }

  for (uint8_t i = 0; i < 4; i++) {
    pinMode(leds[i], OUTPUT);
    digitalWrite(leds[i], LOW); 
  }


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


  #ifdef COMMENTO
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
  #endif // COMMENTO


  //Serial.begin(115200);
  //Serial.printf("\n\nWireless channel da USB_DATA = %d\n\n", usb_data_wireless.channel);
  //Serial.printf("\n\nWireless channel da espnow_wifi_channel = %d\n\n", espnow_wifi_channel);
  
  /*
  Serial.begin(9600);
  Serial.setTimeout(0);
  */
    
  if (animTaskHandleLink != NULL) {
      vTaskDelete(animTaskHandleLink);
      animTaskHandleLink = NULL;
  }

  vTaskDelay(pdMS_TO_TICKS(150));

  //if ((usb_data_wireless.channel == espnow_wifi_channel) && (espnow_wifi_channel == 1))  {
  //SerialWireless.mac_esp_another_card, 6) != 0 && // controlla che sia nei peer la board da cui arrivano i paccheti
  //    memcmp(info->des_addr, peerAddress
  //if ((peerAddress[0] == 0xFF))  {  
  if (esp_now_is_peer_exist(peerAddress) && !memcmp(peerAddress, SerialWireless.mac_esp_another_card, 6)) {
    digitalWrite(leds[0], HIGH);
    digitalWrite(leds[1], HIGH);
    digitalWrite(leds[2], HIGH);
    digitalWrite(leds[3], HIGH);
    vTaskDelay(pdMS_TO_TICKS(1500));

  }


  for (uint8_t i = 0; i < 4; i++) {
    //pinMode(leds[i], OUTPUT);
    digitalWrite(leds[i], LOW); 
  }

  
  if (usb_data_wireless.devicePlayer >= 1 && usb_data_wireless.devicePlayer <= 4) {
    digitalWrite(leds[usb_data_wireless.devicePlayer - 1], HIGH); // accende fisso il led del player corrispondente 1,2,3,4
  } else {
    // Gestione errore: accende il primo e l'ultimo led
    digitalWrite(leds[0], HIGH);
    digitalWrite(leds[3], HIGH);
  }

  //vTaskDelay(pdMS_TO_TICKS(2000)); 

  
  esp_timer_create_args_t timer_args = {
    .callback = &timer_callback_send_repeat,
    .arg = NULL,
    .dispatch_method = ESP_TIMER_TASK,
    .name = "timer_send_repeat"
  };
    
  esp_timer_create(&timer_args, &timer_handle_send_repeat);
  esp_timer_start_periodic(timer_handle_send_repeat, TIME_REPEAT_SEND);
  
  //vTaskDelay(pdMS_TO_TICKS(5000));

}


void loop()
{
  //vTaskDelay(pdMS_TO_TICKS(5000));

  unsigned long millis_current = millis(); 
  
  uint8_t bitMask = 1;
  buttons_state = 0;

  for (uint8_t i = 0; i < NUM_BUTTONS; i++, bitMask <<= 1) {
    if (!buttons[i].debounceTime) {
      buttons[i].currentState = buttons[i].pin >= 0 ? !(bool)digitalRead(buttons[i].pin) : false;
      //gpio_get_level((gpio_num_t)buttons[i].pin);
      
      if (buttons[i].currentState) {
        buttons_state |= bitMask; 
      }
      
      if (buttons[i].currentState != buttons[i].lastState) {
        buttons[i].debounceTime = DEBOUNCE_DELAY;
        buttons[i].lastDebounceTime = millis_current; 
        buttons[i].lastState = buttons[i].currentState;
      }
    }
    else {
      unsigned long aux = millis_current - buttons[i].lastDebounceTime;
      
      if (aux >= buttons[i].debounceTime) {
          buttons[i].debounceTime = 0;
      } else {
          buttons[i].debounceTime -= aux; 
          buttons[i].lastDebounceTime = millis_current;
      }
      
      if (buttons[i].lastState) {
        buttons_state |= bitMask; 
      }  
    }
  }

  //vTaskDelay(pdMS_TO_TICKS(50000));

 
  if ((buttons_state != buttons_last_state) || (send_packet_pedal == true)) {
    // Invia pacchetto con posizione pedali alla lightgun - invia il byte buttons_state
    //vTaskDelay(pdMS_TO_TICKS(5000));
    //buffer_aux[0] = buttons_state;
    
    //SerialWireless.SendPacket((const uint8_t *)buffer_aux, 1, PACKET_TX::PEDAL_TX);
    SerialWireless.SendPacket((const uint8_t *)&buttons_state, 1, PACKET_TX::PEDAL_TX);

    buttons_last_state = buttons_state;
    
    
    esp_timer_restart(timer_handle_send_repeat, TIME_REPEAT_SEND);
     
    send_packet_pedal = false;

    #ifndef COMMENTO
     for (uint8_t i = 0; i < 4; i++) {
    //pinMode(leds[i], OUTPUT);
    digitalWrite(leds[i], LOW); 
  }

    if (buttons_state == 1 || buttons_state == 3) digitalWrite(leds[0], HIGH);
    if (buttons_state == 2 || buttons_state == 3)digitalWrite(leds[3], HIGH);
    vTaskDelay(pdMS_TO_TICKS(50));
    digitalWrite(leds[0], LOW);
    digitalWrite(leds[3], LOW);

    digitalWrite(leds[usb_data_wireless.devicePlayer - 1], HIGH);
   #endif // COMMENTO
  }
  

  //vTaskDelay(pdMS_TO_TICKS(5000));
  
  vTaskDelay(pdMS_TO_TICKS(5));  // fai polling ogni 5ms
}