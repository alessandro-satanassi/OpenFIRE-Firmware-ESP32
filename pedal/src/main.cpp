/*!
 * @file main.cpp
 * @brief main PEDAL
 * @n CPP main PEDAL
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */


#include <Arduino.h>

#include "TinyUSB_Devices.h"

#include "OpenFIRE-PEDAL-version.h"

// ===================================================================================
// HARDWARE STATE & CONFIGURATION
// ===================================================================================

bool display_init = false;

const int8_t leds[4] = {PIN_LED1, 
                        PIN_LED2, 
                        PIN_LED3, 
                        PIN_LED4}; // The GPIOs to which the 4 LEDs are connected / I GPIO ove sono collegati i 4 LED

/*
We maintain the state of all pedals in a single byte (bitmask).
This allows the entire module state to be transmitted using a single radio payload byte,
minimizing network overhead and airtime on ESP-NOW.
/
Manteniamo lo stato di tutti i pedali in un singolo byte (maschera di bit). 
Questo permette di trasmettere l'intero stato del modulo con un solo byte di payload radio, 
minimizzando l'overhead di rete e i tempi di volo (Air Time) su ESP-NOW.
*/
uint8_t buttons_state = 0;
volatile uint8_t buttons_last_state = 0;

// ===================================================================================
// Keep-alive and debouncing management / GESTIONE KEEP-ALIVE E DEBOUNCING
// ===================================================================================

volatile bool send_packet_pedal = false; // true if a packet needs to be sent / true se bisogna spedire un pacchetto
#define DEBOUNCE_DELAY 15 // Debounce time for buttons in ms (perhaps 20ms is safer) / tempo di deboincing per pulsanti in ms (forse 20ms è più sicuro)

/*
When a pedal is held down for an extended period (e.g., to take cover in Time Crisis),
the dongle might lose the connection or interpret the lack of signal as a disconnection.
This timer ensures the active state is periodically re-sent (keep-alive).
NOTE: The value is transmitted every TIME_REPEAT_SEND only if a button is pressed,
to minimize collisions. The value must be specified in microseconds (i.e., ms x 1000).
/
Quando un pedale viene tenuto premuto a lungo (es. per ripararsi in Time Crisis), 
il dongle potrebbe perdere la connessione o interpretare il silenzio radio come disconnessione.
Questo timer garantisce un reinvio periodico (keep-alive) dello stato attivo.
NOTA: il valore è trasmesso ogni TIME_REPEAT_SEND solo se qualche buttone è premuto 
per ridurre al minimo le collisioni. Il valore va specificato in microsecondi quindi ms x 1000
*/
#define TIME_REPEAT_SEND (uint64_t)69000 // time in milliseconds (set to 69ms, perhaps better at 50ms) / tempo in microsendi (impostato a 69ms, forse meglio a 50ms)

#define NUM_BUTTONS 2

/*
The Button struct encapsulates the state machine for software debouncing.
Using this struct allows you to scale the number of pedals without having to
duplicate the logic in the main loop.
/
La struttura Button incapsula la Macchina a Stati per il debouncing software.
L'uso di questa struct permette di scalare il numero di pedali senza dover 
duplicare la logica nel loop principale.
*/
struct Button {
    int8_t pin;
    bool currentState;
    bool lastState;
    uint8_t debounceTime; // 0 means that it does not require debouncing. / 0 vuole dire che non abbisogna di deboucing
    unsigned long lastDebounceTime;
};

Button buttons[NUM_BUTTONS] = {
    {PIN_PEDAL, LOW, LOW, 0, 0}, // Pedal 1 / Pedale 1
    {PIN_PEDAL2, LOW, LOW, 0, 0}  // Pedal 2 / Pedale 2
};


/*
===================================================================================
Visual feedback via the 4 LEDs / FEEDBACK VISIVO CON I 4 LED
===================================================================================
Wireless layer initialization is a blocking operation and can take several seconds.
We leverage the ESP32's dual-core processor by delegating an independent task to
a FreeRTOS thread. This provides immediate visual feedback (Knight Rider-style animation)
to the user, confirming that the pedal is active and searching for the network.
/
L'inizializzazione del layer Wireless è bloccante e può richiedere secondi.
Sfruttiamo il processore dual-core dell'ESP32 delegando un task indipendente a 
un FreeRTOS thread. Questo fornisce un feedback visivo immediato (animazione Supercar)
all'utente, confermando che il pedale è vivo e sta cercando la rete.
*/

void animTaskLink(void *pvParameters) {
   
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

/*
===============================================================================
AUTOMATIC RESEND TIMER (KEEP-ALIVE) / TIMER REINVIO AUTOMATICO (KEEP-ALIVE)
===============================================================================
CALLBACK
This callback comes in a response Hardware Timer interrupt (or task ad
high priority). It must be extremely brief and non-blocking.
We simply set the atomic flag `send_packet_pedal`.
/
CALLBACK
Questa callback viene eseguita in un contesto Hardware Timer interrupt (o task ad 
alta priorità). Deve essere estremamente breve e non bloccante. 
Alziamo semplicemente la flag atomica `send_packet_pedal`.
*/
void timer_callback_send_repeat(void* arg) {
  
  if (buttons_last_state) send_packet_pedal = true;
 
}

esp_timer_handle_t timer_handle_send_repeat;

// ===============================================================================
// MAIN SETUP
// ===============================================================================

// The main show!
void setup() {

  // Pedal Configuration (Inputs with Pull-Up) / Configurazione Pedali (Ingressi con Pull-Up)
  for (uint8_t i = 0; i < NUM_BUTTONS; i++) {
    pinMode(buttons[i].pin, INPUT_PULLUP);
  }

  for (uint8_t i = 0; i < 4; i++) {
    pinMode(leds[i], OUTPUT);
    digitalWrite(leds[i], LOW); 
  }

  // Creation of the asynchronous animation task. / Creazione del task visivo asincrono.
  TaskHandle_t animTaskHandleLink = NULL;  
  xTaskCreatePinnedToCore(
        animTaskLink,          // task function / funzione del task
        "AnimTaskLink",        // name / nome
        4096,                  // stack size
        NULL,                  // parameters / parametri
        1,                     // priority / priorità
        &animTaskHandleLink,   // handle
        APP_CPU_NUM            // core (you can use 0 or 1) / (puoi usare 0 o 1)
  );  
  
  // ====== wireless connection management / gestione connessione wireless ======
  // Blocking phase: The pedal waits here until the ESP-NOW handshake is established. / Fase bloccante: Il pedale attende qui finché non stabilisce l'handshake ESP-NOW.
  SerialWireless.init_wireless();
  SerialWireless.begin();
  SerialWireless.connection_pedal();
  // ====== Wireless management complete... it proceeds only after the device has been paired. / fine gestione wireless .. va avanti solo dopo che si è accoppiato il dispositivo ======
  
  /*
  Connection established. Terminating the animation task, as control of the LEDs is now passed to the system logic.
  /
  Connessione stabilita. Uccidiamo il task visivo poiché ora il controllo dei LED passa alla logica di sistema.
  */
  if (animTaskHandleLink != NULL) {
      vTaskDelete(animTaskHandleLink);
      animTaskHandleLink = NULL;
  }

  vTaskDelay(pdMS_TO_TICKS(150));

  for (uint8_t i = 0; i < 4; i++) {
    digitalWrite(leds[i], LOW); 
  }

  // Visual confirmation of successful mating. / Riscontro visivo di accoppiamento riuscito.
  if (usb_data_wireless.devicePlayer >= 1 && usb_data_wireless.devicePlayer <= 4) {
    digitalWrite(leds[usb_data_wireless.devicePlayer - 1], HIGH); // turns on the corresponding player LED 1,2,3,4 / accende fisso il led del player corrispondente 1,2,3,4
  } else {
    // Error handling: turns on the first and last LEDs. / Gestione errore: accende il primo e l'ultimo led
    digitalWrite(leds[0], HIGH);
    digitalWrite(leds[3], HIGH);
  }
  
  // Hardware timer initialization for keep-alive. / Inizializzazione Timer Hardware per il keep-alive.
  esp_timer_create_args_t timer_args = {
    .callback = &timer_callback_send_repeat,
    .arg = NULL,
    .dispatch_method = ESP_TIMER_TASK,
    .name = "timer_send_repeat"
  };
    
  esp_timer_create(&timer_args, &timer_handle_send_repeat);
  esp_timer_start_periodic(timer_handle_send_repeat, TIME_REPEAT_SEND);
 
}

// ===============================================================================
// MAIN LOOP
// ===============================================================================

void loop()
{
  unsigned long millis_current = millis(); 
  
  uint8_t bitMask = 1;
  buttons_state = 0;

  // Non-blocking pedal scanning and debouncing / SCANSIONE PEDALI E DEBOUNCING NON-BLOCCANTE
  for (uint8_t i = 0; i < NUM_BUTTONS; i++, bitMask <<= 1) {
    if (!buttons[i].debounceTime) {
      // Pedal debounced (or steady state reached). Reading the hardware. / Pedale a riposo (o stato stabile raggiunto). Leggiamo l'hardware.
      buttons[i].currentState = buttons[i].pin >= 0 ? !(bool)digitalRead(buttons[i].pin) : false;
      
      if (buttons[i].currentState) {
        buttons_state |= bitMask; 
      }
      /*
      If there is a physical change, we trigger the debounce window. 
      During this period we will ignore the switch's bounces.
      /
      Se c'è un cambiamento fisico, inneschiamo la finestra di debounce.
      In questo periodo di "cecità" volontaria ignoreremo i rimbalzi meccanici dello switch.
      */
      if (buttons[i].currentState != buttons[i].lastState) {
        buttons[i].debounceTime = DEBOUNCE_DELAY;
        buttons[i].lastDebounceTime = millis_current; 
        buttons[i].lastState = buttons[i].currentState;
      }
    }
    else {
      // Debounce window active. We reduce the timer by scaling the elapsed time delta (`aux`). / Finestra di debounce attiva. Riduciamo il timer scalando il delta temporale trascorso (`aux`).
      unsigned long aux = millis_current - buttons[i].lastDebounceTime;
      
      if (aux >= buttons[i].debounceTime) {
          buttons[i].debounceTime = 0;
      } else {
          buttons[i].debounceTime -= aux; 
          buttons[i].lastDebounceTime = millis_current;
      }
      
      // Store the button state for debounce checking. / Manteniamo lo stato consolidato in attesa che il rimbalzo hardware finisca.
      if (buttons[i].lastState) {
        buttons_state |= bitMask; 
      }  
    }
  }
  
  /*
  EVENT-DRIVEN RADIO TRANSMISSION LOGIC
  We transmit over the network EXCLUSIVELY in two scenarios:
  A) A physical state change has occurred (pedal press/release). 
  B) The hardware timer has raised the keep-alive flag (pedal held down).
  /
  LOGICA EVENT-DRIVEN DI TRASMISSIONE RADIO
  Trasmettiamo sulla rete ESCLUSIVAMENTE in due scenari:
  A) È avvenuto un cambio di stato fisico (salita/discesa di un pedale).
  B) Il timer hardware ha alzato la flag per il keep-alive (pedale mantenuto premuto).
  */
  if ((buttons_state != buttons_last_state) || (send_packet_pedal == true)) {
    // Send packet with pedal position to the lightgun - send the buttons_state byte / Invia pacchetto con posizione pedali alla lightgun - invia il byte buttons_state
 
    SerialWireless.SendPacket((const uint8_t *)&buttons_state, 1, PACKET_TX::PEDAL_TX);

    buttons_last_state = buttons_state;
     
    esp_timer_restart(timer_handle_send_repeat, TIME_REPEAT_SEND);
     
    send_packet_pedal = false;

  }
      
  // Voluntary delay to yield CPU time (e.g. to WiFi/ESP-NOW tasks) and rate-limit this loop. / Rate limiting volontario per cedere CPU al layer WiFi/ESP-NOW sottostante.
  vTaskDelay(pdMS_TO_TICKS(5));  // Approximate 5 ms polling interval. / fai polling ogni 5ms
}