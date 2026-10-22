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

#include <WiFi.h>
#include <esp_sleep.h>
#include <driver/gpio.h>

//#include "OpenFIRE-PEDAL-version.h"

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
// Repeats of the released state still to send: written only by the loop, the timer only reads it. / Ripetizioni dello stato rilasciato ancora da inviare: la scrive solo il loop, il timer la legge soltanto.
#define RELEASE_REPEATS 2
volatile uint8_t release_repeats = 0;
#define DEBOUNCE_DELAY 15 // Debounce time for buttons in ms (perhaps 20ms is safer) / tempo di deboincing per pulsanti in ms (forse 20ms è più sicuro)

/*
When a pedal is held down for an extended period (e.g., to take cover in Time Crisis),
the dongle might lose the connection or interpret the lack of signal as a disconnection.
This timer ensures the active state is periodically re-sent (keep-alive).
NOTE: The value is transmitted every TIME_REPEAT_SEND only if a button is pressed,
to minimize collisions. After a release, the released state is also sent RELEASE_REPEATS more times,
so that a lost release packet is recovered without waiting for the receiver's timeout. The value must be specified in microseconds (i.e., ms x 1000).
/
Quando un pedale viene tenuto premuto a lungo (es. per ripararsi in Time Crisis), 
il dongle potrebbe perdere la connessione o interpretare il silenzio radio come disconnessione.
Questo timer garantisce un reinvio periodico (keep-alive) dello stato attivo.
NOTA: il valore è trasmesso ogni TIME_REPEAT_SEND solo se qualche buttone è premuto 
per ridurre al minimo le collisioni. Dopo un rilascio, lo stato rilasciato viene inviato ancora RELEASE_REPEATS volte,
così un pacchetto di rilascio perso si recupera senza attendere il timeout del ricevitore. Il valore va specificato in microsecondi quindi ms x 1000
*/
#define TIME_REPEAT_SEND (uint64_t)69000 // time in milliseconds (set to 69ms, perhaps better at 50ms) / tempo in microsendi (impostato a 69ms, forse meglio a 50ms)

#define NUM_BUTTONS 2

/*
===================================================================================
POWER SAVING / RISPARMIO ENERGETICO
===================================================================================
If the pedal is left on, it saves the battery.
- Paired with a gun: after SLEEP_MINUTES from the last pedal press or release, the radio switches off and the player LED
  blinks. A press switches the radio back on and the pedal carries on with the same gun, without pairing again.
- Waiting for a gun: the search lasts SLEEP_SEARCH_MINUTES, then the radio switches off and LEDs 1 and 4 blink.
  A press restarts the pedal, which searches again.
0 in SLEEP_MINUTES disables the power saving once paired; 0 in SLEEP_SEARCH_MINUTES makes the search endless.
/
Se il pedale resta acceso, risparmia la batteria.
- Associato a una pistola: dopo SLEEP_MINUTES dall'ultima pressione o rilascio di un pedale, la radio si spegne e lampeggia
  il LED del giocatore. Una pressione riaccende la radio e il pedale continua con la stessa pistola, senza riassociarsi.
- In attesa di una pistola: la ricerca dura SLEEP_SEARCH_MINUTES, poi la radio si spegne e lampeggiano i LED 1 e 4.
  Una pressione riavvia il pedale, che ricomincia a cercare.
0 in SLEEP_MINUTES disattiva il risparmio energetico da associato; 0 in SLEEP_SEARCH_MINUTES rende la ricerca senza fine.
*/
#define SLEEP_MINUTES 15 // 1 // 30          // minutes without presses before power saving / minuti senza pressioni prima del risparmio energetico
#define SLEEP_SEARCH_MINUTES 5 // 1 // 10   // minutes of search for a gun before power saving / minuti di ricerca della pistola prima del risparmio energetico
#define SLEEP_BLINK_PERIOD_MS 1000 // 4000 // one blink every / un lampeggio ogni
#define SLEEP_BLINK_ON_MS 20 // 30       // how long the LEDs stay on in each blink / per quanto restano accesi i LED a ogni lampeggio

/*
CPU frequency: 80 MHz is the lowest the radio (Wi-Fi/ESP-NOW) supports. The pedal only reads two switches and sends a few
packets per second, so it has plenty of margin and saves about 30 mA from the battery compared with 240 MHz.
If start-up or pairing ever misbehave, try 160.
/
Frequenza della CPU: 80 MHz è la minima ammessa dalla radio (Wi-Fi/ESP-NOW). Il pedale legge solo due interruttori e invia
pochi pacchetti al secondo, quindi ha margine in abbondanza e risparmia circa 30 mA dalla batteria rispetto a 240 MHz.
Se avvio o accoppiamento dovessero dare problemi, provare 160.
*/
#define PEDAL_CPU_FREQ_MHZ 80

#define SLEEP_MS ((unsigned long)SLEEP_MINUTES * 60000UL)
#define SLEEP_SEARCH_MS ((unsigned long)SLEEP_SEARCH_MINUTES * 60000UL)

unsigned long last_activity_millis = 0; // last press or release of a pedal / ultima pressione o rilascio di un pedale

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
  
  if (buttons_last_state || release_repeats) send_packet_pedal = true;
 
}

esp_timer_handle_t timer_handle_send_repeat;

// Steady LEDs while connected: the player LED, or the first and last LED if the player is not valid. / LED fissi da connesso: il LED del giocatore, oppure il primo e l'ultimo se il giocatore non è valido.
uint8_t player_leds_mask() {
  if (usb_data_wireless.devicePlayer >= 1 && usb_data_wireless.devicePlayer <= 4) return 1 << (usb_data_wireless.devicePlayer - 1);
  return 0b1001;
}

void set_leds(uint8_t mask) {
  for (uint8_t i = 0; i < 4; i++) {
    if (leds[i] >= 0) digitalWrite(leds[i], (mask >> i) & 1 ? HIGH : LOW);
  }
}

// Light sleep for the given time; returns true if a pedal press woke the pedal. / Sonno leggero per il tempo indicato; restituisce true se a svegliare il pedale è stata una pressione.
bool light_sleep_ms(uint32_t ms) {
  esp_sleep_enable_timer_wakeup((uint64_t)ms * 1000);
  esp_light_sleep_start();
  return esp_sleep_get_wakeup_cause() == ESP_SLEEP_WAKEUP_GPIO;
}

/*
POWER SAVING LOOP / CICLO DI RISPARMIO ENERGETICO
The chip stays in light sleep and wakes up only to blink the LEDs, or when a pedal is pressed.
With a gun paired it returns at the first press; while waiting for a gun a press restarts the pedal.
/
Il chip resta in sonno leggero e si sveglia solo per far lampeggiare i LED, o quando si preme un pedale.
Con una pistola associata ritorna alla prima pressione; in attesa di una pistola una pressione riavvia il pedale.
*/
void power_save(bool paired) {
  if (paired) esp_timer_stop(timer_handle_send_repeat);
  SerialWireless.radio_off_pedal();
  set_leds(0);

  // Pins keep their normal setting during light sleep: LEDs stay lit for the blink, pedals keep their pull-up. / I pin mantengono la configurazione normale nel sonno leggero: i LED restano accesi per il lampeggio, i pedali il pull-up.
  for (uint8_t i = 0; i < 4; i++) {
    if (leds[i] >= 0) gpio_sleep_sel_dis((gpio_num_t)leds[i]);
  }
  for (uint8_t i = 0; i < NUM_BUTTONS; i++) {
    if (buttons[i].pin >= 0) {
      gpio_sleep_sel_dis((gpio_num_t)buttons[i].pin);
      gpio_wakeup_enable((gpio_num_t)buttons[i].pin, GPIO_INTR_LOW_LEVEL); // pressed pedal = pin to ground / pedale premuto = pin a massa
    }
  }
  esp_sleep_enable_gpio_wakeup();

  for (;;) {
    set_leds(paired ? player_leds_mask() : 0b1001); // player LED, or LEDs 1 and 4 while waiting for a gun / LED del giocatore, oppure LED 1 e 4 in attesa di una pistola
    bool pressed = light_sleep_ms(SLEEP_BLINK_ON_MS);
    set_leds(0);
    if (pressed || light_sleep_ms(SLEEP_BLINK_PERIOD_MS - SLEEP_BLINK_ON_MS)) break;
  }

  // Woken by a pedal. / Svegliato da un pedale.
  if (!paired) esp_restart();

  for (uint8_t i = 0; i < NUM_BUTTONS; i++) {
    if (buttons[i].pin >= 0) gpio_wakeup_disable((gpio_num_t)buttons[i].pin);
  }
  esp_sleep_disable_wakeup_source(ESP_SLEEP_WAKEUP_ALL);

  SerialWireless.radio_on_pedal();
  set_leds(player_leds_mask());
  last_activity_millis = millis();
  esp_timer_start_periodic(timer_handle_send_repeat, TIME_REPEAT_SEND);
}

// ===============================================================================
// MAIN SETUP
// ===============================================================================

// The main show!
void setup() {

  #ifdef PEDAL_CPU_FREQ_MHZ
    setCpuFrequencyMhz(PEDAL_CPU_FREQ_MHZ);
  #endif // PEDAL_CPU_FREQ_MHZ

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
  // The Arduino core applies its own Wi-Fi power save at every radio start (default: modem sleep); the pedal wants none, also after waking up. / Il core Arduino applica il suo risparmio Wi-Fi a ogni avvio della radio (predefinito: modem sleep); il pedale non ne vuole, anche dopo il risveglio.
  WiFi.setSleep(WIFI_PS_NONE);
  SerialWireless.init_wireless();
  SerialWireless.begin();
  bool paired = SerialWireless.connection_pedal(SLEEP_SEARCH_MS); // after SLEEP_SEARCH_MINUTES without a gun it gives up / dopo SLEEP_SEARCH_MINUTES senza pistola rinuncia
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

  if (!paired) {
    // No gun found: power saving straight away (LEDs 1 and 4 blink); a press restarts the pedal. / Nessuna pistola trovata: subito risparmio energetico (lampeggiano i LED 1 e 4); una pressione riavvia il pedale.
    power_save(false);
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

  last_activity_millis = millis();
 
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
    // The current state is always sent, never a stored one. / Si invia sempre lo stato corrente, mai uno memorizzato.
    if (buttons_state) release_repeats = 0;                                      // Pressed: no release repeats. / Premuto: nessuna ripetizione del rilascio.
    else if (buttons_last_state) release_repeats = RELEASE_REPEATS;              // Just released: schedule the repeats. / Appena rilasciato: prepara le ripetizioni.
    else if (release_repeats) release_repeats = release_repeats - 1;             // Repeat of the released state. / Ripetizione dello stato rilasciato.
 
    SerialWireless.SendPacket((const uint8_t *)&buttons_state, 1, PACKET_TX::PEDAL_TX);

    buttons_last_state = buttons_state;
     
    esp_timer_restart(timer_handle_send_repeat, TIME_REPEAT_SEND);
     
    send_packet_pedal = false;

  }

  // A pedal held down counts as activity too. / Anche un pedale tenuto premuto conta come attività.
  if (buttons_state || release_repeats) last_activity_millis = millis_current;
  else if ((SLEEP_MS > 0) && ((millis_current - last_activity_millis) >= SLEEP_MS)) power_save(true);
      
  // Voluntary delay to yield CPU time (e.g. to WiFi/ESP-NOW tasks) and rate-limit this loop. / Rate limiting volontario per cedere CPU al layer WiFi/ESP-NOW sottostante.
  vTaskDelay(pdMS_TO_TICKS(5));  // Approximate 5 ms polling interval. / fai polling ogni 5ms
}