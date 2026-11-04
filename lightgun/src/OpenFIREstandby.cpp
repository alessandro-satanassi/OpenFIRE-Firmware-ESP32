/*!
 * @brief Battery ESP-NOW standby / Standby ESP-NOW a batteria.
 * @copyright Alessandro Satanassi, 2026, GNU Lesser General Public License
 *
 * EN: After the standby minutes (see OpenFIREstandby.h) without ESP-NOW transmission attempts
 * the gun turns the radio off and light-sleeps, blinking faintly. The trigger wakes it:
 * same dongle, no new pairing, and the wake-up pull is not sent as a shot.
 * IT: Dopo i minuti di standby (vedi OpenFIREstandby.h) senza tentativi di trasmissione ESP-NOW
 * la pistola spegne la radio e va in light sleep, con un lampeggio tenue. Il grilletto la
 * risveglia: stesso dongle, nessun nuovo abbinamento, e il tiro di risveglio non e' uno sparo.
 *
 * EN: How it works:
 *  1. Begin() creates a one-minute esp_timer, only on battery with ESP-NOW to the dongle.
 *  2. The timer callback only raises standbyRequested, when no transmission was attempted for
 *     the standby minutes (lastRadioTxMs, updated by OpenFIRE_Wireless at every send).
 *  3. The main core sees the flag after GetPosition() and calls Check(), which does all the work:
 *     it pauses the input core, checks that nothing is pending, turns the radio off and sleeps.
 *  4. With DUAL_CORE the input core stays blocked in WaitIfPaused() while the radio is off.
 * IT: Come funziona:
 *  1. Begin() crea un esp_timer di un minuto, solo a batteria con ESP-NOW verso il dongle.
 *  2. La callback del timer alza soltanto standbyRequested, quando non ci sono stati tentativi di
 *     trasmissione per i minuti di standby (lastRadioTxMs, aggiornato da OpenFIRE_Wireless a ogni invio).
 *  3. Il core principale vede il flag dopo GetPosition() e chiama Check(), che fa tutto il lavoro:
 *     ferma il core degli input, controlla che non ci sia nulla in sospeso, spegne la radio e dorme.
 *  4. Con DUAL_CORE il core degli input resta bloccato in WaitIfPaused() finche' la radio e' spenta.
 */
#include "OpenFIREstandby.h"
#if defined(ARDUINO_ARCH_ESP32) && defined(OPENFIRE_WIRELESS_ENABLE)

#include "OpenFIREcommon.h"
#include "OpenFIREFeedback.h"
#include "OpenFIREserial.h"
#include "OpenFIRElights.h"
#include "OpenFIREweb.h"
#include <esp_sleep.h>
#include <esp_timer.h>
#include <driver/gpio.h>

// EN/IT: Defined in OpenFIRE_Wireless.cpp, not declared in its header / definito in OpenFIRE_Wireless.cpp, non dichiarato nel suo header.
extern volatile bool radioFree;
#ifdef DUAL_CORE
// EN/IT: Input task (loop1), created in main.cpp / task degli input (loop1), creato in main.cpp.
extern TaskHandle_t task_loop1;
#endif

volatile bool OF_Standby::standbyRequested = false;

namespace {
constexpr uint64_t CheckPeriodUs = 60ULL * 1000000ULL;
esp_timer_handle_t checkTimer = nullptr;
// EN/IT: Standby sign, a short blink every few seconds / segno di standby, breve lampeggio ogni pochi secondi.
constexpr uint32_t BlinkMs = 40, BlinkPeriodMs = 3000;
// EN/IT: Longest wait for the input core to stop / attesa massima perche' il core dei tasti si fermi.
constexpr uint32_t PauseWaitMs = 100;

#ifdef DUAL_CORE
// EN: Handshake between the two cores. The main core sets pauseInputs; the input core answers with
//     inputsPaused when it is blocked in WaitIfPaused(). Every write to these flags is followed by
//     __sync_synchronize(), a full memory barrier (memw on the ESP32-S3): the write is visible to the
//     other core before anything that follows, and later reads are not moved before it.
// IT: Accordo fra i due core. Il core principale alza pauseInputs; il core degli input risponde con
//     inputsPaused quando e' bloccato in WaitIfPaused(). Ogni scrittura su questi flag e' seguita da
//     __sync_synchronize(), una barriera di memoria completa (memw sull'ESP32-S3): la scrittura e' visibile
//     all'altro core prima di tutto quello che segue, e le letture successive non vengono anticipate.
volatile uint32_t pauseInputs = 0, inputsPaused = 0;
#endif

// EN: Read at every check, so a new saved value applies without restarting. 0 = standby disabled.
//     The saved value is validated only here: anything outside standbyMinutesChoices means the default.
// IT: Letto a ogni controllo, cosi' un nuovo valore salvato vale senza riavvio. 0 = standby disabilitato.
//     Il valore salvato e' validato solo qui: tutto cio' che non e' in standbyMinutesChoices vale il default.
uint32_t standbyMinutes() {
    #ifdef OPENFIRE_STANDBY_MINUTES_FORCE
    return OPENFIRE_STANDBY_MINUTES_FORCE;
    #else
    const uint32_t saved = OF_Prefs::settings[OF_Const::standbyMinutes];
    for(const int choice : OF_Const::standbyMinutesChoices)
        if(saved == uint32_t(choice)) return saved;
    return OF_Const::standbyMinutesDefault; // EN/IT: not in the list / non in lista
    #endif
}

// EN/IT: HID reports prepared but not yet sent / report HID preparati ma non ancora inviati.
bool reportsPending() {
    return TinyUSBDevices.newReport[TinyUSBDevices_::reportMouse] ||
           TinyUSBDevices.newReport[TinyUSBDevices_::reportKeyboard] ||
           TinyUSBDevices.newReport[TinyUSBDevices_::reportGamepad];
}

// EN: Standby only on battery, outside web configuration, in Run mode. Run is checked here too
//     because a dock from the input core can change the mode in the middle of a main loop cycle.
// IT: Standby solo a batteria, fuori dalla configurazione web, in modalita' Run. Il Run si controlla
//     anche qui perche' un dock dal core degli input puo' cambiare modalita' a meta' ciclo principale.
bool eligible() {
    return TinyUSBDevices.onBattery &&
           !OF_WebConfigModeActive &&
           FW_Common::gunMode == FW_Const::GunMode_Run;
}

// EN: Nothing in flight or queued: no transmission in progress (radioFree), empty TX buffer
//     (readIndex/writeIndex), empty ESP-NOW receive queue (_reader/_writer), empty serial write
//     buffer and nothing left to read. The radio task itself is unchanged, as on the pedal.
// IT: Nessun invio in corso o dato in coda: nessuna trasmissione in corso (radioFree), buffer TX vuoto
//     (readIndex/writeIndex), coda di ricezione ESP-NOW vuota (_reader/_writer), buffer di scrittura
//     seriale vuoto e niente da leggere. Il task radio resta invariato, come nel pedale.
bool radioIdle() {
    return radioFree && SerialWireless.readIndex == SerialWireless.writeIndex &&
           SerialWireless._reader == SerialWireless._writer &&
           SerialWireless._readerSerialWrite == SerialWireless._writerSerialWrite && !SerialWireless.available();
}

// EN: Local actuators only: no HID report is sent. Old serial commands must not restart after waking,
//     so the solenoid and rumble queues (Solenoid, SolPulse, Rumble, RumbPulse) are emptied; LED queues stay.
// IT: Solo gli attuatori: nessun report HID. I vecchi comandi seriali non devono ripartire al risveglio,
//     quindi si svuotano le code di solenoide e rumble (Solenoid, SolPulse, Rumble, RumbPulse); le code LED restano.
void stopFeedback() {
    #ifdef MAMEHOOKER
    for(unsigned i = 0; i <= OF_Serial::SerialQueue_RumbPulse; ++i) OF_Serial::serialQueue[i] = false;
    #endif
    OF_FFB::FFBShutdown();
    #ifdef USES_RUMBLE
    // EN: FFBShutdown() leaves the rumble pin on PWM at duty 0, but PWM stops in light sleep:
    //     detach it and hold a steady LOW.
    // IT: FFBShutdown() lascia il pin del rumble sul PWM a duty 0, ma il PWM si ferma nel light sleep:
    //     lo si stacca e si tiene un LOW fisso.
    if(OF_Prefs::pins[OF_Const::rumblePin] >= 0) {
        pinMode(OF_Prefs::pins[OF_Const::rumblePin], OUTPUT);
        digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
    }
    #endif
}

// EN/IT: The SSD1306 keeps its picture while off / l'SSD1306 conserva l'immagine da spento.
void displaySleep(bool sleep) {
    #ifdef USES_DISPLAY
    if(!FW_Common::OLED.display) return;
    #ifdef USE_LOVYAN_GFX
    FW_Common::OLED.display->getPanel()->setSleep(sleep);   // EN/IT: the brightness stays as it is / la luminosita' resta com'e'
    #else
    FW_Common::OLED.display->ssd1306_command(sleep ? SSD1306_DISPLAYOFF : SSD1306_DISPLAYON);
    #endif
    #endif
}

// EN/IT: Trigger read directly from its pin, outside buttons.Poll() / grilletto letto direttamente dal pin, fuori da buttons.Poll().
bool triggerDown() {
    const unsigned i = FW_Const::BtnIdx_Trigger;
    return (digitalRead(LightgunButtons::ButtonDesc[i].pin) ^ ((LightgunButtons::invertMask >> i) & 1U)) == LOW;
}

// EN: Light sleep for ms. Returns true if the trigger woke it, or if the sleep failed: a failure
//     counts as a wake-up, otherwise the blink loop would spin without ever sleeping.
// IT: Light sleep per ms. Restituisce true se l'ha svegliata il grilletto, o se il sonno e' fallito:
//     un errore vale come risveglio, altrimenti il ciclo del lampeggio girerebbe senza mai dormire.
bool sleepFor(uint32_t ms) {
    if(esp_sleep_enable_timer_wakeup(uint64_t(ms) * 1000) != ESP_OK || esp_light_sleep_start() != ESP_OK)
        return true;
    vTaskDelay(1);   // EN/IT: let the idle tasks run between sleeps / lascia girare i task idle fra un sonno e l'altro
    return esp_sleep_get_wakeup_cause() == ESP_SLEEP_WAKEUP_GPIO;
}

// EN: Radio is off: sleep until the trigger, then reconnect and resume the preserved tracking.
//     Light sleep preserves the tracking state, so Run resumes without replacing the last coordinates.
// IT: Radio spenta: dorme fino al grilletto, poi si riconnette e riprende il tracking conservato.
//     Il light sleep conserva lo stato del tracking, quindi il Run riprende senza sostituire le ultime coordinate.
void SleepAndWake() {
    stopFeedback();
    displaySleep(true);
    // EN: Keep the active GPIO settings during light sleep, including actuator gates and LED data lines.
    //     Done here, after stopFeedback() and StandbyBegin() have set their pins.
    // IT: Mantiene la configurazione attiva dei GPIO nel sonno, inclusi gate degli attuatori e linee dati dei LED.
    //     Fatto qui, dopo che stopFeedback() e StandbyBegin() hanno impostato i loro pin.
    for(int pin = 0; pin < GPIO_NUM_MAX; ++pin)
        if(GPIO_IS_VALID_GPIO(pin)) gpio_sleep_sel_dis((gpio_num_t)pin);

    // EN: Wake on the trigger's pressed level (HIGH if the button is inverted). If arming fails
    //     the gun does not sleep and goes straight to the reconnection below.
    // IT: Risveglio sul livello di grilletto premuto (HIGH se il tasto e' invertito). Se l'armamento
    //     fallisce la pistola non dorme e passa subito alla riconnessione qui sotto.
    const gpio_num_t trigger = (gpio_num_t)LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Trigger].pin;
    const bool inverted = (LightgunButtons::invertMask >> FW_Const::BtnIdx_Trigger) & 1U;
    const bool armed = gpio_wakeup_enable(trigger, inverted ? GPIO_INTR_HIGH_LEVEL : GPIO_INTR_LOW_LEVEL) == ESP_OK &&
                       esp_sleep_enable_gpio_wakeup() == ESP_OK;
    // EN/IT: BlinkMs on, then the rest of BlinkPeriodMs off / BlinkMs acceso, poi il resto di BlinkPeriodMs spento.
    while(armed && !triggerDown()) {
        #ifdef LED_ENABLE
        OF_RGB::StandbyBlink(true);
        #endif
        const bool woken = sleepFor(BlinkMs);
        #ifdef LED_ENABLE
        OF_RGB::StandbyBlink(false);
        #endif
        if(woken || sleepFor(BlinkPeriodMs - BlinkMs)) break;
    }
    gpio_wakeup_disable(trigger);
    esp_sleep_disable_wakeup_source(ESP_SLEEP_WAKEUP_GPIO);
    esp_sleep_disable_wakeup_source(ESP_SLEEP_WAKEUP_TIMER);

    // EN: If the dongle does not answer (PC off, dongle restarted) start again as at power-on.
    // IT: Se il dongle non risponde (PC spento, dongle riavviato) riparte come all'accensione.
    if(!SerialWireless.radio_on_gun() || !SerialWireless.connection_gun_at_last_dongle()) ESP.restart();

    // EN: The wake-up pull is not a shot: wait until the trigger has been released for the whole
    //     debounce time. Every reading of "pressed" restarts the count.
    // IT: Il tiro di risveglio non e' uno sparo: attende che il grilletto resti rilasciato per tutto
    //     il tempo di debounce. Ogni lettura di "premuto" fa ripartire il conteggio.
    const uint32_t settle = LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Trigger].debounceTicks;
    for(uint32_t released = millis(); triggerDown() || millis() - released < settle; vTaskDelay(1))
        if(triggerDown()) released = millis();

    displaySleep(false);
}

void Standby() {
    #ifdef DUAL_CORE
    // EN: Stop the input core at the start of its main loop: no inputs processed while the radio is off.
    // IT: Ferma il core degli input all'inizio del suo ciclo principale: nessun input elaborato a radio spenta.
    pauseInputs = 1; __sync_synchronize();
    // EN: Resume runs in its destructor, so the input core restarts on every way out of Standby(),
    //     early returns included. It clears the request before notifying; an early notification
    //     stays pending until the input task waits, and WaitIfPaused() handles it.
    // IT: Resume agisce nel suo distruttore, quindi il core degli input riparte a ogni uscita da Standby(),
    //     compresi i return anticipati. Azzera la richiesta prima di notificare; una notifica anticipata
    //     resta pendente fino all'attesa del task input, e WaitIfPaused() la gestisce.
    struct Resume { ~Resume() { pauseInputs = 0; __sync_synchronize(); xTaskNotifyGive(task_loop1); } } resume;
    // EN: If loop1() is busy (menu, dock, serial processing) it does not reach WaitIfPaused(): give up, the timer retries later.
    // IT: Se loop1() e' occupato (menu, dock, elaborazione seriale) non arriva a WaitIfPaused(): si rinuncia, il timer riprova piu' tardi.
    for(const uint32_t asked = millis(); !inputsPaused; vTaskDelay(1))
        if(millis() - asked >= PauseWaitMs) return;
    __sync_synchronize();   // EN/IT: the reads below come after seeing inputsPaused / le letture qui sotto avvengono dopo aver visto inputsPaused
    #endif
    // EN: A camera update may have queued a report just before this check.
    // IT: Un aggiornamento camera puo' avere preparato un report appena prima di questo controllo.
    if(reportsPending() || !radioIdle()) return;

    #ifdef LED_ENABLE
    OF_RGB::StandbyBegin();
    #endif
    // EN: radioIdle() again: with the input core stopped only radio reception can change it, and the
    //     LED update above takes some time. radio_off_gun() restarts the board if it fails, so radioOff
    //     is false only when the radio was no longer idle.
    // IT: Di nuovo radioIdle(): con il core degli input fermo puo' cambiarlo solo la ricezione radio, e
    //     l'aggiornamento dei LED qui sopra richiede un po' di tempo. radio_off_gun() riavvia la scheda se
    //     fallisce, quindi radioOff e' false solo quando la radio non era piu' libera.
    bool radioOff = radioIdle() && SerialWireless.radio_off_gun();
    // EN: A last packet may arrive during shutdown: resume instead of sleeping.
    // IT: Un ultimo pacchetto puo' arrivare durante lo spegnimento: riprende invece di dormire.
    if(radioOff && !radioIdle()) {
        if(!SerialWireless.radio_on_gun()) ESP.restart();
        radioOff = false;
    }
    if(radioOff) SleepAndWake();
    #ifdef LED_ENABLE
    OF_RGB::StandbyEnd();
    #endif
}
} // namespace

void OF_Standby::Begin() {
    // EN: checkTimer guards against a second call. On USB, or wireless to something other than
    //     the dongle, no timer is created and standby never happens.
    // IT: checkTimer protegge da una seconda chiamata. Su USB, o in wireless verso altro che il
    //     dongle, il timer non viene creato e lo standby non avviene mai.
    if(checkTimer || !TinyUSBDevices.onBattery ||
       TinyUSBDevices.wireless_mode != WIRELESS_MODE::ENABLE_ESP_NOW_TO_DONGLE) return;

    const esp_timer_create_args_t args = {
        // EN: Check inactivity once a minute; the main core handles standby outside this shared timer task.
        //     lastRadioTxMs is copied once because it is volatile and written by the radio code.
        // IT: Controlla l'inattivita' una volta al minuto; il core principale gestisce lo standby fuori dal task condiviso dei timer.
        //     lastRadioTxMs si copia una volta perche' e' volatile e scritto dal codice radio.
        .callback = [](void*) {
            const uint32_t idleMs = standbyMinutes() * 60000UL;
            const uint32_t lastTx = lastRadioTxMs;
            if(idleMs && uint32_t(millis() - lastTx) >= idleMs) standbyRequested = true;
        },
        .arg = nullptr,
        .dispatch_method = ESP_TIMER_TASK,
        .name = "gun_standby",
        .skip_unhandled_events = true
    };
    if(esp_timer_create(&args, &checkTimer) != ESP_OK) return;
    // EN/IT: If the timer does not start, standby stays off until restart / se il timer non parte, lo standby resta spento fino al riavvio.
    if(esp_timer_start_periodic(checkTimer, CheckPeriodUs) != ESP_OK) {
        esp_timer_delete(checkTimer);
        checkTimer = nullptr;
    }
}

void OF_Standby::WaitIfPaused() {
    #ifdef DUAL_CORE
    // EN/IT: In normal play only this read runs / durante il gioco normale si esegue solo questa lettura.
    if(pauseInputs) {
        inputsPaused = 1; __sync_synchronize();
        // EN: Block until resumed, without periodic polling. A stale notification only makes
        //     ulTaskNotifyTake() return once more: the while checks pauseInputs again.
        // IT: Attende la ripresa senza polling periodico. Una notifica vecchia fa solo tornare
        //     ulTaskNotifyTake() una volta in piu': il while ricontrolla pauseInputs.
        while(pauseInputs) ulTaskNotifyTake(pdTRUE, portMAX_DELAY);
        inputsPaused = 0; __sync_synchronize();
    }
    #endif
}

void OF_Standby::Check() {
    // EN: Cleared first: if esp_timer_stop() fails below, Check() does not run again at every frame.
    // IT: Azzerato subito: se esp_timer_stop() qui sotto fallisce, Check() non riparte a ogni frame.
    standbyRequested = false;
    if(!eligible()) {
        lastRadioTxMs = millis();
        return;
    }

    // EN: No timer requests during standby; the one-minute timer resumes only after returning to Run.
    //     It also fails when checkTimer is nullptr after a failed restart: standby then stays off.
    // IT: Nessuna richiesta del timer durante lo standby; il timer di un minuto riparte solo al ritorno nel Run.
    //     Fallisce anche quando checkTimer e' nullptr dopo un riavvio del timer fallito: lo standby resta spento.
    if(esp_timer_stop(checkTimer) != ESP_OK) return;
    Standby();
    // EN: After a wake, or a refused attempt, a whole idle time again. The flag is cleared a second
    //     time to discard a timer callback that was already running when the timer was stopped.
    // IT: Dopo un risveglio o un tentativo rifiutato, di nuovo tutto il tempo. Il flag si azzera una
    //     seconda volta per scartare una callback del timer gia' in esecuzione quando il timer e' stato fermato.
    lastRadioTxMs = millis();
    standbyRequested = false;
    // EN/IT: If the timer does not restart, standby stays off until restart / se il timer non riparte, lo standby resta spento fino al riavvio.
    if(esp_timer_start_periodic(checkTimer, CheckPeriodUs) != ESP_OK) {
        esp_timer_delete(checkTimer);
        checkTimer = nullptr;
    }
}
#endif
