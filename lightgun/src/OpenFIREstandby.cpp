/*!
 * @brief Battery ESP-NOW standby / Standby ESP-NOW a batteria.
 * @copyright Alessandro Satanassi, 2026, GNU Lesser General Public License
 *
 * EN: After OPENFIRE_STANDBY_MINUTES without ESP-NOW transmission attempts
 * the gun turns the radio off and light-sleeps, blinking faintly. The trigger wakes it:
 * same dongle, no new pairing, and the wake-up pull is not sent as a shot.
 * IT: Dopo OPENFIRE_STANDBY_MINUTES senza tentativi di trasmissione ESP-NOW
 * la pistola spegne la radio e va in light sleep, con un lampeggio tenue. Il grilletto la
 * risveglia: stesso dongle, nessun nuovo abbinamento, e il tiro di risveglio non e' uno sparo.
 */
#include "OpenFIREstandby.h"
#if defined(ARDUINO_ARCH_ESP32) && defined(OPENFIRE_WIRELESS_ENABLE) && OPENFIRE_STANDBY_MINUTES > 0

#include "OpenFIREcommon.h"
#include "OpenFIREFeedback.h"
#include "OpenFIREserial.h"
#include "OpenFIRElights.h"
#include "OpenFIREweb.h"
#include <esp_sleep.h>
#include <esp_timer.h>
#include <driver/gpio.h>

extern volatile bool radioFree;

volatile bool OF_Standby::standbyRequested = false;

namespace {
constexpr uint32_t IdleMs = OPENFIRE_STANDBY_MINUTES * 60000UL;
constexpr uint64_t CheckPeriodUs = 60ULL * 1000000ULL;
esp_timer_handle_t checkTimer = nullptr;
// EN/IT: Standby sign, a short blink every few seconds / segno di standby, breve lampeggio ogni pochi secondi.
constexpr uint32_t BlinkMs = 40, BlinkPeriodMs = 3000;
// EN/IT: Longest wait for the input core to stop / attesa massima perche' il core dei tasti si fermi.
constexpr uint32_t PauseWaitMs = 100;

#ifdef DUAL_CORE
volatile uint32_t pauseInputs = 0, inputsPaused = 0;
#endif

bool reportsPending() {
    return TinyUSBDevices.newReport[TinyUSBDevices_::reportMouse] ||
           TinyUSBDevices.newReport[TinyUSBDevices_::reportKeyboard] ||
           TinyUSBDevices.newReport[TinyUSBDevices_::reportGamepad];
}

bool eligible() {
    return TinyUSBDevices.onBattery &&
           !OF_WebConfigModeActive &&
           FW_Common::gunMode == FW_Const::GunMode_Run;
}

// EN: Nothing in flight or queued. The radio task itself is unchanged, as on the pedal.
// IT: Nessun invio in corso o dato in coda. Il task radio resta invariato, come nel pedale.
bool radioIdle() {
    return radioFree && SerialWireless.readIndex == SerialWireless.writeIndex &&
           SerialWireless._reader == SerialWireless._writer &&
           SerialWireless._readerSerialWrite == SerialWireless._writerSerialWrite && !SerialWireless.available();
}

// EN: Local actuators only: no HID report is sent. Old serial commands must not restart after waking.
// IT: Solo gli attuatori: nessun report HID. I vecchi comandi seriali non devono ripartire al risveglio.
void stopFeedback() {
    #ifdef MAMEHOOKER
    for(unsigned i = 0; i <= OF_Serial::SerialQueue_RumbPulse; ++i) OF_Serial::serialQueue[i] = false;
    #endif
    OF_FFB::FFBShutdown();
    #ifdef USES_RUMBLE
    // EN/IT: Steady LOW instead of PWM, which stops in light sleep / LOW fisso invece del PWM, che si ferma nel sonno.
    if(OF_Prefs::pins[OF_Const::rumblePin] >= 0) {
        pinMode(OF_Prefs::pins[OF_Const::rumblePin], OUTPUT);
        digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
    }
    #endif
}

// EN: The SSD1306 keeps its picture while off / l'SSD1306 conserva l'immagine da spento.
void displaySleep(bool sleep) {
    #ifdef USES_DISPLAY
    if(!FW_Common::OLED.display) return;
    #ifdef USE_LOVYAN_GFX
    FW_Common::OLED.display->getPanel()->setSleep(sleep);   // the brightness stays as it is
    #else
    FW_Common::OLED.display->ssd1306_command(sleep ? SSD1306_DISPLAYOFF : SSD1306_DISPLAYON);
    #endif
    #endif
}

bool triggerDown() {
    const unsigned i = FW_Const::BtnIdx_Trigger;
    return (digitalRead(LightgunButtons::ButtonDesc[i].pin) ^ ((LightgunButtons::invertMask >> i) & 1U)) == LOW;
}

// EN/IT: Light sleep for ms; true if the trigger woke it (or the sleep failed) / true se l'ha svegliata il grilletto.
bool sleepFor(uint32_t ms) {
    if(esp_sleep_enable_timer_wakeup(uint64_t(ms) * 1000) != ESP_OK || esp_light_sleep_start() != ESP_OK)
        return true;
    vTaskDelay(1);   // EN/IT: let the idle tasks run between sleeps / lascia girare i task idle fra un sonno e l'altro
    return esp_sleep_get_wakeup_cause() == ESP_SLEEP_WAKEUP_GPIO;
}

// EN: Radio is off: sleep until the trigger, then reconnect and resume the preserved tracking.
// IT: Radio spenta: dorme fino al grilletto, poi si riconnette e riprende il tracking conservato.
void SleepAndWake() {
    stopFeedback();
    displaySleep(true);
    // EN: Keep the active GPIO settings during light sleep, including actuator gates and LED data lines.
    // IT: Mantiene la configurazione attiva dei GPIO nel sonno, inclusi gate degli attuatori e linee dati dei LED.
    for(int pin = 0; pin < GPIO_NUM_MAX; ++pin)
        if(GPIO_IS_VALID_GPIO(pin)) gpio_sleep_sel_dis((gpio_num_t)pin);

    const gpio_num_t trigger = (gpio_num_t)LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Trigger].pin;
    const bool inverted = (LightgunButtons::invertMask >> FW_Const::BtnIdx_Trigger) & 1U;
    const bool armed = gpio_wakeup_enable(trigger, inverted ? GPIO_INTR_HIGH_LEVEL : GPIO_INTR_LOW_LEVEL) == ESP_OK &&
                       esp_sleep_enable_gpio_wakeup() == ESP_OK;
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

    // EN: The wake-up pull is not a shot: wait for the release, stable for the trigger's debounce time.
    // IT: Il tiro di risveglio non e' uno sparo: attende il rilascio, stabile per il debounce del grilletto.
    const uint32_t settle = LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Trigger].debounceTicks;
    for(uint32_t released = millis(); triggerDown() || millis() - released < settle; vTaskDelay(1))
        if(triggerDown()) released = millis();

    // EN: Light sleep preserves tracking state; resume without replacing the last coordinates.
    // IT: Il light sleep conserva lo stato del tracking; riprende senza sostituire le ultime coordinate.

    displaySleep(false);
}

void Standby() {
    #ifdef DUAL_CORE
    // EN: Stop the input core at the start of its main loop: no inputs processed while the radio is off.
    // IT: Ferma il core degli input all'inizio del suo ciclo principale: nessun input elaborato a radio spenta.
    pauseInputs = 1; __sync_synchronize();
    struct Resume { ~Resume() { pauseInputs = 0; __sync_synchronize(); } } resume;
    for(const uint32_t asked = millis(); !inputsPaused; vTaskDelay(1))
        if(millis() - asked >= PauseWaitMs) return;   // EN/IT: busy (menu, dock): later / occupato: piu' tardi
    __sync_synchronize();
    #endif
    // EN: A camera update may have queued a report just before this check.
    // IT: Un aggiornamento camera puo' avere preparato un report appena prima di questo controllo.
    if(reportsPending() || !radioIdle()) return;

    #ifdef LED_ENABLE
    OF_RGB::StandbyBegin();
    #endif
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
    if(checkTimer || !TinyUSBDevices.onBattery ||
       TinyUSBDevices.wireless_mode != WIRELESS_MODE::ENABLE_ESP_NOW_TO_DONGLE) return;

    const esp_timer_create_args_t args = {
        // EN: Check inactivity once a minute; the main core handles standby outside this shared timer task.
        // IT: Controlla l'inattivita' una volta al minuto; il core principale gestisce lo standby fuori dal task condiviso dei timer.
        .callback = [](void*) {
            const uint32_t lastTx = lastRadioTxMs;
            if(uint32_t(millis() - lastTx) >= IdleMs) standbyRequested = true;
        },
        .arg = nullptr,
        .dispatch_method = ESP_TIMER_TASK,
        .name = "gun_standby",
        .skip_unhandled_events = true
    };
    if(esp_timer_create(&args, &checkTimer) != ESP_OK) return;
    if(esp_timer_start_periodic(checkTimer, CheckPeriodUs) != ESP_OK) {
        esp_timer_delete(checkTimer);
        checkTimer = nullptr;
    }
}

void OF_Standby::WaitIfPaused() {
    #ifdef DUAL_CORE
    if(pauseInputs) {
        inputsPaused = 1; __sync_synchronize();
        while(pauseInputs) vTaskDelay(1);
        inputsPaused = 0; __sync_synchronize();
    }
    #endif
}

void OF_Standby::Check() {
    standbyRequested = false;
    if(!eligible()) {
        lastRadioTxMs = millis();
        return;
    }

    // EN: No timer requests during standby; the one-minute timer resumes only after returning to Run.
    // IT: Nessuna richiesta del timer durante lo standby; il timer di un minuto riparte solo al ritorno nel Run.
    if(esp_timer_stop(checkTimer) != ESP_OK) return;
    Standby();
    // EN/IT: After a wake, or a refused attempt, a whole idle time again / dopo un risveglio o un tentativo rifiutato, di nuovo tutto il tempo.
    lastRadioTxMs = millis();
    standbyRequested = false;
    if(esp_timer_start_periodic(checkTimer, CheckPeriodUs) != ESP_OK) {
        esp_timer_delete(checkTimer);
        checkTimer = nullptr;
    }
}
#endif
