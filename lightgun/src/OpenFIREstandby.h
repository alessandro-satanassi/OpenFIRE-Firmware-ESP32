/*!
 * @brief Battery ESP-NOW standby / Standby ESP-NOW a batteria.
 * @copyright Alessandro Satanassi, 2026, GNU Lesser General Public License
 */
#pragma once

// EN: Minutes without ESP-NOW transmission attempts before standby; 0 disables it.
// IT: Minuti senza tentativi di trasmissione ESP-NOW prima dello standby; 0 lo disabilita.
#ifndef OPENFIRE_STANDBY_MINUTES
    #define OPENFIRE_STANDBY_MINUTES 15
#endif

class OF_Standby {
public:
#if defined(ARDUINO_ARCH_ESP32) && defined(OPENFIRE_WIRELESS_ENABLE) && OPENFIRE_STANDBY_MINUTES > 0
    static volatile bool standbyRequested;
    static void Begin(); // EN/IT: one-minute timer, started once in setup / timer di un minuto, avviato una volta nel setup
    static void WaitIfPaused(); // EN/IT: input loop waits here during standby / il ciclo degli input attende qui durante lo standby
    static void Check(); // EN/IT: main core, only when requested by the timer / core principale, solo su richiesta del timer
#else
    static constexpr bool standbyRequested = false;
    static void Begin() {}
    static void WaitIfPaused() {}
    static void Check() {}
#endif
};
