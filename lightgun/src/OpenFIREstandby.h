/*!
 * @brief Battery ESP-NOW standby / Standby ESP-NOW a batteria.
 * @copyright Alessandro Satanassi, 2026, GNU Lesser General Public License
 */
#pragma once

// EN: Minutes without ESP-NOW transmission attempts before standby: saved setting standbyMinutes
//     (list and default in OpenFIREshared.h), 0 = standby disabled.
//     For tests, -D OPENFIRE_STANDBY_MINUTES_FORCE=N in platformio.ini uses N whatever is saved (0 = disabled).
// IT: Minuti senza tentativi di trasmissione ESP-NOW prima dello standby: impostazione salvata standbyMinutes
//     (lista e default in OpenFIREshared.h), 0 = standby disabilitato.
//     Per le prove, -D OPENFIRE_STANDBY_MINUTES_FORCE=N in platformio.ini usa N qualunque cosa sia salvata (0 = disabilitato).
#ifdef OPENFIRE_STANDBY_MINUTES
    #warning "OPENFIRE_STANDBY_MINUTES is ignored: use OPENFIRE_STANDBY_MINUTES_FORCE"
#endif

// EN: Used by main.cpp: Begin() at the end of setup(), Check() in ExecRunMode() after GetPosition()
//     when standbyRequested is set, WaitIfPaused() at the start of every loop1() cycle.
//     Without ESP32 + OPENFIRE_WIRELESS_ENABLE everything below is empty and the compiler removes it.
// IT: Usata da main.cpp: Begin() alla fine di setup(), Check() in ExecRunMode() dopo GetPosition()
//     quando standbyRequested e' alzato, WaitIfPaused() all'inizio di ogni ciclo di loop1().
//     Senza ESP32 + OPENFIRE_WIRELESS_ENABLE tutto quello che segue e' vuoto e il compilatore lo elimina.
class OF_Standby {
public:
#if defined(ARDUINO_ARCH_ESP32) && defined(OPENFIRE_WIRELESS_ENABLE)
    // EN/IT: raised by the timer, cleared by Check() and on entering Run / alzato dal timer, azzerato da Check() e all'ingresso nel Run
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
