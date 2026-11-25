/*!
 * @brief UART link with the SCOPE / Collegamento UART con lo SCOPE.
 * @copyright Alessandro Satanassi, 2026, GNU Lesser General Public License
 */
#pragma once

#include <stdint.h>

// EN: Only on the ESP32-S3 and only if the pins "Scope UART TX" and "Scope UART RX" are set.
//     115200 baud 8N1, same frame in both directions: 0xA5, type, length, data, CRC-8
//     (polynomial 0x9B, initial value 0, no reflection, no final XOR, on type, length and data).
//     Unknown types are ignored. Types 0x01-0x7F from the gun to the SCOPE, 0x81-0xFF from the SCOPE to the gun.
//     Used by main.cpp: Begin() once the USB or the dongle is connected; by OpenFIREcommon.cpp:
//     SendRoll() at every camera frame. SetDisplay() can be called at any time.
//     Without ARDUINO_ARCH_ESP32 everything below is empty and the compiler removes it.
// IT: Solo sulle ESP32-S3 e solo se i pin "Scope UART TX" e "Scope UART RX" sono impostati.
//     115200 baud 8N1, stessa trama nei due sensi: 0xA5, tipo, lunghezza, dati, CRC-8
//     (polinomio 0x9B, valore iniziale 0, senza riflessione, senza XOR finale, su tipo, lunghezza e dati).
//     I tipi sconosciuti si ignorano. Tipi da 0x01 a 0x7F dalla pistola allo SCOPE, da 0x81 a 0xFF dallo SCOPE alla pistola.
//     Usata da main.cpp: Begin() appena collegati USB o dongle; da OpenFIREcommon.cpp: SendRoll() a ogni
//     frame della camera. SetDisplay() si puo' chiamare in qualsiasi momento.
//     Senza ARDUINO_ARCH_ESP32 tutto quello che segue e' vuoto e il compilatore lo elimina.
class OF_Scope {
public:
#ifdef ARDUINO_ARCH_ESP32
    // EN/IT: Hello up to 3 times, every 200 ms; on Present it sends the player / Ciao fino a 3 volte, ogni 200 ms; al Presente manda il giocatore
    static void Begin();
    // EN/IT: roll in radians (Ang() of Square/Diamond) / roll in radianti (Ang() di Square/Diamond)
    static void SendRoll(float angle);
    // EN/IT: display of the SCOPE on or off / display dello SCOPE acceso o spento
    static void SetDisplay(bool on);
    // EN/IT: inputs of the SCOPE, reserved for the future / ingressi dello SCOPE, riservato per il futuro
    static void OnInputs(const uint8_t *data, uint8_t length);
#else
    static void Begin() {}
    static void SendRoll(float) {}
    static void SetDisplay(bool) {}
#endif
};

// EN/IT: true once the SCOPE answered Present / true appena lo SCOPE ha risposto Presente
#ifdef ARDUINO_ARCH_ESP32
extern volatile bool scope_present;
#else
constexpr bool scope_present = false;
#endif
