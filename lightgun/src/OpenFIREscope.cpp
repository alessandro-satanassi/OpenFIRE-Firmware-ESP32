/*!
 * @brief UART link with the SCOPE / Collegamento UART con lo SCOPE.
 * @copyright Alessandro Satanassi, 2026, GNU Lesser General Public License
 *
 * EN: How it works:
 *  1. Begin(), after the USB or dongle connection: if the two pins are set it opens the UART and sends
 *     Hello up to 3 times, every 200 ms. At Present (scope_present = true) it sends the player (3 times).
 *     No answer: the UART is closed and nothing else is sent.
 *  2. SendRoll(), at every camera frame: the roll in degrees, only when it changes and at most every 15 ms,
 *     plus the last one once more 50 ms after the last change. No answer is awaited.
 *  3. The reception does not poll: onReceive() is called by the UART task of the core when bytes arrive.
 * IT: Come funziona:
 *  1. Begin(), dopo il collegamento USB o dongle: se i due pin sono impostati apre la UART e manda
 *     Ciao fino a 3 volte, ogni 200 ms. Al Presente (scope_present = true) manda il giocatore (3 volte).
 *     Nessuna risposta: la UART si chiude e non si manda altro.
 *  2. SendRoll(), a ogni frame della camera: il roll in gradi, solo quando cambia e al massimo ogni 15 ms,
 *     piu' l'ultimo ancora una volta 50 ms dopo l'ultima variazione. Non si aspetta risposta.
 *  3. La ricezione non interroga: onReceive() la chiama il task UART del core quando arrivano byte.
 */
#include "OpenFIREscope.h"
#ifdef ARDUINO_ARCH_ESP32

#include <Arduino.h>
#include <HardwareSerial.h>
#include <math.h>
#include <string.h>
#include "OpenFIREcommon.h"
#include "OpenFIREprefs.h"

volatile bool scope_present = false;

namespace {
constexpr uint32_t Baud = 115200;
constexpr uint8_t  Start = 0xA5;
constexpr uint8_t  ProtocolVersion = 1;

// EN/IT: gun -> SCOPE / pistola -> SCOPE
constexpr uint8_t CmdHello   = 0x01;   // protocol version / versione protocollo
constexpr uint8_t CmdPlayer  = 0x02;   // player 1-4 / giocatore 1-4
constexpr uint8_t CmdRoll    = 0x03;   // int8 degrees -90..+90, positive = clockwise as seen by the shooter / gradi int8 -90..+90, positivo = orario visto da chi spara
constexpr uint8_t CmdDisplay = 0x04;   // 0 off, 1 on / 0 spento, 1 acceso
// EN/IT: SCOPE -> gun / SCOPE -> pistola
constexpr uint8_t MsgPresent = 0x81;   // protocol version, firmware major, minor, patch / versione protocollo, firmware major, minor, patch
constexpr uint8_t MsgInputs  = 0x83;   // reserved: buttons, wheel / riservato: tasti, rotellina

constexpr uint8_t  MaxData = 16;
constexpr uint8_t  HelloTries = 3;
constexpr uint32_t HelloWaitMs = 200;
constexpr uint32_t RollMinMs = 15;
// EN: the last roll sent once more this long after the last change: a lost frame is not kept while the gun is still
// IT: l'ultimo roll rimandato una volta dopo tanto dall'ultima variazione: una trama persa non resta a pistola ferma
constexpr uint32_t RollRepeatMs = 50;
// EN: Player and Display are sent 3 times: a frame spoiled by a disturbance is not lost (the SCOPE ignores the repeats).
// IT: Giocatore e Display si mandano 3 volte: una trama rovinata da un disturbo non va persa (lo SCOPE ignora le ripetizioni).
constexpr uint8_t  RepeatTimes = 3;
constexpr uint32_t PlayerRepeatMs = 10;
// EN: a new degree is sent only 0.1 degrees past the half: no flicker between two values
// IT: un nuovo grado si manda solo 0,1 gradi oltre la meta': niente sfarfallio fra due valori
constexpr float    RollHysteresis = 0.6f;

// EN/IT: not "Serial...": the gun has substitution macros for Serial / non "Serial...": la pistola ha macro di sostituzione per Serial
HardwareSerial scope_uart(1);

int8_t   lastRoll = 0;       // EN/IT: the SCOPE starts from 0 / lo SCOPE parte da 0
uint32_t lastRollMs = 0;
bool     rollRepeat = false; // EN/IT: last roll still to be sent once more / ultimo roll ancora da rimandare una volta

uint8_t Crc8Step(uint8_t crc, uint8_t value)
{
    crc ^= value;
    for(uint8_t bit = 0; bit < 8; ++bit)
        crc = (crc & 0x80) ? (uint8_t)((crc << 1) ^ 0x9B) : (uint8_t)(crc << 1);
    return crc;
}

void SendFrame(uint8_t type, const uint8_t *data, uint8_t length)
{
    uint8_t frame[4 + MaxData];
    frame[0] = Start;
    frame[1] = type;
    frame[2] = length;
    memcpy(&frame[3], data, length);
    uint8_t crc = 0;
    for(uint8_t i = 1; i < 3 + length; ++i) crc = Crc8Step(crc, frame[i]);
    frame[3 + length] = crc;
    scope_uart.write(frame, 4 + length);
}

void HandleFrame(uint8_t type, const uint8_t *data, uint8_t length)
{
    switch(type) {
    case MsgPresent: scope_present = true; break;
    case MsgInputs:  OF_Scope::OnInputs(data, length); break;
    default: break;   // EN/IT: unknown type: ignored / tipo sconosciuto: ignorato
    }
}

// EN: Called by the UART task of the core, not by the loops of the gun. Wrong CRC: frame discarded.
// IT: La chiama il task UART del core, non i loop della pistola. CRC sbagliato: trama scartata.
void OnReceive()
{
    static uint8_t state = 0, type = 0, length = 0, count = 0, crc = 0;
    static uint8_t data[MaxData];

    while(scope_uart.available() > 0) {
        const uint8_t b = (uint8_t)scope_uart.read();
        switch(state) {
        case 0:   // EN/IT: start / inizio
            if(b == Start) { crc = 0; state = 1; }
            break;
        case 1:   // EN/IT: type / tipo
            type = b; crc = Crc8Step(crc, b); state = 2;
            break;
        case 2:   // EN/IT: length / lunghezza
            if(b > MaxData) { state = 0; break; }
            length = b; count = 0; crc = Crc8Step(crc, b); state = length ? 3 : 4;
            break;
        case 3:   // EN/IT: data / dati
            data[count++] = b; crc = Crc8Step(crc, b);
            if(count == length) state = 4;
            break;
        default:  // EN/IT: CRC
            if(b == crc) HandleFrame(type, data, length);
            state = 0;
            break;
        }
    }
}
} // namespace

void OF_Scope::Begin()
{
    const int8_t txPin = OF_Prefs::pins[OF_Const::scope_UART_TX];
    const int8_t rxPin = OF_Prefs::pins[OF_Const::scope_UART_RX];
    if(txPin < 0 || rxPin < 0)
        return;

    scope_uart.onReceive(OnReceive);
    scope_uart.begin(Baud, SERIAL_8N1, rxPin, txPin);

    for(uint8_t i = 0; i < HelloTries && !scope_present; ++i) {
        const uint8_t version = ProtocolVersion;
        SendFrame(CmdHello, &version, 1);
        const uint32_t start = millis();
        while(!scope_present && millis() - start < HelloWaitMs)
            vTaskDelay(1);
    }

    if(!scope_present) {
        scope_uart.end();
        return;
    }

    // EN/IT: player of the gun, as the Start/Select buttons / giocatore della pistola, come i tasti Start/Select
    const char startBtn = FW_Common::playerStartBtn;
    const uint8_t player = (startBtn >= '1' && startBtn <= '4') ? (uint8_t)(startBtn - '0') : 1;
    for(uint8_t i = 0; i < RepeatTimes; ++i) {
        if(i) vTaskDelay(pdMS_TO_TICKS(PlayerRepeatMs));
        SendFrame(CmdPlayer, &player, 1);
    }
}

void OF_Scope::SendRoll(float angle)
{
    if(!scope_present)
        return;

    float degrees = angle * (180.0f / (float)M_PI);
    if(degrees > 90.0f) degrees = 90.0f;
    else if(degrees < -90.0f) degrees = -90.0f;
    else if(!(degrees == degrees)) return;   // EN/IT: NaN

    const uint32_t now = millis();
    if(fabsf(degrees - lastRoll) < RollHysteresis) {
        // EN/IT: no change: only the single repeat of the last one / nessuna variazione: solo l'unica ripetizione dell'ultimo
        if(rollRepeat && now - lastRollMs >= RollRepeatMs) {
            rollRepeat = false;
            const uint8_t data = (uint8_t)lastRoll;
            SendFrame(CmdRoll, &data, 1);
        }
        return;
    }
    if(now - lastRollMs < RollMinMs)
        return;

    lastRoll = (int8_t)lroundf(degrees);
    lastRollMs = now;
    rollRepeat = true;
    const uint8_t data = (uint8_t)lastRoll;
    SendFrame(CmdRoll, &data, 1);
}

void OF_Scope::SetDisplay(bool on)
{
    if(!scope_present)
        return;
    const uint8_t data = on ? 1 : 0;
    for(uint8_t i = 0; i < RepeatTimes; ++i)
        SendFrame(CmdDisplay, &data, 1);
}

// EN: Reserved: inputs of the SCOPE (buttons, wheel), for future features. Runs in the UART task of the core.
// IT: Riservato: ingressi dello SCOPE (tasti, rotellina), per funzioni future. Gira nel task UART del core.
void OF_Scope::OnInputs(const uint8_t *data, uint8_t length)
{
    (void)data;
    (void)length;
}

#endif // ARDUINO_ARCH_ESP32
