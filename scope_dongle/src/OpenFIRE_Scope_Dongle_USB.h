/*!
 * @file OpenFIRE_Scope_Dongle_USB.h
 * @brief USB link between the PC daemon and the scope dongle
 * @n CPP USB link between the PC daemon and the scope dongle
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#pragma once

// ===================================================================================
// USB PC <-> SCOPE DONGLE / USB PC <-> DONGLE DEL MIRINO
// ===================================================================================
// EN: USB Serial/JTAG of the C5, seen by the PC as a serial port. Frame in both directions:
//       0xA5 0x5A, length (16 bit little-endian, only the content), content
//     Content: 0x4D video packet (forwarded to the scope), 0x43 'C' + ASCII control message,
//     0x54 'T' bandwidth test (counted and discarded). No CRC: the USB already checks the data.
//     PC -> dongle: WHOAMI, START, PAUSE, RESCAN;SKIP=36,40
//     dongle -> PC: OF_SCOPE:BOOT / OF_SCOPE:UNPAIRED;FW=..;CH=.. / OF_SCOPE:Pn;FW=..;CH=..
//                   OF_STAT;RX=..;TX=..;DROP=..;QFULL=.. every 500 ms
// IT: USB Serial/JTAG del C5, vista dal PC come porta seriale. Trama nei due sensi:
//       0xA5 0x5A, lunghezza (16 bit little-endian, solo il contenuto), contenuto
//     Contenuto: 0x4D pacchetto video (inoltrato al mirino), 0x43 'C' + messaggio di controllo ASCII,
//     0x54 'T' prova di banda (contata e scartata). Niente CRC: la USB controlla gia' i dati.
//     PC -> dongle: WHOAMI, START, PAUSE, RESCAN;SKIP=36,40
//     dongle -> PC: OF_SCOPE:BOOT / OF_SCOPE:UNPAIRED;FW=..;CH=.. / OF_SCOPE:Pn;FW=..;CH=..
//                   OF_STAT;RX=..;TX=..;DROP=..;QFULL=.. ogni 500 ms

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32) && defined(SCOPE_DONGLE)

// EN: USB diagnostics on the bottom line of the display (temporary): comment out to remove it.
// IT: Diagnostica della USB sull'ultima riga del display (temporanea): commentare per toglierla.
//#define SCOPE_DONGLE_USB_DEBUG

// EN: Starts the USB and its task: call it first in setup(), so the dongle answers WHOAMI from power on
//     (OF_SCOPE:BOOT during the channel search). The task also tells the PC by itself every change of state.
// IT: Avvia la USB e il suo task: chiamarla per prima nel setup(), cosi' il dongle risponde a WHOAMI dall'accensione
//     (OF_SCOPE:BOOT durante la ricerca del canale). Il task comunica anche da solo al PC ogni cambio di stato.
void scopeUsbBegin();

#ifdef SCOPE_DONGLE_USB_DEBUG
// EN: "W<WHOAMI received> F<answers lost> K<flush again> Q<IN FIFO resets> B<USB bus resets> E<IN buffer free>
//      S<IN endpoint state, hex>" (W and F up to 99, K Q B up to 9)
// IT: "W<WHOAMI ricevuti> F<risposte perse> K<flush ripetuti> Q<azzeramenti FIFO IN> B<reset del bus USB> E<buffer IN libero>
//      S<stato dell'endpoint IN, esadecimale>" (W e F fino a 99, K Q B fino a 9)
void scopeUsbDebugText(char *out, size_t size);
#endif // SCOPE_DONGLE_USB_DEBUG

#endif // OPENFIRE_WIRELESS_ENABLE && ARDUINO_ARCH_ESP32 && SCOPE_DONGLE
