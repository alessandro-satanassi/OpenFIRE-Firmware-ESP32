/*!
 * @file OpenFIRE_Scope_Video.h
 * @brief Video of the scope: newest frame, roll rotation and display
 * @n CPP Video of the scope: newest frame, roll rotation and display
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#pragma once

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32) && defined(SCOPE) && defined(USES_DISPLAY) && defined(USE_LOVYAN_GFX)

#define LGFX_USE_V1
#include <LovyanGFX.hpp>
#include "LGFX_128_GC9A01_240x240.hpp"

// ===================================================================================
// VIDEO OF THE SCOPE / VIDEO DEL MIRINO
// ===================================================================================
// EN: Own task: takes the slices of the newest frame from the ring of the library, decodes them, rotates the
//     frame with the roll of the lightgun (scope_roll) and sends it to the display with the DMA.
//     The panel is shared with the state screens of the main: whoever draws takes the lock first.
//     While the video is on screen the lock belongs to the video; at "no signal" it goes back to the main.
// IT: Task proprio: prende dall'anello della libreria le fasce del frame piu' recente, le decodifica, ruota il
//     frame con il roll della lightgun (scope_roll) e lo manda al display con il DMA.
//     Il pannello e' condiviso con le schermate di stato del main: chi disegna prende prima il blocco.
//     Mentre il video e' sullo schermo il blocco e' del video; al "nessun segnale" torna al main.

// EN: Output slices for the DMA (2-6): with more slices the rotation ends earlier and the next decoding starts
//     while the DMA works. Without memory it goes down to 2.
// IT: Fasce di uscita per il DMA (2-6): con piu' fasce la rotazione finisce prima e la decodifica successiva parte
//     mentre il DMA lavora. Senza memoria scende fino a 2.
#define SCOPE_VIDEO_OUT_SLICES   6

#define SCOPE_VIDEO_SHOW_MS      40    // a frame is shown 40 ms after its first slice / un frame si mostra 40 ms dopo la sua prima fascia
#define SCOPE_VIDEO_BLACK_MS     100   // slice missing from every frame for 100 ms -> black / fascia mancante da tutti i frame per 100 ms -> nera
#define SCOPE_VIDEO_NOSIGNAL_MS  500   // no slice for 500 ms and no pause -> "no signal" / nessuna fascia per 500 ms e nessuna pausa -> "nessun segnale"
#define SCOPE_VIDEO_TASK_PRIO    1     // below the UART (2): the roll always arrives / sotto la UART (2): il roll arriva sempre

// EN: Roll: scope_roll in degrees, from -90 to +90, positive = lightgun turned clockwise as seen by the shooter.
//     The image turns the opposite way to stay level. If on the bench it turns the wrong way, uncomment.
// IT: Roll: scope_roll in gradi, da -90 a +90, positivo = lightgun girata in senso orario vista da chi spara.
//     L'immagine gira al contrario per restare in piano. Se al banco gira dalla parte sbagliata, decommentare.
//#define SCOPE_ROLL_INVERT
#define SCOPE_ROLL_MAX           90

// EN: Statistics on the USB serial every second (frames, slices, times). / IT: statistiche sulla seriale USB ogni secondo (frame, fasce, tempi).
//#define SCOPE_VIDEO_DEBUG

// EN: To call once, after the display init and before the task of the state screens.
// IT: Da chiamare una volta, dopo l'avvio del display e prima del task delle schermate di stato.
void scopeVideoBegin(LGFX *display);

// EN: Lock of the panel for the state screens: false if the video is on screen (nothing to draw).
// IT: Blocco del pannello per le schermate di stato: false se il video e' sullo schermo (niente da disegnare).
bool scopePanelTryLock();
void scopePanelUnlock();

// EN: Display on or off as scope_display_on: only the panel, reception and decoding go on. Call with the panel locked.
// IT: Display acceso o spento come scope_display_on: solo il pannello, ricezione e decodifica continuano. Da chiamare col pannello bloccato.
void scopePanelPowerApply();

// EN: Paired, no video on screen, no slice for 500 ms and no valid pause: the main shows "no signal".
// IT: Abbinato, nessun video sullo schermo, nessuna fascia da 500 ms e nessuna pausa valida: il main mostra "nessun segnale".
bool scopeVideoNoSignal();

#endif // OPENFIRE_WIRELESS_ENABLE && ARDUINO_ARCH_ESP32 && SCOPE && USES_DISPLAY && USE_LOVYAN_GFX
