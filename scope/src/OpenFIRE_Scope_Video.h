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

// EN: Output bands of 8 rows for the DMA (2-6): with more bands the rotation ends earlier and the next decoding starts
//     while the DMA works. Beyond 2, a band is taken only if SCOPE_VIDEO_HEAP_MARGIN bytes stay free for the Wi-Fi
//     (its buffers for the packets received are taken while the video runs).
// IT: Fasce di uscita da 8 righe per il DMA (2-6): con piu' fasce la rotazione finisce prima e la decodifica successiva
//     parte mentre il DMA lavora. Oltre 2, una fascia si prende solo se restano liberi SCOPE_VIDEO_HEAP_MARGIN byte per il
//     Wi-Fi (i suoi buffer per i pacchetti ricevuti si prendono mentre il video gira).
#define SCOPE_VIDEO_OUT_SLICES   6
#define SCOPE_VIDEO_HEAP_MARGIN  16384

// EN: Fixed array of the frame in YCbCr (201 coded MCU x 384 bytes), shared with the sprite of the turning logo.
// IT: Array fisso del frame in YCbCr (201 MCU codificati x 384 byte), condiviso con lo sprite del logo che gira.
#define SCOPE_VIDEO_SHARED_BYTES 77184

// EN: Incomplete frame shown 150 ms after its first slice. Normally a frame is shown when it is complete or when the
//     first slice of the next one arrives; this time is only for the last slices lost with no frame after them.
//     It must be longer than the time of a whole frame on the link (the slices arrive spread at the speed of the USB).
// IT: Frame incompleto mostrato 150 ms dopo la sua prima fascia. Di solito un frame si mostra quando e' completo o
//     quando arriva la prima fascia del successivo; questo tempo serve solo per le ultime fasce perse senza un frame
//     dopo. Deve superare il tempo di un frame intero sul collegamento (le fasce arrivano distribuite alla velocita' della USB).
#define SCOPE_VIDEO_SHOW_MS      150
#define SCOPE_VIDEO_BLACK_MS     300   // slice missing from every frame for 300 ms -> black / fascia mancante da tutti i frame per 300 ms -> nera
#define SCOPE_VIDEO_NOSIGNAL_MS  500   // no slice for 500 ms and no pause -> "no signal" / nessuna fascia per 500 ms e nessuna pausa -> "nessun segnale"
#define SCOPE_VIDEO_TASK_PRIO    1     // below the UART (2): the roll always arrives / sotto la UART (2): il roll arriva sempre

// EN: Roll: scope_roll in degrees, from -90 to +90, positive = lightgun turned clockwise as seen by the shooter.
//     The image turns the opposite way to stay level. If on the bench it turns the wrong way, uncomment.
// IT: Roll: scope_roll in gradi, da -90 a +90, positivo = lightgun girata in senso orario vista da chi spara.
//     L'immagine gira al contrario per restare in piano. Se al banco gira dalla parte sbagliata, decommentare.
//#define SCOPE_ROLL_INVERT
#define SCOPE_ROLL_MAX           90
// EN: Steadiness of the image (to tune on the hardware):
//     - ZERO_ZONE: with the image straight, it stays straight while the roll is within +-ZERO_ZONE degrees; after a
//       turn it goes back exactly straight when the roll stays within +-1 degree for 250 ms. 0 = no zone.
//     - HYSTERESIS: elsewhere the image turns only when the roll moves at least this many degrees from the angle shown.
//     - STABLE_MS: a change smaller than 5 degrees must last this long (recoil kicks and tremor peaks are ignored);
//       5 degrees or more is applied at once and then the turn is followed degree by degree. 0 = no filter.
//     ZERO_ZONE 0, HYSTERESIS 1 and STABLE_MS 0 = the image follows every degree.
// IT: Stabilita' dell'immagine (da regolare sull'hardware):
//     - ZERO_ZONE: con l'immagine dritta, resta dritta finche' il roll e' entro +-ZERO_ZONE gradi; dopo una rotazione
//       torna esattamente dritta quando il roll resta entro +-1 grado per 250 ms. 0 = nessuna zona.
//     - HYSTERESIS: altrove l'immagine gira solo quando il roll si sposta di almeno questi gradi dall'angolo mostrato.
//     - STABLE_MS: un cambio sotto i 5 gradi deve durare almeno questo tempo (i colpi del rinculo e i picchi del tremolio
//       sono ignorati); da 5 gradi in su si applica subito e poi la rotazione e' seguita grado per grado. 0 = nessun filtro.
//     ZERO_ZONE 0, HYSTERESIS 1 e STABLE_MS 0 = l'immagine segue ogni grado.
#define SCOPE_ROLL_ZERO_ZONE     3
#define SCOPE_ROLL_HYSTERESIS    3
#define SCOPE_ROLL_STABLE_MS     80

// EN: Statistics on the display, over the logo, every second (memory in KB, frames, slices, times): only while there is
//     no video and the logo is not turning. / IT: statistiche sul display, sopra il logo, ogni secondo (memoria in KB,
//     frame, fasce, tempi): solo mentre non c'e' video e il logo non gira.
#define SCOPE_VIDEO_DEBUG

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

// EN: Memory of the sprite of the turning logo (SCOPE_VIDEO_SHARED_BYTES), the same array of the video frame: NULL
//     while the video has it (try again later). Given back with scopeLogoMemoryGive() after deleteSprite().
//     They set scope_logo_spin_on.
// IT: Memoria dello sprite del logo che gira (SCOPE_VIDEO_SHARED_BYTES), lo stesso array del frame del video: NULL
//     finche' ce l'ha il video (riprovare dopo). Si restituisce con scopeLogoMemoryGive() dopo deleteSprite().
//     Impostano scope_logo_spin_on.
uint8_t *scopeLogoMemoryTake();
void scopeLogoMemoryGive();

// EN: Paired, no video on screen, no slice for 500 ms and no valid pause: the main shows "no signal".
// IT: Abbinato, nessun video sullo schermo, nessuna fascia da 500 ms e nessuna pausa valida: il main mostra "nessun segnale".
bool scopeVideoNoSignal();

#endif // OPENFIRE_WIRELESS_ENABLE && ARDUINO_ARCH_ESP32 && SCOPE && USES_DISPLAY && USE_LOVYAN_GFX
