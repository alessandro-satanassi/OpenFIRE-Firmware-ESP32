<a id="english-version"></a>

[🏠 Back to Home](../../README.md#english-version) / [Lightgun Firmware](../README.md#english-version) / **Operational Manual**

<p align="center">
  <a href="#english-version"><img src="../../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# OpenFIRE - The Enclosed Instruction Book!

*... from the original [OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/blob/OpenFIRE-dev/OpenFIREmain/README.md) project repository:*

## Table of Contents:
 - [IR Emitter Setup](#ir-emitter-setup)
 - [Board Configuration](#board-configuration)
 - [Operations Manual](#operations-manual)
   - [Run Modes](#run-modes)
   - [Default Buttons](#default-buttons)
   - [Default Buttons in Pause Mode](#default-buttons-in-pause-mode-hotkey)
   - [How to Calibrate](#how-to-calibrate)
   - [IR Camera Sensitivity](#ir-camera-sensitivity)
   - [Profiles](#profiles)
   - [Software Toggles](#software-toggles)
   - [Saving Settings to Flash](#saving-settings-to-flash)
   - [Test Mode](#test-mode)
 - [Technical Details & Assorted Errata](#technical-details--assorted-errata)
   - [Serial Handoff (Mame Hooker) Mode](#serial-handoff-mame-hooker-mode)
   - [Change USB ID for Multiple Guns](#change-usb-id-for-multiple-guns)

## IR Emitter setup
The IR emitters can be arranged in either two ways:
 - As dual lightbars, or "square" layout *(vertical rectangle with the width shorter than the height)* (Samco, preferred)
 - At four separate points on the display, or "diamond" layout (Xwiigun and other legacy lightgun systems)

if you're playing on a small PC monitor, you can use 2 Wii sensor bars; one on top of your screen and one below. However, if you're playing on a TV, you should consider building or buying a set of high power black IR LEDs and arranging them like (larger) sensor bars at the top and bottom of the display.

The **OpenFIRE Desktop App** has an alignment assistant that can be used to help align your emitters to the display (by selecting ***Help->Open IR Emitter Alignment Assistant***). 

<table>
  <tr>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_app_001.png">
        <img src="../docs/img/IR_Emitter_app_001.png">
      </a>
    </td>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_app_002.png">
        <img src="../docs/img/IR_Emitter_app_002.png">
      </a>
    </td>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_monitor.jpeg">
        <img src="../docs/img/IR_Emitter_monitor.jpeg">
      </a>
    </td>
  </tr>
</table>

## Board Configuration
The gun is configured through the companion [OpenFIRE Desktop Application](https://github.com/TeamOpenFIRE/OpenFIRE-App), which can be launched even without a microcontroller plugged in and will automatically search for compatible flashed devices. You can easily download the latest version of the App from the **[OpenFIRE ESP32 Tools](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-Tools/?lang=en)** portal.

When the board's COM port is selected in the app, the gun will go into a *Docked* state - this is what allows for real-time configuration of pin mappings, settings, changing between and renaming calibration profiles, and testing button inputs and force feedback devices. The camera will be disabled while in this mode, unless the gun is set to IR test mode.

## First-time Setup
When flashing a new board with OpenFIRE, or after clearing the flash, the first time it's plugged in will prompt for the user to pull the trigger button to start initial calibration - this can be accomplished from the App using any of the *Calibrate Profile* buttons, or pressing the trigger for standalone calibration (see the [How to Calibrate](#how-to-calibrate) section for more information). If your build is using custom pins, or you would like to change any settings at this point, the gun can be docked to the OpenFIRE App and configured prior to starting initial calibration - at least *Trigger* and *Button A* should be mapped and confirmed working in the *Gun Tests* tab.

## Operations Manual
The light gun operates as an absolute positioning mouse (like a stylus!) until the button/combination is pressed to enter pause mode. Alternatively, the gun can be signaled to output using its corresponding HID Gamepad device using a Serial Feedback Distributor program such as MAMEHOOKER - see the [Serial Handoff](#serial-handoff-mame-hooker-mode) section for more info.

Any serial terminal (Arduino IDE's Serial Monitor, *PuTTY,* *screen,* etc.) can be used to see information while the gun is paused and during standalone calibration.

Note that the buttons in pause mode (and to enter pause mode) activate when the last button of the combination releases. This is used to detect and differentiate button combinations vs a single button press.

* Note: At its peak, the mouse position updates at 209Hz, or roughly every ~4.8ms, so it is extremely responsive.

### Run modes
The gun has the following modes of operation:
1. Normal - The mouse position updates from each frame from the IR positioning camera (no averaging)
2. Averaging - The position is calculated from a 2 frame moving average (current + previous position)
3. Averaging2 - The position is calculated from a weighted average of the current frame and 2 previous frames
4. Processing - Test mode for use with the Desktop App (this mode is prevented from being assigned to a profile)

The averaging modes are subtle but do reduce the motion jitter a bit without adding much if any noticeable lag.

> The ESP32 port of OpenFIRE automatically incorporates advanced anti-jitter algorithms, so it is recommended to use **Normal** mode.

### Default Buttons
- Trigger: Left mouse button
- A: Right mouse button (In low buttons mode, Start if pressed offscreen)
- B: Middle mouse button (In low buttons mode, Select if pressed offscreen)
- C/Reload: Mouse button 4/Side Button 1/Back
- Pump Action (Cabela's or alike): Right mouse button
- Start: 1 key
- Select: 5 key
- Up/Down/Left/Right: Keyboard arrow keys
- Pedal Main: Mouse button 4/Side Button 1/Back
- Alt Pedal: Mouse button 5/Side Button 2/Forward
- C + Start: Esc key

Pause mode can be entered by either pressing C + Select by default, pressing the *Home Button* if used in current pin layout, or *holding the trigger plus the A Button with **no IR points in sight*** if hold-to-pause is enabled - pointing the gun towards the ground is recommended here.

#### Default Buttons in Pause mode (Hotkey)
- A, B, Start, Select: select a profile
- Start + Down: Normal gun mode (averaging disabled)
- Start + Up: Normal gun with averaging, switch between the 2 averaging modes (use serial monitor to see the setting)
- B + Down: Decrease IR camera sensitivity (use a serial monitor to see the setting)
- B + Up: Increase IR camera sensitivity (use a serial monitor to see the setting)
- C/Reload: Exit pause mode
- Left: Toggle Rumble *(when no rumble switch is detected)*
- Right: Toggle Solenoid *(when no solenoid switch is detected)*
- Trigger: Begin calibration
- Start + Select: save settings to non-volatile flash storage space

#### Controls for Simple Pause Menu
- A: Navigate Cursor Up
- B: Navigate Cursor Down
- Trigger: Select option
- C: Exit pause mode
  - Holding A or B for half the duration of the hold-to-pause time (so ~2s by default) will also exit the simple pause menu.
Available options in simple pause menu are as follows, from first option to last before rolling back:
* Calibrate current profile (always the initial option)
* Switch profiles (submenu)
  * Select from profile 1-4 using the navigation buttons/trigger to select, or press C to back out.
* Save settings to non-volatile memory
* Toggle Rumble *(when rumble is enabled & no switch is detected)*
* Toggle Solenoid *(when solenoid is enabled & no switch is detected)*
* Send escape key signal to the PC

### How to calibrate
##### These instructions apply to the standalone on-board calibration process; the Calibration screens in the OpenFIRE App has a similar procedure with more info to guide the user through the process.
1. Select the profile to calibrate - either through pressing A/B/Start/Select in the Hotkey Pause Mode, or selecting in the Simple Pause Menu - and pull the trigger to begin calibration. Alternatively, calibration can be initialized from the OpenFIRE App. 
2. Shoot the pointer at center of the screen and press the trigger while keeping a steady aim.
3. The mouse should move to the four edges of the screen; first topmost, bottommost, leftmost, rightmost. Shoot the edge of the screen where the cursor is at.
4. When the cursor returns to the center, shoot the cursor to finish calibration.
5. The new calibration profile will be applied and you'll be able to test the tracking. A good sign of a good calibration is maintaining as close to line-of-sight accuracy as possible when aiming at the screen edges and corners.
   - If the calibration is good, pull the trigger to confirm.
   - If you want to start calibration over, press the A or B button in cali verification to restart from the center point in Step 2.
   - Calibration can be canceled outright by pressing C/Reload at any time, or the A/B buttons any time before cali verification.
Remember to save your calibration and current profile afterwards, either by sending the save signal from the app, pressing Start+Select in Hotkey Pause Mode, or selecting the third "Save Settings" option in the Simple Pause Menu.

### IR Camera Sensitivity
The IR camera sensitivity can be adjusted. It is recommended to adjust the sensitivity as high as possible. If the IR sensitivity is too low then the pointer precision can suffer. However, too high of a sensitivity can cause the camera to pick up unwanted reflections that will cause the pointer to jump around. It is impossible to know which setting will work best since it is dependent on the specific setup. It depends on how bright the IR emitters are, the distance, camera lens, and if shiny surfaces may cause reflections.

A sign that the IR sensitivity is too low is if the pointer moves in noticeable coarse steps, as if it has a low resolution to it. If you have the sensitivity level set to max and you notice this, the IR emitters may not be bright enough.

A sign that the IR sensitivity is too high is if the pointer jumps around erratically. If this happens only while aiming at certain areas of the screen, this is a good indication that a reflection is being detected by the camera. If the sensitivity is at max, step it down to high or minimum. Obviously, the best solution is to eliminate the reflective surface. The Desktop App's Test Mode can help diagnose this problem, since it will visually display the 4 IR points.

### Profiles
The main OpenFIRE builds are configured with 4 calibration profiles available. Each profile has its own calibration data, run mode, and IR camera sensitivity settings. Each profile can be selected from pause mode by pressing the associated button (A/B/Start/Select), or selecting them via the profiles submenu in simple pause menu.

### Software Toggles
Hardware features can be toggled at runtime, even without hardware switches defined!

While in pause mode, the toggles are as follows (color indicating what the board's builtin LED lights up with):
- Left D-Pad: **Rumble Toggle** (Salmon) - Enables/disables the rumble functionality. When enabled, the motor will engage for a short period.
- Right D-Pad: **Solenoid Toggle** (Yellow) - Enables/disables the solenoid force feedback. When enabled, the solenoid will engage for a short period.
These can also be done from the respective setting in the Simple Pause Menu.

The current state of these settings (except Offscreen Button Mode) are saved when committed to, and pulled from flash storage space at boot.

#### Saving Settings to Flash
The calibration data, profile settings, and extended gun options like custom pins mapping and rumble intensity, can be saved in non-volatile memory by pressing Start + Select from the gun itself, or whenever new settings are committed from the Desktop App. The currently selected calibration profile is saved as the default for when the light gun is plugged in - gun settings (pins mapping, force feedback, etc.) applies to *all profiles.*

#### Test Mode
Test Mode lets you visually see the IR points as seen by the camera in the Desktop App's IR Test Mode screen. This is very useful for aligning the camera when building your light gun, and for testing that the camera tracks all 4 points properly, as well as observing possible reflections. The validity of the test points shape (square in default IR layout, diamond in alt IR layout) depends on the current profile used and its IR layout setting.

## Technical Details & Assorted Errata

### Serial Handoff (Mame Hooker) Mode
The gun will automatically hand off control to an instance of Mame Hooker that's connected once a start code has been detected! If available, the onboard LED and any *non-static* external NeoPixels will change to a mid-intensity white to signal serial handoff mode (unless any LED events trigger it to change, which will follow those thereafter).

If you aren't already familiar with Mame Hooker, **you'll need compatible inis for each game you play** and **the gun's COM port should be set to match the player number** (COM1 for P1, COM2 for P2, etc.)! COM port assignment can be done in Windows via the Device Manager, or Linux via settings in the Wine registry of the prefix your game/Mame Hooker is started in. [Consult the wiki page on MAMEHOOKER for more information!](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/wiki/MAMEHOOKER_Documentation_EN) For Linux users wanting to use their gun with native emulators' force feedback (currently MAME, Flycast, or their RetroArch ports), consider trying [QMamehook](https://github.com/SeongGino/QMamehook).

### Change USB ID for Multiple Guns
If you intend to use multiple OpenFIRE guns, you'll want to change what the board reports itself as.

These are known as the **USB Implementer's Forum (USB-IF) identifiers**, and if multiple devices share a common display name and/or Product/Vendor ID, apps like RetroArch and TeknoParrot that read individual mouse devices will get VERY confused.

These parameters can be saved to and loaded from flash storage space if a customized Device Product ID (PID) & Device Name are detected, which can be committed from the Desktop App.

---
### 💬 Questions or Issues?
For technical support and to join the discussion, please refer to the [Community & Support Section](../README.md#community-support-english) in the Main Repository.


---

<p align="center"> 🔸 🔸 🔸 </p>

---

<a id="versione-italiana"></a>

[🏠 Torna alla Home](../../README.md#versione-italiana) / [Lightgun Firmware](../README.md#versione-italiana) / **Manuale Operativo**

<p align="center">
  <a href="#english-version"><img src="../../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# OpenFIRE - Il Manuale di utilizzo!

*... dal sito del progetto originale [OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/blob/OpenFIRE-dev/OpenFIREmain/README.md):*

## Indice:
 - [Configurazione Emettitori IR](#configurazione-emettitori-ir-italiano)
 - [Configurazione della Scheda](#configurazione-della-scheda-italiano)
 - [Prima Configurazione](#prima-configurazione-italiano)
 - [Manuale Operativo](#manuale-operativo-italiano)
   - [Modalità di Funzionamento](#modalità-di-funzionamento-italiano)
   - [Pulsanti Predefiniti](#pulsanti-predefiniti-italiano)
   - [Pulsanti Predefiniti in Modalità Pausa](#pulsanti-predefiniti-in-modalità-pausa-hotkey-italiano)
   - [Controlli per il Menu di Pausa Semplificato](#controlli-per-il-menu-di-pausa-semplificato-italiano)
   - [Come Calibrare](#come-calibrare-italiano)
   - [Sensibilità della Telecamera IR](#sensibilità-della-telecamera-ir-italiano)
   - [Profili](#profili-italiano)
   - [Interruttori Software (Toggle)](#interruttori-software-toggle-italiano)
   - [Salvataggio delle Impostazioni nella Flash](#salvataggio-delle-impostazioni-nella-flash-italiano)
   - [Modalità di Test](#modalità-di-test-italiano)
 - [Dettagli Tecnici e Note Varie](#dettagli-tecnici-e-note-varie-italiano)
   - [Modalità Serial Handoff (Mame Hooker)](#modalità-serial-handoff-mame-hooker-italiano)
   - [Modifica dell'ID USB per Pistole Multiple](#modifica-dell-id-usb-per-pistole-multiple-italiano)


<a id="configurazione-emettitori-ir-italiano"></a>

## Configurazione Emettitori IR
Gli emettitori IR possono essere disposti in due modi:
 - Come doppie barre luminose, o layout "rettangolare" *(rettangolo verticale con base più corta dell'altezza)* (Samco, consigliato)
 - In quattro punti separati sul display, o layout a "diamante" (Xwiigun e altri vecchi sistemi lightgun)

Se giochi su un piccolo monitor per PC, puoi usare 2 barre sensore Wii; una sopra lo schermo e una sotto. Tuttavia, se giochi su una TV, dovresti prendere in considerazione la costruzione o l'acquisto di un set di LED IR neri ad alta potenza e disporli come barre sensore (più grandi) nella parte superiore e inferiore del display.

L'App Desktop OpenFIRE dispone di un assistente di allineamento che può aiutarti ad allineare gli emettitori al display (selezionando *Help -> Open IR Emitter Alignment Assistant*).

<table>
  <tr>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_app_001.png">
        <img src="../docs/img/IR_Emitter_app_001.png">
      </a>
    </td>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_app_002.png">
        <img src="../docs/img/IR_Emitter_app_002.png">
      </a>
    </td>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_monitor.jpeg">
        <img src="../docs/img/IR_Emitter_monitor.jpeg">
      </a>
    </td>
  </tr>
</table>

<a id="configurazione-della-scheda-italiano"></a>

## Configurazione della Scheda
La pistola viene configurata tramite l'[Applicazione Desktop OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-App), che può essere avviata anche senza aver prima collegato un microcontrollore e cercherà automaticamente i dispositivi già flashati compatibili. Puoi scaricare in modo facilitato la versione più recente dell'App dal portale **[OpenFIRE ESP32 Tools](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-Tools/?lang=it)**.

Quando la porta COM della scheda viene selezionata nell'app, la pistola entrerà in uno stato *Docked* (Ancorato) - questo è ciò che consente la configurazione in tempo reale della mappatura dei pin, delle impostazioni, la modifica e rinomina dei profili di calibrazione e il test degli input dei pulsanti e dei dispositivi di force feedback. La telecamera sarà disabilitata in questa modalità, a meno che la pistola non sia impostata sulla modalità di test IR.

<a id="prima-configurazione-italiano"></a>

## Prima Configurazione
Quando si esegue il flashing di una nuova scheda con OpenFIRE, o dopo aver formattato la memoria flash, la prima volta che la si collega verrà richiesto all'utente di premere il grilletto per avviare la calibrazione iniziale. Questo può essere fatto dall'App utilizzando uno qualsiasi dei pulsanti *Calibrate Profile*, oppure premendo il grilletto fisicamente per la calibrazione autonoma (vedi la sezione [Come Calibrare](#come-calibrare-italiano) per maggiori informazioni). Se la tua build utilizza pin personalizzati, o desideri modificare qualsiasi impostazione a questo punto, la pistola può essere collegata all'App OpenFIRE e configurata prima di avviare la calibrazione iniziale - assicurati che almeno il *Trigger* (Grilletto) e il *Button A* siano mappati e confermati come funzionanti nella scheda *Gun Tests*.

<a id="manuale-operativo-italiano"></a>

## Manuale Operativo
La lightgun funziona come un mouse a posizionamento assoluto (come un pennino per tablet!) finché non viene premuto il pulsante/combinazione per entrare in modalità pausa. In alternativa, si può istruire la pistola a inviare output tramite il suo corrispondente dispositivo HID Gamepad utilizzando un programma di distribuzione di feedback seriale come MAMEHOOKER - vedi la sezione [Modalità Serial Handoff](#modalità-serial-handoff-mame-hooker-italiano) per maggiori informazioni.

Qualsiasi terminale seriale (Monitor Seriale dell'IDE Arduino, *PuTTY*, *screen*, ecc.) può essere utilizzato per visualizzare le informazioni mentre la pistola è in pausa e durante la calibrazione autonoma.

Nota che i pulsanti in modalità pausa (e la combinazione per entrare in modalità pausa) si attivano quando viene rilasciato *l'ultimo pulsante* della combinazione. Questo sistema viene utilizzato per rilevare e differenziare le combinazioni di tasti dalla singola pressione.

* Nota: Al suo picco, la posizione del mouse si aggiorna a 209Hz (sia via cavo USB che via wireless), o all'incirca ogni ~4.8ms, risultando estremamente reattiva.

<a id="modalità-di-funzionamento-italiano"></a>

### Modalità di Funzionamento
La pistola ha le seguenti modalità operative:
1. **Normal** - La posizione del mouse si aggiorna a ogni frame ricevuto dalla telecamera di posizionamento IR (nessuna mediazione).
2. **Averaging** - La posizione è calcolata tramite una media mobile su 2 frame (posizione attuale + precedente).
3. **Averaging2** - La posizione è calcolata tramite una media ponderata del frame attuale e dei 2 frame precedenti.
4. **Processing** - Modalità di test da utilizzare con l'App Desktop (non è possibile assegnare questa modalità a un profilo).

Le modalità *Averaging* sono sottili ma riducono un po' il jitter (tremolio) del movimento senza aggiungere lag percettibile.

> Il porting di OpenFIRE per ESP32 aggiunge in automatico algoritmi avanzati anti jitter (tremolio), quindi si consiglia di impostare la modalità **Normal**.

<a id="pulsanti-predefiniti-italiano"></a>

### Pulsanti Predefiniti
- **Trigger (Grilletto):** Tasto sinistro del mouse
- **A:** Tasto destro del mouse (In modalità "low buttons", agisce da Start se premuto puntando fuori dallo schermo)
- **B:** Tasto centrale del mouse (In modalità "low buttons", agisce da Select se premuto puntando fuori dallo schermo)
- **C/Reload:** Tasto mouse 4 / Pulsante laterale 1 / Indietro
- **Pump Action (Ricarica a pompa, es. Cabela's):** Tasto destro del mouse
- **Start:** Tasto 1 della tastiera
- **Select:** Tasto 5 della tastiera
- **Su/Giù/Sinistra/Destra:** Frecce direzionali della tastiera
- **Pedale Principale:** Tasto mouse 4 / Pulsante laterale 1 / Indietro
- **Pedale Secondario (Alt):** Tasto mouse 5 / Pulsante laterale 2 / Avanti
- **C + Start:** Tasto Esc della tastiera

Si può entrare in modalità Pausa premendo **C + Select** (impostazione predefinita), premendo il tasto **Home** (se presente nel layout dei pin), oppure **tenendo premuto il grilletto insieme al pulsante A senza alcun punto IR in vista** se l'opzione *hold-to-pause* è abilitata - in quest'ultimo caso si consiglia di puntare la pistola verso il pavimento.

<a id="pulsanti-predefiniti-in-modalità-pausa-hotkey-italiano"></a>

### Pulsanti Predefiniti in Modalità Pausa (Hotkey)
- **A, B, Start, Select:** Seleziona un profilo.
- **Start + Giù:** Modalità pistola Normal (Averaging disabilitato).
- **Start + Su:** Modalità pistola Normal con Averaging, passa da una modalità di media all'altra (usa il monitor seriale per vedere l'impostazione).
- **B + Giù:** Diminuisci la sensibilità della telecamera IR (usa il monitor seriale per vedere l'impostazione).
- **B + Su:** Aumenta la sensibilità della telecamera IR (usa il monitor seriale per vedere l'impostazione).
- **C/Reload:** Esci dalla modalità pausa.
- **Sinistra:** Attiva/Disattiva Rumble *(se non viene rilevato alcuno switch fisico per il rumble)*.
- **Destra:** Attiva/Disattiva Solenoide *(se non viene rilevato alcuno switch fisico per il solenoide)*.
- **Trigger:** Inizia la calibrazione.
- **Start + Select:** Salva le impostazioni nello spazio di archiviazione flash non volatile.

<a id="controlli-per-il-menu-di-pausa-semplificato-italiano"></a>

#### Controlli per il Menu di Pausa Semplificato
- **A:** Muovi il cursore Su
- **B:** Muovi il cursore Giù
- **Trigger:** Seleziona l'opzione
- **C:** Esci dal menu di pausa
  - *Tenendo premuto A o B per metà della durata del tempo di hold-to-pause (circa ~2 secondi di default) si uscirà anche dal menu di pausa semplice.*
  
Le opzioni disponibili nel menu di pausa semplificato sono le seguenti (dalla prima all'ultima, per poi ricominciare):
* Calibra il profilo corrente (sempre la prima opzione iniziale)
* Cambia profilo (sottomenu)
  * Scegli tra i profili 1-4 usando i pulsanti di navigazione/grilletto per selezionare, o premi C per tornare indietro.
* Salva le impostazioni nella memoria non volatile
* Attiva/Disattiva Rumble *(quando abilitato e senza switch fisico)*
* Attiva/Disattiva Solenoide *(quando abilitato e senza switch fisico)*
* Invia segnale del tasto Esc al PC

<a id="come-calibrare-italiano"></a>

### Come Calibrare
*Queste istruzioni si applicano al processo di calibrazione integrato e autonomo della pistola; le schermate di Calibrazione nell'App OpenFIRE hanno una procedura simile con maggiori informazioni a schermo per guidare l'utente.*

1. Seleziona il profilo da calibrare (tramite A/B/Start/Select nella modalità Hotkey Pause, o selezionandolo nel Menu Pausa Semplificato) e premi il grilletto per iniziare la calibrazione. In alternativa, puoi avviarla dall'App OpenFIRE.
2. Punta il mirino al centro dello schermo e premi il grilletto mantenendo la mira stabile.
3. Il cursore si sposterà sui quattro bordi dello schermo; prima in alto, poi in basso, a sinistra e a destra. Spara sul bordo dello schermo nel punto esatto in cui si trova il cursore.
4. Quando il cursore torna al centro, spara sul cursore per terminare la calibrazione.
5. Il nuovo profilo di calibrazione verrà applicato e potrai testare il tracciamento. Un buon indicatore di una calibrazione corretta è il mantenimento di una precisione il più vicino possibile alla linea di vista (line-of-sight) quando si mira ai bordi e agli angoli dello schermo.
   - Se la calibrazione è buona, premi il grilletto per confermare.
   - Se desideri ricominciare la calibrazione, premi il pulsante A o B nella schermata di verifica per ripartire dal punto centrale (Passo 2).
   - La calibrazione può essere annullata del tutto premendo C/Reload in qualsiasi momento, o i pulsanti A/B in qualsiasi momento prima della verifica finale.

Ricordati di **salvare la calibrazione** e il profilo corrente subito dopo, inviando il segnale di salvataggio dall'app, premendo *Start+Select* nella modalità Hotkey Pause, o scegliendo la terza opzione "Save Settings" nel Menu di Pausa Semplificato.

<a id="sensibilità-della-telecamera-ir-italiano"></a>

### Sensibilità della Telecamera IR
La sensibilità della telecamera IR può essere regolata. Si consiglia di impostarla il più in alto possibile. Se la sensibilità IR è troppo bassa, la precisione del puntatore ne risentirà. Tuttavia, una sensibilità troppo elevata potrebbe far sì che la telecamera rilevi riflessi indesiderati, causando salti improvvisi del puntatore. È impossibile sapere a priori quale impostazione funzionerà meglio, poiché dipende dalle specifiche del tuo setup (luminosità degli emettitori IR, distanza, lente della telecamera ed eventuali superfici lucide che causano riflessi).

Un segno che la sensibilità IR è **troppo bassa** si verifica quando il puntatore si muove in modo "scattoso" e poco fluido, come se avesse una bassa risoluzione. Se noti questo problema nonostante la sensibilità sia impostata al massimo, è probabile che i tuoi emettitori IR non siano abbastanza luminosi.

Un segno che la sensibilità IR è **troppo alta** si verifica quando il puntatore salta in modo irregolare ed erratico. Se ciò accade solo mentre miri a determinate aree dello schermo, è un chiaro indicatore che la telecamera sta rilevando un riflesso. Se la sensibilità è al massimo, riducila su alto o minimo. Ovviamente, la soluzione migliore rimane l'eliminazione della superficie riflettente. La Modalità Test dell'App Desktop può aiutare a diagnosticare questo problema visualizzando a schermo i 4 punti IR.

<a id="profili-italiano"></a>

### Profili
Le build principali di OpenFIRE sono configurate con 4 profili di calibrazione disponibili. Ogni profilo ha i propri dati di calibrazione, la modalità operativa (run mode) e le impostazioni di sensibilità della telecamera IR. Ogni profilo può essere richiamato dalla modalità pausa premendo il pulsante associato (A/B/Start/Select) o selezionandolo dal sottomenu dei profili nel menu di pausa semplificato.

<a id="interruttori-software-toggle-italiano"></a>

### Interruttori Software (Toggle)
Le funzioni hardware possono essere attivate o disattivate in tempo reale, anche se non hai cablato degli interruttori fisici dedicati!

Mentre sei in modalità pausa, i controlli toggle sono i seguenti (il colore indica come si illumina il LED integrato sulla scheda):
- **D-Pad Sinistra: Rumble Toggle** (Salmone) - Abilita/disabilita la funzione rumble (vibrazione). Quando abilitato, il motore si attiverà per un breve periodo come conferma.
- **D-Pad Destra: Solenoid Toggle** (Giallo) - Abilita/disabilita il force feedback del solenoide. Quando abilitato, il solenoide scatterà per un breve periodo come conferma.

Queste operazioni possono essere eseguite anche dalle rispettive opzioni nel Menu di Pausa Semplificato. Lo stato corrente di queste impostazioni (ad eccezione della modalità Offscreen Button) viene salvato nella memoria flash al momento del salvataggio manuale e ricaricato all'avvio.

<a id="salvataggio-delle-impostazioni-nella-flash-italiano"></a>

#### Salvataggio delle Impostazioni nella Flash
I dati di calibrazione, le impostazioni dei profili e le opzioni estese della pistola (come la mappatura personalizzata dei pin e l'intensità del rumble) possono essere salvati nella memoria non volatile premendo **Start + Select** sulla pistola stessa, oppure ogni volta che vengono applicate nuove impostazioni dall'App Desktop. Il profilo di calibrazione attualmente selezionato al momento del salvataggio viene impostato come predefinito per le successive accensioni. Le impostazioni generali della pistola (mappatura pin, force feedback, ecc.) si applicano a *tutti i profili*.

<a id="modalità-di-test-italiano"></a>

#### Modalità di Test
La Modalità di Test ti consente di visualizzare i punti IR esattamente come vengono visti dalla telecamera attraverso la schermata "IR Test Mode" dell'App Desktop. Questo è estremamente utile per allineare la telecamera durante la costruzione della tua lightgun e per verificare che tracci correttamente tutti e 4 i punti, oltre a farti osservare eventuali riflessi fastidiosi. La validità della forma dei punti di test (un quadrato nel layout IR predefinito, un diamante nel layout IR alternativo) dipende dal profilo attualmente in uso e dalle sue impostazioni relative al layout IR.

<a id="dettagli-tecnici-e-note-varie-italiano"></a>

## Dettagli Tecnici e Note Varie

<a id="modalità-serial-handoff-mame-hooker-italiano"></a>

### Modalità Serial Handoff (Mame Hooker)
La pistola passerà automaticamente il controllo a un'istanza in esecuzione di Mame Hooker (o app similari collegate in seriale) non appena rileverà un codice di avvio compatibile! Se disponibile, il LED integrato sulla scheda e gli eventuali NeoPixel esterni *non statici* diventeranno di colore bianco a media intensità per segnalare l'attivazione della modalità Serial Handoff (salvo modifiche in-game che impongano colori diversi).

Se non hai familiarità con Mame Hooker, **avrai bisogno dei file `.ini` compatibili per ogni gioco che utilizzi** e **la porta COM della pistola dovrebbe essere impostata per corrispondere al numero del giocatore** (COM1 per il P1, COM2 per il P2, ecc.)! L'assegnazione della porta COM può essere fatta su Windows tramite "Gestione dispositivi", o su Linux tramite le impostazioni di registro di Wine nel prefisso in cui avvii il gioco/Mame Hooker. [Consulta la pagina wiki su MAMEHOOKER per maggiori informazioni!](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/wiki/MAMEHOOKER_Documentation_IT) Per gli utenti Linux che desiderano utilizzare la pistola con il force feedback nativo degli emulatori (attualmente MAME, Flycast o i relativi port su RetroArch), si consiglia di provare [QMamehook](https://github.com/SeongGino/QMamehook).

<a id="modifica-dell-id-usb-per-pistole-multiple-italiano"></a>

### Modifica dell'ID USB per Pistole Multiple
Se intendi utilizzare più lightgun OpenFIRE contemporaneamente sullo stesso PC, dovrai modificare i dati identificativi della scheda.

Questi sono noti come **identificatori USB Implementer's Forum (USB-IF)**. Se più dispositivi condividono un nome display e/o un ID Prodotto/Venditore (PID/VID) comune, applicazioni come RetroArch e TeknoParrot (che leggono individualmente i singoli mouse) entreranno fortemente in confusione.

Questi parametri possono essere modificati tramite l'App Desktop e salvati/caricati dallo spazio di archiviazione flash non appena vengono rilevati un Device Product ID (PID) e un Nome Dispositivo personalizzati.

---
### 💬 Domande o Problemi?
Per supporto tecnico e per unirti alla community, consulta la [Sezione Community e Supporto](../../README.md#community-support-italiano) nella Home del progetto.