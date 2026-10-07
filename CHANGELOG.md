<a id="english-version"></a>

<p align="center">
  <a href="#english-version"><img src="docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# CHANGELOG

All notable changes to this project will be documented in this file.

## [Unreleased]
### Added
### Changed
### Fixed
### Breaking Changes

---

## [7.0.0]

### Added

* **OpenFIRE ESP32 WebApp**, the new configuration tool: online from a computer browser, through the gun's USB cable or its paired dongle, or offline from the WebApp stored in the gun, over Wi-Fi (phones included) or USB networking. No program to install.
* **Three cameras, one firmware:** DFRobot SEN0158 / Wii Cam, PixArt PAJ7025R2 and the wide-angle PixArt PAJ7025R3, chosen in the WebApp. The firmware corrects the distortion of the R3's wide-angle lens.
* **Startup shortcuts**, held for about 2 seconds while switching on: **Trigger + A** prepares the gun for a firmware update, **B** starts the offline WebApp.
* **Startup mode** (mouse or gamepad) saved in the gun, and a **Wireless pedal** option to skip the pedal search when you do not use one.
* **IR camera test:** each emitter circle shows how large and bright the camera sees its light spot; an emitter that is not seen appears as a dashed circle with a red X. The gun's aim is drawn as a crosshair, and a legend explains every symbol and underlines the layout in use (Square or Diamond).
* **Calibration from the WebApp or the desktop App checks the IR emitters:** the crosshair is green when the camera sees all four emitters well, turns orange when one is weak and red when one is missing. A panel shows each emitter, a target shot is refused while an emitter is missing, and the window explains why. Six dots show the progress, a brief white flash on the crosshair confirms each accepted shot, a reminder asks you to stand centrally without rotating the gun, and a legend explains the colours.
* **Calibration from the gun:** the cursor traces a small circle on each target, so it is easy to spot. With an OLED display, the screen shows where the target is, the step (1/6 to 6/6) and, at the end, how to confirm or repeat.
* **DFRobot/Wii camera:** the IR test and the calibration check also read the brightness of the light spots.
* **Save and Send Settings** pulses while there are unsaved changes (WebApp and desktop App).
* **Analog stick: Invert X Axis and Invert Y Axis** in the **Button Mapping** tab (WebApp and desktop App), for sticks wired or mounted the other way round. They apply to the gamepad stick, the D-pad and the keyboard arrows, and to the stick test; both are off by default, so existing guns behave as before.
* The **Web Flasher** can switch a gun running 7.0.0 into update mode by itself (the serial port changes once: select the new port and retry).
* **Project hub** linking the WebApp, the Web Flasher, the tools and the documentation. Each firmware version opens the WebApp made for it.

### Changed

* Refined IR tracking and faster recovery when some LEDs disappear from view. The Square layout supports a vertical rectangle or a wider one with the LEDs at the screen corners; recalibrate after moving the emitters.
* The pause menu also accepts the Up/Down directional buttons.
* The camera test and calibration screens (WebApp, desktop App and OLED) keep the real proportions of the camera image.
* One firmware file per board, with no separate camera or NoFS/Full files. A clean installation erases the whole flash and writes the same file.
* New configuration protocol and rewritten user documentation. A compatible desktop App for ESP32 7.x remains available on the Tools page.

### Fixed

* Hotkey pause mode: the Left/Right rumble and solenoid toggles now work when no hardware switch is fitted, as in the Simple Pause Menu, and are ignored when a hardware switch is fitted or the rumble/solenoid output is not mapped.
* After the first calibration started from the gun on a new or clean-installed board, the gun returns to normal operation as soon as you confirm the calibration. As before, this first calibration is saved automatically; calibrations started from pause mode must be saved separately.

### Migration from 6.2.1

* A clean installation is recommended; note your settings first, then configure and calibrate again. Erasing removes all saved settings and profiles.
* Use the new WebApp or the compatible ESP32 desktop App; earlier configuration Apps do not work with the new protocol. Game outputs and MAMEHOOKER commands are unchanged.
* Use lightgun, dongle and pedal firmware from the same release. While the gun runs the offline WebApp (B held at startup), its USB serial port is not available.

---

## [6.2.1] - 2026-06-30

### Added
* **Persistent COM Port Assignment (USB Serial Number):** Added a unique USB Serial Number to the lightgun's ESP32 board USB descriptor. This ensures Windows always assigns the same COM port to the same lightgun, regardless of which USB port or USB dongle is used.

---

<a id="versione-italiana"></a>

<p align="center">
  <a href="#english-version"><img src="docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# CRONOLOGIA MODIFICHE

Tutte le modifiche rilevanti apportate a questo progetto saranno documentate in questo file.

## [Non rilasciato]
### Aggiunte
### Modifiche
### Correzioni
### Modifiche Incompatibili

---

## [7.0.0]

### Aggiunte

* **WebApp OpenFIRE ESP32**, il nuovo strumento di configurazione: online dal browser di un computer, tramite il cavo USB della pistola o il suo dongle associato, oppure offline con la WebApp contenuta nella pistola, via Wi-Fi (anche da telefono) o rete USB. Nessun programma da installare.
* **Tre telecamere, un solo firmware:** DFRobot SEN0158 / Wii Cam, PixArt PAJ7025R2 e PixArt PAJ7025R3 grandangolare, da scegliere nella WebApp. Il firmware corregge la distorsione dell'ottica grandangolare della R3.
* **Scorciatoie all'accensione**, da tenere premute per circa 2 secondi: **Grilletto + A** prepara la pistola all'aggiornamento del firmware, **B** avvia la WebApp offline.
* **Modalità all'avvio** (mouse o gamepad) salvata nella pistola, e opzione **Pedale wireless** per saltare la ricerca del pedale quando non lo usi.
* **Test della telecamera IR:** ogni cerchio degli emettitori mostra quanto è grande e luminosa la macchia di luce vista dalla telecamera; un emettitore non visto appare come un cerchio tratteggiato con una X rossa. Il puntamento della pistola è disegnato come un mirino e una legenda spiega ogni simbolo e sottolinea il layout in uso (Square o Diamond).
* **La calibrazione dalla WebApp o dall'App desktop controlla gli emettitori IR:** il mirino è verde quando la telecamera vede bene tutti e quattro gli emettitori, diventa arancione se uno è debole e rosso se ne manca uno. Un riquadro mostra ogni emettitore, il tiro su un bersaglio viene rifiutato finché manca un emettitore e la finestra spiega il motivo. Sei pallini mostrano l'avanzamento, un breve lampo bianco sul mirino conferma ogni tiro accettato, un promemoria ricorda di stare al centro senza ruotare la pistola e una legenda spiega i colori.
* **Calibrazione dalla pistola:** il cursore descrive un piccolo cerchio su ogni bersaglio, così si individua subito. Con un display OLED lo schermo mostra dove si trova il bersaglio, il passo (da 1/6 a 6/6) e, alla fine, come confermare o ripetere.
* **Telecamera DFRobot/Wii:** il test IR e il controllo della calibrazione leggono anche la luminosità delle macchie di luce.
* **Salva e invia impostazioni** pulsa finché ci sono modifiche non salvate (WebApp e App desktop).
* **Stick analogico: Inverti asse X e Inverti asse Y** nella scheda **Mappatura Pulsanti** (WebApp e App desktop), per stick collegati o montati al contrario. Valgono per lo stick del gamepad, il D-pad e le frecce della tastiera, e anche per il test dello stick; sono disattivate di default, quindi le pistole esistenti funzionano come prima.
* Il **Web Flasher** può portare da solo in modalità aggiornamento una pistola con la 7.0.0 (la porta seriale cambia una volta: seleziona la nuova porta e riprova).
* **Portale del progetto** con WebApp, Web Flasher, strumenti e documentazione. Ogni versione del firmware apre la WebApp della stessa versione.

### Modifiche

* Tracciamento IR affinato e recupero più rapido quando alcuni LED escono dal campo visivo. Il layout Square accetta un rettangolo verticale o uno più largo con i LED agli angoli dello schermo; ricalibra dopo aver spostato gli emettitori.
* Il menu di pausa accetta anche i tasti direzionali Su/Giù.
* Le schermate di test della telecamera e di calibrazione (WebApp, App desktop e OLED) mantengono le proporzioni reali dell'immagine della telecamera.
* Un solo file firmware per scheda, senza file separati per telecamera o NoFS/Full. L'installazione pulita cancella l'intera flash e scrive lo stesso file.
* Nuovo protocollo di configurazione e documentazione utente riscritta. Nella pagina Tools resta disponibile un'App desktop compatibile con ESP32 7.x.

### Correzioni

* Modalità pausa Hotkey: i comandi Sinistra/Destra per attivare/disattivare rumble e solenoide ora funzionano quando non è montato un interruttore fisico, come nel Menu di Pausa Semplificato, e vengono ignorati se l'interruttore fisico è montato o se l'uscita rumble/solenoide non è mappata.
* Dopo la prima calibrazione avviata dalla pistola su una scheda nuova o dopo un'installazione pulita, la pistola torna al funzionamento normale appena confermi la calibrazione. Come in precedenza, questa prima calibrazione viene salvata automaticamente; quelle avviate dalla pausa vanno salvate separatamente.

### Passaggio dalla 6.2.1

* È consigliata un'installazione pulita; annota prima le impostazioni, poi configura e calibra di nuovo. La cancellazione elimina tutte le impostazioni e i profili salvati.
* Usa la nuova WebApp o l'App desktop ESP32 compatibile; le App di configurazione precedenti non funzionano con il nuovo protocollo. Uscite di gioco e comandi MAMEHOOKER restano invariati.
* Usa firmware della stessa release per lightgun, dongle e pedale. Mentre la pistola esegue la WebApp offline (B premuto all'avvio), la sua porta seriale USB non è disponibile.

---

## [6.2.1] - 2026-06-30

### Aggiunte
* **Assegnazione fissa della porta COM (USB Serial Number):** Aggiunto un numero di serie USB univoco al descrittore USB della scheda ESP32 della lightgun. Questo garantisce che Windows assegni sempre la stessa porta COM alla stessa lightgun, indipendentemente dalla porta USB o dal dongle USB utilizzato.
