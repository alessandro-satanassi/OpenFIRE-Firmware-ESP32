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

* Official configuration WebApp: online via USB serial (directly or through the paired dongle), and offline from the lightgun over Wi-Fi or USB NCM.
* PAJ7025R2 and wide-angle PAJ7025R3 cameras, selectable alongside DFRobot/Wii in the same board firmware.
* Startup shortcuts: **Trigger + A** for firmware flashing; **B** for the integrated WebApp. Hold for about 2 seconds.
* Saved mouse/gamepad startup mode and an explicit wireless pedal option.
* IR camera test: each emitter circle shows the size and brightness of the light spot measured by the camera; an emitter that is not seen is shown as a dashed circle with a red X.
* Calibration from the WebApp or the desktop App checks the IR emitters: the crosshair turns green when all four are seen well and orange when one is weak, a small panel shows each emitter, a pulsing halo highlights the crosshair, and target shots are refused while an emitter is not seen. The calibration started from the gun works as before.
* DFRobot/Wii camera: the IR test and the calibration check also read the brightness of the light spots, from the camera's full data format.
* **Save and Send Settings** pulses while there are unsaved changes (WebApp and desktop App).
* The Web Flasher can switch a lightgun running 7.0.0 into update mode by itself (the serial port changes once: select the new port and retry).
* Central project hub, version-matched WebApp and Web Flasher.

### Changed

* Refined IR tracking and recovery when some LEDs disappear. Square layout supports vertical or wide rectangles, including screen-corner placement; recalibrate after moving the emitters.
* The pause menu also accepts the Up/Down directional buttons.
* The camera test and calibration views (WebApp, desktop App and OLED) keep the real proportions of the camera image.
* One firmware image per board, not separate camera or NoFS/Full images. Clean installation erases the whole flash before writing the same image.
* Updated configuration communication and user documentation. The compatible ESP32 desktop App remains available from the Tools page.

### Fixed

* Hotkey pause mode: the Left/Right rumble and solenoid toggles now work when no hardware switch is configured, as in the Simple Pause Menu, and are ignored when a hardware switch is fitted or the rumble/solenoid output is not mapped.

### Migration from 6.2.1

* A clean installation is recommended; note your settings first, then configure and calibrate again. Erasing removes all saved settings and profiles.
* Use the new WebApp or a compatible ESP32 desktop App; old configuration Apps are incompatible with the new protocol. Normal game HID outputs and MAMEHOOKER commands remain available.
* Use matching releases for lightgun, dongle and pedal. USB serial is unavailable on the lightgun while its special USB network configuration mode is active.

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

* WebApp ufficiale di configurazione: online tramite seriale USB (direttamente o attraverso il dongle associato), e offline dalla lightgun tramite Wi-Fi o USB NCM.
* Telecamere PAJ7025R2 e PAJ7025R3 grandangolare, selezionabili insieme a DFRobot/Wii nello stesso firmware della scheda.
* Scorciatoie all'avvio: **Grilletto + A** per il flashing; **B** per la WebApp integrata. Tenere premuto per circa 2 secondi.
* Modalità mouse/gamepad salvabile per l'avvio e opzione esplicita per il pedale wireless.
* Test della telecamera IR: ogni cerchio degli emettitori mostra grandezza e luminosità della macchia di luce misurata dalla telecamera; un emettitore non visto è indicato da un cerchio tratteggiato con una X rossa.
* La calibrazione dalla WebApp o dall'App desktop controlla gli emettitori IR: il mirino diventa verde quando tutti e quattro sono visti bene e arancione quando uno è debole, un piccolo riquadro mostra ogni emettitore, un alone pulsante mette in evidenza il mirino e i tiri sui bersagli vengono rifiutati finché un emettitore non è visto. La calibrazione avviata dalla pistola funziona come prima.
* Telecamera DFRobot/Wii: il test IR e il controllo della calibrazione leggono anche la luminosità delle macchie di luce, dal formato dati completo della telecamera.
* **Salva e invia impostazioni** pulsa finché ci sono modifiche non salvate (WebApp e App desktop).
* Il Web Flasher può portare da solo in modalità aggiornamento una lightgun che esegue la 7.0.0 (la porta seriale cambia una volta: selezionare la nuova porta e riprovare).
* Hub centrale del progetto, WebApp abbinata alla versione e Web Flasher.

### Modifiche

* Affinati il tracciamento IR e il recupero quando alcuni LED scompaiono. Il layout Square supporta rettangoli verticali o larghi, anche con LED agli angoli dello schermo; ricalibrare dopo aver spostato gli emettitori.
* Il menu di pausa accetta anche i tasti direzionali Su/Giù.
* Le schermate di test della telecamera e di calibrazione (WebApp, App desktop e OLED) mantengono le proporzioni reali dell'immagine della telecamera.
* Un'unica immagine firmware per scheda, senza varianti separate per telecamera o NoFS/Full. L'installazione pulita cancella l'intera flash prima di scrivere la stessa immagine.
* Aggiornate la comunicazione di configurazione e la documentazione utente. L'App desktop ESP32 compatibile resta disponibile nella pagina Tools.

### Correzioni

* Modalità pausa Hotkey: i comandi Sinistra/Destra per attivare/disattivare rumble e solenoide ora funzionano quando non è configurato un interruttore fisico, come nel Menu di Pausa Semplificato, e vengono ignorati se è presente l'interruttore fisico o se l'uscita rumble/solenoide non è mappata.

### Passaggio dalla 6.2.1

* È consigliata un'installazione pulita; annotare prima le impostazioni, poi configurare e calibrare nuovamente. La cancellazione elimina tutte le impostazioni e i profili salvati.
* Usare la nuova WebApp o un'App desktop ESP32 compatibile; le vecchie App di configurazione non sono compatibili con il nuovo protocollo. Restano disponibili le normali uscite HID di gioco e i comandi MAMEHOOKER.
* Usare release corrispondenti per lightgun, dongle e pedale. La seriale USB della lightgun non è disponibile mentre è attiva la modalità speciale di configurazione con rete USB.

---

## [6.2.1] - 2026-06-30

### Aggiunte
* **Assegnazione fissa della porta COM (USB Serial Number):** Aggiunto un numero di serie USB univoco al descrittore USB della scheda ESP32 della lightgun. Questo garantisce che Windows assegni sempre la stessa porta COM alla stessa lightgun, indipendentemente dalla porta USB o dal dongle USB utilizzato.

