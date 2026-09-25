<a id="english-version"></a>

<p align="center">
  <a href="#english-version"><img src="../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# Lightgun Firmware (ESP32-S3) - NEWS

## OpenFIRE ESP32 - ver 7.0.0 Release Candidate 1 (RC1)

The OpenFIRE core remains practically unchanged and guarantees total compatibility (100%), however, the surrounding ecosystem has undergone a profound technical revision. This 7.0.0 RC1 replaces the current stable version 6.2.1. 

The new features introduced are detailed below.

### Tracking, Sensors, and Calibration
*   **Improved lost LED management:** The tracking system has been refined to handle the loss of IR LED signals more efficiently, a frequent condition when extremely close to the monitor or during fast movements. The camera rotation management up to approximately ±80 degrees remains fully supported.
*   **Free sensor positioning (Square Configuration):** In addition to the standard recommended formation (vertical rectangle with the base smaller than the height), the IR transmitters can now be placed directly at the corners of the monitor. Recognition is automatic, although the ideal "Square" formation remains the preferable configuration for maximum performance (naturally, calibration must be performed again whenever the IR transmitters are repositioned).
*   **Dynamic Multi-Cam Support:** The firmware now supports three different camera modules, dynamically selectable from the Configuration App:
    1.  DFRobot SEN0158 / Wii Cam (I2C connection - 2 wires)
    2.  PixArt PAJ7025R2 (SPI connection - 4 wires)
    3.  PixArt PAJ7025R3 with integrated wide-angle lens (SPI connection - 4 wires)

> **WARNING - IR LED Selection:** 
> For DFRobot / Wii Cam cameras, it is necessary to use 940nm IR LEDs. 
> For PixArt PAJ7025R2 and R3 sensors, the use of 850nm IR LEDs is strictly required. Although PixArt sensors can detect 940nm LEDs if they have high power, they are components specifically designed for the 850nm wavelength. Using different frequencies will significantly degrade the tracking quality and the maximum operational distance.

*   **Realistic visualization:** The camera test and calibration interface (both on the microcontroller display and on the App) has been slightly modified. The image is now more proportioned and faithful to the real geometry.

### Transition to the Integrated WebApp
This release introduces a structural change in the configuration interface. The Qt-based desktop configuration App has been updated (implementing the new protocol, multi-cam support, and the wireless pedal), but for the OpenFIRE ESP32 project, it is now to be considered "Legacy". 

The new official configuration tool is a **WebApp** built from scratch. The interface replicates the Qt App 1:1 in layout and functionality, eliminating the need to download or install any software on the PC.
*   **Online Configuration:** By accessing the dedicated web portal, the system automatically identifies the firmware version installed on the lightgun and launches the correct instance of the WebApp (https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/).
*   **Offline Configuration (Integrated into the Microcontroller):** The entire WebApp package is physically stored within the microcontroller. By booting the device in "Configuration Mode", the microcontroller exposes a Wi-Fi Access Point (IP: 192.168.4.1) or, if connected via USB OTG cable, creates a virtual NCM netUSB network (accessible at IP 192.168.7.1 or by typing "openfire" in the browser).
*   **Operational flexibility:** The configuration WebApp is fully operational even when the lightgun is interfaced to the PC wirelessly via ESP-NOW DONGLE. The usage combinations are almost completely supported.
*   *Development note:* Porting from the Qt environment to the WebApp required a significant investment of resources. The final refinements on the interface texts are underway, but at a functional level, the system is fully operational. The code of the upstream repository (Board) will continue to be monitored to ensure maximum compatibility and for the porting of any future implementations.

### Communication and Flashing
*   **New Serial Protocol:** Communication between the configuration application and the lightgun is managed via a rewritten serial protocol, structured to ensure greater robustness, fault tolerance, and scalability for future updates.
*   **Libraries Update:** The base firmware libraries have been updated to the latest versions.
*   **Flashing via "Magic Key":** The firmware is now able to interpret the flashing command via software, automatically disconnecting the USB OTG. It is no longer necessary to physically access the microcontroller to press the physical Boot/Reset buttons. The procedure will require two attempts (the time necessary for the serial port to re-enumerate post-reboot), but will succeed on the second boot. It is still possible to force the flashing mode by holding down the appropriate button combination for 2 seconds when booting the device.

### General Improvements and Ecosystem
*   **Menu Navigation:** Following a user's suggestion, it is now possible to scroll through the menus in pause mode using the "Up" and "Down" directional buttons as well, in addition to the classic "A" and "B" buttons.
*   Various minor bugfixes and source code optimizations have been implemented.
*   **Central Hub:** A unified web portal has been published that gathers the entire OpenFIRE ESP32 ecosystem, available at: https://alessandro-satanassi.github.io/OpenFIRE-ESP32/

Given the extent and depth of the modifications made to the code, we kindly ask you to thoroughly test this RC1 and report any anomalies or suggestions before the release of the stable version ("Stable"). 

**A clean installation (Erase Flash) is highly recommended**. Although the configuration files of previous versions are theoretically compatible and tend to align after the first saves, starting from a clean memory prevents potential conflicts.

## SPECIAL BOOT MODES

By holding down specific button combinations for at least 2 seconds during boot, the lightgun will enter one of the following special modes:

**'TRIGGER' Button + 'A' Button (Firmware Update Mode)**
The microcontroller prepares for firmware flashing. The message *"Ready for firmware update"* will appear on the OLED display. At this point, you can proceed with the update using the [Web Flasher](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/) or the `esptool` utility (for advanced users).

**'B' Button (Offline Configuration Mode via internal WebApp)**
This mode starts a webserver integrated into the firmware, allowing you to configure the lightgun via a locally hosted WebApp without needing an internet connection. Successful activation is confirmed on the OLED display: the top bar will have inverted colors and show a gear icon in the top right corner.

System behavior varies depending on how the lightgun is connected:

**1. Connection via DONGLE (Wireless)**
*   If booted in this mode without a USB cable, the lightgun will scan the network channels waiting to pair with an available Dongle. This is essential to test "in-game" tracking even in *undocked* mode.
*   **WebApp Access:** Done via the Wi-Fi Access Point exposed directly by the microcontroller, named `"OpenFIRE_Config"` (IP: `192.168.4.1`). By connecting to this network (even from a smartphone), the WebApp should open automatically in your browser.
*   **Functionality:** The system is fully operational. Since the virtual serial port remains active, the use of Mamehooker is also supported.
*   *Technical note:* With the serial port active, you could also use the online WebApp via Web Serial API at [OpenFIRE-ESP32-WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/). However, the "Button B" mode is specifically designed to provide an alternative when no internet connection is available.

**2. Connection via USB OTG cable**
*   **WebApp Access:** In addition to the standard Wi-Fi Access Point (`192.168.4.1`), the system creates a virtual NCM netUSB network. This allows you to access the interface by simply typing `openfire` or the IP `192.168.7.1` in the browser's address bar, without needing to connect to the microcontroller's AP.
*   **Functionality:** Due to technical limitations related to the maximum number of available USB endpoints, the serial port is disabled. The lightgun tracking will work normally, but it will not be possible to use serial-dependent software, such as Mamehooker.

**Exiting Configuration Mode**
Once the parameters are configured and saved, it is highly recommended to reboot the lightgun to return to normal operating mode. Although the peripheral continues to function (except for the serial limitation in USB), a reboot is recommended to:
*   **Free up resources:** by terminating the webserver, all the microcontroller's processing power is dedicated exclusively to the tracking algorithms.
*   **Improve radio stability:** disabling the Wi-Fi Access Point eliminates potential radio interference.

---

<p align="center"> 🔸 🔸 🔸 </p>

---

<a id="versione-italiana"></a>

<p align="center">
  <a href="#english-version"><img src="../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# Lightgun Firmware (ESP32-S3) - NOVITÀ

## OpenFIRE ESP32 - ver 7.0.0 Release Candidate 1 (RC1)

Il core di OpenFIRE rimane praticamente invariato e garantisce una compatibilità totale (100%), tuttavia l'ecosistema circostante è stato sottoposto a una profonda revisione tecnica. Questa 7.0.0 RC1 va a sostituire l'attuale versione stabile 6.2.1. 

Di seguito vengono dettagliate le novità introdotte.

### Tracking, Sensori e Calibrazione
*   **Gestione LED persi migliorata:** Il sistema di tracking è stato affinato per gestire in modo più efficiente la perdita del segnale dei LED IR, una condizione frequente in caso di estrema vicinanza al monitor o di movimenti rapidi. Rimane pienamente supportata la gestione della rotazione della telecamera fino a circa ±80 gradi.
*   **Posizionamento libero dei sensori (Configurazione Square):** Oltre alla formazione standard consigliata (rettangolo verticale con base minore dell'altezza), i trasmettitori IR possono ora essere posizionati direttamente agli angoli del monitor. Il riconoscimento avviene in modo automatico, sebbene la formazione "Square" ideale rimanga la configurazione preferibile per ottenere le massime prestazioni (naturalmente va rieffettuata la calibrazione ogni volta che si riposizionano i trasmettitori IR).
*   **Supporto Multi-Cam Dinamico:** Il firmware ora supporta tre differenti moduli fotocamera, selezionabili dinamicamente dall'App di configurazione:
    1.  DFRobot SEN0158 / Wii Cam (connessione I2C - 2 cavi)
    2.  PixArt PAJ7025R2 (connessione SPI - 4 cavi)
    3.  PixArt PAJ7025R3 con obiettivo grandangolare integrato (connessione SPI - 4 cavi)

> **ATTENZIONE - Scelta dei LED IR:** 
> Per le telecamere DFRobot / Wii Cam è necessario utilizzare LED IR da 940nm. 
> Per i sensori PixArt PAJ7025R2 e R3 è tassativo l'uso di LED IR da 850nm. Sebbene i sensori PixArt riescano a rilevare LED da 940nm se dotati di potenza elevata, si tratta di componenti progettati specificamente per la lunghezza d'onda di 850nm. L'utilizzo di frequenze diverse degraderà in modo significativo la qualità del tracking e la massima distanza operativa.

*   **Visualizzazione realistica:** L'interfaccia di test della telecamera e della calibrazione (sia sul display del microcontrollore che sull'App) è stata leggermente modificata. L'immagine risulta ora più proporzionata e fedele alla geometria reale.

### Transizione alla WebApp Integrata
Questa release introduce un cambiamento strutturale nell'interfaccia di configurazione. L'App di configurazione desktop basata su Qt è stata aggiornata (implementando il nuovo protocollo, il supporto multi-cam e il pedale wireless), ma per il progetto OpenFIRE ESP32 è ora da considerarsi "Legacy". 

Il nuovo strumento ufficiale di configurazione è una **WebApp** sviluppata ex novo. L'interfaccia riprende l'App Qt 1:1 nel layout e nelle funzionalità, eliminando la necessità di scaricare o installare alcun software sul PC.
*   **Configurazione Online:** Accedendo al portale web dedicato, il sistema identifica automaticamente la versione del firmware installata sulla lightgun e avvia l'istanza corretta della WebApp (https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/).
*   **Configurazione Offline (Integrata nel Microcontrollore):** L'intero pacchetto della WebApp è archiviato fisicamente all'interno del microcontrollore. Avviando il dispositivo in "Modalità Configurazione", il microcontrollore espone un Access Point Wi-Fi (IP: 192.168.4.1) oppure, se collegato via cavo USB OTG, crea una rete virtuale NCM netUSB (raggiungibile all'IP 192.168.7.1 o digitando "openfire" nel browser).
*   **Flessibilità operativa:** La WebApp di configurazione risulta pienamente operativa anche quando la lightgun è interfacciata al PC in modalità wireless tramite DONGLE ESP-NOW. Le combinazioni di utilizzo sono quasi totalmente supportate.
*   *Nota di sviluppo:* Il porting dall'ambiente Qt alla WebApp ha richiesto un notevole impiego di risorse. Sono in corso gli ultimi perfezionamenti sui testi dell'interfaccia, ma a livello funzionale il sistema è pienamente operativo. Il codice del repository upstream (Board) continuerà ad essere monitorato per garantire la massima compatibilità e per il porting di eventuali future implementazioni.

### Comunicazione e Flashing
*   **Nuovo Protocollo Seriale:** La comunicazione tra l'applicativo di configurazione e la lightgun è gestita tramite un protocollo seriale riscritto, strutturato per garantire maggiore robustezza, tolleranza agli errori e scalabilità per futuri aggiornamenti.
*   **Aggiornamento Librerie:** Le librerie di base del firmware sono state aggiornate alle versioni più recenti.
*   **Flashing tramite "Magic Key":** Il firmware è ora in grado di interpretare il comando di flashing via software, disconnettendo automaticamente l'USB OTG. Non è più necessario accedere fisicamente al microcontrollore per premere i pulsanti fisici Boot/Reset. La procedura richiederà due tentativi (il tempo necessario alla porta seriale per la renumerazione post-riavvio), ma andrà a buon fine al secondo avvio. Rimane possibile forzare la modalità di flashing tenendo premuta l'apposita combinazione di tasti per 2 secondi all'avvio del dispositivo.

### Miglioramenti Generali ed Ecosistema
*   **Navigazione Menu:** A seguito del suggerimento di un utente, è ora possibile scorrere i menu in modalità pausa utilizzando anche i pulsanti direzionali "Su" e "Giù", in aggiunta ai classici bottoni "A" e "B".
*   Sono stati implementati vari bugfix minori e ottimizzazioni del codice sorgente.
*   **Hub Centrale:** È stato pubblicato un portale web unificato che raccoglie l'intero ecosistema di OpenFIRE ESP32, disponibile all'indirizzo: https://alessandro-satanassi.github.io/OpenFIRE-ESP32/

Considerata l'entità e la profondità delle modifiche apportate al codice, si richiede la cortesia di testare a fondo questa RC1 e di segnalare eventuali anomalie o suggerimenti prima del rilascio della versione stabile ("Stable"). 

**Si raccomanda vivamente di eseguire un'installazione pulita (Erase Flash)**. Sebbene i file di configurazione delle versioni precedenti siano teoricamente compatibili e tendano ad allinearsi dopo i primi salvataggi, partire da una memoria pulita previene potenziali conflitti.

## AVVIO IN MODALITÀ SPECIALI

Tenendo premute specifiche combinazioni di tasti per almeno 2 secondi durante l'avvio, la lightgun entrerà in una delle seguenti modalità speciali:

**Tasto 'TRIGGER' + Tasto 'A' (Modalità Aggiornamento Firmware)**
Il microcontrollore si predispone per la scrittura (flashing) del firmware. Sul display OLED comparirà il messaggio: *"Ready for firmware update"*. A questo punto è possibile procedere all'aggiornamento utilizzando il [Web Flasher](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/) oppure l'utility `esptool` (per gli utenti più esperti).

**Tasto 'B' (Modalità Configurazione Offline via WebApp interna)**
Questa modalità avvia un webserver integrato nel firmware, permettendo di configurare la lightgun tramite una WebApp ospitata localmente, senza alcuna necessità di connessione Internet. L'avvenuta attivazione è confermata sul display OLED: la barra superiore avrà i colori invertiti e mostrerà l'icona di un ingranaggio in alto a destra.

Il comportamento del sistema varia a seconda di come la lightgun è collegata:

**1. Collegamento tramite DONGLE (Wireless)**
*   Se avviata in questa modalità senza cavo USB, la lightgun scansionerà i canali di rete in attesa di agganciarsi a un Dongle disponibile. Questo è fondamentale per poter testare il puntamento "in game" anche in modalità *undocked*.
*   **Accesso alla WebApp:** Avviene tramite l'Access Point Wi-Fi esposto direttamente dal microcontrollore, denominato `"OpenFIRE_Config"` (IP: `192.168.4.1`). Collegandosi a questa rete (anche da smartphone), la WebApp dovrebbe aprirsi in automatico nel browser.
*   **Funzionalità:** Il sistema è pienamente operativo. Poiché la porta seriale virtuale rimane attiva, è supportato anche l'uso di Mamehooker.
*   *Nota tecnica:* Essendo la seriale attiva, potresti usare anche la WebApp online tramite Web Serial API all'indirizzo [OpenFIRE-ESP32-WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/). Tuttavia, la modalità "Tasto B" è pensata proprio per offrire un'alternativa quando non si ha internet a disposizione.

**2. Collegamento tramite cavo USB OTG**
*   **Accesso alla WebApp:** Oltre al normale Access Point Wi-Fi (`192.168.4.1`), il sistema crea una rete virtuale NCM netUSB. Questo permette di accedere all'interfaccia semplicemente digitando `openfire` oppure l'IP `192.168.7.1` nella barra degli indirizzi del browser, senza la necessità di collegarsi all'AP del micro.
*   **Funzionalità:** Per limitazioni tecniche legate al numero massimo di endpoint USB disponibili, la porta seriale viene disabilitata. Il puntamento della lightgun funzionerà regolarmente, ma non sarà possibile utilizzare software dipendenti dalla seriale, come Mamehooker.

**Uscita dalla Modalità di Configurazione**
Una volta configurati e salvati i parametri, è fortemente raccomandato riavviare la lightgun per tornare alla normale modalità operativa. Sebbene la periferica continui a funzionare (salvo il limite della seriale in USB), il riavvio è consigliato per:
*   **Liberare risorse:** terminando il webserver, tutta la potenza di calcolo del microcontrollore viene dedicata esclusivamente agli algoritmi di puntamento.
*   **Migliorare la stabilità radio:** disattivando l'Access Point Wi-Fi si eliminano potenziali interferenze radio.