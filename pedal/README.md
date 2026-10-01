<a id="english-version"></a>

[Back to Home](../README.md#english-version) / **Pedal Firmware**

<p align="center">
  <a href="#english-version"><img src="../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# Pedal Firmware (ESP32-S3)

<p align="center">
  <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases"><img src="https://img.shields.io/github/v/release/alessandro-satanassi/OpenFIRE-Firmware-ESP32?include_prereleases&style=flat-square&color=007ec6&label=latest%20version" alt="Latest Version"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/pedal"><img src="https://img.shields.io/github/languages/top/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&color=success" alt="Top Language"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/pedal"><img src="https://img.shields.io/badge/Platform-PlatformIO-orange?style=flat-square&logo=platformio" alt="PlatformIO"></a> <a href="../README.md#community-support-english"><img src="https://img.shields.io/badge/Discord-Community-5865F2?style=flat-square&logo=discord&logoColor=white" alt="Discord Community"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/LICENSE"><img src="https://img.shields.io/github/license/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square" alt="License"></a>
</p>

<p align="center">
  <img src="docs/img/copertina_pedal.jpg" alt="OpenFIRE Firmware ESP32 pedal" width="100%">
</p>

---
> **Hardware sponsored by [PCBWay](https://www.pcbway.com)**
---

The wireless pedal is designed for arcade "cover shooter" games, such as the *Time Crisis* series, where the pedal is used to take cover and to come back out and shoot. 
Operating in a completely wireless mode, it allows the user to position it freely on the floor, eliminating any annoying physical constraints with the PC or the gun. To meet the strict responsiveness requirements imposed by arcade games, the pedal leverages the power of the ESP-NOW protocol, guaranteeing ultra-low latency. 
The system is totally "Plug & Play" and does not require any additional management software: the Lightgun instantly receives the wireless signal from the pedal and translates it into a standard input, exactly as if it were physically connected via cable.

## Supported Hardware and Wiring

To build the Wireless Pedal, you will need to wire the components to a standard development board. You can use boards like the ultra-compact **Waveshare ESP32-S3-ZERO** (highly recommended due to space constraints inside the pedal chassis), the **Waveshare ESP32-S3-PICO**, or the classic **ESP32-S3-DevKitC-1**.

**Antenna:** keep a small clear area around the board's antenna (the printed antenna at the end of the ESP32-S3 module). Do not route wires over it or right next to it, and do not cover it with metal parts: wires touching the antenna greatly reduce the wireless range and reliability.

### Pedal Components

The firmware manages up to two physical pedals and provides visual feedback via 4 LEDs.

| Component | Type | Description |
| :--- | :---: | :--- |
| **Pedal 1 (Main)** | **Mandatory** | The microswitch or physical pedal button. Without this input, the device has no utility. |
| **Pedal 2 (Secondary)** | *Optional* | A second input for dual-pedal configurations (e.g., for specific arcade titles or extra functions). |
| **4x Status LEDs** | *Optional* | Four standard LEDs to display boot status, searching, and successful connection. |

> [!NOTE]
> **The LEDs and the second pedal are not mandatory for basic operation.** If you decide to omit them, make sure to leave the corresponding contacts disconnected. 

Here are the graphical pinouts of the various boards for use as a Pedal:

| Waveshare ESP32-S3-PICO | Waveshare ESP32-S3-ZERO | ESP32-S3-DevKitC-1 |
| :---: | :---: | :---: |
| <img src="docs/board_scheme/PEDAL-esp32-s3-pico.svg" width="100%" alt="Pedal PICO Pinout"> | <img src="docs/board_scheme/PEDAL-esp32-s3-zero.svg" width="100%" alt="Pedal ZERO Pinout"> | <img src="docs/board_scheme/PEDAL-ESP32S3-Devkit-C.svg" width="100%" alt="Pedal DevKitC Pinout"> |

<br>

<p align="center">
  <b>Practical wiring example on a Breadboard (Waveshare ESP32-S3-ZERO):</b><br><br>
  <img src="docs/img/pedal_bb.png" width="80%" alt="Breadboard wiring scheme for Pedal ZERO">
</p>

---

## Firmware Installation and Flashing

Choose **Pedal**, your exact board and its flash/PSRAM variant in the Web Flasher. Use firmware from the same release as the lightgun (7.0.0). Flash this device through its own USB port; the lightgun’s Trigger + A / B startup shortcuts do not apply to this board. Use its physical BOOT/RESET controls if download mode is needed.

#### WEB FLASHER (Recommended for all users)
The easiest, fastest, and safest way to install or update the firmware. It does not require installing any drivers or external software: it runs entirely within your browser.
* **Requirements:** PC/Mac with Google Chrome, Microsoft Edge, or Opera.

**[LAUNCH OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?lang=en)**

---

#### *FOR ADVANCED USERS:*
If your browser does not support the Web Flasher, or if you prefer to proceed via command line or external tools, you can download and install the specific files.

The firmware loading process (flashing) is identical to the one for the Lightgun and the Dongle. 

### Simplified Procedure with Script
This is the fastest method and does not require the installation of additional software on the computer.

1. Go to the **[Releases](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases)** page.
2. Download the "Simplified Procedure" ZIP for your specific pedal model and operating system (e.g., `OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8-windows-64bit.zip`). 
   > *Warning: make sure to choose the exact file for your board.*
3. Extract the entire content of the ZIP archive into a folder on your PC.
4. Connect the device to the PC via USB port.
5. Run the `flash_firmware` script (on Windows it will be the `.bat` file) and follow the on-screen instructions.

### Troubleshooting
* **Flashing won't start (Connecting...):** Some boards can be reluctant to automatically enter download mode. If the script gets stuck repeating the word `Connecting...`, press and hold the small physical **BOOT** (or `B`) button on the device until the installation begins.
* **Antivirus False Positive (Windows):** The script uses the original `esptool.exe` utility by Espressif. Some antivirus software might block it or flag it as a false positive. The file is 100% safe; you may need to temporarily add it to your exceptions.
* **Manual Installation:** For advanced users, individual `.bin` files are also provided on the Release page to be flashed manually using graphical tools like NodeMCU PyFlasher.

---

## Boot and Synchronization Sequence (Pairing)

**Before pairing:** in the lightgun’s WebApp, enable **Wireless pedal**, leave both wired pedal inputs unmapped, save and restart the gun. A wired pedal mapped to either input prevents wireless pedal discovery. If you do not use a wireless pedal, disable this option to skip the search. The wireless pedal works only with a gun connected through the dongle; when the gun is connected to the computer by USB cable, the pedal is not searched, so use a wired pedal instead.

The Pedal operates completely autonomously and without pairing switches. Pairing with the Lightgun occurs within a specific "time window" at startup.

Here is what happens to establish the connection:

1. **Pedal Power-on:** As soon as it is powered (via battery or power bank), the Pedal enters a "silent listening" state *(if you have connected the 4 LEDs, you will see an animation indicating the waiting state, with the LEDs scanning back and forth like KITT from Knight Rider)*.
2. **Lightgun Search Window:** When you turn on your wireless Lightgun, it first connects to the Dongle on the PC. Immediately after, if the option above is enabled and neither wired pedal input is mapped, it opens a search window of about **10 seconds** dedicated to the pedal.
3. **Pairing:** If the Pedal is powered on and within range during these 10 seconds, the Lightgun detects it, "locks" onto it exclusively, and stops searching *(the LEDs on the pedal will confirm the successful connection, and LED 1, 2, 3, or 4 will remain solidly lit depending on the player assigned to the lightgun)*.
4. **Operation:** From this moment on, physical presses on the pedal are instantly transmitted to the gun, which processes them and sends them to the PC along with the trigger and optical tracking data. The two pedal inputs work like the gun's **Pedal** and **Alt Pedal** buttons: by default mouse button 4 and mouse button 5, remappable in the WebApp's **Button Mapping** tab.
5. **Fast Reconnection:** If for some reason you turn off the gun, but leave the Pedal and the Dongle (both already paired) powered on, as soon as you turn the lightgun back on, the reconnection to the pedal will be instantaneous, skipping the 10-second search.

> [!IMPORTANT]
> **Pedal Restart and New Pairing**
> Every time you turn the Pedal off and back on (or disconnect its power supply), it resets its memory and returns to Step 1, waiting for a connection.
> *Consequently, if the Pedal is restarted, you must necessarily **turn off and on your Lightgun as well** to force it to reopen the famous "10-second window" and allow it to transmit the initial connection signal again.*

---
### Questions or Issues?
For technical support and to join the discussion, please refer to the [Community & Support Section](../README.md#community-support-english) in the Main Repository.

---

<a id="versione-italiana"></a>

[Torna alla Home](../README.md#versione-italiana) / **Pedal Firmware**

<p align="center">
  <a href="#english-version"><img src="../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# Pedal Firmware (ESP32-S3)

<p align="center">
  <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases"><img src="https://img.shields.io/github/v/release/alessandro-satanassi/OpenFIRE-Firmware-ESP32?include_prereleases&style=flat-square&color=007ec6&label=ultima%20versione" alt="Ultima Versione"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/pedal"><img src="https://img.shields.io/github/languages/top/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&color=success" alt="Linguaggio Principale"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/pedal"><img src="https://img.shields.io/badge/Platform-PlatformIO-orange?style=flat-square&logo=platformio" alt="PlatformIO"></a> <a href="../README.md#community-support-italiano"><img src="https://img.shields.io/badge/Discord-Community-5865F2?style=flat-square&logo=discord&logoColor=white" alt="Discord Community"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/LICENSE"><img src="https://img.shields.io/github/license/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&label=licenza" alt="Licenza"></a>
</p>

<p align="center">
  <img src="docs/img/copertina_pedal.jpg" alt="OpenFIRE Firmware ESP32 pedal" width="100%">
</p>

---
> **Hardware sponsored by [PCBWay](https://www.pcbway.com)**
---

Il pedale wireless è pensato per i giochi arcade di tipo "cover shooter", come la serie *Time Crisis*, in cui il pedale serve a mettersi al riparo e a uscire per sparare. 
Operando in modalità completamente senza fili, consente all'utente di posizionarlo liberamente sul pavimento, eliminando qualsiasi fastidioso vincolo fisico con il PC o con la pistola. Per soddisfare i severi requisiti di reattività imposti dai giochi arcade, il pedale sfrutta la potenza del protocollo ESP-NOW, garantendo latenze bassissime. 
Il sistema è totalmente "Plug & Play" e non richiede alcun software di gestione aggiuntivo: la Lightgun riceve istantaneamente il segnale wireless dal pedale e lo traduce in un input standard, esattamente come se quest'ultimo fosse collegato fisicamente via cavo.

## Hardware Supportato e Cablaggio

Per realizzare il Pedale Wireless, dovrai cablare i componenti su una scheda di sviluppo standard. Puoi utilizzare board come la compattissima **Waveshare ESP32-S3-ZERO** (la più consigliata per motivi di spazio all'interno dello chassis del pedale), la **Waveshare ESP32-S3-PICO** o la classica **ESP32-S3-DevKitC-1**.

**Antenna:** lascia libera una piccola area intorno all'antenna della scheda (l'antenna stampata all'estremità del modulo ESP32-S3). Non far passare fili sopra o a ridosso e non coprirla con parti metalliche: fili che toccano l'antenna riducono molto la portata e l'affidabilità del collegamento wireless.

### Componenti del Pedale

Il firmware gestisce fino a due pedali fisici e fornisce un feedback visivo tramite 4 LED.

| Componente | Tipologia | Descrizione |
| :--- | :---: | :--- |
| **Pedale 1 (Principale)** | **Obbligatorio** | Il microswitch o pulsante fisico del pedale. Senza questo ingresso, il dispositivo non ha utilità. |
| **Pedale 2 (Secondario)** | *Opzionale* | Un secondo ingresso per configurazioni a doppio pedale (es. per titoli arcade specifici o funzioni extra). |
| **4x LED di Stato** | *Opzionali* | Quattro LED standard per visualizzare lo stato dell'avvio, la ricerca e l'avvenuta connessione. |

> [!NOTE]
> **I LED e il secondo pedale non sono obbligatori per il funzionamento di base.** Se decidi di ometterli, assicurati comunque di non collegare i relativi contatti. 

Ecco gli schemi grafici delle varie schede per l'utilizzo come Pedale:

| Waveshare ESP32-S3-PICO | Waveshare ESP32-S3-ZERO | ESP32-S3-DevKitC-1 |
| :---: | :---: | :---: |
| <img src="docs/board_scheme/PEDAL-esp32-s3-pico.svg" width="100%" alt="Pinout Pedal PICO"> | <img src="docs/board_scheme/PEDAL-esp32-s3-zero.svg" width="100%" alt="Pinout Pedal ZERO"> | <img src="docs/board_scheme/PEDAL-ESP32S3-Devkit-C.svg" width="100%" alt="Pinout Pedal DevKitC"> |

<br>

<p align="center">
  <b>Esempio pratico di cablaggio su Breadboard (Waveshare ESP32-S3-ZERO):</b><br><br>
  <img src="docs/img/pedal_bb.png" width="80%" alt="Schema di collegamento su breadboard per Pedal ZERO">
</p>

---

## Installazione e Flashing del Firmware

Nel Web Flasher scegli **Pedal**, la scheda esatta e la sua variante flash/PSRAM. Usa il firmware della stessa release della lightgun (7.0.0). Installa il firmware dalla porta USB del dispositivo stesso: le combinazioni Grilletto + A / B all’avvio della lightgun non si applicano a questa scheda. Se serve la modalità download, usa i suoi comandi fisici BOOT/RESET.

####  WEB FLASHER (Consigliato per qualsiasi utente)
Il modo più semplice, veloce e sicuro per installare o aggiornare il firmware. Non richiede l'installazione di driver o software esterni: viene eseguito interamente dal tuo browser.
* **Requisiti:** PC/Mac con Google Chrome, Microsoft Edge o Opera.

**[AVVIA OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?lang=it)**

---

#### *PER UTENTI ESPERTI:*
Se il tuo browser non supporta il Web Flasher, o se preferisci procedere tramite riga di comando o tool esterni, puoi scaricare ed installare i file specifici.

Il processo di caricamento del firmware (flashing) è identico a quello previsto per la Lightgun e per il Dongle. 

### Procedura Semplificata con Script
Questo è il metodo più veloce e non richiede l'installazione di software aggiuntivi sul computer.

1. Vai alla pagina delle **[Releases](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases)**.
2. Scarica lo ZIP "Procedura Semplificata" relativo al tuo specifico modello di pedale e al tuo sistema operativo (es. `OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8-windows-64bit.zip`). 
   > *Attenzione: assicurati di scegliere il file esatto per la tua scheda.*
3. Estrai l'intero contenuto dell'archivio ZIP in una cartella sul tuo PC.
4. Collega il dispositivo al PC tramite porta USB.
5. Esegui lo script `flash_firmware` (su Windows sarà il file `.bat`) e segui le istruzioni a schermo.

### Risoluzione dei Problemi (Troubleshooting)
* **Il Flashing non parte (Connecting...):** Alcune schede possono essere restie a entrare automaticamente in modalità download. Se lo script si blocca ripetendo la scritta `Connecting...`, tieni premuto il piccolo pulsante fisico **BOOT** (o `B`) presente sul dispositivo finché l'installazione non inizia.
* **Falso Positivo Antivirus (Windows):** Lo script utilizza l'utility `esptool.exe` originale di Espressif. Alcuni software antivirus potrebbero bloccarlo o segnalarlo come falso positivo. Il file è sicuro al 100%; potresti doverlo aggiungere momentaneamente alle eccezioni.
* **Installazione Manuale:** Per gli utenti avanzati, nella pagina delle Release sono forniti anche i singoli file `.bin` da flashare manualmente utilizzando tool grafici come NodeMCU PyFlasher.

---

## Sequenza di Avvio e Sincronizzazione (Pairing)

**Prima dell’associazione:** nella WebApp della lightgun abilita **Pedale wireless**, lascia non assegnati entrambi gli ingressi dei pedali cablati, salva e riavvia la pistola. Un pedale cablato mappato su uno dei due ingressi impedisce la ricerca del pedale wireless. Se non usi un pedale wireless, disabilita l’opzione per saltare la ricerca. Il pedale wireless funziona solo con una pistola collegata tramite il dongle; quando la pistola è collegata al computer con il cavo USB il pedale non viene cercato, quindi usa invece un pedale cablato.

Il Pedale opera in totale autonomia e senza interruttori di pairing. L'associazione con la Lightgun avviene in una specifica "finestra temporale" all'avvio.

Ecco cosa succede per instaurare la connessione:

1. **Accensione del Pedale:** Appena alimentato (tramite batteria o power bank), il Pedale si mette in stato di "ascolto silenzioso" *(se hai collegato i 4 LED, vedrai un'animazione che indica lo stato di attesa, con i LED che scorrono avanti e indietro come KITT di Supercar)*.
2. **Finestra di Ricerca della Lightgun:** Quando accendi la tua Lightgun wireless, quest'ultima si collega prima al Dongle sul PC. Immediatamente dopo, se l’opzione indicata sopra è abilitata e nessun ingresso del pedale cablato è mappato, apre una finestra di ricerca di circa **10 secondi** dedicata al pedale.
3. **Associazione (Pairing):** Se il Pedale è acceso e nel raggio d'azione durante questi 10 secondi, la Lightgun lo rileva, si "aggancia" in via esclusiva e interrompe la ricerca *(i LED sul pedale confermeranno l'avvenuta connessione e rimarrà acceso in modo fisso il LED 1, 2, 3 o 4 a seconda del player assegnato alla lightgun)*.
4. **Operatività:** Da questo momento, le pressioni fisiche sul pedale vengono trasmesse istantaneamente alla pistola, che le elabora e le invia al PC insieme al grilletto e al tracciamento ottico. I due ingressi del pedale funzionano come i pulsanti **Pedale** e **Pedale alternativo** della pistola: in modo predefinito tasto mouse 4 e tasto mouse 5, rimappabili nella scheda **Mappatura Pulsanti** della WebApp.
5. **Riconnessione Veloce:** Se per qualche motivo spegni la pistola, ma lasci il Pedale e il Dongle (entrambi già accoppiati) accesi, non appena riaccenderai la lightgun la riconnessione al pedale sarà istantanea, saltando i 10 secondi di ricerca.

> [!IMPORTANT]
> **Riavvio del Pedale e nuovo Pairing**
> Ogni volta che spegni e riaccendi il Pedale (o ne scolleghi l'alimentazione), questo resetta la sua memoria e torna al Punto 1, mettendosi in attesa.
> *Di conseguenza, in caso di riavvio del Pedale, dovrai necessariamente **spegnere e riaccendere anche la tua Lightgun** per forzarla a riaprire la famosa "finestra di 10 secondi" e permetterle di trasmettere di nuovo il segnale iniziale di connessione.*

---
### Domande o Problemi?
Per supporto tecnico e per unirti alla community, consulta la [Sezione Community e Supporto](../README.md#community-support-italiano) nella Home del progetto.
