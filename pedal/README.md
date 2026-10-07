<a id="english-version"></a>

[Back to Home](../README.md#english-version) / **Pedal Firmware**

<p align="center">
  <a href="#english-version"><img src="../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# Pedal Firmware (ESP32-S3)

<p align="center">
  <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases"><img src="https://img.shields.io/github/v/release/alessandro-satanassi/OpenFIRE-Firmware-ESP32?include_prereleases&style=flat-square&color=007ec6&label=latest%20version" alt="Latest Version"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/pedal"><img src="https://img.shields.io/github/languages/top/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&color=success" alt="Top Language"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/pedal"><img src="https://img.shields.io/badge/Platform-PlatformIO-orange?style=flat-square&logo=platformio" alt="PlatformIO"></a> <a href="../README.md#community-support-english"><img src="https://img.shields.io/badge/Discord-Community-5865F2?style=flat-square&logo=discord&logoColor=white" alt="Discord Community"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/LICENSE"><img src="https://img.shields.io/github/license/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square" alt="License"></a>
</p>

<p align="center">
  <img src="docs/img/copertina_pedal.jpg" alt="OpenFIRE Firmware ESP32 pedal" width="100%">
</p>

---
> **Hardware sponsored by [PCBWay](https://www.pcbway.com)**
---

The wireless pedal is made for arcade cover shooters such as the *Time Crisis* series, where you press the pedal to come out and shoot and release it to take cover. With no cable you can place it wherever you like on the floor, and the ESP-NOW radio keeps the response as fast as these games demand. There is nothing to configure on the computer: the gun receives the pedal presses and passes them on as if the pedal were wired to it.

## Supported Hardware and Wiring

The pedal is built on a standard development board. The very small **Waveshare ESP32-S3-ZERO** is the best fit inside a pedal; the **Waveshare ESP32-S3-PICO** and the classic **ESP32-S3-DevKitC-1** work too.

**Antenna:** keep a small clear area around the board's antenna (the printed antenna at the end of the ESP32-S3 module). Do not route wires over it or right next to it, and do not cover it with metal: wires touching the antenna greatly reduce the wireless range and reliability.

### Pedal Components

The firmware handles up to two pedals and shows its status on 4 LEDs.

| Component | Type | Description |
| :--- | :---: | :--- |
| **Pedal 1 (Main)** | **Mandatory** | The switch operated by the main pedal. |
| **Pedal 2 (Secondary)** | *Optional* | A second input, for games that use two pedals or for an extra function. |
| **4x Status LEDs** | *Optional* | Four standard LEDs showing startup, search and connection, then the player number. |

> [!NOTE]
> **The LEDs and the second pedal are optional.** If you leave them out, simply leave their pins unconnected.

Pinouts of the supported boards for use as a pedal:

| Waveshare ESP32-S3-PICO | Waveshare ESP32-S3-ZERO | ESP32-S3-DevKitC-1 |
| :---: | :---: | :---: |
| <img src="docs/board_scheme/PEDAL-esp32-s3-pico.svg" width="100%" alt="Pedal PICO Pinout"> | <img src="docs/board_scheme/PEDAL-esp32-s3-zero.svg" width="100%" alt="Pedal ZERO Pinout"> | <img src="docs/board_scheme/PEDAL-ESP32S3-Devkit-C.svg" width="100%" alt="Pedal DevKitC Pinout"> |

<br>

<p align="center">
  <b>Wiring example on a breadboard (Waveshare ESP32-S3-ZERO):</b><br><br>
  <img src="docs/img/pedal_bb.png" width="80%" alt="Breadboard wiring scheme for Pedal ZERO">
</p>

---

## Firmware Installation and Flashing

In the Web Flasher choose **Pedal**, your exact board and its flash/PSRAM variant. Use firmware from the same release as the gun (7.0.0). The pedal is programmed through its own USB port; the gun's Trigger + A and B startup shortcuts do not apply to it, so if it does not enter update mode by itself, use its BOOT/RESET buttons.

#### WEB FLASHER (Recommended for all users)
The easiest, fastest and safest way to install or update the firmware. Nothing to install: it runs in your browser.
* **Requirements:** a computer with Google Chrome, Microsoft Edge or Opera.

**[LAUNCH OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?lang=en)**

---

#### *FOR ADVANCED USERS:*
If your browser does not support the Web Flasher, or you prefer the command line, you can install the firmware files yourself, exactly as for the gun and the dongle.

### Simplified Procedure with Script
The fastest way, with nothing else to install on the computer.

1. Go to the **[Releases](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases)** page.
2. Download the "Simplified Procedure" ZIP for your board and operating system (e.g. `OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8-windows-64bit.zip`).
   > *Make sure you choose the file of your exact board.*
3. Extract the whole ZIP into a folder on your computer.
4. Connect the pedal board to the computer with a USB cable.
5. Run the `flash_firmware` script (`.bat` on Windows, `.sh` on Linux and macOS) and follow the instructions on screen.

### Troubleshooting
* **The installation does not start (`Connecting...`):** some boards do not enter update mode by themselves. If the script keeps repeating `Connecting...`, hold the small **BOOT** (or `B`) button on the board until the installation starts.
* **Antivirus false positive (Windows):** the script uses Espressif's original `esptool.exe`, which some antivirus programs may block or flag. Download the package from the official Releases page and check its origin before allowing it to run; do not disable your antivirus globally.
* **Manual installation:** the single `.bin` file of each board is also on the Releases page, for esptool or graphical tools such as NodeMCU PyFlasher (address `0x0`).

---

## Boot and Synchronization Sequence (Pairing)

**Before pairing:** in the WebApp, enable **Wireless pedal** in the gun's settings, leave both wired pedal inputs unmapped, save and restart the gun. A wired pedal mapped to either input prevents the wireless pedal from being found. The wireless pedal works when the gun plays through the dongle; with the gun connected to the computer by USB cable it is not searched, so use a wired pedal instead.

The pedal has no pairing button: it pairs with the gun during a short window when the gun starts.

1. **Pedal on:** as soon as it is powered (battery or power bank), the pedal waits for a gun *(with the 4 LEDs connected, a light runs back and forth, like KITT from Knight Rider)*.
2. **Search window:** when you switch on the gun, it first connects to the dongle; right after that, with the option above enabled and no wired pedal mapped, it looks for a pedal for about **10 seconds**.
3. **Pairing:** if the pedal is on and in range during those 10 seconds, the gun finds it, pairs with it exclusively and stops searching *(the LEDs confirm the connection, then LED 1, 2, 3 or 4 stays on, showing the gun's player number)*.
4. **Play:** from now on every press reaches the gun instantly, and the gun sends it to the computer together with the trigger and the aim. The two pedal inputs work like the gun's **Pedal** and **Alt Pedal** buttons: mouse button 4 and mouse button 5 by default, remappable in the WebApp's **Button Mapping** tab.
5. **Fast reconnection:** if you switch the gun off while the pedal and the dongle stay on, the pedal reconnects almost instantly when you switch the gun back on, without the 10-second search.

> [!IMPORTANT]
> **Pedal restart and new pairing**
>
> Every time the pedal is switched off and on again (or loses power), it goes back to step 1 and waits for a gun.
> *So after restarting the pedal, **switch the gun off and on again** as well, so that it opens the 10-second search again.*

---
### Questions or Issues?
For technical support and to join the discussion, see the [Community & Support section](../README.md#community-support-english) on the project's home page.

---

<a id="versione-italiana"></a>

[Torna alla Home](../README.md#versione-italiana) / **Pedal Firmware**

<p align="center">
  <a href="#english-version"><img src="../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# Pedal Firmware (ESP32-S3)

<p align="center">
  <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases"><img src="https://img.shields.io/github/v/release/alessandro-satanassi/OpenFIRE-Firmware-ESP32?include_prereleases&style=flat-square&color=007ec6&label=ultima%20versione" alt="Ultima Versione"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/pedal"><img src="https://img.shields.io/github/languages/top/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&color=success" alt="Linguaggio Principale"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/pedal"><img src="https://img.shields.io/badge/Platform-PlatformIO-orange?style=flat-square&logo=platformio" alt="PlatformIO"></a> <a href="../README.md#community-support-italiano"><img src="https://img.shields.io/badge/Discord-Community-5865F2?style=flat-square&logo=discord&logoColor=white" alt="Discord Community"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/LICENSE"><img src="https://img.shields.io/github/license/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&label=licenza" alt="Licenza"></a>
</p>

<p align="center">
  <img src="docs/img/copertina_pedal.jpg" alt="OpenFIRE Firmware ESP32 pedal" width="100%">
</p>

---
> **Hardware sponsored by [PCBWay](https://www.pcbway.com)**
---

Il pedale wireless è pensato per i cover shooter arcade come la serie *Time Crisis*, dove premi il pedale per uscire allo scoperto e sparare e lo rilasci per ripararti. Senza cavo puoi metterlo dove vuoi sul pavimento, e la radio ESP-NOW mantiene la risposta rapida quanto questi giochi richiedono. Sul computer non c'è nulla da configurare: la pistola riceve le pressioni del pedale e le inoltra come se il pedale fosse collegato a lei via cavo.

## Hardware Supportato e Cablaggio

Il pedale si costruisce su una scheda di sviluppo standard. La piccolissima **Waveshare ESP32-S3-ZERO** è la più adatta da inserire in un pedale; funzionano anche la **Waveshare ESP32-S3-PICO** e la classica **ESP32-S3-DevKitC-1**.

**Antenna:** lascia libera una piccola area intorno all'antenna della scheda (l'antenna stampata all'estremità del modulo ESP32-S3). Non far passare fili sopra o a ridosso e non coprirla con metallo: fili che toccano l'antenna riducono molto la portata e l'affidabilità del collegamento wireless.

### Componenti del Pedale

Il firmware gestisce fino a due pedali e mostra il suo stato su 4 LED.

| Componente | Tipo | Descrizione |
| :--- | :---: | :--- |
| **Pedale 1 (Principale)** | **Obbligatorio** | L'interruttore azionato dal pedale principale. |
| **Pedale 2 (Secondario)** | *Opzionale* | Un secondo ingresso, per i giochi che usano due pedali o per una funzione in più. |
| **4x LED di Stato** | *Opzionali* | Quattro LED standard che mostrano avvio, ricerca e connessione, poi il numero del giocatore. |

> [!NOTE]
> **I LED e il secondo pedale sono facoltativi.** Se non li monti, lascia semplicemente scollegati i relativi pin.

Pinout delle schede supportate per l'uso come pedale:

| Waveshare ESP32-S3-PICO | Waveshare ESP32-S3-ZERO | ESP32-S3-DevKitC-1 |
| :---: | :---: | :---: |
| <img src="docs/board_scheme/PEDAL-esp32-s3-pico.svg" width="100%" alt="Pinout Pedale PICO"> | <img src="docs/board_scheme/PEDAL-esp32-s3-zero.svg" width="100%" alt="Pinout Pedale ZERO"> | <img src="docs/board_scheme/PEDAL-ESP32S3-Devkit-C.svg" width="100%" alt="Pinout Pedale DevKitC"> |

<br>

<p align="center">
  <b>Esempio di cablaggio su breadboard (Waveshare ESP32-S3-ZERO):</b><br><br>
  <img src="docs/img/pedal_bb.png" width="80%" alt="Schema di cablaggio su breadboard per Pedale ZERO">
</p>

---

## Installazione e Flashing del Firmware

Nel Web Flasher scegli **Pedal**, la tua scheda esatta e la sua variante flash/PSRAM. Usa firmware della stessa release della pistola (7.0.0). Il pedale si programma tramite la sua porta USB; le scorciatoie all'avvio della pistola (Grilletto + A e B) non valgono per il pedale, quindi se non entra da solo in modalità aggiornamento usa i suoi pulsanti BOOT/RESET.

#### WEB FLASHER (Consigliato per qualsiasi utente)
Il modo più semplice, veloce e sicuro per installare o aggiornare il firmware. Niente da installare: funziona nel browser.
* **Requisiti:** un computer con Google Chrome, Microsoft Edge o Opera.

**[AVVIA OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?lang=it)**

---

#### *PER UTENTI ESPERTI:*
Se il tuo browser non supporta il Web Flasher, o preferisci la riga di comando, puoi installare manualmente i file del firmware, esattamente come per la pistola e il dongle.

### Procedura Semplificata con Script
Il metodo più veloce, senza nient'altro da installare sul computer.

1. Vai alla pagina delle **[Releases](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases)**.
2. Scarica lo ZIP "Procedura Semplificata" per la tua scheda e il tuo sistema operativo (es. `OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8-windows-64bit.zip`).
   > *Assicurati di scegliere il file della tua scheda esatta.*
3. Estrai tutto lo ZIP in una cartella del computer.
4. Collega la scheda del pedale al computer con un cavo USB.
5. Esegui lo script `flash_firmware` (`.bat` su Windows, `.sh` su Linux e macOS) e segui le istruzioni a schermo.

### Risoluzione dei Problemi
* **L'installazione non parte (`Connecting...`):** alcune schede non entrano da sole in modalità aggiornamento. Se lo script continua a ripetere `Connecting...`, tieni premuto il piccolo pulsante **BOOT** (o `B`) della scheda finché l'installazione non parte.
* **Falso positivo dell'antivirus (Windows):** lo script usa l'`esptool.exe` originale di Espressif, che alcuni antivirus possono bloccare o segnalare. Scarica il pacchetto dalla pagina Releases ufficiale e verificane la provenienza prima di consentirne l'esecuzione; non disattivare l'antivirus globalmente.
* **Installazione manuale:** nella pagina Releases c'è anche il singolo file `.bin` di ogni scheda, per esptool o programmi grafici come NodeMCU PyFlasher (indirizzo `0x0`).

---

## Sequenza di Avvio e Sincronizzazione (Pairing)

**Prima dell'associazione:** nella WebApp attiva **Pedale wireless** nelle impostazioni della pistola, lascia non assegnati entrambi gli ingressi dei pedali cablati, salva e riavvia la pistola. Un pedale cablato assegnato a uno dei due ingressi impedisce di trovare il pedale wireless. Il pedale wireless funziona quando la pistola gioca tramite il dongle; con la pistola collegata al computer via cavo USB non viene cercato, quindi usa un pedale cablato.

Il pedale non ha un pulsante di associazione: si associa alla pistola durante una breve finestra all'accensione della pistola.

1. **Pedale acceso:** appena riceve alimentazione (batteria o power bank), il pedale attende una pistola *(con i 4 LED collegati, una luce scorre avanti e indietro, come KITT di Supercar)*.
2. **Finestra di ricerca:** quando accendi la pistola, questa si collega prima al dongle; subito dopo, con l'opzione sopra attiva e nessun pedale cablato assegnato, cerca un pedale per circa **10 secondi**.
3. **Associazione:** se il pedale è acceso e a portata durante quei 10 secondi, la pistola lo trova, vi si associa in modo esclusivo e smette di cercare *(i LED confermano la connessione, poi resta acceso il LED 1, 2, 3 o 4, che indica il numero del giocatore della pistola)*.
4. **Gioco:** da questo momento ogni pressione arriva all'istante alla pistola, che la invia al computer insieme al grilletto e alla mira. I due ingressi del pedale funzionano come i pulsanti **Pedale** e **Pedale alternativo** della pistola: per impostazione predefinita tasto mouse 4 e tasto mouse 5, rimappabili nella scheda **Mappatura Pulsanti** della WebApp.
5. **Riconnessione rapida:** se spegni la pistola mentre pedale e dongle restano accesi, alla riaccensione della pistola il pedale si ricollega quasi all'istante, senza i 10 secondi di ricerca.

> [!IMPORTANT]
> **Riavvio del pedale e nuova associazione**
>
> Ogni volta che il pedale viene spento e riacceso (o perde l'alimentazione), torna al passo 1 e attende una pistola.
> *Quindi, dopo aver riavviato il pedale, **spegni e riaccendi anche la pistola**, perché riapra i 10 secondi di ricerca.*

---
### Domande o Problemi?
Per supporto tecnico e per unirti alla community, consulta la [sezione Community e Supporto](../README.md#community-support-italiano) nella Home del progetto.
