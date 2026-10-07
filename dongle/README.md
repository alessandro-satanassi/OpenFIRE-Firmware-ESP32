<a id="english-version"></a>

[Back to Home](../README.md#english-version) / **Dongle Firmware**

<p align="center">
  <a href="#english-version"><img src="../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# Dongle Firmware (ESP32-S3)

<p align="center">
  <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases"><img src="https://img.shields.io/github/v/release/alessandro-satanassi/OpenFIRE-Firmware-ESP32?include_prereleases&style=flat-square&color=007ec6&label=latest%20version" alt="Latest Version"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/dongle"><img src="https://img.shields.io/github/languages/top/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&color=success" alt="Top Language"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/dongle"><img src="https://img.shields.io/badge/Platform-PlatformIO-orange?style=flat-square&logo=platformio" alt="PlatformIO"></a> <a href="../README.md#community-support-english"><img src="https://img.shields.io/badge/Discord-Community-5865F2?style=flat-square&logo=discord&logoColor=white" alt="Discord Community"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/LICENSE"><img src="https://img.shields.io/github/license/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square" alt="License"></a>
</p>

<p align="center">
  <img src="docs/img/copertina_dongle.jpg" alt="OpenFIRE Firmware ESP32 dongle" width="100%">
</p>

---
> **Hardware sponsored by [PCBWay](https://www.pcbway.com)**
---

The dongle is what lets you play without cables. Plug it into a USB port of the computer, switch on the gun, and the two find each other by themselves over the ESP32's ESP-NOW radio. Data flows both ways, force feedback and configuration included, so the computer sees exactly the same gun as with the USB cable: same devices, same 209 Hz update rate, no drivers.

## Supported Hardware and Wiring

There are two ways to get a dongle: a ready-made USB stick (the quickest, recommended choice) or a standard development board wired by you.

### 1. "Plug & Play" Solutions (No soldering)

These sticks have a USB connector on the board itself: install the firmware and plug them in. No soldering, no wiring.

* **LILYGO T-Dongle-S3:** an excellent choice. Compact, in a plastic case, with a small colour display that shows the connection status and the project logo.
* **ESP32-S3 Pocket Dongle S3:** another excellent ready-made option, also with a display, ready to plug into the computer or a USB hub.

<table width="100%">
  <tr>
    <th width="33%" align="center">LILYGO T-Dongle-S3</th>
    <th width="33%" align="center">ESP32-S3 Pocket Dongle S3</th>
    <th width="33%" align="center">Installation Example</th>
  </tr>
  <tr>
    <td align="center"><img src="docs/board_scheme/LILYGO-T-Dongle-S3-ESP32-S3.svg" width="100%" alt="LILYGO T-Dongle-S3"></td>
    <td align="center"><img src="docs/board_scheme/esp32-s3-pocket-dongle-s3.svg" width="100%" alt="Pocket Dongle S3"></td>
    <td align="center"><img src="docs/img/dongle_on_PC.jpg" width="100%" alt="Dongle on PC example"></td>
  </tr>
</table>

### 2. "Do-It-Yourself" Solution (Wiring on standard development boards)
To use a standard board, or to fit the receiver in a case of your own, you can use the **Waveshare ESP32-S3-PICO**, the very small **Waveshare ESP32-S3-ZERO** or the classic **ESP32-S3-DevKitC-1**.

These boards connect to the computer with a USB cable, the same one used to program them, which carries both power and data (USB OTG port).

**Antenna:** keep a small clear area around the board's antenna (the printed antenna at the end of the ESP32-S3 module). Do not route wires over it or right next to it, and do not cover it with metal: wires touching the antenna greatly reduce the wireless range and reliability.

| Optional Component | Image |
| :--- | :---: |
| To show information, you can add a **160x80 IPS colour LCD - 0.96 inches - ST7735**. | <img src="docs/img/IPS_0_96_pollici_TFT_display_LCD_ST7735_80x160.png" width="160" alt="ST7735 0.96 inch Display"> |

> [!NOTE]
> **The display is optional but highly recommended: it shows the assigned player, the name of the connected gun and the connection status at a glance.**

Pinouts of the supported boards for use as a dongle:

| Waveshare ESP32-S3-PICO | Waveshare ESP32-S3-ZERO | ESP32-S3-DevKitC-1 |
| :---: | :---: | :---: |
| <img src="docs/board_scheme/DONGLE-esp32-s3-pico.svg" width="100%" alt="PICO Dongle Pinout"> | <img src="docs/board_scheme/DONGLE-esp32-s3-zero.svg" width="100%" alt="ZERO Dongle Pinout"> | <img src="docs/board_scheme/DONGLE-ESP32S3-Devkit-C.svg" width="100%" alt="DevKitC Dongle Pinout"> |

<br>

<p align="center">
  <b>Wiring example on a breadboard (Waveshare ESP32-S3-PICO):</b><br><br>
  <img src="docs/img/dongle_PICO_bb.png" width="80%" alt="Breadboard wiring scheme for PICO Dongle">
</p>

---

## Firmware Installation and Flashing

In the Web Flasher choose **Dongle**, your exact board and its flash/PSRAM variant. Use firmware from the same release as the gun (7.0.0). The dongle is programmed through its own USB port; the gun's Trigger + A and B startup shortcuts do not apply to it, so if it does not enter update mode by itself, use its BOOT/RESET buttons.

#### WEB FLASHER (Recommended for all users)
The easiest, fastest and safest way to install or update the firmware. Nothing to install: it runs in your browser.
* **Requirements:** a computer with Google Chrome, Microsoft Edge or Opera.

**[LAUNCH OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?lang=en)**

---

#### *FOR ADVANCED USERS:*
If your browser does not support the Web Flasher, or you prefer the command line, you can install the firmware files yourself, exactly as for the gun.

### Simplified Procedure with Script
The fastest way, with nothing else to install on the computer.

1. Go to the **[Releases](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases)** page.
2. Download the "Simplified Procedure" ZIP for your receiver and operating system (e.g. `OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-windows-64bit.zip`).
   > *Make sure you choose the file of your exact board.*
3. Extract the whole ZIP into a folder on your computer.
4. Plug the dongle into the computer.
5. Run the `flash_firmware` script (`.bat` on Windows, `.sh` on Linux and macOS) and follow the instructions on screen.

### Troubleshooting
* **The installation does not start (`Connecting...`):** some boards and sticks do not enter update mode by themselves. If the script keeps repeating `Connecting...`, hold the small **BOOT** (or `B`) button on the device until the installation starts.
* **Antivirus false positive (Windows):** the script uses Espressif's original `esptool.exe`, which some antivirus programs may block or flag. Download the package from the official Releases page and check its origin before allowing it to run; do not disable your antivirus globally.
* **Manual installation:** the single `.bin` file of each board is also on the Releases page, for esptool or graphical tools such as NodeMCU PyFlasher (address `0x0`).

---

## Boot and Synchronization Sequence (Pairing)

**Configuring the gun through the dongle:** once the gun is paired, open the [WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/?lang=en) in Chrome or Edge on a computer and choose the dongle's serial port: you configure the gun wirelessly, as if it were connected by cable. The firmware, however, is always installed through each device's own USB port. Close other programs using the same port first.

The dongle needs no buttons and no setup: it does everything by itself. When you plug it into the computer:

1. **Channel selection:** the dongle checks the radio channels around it and picks the one with the least interference (about 15 seconds; with a display you can follow the phases on screen).
2. **Waiting:** it then waits for a gun *(with a display you see a search or wait icon)*.
3. **Pairing:** when you switch on the gun without its USB cable, it announces itself; the first free dongle answers and the two pair exclusively.
4. **Play:** from now on the dongle turns the gun's data into standard mouse, keyboard and gamepad inputs, and carries the serial port too. The computer handles the gun exactly like a wired controller.
5. **Fast reconnection:** if you switch the gun off while the dongle stays plugged in, it reconnects almost instantly when you switch it back on.

> [!TIP]
> **Recommended order:** 1) plug in the dongle and wait about 15 seconds, until it is waiting for a gun; 2) switch on the wireless pedal, if you use one; 3) switch on the gun, without its USB cable connected to the computer. With more guns (up to four), use one dongle per gun: see [several guns](../lightgun/src/README.md#multiple-guns-and-multiplayer).

> [!IMPORTANT]
> **Dongle restart and new pairing**
>
> Every time the dongle is unplugged and plugged in again, it starts over from step 1, ready to pair with any gun (the same one too, with a new pairing).
> *So after restarting the dongle, **switch the gun off and on again** as well, so that it looks for a dongle again.*

---
### Questions or Issues?
For technical support and to join the discussion, see the [Community & Support section](../README.md#community-support-english) on the project's home page.

---

<a id="versione-italiana"></a>

[Torna alla Home](../README.md#versione-italiana) / **Dongle Firmware**

<p align="center">
  <a href="#english-version"><img src="../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# Dongle Firmware (ESP32-S3)

<p align="center">
  <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases"><img src="https://img.shields.io/github/v/release/alessandro-satanassi/OpenFIRE-Firmware-ESP32?include_prereleases&style=flat-square&color=007ec6&label=ultima%20versione" alt="Ultima Versione"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/dongle"><img src="https://img.shields.io/github/languages/top/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&color=success" alt="Linguaggio Principale"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/dongle"><img src="https://img.shields.io/badge/Platform-PlatformIO-orange?style=flat-square&logo=platformio" alt="PlatformIO"></a> <a href="../README.md#community-support-italiano"><img src="https://img.shields.io/badge/Discord-Community-5865F2?style=flat-square&logo=discord&logoColor=white" alt="Discord Community"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/LICENSE"><img src="https://img.shields.io/github/license/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&label=licenza" alt="Licenza"></a>
</p>

<p align="center">
  <img src="docs/img/copertina_dongle.jpg" alt="OpenFIRE Firmware ESP32 dongle" width="100%">
</p>

---
> **Hardware sponsored by [PCBWay](https://www.pcbway.com)**
---

Il dongle è ciò che permette di giocare senza fili. Inseriscilo in una porta USB del computer, accendi la pistola, e i due si trovano da soli tramite la radio ESP-NOW dell'ESP32. I dati viaggiano in entrambe le direzioni, force feedback e configurazione compresi, quindi il computer vede esattamente la stessa pistola che vedrebbe con il cavo USB: stessi dispositivi, stesso aggiornamento a 209 Hz, nessun driver.

## Hardware Supportato e Cablaggio

Ci sono due modi per avere un dongle: una chiavetta USB già pronta (la scelta più rapida e consigliata) oppure una scheda di sviluppo standard cablata da te.

### 1. Soluzioni "Plug & Play" (Senza saldature)

Queste chiavette hanno il connettore USB direttamente sulla scheda: installi il firmware e le inserisci. Niente saldature, niente cablaggi.

* **LILYGO T-Dongle-S3:** un'ottima scelta. Compatta, con involucro in plastica e un piccolo display a colori che mostra lo stato della connessione e il logo del progetto.
* **ESP32-S3 Pocket Dongle S3:** un'altra ottima soluzione già pronta, anch'essa con display, da inserire nel computer o in un hub USB.

<table width="100%">
  <tr>
    <th width="33%" align="center">LILYGO T-Dongle-S3</th>
    <th width="33%" align="center">ESP32-S3 Pocket Dongle S3</th>
    <th width="33%" align="center">Esempio di Installazione</th>
  </tr>
  <tr>
    <td align="center"><img src="docs/board_scheme/LILYGO-T-Dongle-S3-ESP32-S3.svg" width="100%" alt="LILYGO T-Dongle-S3"></td>
    <td align="center"><img src="docs/board_scheme/esp32-s3-pocket-dongle-s3.svg" width="100%" alt="Pocket Dongle S3"></td>
    <td align="center"><img src="docs/img/dongle_on_PC.jpg" width="100%" alt="Esempio dongle su PC"></td>
  </tr>
</table>

### 2. Soluzione "Fai-Da-Te" (Cablaggio su schede di sviluppo standard)
Per usare una scheda standard, o per inserire il ricevitore in un involucro tuo, puoi usare la **Waveshare ESP32-S3-PICO**, la piccolissima **Waveshare ESP32-S3-ZERO** o la classica **ESP32-S3-DevKitC-1**.

Queste schede si collegano al computer con un cavo USB, lo stesso usato per programmarle, che porta sia l'alimentazione sia i dati (porta USB OTG).

**Antenna:** lascia libera una piccola area intorno all'antenna della scheda (l'antenna stampata all'estremità del modulo ESP32-S3). Non far passare fili sopra o a ridosso e non coprirla con metallo: fili che toccano l'antenna riducono molto la portata e l'affidabilità del collegamento wireless.

| Componente Opzionale | Immagine |
| :--- | :---: |
| Per mostrare le informazioni puoi aggiungere un **display LCD IPS a colori 160x80 - 0,96 pollici - ST7735**. | <img src="docs/img/IPS_0_96_pollici_TFT_display_LCD_ST7735_80x160.png" width="160" alt="Display ST7735 0,96 pollici"> |

> [!NOTE]
> **Il display è facoltativo ma vivamente consigliato: mostra a colpo d'occhio il giocatore assegnato, il nome della pistola collegata e lo stato della connessione.**

Pinout delle schede supportate per l'uso come dongle:

| Waveshare ESP32-S3-PICO | Waveshare ESP32-S3-ZERO | ESP32-S3-DevKitC-1 |
| :---: | :---: | :---: |
| <img src="docs/board_scheme/DONGLE-esp32-s3-pico.svg" width="100%" alt="Pinout Dongle PICO"> | <img src="docs/board_scheme/DONGLE-esp32-s3-zero.svg" width="100%" alt="Pinout Dongle ZERO"> | <img src="docs/board_scheme/DONGLE-ESP32S3-Devkit-C.svg" width="100%" alt="Pinout Dongle DevKitC"> |

<br>

<p align="center">
  <b>Esempio di cablaggio su breadboard (Waveshare ESP32-S3-PICO):</b><br><br>
  <img src="docs/img/dongle_PICO_bb.png" width="80%" alt="Schema di cablaggio su breadboard per Dongle PICO">
</p>

---

## Installazione e Flashing del Firmware

Nel Web Flasher scegli **Dongle**, la tua scheda esatta e la sua variante flash/PSRAM. Usa firmware della stessa release della pistola (7.0.0). Il dongle si programma tramite la sua porta USB; le scorciatoie all'avvio della pistola (Grilletto + A e B) non valgono per il dongle, quindi se non entra da solo in modalità aggiornamento usa i suoi pulsanti BOOT/RESET.

#### WEB FLASHER (Consigliato per qualsiasi utente)
Il modo più semplice, veloce e sicuro per installare o aggiornare il firmware. Niente da installare: funziona nel browser.
* **Requisiti:** un computer con Google Chrome, Microsoft Edge o Opera.

**[AVVIA OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?lang=it)**

---

#### *PER UTENTI ESPERTI:*
Se il tuo browser non supporta il Web Flasher, o preferisci la riga di comando, puoi installare manualmente i file del firmware, esattamente come per la pistola.

### Procedura Semplificata con Script
Il metodo più veloce, senza nient'altro da installare sul computer.

1. Vai alla pagina delle **[Releases](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases)**.
2. Scarica lo ZIP "Procedura Semplificata" per il tuo ricevitore e il tuo sistema operativo (es. `OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-windows-64bit.zip`).
   > *Assicurati di scegliere il file della tua scheda esatta.*
3. Estrai tutto lo ZIP in una cartella del computer.
4. Inserisci il dongle nel computer.
5. Esegui lo script `flash_firmware` (`.bat` su Windows, `.sh` su Linux e macOS) e segui le istruzioni a schermo.

### Risoluzione dei Problemi
* **L'installazione non parte (`Connecting...`):** alcune schede e chiavette non entrano da sole in modalità aggiornamento. Se lo script continua a ripetere `Connecting...`, tieni premuto il piccolo pulsante **BOOT** (o `B`) del dispositivo finché l'installazione non parte.
* **Falso positivo dell'antivirus (Windows):** lo script usa l'`esptool.exe` originale di Espressif, che alcuni antivirus possono bloccare o segnalare. Scarica il pacchetto dalla pagina Releases ufficiale e verificane la provenienza prima di consentirne l'esecuzione; non disattivare l'antivirus globalmente.
* **Installazione manuale:** nella pagina Releases c'è anche il singolo file `.bin` di ogni scheda, per esptool o programmi grafici come NodeMCU PyFlasher (indirizzo `0x0`).

---

## Sequenza di Avvio e Sincronizzazione (Pairing)

**Configurare la pistola tramite il dongle:** con la pistola associata, apri la [WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/?lang=it) con Chrome o Edge su computer e scegli la porta seriale del dongle: configuri la pistola senza fili, come se fosse collegata via cavo. Il firmware invece si installa sempre dalla porta USB di ciascun dispositivo. Chiudi prima gli altri programmi che usano la stessa porta.

Il dongle non richiede pulsanti né impostazioni: fa tutto da solo. Quando lo inserisci nel computer:

1. **Scelta del canale:** il dongle controlla i canali radio intorno a sé e sceglie quello con meno interferenze (circa 15 secondi; con un display puoi seguire le fasi a schermo).
2. **Attesa:** poi resta in attesa di una pistola *(con un display vedi un'icona di ricerca o di attesa)*.
3. **Associazione:** quando accendi la pistola senza il cavo USB, questa si annuncia; il primo dongle libero risponde e i due si associano in modo esclusivo.
4. **Gioco:** da questo momento il dongle trasforma i dati della pistola in normali comandi di mouse, tastiera e gamepad, e trasporta anche la porta seriale. Il computer gestisce la pistola esattamente come un controller cablato.
5. **Riconnessione rapida:** se spegni la pistola lasciando il dongle inserito, alla riaccensione si ricollega quasi all'istante.

> [!TIP]
> **Ordine di accensione consigliato:** 1) inserisci il dongle e attendi circa 15 secondi, finché è in attesa di una pistola; 2) accendi il pedale wireless, se lo usi; 3) accendi la pistola, senza il cavo USB collegato al computer. Con più pistole (fino a quattro) usa un dongle per ogni pistola: vedi [più pistole](../lightgun/src/README.md#modifica-dell-id-usb-per-pistole-multiple-italiano).

> [!IMPORTANT]
> **Riavvio del dongle e nuova associazione**
>
> Ogni volta che il dongle viene scollegato e ricollegato, riparte dal passo 1, pronto ad associarsi a qualsiasi pistola (anche la stessa, con una nuova associazione).
> *Quindi, dopo aver riavviato il dongle, **spegni e riaccendi anche la pistola**, perché torni a cercare un dongle.*

---
### Domande o Problemi?
Per supporto tecnico e per unirti alla community, consulta la [sezione Community e Supporto](../README.md#community-support-italiano) nella Home del progetto.
