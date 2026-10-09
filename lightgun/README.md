<a id="english-version"></a>

[Back to Home](../README.md#english-version) / **Lightgun Firmware**

<p align="center">
  <a href="#english-version"><img src="../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# Lightgun Firmware (ESP32-S3)

<p align="center">
  <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases"><img src="https://img.shields.io/github/v/release/alessandro-satanassi/OpenFIRE-Firmware-ESP32?include_prereleases&style=flat-square&color=007ec6&label=latest%20version" alt="Latest Version"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/lightgun"><img src="https://img.shields.io/github/languages/top/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&color=success" alt="Top Language"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/lightgun"><img src="https://img.shields.io/badge/Platform-PlatformIO-orange?style=flat-square&logo=platformio" alt="PlatformIO"></a> <a href="../README.md#community-support-english"><img src="https://img.shields.io/badge/Discord-Community-5865F2?style=flat-square&logo=discord&logoColor=white" alt="Discord Community"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/LICENSE"><img src="https://img.shields.io/github/license/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square" alt="License"></a>
</p>

<p align="center">
  <img src="docs/img/cornice_lightgun.jpg" alt="OpenFIRE Firmware ESP32 lightgun" width="100%">
</p>

<p align="center">
  <b><a href="src/README.md#english-version">OPERATIONAL MANUAL & USER INSTRUCTIONS</a></b>
</p>

---
> **Hardware sponsored by [PCBWay](https://www.pcbway.com)**
---

Everything you need to build the lightgun and install its firmware: the parts, how to wire them, which board to choose and what happens when you switch it on. Once it is built, the [operational manual](src/README.md#english-version) takes you through configuration, calibration and play.

---

## Hardware Requirements

A working gun needs very little: a board, a camera, four IR emitters and a trigger. From there you can add buttons, recoil, rumble, lights and a display, up to a complete arcade gun.

The picture below shows the reference build (from the PICON-AS project), with the essential parts and the optional modules supported by the firmware.

<p align="center">
  <img src="docs/img/picon_AS_opzioni_refactory.png" alt="Lightgun Hardware Architecture Reference" width="100%">
</p>

### Mandatory Components (Essentials)
To aim and shoot you need:

* **Microcontroller:** an **ESP32-S3** board. Recommended and tested boards:
  * *ESP32-S3-WROOM1-DevKitC-1*
  * *Waveshare ESP32-S3-PICO*
  * *Waveshare ESP32-S3-ZERO (Mini)*
  * **Antenna:** keep a small clear area around the board's antenna (the printed antenna at the end of the ESP32-S3 module). Do not route wires over it or right next to it, and do not cover it with metal: wires touching the antenna greatly reduce the wireless range and reliability.
* **Optical Sensor:** a **DFRobot SEN0158 / Wii Cam**, a **PixArt PAJ7025R2** or a **PixArt PAJ7025R3** (wide-angle lens). The same firmware works with all three: select yours in **Gun Settings → CAMERA Model**, assign its pins in **Board Layout**, then save. Calibrate again after changing the camera or lens.
  * DFRobot / Wii Cam connects over **I2C** (**Camera SDA**, **Camera SCL**). A bare Wii camera, without the DFRobot board, also needs a clock signal: the gun generates it on the pin you assign to **Wii Cam Clock**.
  * PAJ7025R2 / R3 connect over **SPI** (**Camera SPI MISO**, **MOSI**, **SCK** and **CS**), besides power and ground; follow the wiring and supply of your camera module.
  * **Which one to choose?** The DFRobot SEN0158 is no longer made; a Wii camera can take its place, with the same **DFRobot SEN0158 / Wii** setting (the default after a clean installation). The PixArt PAJ7025R2 and R3 are the newer, high-performance cameras: they need SPI wiring and 850 nm emitters, and they also measure the brightness of each IR point, shown in the WebApp IR test. The R3's wide-angle lens suits playing close to the screen, though by itself it does not give more accuracy or range. Purchase links and wiring for the supported cameras are on the [PICON-AS site](https://alessandro-satanassi.github.io/OpenFIRE-PICON-AS-ESP32/) (currently in Italian).
* **IR Emitters:** 4 infrared LEDs, **940 nm for DFRobot / Wii Cam** and **850 nm for PixArt PAJ7025R2 / R3**. Wii sensor bars suit the DFRobot/Wii camera only. Give your LEDs a suitable supply and current limiting. Example: the OSRAM SFH 4547 (940 nm) suits DFRobot/Wii, not the PixArt cameras.
* **Primary Input:** at least 1 switch, used as the trigger.

### Highly Recommended Components
The gun works with the trigger alone, but most lightgun games also need:
* **A and B Buttons:** 2 more switches, for reloading and secondary actions.

### Optional Modules:
Add whatever your build needs:

**Switches (Buttons):**
Besides the trigger, the firmware handles up to 13 buttons, all assignable in the OpenFIRE ESP32 WebApp (see the pinout pictures below for the defaults):
* **C Button:** one more button alongside A and B.
* **D-Pad Module:** a 5-way navigation switch (Up, Down, Left, Right, and Centre, used as the Home button that enters pause mode).
* **System Buttons:** 2 buttons for *Start* and *Select*.
* **Pump Switch:** a switch for pump-action reloading.
* **Pedal Button:** a switch for the main pedal *(do not map it if you use the wireless pedal)*.
* **Alt Pedal Button:** a switch for a second pedal, for games that use one *(do not map it if you use the wireless pedal)*.

**Additional Controls:**
* **Analog Joystick:** a 2-axis analog joystick module.

**Force Feedback, Haptics, and Hardware Switches:**
* **Solenoid (12V-24V):** any solenoid with its MOSFET driver board *(needs a separate, adjustable 12-24V power supply)*.
  * *For fully wireless builds:* with 12V solenoids only, you can use a high-discharge (at least 10A) 3.7V Li-ion battery with a powerful booster board such as the **MP3429** (3.7V -> 12V). Both 18650 and 21700 cells work, but a **21700** is strongly recommended for a reasonable battery life.
  * *Beware of interference (EMI):* solenoids can cause USB or radio disconnections if the wiring is too thin. Use solenoid power wires of at least **24 AWG**.
* **Temperature Sensor (TMP36):** placed against the solenoid, it protects it from overheating during long sessions.
  * *Wiring tip:* give the sensor its own ground (GND) wire straight to a GND pin of the board, not shared with the camera or other modules. Current flowing through a shared ground wire makes the sensor read too high (for example about 90 °C instead of 20 °C). A 0.1 µF capacitor between its VCC and GND pins, close to the sensor, makes the reading steadier.
* **Rumble Motor (5V):** any gamepad vibration motor with a basic driver board.
* **2-way SPDT Switches:** to turn Rumble, Solenoid or Rapid Fire on and off with a physical switch. Wire each switch between its pin and GND, then assign the pin (**Rumble Switch**, **Solenoid Switch** or **Autofire Switch**) in **Board Layout**: closed is on, open is off, and the switch has priority over the pause menu. Without switches, rumble and solenoid are controlled from the pause menu.

**Lighting and Display:**
* **NeoPixel:** WS2812B NeoPixel modules for lighting that reacts to the game in real time.
* **RGB LED:** standard 4-pin RGB LEDs, for the same purpose.
* **OLED Display:** a small 128x64 **SSD1306** I2C screen (4 pins) for the menus, the connection status and the life/ammo counters, and to guide you through calibration. It is disabled by default: enable it in the WebApp (**Gun Settings → I2C Peripherals**), as explained in the [manual](src/README.md#camera-display-and-startup-settings).

### Resources and Assembly Guides
Purchase links for the components and step-by-step assembly instructions are in my hardware project **[PICON-AS](https://alessandro-satanassi.github.io/OpenFIRE-PICON-AS-ESP32/)**.

The OpenFIRE project's illustrated guide to wiring the components is also included. *The pictures show a Raspberry Pi Pico (RP2040), but the wiring is the same on any ESP32-S3 pin you assign.*
* **[OpenFIRE Hardware Guide](docs/guide/OpenFIRE-Hardware-Guide_1.2.pdf)**

---

## Default Board Pinouts

The pictures below show the default pin assignment of each board.
> **Note:** the functions of the available pins can be reassigned. If it is easier to solder a component to another suitable pin, change it in the **Board Layout** tab of the **OpenFIRE ESP32 WebApp**; leave pins reserved for USB, startup or other board functions alone.

| ESP32-S3-DevKitC-1 | Waveshare ESP32-S3-PICO | Waveshare ESP32-S3-ZERO |
| :---: | :---: | :---: |
| <img src="docs/board_scheme/LIGHTGUN-ESP32S3-Devkit-C.svg" width="100%" alt="Pinout DevKitC-1"> | <img src="docs/board_scheme/LIGHTGUN-esp32-s3-pico.svg" width="100%" alt="Pinout S3-PICO"> | <img src="docs/board_scheme/LIGHTGUN-esp32-s3-zero.svg" width="100%" alt="Pinout S3-ZERO"> |

---

## Firmware Installation and Flashing

#### WEB FLASHER (Recommended for all users)
The easiest, fastest and safest way to install or update the firmware. Nothing to install: it runs in your browser.
* **Requirements:** a computer with Google Chrome, Microsoft Edge or Opera.

**[LAUNCH OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?lang=en)**

Choose **Lightgun**, your board and its exact variant, then **Base Update** to keep your settings or **Clean Install** to erase everything (recommended for a new board or when coming from 6.2.1).

<a id="board-variant"></a>

> **Which board variant do I have?** Choose the exact board and variant, for example **DevKitC-1 N16R8** or **N8R2**, **Waveshare ZERO N8R8** or **N4R2**. The code is printed on the metal shield of the module or on the chip (for example *ESP32-S3-WROOM-1-N16R8*, or *ESP32-S3FH4R2* on a Waveshare ZERO) and is given in the seller's description: the number after **N** (or **H**) is the flash size in MB, the one after **R** the PSRAM size in MB. If the board does not start, or keeps restarting after the installation, install again choosing the correct variant, using the board's BOOT/RESET buttons if needed.

---

#### *FOR ADVANCED USERS:*
If your browser does not support the Web Flasher, or you prefer the command line, you can install the firmware files yourself. The ESP32-S3 is programmed with binary (`.bin`) files over its USB port; ready-made packages for every supported board are on the **[Releases](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases)** page.

### Method 1: Simplified Script Procedure
The fastest way, with nothing else to install. Packages are available for **Windows, Linux and macOS**.

1. On the **[Releases](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases)** page, download the "Simplified Procedure" ZIP for your board and operating system (e.g. `OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip`).
2. Extract the whole ZIP into a folder on your computer.
3. Connect the board with a USB cable *(a data cable, not a charge-only one)*.
4. Run the `flash_firmware` script (`.bat` on Windows, `.sh` on Linux and macOS).
5. The script finds the serial port by itself and guides you through the installation.

> **Normal update or clean installation**
>
> There is **one file per board**, whatever the camera.
> * **Normal update:** writes the firmware and keeps the saved settings and calibrations. Press **Enter** in the script.
> * **Clean installation:** erases the whole flash, then writes the **same firmware**. Type **1** in the script. **All settings and calibrations are deleted**; defaults are created at the first start.
>
> Coming from 6.2.1, a clean installation is recommended. Note your settings first, then configure and calibrate again.

### Method 2: Manual Installation
With **esptool** or a graphical tool such as **NodeMCU PyFlasher**, write the single `.bin` file from the Releases page at address `0x0`. For a clean installation with esptool 5.x add `--erase-all` to `write-flash`; leave it out for a normal update. Installation always goes through the board's USB port, not over Wi-Fi or the wireless dongle.

---

### Troubleshooting

* **Entering update mode:** on a gun running 7.0.0, hold **Trigger + A while switching on, for about 2 seconds**; the OLED, if fitted, shows **Ready for firmware update**. Connect the gun's own USB OTG port, not the dongle. If the gun is running normally, the Web Flasher can also restart it into update mode by itself: the port disappears and comes back, so wait, start again and select the new port. For a first installation, older firmware or a recovery, use the board's **BOOT/RESET** buttons. The board's BOOT button is not the gun's B button.
* **Antivirus false positive (Windows):** the script uses Espressif's original `esptool.exe`, which some antivirus programs may block or flag. Download the package from the official Releases page and check its origin before allowing it to run; do not disable your antivirus globally.
* **After the installation:** open the **[OpenFIRE ESP32 WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/?lang=en)**, or use the offline WebApp described below.

<a id="manual-update-mode"></a>

**BOOT/RESET, for a first installation or recovery:** on an ESP32-S3 board fitted with both buttons, hold **BOOT**, press and release **RESET** (sometimes labelled **EN**), then release **BOOT**. Choose the port that appears and install the firmware. This is the [manual download-mode procedure](https://docs.espressif.com/projects/esptool/en/latest/esp32s3/advanced-topics/boot-mode-selection.html#manual-bootloader) described by Espressif; if your board has different buttons, follow its own instructions. After installation, reset or switch the board off and on again without holding BOOT.

## Configuration with the WebApp

* **Online:** open the [WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/?lang=en) in Chrome or Edge on a computer, then select the lightgun or its paired dongle.
* **Offline:** switch on the gun holding **B** for about 2 seconds. Join **OpenFIRE_Config**, accept the network without Internet, then open **http://openfire.local** or **http://192.168.4.1** in your normal browser. With a USB data cable you can instead open **http://192.168.7.1**, if the computer supports USB networking (NCM).
* **Save and restart:** save, close the WebApp, then restart the gun without holding any buttons to return to normal play. While the offline WebApp is running, serial programs such as MAMEHOOKER cannot use the gun's USB port.

The [operational manual](src/README.md#configuration-with-the-webapp) explains every step, the camera choice, the startup buttons and the limits of each connection.

## Boot Sequence and Connectivity

Here is what the gun does at a normal start, without any button combination:

### 1. Automatic Joystick Calibration
If an analog joystick is fitted, the gun reads its centre position as soon as it is powered.
> [!IMPORTANT]
> This takes about two seconds: **do not touch the joystick** meanwhile.

### 2. Connectivity Management (USB vs Wireless)
The gun chooses how to talk to the computer by itself:
* **USB cable:** if the cable is connected to the computer at startup, the gun plays wired and switches the radio off to save power. (The B startup mode turns on the configuration network instead.)
* **Wireless (ESP-NOW):** without the cable, the wireless sequence starts:
  1. The **dongle**, already plugged into the computer, picks the radio channel with the least interference and waits.
  2. The **gun** announces itself on every channel until a free dongle answers. *(If you plug in the USB cable meanwhile, the gun stops searching and plays wired.)*
  3. The two **pair**. From now on the computer sees the gun exactly as if it were connected by cable, with the same performance.

### 3. Wireless Pedal Search (Optional)
Right after pairing with the dongle (with **Wireless pedal** enabled in Gun Settings and no wired pedal pin assigned), the gun looks for a free wireless pedal for **10 seconds**.
* If it finds one, they pair, and the pedal LEDs show the player number.
* If you do not use a wireless pedal, disable **Wireless pedal** and save, to skip the search.
* The pedal is searched only when the gun plays through the dongle; with the gun connected by USB cable, use a wired pedal.
* If no pedal answers within 10 seconds, the start completes normally and you are ready to play.

### 4. Fast Reconnection
If only the gun is switched off and on again (for example when the battery runs low), it first looks for the dongle and pedal it was paired with, and reconnects almost instantly.
* *Note:* this works only if the dongle and the pedal stayed on in the meantime. If they were restarted, they wait for a new pairing and the gun performs a full search.

### 5. Connection Status (Visual Feedback)
With an OLED display, an icon at the top shows the connection:
* **Wi-Fi icon:** connected through the wireless dongle.
* **USB icon:** connected by cable.

---

## Operational Manual and User Instructions

> [!IMPORTANT]
> **Gun built and firmware installed? Great job, you are almost there!**
>
> The gun has a pause menu, button combinations and a calibration procedure worth knowing before your first game.

The **Operational Manual** explains:

* **Calibration:** how to calibrate so the pointer lands exactly where you aim.
* **Controls and shortcuts:** the default buttons and how to enter pause mode.
* **Profiles:** how to switch profiles and save your settings in the gun.
* **Adjustments on the fly:** switching rumble and solenoid on and off, and adjusting the IR camera sensitivity, straight from the gun.

**[CLICK HERE TO READ THE FULL OPERATIONAL MANUAL](src/README.md#english-version)**

---

## Technical Information for Developers

This guide is for people who build and use the gun. To modify the firmware or contribute to the project, see the [build guide](docs/COMPILING.md): the firmware is built with **PlatformIO**, using the configuration included in this repository.

---
### Questions or Issues?
For technical support and to join the discussion, see the [Community & Support section](../README.md#community-support-english) on the project's home page.

---

<a id="versione-italiana"></a>

[Torna alla Home](../README.md#versione-italiana) / **Lightgun Firmware**

<p align="center">
  <a href="#english-version"><img src="../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# Lightgun Firmware (ESP32-S3)

<p align="center">
  <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases"><img src="https://img.shields.io/github/v/release/alessandro-satanassi/OpenFIRE-Firmware-ESP32?include_prereleases&style=flat-square&color=007ec6&label=ultima%20versione" alt="Ultima Versione"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/lightgun"><img src="https://img.shields.io/github/languages/top/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&color=success" alt="Linguaggio Principale"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/tree/main/lightgun"><img src="https://img.shields.io/badge/Platform-PlatformIO-orange?style=flat-square&logo=platformio" alt="PlatformIO"></a> <a href="../README.md#community-support-italiano"><img src="https://img.shields.io/badge/Discord-Community-5865F2?style=flat-square&logo=discord&logoColor=white" alt="Discord Community"></a> <a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/LICENSE"><img src="https://img.shields.io/github/license/alessandro-satanassi/OpenFIRE-Firmware-ESP32?style=flat-square&label=licenza" alt="Licenza"></a>
</p>

<p align="center">
  <img src="docs/img/cornice_lightgun.jpg" alt="OpenFIRE Firmware ESP32 lightgun" width="100%">
</p>

<p align="center">
  <b><a href="src/README.md#versione-italiana">MANUALE OPERATIVO E ISTRUZIONI PER L'USO</a></b>
</p>

---
> **Hardware sponsored by [PCBWay](https://www.pcbway.com)**
---

Tutto quello che serve per costruire la lightgun e installarne il firmware: i componenti, come collegarli, quale scheda scegliere e cosa succede quando la accendi. Una volta costruita, il [manuale operativo](src/README.md#versione-italiana) ti guida in configurazione, calibrazione e gioco.

---

## Requisiti Hardware

Per una pistola funzionante serve pochissimo: una scheda, una telecamera, quattro emettitori IR e un grilletto. Da lì puoi aggiungere pulsanti, rinculo, vibrazione, luci e display, fino a una pistola arcade completa.

L'immagine qui sotto mostra la costruzione di riferimento (dal progetto PICON-AS), con i componenti essenziali e i moduli opzionali supportati dal firmware.

<p align="center">
  <img src="docs/img/picon_AS_opzioni_refactory.png" alt="Schema Architettura Hardware Lightgun" width="100%">
</p>

### Componenti Obbligatori (Essenziali)
Per mirare e sparare servono:

* **Microcontrollore:** una scheda **ESP32-S3**. Schede consigliate e testate:
  * *ESP32-S3-WROOM1-DevKitC-1*
  * *Waveshare ESP32-S3-PICO*
  * *Waveshare ESP32-S3-ZERO (Mini)*
  * **Antenna:** lascia libera una piccola area intorno all'antenna della scheda (l'antenna stampata all'estremità del modulo ESP32-S3). Non far passare fili sopra o a ridosso e non coprirla con metallo: fili che toccano l'antenna riducono molto la portata e l'affidabilità del collegamento wireless.
* **Sensore Ottico:** una **DFRobot SEN0158 / Wii Cam**, una **PixArt PAJ7025R2** o una **PixArt PAJ7025R3** (ottica grandangolare). Lo stesso firmware funziona con tutte e tre: seleziona la tua in **Impostazioni Gun → Modello TELECAMERA**, assegna i suoi pin nel **Layout Scheda**, poi salva. Ripeti la calibrazione dopo aver cambiato telecamera o lente.
  * DFRobot / Wii Cam si collega in **I2C** (**SDA telecamera**, **SCL telecamera**). Una telecamera Wii senza la scheda DFRobot richiede anche un segnale di clock: la pistola lo genera sul pin a cui assegni **Clock Wii Cam**.
  * PAJ7025R2 / R3 si collegano in **SPI** (**SPI MISO**, **MOSI**, **SCK** e **CS telecamera**), oltre ad alimentazione e massa; rispetta cablaggio e alimentazione del tuo modulo telecamera.
  * **Quale scegliere?** La DFRobot SEN0158 non è più in produzione; al suo posto si può usare una telecamera Wii, con la stessa impostazione **DFRobot SEN0158 / Wii** (quella predefinita dopo un'installazione pulita). Le PixArt PAJ7025R2 e R3 sono le telecamere più recenti e prestanti: richiedono il collegamento SPI ed emettitori da 850 nm, e misurano anche la luminosità di ogni punto IR, mostrata nel test IR della WebApp. L'ottica grandangolare della R3 è adatta a chi gioca vicino allo schermo, anche se da sola non dà più precisione o portata. Link per l'acquisto e collegamenti delle telecamere supportate sono sul [sito PICON-AS](https://alessandro-satanassi.github.io/OpenFIRE-PICON-AS-ESP32/).
* **Emettitori IR:** 4 LED infrarossi, **da 940 nm per DFRobot / Wii Cam** e **da 850 nm per PixArt PAJ7025R2 / R3**. Le barre sensore Wii vanno bene solo con la telecamera DFRobot/Wii. Dai ai LED un'alimentazione e una limitazione di corrente adatte. Esempio: l'OSRAM SFH 4547 (940 nm) è adatto a DFRobot/Wii, non alle telecamere PixArt.
* **Input Primario:** almeno 1 interruttore, usato come grilletto.

### Componenti Altamente Consigliati
La pistola funziona anche con il solo grilletto, ma la maggior parte dei giochi per lightgun richiede anche:
* **Pulsanti A e B:** 2 interruttori in più, per ricaricare e per le azioni secondarie.

### Moduli Opzionali:
Aggiungi quello che serve alla tua costruzione:

**Interruttori (Pulsanti):**
Oltre al grilletto il firmware gestisce fino a 13 pulsanti, tutti assegnabili nella WebApp OpenFIRE ESP32 (vedi le immagini dei pinout più sotto per quelli predefiniti):
* **Pulsante C:** un pulsante in più oltre ad A e B.
* **Modulo D-Pad:** un interruttore di navigazione a 5 vie (Su, Giù, Sinistra, Destra, e Centro, usato come pulsante Home che entra in modalità pausa).
* **Pulsanti di Sistema:** 2 pulsanti per *Start* e *Select*.
* **Switch a Pompa:** un interruttore per la ricarica a pompa.
* **Pulsante Pedale:** un interruttore per il pedale principale *(non assegnarlo se usi il pedale wireless)*.
* **Pulsante Pedale alternativo:** un interruttore per un secondo pedale, per i giochi che lo usano *(non assegnarlo se usi il pedale wireless)*.

**Controlli Aggiuntivi:**
* **Joystick Analogico:** un modulo joystick analogico a 2 assi.

**Force feedback, feedback aptico e interruttori fisici:**
* **Solenoide (12V-24V):** qualsiasi solenoide con la sua scheda driver MOSFET *(richiede un alimentatore separato regolabile da 12-24V)*.
  * *Per costruzioni completamente senza fili:* solo con solenoidi a 12 V puoi usare una batteria Li-ion da 3,7 V ad alta scarica (almeno 10 A) con una scheda booster potente come l'**MP3429** (3,7 V -> 12 V). Vanno bene sia le 18650 sia le 21700, ma per una durata ragionevole è vivamente consigliata una **21700**.
  * *Attenzione alle interferenze (EMI):* i solenoidi possono causare disconnessioni USB o radio se il cablaggio è troppo sottile. Usa per il solenoide fili di alimentazione di almeno **24 AWG**.
* **Sensore di Temperatura (TMP36):** posizionato a contatto con il solenoide, lo protegge dal surriscaldamento durante le sessioni lunghe.
  * *Consiglio di cablaggio:* collega il GND del sensore con un filo dedicato direttamente a un pin GND della scheda, non condiviso con la telecamera o altri moduli. La corrente che scorre in un filo di massa in comune fa leggere al sensore valori troppo alti (ad esempio circa 90 °C invece di 20 °C). Un condensatore da 0,1 µF tra i suoi pin VCC e GND, vicino al sensore, rende la lettura più stabile.
* **Motore Rumble (5V):** qualsiasi motorino di vibrazione per gamepad con una semplice scheda driver.
* **Interruttori SPDT a 2 vie:** per accendere e spegnere Rumble, Solenoide o Fuoco Rapido con un interruttore fisico. Collega ogni interruttore tra il suo pin e GND, poi assegna il pin (**Interruttore Rumble**, **Interruttore Solenoide** o **Interruttore Autofire**) in **Layout Scheda**: chiuso è acceso, aperto è spento, e l'interruttore ha la precedenza sul menu di pausa. Senza interruttori, rumble e solenoide si comandano dal menu di pausa.

**Illuminazione e Display:**
* **NeoPixel:** moduli NeoPixel WS2812B per luci che reagiscono al gioco in tempo reale.
* **LED RGB:** classici LED RGB a 4 pin, per lo stesso scopo.
* **Display OLED:** un piccolo schermo I2C **SSD1306** da 128x64 (4 pin) per i menu, lo stato della connessione e i contatori di vite/munizioni, e per guidarti nella calibrazione. È disattivato per impostazione predefinita: attivalo nella WebApp (**Impostazioni Gun → Periferiche I2C**), come spiegato nel [manuale](src/README.md#telecamera-display-e-impostazioni-di-avvio).

### Risorse e Guide all'Assemblaggio
Link per l'acquisto dei componenti e istruzioni di montaggio passo per passo sono nel mio progetto hardware **[PICON-AS](https://alessandro-satanassi.github.io/OpenFIRE-PICON-AS-ESP32/)**.

È inclusa anche la guida illustrata del progetto OpenFIRE per il collegamento dei componenti. *Le immagini mostrano un Raspberry Pi Pico (RP2040), ma i collegamenti sono gli stessi su qualsiasi pin dell'ESP32-S3 che assegni.*
* **[Guida Componenti OpenFIRE](docs/guide/OpenFIRE-Hardware-Guide_1.2.pdf)**

---

## Pinout Predefiniti delle Schede

Le immagini qui sotto mostrano l'assegnazione predefinita dei pin di ogni scheda.
> **Nota:** puoi riassegnare le funzioni dei pin disponibili. Se è più comodo saldare un componente su un altro pin adatto, cambialo nella sezione **Layout Scheda** della **WebApp OpenFIRE ESP32**; lascia liberi i pin riservati a USB, avvio o altre funzioni della scheda.

| ESP32-S3-DevKitC-1 | Waveshare ESP32-S3-PICO | Waveshare ESP32-S3-ZERO |
| :---: | :---: | :---: |
| <img src="docs/board_scheme/LIGHTGUN-ESP32S3-Devkit-C.svg" width="100%" alt="Pinout DevKitC-1"> | <img src="docs/board_scheme/LIGHTGUN-esp32-s3-pico.svg" width="100%" alt="Pinout S3-PICO"> | <img src="docs/board_scheme/LIGHTGUN-esp32-s3-zero.svg" width="100%" alt="Pinout S3-ZERO"> |

---

## Installazione e Flashing del Firmware

#### WEB FLASHER (Consigliato per qualsiasi utente)
Il modo più semplice, veloce e sicuro per installare o aggiornare il firmware. Niente da installare: funziona nel browser.
* **Requisiti:** un computer con Google Chrome, Microsoft Edge o Opera.

**[AVVIA OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?lang=it)**

Scegli **Lightgun**, la tua scheda e la sua variante esatta, poi **Aggiornamento Base** per conservare le impostazioni oppure **Installazione Pulita** per cancellare tutto (consigliata per una scheda nuova o passando dalla 6.2.1).

<a id="variante-scheda"></a>

> **Quale variante di scheda ho?** Scegli la scheda e la variante esatte, ad esempio **DevKitC-1 N16R8** o **N8R2**, **Waveshare ZERO N8R8** o **N4R2**. La sigla è stampata sulla schermatura metallica del modulo o sul chip (ad esempio *ESP32-S3-WROOM-1-N16R8*, oppure *ESP32-S3FH4R2* su una Waveshare ZERO) ed è indicata nella descrizione del venditore: il numero dopo **N** (o **H**) è la memoria flash in MB, quello dopo **R** la PSRAM in MB. Se dopo l'installazione la scheda non si avvia o si riavvia di continuo, ripeti l'installazione scegliendo la variante corretta, usando se necessario i pulsanti BOOT/RESET della scheda.

---

#### *PER UTENTI ESPERTI:*
Se il tuo browser non supporta il Web Flasher, o preferisci la riga di comando, puoi installare manualmente i file del firmware. L'ESP32-S3 si programma con file binari (`.bin`) attraverso la sua porta USB; nella pagina delle **[Releases](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases)** trovi pacchetti pronti per ogni scheda supportata.

### Metodo 1: Procedura Semplificata con Script
Il metodo più veloce, senza nient'altro da installare. I pacchetti sono disponibili per **Windows, Linux e macOS**.

1. Nella pagina delle **[Releases](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases)** scarica lo ZIP "Procedura Semplificata" per la tua scheda e il tuo sistema operativo (es. `OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip`).
2. Estrai tutto lo ZIP in una cartella del computer.
3. Collega la scheda con un cavo USB *(un cavo dati, non solo di ricarica)*.
4. Esegui lo script `flash_firmware` (`.bat` su Windows, `.sh` su Linux e macOS).
5. Lo script trova da solo la porta seriale e ti guida nell'installazione.

> **Aggiornamento normale oppure installazione pulita**
>
> C'è **un solo file per scheda**, qualunque sia la telecamera.
> * **Aggiornamento normale:** scrive il firmware e conserva impostazioni e calibrazioni salvate. Premi **Invio** nello script.
> * **Installazione pulita:** cancella tutta la flash, poi scrive lo **stesso firmware**. Digita **1** nello script. **Tutte le impostazioni e calibrazioni vengono eliminate**; i valori predefiniti vengono creati al primo avvio.
>
> Passando dalla 6.2.1 è consigliata l'installazione pulita. Annota prima le impostazioni, poi configura e calibra di nuovo.

### Metodo 2: Installazione Manuale
Con **esptool** o un programma grafico come **NodeMCU PyFlasher**, scrivi il singolo file `.bin` della pagina Releases all'indirizzo `0x0`. Per un'installazione pulita con esptool 5.x aggiungi `--erase-all` a `write-flash`; omettilo per un aggiornamento normale. L'installazione passa sempre dalla porta USB della scheda, non dal Wi-Fi o dal dongle wireless.

---

### Risoluzione dei Problemi (Troubleshooting)

* **Ingresso in modalità aggiornamento:** su una pistola con la 7.0.0 tieni premuti **Grilletto + A all'accensione per circa 2 secondi**; l'OLED, se presente, mostra **Ready for firmware update**. Collega la porta USB OTG della pistola, non il dongle. Se la pistola è in funzionamento normale, il Web Flasher può anche riavviarla da solo in modalità aggiornamento: la porta scompare e ricompare, quindi attendi, ricomincia e seleziona la nuova porta. Per la prima installazione, firmware precedenti o un recupero, usa i pulsanti **BOOT/RESET** della scheda. Il BOOT della scheda non è il pulsante B della pistola.
* **Falso positivo dell'antivirus (Windows):** lo script usa l'`esptool.exe` originale di Espressif, che alcuni antivirus possono bloccare o segnalare. Scarica il pacchetto dalla pagina Releases ufficiale e verificane la provenienza prima di consentirne l'esecuzione; non disattivare l'antivirus globalmente.
* **Dopo l'installazione:** apri la **[WebApp OpenFIRE ESP32](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/?lang=it)**, oppure usa la WebApp offline descritta sotto.

<a id="modalita-aggiornamento-manuale"></a>

**BOOT/RESET, per la prima installazione o un recupero:** su una scheda ESP32-S3 con entrambi i pulsanti, tieni premuto **BOOT**, premi e rilascia **RESET** (a volte indicato come **EN**), poi rilascia **BOOT**. Seleziona la porta che compare e installa il firmware. È la [procedura di ingresso manuale in modalità download](https://docs.espressif.com/projects/esptool/en/latest/esp32s3/advanced-topics/boot-mode-selection.html#manual-bootloader) descritta da Espressif; se la tua scheda ha pulsanti diversi, segui le sue istruzioni. Dopo l'installazione riavvia o spegni e riaccendi la scheda senza tenere premuto BOOT.

## Configurazione con la WebApp

* **Online:** apri la [WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/?lang=it) con Chrome o Edge su computer, poi seleziona la lightgun o il dongle associato.
* **Offline:** accendi la pistola tenendo premuto **B** per circa 2 secondi. Collegati a **OpenFIRE_Config**, accetta la rete senza Internet, poi apri **http://openfire.local** oppure **http://192.168.4.1** nel browser normale. Con un cavo dati USB puoi invece aprire **http://192.168.7.1**, se il computer supporta la rete via USB (NCM).
* **Salva e riavvia:** salva, chiudi la WebApp, poi riavvia la pistola senza tenere premuti pulsanti per tornare al gioco normale. Mentre la WebApp offline è attiva, i programmi seriali come MAMEHOOKER non possono usare la porta USB della pistola.

Il [manuale operativo](src/README.md#configurazione-con-la-webapp) spiega ogni passaggio, la scelta della telecamera, i pulsanti all'avvio e i limiti di ogni collegamento.

## Sequenza di Avvio e Connettività

Ecco cosa fa la pistola a un avvio normale, senza combinazioni di pulsanti:

### 1. Calibrazione Automatica del Joystick
Se è montato un joystick analogico, la pistola ne legge la posizione centrale appena riceve alimentazione.
> [!IMPORTANT]
> Ci vogliono circa due secondi: nel frattempo **non toccare il joystick**.

### 2. Gestione della Connettività (USB vs Wireless)
La pistola sceglie da sola come comunicare con il computer:
* **Cavo USB:** se all'avvio il cavo è collegato al computer, la pistola gioca via cavo e spegne la radio per risparmiare energia. (L'avvio con B accende invece la rete di configurazione.)
* **Senza fili (ESP-NOW):** senza il cavo parte la sequenza wireless:
  1. Il **dongle**, già inserito nel computer, sceglie il canale radio con meno interferenze e si mette in attesa.
  2. La **pistola** si annuncia su tutti i canali finché un dongle libero non risponde. *(Se nel frattempo colleghi il cavo USB, la pistola smette di cercare e gioca via cavo.)*
  3. I due si **associano**. Da questo momento il computer vede la pistola esattamente come se fosse collegata via cavo, con le stesse prestazioni.

### 3. Ricerca del Pedale Wireless (Opzionale)
Subito dopo l'associazione con il dongle (con **Pedale wireless** attivo in **Impostazioni Gun** e nessun pin di pedale cablato assegnato), la pistola cerca un pedale wireless libero per **10 secondi**.
* Se lo trova si associano, e i LED del pedale mostrano il numero del giocatore.
* Se non usi un pedale wireless, disattiva **Pedale wireless** e salva, per saltare la ricerca.
* Il pedale viene cercato solo quando la pistola gioca tramite il dongle; con la pistola collegata via cavo USB, usa un pedale cablato.
* Se nessun pedale risponde entro 10 secondi, l'avvio si conclude normalmente e sei pronto a giocare.

### 4. Riconnessione Rapida
Se viene spenta e riaccesa solo la pistola (ad esempio quando la batteria si scarica), cerca per primi il dongle e il pedale a cui era associata e si ricollega quasi all'istante.
* *Nota:* funziona solo se nel frattempo dongle e pedale sono rimasti accesi. Se sono stati riavviati, attendono una nuova associazione e la pistola esegue una ricerca completa.

### 5. Stato della Connessione (Feedback Visivo)
Con un display OLED, un'icona in alto indica il collegamento:
* **Icona Wi-Fi:** collegata tramite il dongle wireless.
* **Icona USB:** collegata via cavo.

---

## Manuale Operativo e Istruzioni d'Uso

> [!IMPORTANT]
> **Pistola costruita e firmware installato? Ottimo lavoro, ci sei quasi!**
>
> La pistola ha un menu di pausa, combinazioni di pulsanti e una procedura di calibrazione che vale la pena conoscere prima di iniziare a giocare.

Il **Manuale Operativo** spiega:

* **Calibrazione:** come calibrare perché il puntatore arrivi esattamente dove miri.
* **Comandi e scorciatoie:** i pulsanti predefiniti e come entrare in modalità pausa.
* **Profili:** come cambiare profilo e salvare le impostazioni nella pistola.
* **Regolazioni al volo:** accendere e spegnere rumble e solenoide e regolare la sensibilità della telecamera IR direttamente dalla pistola.

**[CLICCA QUI PER LEGGERE IL MANUALE OPERATIVO COMPLETO](src/README.md#versione-italiana)**

---

## Informazioni Tecniche per gli Sviluppatori

Questa guida è per chi costruisce e usa la pistola. Per modificare il firmware o contribuire al progetto, consulta la [guida alla compilazione](docs/COMPILING.md): il firmware si compila con **PlatformIO**, usando la configurazione inclusa in questo repository.

---
### Domande o Problemi?
Per supporto tecnico e per unirti alla community, consulta la [sezione Community e Supporto](../README.md#community-support-italiano) nella Home del progetto.
