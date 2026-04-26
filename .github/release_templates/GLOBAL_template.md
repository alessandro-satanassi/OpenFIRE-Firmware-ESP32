# GLOBAL RELEASE - FIRMWARE **LIGHTGUN** **DONGLE** **PEDAL** - OpenFIRE ESP32


This release is updated to commit **8b651a2** (April 19, 2026) of the original [OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware) firmware.

> [!WARNING]
> **NEW! PEDAL FIRMWARE NOTICE:** The **PEDAL** firmware is implemented and working well, but at the moment I have only preconfigured the WAVESHARE_ESP32_S3_ZERO_N4R2 (Esp32s3 Mini) board..


### Supported Boards (ESP32-S3)

***Lightgun:***
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **WAVESHARE_ESP32_S3_PICO**

***Dongle:***
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **WAVESHARE_ESP32_S3_PICO**
- **LILYGO_T_DONGLE_S3** *(this one has an integrated DISPLAY)* 

***Pedal:***
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **WAVESHARE_ESP32_S3_PICO**
- **WAVESHARE_ESP32_S3_ZERO_N4R2**  *(Mini)*

---

### INSTALLATION

#### *OPTION 1: Manual Installation (.bin Files)*
This option is for those who wish to flash the firmware manually using the official [esptool](https://github.com/espressif/esptool) utility or the [NodeMCU PyFlasher](https://github.com/marcelstoer/nodemcu-pyflasher) graphical interface.

#### ***Lightgun:***

- **1. Base Version (NoFS):** Updates only the application part. **Does not overwrite** existing calibrations or configurations. If the filesystem is missing, it will still be created on first boot.
  - [OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
  - [OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin)

- **2. Full Version (Full):** Installs both firmware and factory filesystem. **Warning:** this procedure formats the microcontroller and deletes all previous data or settings.
  - [OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
  - [OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin)

#### ***Dongle:***

- [OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO.bin)
- [OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3.bin)

#### ***Pedal:***

- [OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO.bin)
- [OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2.bin)

---

#### *OPTION 2: Simplified Procedure (Ready-made Packages)*
As an alternative to the manual procedure, you can use these packages which already include the firmware files and tools for a guided installation.
Extract the entire contents of the ZIP into a folder on your PC. Then, run the **"flash_firmware"** script.
*Note for Windows users: the esptool.exe file might trigger an antivirus false positive; the file is safe and is extracted from official Espressif sources.*

#### ***Lightgun:***
The script will guide you through the installation, allowing you to choose between Base or Full versions, and will automatically search for the serial port.
* **Windows (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
* **Linux (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)

#### ***Dongle:***
The script will guide you through the installation and automatically search for the serial port.
* **Windows (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
  - [Download for LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-windows-64bit.zip)
* **Linux (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
  - [Download for LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)
  - [Download for LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-macos-64bit.zip)

#### ***Pedal:***
The script will guide you through the installation and automatically search for the serial port.
* **Windows (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-windows-64bit.zip)
* **Linux (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-macos-64bit.zip)

---

#### Compatibility
Compatibility between **lightgun**, **dongle**, and **pedal** is guaranteed using the firmware listed above.

To configure all lightgun features (pins, buttons, name, player, etc.), the original **OpenFIRE App** (commit **d0f3334** of April 19, 2026) **Release 3.0.3 - Tokinomiya Turbo: Hyper Fighting** is required:
[Download OpenFIRE App Release 3.0.3](https://github.com/TeamOpenFIRE/OpenFIRE-App/releases/tag/v3.0.3)

---
---

# GLOBAL RELEASE - FIRMWARE **LIGHTGUN** **DONGLE** **PEDAL** - OpenFIRE ESP32


Questa release è aggiornata al commit **8b651a2** (19 aprile 2026) del firmware originale [OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware).

> [!WARNING]
> **NOVITA! AVVISO SUL FIRMWARE PEDAL:** Il firmware **PEDAL** è implementato e ben funzionante, ma al momento ho preconfigurato solo la board WAVESHARE_ESP32_S3_ZERO_N4R2 (Esp32s3 Mini).

---

### Schede Supportate (ESP32-S3)

***Lightgun:***
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **WAVESHARE_ESP32_S3_PICO**

***Dongle:***
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **WAVESHARE_ESP32_S3_PICO**
- **LILYGO_T_DONGLE_S3** *(questo ha il DISPLAY integrato)* 

***Pedal:***
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **WAVESHARE_ESP32_S3_PICO**
- **WAVESHARE_ESP32_S3_ZERO_N4R2**  *(Mini)*

---

### INSTALLAZIONE

#### *OPZIONE 1: Installazione Manuale (File .bin)*
Questa opzione è dedicata a chi vuole caricare il firmware manualmente utilizzando l'utility ufficiale [esptool](https://github.com/espressif/esptool) o l'interfaccia grafica [NodeMCU PyFlasher](https://github.com/marcelstoer/nodemcu-pyflasher).

#### ***Lightgun:***

- **1. Versione Base (NoFS):** Aggiorna solo la parte applicativa. **Non sovrascrive** calibrazioni o configurazioni esistenti. Se il filesystem è assente, verrà comunque creato al primo avvio.
  - [OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
  - [OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin)

- **2. Versione Completa (Full):** Installa firmware e filesystem di fabbrica. **Attenzione:** questa procedura formatta il microcontrollore e cancella ogni dato o impostazione precedente.
  - [OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
  - [OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin)

#### ***Dongle:***

- [OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO.bin)
- [OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3.bin)

#### ***Pedal:***

- [OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO.bin)
- [OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2.bin)

---

#### *OPZIONE 2: Procedura Semplificata (Pacchetti pronti)*
In alternativa alla procedura manuale, puoi utilizzare questi pacchetti che includono già i file firmware e gli strumenti per l'installazione guidata.
Estrai l'intero contenuto dello ZIP in una cartella sul tuo PC. Successivamente, avvia lo script **"flash_firmware"**.
*Nota per utenti Windows: il file esptool.exe potrebbe generare un falso positivo dell'antivirus; il file è sicuro ed è estratto dai sorgenti originali Espressif.*

#### ***Lightgun:***
Lo script ti guiderà nell'installazione permettendoti di scegliere tra versione Base o Full e cercherà automaticamente la porta seriale.
* **Windows (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
* **Linux (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)

#### ***Dongle:***
Lo script ti guiderà nell'installazione e cercherà automaticamente la porta seriale.
* **Windows (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
  - [Download per LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-windows-64bit.zip)
* **Linux (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
  - [Download per LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)
  - [Download per LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-macos-64bit.zip)

#### ***Pedal:***
Lo script ti guiderà nell'installazione e cercherà automaticamente la porta seriale.
* **Windows (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-windows-64bit.zip)
* **Linux (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-macos-64bit.zip)


---

#### Compatibilità
La compatibilità tra **lightgun**, **dongle** e **pedal** è garantita utilizzando i firmware sopra elencati.

Per configurare la lightgun (pin, pulsanti, nome, player, ecc.) è richiesta l'**App OpenFIRE** originale (commit **d0f3334** del 19 aprile 2026) **Release 3.0.3 - Tokinomiya Turbo: Hyper Fighting**:
[Scarica OpenFIRE App Release 3.0.3](https://github.com/TeamOpenFIRE/OpenFIRE-App/releases/tag/v3.0.3)

---
