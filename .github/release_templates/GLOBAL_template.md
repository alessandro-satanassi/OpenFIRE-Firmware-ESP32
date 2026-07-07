<a id="english-version"></a>

<p align="center">
  <a href="#english-version"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

## Release Notes

**[OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware) Core:** Aligned to commit `8b651a2` of April 19, 2026 (version 6.2 - Long Bridge)

### Quick Install
Update the firmware directly from your browser via the WebFlasher: [Launch WebFlasher](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher)

---

### New Features
* *No new features introduced*

### Bug Fixes
* *No bugs fixed*

### Documentation
* *No documentation updates*

### Chores
* *No minor changes*

### Known Issues
* *No known issues*

---

### Supported Boards (ESP32-S3)

<table>
  <tr>
    <td align="center"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/board_scheme/ESP32S3-Devkit-C.svg" width="100%"></td>
    <td align="center"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/board_scheme/esp32-s3-pico.svg" width="100%"></td>
    <td align="center"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/board_scheme/esp32-s3-zero.svg" width="100%"></td>
    <td align="center"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/board_scheme/LILYGO-T-Dongle-S3-ESP32-S3.svg" width="100%"></td>
    <td align="center"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/board_scheme/esp32-s3-pocket-dongle-s3.svg" width="100%"></td>
  </tr>
</table>

***Lightgun:***
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **ESP32_S3_WROOM1_DevKitC_1_N8R2**
- **WAVESHARE_ESP32_S3_PICO**
- **WAVESHARE_ESP32_S3_ZERO_N8R8** *(Mini)*
- **WAVESHARE_ESP32_S3_ZERO_N4R2** *(Mini)*

***Dongle:***
- **LILYGO_T_DONGLE_S3** *(this one has an integrated DISPLAY)*
- **GNPE_POCKET_DONGLE_S3_N16R8** *(this one has an integrated DISPLAY)*
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **ESP32_S3_WROOM1_DevKitC_1_N8R2**
- **WAVESHARE_ESP32_S3_PICO**
- **WAVESHARE_ESP32_S3_ZERO_N8R8** *(Mini)*
- **WAVESHARE_ESP32_S3_ZERO_N4R2** *(Mini)*

***Pedal:***
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **ESP32_S3_WROOM1_DevKitC_1_N8R2**
- **WAVESHARE_ESP32_S3_PICO**
- **WAVESHARE_ESP32_S3_ZERO_N8R8** *(Mini)*
- **WAVESHARE_ESP32_S3_ZERO_N4R2** *(Mini)*

---

### INSTALLATION

#### WEB FLASHER (Recommended for all users)
The easiest, fastest, and safest way to install or update the firmware. It does not require installing any drivers or external software: it runs entirely within your browser.
* **Requirements:** PC/Mac with Google Chrome, Microsoft Edge, or Opera.

👉 **[LAUNCH OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?tag=<TAG>&lang=en)**

---

#### *FOR ADVANCED USERS:*
If your browser does not support the Web Flasher, or if you prefer to proceed via command line or external tools, you can download the specific files:

#### *OPTION 1: Manual Installation (.bin Files)*
This option is for those who wish to flash the firmware manually using the official [esptool](https://github.com/espressif/esptool) utility or the [NodeMCU PyFlasher](https://github.com/marcelstoer/nodemcu-pyflasher) graphical interface.

#### ***Lightgun:***

- **1. Base Version (NoFS):** Updates only the application part. **Does not overwrite** existing calibrations or configurations. If the filesystem is missing, it will still be created on first boot.
  - [OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
  - [OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin)
  - [OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin)
  - [OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_ZERO_N8R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_ZERO_N8R8.bin)
  - [OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_ZERO_N4R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_ZERO_N4R2.bin)

- **2. Full Version (Full):** Installs both firmware and factory filesystem. **Warning:** this procedure formats the microcontroller and deletes all previous data or settings.
  - [OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
  - [OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin)  
  - [OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin)
  - [OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N8R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N8R8.bin)
  - [OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N4R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N4R2.bin)

#### ***Dongle:***

- [OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3.bin)
- [OpenFIRE-DONGLE-GNPE_POCKET_DONGLE_S3_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-GNPE_POCKET_DONGLE_S3_N16R8.bin)
- [OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin)
- [OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO.bin)
- [OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N8R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N8R8.bin)
- [OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N4R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N4R2.bin)

#### ***Pedal:***

- [OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin)
- [OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO.bin)
- [OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8.bin)
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
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N8R2-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N8R8-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N4R2-windows-64bit.zip)
* **Linux (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N8R2-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N8R8-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N4R2-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N8R2-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N8R8-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N4R2-macos-64bit.zip)

#### ***Dongle:***
The script will guide you through the installation and automatically search for the serial port.
* **Windows (64bit)**
  - [Download for LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-windows-64bit.zip)
  - [Download for GNPE_POCKET_DONGLE_S3_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-GNPE_POCKET_DONGLE_S3_N16R8-windows-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N8R2-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N8R8-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N4R2-windows-64bit.zip)
* **Linux (64bit)**
  - [Download for LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-linux-amd-64bit.zip)
  - [Download for GNPE_POCKET_DONGLE_S3_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-GNPE_POCKET_DONGLE_S3_N16R8-linux-amd-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N8R2-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N8R8-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N4R2-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download for LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-macos-64bit.zip)
  - [Download for GNPE_POCKET_DONGLE_S3_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-GNPE_POCKET_DONGLE_S3_N16R8-macos-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N8R2-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N8R8-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N4R2-macos-64bit.zip)

#### ***Pedal:***
The script will guide you through the installation and automatically search for the serial port.
* **Windows (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N8R2-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-windows-64bit.zip)
* **Linux (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N8R2-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N8R2-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-macos-64bit.zip)

---

#### Compatibility
Compatibility between **lightgun**, **dongle**, and **pedal** is guaranteed using the firmware listed above.

To configure the lightgun (pins, buttons, calibration, name, player, etc.) the **[OpenFIRE App](https://github.com/TeamOpenFIRE/OpenFIRE-App)** is required. 
You can easily download the most recent pre-compiled binaries for your operating system (Windows, Linux, ARM) directly from our Tools portal:
👉 **[Download OpenFIRE App from ESP32 Ecosystem](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-Tools/?lang=en)**

---

<p align="center"> 🔸 🔸 🔸 </p>

---

<a id="versione-italiana"></a>

<p align="center">
  <a href="#english-version"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

## Note di rilascio

**Core di [OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware):** Allineato al commit `8b651a2` del 19 aprile 2026 (versione 6.2 - Long Bridge)

### Installazione Rapida
Aggiorna il firmware direttamente dal browser tramite il WebFlasher: [Avvia WebFlasher](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher)

---

### Nuove funzionalità
* *Nessuna nuova funzionalità introdotta*

### Correzione di bug
* *Nessun bug corretto*

### Documentazione
* *Nessun aggiornamento alla documentazione*

### Modifiche minori
* *Nessuna modifica minore*

### Problemi noti
* *Nessun problema noto*

---

### Schede Supportate (ESP32-S3)

<table>
  <tr>
    <td align="center"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/board_scheme/ESP32S3-Devkit-C.svg" width="100%"></td>
    <td align="center"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/board_scheme/esp32-s3-pico.svg" width="100%"></td>
    <td align="center"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/board_scheme/esp32-s3-zero.svg" width="100%"></td>
    <td align="center"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/board_scheme/LILYGO-T-Dongle-S3-ESP32-S3.svg" width="100%"></td>
    <td align="center"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/board_scheme/esp32-s3-pocket-dongle-s3.svg" width="100%"></td>
  </tr>
</table>

***Lightgun:***
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **ESP32_S3_WROOM1_DevKitC_1_N8R2**
- **WAVESHARE_ESP32_S3_PICO**
- **WAVESHARE_ESP32_S3_ZERO_N8R8** *(Mini)*
- **WAVESHARE_ESP32_S3_ZERO_N4R2** *(Mini)*

***Dongle:***
- **LILYGO_T_DONGLE_S3** *(questo ha il DISPLAY integrato)*
- **GNPE_POCKET_DONGLE_S3_N16R8** *(questo ha il DISPLAY integrato)*
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **ESP32_S3_WROOM1_DevKitC_1_N8R2**
- **WAVESHARE_ESP32_S3_PICO**
- **WAVESHARE_ESP32_S3_ZERO_N8R8** *(Mini)*
- **WAVESHARE_ESP32_S3_ZERO_N4R2** *(Mini)*

***Pedal:***
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **ESP32_S3_WROOM1_DevKitC_1_N8R2**
- **WAVESHARE_ESP32_S3_PICO**
- **WAVESHARE_ESP32_S3_ZERO_N8R8** *(Mini)*
- **WAVESHARE_ESP32_S3_ZERO_N4R2** *(Mini)*

---

### INSTALLAZIONE

####  WEB FLASHER (Consigliato per qualsiasi utente)
Il modo più semplice, veloce e sicuro per installare o aggiornare il firmware. Non richiede l'installazione di driver o software esterni: viene eseguito interamente dal tuo browser.
* **Requisiti:** PC/Mac con Google Chrome, Microsoft Edge o Opera.

👉 **[AVVIA OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?tag=<TAG>&lang=it)**

---

#### *PER UTENTI ESPERTI:*
Se il tuo browser non supporta il Web Flasher, o se preferisci procedere tramite riga di comando o tool esterni, puoi scaricare i file specifici:

#### *OPZIONE 1: Installazione Manuale (File .bin)*
Questa opzione è dedicata a chi vuole caricare il firmware manualmente utilizzando l'utility ufficiale [esptool](https://github.com/espressif/esptool) o l'interfaccia grafica [NodeMCU PyFlasher](https://github.com/marcelstoer/nodemcu-pyflasher).

#### ***Lightgun:***

- **1. Versione Base (NoFS):** Aggiorna solo la parte applicativa. **Non sovrascrive** calibrazioni o configurazioni esistenti. Se il filesystem è assente, verrà comunque creato al primo avvio.
  - [OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
  - [OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin)
  - [OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin)
  - [OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_ZERO_N8R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_ZERO_N8R8.bin)
  - [OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_ZERO_N4R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_ZERO_N4R2.bin)

- **2. Versione Completa (Full):** Installa firmware e filesystem di fabbrica. **Attenzione:** questa procedura formatta il microcontrollore e cancella ogni dato o impostazione precedente.
  - [OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
  - [OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin)  
  - [OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin)
  - [OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N8R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N8R8.bin)
  - [OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N4R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N4R2.bin)

#### ***Dongle:***

- [OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3.bin)
- [OpenFIRE-DONGLE-GNPE_POCKET_DONGLE_S3_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-GNPE_POCKET_DONGLE_S3_N16R8.bin)
- [OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin)
- [OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO.bin)
- [OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N8R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N8R8.bin)
- [OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N4R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N4R2.bin)

#### ***Pedal:***

- [OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N8R2.bin)
- [OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO.bin)
- [OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8.bin)
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
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N8R2-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N8R8-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N4R2-windows-64bit.zip)
* **Linux (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N8R2-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N8R8-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N4R2-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N8R2-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N8R8-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_ZERO_N4R2-macos-64bit.zip)

#### ***Dongle:***
Lo script ti guiderà nell'installazione e cercherà automaticamente la porta seriale.
* **Windows (64bit)**
  - [Download per LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-windows-64bit.zip)
  - [Download per GNPE_POCKET_DONGLE_S3_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-GNPE_POCKET_DONGLE_S3_N16R8-windows-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N8R2-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N8R8-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N4R2-windows-64bit.zip)
* **Linux (64bit)**
  - [Download per LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-linux-amd-64bit.zip)
  - [Download per GNPE_POCKET_DONGLE_S3_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-GNPE_POCKET_DONGLE_S3_N16R8-linux-amd-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N8R2-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N8R8-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N4R2-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download per LILYGO_T_DONGLE_S3](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-LILYGO_T_DONGLE_S3-macos-64bit.zip)
  - [Download per GNPE_POCKET_DONGLE_S3_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-GNPE_POCKET_DONGLE_S3_N16R8-macos-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-ESP32_S3_WROOM1_DevKitC_1_N8R2-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N8R8-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-DONGLE-WAVESHARE_ESP32_S3_ZERO_N4R2-macos-64bit.zip)

#### ***Pedal:***
Lo script ti guiderà nell'installazione e cercherà automaticamente la porta seriale.
* **Windows (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N8R2-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-windows-64bit.zip)
* **Linux (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N8R2-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N8R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-ESP32_S3_WROOM1_DevKitC_1_N8R2-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N8R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N8R8-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_ZERO_N4R2](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-PEDAL-WAVESHARE_ESP32_S3_ZERO_N4R2-macos-64bit.zip)

---

#### Compatibilità
La compatibilità tra **lightgun**, **dongle** e **pedal** è garantita utilizzando i firmware sopra elencati.

Per configurare la lightgun (pin, pulsanti, calibrazione, nome, player, ecc.) è necessaria l'**[App OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-App)**. 
Puoi scaricare facilmente gli eseguibili precompilati più recenti per il tuo sistema operativo (Windows, Linux, ARM) direttamente dal nostro portale Tools:
👉 **[Scarica OpenFIRE App dall'Ecosistema ESP32](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-Tools/?lang=it)**

---
