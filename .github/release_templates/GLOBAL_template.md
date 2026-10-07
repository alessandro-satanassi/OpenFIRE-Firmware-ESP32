<a id="english-version"></a>

<p align="center">
  <a href="#english-version"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

## Release Notes

<details>
<summary>Original OpenFIRE core reference</summary>

This ESP32 release is based on the [OpenFIRE core](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware) at commit `8b651a2` of April 19, 2026 (version 6.2 - Long Bridge). This is the upstream reference, not the version to choose for the ESP32 devices below.

</details>

### Quick Install
Install or update the firmware directly from your browser with the Web Flasher: [Launch Web Flasher](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?tag=<TAG>&lang=en)

New to OpenFIRE? Start from [Getting Started](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32#getting-started) · [Common Problems](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/lightgun/src/README.md#common-problems) · [Full changelog](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/CHANGELOG.md)

---

### New Features

* Configuration with the new [WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/?lang=en): online from Chrome or Edge, through the gun's USB cable or its paired dongle, or offline from the WebApp stored in the gun, over Wi-Fi (phones included) or USB networking.
* DFRobot/Wii, PAJ7025R2 and PAJ7025R3 cameras in the same firmware, chosen in the WebApp. Use 940 nm emitters for DFRobot/Wii and 850 nm for PixArt.
* Hold **B** while switching on the gun, for about 2 seconds, to start the offline WebApp: join **OpenFIRE_Config**, accept the network without Internet, then open **http://openfire.local/** or **http://192.168.4.1/** in your normal browser. Over USB, use **http://192.168.7.1/** on computers that support it.
* Hold **Trigger + A** while switching on, for about 2 seconds, to put a gun already running 7.0.0 into firmware update mode. For a first installation or a recovery, use the board's BOOT/RESET buttons. Install through the gun's own USB OTG port, not through the dongle.
* Calibration from the WebApp or the desktop App checks the IR emitters at every target: green crosshair when all four are seen well, orange when one is weak, red (shot refused) when one is missing, with a panel showing each emitter, progress dots and a legend.
* Calibration from the gun: the cursor traces a small circle on each target, and an OLED display shows where the target is, the step and how to confirm.
* Clearer IR camera test: each emitter circle shows the size and brightness of the light spot seen by the camera; emitters not seen are marked with a red X, and a legend explains every symbol.
* Saved mouse/gamepad startup mode, wireless pedal option and Up/Down navigation in the pause menu.

### Improvements and Fixes

* Refined IR tracking and faster recovery of temporarily hidden LEDs; the Square layout supports vertical and wide rectangles.
* Hotkey pause mode: the rumble/solenoid toggles now work when no hardware switch is fitted.
* After the first calibration started from the gun on a new board, the gun returns to normal operation once you confirm it. As before, this first calibration is saved automatically; calibrations started from pause mode must be saved separately.
* Rewritten bilingual user guides: see the [operating manual](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/lightgun/src/README.md#english-version).

### Installation and Usage Notes

* One firmware file per board, whatever the camera. A clean installation uses the same file but erases all settings and calibrations.
* A clean installation is recommended when upgrading from 6.2.1: note your settings first, then configure and calibrate again.
* Use lightgun, dongle and pedal firmware from this same release.
* Starting the offline WebApp replaces the gun's USB serial port with USB networking, even if you configure over Wi-Fi. While a configuration App is connected, normal play is suspended. Save and restart normally after configuring.
* The Web Flasher can also restart a gun running 7.0.0 into update mode by itself. After this restart the serial port changes: select the new port and retry.

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

**[LAUNCH OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?tag=<TAG>&lang=en)**

---

#### *FOR ADVANCED USERS:*
If your browser does not support the Web Flasher, or if you prefer to proceed via command line or external tools, you can download the specific files:

#### *OPTION 1: Manual Installation (.bin Files)*
This option is for those who wish to flash the firmware manually using the official [esptool](https://github.com/espressif/esptool) utility or the [NodeMCU PyFlasher](https://github.com/marcelstoer/nodemcu-pyflasher) graphical interface.

#### ***Lightgun:***

A single firmware file is provided for each board. 
**Normal Update:** Flash the selected `.bin` file to address `0x0000`. This will not overwrite your existing calibrations or configurations. If the filesystem is missing, it will be automatically created on the first boot using default settings.
**Clean Install:** if you want a 'clean' installation, use the `--erase-all` option alongside the esptool `write-flash` command to erase the entire flash memory before writing the firmware. **Warning: with a 'clean' installation, all previous data on the flash memory, including settings and calibrations, will be permanently deleted.**

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
*Note for Windows users: the original Espressif esptool.exe may trigger an antivirus alert. Download this package from the official Releases page and check its origin before allowing it to run; do not disable your antivirus globally.*

#### ***Lightgun:***
The script will guide you through the installation and automatically search for the serial port. During the guided installation, press **Enter** for a normal update, or type **1** for a 'clean' install which **deletes all data on the flash memory, including settings and calibrations**.
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

Configure the lightgun with the official **[WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/?lang=en)** using Chrome or Edge on a computer, or with the integrated offline WebApp. Configuration Apps for earlier firmware are not compatible with the new configuration protocol.

---

<a id="versione-italiana"></a>

<p align="center">
  <a href="#english-version"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/raw/main/docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

## Note di rilascio

<details>
<summary>Riferimento al core OpenFIRE originale</summary>

Questa release ESP32 si basa sul [core OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware) al commit `8b651a2` del 19 aprile 2026 (versione 6.2 - Long Bridge). È il riferimento al progetto originale, non la versione da scegliere per i dispositivi ESP32 elencati sotto.

</details>

### Installazione Rapida
Installa o aggiorna il firmware direttamente dal browser con il Web Flasher: [Avvia Web Flasher](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?tag=<TAG>&lang=it)

Non conosci ancora OpenFIRE? Parti dai [Primi passi](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32#primi-passi) · [Problemi comuni](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/lightgun/src/README.md#problemi-comuni-italiano) · [Cronologia modifiche completa](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/CHANGELOG.md#versione-italiana)

---

### Nuove funzionalità

* Configurazione con la nuova [WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/?lang=it): online da Chrome o Edge, tramite il cavo USB della pistola o il suo dongle associato, oppure offline con la WebApp contenuta nella pistola, via Wi-Fi (anche da telefono) o rete USB.
* Telecamere DFRobot/Wii, PAJ7025R2 e PAJ7025R3 nello stesso firmware, da scegliere nella WebApp. Usa emettitori da 940 nm per DFRobot/Wii e da 850 nm per PixArt.
* Tieni premuto **B** all'accensione della pistola, per circa 2 secondi, per avviare la WebApp offline: collegati a **OpenFIRE_Config**, accetta la rete senza Internet, poi apri **http://openfire.local/** o **http://192.168.4.1/** nel browser normale. Via USB usa **http://192.168.7.1/** sui computer che la supportano.
* Tieni premuti **Grilletto + A** all'accensione, per circa 2 secondi, per portare in modalità aggiornamento firmware una pistola che esegue già la 7.0.0. Per la prima installazione o un recupero usa i pulsanti BOOT/RESET della scheda. Installa dalla porta USB OTG della pistola, non tramite il dongle.
* La calibrazione dalla WebApp o dall'App desktop controlla gli emettitori IR a ogni bersaglio: mirino verde quando tutti e quattro sono visti bene, arancione quando uno è debole, rosso (tiro rifiutato) quando ne manca uno, con un riquadro che mostra ogni emettitore, pallini di avanzamento e una legenda.
* Calibrazione dalla pistola: il cursore descrive un piccolo cerchio su ogni bersaglio, e un display OLED mostra dove si trova il bersaglio, il passo e come confermare.
* Test della telecamera IR più chiaro: ogni cerchio degli emettitori mostra grandezza e luminosità della macchia di luce vista dalla telecamera; gli emettitori non visti sono segnati con una X rossa e una legenda spiega ogni simbolo.
* Modalità mouse/gamepad salvabile per l'avvio, opzione pedale wireless e navigazione Su/Giù nel menu di pausa.

### Miglioramenti e correzioni

* Tracciamento IR affinato e recupero più rapido dei LED temporaneamente nascosti; il layout Square supporta rettangoli verticali e larghi.
* Modalità pausa Hotkey: i comandi per rumble e solenoide ora funzionano quando non è montato un interruttore fisico.
* Dopo la prima calibrazione avviata dalla pistola su una scheda nuova, la pistola torna al funzionamento normale appena la confermi. Come in precedenza, questa prima calibrazione viene salvata automaticamente; quelle avviate dalla pausa vanno salvate separatamente.
* Guide utente bilingui riscritte: consulta il [manuale operativo](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/lightgun/src/README.md#versione-italiana).

### Note di installazione e utilizzo

* Un solo file firmware per scheda, qualunque sia la telecamera. L'installazione pulita usa lo stesso file ma elimina tutte le impostazioni e calibrazioni.
* Passando dalla 6.2.1 è consigliata un'installazione pulita: annota prima le impostazioni, poi configura e calibra di nuovo.
* Usa firmware di questa stessa release per lightgun, dongle e pedale.
* L'avvio della WebApp offline sostituisce la porta seriale USB della pistola con la rete USB, anche se configuri via Wi-Fi. Mentre un'App di configurazione è collegata, il gioco normale è sospeso. Al termine della configurazione salva e riavvia normalmente.
* Il Web Flasher può anche riavviare da solo in modalità aggiornamento una pistola che esegue la 7.0.0. Dopo questo riavvio la porta seriale cambia: seleziona la nuova porta e riprova.

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

**[AVVIA OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?tag=<TAG>&lang=it)**

---

#### *PER UTENTI ESPERTI:*
Se il tuo browser non supporta il Web Flasher, o se preferisci procedere tramite riga di comando o tool esterni, puoi scaricare i file specifici:

#### *OPZIONE 1: Installazione Manuale (File .bin)*
Questa opzione è dedicata a chi vuole caricare il firmware manualmente utilizzando l'utility ufficiale [esptool](https://github.com/espressif/esptool) o l'interfaccia grafica [NodeMCU PyFlasher](https://github.com/marcelstoer/nodemcu-pyflasher).

#### ***Lightgun:***

Viene fornito un singolo file firmware per ogni scheda. 
**Aggiornamento normale:** Scrivi il file `.bin` scelto all'indirizzo `0x0000`. Non sovrascriverà calibrazioni o configurazioni esistenti. Se il filesystem è assente, verrà creato automaticamente al primo avvio utilizzando le impostazioni predefinite.
**Installazione pulita:** se vuoi un'installazione 'pulita', usa l'opzione `--erase-all` assieme al comando `write-flash` di esptool per cancellare l'intera flash prima di scrivere il firmware. **Attenzione: con l'installazione 'pulita' tutti i dati precedenti nella flash, comprese impostazioni e calibrazioni, verranno eliminati definitivamente.**

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
*Nota per utenti Windows: l'esptool.exe originale di Espressif può generare un avviso dell'antivirus. Scarica questo pacchetto dalla pagina Releases ufficiale e verificane la provenienza prima di consentirne l'esecuzione; non disattivare l'antivirus globalmente.*

#### ***Lightgun:***
Lo script ti guiderà nell'installazione e cercherà automaticamente la porta seriale. Durante l'installazione guidata premi **Invio** per un aggiornamento normale, oppure digita **1** per un'installazione 'pulita' che **cancella tutti i dati nella flash, comprese impostazioni e calibrazioni**.
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

Configura la lightgun con la **[WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/?lang=it)** ufficiale usando Chrome o Edge su computer, oppure con la WebApp offline integrata. Le App di configurazione dei firmware precedenti non sono compatibili con il nuovo protocollo di configurazione.

---
