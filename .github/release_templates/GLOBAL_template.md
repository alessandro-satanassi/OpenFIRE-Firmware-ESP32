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

* **Aim test (WebApp and desktop App):** a new **Open Aim Test...** button in **Gun Tests**, next to **Open IR Camera Tester...**, opens a fullscreen grid with nine targets to check the aim, the buttons and the recoil together. The crosshair follows the aim of the current profile and turns orange when only two or three emitters are seen; the trigger leaves a hit mark that fades out and fires the solenoid, or the rumble when the solenoid is disabled; the other buttons show their name. **H** hides the key line, **F** turns the recoil off and on, **Esc** closes. See [Test Mode](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/lightgun/src/README.md#test-mode).
* **Polling rate test (WebApp and desktop App):** measures how often the gun really sends its aim to the computer, with the USB cable or through the dongle. In the WebApp use **Measure Polling Rate** in the menu bar (also on the home page of the published WebApp), in the desktop App the button next to the COM port. With the gun not connected to the App and in mouse output, keep moving it: the test shows the current rate, the average, maximum and minimum and a graph. The expected value is about 209 Hz, the frame rate of the camera. Nothing is sent to the gun.
* **Wireless transmit power of the lightgun:** in **Gun Settings → Input / Output**, shown with **View > Show Unsafe Settings**, choose the power with which the gun transmits to the dongle and the pedal: one of the levels the ESP32-S3 radio really uses, from 2 to 20 dBm. The default, 16.5 dBm, is the power used so far. It takes effect at the next start of the gun. See [Operations Manual](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/lightgun/src/README.md#operations-manual).

### Improvements and Fixes

* **Data from the computer through the dongle no longer held back:** when the radio queue was full, part of the data sent by the desktop App or by MAMEHooker to a wireless gun could stay waiting until the next data from the computer. Now the rest is sent as soon as there is room.
* **Keyboard and gamepad over the dongle:** with the unified wireless packet, keyboard and gamepad reports now send the whole state of mouse, keyboard and gamepad together, as the mouse already did.
* **Damaged wireless packets are discarded:** a corrupted or foreign packet is now dropped before it is decoded, instead of writing past its end. This applies to the lightgun, the dongle and the pedal.
* This release includes everything introduced in 7.1.1 (tolerant calibration, more precise aim with the PAJ7025R3, longer solenoid timings) and in 7.1.0 (calibration with the four corners, parallax correction, guided lens refinement, IR LED placement assistant, wireless power saving and more). See the [7.1.1](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/tag/v7.1.1) and [7.1.0](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/tag/v7.1.0) release notes.
* Full list of changes in the [changelog](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/CHANGELOG.md#english-version).

### Installation and Usage Notes

* **Updating from 7.1.1 or 7.1.0:** a normal update keeps settings, pins and calibration profiles, and there is no need to calibrate again. Update the lightgun, the dongle and the pedal: in this release all three firmware change.
* The buttons and *Aim unavailable* in the aim test and the wireless transmit power need this lightgun firmware. The polling rate test works with any firmware.
* **Updating from 7.0.0:** a normal update keeps settings, pins and calibration profiles. Calibrate each profile again to get the corner refinement and, if you use it, the parallax correction.
* **Lightgun on a 4 MB board (Waveshare ESP32-S3 Zero N4R2), coming from 7.0.0:** 7.1.0 changed the flash partitions, so do a clean installation: settings, pins and calibration profiles go back to their defaults. Note them down before updating. The 4 MB dongle and pedal are not affected.
* **Coming from 6.2.1:** a clean installation is recommended; note your settings first, then configure and calibrate again. Layout files exported with 6.2.1 cannot be imported: set the pins again and export a new file.
* Use lightgun, dongle and pedal firmware from this same release.
* Lightgun wireless power saving is on by default (15 minutes, only on battery through the dongle): change it or disable it in **Gun Settings → Input / Output**.
* One firmware file per board, whatever the camera. A clean installation uses the same file but erases all settings and calibrations.
* To update a gun already running 7.0.0 or later, hold **Trigger + A** while switching on, for about 2 seconds, or let the Web Flasher restart it into update mode (the serial port then changes: select the new port and retry). For a first installation or a recovery, use the board's BOOT/RESET buttons. Install through the gun's own USB OTG port, not through the dongle.
* If the temperature sensor pin is assigned but no TMP36 is fitted, the solenoid does not fire. The default layouts of the ESP32-S3 DevKitC-1 and the Waveshare ESP32-S3 Pico assign this pin: without the sensor, enable **Use Custom Pins** in **Board Layout** and set **Temperature Sensor** to **Unmapped**.
* Starting the offline WebApp replaces the gun's USB serial port with USB networking, even if you configure over Wi-Fi. While a configuration App is connected, normal play is suspended. Save and restart normally after configuring.

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
The recommended way to install or update the firmware. It does not require installing any drivers or external software: it runs entirely within your browser.
* **Requirements:** PC/Mac with Google Chrome, Microsoft Edge, or Opera.

**[LAUNCH OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?tag=<TAG>&lang=en)**

---

#### *FOR ADVANCED USERS:*
If your browser does not support the Web Flasher, or if you prefer to proceed via command line or external tools, you can download the specific files:

#### *OPTION 1: Manual Installation (.bin Files)*
This option is for those who wish to flash the firmware manually using the official [esptool](https://github.com/espressif/esptool) utility or the [NodeMCU PyFlasher](https://github.com/marcelstoer/nodemcu-pyflasher) graphical interface.

#### ***Lightgun:***

A single firmware file is provided for each board. 
**Normal Update:** Flash the selected `.bin` file to address `0x0000`. This will not overwrite your existing calibrations or configurations (except the 4 MB Waveshare ESP32-S3 Zero N4R2 lightgun: see the installation notes above). If the filesystem is missing, it will be automatically created on the first boot using default settings.
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

* **Prova puntamento (WebApp e App desktop):** un nuovo pulsante **Apri Prova puntamento...** in **Gun Tests**, accanto ad **Apri Tester telecamera IR...**, apre una griglia a schermo intero con nove bersagli per controllare insieme la mira, i pulsanti e il rinculo. Il mirino segue la mira del profilo corrente e diventa arancione quando si vedono solo due o tre emettitori; il grilletto lascia un segno che svanisce e aziona il solenoide, oppure il rumble quando il solenoide è disattivato; gli altri pulsanti mostrano il loro nome. **H** nasconde la riga dei tasti, **F** disattiva e riattiva il rinculo, **Esc** chiude. Vedi [Modalità di Test](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/lightgun/src/README.md#modalità-di-test-italiano).
* **Misura polling rate (WebApp e App desktop):** misura quanto spesso la pistola invia davvero la mira al computer, con il cavo USB o tramite il dongle. Nella WebApp usa **Misura polling rate** nella barra dei menù (anche nella pagina iniziale della WebApp pubblicata), nell'App desktop il pulsante accanto alla Porta COM. Con la pistola non collegata all'App e in uscita mouse, continua a muoverla: il test mostra la frequenza attuale, la media, il massimo e il minimo e un grafico. Il valore atteso è circa 209 Hz, la frequenza della camera. Alla pistola non viene inviato nulla.
* **Potenza di trasmissione wireless della lightgun:** in **Impostazioni Gun → Input / Output**, visibile con **Visualizza > Mostra Impostazioni Non Sicure**, scegli la potenza con cui la pistola trasmette al dongle e al pedale: uno dei livelli che la radio dell'ESP32-S3 usa davvero, da 2 a 20 dBm. Il predefinito, 16,5 dBm, è la potenza usata finora. Vale dal prossimo avvio della pistola. Vedi [Manuale Operativo](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/lightgun/src/README.md#manuale-operativo-italiano).

### Miglioramenti e correzioni

* **I dati dal computer tramite il dongle non restano più indietro:** quando la coda radio era piena, una parte dei dati inviati dall'App desktop o da MAMEHooker a una pistola wireless poteva restare in attesa fino ai dati successivi del computer. Ora il resto viene inviato appena c'è spazio.
* **Tastiera e gamepad tramite il dongle:** con il pacchetto wireless unificato, anche i report di tastiera e gamepad inviano ora insieme lo stato intero di mouse, tastiera e gamepad, come già faceva il mouse.
* **I pacchetti wireless danneggiati vengono scartati:** un pacchetto corrotto o estraneo ora viene eliminato prima di essere decodificato, invece di scrivere oltre la sua fine. Vale per lightgun, dongle e pedale.
* Questa release comprende tutto ciò che è stato introdotto con la 7.1.1 (calibrazione tollerante, mira più precisa con la PAJ7025R3, tempi del solenoide più lunghi) e con la 7.1.0 (calibrazione con i quattro angoli, correzione parallasse, affinamento guidato della lente, assistente per il posizionamento dei LED IR, risparmio energetico wireless e altro). Vedi le note di rilascio della [7.1.1](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/tag/v7.1.1) e della [7.1.0](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/tag/v7.1.0).
* Elenco completo delle modifiche nella [cronologia modifiche](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/blob/main/CHANGELOG.md#versione-italiana).

### Note di installazione e utilizzo

* **Aggiornamento dalla 7.1.1 o dalla 7.1.0:** un aggiornamento normale mantiene impostazioni, pin e profili di calibrazione, e non serve calibrare di nuovo. Aggiorna lightgun, dongle e pedale: in questa release cambiano tutti e tre i firmware.
* I pulsanti e *Puntamento non disponibile* nella prova puntamento e la potenza di trasmissione wireless richiedono questo firmware della lightgun. La misura del polling rate funziona con qualsiasi firmware.
* **Aggiornamento dalla 7.0.0:** un aggiornamento normale mantiene impostazioni, pin e profili di calibrazione. Calibra di nuovo ogni profilo per avere l'affinamento con gli angoli e, se la usi, la correzione parallasse.
* **Lightgun su scheda da 4 MB (Waveshare ESP32-S3 Zero N4R2), se arrivi dalla 7.0.0:** la 7.1.0 ha cambiato le partizioni della flash, quindi fai un'installazione pulita: impostazioni, pin e profili di calibrazione tornano ai valori predefiniti. Annotali prima di aggiornare. Dongle e pedale da 4 MB non sono interessati.
* **Se arrivi dalla 6.2.1:** è consigliata un'installazione pulita; annota prima le impostazioni, poi configura e calibra di nuovo. I file di layout esportati con la 6.2.1 non si possono importare: reimposta i pin ed esporta un nuovo file.
* Usa firmware di questa stessa release per lightgun, dongle e pedale.
* Il risparmio energetico wireless della lightgun è attivo per impostazione predefinita (15 minuti, solo a batteria tramite il dongle): cambialo o disattivalo in **Impostazioni Gun → Input / Output**.
* Un solo file firmware per scheda, qualunque sia la telecamera. L'installazione pulita usa lo stesso file ma elimina tutte le impostazioni e calibrazioni.
* Per aggiornare una pistola che esegue già la 7.0.0 o successiva, tieni premuti **Grilletto + A** all'accensione, per circa 2 secondi, oppure lascia che il Web Flasher la riavvii in modalità aggiornamento (la porta seriale poi cambia: seleziona la nuova porta e riprova). Per la prima installazione o un recupero usa i pulsanti BOOT/RESET della scheda. Installa dalla porta USB OTG della pistola, non tramite il dongle.
* Se il pin del sensore di temperatura è assegnato ma il TMP36 non è montato, il solenoide non scatta. I layout predefiniti di ESP32-S3 DevKitC-1 e Waveshare ESP32-S3 Pico assegnano questo pin: senza il sensore, attiva **Usa pin personalizzati** in **Layout Scheda** e imposta **Sensore di temperatura** su **Non assegnato**.
* L'avvio della WebApp offline sostituisce la porta seriale USB della pistola con la rete USB, anche se configuri via Wi-Fi. Mentre un'App di configurazione è collegata, il gioco normale è sospeso. Al termine della configurazione salva e riavvia normalmente.

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
Il modo consigliato per installare o aggiornare il firmware. Non richiede l'installazione di driver o software esterni: viene eseguito interamente dal tuo browser.
* **Requisiti:** PC/Mac con Google Chrome, Microsoft Edge o Opera.

**[AVVIA OPENFIRE ESP32 WEB FLASHER](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?tag=<TAG>&lang=it)**

---

#### *PER UTENTI ESPERTI:*
Se il tuo browser non supporta il Web Flasher, o se preferisci procedere tramite riga di comando o tool esterni, puoi scaricare i file specifici:

#### *OPZIONE 1: Installazione Manuale (File .bin)*
Questa opzione è dedicata a chi vuole caricare il firmware manualmente utilizzando l'utility ufficiale [esptool](https://github.com/espressif/esptool) o l'interfaccia grafica [NodeMCU PyFlasher](https://github.com/marcelstoer/nodemcu-pyflasher).

#### ***Lightgun:***

Viene fornito un singolo file firmware per ogni scheda. 
**Aggiornamento normale:** Scrivi il file `.bin` scelto all'indirizzo `0x0000`. Non sovrascriverà calibrazioni o configurazioni esistenti (tranne la lightgun Waveshare ESP32-S3 Zero N4R2 da 4 MB: vedi le note di installazione sopra). Se il filesystem è assente, verrà creato automaticamente al primo avvio utilizzando le impostazioni predefinite.
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
