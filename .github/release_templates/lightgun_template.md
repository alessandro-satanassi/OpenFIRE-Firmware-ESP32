# **LIGHTGUN** FIRMWARE - OpenFIRE ESP32


This release is updated to commit **80166a6** (January 16, 2026) of the original [OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware) firmware.

### Supported Boards (ESP32-S3)
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **WAVESHARE_ESP32_S3_PICO**

---

### INSTALLATION

#### *OPTION 1: Manual Installation (.bin Files)*
This option is for those who wish to flash the firmware manually using the official [esptool](https://github.com/espressif/esptool) utility or the [NodeMCU PyFlasher](https://github.com/marcelstoer/nodemcu-pyflasher) graphical interface.

**1. Base Version (NoFS):** Updates only the application part. **Does not overwrite** existing calibrations or configurations. If the filesystem is missing, it will still be created on first boot.
- [OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin)

**2. Full Version (Full):** Installs both firmware and factory filesystem. **Warning:** this procedure formats the microcontroller and deletes all previous data or settings.
- [OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin)

---

#### *OPTION 2: Simplified Procedure (Ready-made Packages)*
As an alternative to the manual procedure, you can use these packages which already include the firmware files and tools for a guided installation.
Extract the entire contents of the ZIP into a folder on your PC. Then, run the **"flash_firmware"** script. The script will allow you to choose between Base or Full versions and will automatically search for the serial port.
*Note for Windows users: the esptool.exe file might trigger an antivirus false positive; the file is safe and is extracted from official Espressif sources.*

* **Windows (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
* **Linux (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download for ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download for WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)

---

#### Compatibility
This version of the **Lightgun** firmware is compatible with:
* **Dongle** firmware version: ______ 
* **Pedal** firmware version: ______ 

To configure all lightgun features (pins, buttons, name, player, etc.), the original **OpenFIRE App** (commit **35e91f0** of January 16, 2026) **or later versions** is required:
[Download for Linux](https://github.com/TeamOpenFIRE/OpenFIRE-App/actions/runs/21079466755) | [Download for Windows](https://github.com/TeamOpenFIRE/OpenFIRE-App/actions/runs/21079466753)

---
---

# FIRMWARE **LIGHTGUN** - OpenFIRE ESP32


Questa release è aggiornata al commit **80166a6** (16 gennaio 2026) del firmware originale [OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware).

### Schede Supportate (ESP32-S3)
- **ESP32_S3_WROOM1_DevKitC_1_N16R8**
- **WAVESHARE_ESP32_S3_PICO**

---

### INSTALLAZIONE

#### *OPZIONE 1: Installazione Manuale (File .bin)*
Questa opzione è dedicata a chi vuole caricare il firmware manualmente utilizzando l'utility ufficiale [esptool](https://github.com/espressif/esptool) o l'interfaccia grafica [NodeMCU PyFlasher](https://github.com/marcelstoer/nodemcu-pyflasher).

**1. Versione Base (NoFS):** Aggiorna solo la parte applicativa. **Non sovrascrive** calibrazioni o configurazioni esistenti. Se il filesystem è assente, verrà comunque creato al primo avvio.
- [OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-NoFS-WAVESHARE_ESP32_S3_PICO.bin)

**2. Versione Completa (Full):** Installa firmware e filesystem di fabbrica. **Attenzione:** questa procedura formatta il microcontrollore e cancella ogni dato o impostazione precedente.
- [OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO.bin)

---

#### *OPZIONE 2: Procedura Semplificata (Pacchetti pronti)*
In alternativa alla procedura manuale, puoi utilizzare questi pacchetti che includono già i file firmware e gli strumenti per l'installazione guidata.
Estrai l'intero contenuto dello ZIP in una cartella sul tuo PC. Successivamente, avvia lo script **"flash_firmware"**. Lo script ti permetterà di scegliere tra versione Base o Full e cercherà automaticamente la porta seriale.
*Nota per utenti Windows: il file esptool.exe potrebbe generare un falso positivo dell'antivirus; il file è sicuro ed è estratto dai sorgenti originali Espressif.*

* **Windows (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)
* **Linux (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)
* **MacOS (64bit)**
  - [Download per ESP32_S3_WROOM1_DevKitC_1_N16R8](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
  - [Download per WAVESHARE_ESP32_S3_PICO](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-LIGHTGUN-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)

---

#### Compatibilità
Questa versione del firmware **Lightgun** è compatibile con:
* Firmware **Dongle** versione: ______ 
* Firmware **Pedale** versione: ______ 

Per configurare la lightgun (pin, pulsanti, nome, player, ecc.) è richiesta l'**App OpenFIRE** originale (commit **35e91f0** del 16 gennaio 2026) **o versioni successive**:
[Scarica per Linux](https://github.com/TeamOpenFIRE/OpenFIRE-App/actions/runs/21079466755) | [Scarica per Windows](https://github.com/TeamOpenFIRE/OpenFIRE-App/actions/runs/21079466753)

---
