============== ENGLISH TRANSLATION ===============

This release is updated to commit 80166a6 (16 jan 2026) of the original firmware https://github.com/TeamOpenFIRE/OpenFIRE-Firmware
At the moment two Esp32-S3 boards are supported:

- **_ESP32 S3 WROOM-1 DevKitC-1 N16R8_**
- **_WAVESHARE ESP32 S3 PICO_**

The firmwares for both micros are provided below: I have prepared a full filesystem version that will overwrite any previous settings and a base (indicated with NoFS) that will not overwrite the previous calibration/configuration settings, but in any case if the filesystem is not yet present in the flash it will be created automatically.
FIRMWARE Base version (updates only the program or installs it if it is missing):

- [OpenFIREfw-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIREfw-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIREfw-NoFS-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIREfw-NoFS-WAVESHARE_ESP32_S3_PICO.bin)

FIRMWARE Complete version of filesystem (deletes all micro data):

- [OpenFIREfw-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIREfw-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIREfw-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIREfw-WAVESHARE_ESP32_S3_PICO.bin)

To upload the firmware on the micro you can use the command line utility "esptool" https://github.com/espressif/esptool or I recommend the following simple program with GUI interface https://github.com/marcelstoer/nodemcu-pyflasherhttps://github.com/marcelstoer/nodemcu-pyflasher

If you want things even simpler and more automatic, I have prepared a .zip package, one for each operating system (windows, linux, macos), that you can download to your computer. Just run the "flash_firmware" script that you will find inside the unzipped directory; the script will ask you if you want to install the basic or filesystem version of the firmware and the serial port where the micro is connected, but if you do not know it, it will automatically search for it. N.B. (the version of esptool.exe for Windows is identified as if it had a virus, but the file is clean, downloaded from the original archive).

Packages . zip for Windows:

- [OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
- [OpenFIRE-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)

Packages .zip for Linux:

- [OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
- [OpenFIRE-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)

Packages .zip for Macos:

- [OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
- [OpenFIRE-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)

For those who want to use the lightgun in wireless mode espnow using the wireless DONGLE, I indicate below the repository where to download the firmware for the DONGLE: https://github.com/alessandro-satanassi/OpenFIRE-DONGLE-ESP32/releases/tag/<TAG>
This firmware version is compatible with the original OpenFIRE App updated to commit 35e91f0 (16 Jan 2026) 
https://github.com/TeamOpenFIRE/OpenFIRE-App/actions/runs/21079466755 (linux)
https://github.com/TeamOpenFIRE/OpenFIRE-App/actions/runs/21079466753 (windows)

============== TRADUZIONE ITALIANA ===============

Questa release è aggiornata al commit 80166a6 (16 gennaio 2026) del firmware originale https://github.com/TeamOpenFIRE/OpenFIRE-Firmware
Al momento sono supportate due schede Esp32-S3:

- **_ESP32 S3 WROOM-1 DevKitC-1 N16R8_**
- **_WAVESHARE ESP32 S3 PICO_**

Vengono di seguito forniti i firmware di entrambi i micro: ho preparato una versione completa di filesystem che sovrascriverà qualsiasi precedente impostazione ed una base (indicata con NoFS) che non sovrascriverà le precedenti impostazioni di calibrazione/configurazione, ma comunque se il filesystem non è  ancora presente nella flash verrà creato in automatico.
FIRMWARE Versione base (aggiorna solo il programma o lo installa se assente):

- [OpenFIREfw-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIREfw-NoFS-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIREfw-NoFS-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIREfw-NoFS-WAVESHARE_ESP32_S3_PICO.bin)

FIRMWARE Versione completa di filesystem (cancella tutti i dati del micro):

- [OpenFIREfw-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIREfw-ESP32_S3_WROOM1_DevKitC_1_N16R8.bin)
- [OpenFIREfw-WAVESHARE_ESP32_S3_PICO.bin](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIREfw-WAVESHARE_ESP32_S3_PICO.bin)

Per fare l'upload del firmware sul micro potete usare  l'utility a riga di comando "esptool" https://github.com/espressif/esptool  oppure vi segnalo  il seguente  programmino semplice semplice con interfaccia GUI https://github.com/marcelstoer/nodemcu-pyflasherhttps://github.com/marcelstoer/nodemcu-pyflasher

Se volete le cose ancor più semplici ed automatiche, ho preparato un pacchetto .zip, uno per ogni sistema operativo (windows, linux, macos), che potete scaricare sul vostro computer. Basta eseguire lo script "flash_firmware" che troverete dentro alla directory scompattata; lo script vi chiederà se volete installare la versione base o con filesystem del firmware e la porta seriale ove è collegato il micro, ma se non la conoscete  la ricerca in automatico. N.B. (la versione di esptool.exe per windows viene identificato come se avesse un virus, ma il file è pulito, scaricato dall'archivio originale).

Pacchetti . zip Windows:

- [OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-windows-64bit.zip)
- [OpenFIRE-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-WAVESHARE_ESP32_S3_PICO-windows-64bit.zip)

Pacchetti .zip per Linux:

- [OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-linux-amd-64bit.zip)
- [OpenFIRE-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-WAVESHARE_ESP32_S3_PICO-linux-amd-64bit.zip)

Pacchetti .zip per Macos:

- [OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-ESP32_S3_WROOM1_DevKitC_1_N16R8-macos-64bit.zip)
- [OpenFIRE-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/releases/download/<TAG>/OpenFIRE-WAVESHARE_ESP32_S3_PICO-macos-64bit.zip)

Per chi vuole utilizzare la lightgun in modalità wireless espnow utilizzando il DONGLE wireless, indico di seguito il repository ove scaricare il firmware per il DONGLE: https://github.com/alessandro-satanassi/OpenFIRE-DONGLE-ESP32/releases/tag/<TAG>
Questa versione del firmware è compatibile con l'App OpenFIRE originale aggiornata al commit 35e91f0 (16 gennaio 2026)  
https://github.com/TeamOpenFIRE/OpenFIRE-App/actions/runs/21079466755 (linux)
https://github.com/TeamOpenFIRE/OpenFIRE-App/actions/runs/21079466753 (windows)
