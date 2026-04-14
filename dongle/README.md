## OpenFIRE-DONGLE for EPS32


![foto_dongle](https://github.com/user-attachments/assets/a38d7390-f7e7-42d7-844c-ad8388670f6d)



https://github.com/user-attachments/assets/70a193a0-686e-4da3-8d30-f0067de6b63c



Questo repository è stato creato per essere utilizzato insieme al codice del progetto OpenFIRE-Firmware-ESP32, un porting del progetto originale 'OpenFIRE-firmware' del TeamOpenFIRE, adattato per funzionare sul microcontrollore ESP32S3.
Il progetto, sviluppato utilizzando PlatformIO, rappresenta il firmware per un ESP32S3 da usare come dongle collegato al PC, al fine di abilitare una connessione wireless tramite il protocollo ESP-NOW di ESP32. Questo dongle è progettato per essere usato in combinazione con il firmware 'OpenFIRE-Firmware-ESP32', da installare sulla lightgun.
Il codice è strutturato per rilevare automaticamente le lightgun e configurarsi in modo autonomo.
La trasmissione tra lightgun e dongle è bidirezionale, consentendo di utilizzare la lightgun come se fosse connessa direttamente al PC via USB.
Il PC non rileva alcuna differenza tra una connessione diretta tramite USB e una connessione wireless tramite dongle.
Desidero esprimere la mia sincera gratitudine al TeamOpenFIRE per la creazione del progetto 'OpenFIRE': tutti i meriti e la mia piena riconoscenza vanno a loro.
Una piccola parte del codice è stata estratta e modificata dal repository "SerialTransfer": https://github.com/PowerBroker2/SerialTransfer.git.
Desidero pertanto ringraziare anche l'autore di tale codice, PowerBroker2, per il suo contributo.


This repository was created to be used together with the code from the OpenFIRE-Firmware-ESP32 project, a porting of the original 'OpenFIRE-firmware' project by TeamOpenFIRE, adapted to work on the ESP32S3 microcontroller.
The project, developed using PlatformIO, serves as the firmware for an ESP32S3 to be used as a dongle connected to the PC, enabling wireless connectivity through the ESP-NOW protocol of ESP32. This dongle is designed to be used in combination with the firmware 'OpenFIRE-Firmware-ESP32', which is installed on the lightgun.
The code is structured to automatically detect lightguns and configure itself independently.
The communication between the lightgun and the dongle is bidirectional, allowing the lightgun to function as if it were directly connected to the PC via USB.
The PC does not notice any difference between a direct USB connection and a wireless connection through the dongle.
I would like to express my sincere gratitude to TeamOpenFIRE for creating the 'OpenFIRE' project; all credit and my deepest appreciation go to them.
A small portion of the code was taken and modified from the repository "SerialTransfer": https://github.com/PowerBroker2/SerialTransfer.git.
I would therefore like to thank the author of this code, PowerBroker2, for their contribution.

