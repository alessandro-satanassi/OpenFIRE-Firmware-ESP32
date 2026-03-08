# OpenFIRE-Firmware-ESP32

> 🛠️ **Hardware sponsored by [PCBWay](https://www.pcbway.com)**
---

## ATTENZIONE

Per clonare correttamente il repository e compilarlo devi usare git clone col parametro --recursive, nel seguente modo:
git clone --recursive https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32.git

## WARNING

To properly clone the repository and compile it you must use git clone with the --recursive parameter, as follows:
git clone --recursive https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32.git


## ... porting OpenFIRE-firmware for EPS32

Questo repository è un porting del progetto originale 'OpenFIRE-firmware' del TeamOpenFIRE, adattato per funzionare sul microcontrollore ESP32S3.
Il progetto è stato sviluppato utilizzando PlatformIO e, salvo alcuni adattamenti e piccole modifiche, il codice rimane sostanzialmente fedele all'originale del TeamOpenFIRE.
Ho inoltre implementato la connessione wireless utilizzando il protocollo ESP-NOW per ESP32. Per farlo, ho sviluppato un dongle ESP32S3 da collegare al PC, il cui codice è disponibile nel mio repository 'OpenFIRE-DONGLE-ESP32'.
Grazie a questa soluzione, il PC non rileva alcuna differenza tra la connessione diretta via USB e quella wireless tramite dongle.
Il codice di questo repository funziona anche sul microcontrollore RP2040, ma in questo caso supporta solo il collegamento diretto tramite USB, senza connessione wireless.
Ogni volta che verranno apportate modifiche al progetto originale 'OpenFIRE-firmware' del TeamOpenFIRE, aggiornerò di conseguenza anche il codice di questo repository.
Ringrazio di cuore il TeamOpenFIRE per la creazione del progetto 'OpenFIRE-firmware': a loro vanno tutti i meriti e la mia piena gratitudine.
Questo è semplicemente un adattamento per il funzionamento su ESP32S3, con l'aggiunta della connessione wireless tramite ESP-NOW.

Questa è la sequenza di avvio del firmware:
1  - calibrazione automatica del joystick analogico (dura un paio di secondo durante i quali non devi toccare il joystick analogico);
2A - se il cavo USB è collegato al PC, la lightgun si connette direttamente via cavo USB
2B - se il cavo USB NON è collegato al PC, la lightgun cerca il canale wi-fi con minore interferenze e più libero per un trasmissione ottimale.
      Quindi si pone in ascolto su tale canale fino a quando non trova la presenza di un 'dongle', ma se nel frattempo viene collegato un cavo USB, si connette via cavo USB.
      Se rileva un 'dongle', prova a collegarsi ed associarsi a tale 'dongle' e se ci riesce si collega al PC tramite tale 'dongle', esattamente come farebbe con cavo USB ma via wifi. Il PC non nota alcuna differenza tra connessione via cavo USB o tramite dongle 'wifi', funzina tutta esattamente allo stesso modo.
      Nel caso in cui la lightgun si dovesse spegnere, al riavvio tenterà subito di riconnettersi all'ultimo 'dongle' associato e se lo trova, la connessione è immediata, altrimenti rinizia la ricerca del canale wifi più libero e si pone alla ricerca di un 'dongle' attivo libero.
3  - la lightgun si collega al PC, se è coolegata via wifi viene mostrata un'icona del wi-fi, altrimenti una incona di un cavo USB



This repository is a porting of the original 'OpenFIRE-firmware' project by TeamOpenFIRE, adapted to work on the ESP32S3 microcontroller.
The project was developed using PlatformIO and, apart from some adaptations and small adjustments, the code remains essentially faithful to the original by TeamOpenFIRE.
I have also implemented wireless connectivity using the ESP-NOW protocol for ESP32. To achieve this, I developed an ESP32S3 dongle to connect to the PC, and its code is available in my repository 'OpenFIRE-DONGLE-ESP32'.
Thanks to this solution, the PC does not detect any difference between a direct USB connection and a wireless connection via dongle.
The code in this repository also works on the RP2040 microcontroller, but in this case, it only supports direct USB connection without wireless connectivity.
Whenever there are changes to the original 'OpenFIRE-firmware' project by TeamOpenFIRE, I will adapt the code in this repository accordingly.
I sincerely thank TeamOpenFIRE for creating the 'OpenFIRE-firmware' project; all credit and gratitude go to them for their work.
This is simply an adaptation to make it work on ESP32S3, with the addition of wireless connectivity via ESP-NOW.

This is the firmware boot sequence:
1  - automatic calibration of the analog joystick (this lasts a couple of seconds, during which you must not touch the analog joystick);
2A - if the USB cable is connected to the PC, the lightgun connects directly via the USB cable.
2B - if the USB cable is NOT connected to the PC, the lightgun searches for the Wi-Fi channel with the least interference and the clearest signal for optimal transmission.
      It then listens on that channel until it detects the presence of a dongle, but if a USB cable is connected in the meantime, it connects via the USB cable.
      If it detects a dongle, it attempts to connect and pair with that dongle. If successful, it connects to the PC via that dongle, just as it would with a USB cable, but via Wi-Fi. The PC doesn't notice any difference between connecting via USB cable or via Wi-Fi dongle; everything works exactly the same.
      If the lightgun shuts down, upon restarting it will immediately attempt to reconnect to the last paired dongle. If it finds one, it connects immediately. Otherwise, it restarts the search for the freest Wi-Fi channel and searches for a free active dongle.
3 - The lightgun connects to the PC. If it is connected via Wi-Fi, a Wi-Fi icon is displayed; otherwise, a USB cable icon is displayed.


## 🤝 Sponsorship & Support

A special thanks to **[PCBWay](https://www.pcbway.com)** for sponsoring the hardware development of this project. Their professional PCB manufacturing service has been fundamental in transforming our schematics into reliable, high-quality physical boards.

We chose PCBWay for their:
* **Manufacturing Precision:** Excellent solder mask and silkscreen quality even on dense designs.
* **Reliability:** Consistent build quality across different batches.
* **Fast Prototyping:** Quick turnaround times that significantly accelerated our testing phase.

<p align="left">
  <a href="https://www.pcbway.com">
    <img src="docs/img/pcbway-logo.png" alt="PCBWay - PCB Prototype & Fabrication" width="200">
  </a>
</p>


![Waveshare-esp32-s3-pico](https://github.com/user-attachments/assets/5f7bf9ae-6ab5-4240-b930-d8cf20cb1c75)

![YD-esp32-s3-wroom1-DevKitC-1-N16R8](https://github.com/user-attachments/assets/6e865f2a-d90f-4dd0-9b57-c992bfd7377f)


## ... segue il link della la pagina originale del progetto
## ... follow the link to the original project page

### OpenFIRE - The Open *Four Infa-Red Emitter* Light Gun System [https://github.com/TeamOpenFIRE/OpenFIRE-Firmware](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware)

![BannerDark](docs/of_bannerLoD.png#gh-dark-mode-only)![BannerLight](docs/of_bannerDoL.png#gh-light-mode-only)

