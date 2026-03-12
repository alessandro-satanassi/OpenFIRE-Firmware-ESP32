# OpenFIRE-Firmware-ESP32

> 🛠️ **Hardware sponsored by [PCBWay](https://www.pcbway.com)**
---

---
### 🔫 PICON-AS: Wireless Lightgun Hardware
Se vuoi costruire una lightgun wireless funzionante completamente a batteria basata su questo firmware, dai un'occhiata al progetto **PICON-AS**. Trovi tutte le istruzioni per l'hardware e il montaggio a questo indirizzo: [Istruzioni PICON-AS](https://alessandro-satanassi.github.io/OpenFIRE-PICON-AS-ESP32/)
Nota: Il sito è già fruibile e i contenuti tecnici sono corretti; stiamo solo ultimando la stesura dettagliata di alcune sezioni delle istruzioni.
---
If you want to build a fully battery-powered wireless lightgun based on this firmware, check out the **PICON-AS** project. You can find all the hardware and assembly instructions at this link: [PICON-AS Instructions](https://alessandro-satanassi.github.io/OpenFIRE-PICON-AS-ESP32/)
Note: The website is already functional and all technical data is correct; we are currently finalizing the detailed step-by-step instructions.
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

## Sequenza di avvio del firmware

1. **Calibrazione automatica del joystick analogico**
   > ⚠️ **IMPORTANTE:** La procedura dura circa due secondi. Durante questa fase **non toccare il joystick** per garantire una calibrazione corretta.

2. **Gestione della connettività**
    * **2A - Connessione via Cavo:** Se il cavo USB è collegato al PC, la lightgun si connette direttamente in modalità cablata.
    * **2B - Connessione Wireless (USB scollegato):**
        * La lightgun scansiona l'ambiente e seleziona il canale Wi-Fi con minori interferenze per una trasmissione ottimale.
        * Rimane in ascolto di un **dongle**. Se durante l'attesa viene collegato un cavo USB, la connessione passa immediatamente su cavo.
        * Una volta rilevato un dongle libero, esegue l'associazione. Il PC gestirà la periferica esattamente come se fosse collegata via cavo (nessuna differenza di funzionamento).
        * **Riconnessione automatica:** In caso di spegnimento, al riavvio cercherà prioritariamente l'ultimo dongle associato per una connessione istantanea. Se non lo trova, farà una nuova scansione dei canali.

3. **Stato della connessione**
   Una volta stabilito il collegamento con il PC, l'interfaccia mostrerà lo stato attuale:
   * 🛜 **Icona Wi-Fi**: Connessione tramite dongle wireless.
   * 🔌 **Icona USB**: Connessione tramite cavo fisico.


This repository is a porting of the original 'OpenFIRE-firmware' project by TeamOpenFIRE, adapted to work on the ESP32S3 microcontroller.
The project was developed using PlatformIO and, apart from some adaptations and small adjustments, the code remains essentially faithful to the original by TeamOpenFIRE.
I have also implemented wireless connectivity using the ESP-NOW protocol for ESP32. To achieve this, I developed an ESP32S3 dongle to connect to the PC, and its code is available in my repository 'OpenFIRE-DONGLE-ESP32'.
Thanks to this solution, the PC does not detect any difference between a direct USB connection and a wireless connection via dongle.
The code in this repository also works on the RP2040 microcontroller, but in this case, it only supports direct USB connection without wireless connectivity.
Whenever there are changes to the original 'OpenFIRE-firmware' project by TeamOpenFIRE, I will adapt the code in this repository accordingly.
I sincerely thank TeamOpenFIRE for creating the 'OpenFIRE-firmware' project; all credit and gratitude go to them for their work.
This is simply an adaptation to make it work on ESP32S3, with the addition of wireless connectivity via ESP-NOW.

## Firmware Boot Sequence

1. **Automatic Joystick Calibration**
   > ⚠️ **IMPORTANT:** The process takes about two seconds. **Do not touch the analog joystick** during this phase to ensure a correct calibration.

2. **Connectivity Management**
    * **2A - Wired Connection:** If the USB cable is connected to the PC, the lightgun will connect directly via USB.
    * **2B - Wireless Connection (USB disconnected):**
        * The lightgun scans for the Wi-Fi channel with the lowest interference for optimal transmission.
        * It listens for an available **dongle**. If a USB cable is connected during this stage, it will switch to wired mode immediately.
        * Once a dongle is detected and paired, the connection is established wirelessly. The PC will treat the device exactly like a standard USB connection.
        * **Auto-Reconnect:** If the lightgun powers off, it will prioritize reconnecting to the last paired dongle for an instant link. If not found, it restarts the channel scan and dongle search.

3. **Connection Status**
   Once the connection to the PC is established, the interface will display the current status:
   * 🛜 **Wi-Fi Icon**: Connected via wireless dongle.
   * 🔌 **USB Icon**: Connected via physical cable.


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

