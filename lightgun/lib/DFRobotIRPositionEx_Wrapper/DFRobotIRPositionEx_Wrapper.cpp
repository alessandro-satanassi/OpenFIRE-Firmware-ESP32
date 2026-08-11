#ifdef PAJ7025_CAM


#include "DFRobotIRPositionEx_Wrapper.h"
#include "OpenFIREConst.h"

#ifdef ARDUINO_ARCH_ESP32  // [ESP32_PORT]
    #define delay(ms) vTaskDelay(pdMS_TO_TICKS(ms))                    
#endif //ARDUINO_ARCH_ESP32


DFRobotIRPositionEx::DFRobotIRPositionEx(SPIClass* spiPort, int8_t csPin) : _spi(spiPort), _csPin(csPin) {
    current_format = DataFormat_Basic;
    seenFlags = 0;
    for(int i=0; i<4; i++){
        positionX[i] = 1023; // Valore fittizio tipico per 'vuoto'
        positionY[i] = 1023;
        unpackedSizes[i] = 15;
    }
}

DFRobotIRPositionEx::~DFRobotIRPositionEx() {
}

bool DFRobotIRPositionEx::begin(uint32_t clock, DataFormat_e format, Sensitivity_e sensitivity) {
    // Inizializza la PAJ7025
    if (!cam.begin(_spi, (uint8_t)_csPin), PAJ7025_SPI_CLOCK_2MHZ) {
        return false;
    }
    
    // ===========================================================
    // configurazione personalizzata di fps. resolution, ecc.ecc.
    // =========================================================== 
    
    // imposta frame rate
    cam.setFrameRate(200); // 200 fps

    // imposta velocità otturatore  (provvisoriamente per led 940nm provare 1000 - la CAM lavora bene con led 850nm)
    cam.setExposure(300); // ideale tra 200 e 400 -- più è basso meglio è per evitare effetto scia o accecamento, ma se troppo basso si rischia di oscurare i led
    
    // imposta guadagno
    //cam.setGain(0x10, 0x00); // minima x2 (default)
    ////cam.setGain(0x10, 0x02); // alta   x4
    ////cam.setGain(0x10, 0x03); // massima x8
    // Converti la sensitivity nel gain/dsp della PAJ7025
    sensitivityLevel(sensitivity);

    // ========= impostazioni DSP ====================
    
    // area_min (default 0x03 - 3)
    // brightness_th (defaul 0x97 - 151)
    // area_max (default 0x2585 - 9605)
    // noise_th (default 0x0A - 10)
    //cam.setDSP(2, 130, 150, 40);
    
    // risoluzione massima 4095x4095 come richiesto per il calcolo mouse
    cam.setResolution(CamMaxX, CamMaxY);

    // imposta massimo numero oggetti relevabili  1 - 16
    cam.setMaxObjects(4); // Limita ai 4 blob che ti interessano

    // ======= fine impostazioni DSP ==================

    // Imposta il formato
    dataFormat(format);

    return true;
}

void DFRobotIRPositionEx::dataFormat(DataFormat_e format) {
    current_format = format;
}

void DFRobotIRPositionEx::sensitivityLevel(Sensitivity_e sensitivity) {
    // Configura dinamicamente il guadagno (Gain) e i filtri di elaborazione interna (DSP)
    // del sensore PAJ7025 in base al livello di sensibilità selezionato dall'utente.
    // L'ordine degli 'if' garantisce la priorità dei profili specifici rispetto ai fallback.

    if (sensitivity == Sensitivity_Min) {
        // Guadagno minimo di fabbrica (2X) - Segnale pulito, zero rumore termico.
        cam.setGain(0x10, 0x00); 
        
        // DSP ottimizzato per la "Golden Zone" (distanza di gioco ideale: 1.0m - 1.5m).
        // - area_min (2): Scarta il rumore di fondo ma intercetta i pixel vicini.
        // - brightness_th (130): Taglia i riflessi deboli della stanza.
        // - area_max (150): Muro anti-finestra/lampade ma sicuro per il blooming ravvicinato.
        // - noise_th (40): "Incolla" efficacemente i 3 LED del modulo in un unico blob solido.
        cam.setDSP(2, 130, 150, 40); 
    }
    else if (sensitivity == Sensitivity_High) {
        // Guadagno intermedio (4X) - Compensa distanze superiori o ambienti più ampi.
        cam.setGain(0x10, 0x02); 
        
        // DSP adattato per tollerare l'ingrossamento del blob (blooming) dovuto all'amplificazione:
        // - area_max (200): Allargato per evitare che il sensore scarti i LED più luminosi/vicini.
        // - noise_th (50): Finestra di fusione più ampia per gestire la sfumatura estesa dei 3 LED.
        cam.setDSP(2, 130, 200, 50); 
    }
    else if (sensitivity == Sensitivity_Max) {
        // Guadagno massimo fisico (8X) - Per distanze estreme o condizioni critiche.
        cam.setGain(0x10, 0x03); 
        
        // DSP spinto ai limiti operativi del silicio:
        // - area_min (1): Permette di agganciare anche un singolo pixel debolissimo a grande distanza.
        // - brightness_th (150): Alzato per fare da diga contro il rumore termico amplificato 8 volte.
        // - area_max (300): Tolleranza massima per macchie di luce molto estese.
        // - noise_th (60): Compensa la forte dispersione ottica generata dall'alta amplificazione.
        cam.setDSP(1, 150, 300, 60); 
    }
    else if (sensitivity == Sensitivity_Default) {
        // Profilo standard di sicurezza (uguale a Min)
        cam.setGain(0x10, 0x00); 
        cam.setDSP(2, 130, 150, 40);       
    }
    else {
        // Fallback di sicurezza per valori numerici fuori dall'enum
        cam.setGain(0x10, 0x00); 
        cam.setDSP(2, 130, 150, 40);
    }
}

void DFRobotIRPositionEx::readAndUnpack(bool updateSeen) {
    PAJ7025_Object objs[4]; // 16
    
    // Mappa DataFormat_Basic a PAJ7025_FORMAT_2_96_BYTE per la massima velocità
    // Mappa DataFormat_Extended a PAJ7025_FORMAT_1_256_BYTE per avere le size (area)
    uint8_t paj_format = (current_format == DataFormat_Extended) ? PAJ7025_FORMAT_1_256_BYTE : PAJ7025_FORMAT_2_96_BYTE;
    
    int valid = cam.readData(objs, paj_format);  // legge 4 record
    
    unsigned int tempSeenFlags = 0;
    
    int validCount = 0;
    
    // Scandagliamo l'intero array da 16 oggetti della PAJ7025 per cercare i blob validi sparsi - ho impostato solo 4 per ottimizzare
    for(int i = 0; i < 4 && validCount < 4; i++) {
        if(objs[i].is_valid) {
            positionX[validCount] = objs[i].cx;
            positionY[validCount] = objs[i].cy;
            // Se formato extended, popola la size. (La vecchia DFRobot dava size da 0 a 15, la scaliamo o diamo l'area diretta)
            if(current_format == DataFormat_Extended) {
                unpackedSizes[validCount] = (objs[i].area > 15) ? 15 : objs[i].area; 
            } else {
                unpackedSizes[validCount] = 15; // default vuoto per compatibilità DFRobot
            }
            tempSeenFlags |= (1 << validCount);
            validCount++;
        }
    }
    
    // Riempi i rimanenti slot (da validCount a 3) con valori "vuoti fittizi"
    for(int i = validCount; i < 4; i++) {
        positionX[i] = 1023;
        positionY[i] = 1023;
        unpackedSizes[i] = 15;
    }
    
    if (updateSeen) {
        seenFlags = tempSeenFlags;
    }
}

void DFRobotIRPositionEx::requestPositionExtended() {
    // SPI è immediato, non c'è asincronismo come I2C, nulla da pre-richiedere
}

void DFRobotIRPositionEx::requestPositionBasic() {
    // Nulla da fare
}

bool DFRobotIRPositionEx::availableExtended() {
    readAndUnpack(true);
    return true; // SPI è sempre subito disponibile
}

bool DFRobotIRPositionEx::availableExtendedNoSeen() {
    readAndUnpack(false);
    return true;
}

bool DFRobotIRPositionEx::availableBasic() {
    readAndUnpack(true);
    return true;
}

bool DFRobotIRPositionEx::availableBasicNoSeen() {
    readAndUnpack(false);
    return true;
}

int DFRobotIRPositionEx::basicAtomic(DFRobotIRPositionEx::Retry_e retry) {
    // SPI non soffre di corruzione asincrona del buffer I2C. La lettura è sempre atomica.
    readAndUnpack(true);
    return Error_Success;
}

int DFRobotIRPositionEx::extendedAtomic(DFRobotIRPositionEx::Retry_e retries) {
    readAndUnpack(true);
    return Error_Success;
}


#endif //PAJ7025_CAM