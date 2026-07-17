#ifdef PAJ7025_CAM


#include "DFRobotIRPositionEx_Wrapper.h"

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
    if (!cam.begin(_spi, (uint8_t)_csPin)) {
        return false;
    }
    
    // Imposta la risoluzione massima 4095x4095 come richiesto per il calcolo mouse
    cam.setResolution(4095, 4095);
    
    // Converti la sensitivity nel gain/dsp della PAJ7025
    sensitivityLevel(sensitivity);
    
    // Imposta il formato
    dataFormat(format);

    return true;
}

void DFRobotIRPositionEx::dataFormat(DataFormat_e format) {
    current_format = format;
}

void DFRobotIRPositionEx::sensitivityLevel(Sensitivity_e sensitivity) {
    // Adatta il gain e il DSP in base alla sensibilità desiderata.
    // L'ordine degli if garantisce che le impostazioni specifiche (Min, High, Max)
    // abbiano sempre la priorità se condividono lo stesso valore numerico di Default.

    if (sensitivity == Sensitivity_Min) {
        cam.setGain(5, 0);
        cam.setDSP(5, 0xA0); // Area minima alta, luminosità richiesta alta
    }
    else if (sensitivity == Sensitivity_High) {
        cam.setGain(20, 0);
        cam.setDSP(2, 0x70);
    }
    else if (sensitivity == Sensitivity_Max) {
        cam.setGain(31, 0);  // Max Gain
        cam.setDSP(1, 0x50); // Prendi anche i riflessi più flebili
    }
    // Viene valutato SOLO SE il valore non è già stato intercettato da Min, High o Max
    else if (sensitivity == Sensitivity_Default) {
        cam.setGain(15, 0);
        cam.setDSP(3, 0x97); // Valori bilanciati per un Default "unico"
    }
    else {
        // Fallback di sicurezza per valori numerici fuori dall'enum
        cam.setGain(15, 0);
        cam.setDSP(3, 0x97);
    }
}


void DFRobotIRPositionEx::readAndUnpack(bool updateSeen) {
    PAJ7025_Object objs[16];
    
    // Mappa DataFormat_Basic a PAJ7025_FORMAT_2_96_BYTE per la massima velocità
    // Mappa DataFormat_Extended a PAJ7025_FORMAT_1_256_BYTE per avere le size (area)
    uint8_t paj_format = (current_format == DataFormat_Extended) ? PAJ7025_FORMAT_1_256_BYTE : PAJ7025_FORMAT_2_96_BYTE;
    
    int valid = cam.readData(objs, paj_format);
    
    unsigned int tempSeenFlags = 0;
    
    int validCount = 0;
    
    // Scandagliamo l'intero array da 16 oggetti della PAJ7025 per cercare i blob validi sparsi
    for(int i = 0; i < 16 && validCount < 4; i++) {
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