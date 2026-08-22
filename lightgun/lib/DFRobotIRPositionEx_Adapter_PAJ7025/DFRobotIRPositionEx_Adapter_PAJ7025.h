#ifdef PAJ7025_CAM

/*!
 * @file DFRobotIRPositionEx_Wrapper.h
 * @brief Wrapper PAJ7025 -> DFRobotIRPositionEx
 * @details Questa libreria offre la stessa interfaccia della vecchia libreria DFRobot
 * ma usa internamente il sensore PAJ7025 su bus SPI (4 punti configurati, 4095x4095).
 */

#ifndef DFRobotIRPositionEx_Adapter_PAJ7025_h
#define DFRobotIRPositionEx_Adapter_PAJ7025_h

#include <stdint.h>
#include <SPI.h>
#include <PixArt_PAJ7025.h> // Richiede l'installazione della nuova libreria PAJ7025

class DFRobotIRPositionEx {
private:
    PAJ7025 cam;
    SPIClass* _spi;
    int8_t _csPin;
    uint8_t current_format;

    int positionX[4];
    int positionY[4];
    int unpackedSizes[4];
    unsigned int seenFlags;

    void readAndUnpack(bool updateSeen);

public:
    /*!
    * @brief Data format
    */
    enum DataFormat_e {
        DataFormat_Basic = 0,       ///< Usa PAJ7025 Format 2 (ultra veloce)
        DataFormat_Extended = 1     ///< Usa PAJ7025 Format 1 (per leggere anche le size)
    };

    /*!
    * @brief Camera sensitivity.
    */
    enum Sensitivity_e {
        Sensitivity_Min = 0,
        Sensitivity_Default = 0,
        Sensitivity_High = 1,
        Sensitivity_Max = 2 
    };

    /*!
    * @brief Error codes.
    */
    enum Errors_e {
        Error_SuccessMismatch = 1,  
        Error_Success = 0,        
        Error_IICerror = -1,      
        Error_DataMismatch = -2,  
    };

    /*!
    * @brief Retry options (Nel PAJ7025 SPI non ci sono quasi mai errori, ma manteniamo l'enum per compatibilità)
    */
    enum Retry_e {
        Retry_0 = 0,    
        Retry_0s = 1,   
        Retry_1 = 2,    
        Retry_1s = 3,   
        Retry_2 = 4,    
        Retry_2s = 5    
    };
    
    /*!
    * @brief Costruttore SPI (Sostituisce il vecchio TwoWire&)
    * @param spiPort Porta SPI (&SPI, &SPI1)
    * @param csPin Pin Chip Select. Deve essere sempre specificato.
    */
    DFRobotIRPositionEx(SPIClass* spiPort, int8_t csPin);
  
    ~DFRobotIRPositionEx();
  
    /*!
    * @brief Inizializza il sensore PAJ7025. Imposta la risoluzione massima a 4095x4095.
    * @param clock Mantenuto per compatibilità con DFRobotIRPositionEx; il PAJ7025 usa SPI a 2 MHz.
    */
    bool begin(uint32_t clock = 2000000 /*2Mhz*/, DataFormat_e format = DataFormat_Basic, Sensitivity_e sensitivity = Sensitivity_Default);

    void dataFormat(DataFormat_e format);

    void sensitivityLevel(Sensitivity_e sensitivity);

    void requestPositionExtended();
    void requestPositionBasic();

    bool availableExtended();
    bool availableExtendedNoSeen();
    bool availableBasic();
    bool availableBasicNoSeen();

    int basicAtomic(DFRobotIRPositionEx::Retry_e retry = DFRobotIRPositionEx::Retry_1s);
    int extendedAtomic(DFRobotIRPositionEx::Retry_e retries = DFRobotIRPositionEx::Retry_1s);

    int x(int index) const { return positionX[index]; }
    int y(int index) const { return positionY[index]; }
    int size(int index) const { return unpackedSizes[index]; }

    const int* xPositions() const { return positionX; }
    const int* yPositions() const { return positionY; }
    const int* sizes() const { return unpackedSizes; }

    unsigned int seen() const { return seenFlags; }
};

#endif // DFRobotIRPositionEx_Adapter_PAJ7025_h

#endif //PAJ7025_CAM