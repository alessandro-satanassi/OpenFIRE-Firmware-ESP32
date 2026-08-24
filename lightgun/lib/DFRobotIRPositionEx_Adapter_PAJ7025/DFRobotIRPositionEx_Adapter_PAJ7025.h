#ifdef PAJ7025_CAM

/*!
 * @file DFRobotIRPositionEx_Adapter_PAJ7025.h
 * @brief Adapter for using the PixArt PAJ7025 camera with OpenFIRE
 * @n Provides a DFRobotIRPositionEx-compatible interface for the PixArt PAJ7025 driver
 *
 * @copyright Alessandro Satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author Alessandro Satanassi
 * @version V1.0
 * @date 2026
 */

#ifndef DFRobotIRPositionEx_Adapter_PAJ7025_h
#define DFRobotIRPositionEx_Adapter_PAJ7025_h

#include <stdint.h>
#include <SPI.h>
#include <PixArt_PAJ7025.h>

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
    enum DataFormat_e {
        DataFormat_Basic = 0,
        DataFormat_Extended = 1
    };

    enum Sensitivity_e {
        Sensitivity_Min = 0,
        Sensitivity_Default = 0,
        Sensitivity_High = 1,
        Sensitivity_Max = 2
    };

    enum Errors_e {
        Error_SuccessMismatch = 1,
        Error_Success = 0,
        Error_IICerror = -1,
        Error_DataMismatch = -2,
    };

    enum Retry_e {
        Retry_0 = 0,
        Retry_0s = 1,
        Retry_1 = 2,
        Retry_1s = 3,
        Retry_2 = 4,
        Retry_2s = 5
    };

    DFRobotIRPositionEx(SPIClass* spiPort, int8_t csPin);

    bool begin(uint32_t clock = 2000000, DataFormat_e format = DataFormat_Basic, Sensitivity_e sensitivity = Sensitivity_Default);

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

#endif

#endif
