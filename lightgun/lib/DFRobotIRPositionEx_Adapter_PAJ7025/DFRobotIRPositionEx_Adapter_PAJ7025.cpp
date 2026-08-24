#ifdef PAJ7025_CAM

/*!
 * @file DFRobotIRPositionEx_Adapter_PAJ7025.cpp
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

#include "DFRobotIRPositionEx_Adapter_PAJ7025.h"
#include "OpenFIREConst.h"

DFRobotIRPositionEx::DFRobotIRPositionEx(SPIClass* spiPort, int8_t csPin) : _spi(spiPort), _csPin(csPin) {
    current_format = DataFormat_Basic;
    seenFlags = 0;
}

bool DFRobotIRPositionEx::begin(uint32_t clock, DataFormat_e format, Sensitivity_e sensitivity) {
    if (_spi == nullptr || _csPin < 0) {
        return false;
    }

    if (!cam.begin(_spi, (uint8_t)_csPin, clock)) {
        return false;
    }

    cam.setFrameRate(CamFPS);
    cam.setExposure(300);
    sensitivityLevel(sensitivity);
    cam.setResolution(CamMaxX, CamMaxY);
    dataFormat(format);

    return true;
}

void DFRobotIRPositionEx::dataFormat(DataFormat_e format) {
    current_format = format;
}

void DFRobotIRPositionEx::sensitivityLevel(Sensitivity_e sensitivity) {
    if (sensitivity > Sensitivity_Max) {
        sensitivity = Sensitivity_Max;
    }

    if (sensitivity == Sensitivity_Min) {
        cam.setGain(0x10, 0x00);
        cam.setDSP(2, 130, 150, 40);
    }
    else if (sensitivity == Sensitivity_High) {
        cam.setGain(0x10, 0x02);
        cam.setDSP(2, 130, 200, 50);
    }
    else if (sensitivity == Sensitivity_Max) {
        cam.setGain(0x10, 0x03);
        cam.setDSP(1, 150, 300, 60);
    }
    else {
        cam.setGain(0x10, 0x00);
        cam.setDSP(2, 130, 150, 40);
    }
}

void DFRobotIRPositionEx::readAndUnpack(bool updateSeen) {
    PAJ7025_Object objs[4];

    uint8_t paj_format = (current_format == DataFormat_Extended) ? PAJ7025_FORMAT_EXTENDED : PAJ7025_FORMAT_BASIC;

    cam.readData(objs, paj_format);

    int validCount = 0;

    for(int i = 0; i < 4; i++) {
        if(objs[i].is_valid) {
            positionX[validCount] = objs[i].cx;
            positionY[validCount] = objs[i].cy;
            if(current_format == DataFormat_Extended) {
                unpackedSizes[validCount] = (objs[i].area > 15) ? 15 : objs[i].area;
            } else {
                unpackedSizes[validCount] = 15;
            }
            validCount++;
        }
    }

    for(int i = validCount; i < 4; i++) {
        positionX[i] = CamResX / 2;
        positionY[i] = CamResY / 2;
        unpackedSizes[i] = 15;
    }

    if (updateSeen) {
        seenFlags = (1U << validCount) - 1U;
    }
}

void DFRobotIRPositionEx::requestPositionExtended() {
}

void DFRobotIRPositionEx::requestPositionBasic() {
}

bool DFRobotIRPositionEx::availableExtended() {
    readAndUnpack(true);
    return true;
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
    (void)retry;
    readAndUnpack(true);
    return Error_Success;
}

int DFRobotIRPositionEx::extendedAtomic(DFRobotIRPositionEx::Retry_e retries) {
    (void)retries;
    readAndUnpack(true);
    return Error_Success;
}

#endif
