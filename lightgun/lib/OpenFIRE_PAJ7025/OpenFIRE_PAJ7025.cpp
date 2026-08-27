#ifdef PAJ7025_CAM

/*!
 * @file OpenFIRE_PAJ7025.cpp
 * @brief Internal OpenFIRE backend adapter for the PixArt PAJ7025 camera.
 *
 * @copyright Alessandro Satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author Alessandro Satanassi
 * @version V1.0
 * @date 2026
 */

#include <string.h>
#include "OpenFIRE_PAJ7025.h"

OpenFIRE_PAJ7025::OpenFIRE_PAJ7025(SPIClass* spiPort, int8_t csPin, const CameraProfile& profile)
    : _spi(spiPort), _csPin(csPin), _profile(&profile), seenFlags(0) {
    memset(extendedObjects, 0, sizeof(extendedObjects));
}

bool OpenFIRE_PAJ7025::begin(uint32_t clock, uint8_t sensitivity) {
    if (_spi == nullptr || _csPin < 0 || _profile == nullptr) {
        return false;
    }

    if (!cam.begin(_spi, (uint8_t)_csPin, clock)) {
        return false;
    }

    cam.setFrameRate(_profile->fps);
    cam.setExposure(300);
    sensitivityLevel(sensitivity);
    cam.setResolution((uint16_t)_profile->camMaxX, (uint16_t)_profile->camMaxY);

    return true;
}

void OpenFIRE_PAJ7025::sensitivityLevel(uint8_t sensitivity) {
    // The public OpenFIRECamera API clamps sensitivity to 0..2 before this
    // internal backend is called.
    if (sensitivity == 0U) {
        cam.setGain(0x10, 0x00);
        cam.setDSP(2, 130, 150, 40);
    }
    else if (sensitivity == 1U) {
        cam.setGain(0x10, 0x02);
        cam.setDSP(2, 130, 200, 50);
    }
    else {
        cam.setGain(0x10, 0x03);
        cam.setDSP(1, 150, 300, 60);
    }
}

int OpenFIRE_PAJ7025::readBasicAndUnpack() {
    PAJ7025_Object rawObjects[4];
    cam.readData(rawObjects, PAJ7025_FORMAT_BASIC);

    int validCount = 0;

    for (int i = 0; i < 4; i++) {
        if (rawObjects[i].is_valid) {
            positionX[validCount] = rawObjects[i].cx;
            positionY[validCount] = rawObjects[i].cy;
            unpackedSizes[validCount] = 15;
            validCount++;
        }
    }

    for (int i = validCount; i < 4; i++) {
        positionX[i] = _profile->camResX / 2;
        positionY[i] = _profile->camResY / 2;
        unpackedSizes[i] = 15;
    }

    seenFlags = (1U << validCount) - 1U;
    return 0;
}

int OpenFIRE_PAJ7025::readExtendedAndUnpack() {
    PAJ7025_Object rawObjects[4];
    cam.readData(rawObjects, PAJ7025_FORMAT_EXTENDED);

    memset(extendedObjects, 0, sizeof(extendedObjects));
    int validCount = 0;

    for (int i = 0; i < 4; i++) {
        if (rawObjects[i].is_valid) {
            positionX[validCount] = rawObjects[i].cx;
            positionY[validCount] = rawObjects[i].cy;
            unpackedSizes[validCount] = (rawObjects[i].area > 15U) ? 15 : (int)rawObjects[i].area;
            extendedObjects[validCount] = rawObjects[i];
            validCount++;
        }
    }

    for (int i = validCount; i < 4; i++) {
        positionX[i] = _profile->camResX / 2;
        positionY[i] = _profile->camResY / 2;
        unpackedSizes[i] = 15;
    }

    seenFlags = (1U << validCount) - 1U;
    return 0;
}

int OpenFIRE_PAJ7025::basicAtomic() {
    return readBasicAndUnpack();
}

int OpenFIRE_PAJ7025::extendedAtomic() {
    return readExtendedAndUnpack();
}

#endif // PAJ7025_CAM
