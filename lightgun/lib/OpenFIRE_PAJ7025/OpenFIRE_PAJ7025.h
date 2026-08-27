#ifdef PAJ7025_CAM

/*!
 * @file OpenFIRE_PAJ7025.h
 * @brief Internal OpenFIRE backend adapter for the PixArt PAJ7025 camera.
 *
 * @copyright Alessandro Satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author Alessandro Satanassi
 * @version V1.0
 * @date 2026
 */

#ifndef OPENFIRE_PAJ7025_H
#define OPENFIRE_PAJ7025_H

#include <stdint.h>
#include <SPI.h>
#include <PixArt_PAJ7025.h>
#include <OpenFIRECameraProfile.h>

class OpenFIRE_PAJ7025 {
private:
    PAJ7025 cam;
    SPIClass* _spi;
    int8_t _csPin;
    const CameraProfile* _profile;

    int positionX[4];
    int positionY[4];
    int unpackedSizes[4];
    PAJ7025_Object extendedObjects[4];
    unsigned int seenFlags;

    int readBasicAndUnpack();
    int readExtendedAndUnpack();

public:
    OpenFIRE_PAJ7025(SPIClass* spiPort, int8_t csPin, const CameraProfile& profile);

    bool begin(uint32_t clock, uint8_t sensitivity);
    void sensitivityLevel(uint8_t sensitivity);

    int basicAtomic();
    int extendedAtomic();

    int x(int index) const { return positionX[index]; }
    int y(int index) const { return positionY[index]; }
    int size(int index) const { return unpackedSizes[index]; }
    const int* xPositions() const { return positionX; }
    const int* yPositions() const { return positionY; }
    const int* sizes() const { return unpackedSizes; }
    const PAJ7025_Object* objects() const { return extendedObjects; }
    unsigned int seen() const { return seenFlags; }
};

#endif // OPENFIRE_PAJ7025_H
#endif // PAJ7025_CAM
