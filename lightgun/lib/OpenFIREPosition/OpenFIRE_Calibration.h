/*!
 * @file OpenFIRE_Calibration.h
 * @brief Six-shot calibration; no work or retained sample buffer in Run mode.
 * @copyright Alessandro Satanassi, 2026, GNU Lesser General Public License
 */
#ifndef OPENFIRE_CALIBRATION_H
#define OPENFIRE_CALIBRATION_H

#include <array>
#include <stdint.h>
#include "OpenFIRECameraProfile.h"

namespace OpenFIRECalibration {

// Existing profile representation: top, bottom, left and right offsets.
struct Parameters {
    float sourceX = 0, sourceY = 0;
    float ledLeft = 0, ledRight = 0;
    int32_t offsets[4] = {};
};

// Mean LED geometry before one trigger pull (tracker coordinates).
struct Measurement {
    std::array<float, 8> mean{};
    bool valid = false;
};

// Pure calibration solve. Failure leaves result exactly equal to original.
bool Refine(const CameraProfile& profile, bool diamond,
            const std::array<Measurement, 6>& shots,
            const Parameters& original, Parameters& result);

// Allocate this session on the heap, only while calibrating (about 10 KiB).
class Session {
public:
    Session(const CameraProfile& profile, bool diamond);
    void reset();
    void invalidate();
    void observe(const int32_t* xy, uint32_t nowMicros);
    uint32_t frameCount() const { return frames; }
    bool capture(uint8_t shot, uint32_t nowMicros);
    bool refine(const Parameters& original, Parameters& result) const;

private:
    CameraProfile profile;
    bool diamond;
    uint32_t maxGap, maxAge;
    uint32_t lastStamp = 0;
    uint32_t frames = 0;
    uint16_t next = 0, count = 0;
    int32_t history[256][8];
    uint32_t stamps[256];
    std::array<Measurement, 6> shots{};
};

} // namespace OpenFIRECalibration
#endif
