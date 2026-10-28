/*!
 * @file OpenFIRE_Calibration.h
 * @brief Calibration from the target shots; no work or retained sample buffer in Run mode.
 * @copyright Alessandro Satanassi, 2026, GNU Lesser General Public License
 */
#ifndef OPENFIRE_CALIBRATION_H
#define OPENFIRE_CALIBRATION_H

#include <array>
#include <stdint.h>
#include "OpenFIRECameraProfile.h"

namespace OpenFIRECalibration {

// Shots in calibration order: centre, top, bottom, left, right, centre, then the optional
// corners top-right, bottom-right, bottom-left, top-left. A shot not taken stays invalid.
constexpr unsigned ShotSlots = 10;

// LEDs not exactly on a rectangle aligned with the screen (skew and keystone), measured with
// the four corners. All zero = no correction, as in profiles calibrated without the corners.
struct Shape {
    float skewX = 0, skewY = 0, keystoneX = 0, keystoneY = 0;
};

// Existing profile representation: top, bottom, left and right offsets, then the LED shape.
struct Parameters {
    float sourceX = 0, sourceY = 0;
    float ledLeft = 0, ledRight = 0;
    int32_t offsets[4] = {};
    Shape shape;
};

// Mean LED geometry before one trigger pull (tracker coordinates).
struct Measurement {
    std::array<float, 8> mean{};
    bool valid = false;
};

// Pure calibration solve. Failure leaves result exactly equal to original.
bool Refine(const CameraProfile& profile, bool diamond,
            const std::array<Measurement, ShotSlots>& shots,
            const Parameters& original, Parameters& result);

// Run mode: Perspective output (x, y) corrected for the LED shape, before the edge offsets.
// The result stays within +-limit; nothing changes with a zero shape.
void ApplyShape(const Shape& shape, int32_t limit, int32_t& x, int32_t& y);

// Run mode: Perspective output (x, y) corrected for the camera mounted below the sight line.
// offset = that distance / vertical distance between the LEDs. A shift fixed in the gun moves
// the aim by offset * (D / D0 - 1) LED heights away from the calibration distance D0, with
// D / D0 = sqrt(calibrationArea / area) (OpenFIRE_Perspective::area(), same units for both).
// angle = tilt of the gun seen by the tracker (Ang(), radians): the shift turns with the gun.
// span = LED height in Perspective units. Nothing changes with a zero offset or area.
void ApplyParallax(float offset, float calibrationArea, float area, float angle, int32_t span, int32_t limit,
                   int32_t& x, int32_t& y);

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
    std::array<Measurement, ShotSlots> shots{};
};

} // namespace OpenFIRECalibration
#endif
