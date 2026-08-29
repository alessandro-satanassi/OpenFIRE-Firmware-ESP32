/*!
 * @file OpenFIRECameraProfile.h
 * @brief Runtime camera model/profile definitions for OpenFIRE.
 *
 * @copyright Alessandro Satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author Alessandro Satanassi
 * @version V1.0
 * @date 2026
 */

#ifndef OPENFIRE_CAMERA_PROFILE_H
#define OPENFIRE_CAMERA_PROFILE_H

#include <stdint.h>
#include <OpenFIREConst.h>
#include "../../src/boards/OpenFIREshared.h"


using CameraModel = decltype(OF_Const::cameraModel_e);

// ============================================================================
// CAMERA PROFILE CONFIGURATION
//
// To add or configure a camera:
//   1. Add the camera model to OpenFIREshared.h.
//   2. Add or modify its profile in the "SUPPORTED CAMERA PROFILES" section.
//
// Camera-specific values are entered only in the profile definitions below.
// All derived values are calculated automatically by MakeProfile().
// ============================================================================

// ============================================================================
// DEFAULT CAMERA PROFILE GEOMETRY
//
// These values are used only to generate the initial/default IR geometry
// for each camera profile.
// ============================================================================

// 5% margin on each side -> 90% usable camera area.
constexpr float IR_SQUARE_USABLE_FACTOR = 0.90f;

// 5% margin on each side -> 90% usable camera area.
constexpr float IR_DIAMOND_USABLE_FACTOR = 0.90f;

// Default Diamond bounding-box width / height ratio follows the monitor.
constexpr float IR_DIAMOND_MONITOR_RATIO =
    (float)res_x / (float)res_y;

// ============================================================================
// CAMERA BUS CLOCKS
// Communication speed used by each camera backend.
// ============================================================================

constexpr uint32_t DFROBOT_I2C_CLOCK = 1000000U; // 1 MHz I2C
constexpr uint32_t PAJ7025_SPI_CLOCK = 2000000U; // 2 MHz SPI

// ============================================================================
// CAMERA PROFILE DATA
// ============================================================================

struct CameraProfile {
    // Camera model
    CameraModel model;

    // Camera coordinate space
    // Virtual/output coordinate resolution provided by the camera.
    int camResX;
    int camResY;
    int camMaxX;
    int camMaxY;

    // OpenFIRE internal coordinate space
    // The shift increases the coordinate precision used by OpenFIRE.
    // Internal coordinates = virtual/output coordinates << camToMouseShift.
    uint8_t camToMouseShift;
    int camToMouseMult;

    int mouseResX;
    int mouseResY;
    int mouseMaxX;
    int mouseMaxY;

    // Real sensor characteristics
    int sensorResX;
    int sensorResY;
    float noiseFactor;

    // Camera operating parameters
    uint16_t fps;
    uint32_t busClock;

    // Lens correction
    float lensRadialK1;
    float lensRadialK2;

    // Default Square geometry - automatically calculated
    int squareTLX;
    int squareTLY;
    int squareTRX;
    int squareTRY;
    int squareBLX;
    int squareBLY;
    int squareBRX;
    int squareBRY;

    // Default Diamond geometry - automatically calculated
    int diamondTCX;
    int diamondTCY;
    int diamondRCX;
    int diamondRCY;
    int diamondBCX;
    int diamondBCY;
    int diamondLCX;
    int diamondLCY;

    // RAW/Test display geometry - automatically calculated
    int testWidth;
    int testHeight;
    int testOffsetX;
    int testOffsetY;
};

namespace OpenFIRE_CameraProfiles {

// ============================================================================
// CAMERA PROFILE BUILDER
//
// Generates all derived values from the camera-specific parameters.
// MakeProfile() is constexpr, so profiles are calculated at compile time.
// Normally this section should not need to be modified when tuning a camera.
// ============================================================================

constexpr CameraProfile MakeProfile(CameraModel model,
                                    int camResX,
                                    int camResY,
                                    uint8_t camToMouseShift,
                                    int sensorResX,
                                    int sensorResY,
                                    float noiseFactor,
                                    uint16_t fps,
                                    uint32_t busClock,
                                    float lensRadialK1,
                                    float lensRadialK2) {
    // OpenFIRE internal coordinate space
    const int camToMouseMult = 1 << camToMouseShift;
    const int mouseResX = camResX * camToMouseMult;
    const int mouseResY = camResY * camToMouseMult;

    // Default Square geometry
    const float squareUsableWidth = (float)camResX * IR_SQUARE_USABLE_FACTOR;
    const float squareUsableHeight = (float)camResY * IR_SQUARE_USABLE_FACTOR;
    const float squareHeightFromWidth = squareUsableWidth / IR_SQUARE_BASE_HEIGHT_RATIO;
    const float squareHeight = (squareHeightFromWidth < squareUsableHeight) ? squareHeightFromWidth : squareUsableHeight;
    const float squareWidth = squareHeight * IR_SQUARE_BASE_HEIGHT_RATIO;
    const float squareCenterX = (float)camResX * 0.5f;
    const float squareCenterY = (float)camResY * 0.5f;

    const int squareTLX = (int)(squareCenterX - squareWidth * 0.5f);
    const int squareTLY = (int)(squareCenterY - squareHeight * 0.5f);
    const int squareTRX = (int)(squareCenterX + squareWidth * 0.5f);
    const int squareTRY = squareTLY;
    const int squareBLX = squareTLX;
    const int squareBLY = (int)(squareCenterY + squareHeight * 0.5f);
    const int squareBRX = squareTRX;
    const int squareBRY = squareBLY;

    // Default Diamond geometry
    const float diamondUsableWidth = (float)camResX * IR_DIAMOND_USABLE_FACTOR;
    const float diamondUsableHeight = (float)camResY * IR_DIAMOND_USABLE_FACTOR;
    const float diamondHeightFromWidth = diamondUsableWidth / IR_DIAMOND_MONITOR_RATIO;
    const float diamondHeight = (diamondHeightFromWidth < diamondUsableHeight) ? diamondHeightFromWidth : diamondUsableHeight;
    const float diamondWidth = diamondHeight * IR_DIAMOND_MONITOR_RATIO;
    const float diamondCenterX = (float)camResX * 0.5f;
    const float diamondCenterY = (float)camResY * 0.5f;

    const int diamondTCX = (int)diamondCenterX;
    const int diamondTCY = (int)(diamondCenterY - diamondHeight * 0.5f);
    const int diamondRCX = (int)(diamondCenterX + diamondWidth * 0.5f);
    const int diamondRCY = (int)diamondCenterY;
    const int diamondBCX = (int)diamondCenterX;
    const int diamondBCY = (int)(diamondCenterY + diamondHeight * 0.5f);
    const int diamondLCX = (int)(diamondCenterX - diamondWidth * 0.5f);
    const int diamondLCY = (int)diamondCenterY;

    // RAW/Test display geometry
    constexpr int screenResX = res_x >> 2;
    constexpr int screenResY = res_y >> 2;
    const bool widthLimited = (mouseResX * screenResY) > (mouseResY * screenResX);
    const int testWidth = widthLimited ? screenResX : (mouseResX * screenResY) / mouseResY;
    const int testHeight = widthLimited ? (mouseResY * screenResX) / mouseResX : screenResY;

    return {
        model,

        // Camera coordinate space
        camResX,
        camResY,
        camResX - 1,
        camResY - 1,

        // OpenFIRE internal coordinate space
        camToMouseShift,
        camToMouseMult,
        mouseResX,
        mouseResY,
        mouseResX - 1,
        mouseResY - 1,

        // Real sensor characteristics
        sensorResX,
        sensorResY,
        noiseFactor,

        // Camera operating parameters
        fps,
        busClock,

        // Lens correction
        lensRadialK1,
        lensRadialK2,

        // Default Square geometry
        squareTLX,
        squareTLY,
        squareTRX,
        squareTRY,
        squareBLX,
        squareBLY,
        squareBRX,
        squareBRY,

        // Default Diamond geometry
        diamondTCX,
        diamondTCY,
        diamondRCX,
        diamondRCY,
        diamondBCX,
        diamondBCY,
        diamondLCX,
        diamondLCY,

        // RAW/Test display geometry
        testWidth,
        testHeight,
        (screenResX - testWidth) / 2,
        (screenResY - testHeight) / 2
    };
}

// ============================================================================
// SUPPORTED CAMERA PROFILES
//
// THIS IS THE MAIN SECTION TO EDIT WHEN ADDING OR TUNING A CAMERA.
//
// Parameters:
//   model               Camera model defined in OpenFIREshared.h
//   camResX/Y           Virtual/output coordinate resolution
//   camToMouseShift     Additional OpenFIRE precision shift
//   sensorResX/Y        Real sensor resolution
//   noiseFactor         Sensor noise factor used by tracking/filtering
//   fps                 Camera frame rate
//   busClock            I2C/SPI communication clock
//   lensRadialK1/K2     Radial lens correction coefficients
// ============================================================================

// ---------------------------------------------------------------------------
// DFRobot SEN0158 / WiiCam
// ---------------------------------------------------------------------------
static constexpr CameraProfile DFRobot_SEN0158 = MakeProfile(
    OF_Const::DFRobot_SEN0158,    // Camera model
    1024,                         // Virtual/output coordinate resolution X
    768,                          // Virtual/output coordinate resolution Y
    2,                            // OpenFIRE internal precision shift (virtual coordinates << shift)
    128,                          // Real sensor resolution X
    96,                           // Real sensor resolution Y
    1.0f,                         // Sensor noise factor
    209,                          // Camera FPS
    DFROBOT_I2C_CLOCK,            // Camera bus clock
    0.006f,                       // Lens radial correction K1
    0.0f                          // Lens radial correction K2
);

// ---------------------------------------------------------------------------
// PixArt PAJ7025R2
// ---------------------------------------------------------------------------
static constexpr CameraProfile PixArt_PAJ7025R2 = MakeProfile(
    OF_Const::PixArt_PAJ7025R2,   // Camera model
    4096,                         // Virtual/output coordinate resolution X
    4096,                         // Virtual/output coordinate resolution Y
    0,                            // OpenFIRE internal precision shift (virtual coordinates << shift)
    98,                           // Real sensor resolution X
    98,                           // Real sensor resolution Y
    1.0f,                         // Sensor noise factor
    209,                          // Camera FPS
    PAJ7025_SPI_CLOCK,            // Camera bus clock
    0.0f,                         // Lens radial correction K1
    0.0f                          // Lens radial correction K2
);

// ---------------------------------------------------------------------------
// PixArt PAJ7025R3
// ---------------------------------------------------------------------------
static constexpr CameraProfile PixArt_PAJ7025R3 = MakeProfile(
    OF_Const::PixArt_PAJ7025R3,   // Camera model
    4096,                         // Virtual/output coordinate resolution X
    4096,                         // Virtual/output coordinate resolution Y
    0,                            // OpenFIRE internal precision shift (virtual coordinates << shift)
    98,                           // Real sensor resolution X
    98,                           // Real sensor resolution Y
    1.0f,                         // Sensor noise factor
    209,                          // Camera FPS
    PAJ7025_SPI_CLOCK,            // Camera bus clock
    0.0f,                         // Lens radial correction K1
    0.0f                          // Lens radial correction K2
);

/*
constexpr const CameraProfile* Find(CameraModel model)
{
    switch (model) {
        case OF_Const::DFRobot_SEN0158:
            return &DFRobot_SEN0158;

        case OF_Const::PixArt_PAJ7025R2:
            return &PixArt_PAJ7025R2;

        case OF_Const::PixArt_PAJ7025R3:
            return &PixArt_PAJ7025R3;

        default:
            return nullptr;
    }
}
*/
} // namespace OpenFIRE_CameraProfiles

#endif // OPENFIRE_CAMERA_PROFILE_H
