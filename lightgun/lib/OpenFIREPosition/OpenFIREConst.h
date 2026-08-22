/*!
 * @file OpenFIREConst.h
 * @brief Global constants for the OpenFIRE light gun.
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V2.0
 * @date 2026
 *
 * I thank you for producing the first original code:
 *
 * @copyright Mike Lynch, 2021
 * @copyright GNU Lesser General Public License
 *
 * @author Mike Lynch
 * @version V1.0
 * @date 2021
 */

#ifndef _OPENFIRECONST_H_
#define _OPENFIRECONST_H_

#include <stdint.h>

// ============================================================================
// CAMERA SELECTION VALIDATION
// ============================================================================

// Exactly one physical camera model must be selected.
#if (defined(CAMERA_DFROBOT_SEN0158) + \
     defined(CAMERA_PIXART_PAJ7025R2) + \
     defined(CAMERA_PIXART_PAJ7025R3)) != 1
    #error "Select exactly one camera model"
#endif

// PAJ7025_CAM identifies the whole PAJ7025 family.
// A PAJ7025 build must also identify the specific R2/R3 model.
#if defined(PAJ7025_CAM) && \
    !defined(CAMERA_PIXART_PAJ7025R2) && \
    !defined(CAMERA_PIXART_PAJ7025R3)
    #error "PAJ7025_CAM requires a PAJ7025 camera model"
#endif

// Conversely, a PAJ7025 model must always belong to the PAJ7025_CAM family.
#if (defined(CAMERA_PIXART_PAJ7025R2) || defined(CAMERA_PIXART_PAJ7025R3)) && \
    !defined(PAJ7025_CAM)
    #error "PAJ7025 camera models require PAJ7025_CAM"
#endif

#ifdef PAJ7025_CAM
    // Configurazione hardware SPI // clock max 14 Mhz
    #define PAJ7025_SPI_CLOCK_2MHZ   2000000
    #define PAJ7025_SPI_CLOCK_4MHZ   4000000
    #define PAJ7025_SPI_CLOCK_8MHZ   8000000
    #define PAJ7025_SPI_CLOCK_14MHZ 14000000  // MAX

    // PixArt PAJ7025 IR positioning camera resolution
    constexpr int CamResX = 4096; //1024; // 2940; //1024; //4096;
    constexpr int CamResY = 4096; //768; // 2940; // 768;  //4096;

    // shift amount for extra precision for the maths
    // since the median is an average of 4 values, use 2 more bits
    constexpr int CamToMouseShift = 0; //2; // 0; // l'uscita è già abbastanza alta non serve aumentarla

    constexpr int CamSensorResX = 98;
    constexpr int CamSensorResY = 98;

    //constexpr float CamNoiseFactor = 1.0f; // 0.765625f; // 1.0f; // 0.765625f;
    constexpr int CamFPS = 209;

    constexpr uint32_t IR_CAMERA_BUS_CLOCK = PAJ7025_SPI_CLOCK_2MHZ;

    #if defined(CAMERA_PIXART_PAJ7025R2)
        constexpr float CamNoiseFactor = 1.0f; // 0.765625f; // 1.0f; // 0.765625f;
        constexpr float CamLensRadialK1 = 0.0f; // 0.006f;
        constexpr float CamLensRadialK2 = 0.0f;
    #elif defined(CAMERA_PIXART_PAJ7025R3)
        constexpr float CamNoiseFactor = 1.0f; // 0.765625f; // 1.0f; // 0.765625f;
        constexpr float CamLensRadialK1 = 0.0f; // 0.006f;
        constexpr float CamLensRadialK2 = 0.0f;
    #endif // CAMERA_PIXART_PAJ7025R2 / CAMERA_PIXART_PAJ7025R3

#elif CAMERA_DFROBOT_SEN0158 // DF ROBOT / WII_CAM
    // DFRobot IR camera IIC clock
    // even with cheap clips and the full length IR cam cable 1MHz is fine
    #define DFROBOT_IR_IIC_CLOCK 1000000

    // DFRobot IR positioning camera resolution
    constexpr int CamResX = 1024;
    constexpr int CamResY = 768;

    // shift amount for extra precision for the maths
    // since the median is an average of 4 values, use 2 more bits
    constexpr int CamToMouseShift = 2;

    constexpr int CamSensorResX = 128;
    constexpr int CamSensorResY = 96;

    constexpr float CamNoiseFactor = 1.0f;
    constexpr int CamFPS = 209;

    constexpr uint32_t IR_CAMERA_BUS_CLOCK = DFROBOT_IR_IIC_CLOCK;

    constexpr float CamLensRadialK1 = 0.006f;
    constexpr float CamLensRadialK2 = 0.0f;

#endif // PAJ7025_CAM // DF ROBOT / WII_CAM

// IR positioning camera maximum X and Y
constexpr int CamMaxX = CamResX - 1;
constexpr int CamMaxY = CamResY - 1;

// multiplier to convert IR camera position to mouse position
constexpr int CamToMouseMult = 1 << CamToMouseShift;

// mouse resolution
constexpr int MouseResX = CamResX * CamToMouseMult;
constexpr int MouseResY = CamResY * CamToMouseMult;

// Mouse position maximum X and Y
constexpr int MouseMaxX = MouseResX - 1;
constexpr int MouseMaxY = MouseResY - 1;

// Perspective code stuff
constexpr int res_x = 1920 << 2;
constexpr int res_y = 1080 << 2;

/////////////////////////////////////////////////////////////////////////////////////////////


// ============================================================================
// DEFAULT IR GEOMETRY
// ============================================================================
//
// SQUARE:
//   A = TL (Top Left)
//   B = TR (Top Right)
//   C = BL (Bottom Left)
//   D = BR (Bottom Right)
//
//   Physical base / height ratio = 0.711
//
// DIAMOND:
//   A = TC (Top Center)
//   B = RC (Right Center)
//   C = BC (Bottom Center)
//   D = LC (Left Center)
//
//   Physical bounding-box width / height ratio = monitor aspect ratio
//
// Both geometries reserve 5% of the camera FOV on every side.
// Therefore, 90% of the camera width and height is usable.
//
// All geometry calculations are constexpr and are resolved at compile time.
// ============================================================================


// ============================================================================
// SQUARE
// ============================================================================

// 5% margin on each side -> 90% usable camera area.
constexpr float IR_SQUARE_USABLE_FACTOR = 0.90f;

// Physical LED rectangle base / height ratio.
// Derived from the Qt alignment geometry.
constexpr float IR_SQUARE_BASE_HEIGHT_RATIO = 0.711f;

// Usable camera dimensions.
constexpr float IR_SQUARE_USABLE_WIDTH = (float)CamResX * IR_SQUARE_USABLE_FACTOR;

constexpr float IR_SQUARE_USABLE_HEIGHT = (float)CamResY * IR_SQUARE_USABLE_FACTOR;

// Maximum possible rectangle height allowed by camera width.
constexpr float IR_SQUARE_HEIGHT_FROM_WIDTH = IR_SQUARE_USABLE_WIDTH / IR_SQUARE_BASE_HEIGHT_RATIO;

// Select the limiting dimension.
constexpr float IR_SQUARE_HEIGHT = (IR_SQUARE_HEIGHT_FROM_WIDTH < IR_SQUARE_USABLE_HEIGHT) ? IR_SQUARE_HEIGHT_FROM_WIDTH : IR_SQUARE_USABLE_HEIGHT;

// Preserve physical 0.711 ratio.
constexpr float IR_SQUARE_WIDTH = IR_SQUARE_HEIGHT * IR_SQUARE_BASE_HEIGHT_RATIO;

// Camera center.
constexpr float IR_SQUARE_CENTER_X = (float)CamResX * 0.5f;

constexpr float IR_SQUARE_CENTER_Y = (float)CamResY * 0.5f;

// --------------------------------------------------------------------------
// SQUARE LED DEFAULT COORDINATES
// --------------------------------------------------------------------------

constexpr int IR_LED_SQUARE_TL_X = (int)(IR_SQUARE_CENTER_X - IR_SQUARE_WIDTH * 0.5f);
constexpr int IR_LED_SQUARE_TL_Y = (int)(IR_SQUARE_CENTER_Y - IR_SQUARE_HEIGHT * 0.5f);
constexpr int IR_LED_SQUARE_TR_X = (int)(IR_SQUARE_CENTER_X + IR_SQUARE_WIDTH * 0.5f);
constexpr int IR_LED_SQUARE_TR_Y = IR_LED_SQUARE_TL_Y;
constexpr int IR_LED_SQUARE_BL_X = IR_LED_SQUARE_TL_X;
constexpr int IR_LED_SQUARE_BL_Y = (int)(IR_SQUARE_CENTER_Y + IR_SQUARE_HEIGHT * 0.5f);
constexpr int IR_LED_SQUARE_BR_X = IR_LED_SQUARE_TR_X;
constexpr int IR_LED_SQUARE_BR_Y = IR_LED_SQUARE_BL_Y;

// ============================================================================
// DIAMOND
// ============================================================================

// 5% margin on each side -> 90% usable camera area.
constexpr float IR_DIAMOND_USABLE_FACTOR = 0.90f;

// The bounding box of the physical Diamond has the same aspect ratio
// as the monitor.
//
// res_x and res_y contain the monitor resolution multiplied by the
// fixed-point precision factor (currently << 2).
//
// The common scaling factor cancels out in the ratio.
constexpr float IR_DIAMOND_MONITOR_RATIO = (float)res_x / (float)res_y;

// Usable camera dimensions.
constexpr float IR_DIAMOND_USABLE_WIDTH = (float)CamResX * IR_DIAMOND_USABLE_FACTOR;

constexpr float IR_DIAMOND_USABLE_HEIGHT = (float)CamResY * IR_DIAMOND_USABLE_FACTOR;

// Maximum possible Diamond height allowed by camera width.
constexpr float IR_DIAMOND_HEIGHT_FROM_WIDTH = IR_DIAMOND_USABLE_WIDTH / IR_DIAMOND_MONITOR_RATIO;

// Select the limiting dimension.
constexpr float IR_DIAMOND_HEIGHT = (IR_DIAMOND_HEIGHT_FROM_WIDTH < IR_DIAMOND_USABLE_HEIGHT) ? IR_DIAMOND_HEIGHT_FROM_WIDTH : IR_DIAMOND_USABLE_HEIGHT;

// Preserve monitor aspect ratio.
constexpr float IR_DIAMOND_WIDTH = IR_DIAMOND_HEIGHT * IR_DIAMOND_MONITOR_RATIO;

// Camera center.
constexpr float IR_DIAMOND_CENTER_X = (float)CamResX * 0.5f;

constexpr float IR_DIAMOND_CENTER_Y = (float)CamResY * 0.5f;

// --------------------------------------------------------------------------
// DIAMOND LED DEFAULT COORDINATES
// --------------------------------------------------------------------------

// A = TC (Top Center)
constexpr int IR_LED_DIAMOND_TC_X = (int)IR_DIAMOND_CENTER_X;
constexpr int IR_LED_DIAMOND_TC_Y = (int)(IR_DIAMOND_CENTER_Y - IR_DIAMOND_HEIGHT * 0.5f);
// B = RC (Right Center)
constexpr int IR_LED_DIAMOND_RC_X = (int)(IR_DIAMOND_CENTER_X + IR_DIAMOND_WIDTH * 0.5f);
constexpr int IR_LED_DIAMOND_RC_Y = (int)IR_DIAMOND_CENTER_Y;
// C = BC (Bottom Center)
constexpr int IR_LED_DIAMOND_BC_X = (int)IR_DIAMOND_CENTER_X;
constexpr int IR_LED_DIAMOND_BC_Y = (int)(IR_DIAMOND_CENTER_Y + IR_DIAMOND_HEIGHT * 0.5f);
// D = LC (Left Center)
constexpr int IR_LED_DIAMOND_LC_X = (int)(IR_DIAMOND_CENTER_X - IR_DIAMOND_WIDTH * 0.5f);
constexpr int IR_LED_DIAMOND_LC_Y = (int)IR_DIAMOND_CENTER_Y;


/////////////////////////////

#endif // _OPENFIRECONST_H_
