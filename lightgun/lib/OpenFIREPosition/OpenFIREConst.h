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
// OPENFIRE GLOBAL COORDINATE SPACE
// ============================================================================

// Camera-independent perspective coordinate space.
constexpr int res_x = 1920 << 2;
constexpr int res_y = 1080 << 2;


// ============================================================================
// GLOBAL IR GEOMETRY
// ============================================================================

// Physical Square LED rectangle base / height ratio.
// Used by Square tracking, profile defaults and calibration.
constexpr float IR_SQUARE_BASE_HEIGHT_RATIO = 0.711f;

#endif // _OPENFIRECONST_H_