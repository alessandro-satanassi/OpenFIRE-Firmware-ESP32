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

#ifdef PAJ7025_CAM
    // PixArt PAH7025R2 IR positioning camera resolution
    constexpr int CamResX = 1024; //4096;
    constexpr int CamResY = 768;  //4096;
#else // DF ROBOT
    // DFRobot IR positioning camera resolution
    constexpr int CamResX = 1024;
    constexpr int CamResY = 768;
#endif // PAJ7025_CAM

// DFRobot IR positioning camera maximum X and Y
constexpr int CamMaxX = CamResX - 1;
constexpr int CamMaxY = CamResY - 1;

#ifdef PAJ7025_CAM
    // shift amount for extra precision for the maths
    // since the median is an average of 4 values, use 2 more bits
    constexpr int CamToMouseShift = 2; // 0;
#else  // DF ROBOT
    // shift amount for extra precision for the maths
    // since the median is an average of 4 values, use 2 more bits
    constexpr int CamToMouseShift = 2;
#endif // PAJ7025_CAM

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

#endif // _OPENFIRECONST_H_
