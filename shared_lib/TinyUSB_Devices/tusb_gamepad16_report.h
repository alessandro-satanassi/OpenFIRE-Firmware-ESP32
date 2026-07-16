/*!
 * @file tusb_gamepad16_report.h
 * @brief tusb_gamepad16_report
 * @n CPP tusb_gamepad16_report
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */


#pragma once

// ===================================================================================
// DATA STRUCTURE: PAYLOAD GAMEPAD 16-BIT
// ===================================================================================
// __attribute__((packed)) is mandatory: it tells the C++ compiler not to insert
// padding bytes (empty space) between variables for memory alignment.
// The absence of padding ensures that the struct matches exactly
// the byte order required by the USB HID protocol for the 16-bit Gamepad.
// /
// __attribute__((packed)) è obbligatorio: dice al compilatore C++ di non inserire 
// byte di padding (spazio vuoto) tra le variabili per allinearle in memoria. 
// L'assenza di padding garantisce che la struct corrisponda millimetricamente 
// al byte-order richiesto dal protocollo HID USB per il Gamepad 16.

// HID Gamepad Protocol Report.
typedef struct __attribute__((packed)) {
    int16_t x;         ///< Delta x  movement of left analog-stick
    int16_t y;         ///< Delta y  movement of left analog-stick
    int16_t z;         ///< Delta z  movement of right analog-joystick
    int16_t rz;        ///< Delta Rz movement of right analog-joystick
    int16_t rx;        ///< Delta Rx movement of analog left trigger
    int16_t ry;        ///< Delta Ry movement of analog right trigger
    uint8_t hat;       ///< Buttons mask for currently pressed buttons in the DPad/hat
    uint32_t buttons;  ///< Buttons mask for currently pressed buttons
} hid_gamepad16_report_t;