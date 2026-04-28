/*
 * This module simulates the standard Arduino "Mouse.h" and
 * "Keyboard.h" API for use with the TinyUSB HID API. Instead of doing
 *  #include <HID.h>
 *  #include <Mouse.h>
 *  #include <Keyboard.h>
 *  
 *  Simply do
 *  
 *  #include <TinyUSB_Mouse_Keyboard.h>
 *  
 *  and this module will automatically select whether or not to use the
 *  standard Arduino mouse and keyboard API or the TinyUSB API. We had to
 *  combine them into a single library because of the way TinyUSB handles
 *  descriptors.
 *  
 *  For details on Arduino Mouse.h see
 *   https://www.arduino.cc/reference/en/language/functions/usb/mouse/
 *  For details on Arduino Keyboard.h see
 *   https://www.arduino.cc/reference/en/language/functions/usb/keyboard/
 *
 *  NOTE: This code is derived from the standard Arduino Mouse.h, Mouse.cpp,
 *    Keyboard.h, and Keyboard.cpp code. The copyright on that original code
 *    is as follows.
 *   
 *  Copyright (c) 2015, Arduino LLC
 *  Original code (pre-library): Copyright (c) 2011, Peter Barrett
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#pragma once

#ifndef _TINYUSB_DEVICES_H_
#define _TINYUSB_DEVICES_H_

#include <Arduino.h>
#include <Adafruit_TinyUSB.h> // 696969 inserito per le define di ambiente in particolare collegate a lighgunButtons.cpp
#include "tusb_gamepad16.h"   // 696969 inserito per gamepad16

#ifdef OPENFIRE_WIRELESS_ENABLE
  #ifdef ARDUINO_ARCH_ESP32
    #include "../OpenFIRE_Wireless/ESP32/OpenFIRE_Wireless.h"
  #else // rp2040
    #include "../OpenFIRE_Wireless/RP2040/OpenFIRE_Wireless.h"
  #endif
#endif //OPENFIRE_WIRELESS_ENABLE

/*****************************
 *   GLOBAL SECTION
 *****************************/
#ifndef _GLOBAL_H_
#define _GLOBAL_H_


class TinyUSBDevices_ {
public:
  /// @brief Constructor
  TinyUSBDevices_() {};

  /// @brief Sends and initializes the array of USB devices to the connected host
  void begin(int polRate);

  /// @brief Sends and initializes the array of BT devices for paired devices
  void beginBT(const char *localName, const char *hidName);

  /// @brief Whether the device is running in USB mode (false) or Bluetooth (true)
  bool onBattery = false;

  enum reportChannels_e {
    reportMouse = 0,
    reportKeyboard,
    reportGamepad
  };

  /// @brief Array of which of the three devices have new data that should be reported.
  bool newReport[3] = { false, false, false };

  uint8_t wireless_mode = 0; // 0 = nessuna connessione wireless altro valore connessione // qualsiasi altro valore deve essere diverso da zero
  bool is_pedal_wireless = false; // true se connesso un pedale wifi (può essere true solo se i pedali non sono abbinati a qualche pin del micro)
  volatile uint8_t pedals_wireless_state = 0;
};
extern TinyUSBDevices_ TinyUSBDevices;

//=====================================

#ifdef OPENFIRE_USE_ESPNOW_UNIFIED_PACKET
void report_all_MOUSE_KEY_PAD_TX_wifi_espnow(void);

#endif // OPENFIRE_USE_ESPNOW_UNIFIED_PACKET

// =====================================



enum HID_RID_e{
  HID_RID_KEYBOARD = 1,
  HID_RID_MOUSE,   //2  
  HID_RID_GAMEPAD // 3  
};

uint8_t const desc_hid_report[] = {
  TUD_HID_REPORT_DESC_KEYBOARD(HID_REPORT_ID(HID_RID_e::HID_RID_KEYBOARD)),
  TUD_HID_REPORT_DESC_ABSMOUSE(HID_REPORT_ID(HID_RID_e::HID_RID_MOUSE)),
  TUD_HID_REPORT_DESC_GAMEPAD16(HID_REPORT_ID(HID_RID_e::HID_RID_GAMEPAD))
};

#endif // _GLOBAL_H_


/*****************************
 *   MOUSE SECTION
 *****************************/ 
#ifndef _ABSMOUSE5_H_
#define _ABSMOUSE5_H_

//#include <stdint.h>

// OpenFIRE define           TinyUSB define in hid.c     
#define MOUSE_LEFT    hid_mouse_button_bm_t::MOUSE_BUTTON_LEFT     // 0x01    // 696969 definito in TinyUSB in hid.h
#define MOUSE_RIGHT   hid_mouse_button_bm_t::MOUSE_BUTTON_RIGHT    // 0x02    // 696969 definito in TinyUSB in hid.h
#define MOUSE_MIDDLE  hid_mouse_button_bm_t::MOUSE_BUTTON_MIDDLE   // 0x04    // 696969 definito in TinyUSB in hid.h
#define MOUSE_BUTTON4 hid_mouse_button_bm_t::MOUSE_BUTTON_BACKWARD // 0x08    // 696969 definito in TinyUSB in hid.h
#define MOUSE_BUTTON5 hid_mouse_button_bm_t::MOUSE_BUTTON_FORWARD  // 0x10    // 696969 definito in TinyUSB in hid.h

// 5 button absolute mouse
class AbsMouse5_
{
  //private:
  public:
  hid_abs_mouse_report_t absmouse5Report = {0,0,0,0,0};
  
  public:
  AbsMouse5_();  // si puo' togliere non serve a nulla
	void report(void);
  void move(int16_t x, int16_t y);
  void move_wheel_pan(int8_t wheel, int8_t pan); // NON INDISPENSABILE AGGIUNTA DA ME
  void press(uint8_t b = hid_mouse_button_bm_t::MOUSE_BUTTON_LEFT); // solo un bottone per volta
	void release(uint8_t b = hid_mouse_button_bm_t::MOUSE_BUTTON_LEFT); // solo un bottone per volta
	void releaseAll() { release(0x1f); }
  void click(uint8_t b = hid_mouse_button_bm_t::MOUSE_BUTTON_LEFT); // clicca tutti i pulsanti come impostati e poi li rilascia subito - valutare se tenerlo // NON INDISPENSABILE AGGIUNTA DA ME
  void buttons(uint8_t b); // imposta tutti i bottoni insieme senza farlo una per volta o li disabilita - valutare se tenerlo // NON INDISPENSABILE AGGIUNTA DA ME
  bool isPressed(uint8_t b = hid_mouse_button_bm_t::MOUSE_BUTTON_LEFT);   // NON INDISPENSABILE AGGIUNTA DA ME
};


// global singleton
extern AbsMouse5_ AbsMouse5;

#endif // _ABSMOUSE5_H_

/******************************
 *    KEYBOARD SECTION
 ******************************/
#ifndef _KEYBOARD_H_
#define _KEYBOARD_H_

 //  Keyboard codes
  //  Note these are different in some respects to the TinyUSB codes but 
  //  are compatible with Arduino Keyboard.h API
  
  #define KEY_LEFT_CTRL   0x80
  #define KEY_LEFT_SHIFT  0x81
  #define KEY_LEFT_ALT    0x82
  #define KEY_LEFT_GUI    0x83
  #define KEY_RIGHT_CTRL  0x84
  #define KEY_RIGHT_SHIFT 0x85
  #define KEY_RIGHT_ALT   0x86
  #define KEY_RIGHT_GUI   0x87
  
  #define KEY_UP_ARROW    0xDA
  #define KEY_DOWN_ARROW  0xD9
  #define KEY_LEFT_ARROW  0xD8
  #define KEY_RIGHT_ARROW 0xD7
  #define KEY_BACKSPACE   0xB2
  #define KEY_TAB         0xB3
  #define KEY_RETURN      0xB0
  #define KEY_ESC         0xB1  // 696969 non corrisponde a HID_KEY_ESCAPE di TinyUSB, boh ???
  #define KEY_INSERT      0xD1
  #define KEY_DELETE      0xD4
  #define KEY_PAGE_UP     0xD3
  #define KEY_PAGE_DOWN   0xD6
  #define KEY_HOME        0xD2
  #define KEY_END         0xD5
  #define KEY_CAPS_LOCK   0xC1
  #define KEY_F1          0xC2
  #define KEY_F2          0xC3
  #define KEY_F3          0xC4
  #define KEY_F4          0xC5
  #define KEY_F5          0xC6
  #define KEY_F6          0xC7
  #define KEY_F7          0xC8
  #define KEY_F8          0xC9
  #define KEY_F9          0xCA
  #define KEY_F10         0xCB
  #define KEY_F11         0xCC
  #define KEY_F12         0xCD
  #define KEY_F13         0xF0
  #define KEY_F14         0xF1
  #define KEY_F15         0xF2
  #define KEY_F16         0xF3
  #define KEY_F17         0xF4
  #define KEY_F18         0xF5
  #define KEY_F19         0xF6
  #define KEY_F20         0xF7
  #define KEY_F21         0xF8
  #define KEY_F22         0xF9
  #define KEY_F23         0xFA
  #define KEY_F24         0xFB
  
  #define SHIFT 0x80 
  #define ALT_GR 0xc0 
  #define ISO_KEY 0x64
  #define ISO_REPLACEMENT 0x32

  const uint8_t KeyboardLayout_it_IT[128] =
  {
    0x00,          // NUL
    0x00,          // SOH
    0x00,          // STX
    0x00,          // ETX
    0x00,          // EOT
    0x00,          // ENQ
    0x00,          // ACK
    0x00,          // BEL
    0x2a,          // BS  Backspace
    0x2b,          // TAB Tab
    0x28,          // LF  Enter
    0x00,          // VT
    0x00,          // FF
    0x00,          // CR
    0x00,          // SO
    0x00,          // SI
    0x00,          // DEL
    0x00,          // DC1
    0x00,          // DC2
    0x00,          // DC3
    0x00,          // DC4
    0x00,          // NAK
    0x00,          // SYN
    0x00,          // ETB
    0x00,          // CAN
    0x00,          // EM
    0x00,          // SUB
    0x00,          // ESC
    0x00,          // FS
    0x00,          // GS
    0x00,          // RS
    0x00,          // US
  
    0x2c,          // ' '
    0x1e|SHIFT,    // !
    0x1f|SHIFT,    // "
    0x34|ALT_GR,   // #
    0x21|SHIFT,    // $
    0x22|SHIFT,    // %
    0x23|SHIFT,    // &
    0x2d,          // '
    0x25|SHIFT,    // (
    0x26|SHIFT,    // )
    0x30|SHIFT,    // *
    0x30,          // +
    0x36,          // ,
    0x38,          // -
    0x37,          // .
    0x24|SHIFT,    // /
    0x27,          // 0
    0x1e,          // 1
    0x1f,          // 2
    0x20,          // 3
    0x21,          // 4
    0x22,          // 5
    0x23,          // 6
    0x24,          // 7
    0x25,          // 8
    0x26,          // 9
    0x37|SHIFT,    // :
    0x36|SHIFT,    // ;
    0x32,          // <
    0x27|SHIFT,    // =
    0x32|SHIFT,    // >
    0x2d|SHIFT,    // ?
    0x33|ALT_GR,   // @
    0x04|SHIFT,    // A
    0x05|SHIFT,    // B
    0x06|SHIFT,    // C
    0x07|SHIFT,    // D
    0x08|SHIFT,    // E
    0x09|SHIFT,    // F
    0x0a|SHIFT,    // G
    0x0b|SHIFT,    // H
    0x0c|SHIFT,    // I
    0x0d|SHIFT,    // J
    0x0e|SHIFT,    // K
    0x0f|SHIFT,    // L
    0x10|SHIFT,    // M
    0x11|SHIFT,    // N
    0x12|SHIFT,    // O
    0x13|SHIFT,    // P
    0x14|SHIFT,    // Q
    0x15|SHIFT,    // R
    0x16|SHIFT,    // S
    0x17|SHIFT,    // T
    0x18|SHIFT,    // U
    0x19|SHIFT,    // V
    0x1a|SHIFT,    // W
    0x1b|SHIFT,    // X
    0x1c|SHIFT,    // Y
    0x1d|SHIFT,    // Z
    0x2f|ALT_GR,   // [
    0x35,          // bslash
    0x30|ALT_GR,   // ]
    0x2e|SHIFT,    // ^
    0x38|SHIFT,    // _
    0x00,          // `  not in this layout
    0x04,          // a
    0x05,          // b
    0x06,          // c
    0x07,          // d
    0x08,          // e
    0x09,          // f
    0x0a,          // g
    0x0b,          // h
    0x0c,          // i
    0x0d,          // j
    0x0e,          // k
    0x0f,          // l
    0x10,          // m
    0x11,          // n
    0x12,          // o
    0x13,          // p
    0x14,          // q
    0x15,          // r
    0x16,          // s
    0x17,          // t
    0x18,          // u
    0x19,          // v
    0x1a,          // w
    0x1b,          // x
    0x1c,          // y
    0x1d,          // z
    0x00,          // {  not supported (requires AltGr+Shift)
    0x35|SHIFT,    // |
    0x00,          // }  not supported (requires AltGr+Shift)
    0x00,          // ~  not in this layout
    0x00           // DEL
  };
  
class Keyboard_ : public Print
{
  //private:
  public:
    hid_keyboard_report_t _keyReport = {0,0,{0,0,0,0,0,0}};

  public:
    Keyboard_(void);
    void report(void);
    size_t write(uint8_t k); // mai usata pa serve per print
    size_t write(const uint8_t *buffer, size_t size); // mai usata ma serve per print
    size_t press(uint8_t k);
    size_t release(uint8_t k); // usata ma non ne vedo il senso
    void releaseAll(void); // usata
};

extern Keyboard_ Keyboard;

#endif // _KEYBOARD_H_

/*****************************
 *   GAMEPAD SECTION
 *****************************/
#ifndef _GAMEPAD_H_
#define _GAMEPAD_H_

// Using the Adafruit Gamepad desc as the basis, but with 4 16-bit axis and a 16-bit 15 button report
// OpenFIRE define       TinyUSB define in hid.h    
#define PAD_A      0    // GAMEPAD_BUTTON_A  GAMEPAD_BUTTON_SOUTH
#define PAD_B      1    // GAMEPAD_BUTTON_B  GAMEPAD_BUTTON_EAST
#define PAD_C      2    // GAMEPAD_BUTTON_C
#define PAD_X      3    // GAMEPAD_BUTTON_X  GAMEPAD_BUTTON_NORTH 
#define PAD_Y      4    // GAMEPAD_BUTTON_Y  GAMEPAD_BUTTON_WEST
#define PAD_Z      5    // GAMEPAD_BUTTON_Z
#define PAD_LB     6    // GAMEPAD_BUTTON_TL
#define PAD_RB     7    // GAMEPAD_BUTTON_TR
#define PAD_LT     8    // GAMEPAD_BUTTON_TL2
#define PAD_RT     9    // GAMEPAD_BUTTON_TR2
#define PAD_SELECT 10   // GAMEPAD_BUTTON_SELECT
#define PAD_START  11   // GAMEPAD_BUTTON_START
#define PAD_HOME   12   // GAMEPAD_BUTTON_MODE
#define PAD_LS     13   // GAMEPAD_BUTTON_THUMBL
#define PAD_RS     14   // GAMEPAD_BUTTON_THUMBR   // qui finiscono i 15 bottoni standard
#define PAD_15     15   // GAMEPAD_BUTTON_15      // TUTTI QUELLI DA ORA IN SEGUITO SONO bottoni CUSTOM
#define PAD_16     16   // GAMEPAD_BUTTON_16
#define PAD_17     17   // GAMEPAD_BUTTON_17
#define PAD_18     18   // GAMEPAD_BUTTON_18
#define PAD_19     19   // GAMEPAD_BUTTON_19
#define PAD_20     20   // GAMEPAD_BUTTON_20
#define PAD_21     21   // GAMEPAD_BUTTON_21
#define PAD_22     22   // GAMEPAD_BUTTON_22
#define PAD_23     23   // GAMEPAD_BUTTON_23
#define PAD_24     24   // GAMEPAD_BUTTON_24
#define PAD_25     25   // GAMEPAD_BUTTON_25
#define PAD_26     26   // GAMEPAD_BUTTON_26
#define PAD_27     27   // GAMEPAD_BUTTON_27
#define PAD_28     28   // GAMEPAD_BUTTON_28
#define PAD_29     29   // GAMEPAD_BUTTON_29
#define PAD_30     30   // GAMEPAD_BUTTON_30
#define PAD_31     31   // GAMEPAD_BUTTON_31
// I CODICI DEI TASTI PER IL DPAD (HAT) INZIAZIONO DA 32 OVVERO DOPO I 32 (0-31) CODICI PER I BOTTONI
#define PAD_UP     32   // da qui iniziano i codici per il DPad/hat
#define PAD_DOWN   33   //       "                 "
#define PAD_LEFT   34   //       "                 "
#define PAD_RIGHT  35   //       "                 "

#define START_DPAD_KEY PAD_UP  // utilizzato come riferimento per posizione per file lightgunButtons.cpp

class Gamepad16_ {
  //private:
  public:
    hid_gamepad16_report_t gamepad16Report = {0,0,0,0,0,0,0,0};
    uint16_t _x = 2048; // A META' CONSIDERANDO IL RANGE DA 0 A 4095
    uint16_t _y = 2048; // A META' CONSIDERANDO IL RANGE DA 0 A 4095
  public:
    Gamepad16_(void);
    void moveCam(uint16_t origX, uint16_t origY);
    void moveStick(uint16_t origX, uint16_t origY);
    void press(uint8_t buttonNum);
    void release(uint8_t buttonNum);
    void padUpdate(uint8_t padMask);
    void report(void);
    void releaseAll(void);
    bool stickRight;
  };
  extern Gamepad16_ Gamepad16;

#endif // _GAMEPAD_H_

#endif // _TINYUSB_DEVICES_H_
