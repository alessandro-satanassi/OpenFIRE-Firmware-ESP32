/*!
 * @file OpenFIREmain.h
 * @brief OpenFIRE main control program.
 *
 * @copyright Samco, https://github.com/samuelballantyne, June 2020
 * @copyright Mike Lynch, July 2021
 * @copyright That One Seong, https://github.com/SeongGino, 2024
 * @copyright GNU Lesser General Public License
 *
 * @author [Sam Ballantyne](samuelballantyne@hotmail.com)
 * @author Mike Lynch
 * @author [That One Seong](SeongsSeongs@gmail.com)
 * @date 2025
 */

#ifndef _OPENFIREMAIN_H_
#define _OPENFIREMAIN_H_

#include <Arduino.h>
#include <Wire.h>
// include TinyUSB or HID depending on USB stack option
#if defined(USE_TINYUSB)
    #include <Adafruit_TinyUSB.h>
#elif defined(CFG_TUSB_MCU)
    #error Incompatible USB stack. Use Adafruit TinyUSB.
#endif

#include <DFRobotIRPositionEx.h>

#include <OpenFIREBoard.h>
#include "OpenFIREDefines.h"
#include "OpenFIREcommon.h"

#ifdef ARDUINO_ARCH_RP2040
  #include <hardware/pwm.h>
  #include <hardware/irq.h>
  // declare PWM ISR
  void rp2040pwmIrq(void);
#elif defined(ARDUINO_ARCH_ESP32)  // 696969 per ESP32
  hw_timer_t *My_timer = NULL;
  void ARDUINO_ISR_ATTR esp32s3pwmIrq(void);
#endif

#define POLL_RATE 1

// TinyUSB devices interface object that's initialized in MainCoreSetup
// TinyUSBDevices_ TUSBDeviceSetup; // 696969 tolto non serve

// Selector for which profile in the profile selector of the simple pause menu you're picking.
uint8_t profileModeSelection = 0;
// Flag to tell if we're in the profile selector submenu of the simple pause menu.
bool pauseModeSelectingProfile = false;

// Timestamp of when we started holding a buttons combo.
unsigned long pauseHoldStartstamp;
bool pauseHoldStarted = false;
bool pauseExitHoldStarted = false;

// Timestamp of last USB packet update.
unsigned long lastUSBpoll = 0;

uint32_t fifoData = 0;

// ============ VARIABILI e COSTANTI AGGIUNTE DA ME ========= 696969

Stream* Serial_OpenFIRE_Stream; // SERVE PER GESTIRE LE SERIALE WIRELESS

// == per gestione Analog Stick ==
#define ANALOG_STICK_MIN_X 0 // al momento non serve
#define ANALOG_STICK_MAX_X 4095 // al momento non serve
#define ANALOG_STICK_MIN_Y 0 // al momento non serve
#define ANALOG_STICK_MAX_Y 4095 // al momento non serve
#define ANALOG_STICK_CENTER_X 2048 // serve per inviare dati nel joystic simulato valore tra 0 a 4095
#define ANALOG_STICK_CENTER_Y 2048 // serve per inviare dati nel joystic simulato valore tra 0 a 4095
#if defined(ARDUINO_ARCH_RP2040)
    uint16_t ANALOG_STICK_DEADZONE_X_MIN = 1900;  // impostati gli stessi valori originali di OpenFIRE 
    uint16_t ANALOG_STICK_DEADZONE_X_MAX = 2200;  // impostati gli stessi valori originali di OpenFIRE 
    uint16_t ANALOG_STICK_DEADZONE_Y_MIN = 1900;  // impostati gli stessi valori originali di OpenFIRE 
    uint16_t ANALOG_STICK_DEADZONE_Y_MAX = 2200;  // impostati gli stessi valori originali di OpenFIRE
    uint16_t ANALOG_STICK_DEADZONE_X_CENTER = ANALOG_STICK_CENTER_X;  // serve per inviare dati nel joystic simulato valore tra 0 a 4095
    uint16_t ANALOG_STICK_DEADZONE_Y_CENTER = ANALOG_STICK_CENTER_Y;  // serve per inviare dati nel joystic simulato valore tra 0 a 4095
#elif defined(ARDUINO_ARCH_ESP32)
    uint16_t ANALOG_STICK_DEADZONE_X_MIN = ANALOG_STICK_MAX_X;  // serve per leggere i valori del joystick collegato al micro (deadzone rispetto al centro, lo calcola in fate di setup)  
    uint16_t ANALOG_STICK_DEADZONE_X_MAX = ANALOG_STICK_MIN_X;  // serve per leggere i valori del joystick collegato al micro (deadzone rispetto al centro, lo calcola in fate di setup)
    uint16_t ANALOG_STICK_DEADZONE_Y_MIN = ANALOG_STICK_MAX_Y;  // serve per leggere i valori del joystick collegato al micro (deadzone rispetto al centro, lo calcola in fate di setup)
    uint16_t ANALOG_STICK_DEADZONE_Y_MAX = ANALOG_STICK_MIN_Y;  // serve per leggere i valori del joystick collegato al micro (deadzone rispetto al centro, lo calcola in fate di setup)
    uint16_t ANALOG_STICK_DEADZONE_X_CENTER = ANALOG_STICK_CENTER_X;  // serve per inviare dati nel joystic simulato valore tra 0 a 4095
    uint16_t ANALOG_STICK_DEADZONE_Y_CENTER = ANALOG_STICK_CENTER_Y;  // serve per inviare dati nel joystic simulato valore tra 0 a 4095
#endif
// == FINE per gestione Analog Stick ==


// ============ DEFINIZIONE DELLE FUNZIONI ================== 696969 per platformio - in arduino IDE si possono anche nonn definire prima
void startIrCamTimer(const int &frequencyHz);
void ExecGunModeDocked();
void TriggerFire();
void TriggerNotFire();
void TriggerFireSimple();
void TriggerNotFireSimple();
void AnalogStickPoll();
void SendEscapeKey();
void SetProfileSelection(const bool &isIncrement);
void SetPauseModeSelection(const bool &isIncrement);
void RumbleToggle();
void SolenoidToggle();
void IncreaseIrSensitivity(const uint32_t &sens);
void DecreaseIrSensitivity(const uint32_t &sens);
void OffscreenToggle();
void AutofireSpeedToggle();
void SelectCalProfileFromBtnMask(const uint32_t &mask);
void ExecRunMode();

#ifdef ARDUINO_ARCH_RP2040
  void rp2040EnablePWMTimer(const unsigned int &slice_num, const unsigned int &frequency);
#endif

// ==========================================================

#endif // _SAMCOENHANCED_H_