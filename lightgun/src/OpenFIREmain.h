/*!
 * @file OpenFIREmain.h
 * @brief OpenFIREmain.h
 * @n CPP OpenFIREmain.h
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
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

// = [ESP32_PORT] ===========================================================================================
#ifdef PAJ7025_CAM
    #include <DFRobotIRPositionEx_Adapter_PAJ7025.h>
#elif CAMERA_DFROBOT_SEN0158
    #include <DFRobotIRPositionEx.h>
#endif // PAJ7025_CAM
// = [ESP32_PORT] ===========================================================================================

#include <OpenFIREBoard.h>
#include "OpenFIREDefines.h"
#include "OpenFIREcommon.h"

#ifdef ARDUINO_ARCH_RP2040
  #include <hardware/pwm.h>
  #include <hardware/irq.h>
  // declare PWM ISR
  void rp2040pwmIrq(void);
#elif defined(ARDUINO_ARCH_ESP32)  // [ESP32_PORT] per ESP32
  hw_timer_t *My_timer = NULL;
  void ARDUINO_ISR_ATTR esp32s3pwmIrq(void);
#endif

#define POLL_RATE 1

// TinyUSB devices interface object that's initialized in MainCoreSetup
// TinyUSBDevices_ TUSBDeviceSetup; // [ESP32_PORT] tolto non serve

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

// ============ VARIABLES AND CONSTANTS ADDED / VARIABILI e COSTANTI AGGIUNTE  ========= [ESP32_PORT]

Stream* Serial_OpenFIRE_Stream; // USED TO HANDLE THE WIRELESS SERIAL / SERVE PER GESTIRE LE SERIALE WIRELESS

// == analog stick handling / per gestione Analog Stick ==
#define ANALOG_STICK_MIN_X 0    // minimum X value / valore minimo X 
#define ANALOG_STICK_MAX_X 4095 // maximum X value / valore massimo X
#define ANALOG_STICK_MIN_Y 0    // minimum Y value / valore minimo Y
#define ANALOG_STICK_MAX_Y 4095 // maximum Y value / valore massimo Y
#define ANALOG_STICK_CENTER_X 2048 // center - used to send data in the simulated joystick, value between 0 and 4095 / centro - serve per inviare dati nel joystic simulato valore tra 0 a 4095
#define ANALOG_STICK_CENTER_Y 2048 // center - used to send data in the simulated joystick, value between 0 and 4095 / centro - serve per inviare dati nel joystic simulato valore tra 0 a 4095
#if defined(ARDUINO_ARCH_RP2040)
    uint16_t ANALOG_STICK_DEADZONE_X_MIN = 1900;  // set to the same original OpenFIRE values / impostati gli stessi valori originali di OpenFIRE 
    uint16_t ANALOG_STICK_DEADZONE_X_MAX = 2200;  // set to the same original OpenFIRE values / impostati gli stessi valori originali di OpenFIRE 
    uint16_t ANALOG_STICK_DEADZONE_Y_MIN = 1900;  // set to the same original OpenFIRE values / impostati gli stessi valori originali di OpenFIRE 
    uint16_t ANALOG_STICK_DEADZONE_Y_MAX = 2200;  // set to the same original OpenFIRE values / impostati gli stessi valori originali di OpenFIRE
    uint16_t ANALOG_STICK_DEADZONE_X_CENTER = ANALOG_STICK_CENTER_X;  // used to send data in the simulated joystick, value between 0 and 4095 / serve per inviare dati nel joystic simulato valore tra 0 a 4095
    uint16_t ANALOG_STICK_DEADZONE_Y_CENTER = ANALOG_STICK_CENTER_Y;  // used to send data in the simulated joystick, value between 0 and 4095 / serve per inviare dati nel joystic simulato valore tra 0 a 4095
#elif defined(ARDUINO_ARCH_ESP32)
    uint16_t ANALOG_STICK_DEADZONE_X_MIN = ANALOG_STICK_MAX_X;  // used to read the values of the joystick connected to the MCU (deadzone relative to the center, calculated during setup) / serve per leggere i valori del joystick collegato al micro (deadzone rispetto al centro, lo calcola in fate di setup)  
    uint16_t ANALOG_STICK_DEADZONE_X_MAX = ANALOG_STICK_MIN_X;  // used to read the values of the joystick connected to the MCU (deadzone relative to the center, calculated during setup) / serve per leggere i valori del joystick collegato al micro (deadzone rispetto al centro, lo calcola in fate di setup)
    uint16_t ANALOG_STICK_DEADZONE_Y_MIN = ANALOG_STICK_MAX_Y;  // used to read the values of the joystick connected to the MCU (deadzone relative to the center, calculated during setup) / serve per leggere i valori del joystick collegato al micro (deadzone rispetto al centro, lo calcola in fate di setup)
    uint16_t ANALOG_STICK_DEADZONE_Y_MAX = ANALOG_STICK_MIN_Y;  // used to read the values of the joystick connected to the MCU (deadzone relative to the center, calculated during setup) / serve per leggere i valori del joystick collegato al micro (deadzone rispetto al centro, lo calcola in fate di setup)
    uint16_t ANALOG_STICK_DEADZONE_X_CENTER = ANALOG_STICK_CENTER_X;  // used to send data in the simulated joystick, value between 0 and 4095 / serve per inviare dati nel joystic simulato valore tra 0 a 4095
    uint16_t ANALOG_STICK_DEADZONE_Y_CENTER = ANALOG_STICK_CENTER_Y;  // used to send data in the simulated joystick, value between 0 and 4095 / serve per inviare dati nel joystic simulato valore tra 0 a 4095
#endif
// == END analog stick handling / FINE per gestione Analog Stick ==


// ============ FUNCTION DEFINITIONS / DEFINIZIONE DELLE FUNZIONI ================== [ESP32_PORT]
// for PlatformIO - in Arduino IDE they don’t need to be defined beforehand / per platformio - in arduino IDE si possono anche non definire prima
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