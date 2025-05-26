 /*!
 * @file OpenFIREcommon.h
 * @brief Shared methods used throughout the OpenFIRE project.
 *
 * @copyright That One Seong, 2025
 * @copyright GNU Lesser General Public License
 */ 

#include <Arduino.h>
#include <Wire.h>
#ifdef ARDUINO_ARCH_RP2040
  // for RP2040 Wii Clock Gen
  #include <pico/stdlib.h>
  #include <hardware/clocks.h>
  #include <hardware/pwm.h>
#endif // ARDUINO_ARCH_RP2040

#include "OpenFIREcommon.h"
#include "OpenFIREFeedback.h"
#include "OpenFIRElights.h"
#include "OpenFIREserial.h"

// ============ 696969 ========== redifinizione di Serial per gestire le connessione wireless seriali ========
#ifdef OPENFIRE_WIRELESS_ENABLE
    extern Stream* Serial_OpenFIRE_Stream;
    #ifdef Serial
        #define AUX_SERIAL Serial
        #undef Serial
    #endif
    #define Serial (*Serial_OpenFIRE_Stream)
#endif // OPENFIRE_WIRELESS_ENABLE
// ============ 696969 ===== fine redifinizione di Serial per gestire le connessione wireless seriali ========

#ifdef ARDUINO_ARCH_ESP32
    ESP32FIFO esp32_fifo(8);
#endif // ARDUINO_ARCH_ESP32



// button object instance (defined in OpenFIREcommon.h/OpenFIREprefs.h)
LightgunButtons FW_Common::buttons(lgbData, ButtonCount);

void FW_Common::FeedbackSet()
{
    #ifdef USES_RUMBLE
        if(OF_Prefs::pins[OF_Const::rumblePin] >= 0)
            pinMode(OF_Prefs::pins[OF_Const::rumblePin], OUTPUT);
        else OF_Prefs::toggles[OF_Const::rumble] = false;
    #endif // USES_RUMBLE

    #ifdef USES_SOLENOID
        if(OF_Prefs::pins[OF_Const::solenoidPin] >= 0)
            pinMode(OF_Prefs::pins[OF_Const::solenoidPin], OUTPUT);
        else OF_Prefs::toggles[OF_Const::solenoid] = false;
    #endif // USES_SOLENOID

    #ifdef USES_SWITCHES
        #ifdef USES_RUMBLE
            if(OF_Prefs::pins[OF_Const::rumbleSwitch] >= 0)
                pinMode(OF_Prefs::pins[OF_Const::rumbleSwitch], INPUT_PULLUP);
        #endif // USES_RUMBLE

        #ifdef USES_SOLENOID
            if(OF_Prefs::pins[OF_Const::solenoidSwitch] >= 0)
                pinMode(OF_Prefs::pins[OF_Const::solenoidSwitch], INPUT_PULLUP);
        #endif // USES_SOLENOID

        if(OF_Prefs::pins[OF_Const::autofireSwitch] >= 0)
            pinMode(OF_Prefs::pins[OF_Const::autofireSwitch], INPUT_PULLUP);
    #endif // USES_SWITCHES

    #ifdef USES_ANALOG
        analogReadResolution(12);
        #ifdef USES_TEMP
        if(OF_Prefs::pins[OF_Const::analogX] >= 0 && OF_Prefs::pins[OF_Const::analogY] >= 0 &&
           OF_Prefs::pins[OF_Const::analogX] != OF_Prefs::pins[OF_Const::analogY] &&
           OF_Prefs::pins[OF_Const::analogX] != OF_Prefs::pins[OF_Const::tempPin] &&
           OF_Prefs::pins[OF_Const::analogY] != OF_Prefs::pins[OF_Const::tempPin])
        #else
        if(OF_Prefs::pins[OF_Const::analogX] >= 0 && OF_Prefs::pins[OF_Const::analogY] >= 0 &&
           OF_Prefs::pins[OF_Const::analogX] != OF_Prefs::pins[OF_Const::analogY])
        #endif // USES_TEMP
            //pinMode(analogPinX, INPUT);
            //pinMode(analogPinY, INPUT);
            analogIsValid = true;
        else analogIsValid = false;
    #endif // USES_ANALOG

    #if defined(LED_ENABLE) && defined(FOURPIN_LED)
    if(OF_Prefs::pins[OF_Const::ledR] < 0 || OF_Prefs::pins[OF_Const::ledG] < 0 || OF_Prefs::pins[OF_Const::ledB] < 0)
        ledIsValid = false;
    else {
        pinMode(OF_Prefs::pins[OF_Const::ledR], OUTPUT);
        pinMode(OF_Prefs::pins[OF_Const::ledG], OUTPUT);
        pinMode(OF_Prefs::pins[OF_Const::ledB], OUTPUT);
        ledIsValid = true;
    }
    #endif // FOURPIN_LED

    #ifdef CUSTOM_NEOPIXEL
    if(OF_Prefs::pins[OF_Const::neoPixel] >= 0)
        OF_RGB::InitExternPixel(OF_Prefs::pins[OF_Const::neoPixel]);
    #endif // CUSTOM_NEOPIXEL

    #ifdef ARDUINO_ARCH_ESP32
    if (OF_Prefs::pins[OF_Const::periphSCL] >= 0 && OF_Prefs::pins[OF_Const::periphSDA] >= 0) {
    #else //rp2040
    if(OF_Prefs::pins[OF_Const::periphSCL] >= 0 && OF_Prefs::pins[OF_Const::periphSDA] >= 0 &&
       bitRead(OF_Prefs::pins[OF_Const::camSCL], 1) != bitRead(OF_Prefs::pins[OF_Const::periphSCL], 1) &&
       bitRead(OF_Prefs::pins[OF_Const::camSDA], 1) != bitRead(OF_Prefs::pins[OF_Const::periphSDA], 1)) {
    #endif
    #ifdef USES_DISPLAY
        // wrapper will manage display validity
        // check it's not using the camera's I2C line
        if(OF_Prefs::toggles[OF_Const::i2cOLED]) {
            if(!OLED.Begin()) { if(OLED.display != nullptr) delete OLED.display; }
        }
    #endif // USES_DISPLAY
    }
}

void FW_Common::PinsReset()
{
    if(dfrIRPos != nullptr) {
        delete dfrIRPos;
        dfrIRPos = nullptr;
    }

    #ifdef USES_RUMBLE
        if(OF_Prefs::pins[OF_Const::rumblePin] >= 0)
            pinMode(OF_Prefs::pins[OF_Const::rumblePin], INPUT);
    #endif // USES_RUMBLE

    #ifdef USES_SOLENOID
        if(OF_Prefs::pins[OF_Const::solenoidPin] >= 0)
            pinMode(OF_Prefs::pins[OF_Const::solenoidPin], INPUT);
    #endif // USES_SOLENOID

    #ifdef USES_SWITCHES
        #ifdef USES_RUMBLE
            if(OF_Prefs::pins[OF_Const::rumbleSwitch] >= 0)
                pinMode(OF_Prefs::pins[OF_Const::rumbleSwitch], INPUT);
        #endif // USES_RUMBLE

        #ifdef USES_SOLENOID
            if(OF_Prefs::pins[OF_Const::solenoidSwitch] >= 0)
                pinMode(OF_Prefs::pins[OF_Const::solenoidSwitch], INPUT);
        #endif // USES_SOLENOID

        if(OF_Prefs::pins[OF_Const::autofireSwitch] >= 0)
            pinMode(OF_Prefs::pins[OF_Const::autofireSwitch], INPUT);
    #endif // USES_SWITCHES

    #ifdef LED_ENABLE
        OF_RGB::LedOff();

        #ifdef FOURPIN_LED
            if(ledIsValid) {
                pinMode(OF_Prefs::pins[OF_Const::ledR], INPUT);
                pinMode(OF_Prefs::pins[OF_Const::ledG], INPUT);
                pinMode(OF_Prefs::pins[OF_Const::ledB], INPUT);
            }
        #endif // FOURPIN_LED

        #ifdef CUSTOM_NEOPIXEL
            if(OF_RGB::externPixel != nullptr) {
                OF_RGB::externPixel->clear();
                delete OF_RGB::externPixel;
                OF_RGB::externPixel = nullptr;
            }
        #endif // CUSTOM_NEOPIXEL
    #endif // LED_ENABLE

    #ifdef USES_DISPLAY
        if(OLED.display != nullptr)
            OLED.Stop();
    #endif // USES_DISPLAY
}

void FW_Common::CameraSet()
{
    #ifdef ARDUINO_ARCH_ESP32
        Wire.setPins(OF_Prefs::pins[OF_Const::camSDA], OF_Prefs::pins[OF_Const::camSCL]); // MODIFICATO 696969 per ESP32
        dfrIRPos = new DFRobotIRPositionEx(Wire);
    #else // rp2040   
    // Sanity check: which channel do these pins correlate to?
    if(bitRead(OF_Prefs::pins[OF_Const::camSCL], 1) && bitRead(OF_Prefs::pins[OF_Const::camSDA], 1)) {
        // I2C1
        if(bitRead(OF_Prefs::pins[OF_Const::camSCL], 0) && !bitRead(OF_Prefs::pins[OF_Const::camSDA], 0)) {
            // SDA/SCL are indeed on verified correct pins
            Wire1.setSDA(OF_Prefs::pins[OF_Const::camSDA]);
            Wire1.setSCL(OF_Prefs::pins[OF_Const::camSCL]);
            dfrIRPos = new DFRobotIRPositionEx(Wire1);
        }
        
    } else if(!bitRead(OF_Prefs::pins[OF_Const::camSCL], 1) && !bitRead(OF_Prefs::pins[OF_Const::camSDA], 1)) {
        // I2C0
        if(bitRead(OF_Prefs::pins[OF_Const::camSCL], 0) && !bitRead(OF_Prefs::pins[OF_Const::camSDA], 0)) {
            // SDA/SCL are indeed on verified correct pins
            Wire.setSDA(OF_Prefs::pins[OF_Const::camSDA]);
            Wire.setSCL(OF_Prefs::pins[OF_Const::camSCL]);
            dfrIRPos = new DFRobotIRPositionEx(Wire);
        }
    }
    #endif

    #ifdef ARDUINO_ARCH_RP2040
    if(OF_Prefs::pins[OF_Const::wiiClockGen] > -1) {
        set_sys_clock_khz(125000, true);
        gpio_set_function(OF_Prefs::pins[OF_Const::wiiClockGen], GPIO_FUNC_PWM);
        uint slice_num = pwm_gpio_to_slice_num(OF_Prefs::pins[OF_Const::wiiClockGen]);
        pwm_set_clkdiv(slice_num, 1.0f);  // for sys_clock = 125MHz
        //pwm_set_clkdiv(slice_num, 3.0f); // for sys_clock = 150MHz
        pwm_set_wrap(slice_num, 4);
        pwm_set_chan_level(slice_num, pwm_gpio_to_channel(OF_Prefs::pins[OF_Const::wiiClockGen]), 2);
        pwm_set_enabled(slice_num, true);
    } else {
        set_sys_clock_khz(133000, true);
        pwm_set_enabled(pwm_gpio_to_slice_num(OF_Prefs::pins[OF_Const::wiiClockGen]), false);
        gpio_set_function(OF_Prefs::pins[OF_Const::wiiClockGen], GPIO_FUNC_SIO);
        gpio_put(OF_Prefs::pins[OF_Const::wiiClockGen], 0);
        gpio_set_dir(OF_Prefs::pins[OF_Const::wiiClockGen], GPIO_IN);
    }
    #endif // ARDUINO_ARCH_RP2040

    // Start IR Camera with basic data format
    if(dfrIRPos != nullptr) {
        if(!dfrIRPos->begin(DFROBOT_IR_IIC_CLOCK, DFRobotIRPositionEx::DataFormat_Basic, (DFRobotIRPositionEx::Sensitivity_e)OF_Prefs::profiles[OF_Prefs::currentProfile].irSens)) {
            delete dfrIRPos;
            dfrIRPos = nullptr;
            PrintIrError();
        } else camNotAvailable = false;
    } else PrintIrError();
}

void FW_Common::SetMode(const FW_Const::GunMode_e &newMode)
{
    if(gunMode == newMode)
        return;
    
    // exit current mode
    switch(gunMode) {
    case FW_Const::GunMode_Run:
        stateFlags |= FW_Const::StateFlag_PrintPreferences;
        // MAKE SURE EVERYTHING IS DISENGAGED:
        OF_FFB::FFBShutdown();
        buttons.ReleaseAll();
        buttons.ReportDisable();
        break;
    case FW_Const::GunMode_Pause:
        break;
    case FW_Const::GunMode_Docked:
        if(newMode != FW_Const::GunMode_Calibration)
            Serial.println("Undocking.");
        break;
    }
    
    // enter new mode
    gunMode = newMode;
    switch(newMode) {
    case FW_Const::GunMode_Run:
        // begin run mode with all 4 points seen
        lastSeen = 0x0F;

        #ifdef USES_DISPLAY
            if(OLED.serialDisplayType == ExtDisplay::ScreenSerial_Both)
                OLED.ScreenModeChange(ExtDisplay::Screen_Mamehook_Dual);
            else if(OF_Serial::serialMode)
                OLED.ScreenModeChange(ExtDisplay::Screen_Mamehook_Single, buttons.analogOutput);
            else OLED.ScreenModeChange(ExtDisplay::Screen_Normal, buttons.analogOutput);

            OLED.TopPanelUpdate("Prof: ", OF_Prefs::profiles[OF_Prefs::currentProfile].name);
        #endif // USES_DISPLAY

        break;
    case FW_Const::GunMode_Calibration:
        #ifdef USES_DISPLAY
            OLED.ScreenModeChange(ExtDisplay::Screen_Calibrating);
            OLED.TopPanelUpdate("Cali: ", OF_Prefs::profiles[OF_Prefs::currentProfile].name);
        #endif // USES_DISPLAY
        break;
    case FW_Const::GunMode_Pause:
        stateFlags |= FW_Const::StateFlag_SavePreferencesEn | FW_Const::StateFlag_PrintSelectedProfile;
        pauseModeSelection = FW_Const::PauseMode_Calibrate;

        #ifdef USES_DISPLAY
          OLED.ScreenModeChange(ExtDisplay::Screen_Pause);
          OLED.TopPanelUpdate("Using ", OF_Prefs::profiles[OF_Prefs::currentProfile].name);

          if(OF_Prefs::toggles[OF_Const::simplePause]) 
              OLED.PauseListUpdate(pauseModeSelection);
          else OLED.PauseScreenShow(OF_Prefs::currentProfile, OF_Prefs::profiles[0].name, OF_Prefs::profiles[1].name, OF_Prefs::profiles[2].name, OF_Prefs::profiles[3].name);
        #endif // USES_DISPLAY

        break;
    case FW_Const::GunMode_Docked:
        stateFlags |= FW_Const::StateFlag_SavePreferencesEn;

        #ifdef USES_DISPLAY
            OLED.ScreenModeChange(ExtDisplay::Screen_Docked);
        #endif // USES_DISPLAY

        break;
    }

    #ifdef LED_ENABLE
        SetLedColorFromMode();
    #endif // LED_ENABLE
}

void FW_Common::SetRunMode(const FW_Const::RunMode_e &newMode)
{
    if(newMode >= FW_Const::RunMode_Count)
        return;

    // block Processing/test modes being applied to a profile
    if(newMode <= FW_Const::RunMode_ProfileMax && OF_Prefs::profiles[OF_Prefs::currentProfile].runMode != newMode) {
        OF_Prefs::profiles[OF_Prefs::currentProfile].runMode = newMode;
        stateFlags |= FW_Const::StateFlag_SavePreferencesEn;
    }
    
    if(runMode != newMode) {
        runMode = newMode;
        //if(!(stateFlags & FW_Const::StateFlag_PrintSelectedProfile))
            //PrintRunMode();
    }
}

// Dedicated calibration method
void FW_Common::ExecCalMode(const bool &fromDesktop)
{
    buttons.ReportDisable();

    uint8_t calStage = 0;
    char buf[6];
    buf[0] = OF_Const::sCaliInfoUpd;

    // hold values in a buffer till calibration is complete
    int topOffset;
    int bottomOffset;
    int leftOffset;
    int rightOffset;

    // backup current values in case the user cancels
    int _topOffset = OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset;
    int _bottomOffset = OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset;
    int _leftOffset = OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset;
    int _rightOffset = OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset;
    float _TLled = OF_Prefs::profiles[OF_Prefs::currentProfile].TLled;
    float _TRled = OF_Prefs::profiles[OF_Prefs::currentProfile].TRled;
    float _adjX = OF_Prefs::profiles[OF_Prefs::currentProfile].adjX;
    float _adjY = OF_Prefs::profiles[OF_Prefs::currentProfile].adjY;

    // set current values to factory defaults
    OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset = 0;
    OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset = 0;
    OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset = 0;
    OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset = 0;

    // Force center mouse to center
    AbsMouse5.move(32768/2, 32768/2);
    AbsMouse5.report();

    // Initialize current mouse positions (local variables)
    int32_t mouseCurrentX = 32768 / 2;
    int32_t mouseCurrentY = 32768 / 2;

    // Initialize variables for incremental movement
    int32_t mouseTargetX = mouseCurrentX;
    int32_t mouseTargetY = mouseCurrentY;
    bool mouseMoving = false;

    // Jack in, CaliMan, execute!!!
    SetMode(FW_Const::GunMode_Calibration);
    Serial.printf("%c%c", OF_Const::sCaliStageUpd, FW_Const::Cali_Init);

    while(gunMode == FW_Const::GunMode_Calibration) {
        buttons.Poll(1);

        if(irPosUpdateTick) {
            irPosUpdateTick = 0;
            GetPosition();
        }

        // Handle incremental mouse movement
        if (mouseMoving) {
            int32_t deltaX = mouseTargetX - mouseCurrentX;
            int32_t deltaY = mouseTargetY - mouseCurrentY;
            int32_t stepX = 30;
            int32_t stepY = 30;

            if (abs(deltaX) < stepX) stepX = abs(deltaX);
            if (abs(deltaY) < stepY) stepY = abs(deltaY);

            if (deltaX != 0)
                mouseCurrentX += (deltaX > 0) ? stepX : -stepX;

            if (deltaY != 0)
                mouseCurrentY += (deltaY > 0) ? stepY : -stepY;

            AbsMouse5.move(mouseCurrentX, mouseCurrentY);
            AbsMouse5.report();

            if (mouseCurrentX == mouseTargetX && mouseCurrentY == mouseTargetY) {
                mouseMoving = false;
                delay(5);  // Optional small delay
            }
        }

        // Handle button presses and calibration stages
        if((buttons.pressedReleased & (FW_Const::ExitPauseModeBtnMask | FW_Const::ExitPauseModeHoldBtnMask) || Serial.read() == OF_Const::serialTerminator) && !justBooted) {
            Serial.printf("%c%c", OF_Const::sCaliStageUpd, FW_Const::Cali_Verify+1);
            Serial.flush();

            // Reapplying backed up data
            OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset = _topOffset;
            OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset = _bottomOffset;
            OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset = _leftOffset;
            OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset = _rightOffset;
            OF_Prefs::profiles[OF_Prefs::currentProfile].TLled = _TLled;
            OF_Prefs::profiles[OF_Prefs::currentProfile].TRled = _TRled;
            OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = _adjX;
            OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = _adjY;

            // Re-print the profile
            stateFlags |= FW_Const::StateFlag_PrintSelectedProfile;

            // Exit back to docked mode or run mode, depending on if pinged from Desktop App
            if(fromDesktop)
                SetMode(FW_Const::GunMode_Docked);
            else SetMode(FW_Const::GunMode_Run);

            return;
        } else if(buttons.pressed == FW_Const::BtnMask_Trigger && !mouseMoving) {
            Serial.printf("%c%c", OF_Const::sCaliStageUpd, ++calStage);
            // Ensure our messages go through, or else the HID reports eat UART.
            Serial.flush();

            switch(calStage) {
                case FW_Const::Cali_Init:
                    // Initial state, nothing to do (but center cursor for desktop use)
                    if(fromDesktop) {
                        AbsMouse5.move(32768/2, 32768/2);
                        AbsMouse5.report();
                    }
                    break;
                case FW_Const::Cali_Top:
                    // Reset Offsets
                    topOffset = 0;
                    bottomOffset = 0;
                    leftOffset = 0;
                    rightOffset = 0;

                    // Set Cam center offsets
                    if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) {
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = (OpenFIREdiamond.testMedianX() - (512 << 2)) * cos(OpenFIREdiamond.Ang()) -
                                                                     (OpenFIREdiamond.testMedianY() - (384 << 2)) * sin(OpenFIREdiamond.Ang()) + (512 << 2);
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = (OpenFIREdiamond.testMedianX() - (512 << 2)) * sin(OpenFIREdiamond.Ang()) +
                                                                     (OpenFIREdiamond.testMedianY() - (384 << 2)) * cos(OpenFIREdiamond.Ang()) + (384 << 2);
                    } else {
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = (OpenFIREsquare.testMedianX() - (512 << 2)) * cos(OpenFIREsquare.Ang()) -
                                                                     (OpenFIREsquare.testMedianY() - (384 << 2)) * sin(OpenFIREsquare.Ang()) + (512 << 2);
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = (OpenFIREsquare.testMedianX() - (512 << 2)) * sin(OpenFIREsquare.Ang()) +
                                                                     (OpenFIREsquare.testMedianY() - (384 << 2)) * cos(OpenFIREsquare.Ang()) + (384 << 2);
                        // Work out LED locations by assuming height is 100%
                        OF_Prefs::profiles[OF_Prefs::currentProfile].TLled = (res_x / 2) - ((OpenFIREsquare.W() * (res_y  / OpenFIREsquare.H())) / 2);
                        OF_Prefs::profiles[OF_Prefs::currentProfile].TRled = (res_x / 2) + ((OpenFIREsquare.W() * (res_y  / OpenFIREsquare.H())) / 2);
                    }

                    // Update Cam centre in perspective library
                    OpenFIREper.source(OF_Prefs::profiles[OF_Prefs::currentProfile].adjX, OF_Prefs::profiles[OF_Prefs::currentProfile].adjY);
                    OpenFIREper.deinit(0);

                    // Set mouse movement to top position
                    if(!fromDesktop) {
                        mouseTargetX = 32768 / 2;
                        mouseTargetY = 0;
                        mouseMoving = true;
                    }
                    break;
                case FW_Const::Cali_Bottom:
                    // Set Offset buffer
                    topOffset = mouseY;

                    buf[1] = 1;
                    memcpy(&buf[2], &topOffset, sizeof(int));
                    Serial.write(buf, sizeof(buf));
                    Serial.flush();

                    // Set mouse movement to bottom position
                    if(!fromDesktop) {
                        mouseTargetX = 32768 / 2;
                        mouseTargetY = 32767;
                        mouseMoving = true;
                    }
                    break;
                case FW_Const::Cali_Left:
                    // Set Offset buffer
                    bottomOffset = (res_y - mouseY);

                    buf[1] = 2;
                    memcpy(&buf[2], &bottomOffset, sizeof(int));
                    Serial.write(buf, sizeof(buf));
                    Serial.flush();

                    // Set mouse movement to left position
                    if(!fromDesktop) {
                        mouseTargetX = 0;
                        mouseTargetY = 32768 / 2;
                        mouseMoving = true;
                    }
                    break;
                case FW_Const::Cali_Right:
                    // Set Offset buffer
                    leftOffset = mouseX;

                    buf[1] = 3;
                    memcpy(&buf[2], &leftOffset, sizeof(int));
                    Serial.write(buf, sizeof(buf));
                    Serial.flush();

                    // Set mouse movement to right position
                    if(!fromDesktop) {
                        mouseTargetX = 32767;
                        mouseTargetY = 32768 / 2;
                        mouseMoving = true;
                    }
                    break;
                case FW_Const::Cali_Center:
                    // Set Offset buffer
                    rightOffset = (res_x - mouseX);

                    buf[1] = 4;
                    memcpy(&buf[2], &rightOffset, sizeof(int));
                    Serial.write(buf, sizeof(buf));
                    Serial.flush();

                    // Save Offset buffer to profile
                    OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset = topOffset;
                    OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset = bottomOffset;
                    OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset = leftOffset;
                    OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset = rightOffset;

                    // Move back to center calibration point
                    if(!fromDesktop) {
                        mouseTargetX = 32768 / 2;
                        mouseTargetY = 32768 / 2;
                        mouseMoving = true;
                    }
                    break;
                case FW_Const::Cali_Verify:
                    // Apply new Cam center offsets with Offsets applied
                    if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) {
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = (OpenFIREdiamond.testMedianX() - (512 << 2)) * cos(OpenFIREdiamond.Ang()) -
                                                                     (OpenFIREdiamond.testMedianY() - (384 << 2)) * sin(OpenFIREdiamond.Ang()) + (512 << 2);
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = (OpenFIREdiamond.testMedianX() - (512 << 2)) * sin(OpenFIREdiamond.Ang()) +
                                                                     (OpenFIREdiamond.testMedianY() - (384 << 2)) * cos(OpenFIREdiamond.Ang()) + (384 << 2);
                    } else {
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = (OpenFIREsquare.testMedianX() - (512 << 2)) * cos(OpenFIREsquare.Ang()) -
                                                                     (OpenFIREsquare.testMedianY() - (384 << 2)) * sin(OpenFIREsquare.Ang()) + (512 << 2);
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = (OpenFIREsquare.testMedianX() - (512 << 2)) * sin(OpenFIREsquare.Ang()) +
                                                                     (OpenFIREsquare.testMedianY() - (384 << 2)) * cos(OpenFIREsquare.Ang()) + (384 << 2);
                    }

                    buf[1] = 5;
                    memcpy(&buf[2], &OF_Prefs::profiles[OF_Prefs::currentProfile].TLled, sizeof(float));
                    Serial.write(buf, sizeof(buf));
                    Serial.flush();

                    buf[1] = 6;
                    memcpy(&buf[2], &OF_Prefs::profiles[OF_Prefs::currentProfile].TRled, sizeof(float));
                    Serial.write(buf, sizeof(buf));
                    Serial.flush();

                    // Update Cam centre in perspective library
                    OpenFIREper.source(OF_Prefs::profiles[OF_Prefs::currentProfile].adjX, OF_Prefs::profiles[OF_Prefs::currentProfile].adjY);
                    OpenFIREper.deinit(0);

                    // Let the user test.
                    SetMode(FW_Const::GunMode_Verification);
                    while(gunMode == FW_Const::GunMode_Verification) {
                        buttons.Poll();

                        if(irPosUpdateTick) {
                            irPosUpdateTick = 0;
                            GetPosition();
                        }

                        // If it's good, move onto calibration finish.
                        if(buttons.pressed == FW_Const::BtnMask_Trigger) {
                            calStage++;
                            // Stay in Verification Mode; the code outside of the calibration loop will catch us.
                            break;
                        // Press A/B to restart calibration for current profile
                        } else if(buttons.pressedReleased & FW_Const::ExitPauseModeHoldBtnMask) {
                            calStage = 0;
                            Serial.printf("%c%c", OF_Const::sCaliStageUpd, FW_Const::Cali_Init);
                            Serial.flush();

                            // (Re)set current values to factory defaults
                            OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset = 0;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset = 0;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset = 0;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset = 0;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = 512 << 2;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = 384 << 2;
                            SetMode(FW_Const::GunMode_Calibration);
                            AbsMouse5.move(32768/2, 32768/2);
                            AbsMouse5.report();
                        // Press C/Home to exit without committing new calibration values
                        } else if(buttons.pressedReleased & FW_Const::ExitPauseModeBtnMask && !justBooted) {
                            Serial.printf("%c%c", OF_Const::sCaliStageUpd, FW_Const::Cali_Verify+1);
                            Serial.flush();

                            // Reapply backed-up data
                            OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset = _topOffset;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset = _bottomOffset;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset = _leftOffset;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset = _rightOffset;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].TLled = _TLled;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].TRled = _TRled;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = _adjX;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = _adjY;

                            // Re-print the profile
                            stateFlags |= FW_Const::StateFlag_PrintSelectedProfile;

                            // Re-apply the calibration stored in the profile
                            if(fromDesktop)
                                SetMode(FW_Const::GunMode_Docked);
                            else SetMode(FW_Const::GunMode_Run);
                            return;
                        }
                    }
                    break;
                default:
                    break;
            }
        }
    }

    // Break calibration
    if(justBooted) {
        // If this is an initial calibration, save it immediately!
        stateFlags |= FW_Const::StateFlag_SavePreferencesEn;
        SavePreferences();
        if(fromDesktop)
            SetMode(FW_Const::GunMode_Docked);
    } else if(fromDesktop) {
        SetMode(FW_Const::GunMode_Docked);
    } else SetMode(FW_Const::GunMode_Run);

    #ifdef USES_RUMBLE
        if(OF_Prefs::toggles[OF_Const::rumble]) {
            analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength]);
            delay(80);
            #ifdef ARDUINO_ARCH_ESP32
                analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0);  // 696969 per EPS32
            #else // rp2040
            digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
            #endif
            delay(50);
            analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength]);
            delay(125);
            #ifdef ARDUINO_ARCH_ESP32
                analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0);  // 696969 per ESP32
            #else // rp2040            
            digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
            #endif
        }
    #endif // USES_RUMBLE

    Serial.printf("%c%c", OF_Const::sCaliStageUpd, FW_Const::Cali_Verify+1);
    Serial.flush();
}

void FW_Common::GetPosition()
{
    if(dfrIRPos != nullptr) {
        int error = dfrIRPos->basicAtomic(DFRobotIRPositionEx::Retry_2);
        if(error == DFRobotIRPositionEx::Error_Success) {
            // if diamond layout, or square
            if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) {
                OpenFIREdiamond.begin(dfrIRPos->xPositions(), dfrIRPos->yPositions(), dfrIRPos->seen());
                OpenFIREper.warp(OpenFIREdiamond.X(0), OpenFIREdiamond.Y(0),
                                OpenFIREdiamond.X(1), OpenFIREdiamond.Y(1),
                                OpenFIREdiamond.X(2), OpenFIREdiamond.Y(2),
                                OpenFIREdiamond.X(3), OpenFIREdiamond.Y(3),
                                res_x / 2, 0, 0,
                                res_y / 2, res_x / 2,
                                res_y, res_x, res_y / 2);
            } else {
                OpenFIREsquare.begin(dfrIRPos->xPositions(), dfrIRPos->yPositions(), dfrIRPos->seen());
                OpenFIREper.warp(OpenFIREsquare.X(0), OpenFIREsquare.Y(0),
                                OpenFIREsquare.X(1), OpenFIREsquare.Y(1),
                                OpenFIREsquare.X(2), OpenFIREsquare.Y(2),
                                OpenFIREsquare.X(3), OpenFIREsquare.Y(3),
                                OF_Prefs::profiles[OF_Prefs::currentProfile].TLled, 0,
                                OF_Prefs::profiles[OF_Prefs::currentProfile].TRled, 0,
                                OF_Prefs::profiles[OF_Prefs::currentProfile].TLled, res_y,
                                OF_Prefs::profiles[OF_Prefs::currentProfile].TRled, res_y);
            }

            #ifdef CAM_SIMPLE_KALMAN_FILTER
            int aux_getX = OpenFIREper.getX();
            int aux_getY = OpenFIREper.getY();
            Kalman_filter(aux_getX, aux_getY);
            mouseX = map(aux_getX, 0, res_x, (0 - OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset), (res_x + OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset));                 
            mouseY = map(aux_getY, 0, res_y, (0 - OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset), (res_y + OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset));          
            #else
            // Output mapped to screen resolution because offsets are measured in pixels
            mouseX = map(OpenFIREper.getX(), 0, res_x, (0 - OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset), (res_x + OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset));                 
            mouseY = map(OpenFIREper.getY(), 0, res_y, (0 - OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset), (res_y + OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset));
            #endif // CAM_SIMPLE_KALMAN_FILTER


            // =============== 696969 ============== filter ===============================
            /*
            #ifdef CAM_SIMPLE_KALMAN_FILTER
            Kalman_filter(mouseX, mouseY);
            #endif // CAM_SIMPLE_KALMAN_FILTER
            */
            // =============== 696969 ============== filter ===============================

            switch(runMode) {
                case FW_Const::RunMode_Average:
                    // 2 position moving average
                    moveIndex ^= 1;
                    moveXAxisArr[moveIndex] = mouseX;
                    moveYAxisArr[moveIndex] = mouseY;
                    mouseX = (moveXAxisArr[0] + moveXAxisArr[1]) / 2;
                    mouseY = (moveYAxisArr[0] + moveYAxisArr[1]) / 2;
                    break;
                case FW_Const::RunMode_Average2:
                    // weighted average of current position and previous 2
                    if(moveIndex < 2)
                        ++moveIndex;
                    else moveIndex = 0;

                    moveXAxisArr[moveIndex] = mouseX;
                    moveYAxisArr[moveIndex] = mouseY;
                    mouseX = (mouseX + moveXAxisArr[0] + moveXAxisArr[1] + moveXAxisArr[2]) / 4;
                    mouseY = (mouseY + moveYAxisArr[0] + moveYAxisArr[1] + moveYAxisArr[2]) / 4;
                    break;
                default:
                    break;
                }

            // Constrain that bisch so negatives don't cause underflow
            int32_t conMoveX = constrain(mouseX, 0, res_x);
            int32_t conMoveY = constrain(mouseY, 0, res_y);

            if(gunMode == FW_Const::GunMode_Run) {
                UpdateLastSeen();

                if(OF_Serial::serialARcorrection) switch(OF_Prefs::profiles[OF_Prefs::currentProfile].aspectRatio) {
                    case OF_Const::ar16_9:
                        conMoveX = map(conMoveX, 966, 6720, 0, 32767);
                        conMoveX = constrain(conMoveX, 0, 32767);
                        conMoveY = map(conMoveY, 0, res_y, 0, 32767);
                        break;
                    case OF_Const::ar16_10:
                        conMoveX = map(conMoveX, 655, 7048, 0, 32767);
                        conMoveX = constrain(conMoveX, 0, 32767);
                        conMoveY = map(conMoveY, 0, res_y, 0, 32767);
                        break;
                    case OF_Const::ar3_2:
                        conMoveX = map(conMoveX, 438, 7264, 0, 32767);
                        conMoveX = constrain(conMoveX, 0, 32767);
                        conMoveY = map(conMoveY, 0, res_y, 0, 32767);
                        break;
                    case OF_Const::ar5_4:
                        conMoveX = map(conMoveX, 0, res_x, 0, 32767);
                        conMoveY = map(conMoveY, 148, 4182, 0, 32767);
                        conMoveY = constrain(conMoveY, 0, 32767);
                        break;
                    case OF_Const::ar4_3:
                    default:
                        // Output mapped to Mouse resolution
                        conMoveX = map(conMoveX, 0, res_x, 0, 32767);
                        conMoveY = map(conMoveY, 0, res_y, 0, 32767);
                        break;
                } else {
                    // Output mapped to Mouse resolution
                    conMoveX = map(conMoveX, 0, res_x, 0, 32767);
                    conMoveY = map(conMoveY, 0, res_y, 0, 32767);
                }

                bool offXAxis = false;
                bool offYAxis = false;

                if(conMoveX == 0 || conMoveX == 32767)
                    offXAxis = true;
                
                if(conMoveY == 0 || conMoveY == 32767)
                    offYAxis = true;

                if(offXAxis || offYAxis)
                     buttons.offScreen = true;
                else buttons.offScreen = false;

                // ============ 696969 ============================
                /*
                #ifdef CAM_SIMPLE_KALMAN_FILTER
                Kalman_filter(conMoveX, conMoveY);
                #endif // CAM_SIMPLE_KALMAN_FILTER
                */
                // ============ 696969 ============================               
                
                if(buttons.analogOutput)
                     Gamepad16.moveCam(conMoveX, conMoveY);
                else AbsMouse5.move(conMoveX, conMoveY);

            } else if(gunMode == FW_Const::GunMode_Verification) {
                
                // ============ 696969 ============================
                /*
                #ifdef CAM_SIMPLE_KALMAN_FILTER
                Kalman_filter(conMoveX, conMoveY);
                #endif // CAM_SIMPLE_KALMAN_FILTER
                */
                // ============ 696969 ============================
                
                // Output mapped to Mouse resolution
                conMoveX = map(conMoveX, 0, res_x, 0, 32767);
                conMoveY = map(conMoveY, 0, res_y, 0, 32767);

                AbsMouse5.move(conMoveX, conMoveY);
                AbsMouse5.report();
            } else {
                if(millis() - testLastStamp > 50) {
                    testLastStamp = millis();
                    // RAW Camera Output mapped to screen res (1920x1080)
                    int rawX[4];
                    int rawY[4];
                    // RAW Output for viewing in processing sketch mapped to 1920x1080 screen resolution
                    for (int i = 0; i < 4; i++) {
                        if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) {
                            rawX[i] = map(OpenFIREdiamond.X(i), 0, 1023 << 2, 1920, 0);
                            rawY[i] = map(OpenFIREdiamond.Y(i), 0, 768 << 2, 0, 1080);
                        } else {
                            rawX[i] = map(OpenFIREsquare.X(i), 0, 1023 << 2, 0, 1920);
                            rawY[i] = map(OpenFIREsquare.Y(i), 0, 768 << 2, 0, 1080);
                        }
                    }

                    if(runMode == FW_Const::RunMode_Processing) {
                        int mouseXscaled = mouseX / 4;
                        int mouseYscaled = mouseY / 4;

                        if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) {
                            int testMedianX = map(OpenFIREdiamond.testMedianX(), 0, 1023 << 2, 1920, 0);
                            int testMedianY = map(OpenFIREdiamond.testMedianY(), 0, 768 << 2, 0, 1080);
                            char buf[49];
                            buf[0] = OF_Const::sTestCoords;
                            memcpy(&buf[1],  &rawX[0],      sizeof(int));
                            memcpy(&buf[5],  &rawY[0],      sizeof(int));
                            memcpy(&buf[9],  &rawX[1],      sizeof(int));
                            memcpy(&buf[13], &rawY[1],      sizeof(int));
                            memcpy(&buf[17], &rawX[2],      sizeof(int));
                            memcpy(&buf[21], &rawY[2],      sizeof(int));
                            memcpy(&buf[25], &rawX[3],      sizeof(int));
                            memcpy(&buf[29], &rawY[3],      sizeof(int));
                            memcpy(&buf[33], &mouseXscaled, sizeof(int));
                            memcpy(&buf[37], &mouseYscaled, sizeof(int));
                            memcpy(&buf[41], &testMedianX,  sizeof(int));
                            memcpy(&buf[45], &testMedianY,  sizeof(int));
                            Serial.write(buf, sizeof(buf));
                        } else {
                            int testMedianX = map(OpenFIREsquare.testMedianX(), 0, 1023 << 2, 0, 1920);
                            int testMedianY = map(OpenFIREsquare.testMedianY(), 0, 768 << 2, 0, 1080);
                            char buf[49];
                            buf[0] = OF_Const::sTestCoords;
                            memcpy(&buf[1],  &rawX[0],      sizeof(int));
                            memcpy(&buf[5],  &rawY[0],      sizeof(int));
                            memcpy(&buf[9],  &rawX[1],      sizeof(int));
                            memcpy(&buf[13], &rawY[1],      sizeof(int));
                            memcpy(&buf[17], &rawX[2],      sizeof(int));
                            memcpy(&buf[21], &rawY[2],      sizeof(int));
                            memcpy(&buf[25], &rawX[3],      sizeof(int));
                            memcpy(&buf[29], &rawY[3],      sizeof(int));
                            memcpy(&buf[33], &mouseXscaled, sizeof(int));
                            memcpy(&buf[37], &mouseYscaled, sizeof(int));
                            memcpy(&buf[41], &testMedianX,  sizeof(int));
                            memcpy(&buf[45], &testMedianY,  sizeof(int));
                            Serial.write(buf, sizeof(buf));
                        }
                    }

                    #ifdef USES_DISPLAY
                        OLED.DrawVisibleIR(rawX, rawY);
                    #endif // USES_DISPLAY
                }
            }
        } else if(error != DFRobotIRPositionEx::Error_DataMismatch)
            PrintIrError();
    } else  PrintIrError();
}

// ============== 696969 =====================================================================
#ifdef CAM_SIMPLE_KALMAN_FILTER   
void FW_Common::Kalman_filter(int& mouseX,int& mouseY)
{
// --- Blocco di codice per il filtro di Kalman ---
// Inserisci questo blocco direttamente nel tuo loop principale, dopo aver letto
// le variabili raw `mouseX` e `mouseY`.
// Al termine di questo blocco, `mouseX` e `mouseY` conterranno le coordinate filtrate.
//
// Le variabili `mouseX` e `mouseY` devono essere dichiarate (es. `int mouseX; int mouseY;`)
// e contenere le coordinate raw *prima* di questo blocco.
/*
Esempio di utilizzo (nel tuo loop principale):

int mouseX = raw_x_from_sensor; // Leggi la coordinata X grezza dal sensore
int mouseY = raw_y_from_sensor; // Leggi la coordinata Y grezza dal sensore

// --- INIZIO BLOCCO FILTRO DI KALMAN ---
*/
unsigned long currentMillis = millis(); // Tempo attuale in millisecondi.
float dt; // Delta tempo (tempo trascorso) in secondi.

// Calcolo del delta tempo (dt). Gestisce il primo avvio per evitare un dt nullo o troppo grande.
if (lastFilterExecutionTime == 0) {
    dt = 0.01f; // Valore di fallback sicuro per il primo ciclo.
} else {
    dt = (float)(currentMillis - lastFilterExecutionTime) / 1000.0f;
}
lastFilterExecutionTime = currentMillis; // Aggiorna il timestamp per la prossima esecuzione.

// Protezione robusta per dt: evita divisioni per zero o valori irrealistici.
if (dt <= 0.0f || dt > 1.0f) {
    dt = 0.01f;
}

float dt_sq = dt * dt; // Pre-calcolo di dt*dt per efficienza nei calcoli.

// --- Logica di rilevamento del "movimento brusco volontario" ---
// Calcola la "magnitudine al quadrato del salto" tra la misura raw e la posizione stimata attuale.
float deltaX_raw_est = (float)mouseX - estimatedX;
float deltaY_raw_est = (float)mouseY - estimatedY;
float currentResidualMagnitudeSq = deltaX_raw_est * deltaX_raw_est + deltaY_raw_est * deltaY_raw_est;

// Soluzione precedente (con sqrt()):
// float currentResidualMagnitude = sqrt(deltaX_raw_est * deltaX_raw_est + deltaY_raw_est * deltaY_raw_est);

// Inizializza i valori dinamici di R (rumore di misurazione) e Q (incertezza del modello) con i valori base.
float current_R_val = BASE_R;
float current_Q_val = BASE_Q;

// Se il quadrato del "salto" è significativo (supera il quadrato della soglia), aumentiamo Q e R.
// Questa ottimizzazione evita una costosa operazione sqrt() per ogni ciclo.
if (currentResidualMagnitudeSq > BRUSQUE_MOVE_THRESHOLD_SQ) { // Confronto con il quadrato della soglia
// Soluzione precedente (con sqrt()):
// if (currentResidualMagnitude > BRUSQUE_MOVE_THRESHOLD) {
    current_R_val *= BRUSQUE_MOVE_R_MULTIPLIER;
    current_Q_val *= BRUSQUE_MOVE_Q_MULTIPLIER;
}

// --- Adattamento di Q e R in base alla velocità complessiva del mouse ---
// Calcola la magnitudine al quadrato della velocità stimata complessiva per evitare la sqrt().
float totalVelocityMagnitudeSq = velocityX * velocityX + velocityY * velocityY;
// La radice quadrata è ancora necessaria qui se il fattore di adattamento dipende dalla magnitudine lineare.
float totalVelocityMagnitude = sqrt(totalVelocityMagnitudeSq);

// Soluzione precedente (calcolava sqrt() due volte se usata per la soglia):
// float totalVelocityMagnitude = sqrt(velocityX * velocityX + velocityY * velocityY);


// Aumenta Q e R proporzionalmente alla velocità. Ciò permette al filtro di essere più reattivo
// durante i movimenti rapidi e più stabile quando il mouse è quasi fermo.
current_Q_val *= (1.0f + totalVelocityMagnitude * VELOCITY_Q_FACTOR);
current_R_val *= (1.0f + totalVelocityMagnitude * VELOCITY_R_FACTOR);

// Applica limiti massimi per Q e R per prevenire instabilità o valori eccessivi.
current_Q_val = fmin(current_Q_val, Q_MAX_ADAPTIVE);
current_R_val = fmin(current_R_val, R_MAX_ADAPTIVE);

// --- Applicazione del filtro di Kalman per l'asse X ---
// 1. Predizione dello stato per X: stima la prossima posizione e velocità di X.
float predictedPosX = estimatedX + (velocityX * dt) + (0.5f * accelerationX * dt_sq);
float predictedVelX = velocityX + (accelerationX * dt);

// 2. Predizione dell'incertezza (P_x): aumenta l'incertezza per il rumore del modello.
P_x = P_x + current_Q_val;
P_x = constrain(P_x, P_MIN, P_MAX); // Limita P_x per stabilità.

// 3. Calcolo del Guadagno di Kalman (K_x): quanto ci fidiamo della misura corrente.
float K_x = P_x / (P_x + current_R_val);

// 4. Correzione con la misura per X: aggiorna le stime usando la misura raw e K_x.
float residualX = (float)mouseX - predictedPosX;
estimatedX = predictedPosX + K_x * residualX;
velocityX = predictedVelX + K_x * (residualX / dt);

// 5. Aggiornamento dell'accelerazione per X: adattiva per rispondere meglio ai cambi di direzione.
float adaptiveMaxAccelerationX = MAX_ACCELERATION_BASE * (0.5f + 0.5f * (1.0f / (1.0f + constrain(fabs(velocityX), 0.0f, MAX_VELOCITY) * 0.01f)));
accelerationX = constrain(accelerationX + (K_x * 2.0f * residualX) / dt_sq, -adaptiveMaxAccelerationX, adaptiveMaxAccelerationX);

// 6. Limitazione della velocità per X: evita che la velocità stimata cresca indefinitamente.
velocityX = constrain(velocityX, -MAX_VELOCITY, MAX_VELOCITY);

// 7. Aggiornamento della covarianza (P_x): riduce l'incertezza dopo la correzione.
P_x = P_x * (1.0f - K_x);
P_x = constrain(P_x, P_MIN, P_MAX); // Mantiene P_x entro limiti.


// --- Applicazione del filtro di Kalman per l'asse Y ---
// La logica è identica a quella dell'asse X, ma applicata alle variabili Y.
// 1. Predizione dello stato per Y
float predictedPosY = estimatedY + (velocityY * dt) + (0.5f * accelerationY * dt_sq);
float predictedVelY = velocityY + (accelerationY * dt);

// 2. Predizione dell'incertezza (P_y)
P_y = P_y + current_Q_val;
P_y = constrain(P_y, P_MIN, P_MAX);

// 3. Calcolo del Guadagno di Kalman (K_y)
float K_y = P_y / (P_y + current_R_val);

// 4. Correzione con la misura per Y
float residualY = (float)mouseY - predictedPosY;
estimatedY = predictedPosY + K_y * residualY;
velocityY = predictedVelY + K_y * (residualY / dt);

// 5. Aggiornamento dell'accelerazione per Y
float adaptiveMaxAccelerationY = MAX_ACCELERATION_BASE * (0.5f + 0.5f * (1.0f / (1.0f + constrain(fabs(velocityY), 0.0f, MAX_VELOCITY) * 0.01f)));
accelerationY = constrain(accelerationY + (K_y * 2.0f * residualY) / dt_sq, -adaptiveMaxAccelerationY, adaptiveMaxAccelerationY);

// 6. Limitazione della velocità per Y
velocityY = constrain(velocityY, -MAX_VELOCITY, MAX_VELOCITY);

// 7. Aggiornamento della covarianza (P_y)
P_y = P_y * (1.0f - K_y);
P_y = constrain(P_y, P_MIN, P_MAX);

// --- AGGIORNAMENTO DELLE VARIABILI mouseX E mouseY CON I VALORI FILTRATI ---
///////////mouseX = (int)estimatedX;
///////////mouseY = (int)estimatedY;
/*
// --- FINE BLOCCO FILTRO DI KALMAN ---
*/

#ifdef COMMENTO
Serial.printf("RawX:%d RawY:%d FiltX:%.2f FiltY:%.2f VelX:%.1f VelY:%.1f AccX:%.1f AccY:%.1f Px:%.3f Py:%.3f R_eff:%.3f Q_eff:%.3f ResMagSq:%.1f TotalVelMag:%.1f\n",
              mouseX, mouseY, estimatedX, estimatedY, velocityX, velocityY, accelerationX, accelerationY, P_x, P_y, current_R_val, current_Q_val, currentResidualMagnitudeSq, totalVelocityMagnitude);
#endif
}
#endif //CAM_SIMPLE_KALMAN_FILTER   
// ============== 696969 =====================================================================

void FW_Common::PrintIrError()
{
    // set flag to warn desktop app when docking
    if(!camNotAvailable)
        camNotAvailable = true;

    if(dockedSaving) {
        char buf[2] = { OF_Const::sError, OF_Const::sErrCam };
        Serial.write(buf, 2);
    } else if(millis() - camWarningTimestamp > CAM_WARNING_INTERVAL) {
        Serial.println("CAMERROR: Not available");
        camWarningTimestamp = millis();
    }
}

void FW_Common::UpdateLastSeen()
{
    if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) {
        if(lastSeen != OpenFIREdiamond.seen()) {
            #ifdef MAMEHOOKER
            if(!OF_Serial::serialMode)
            #endif // MAMEHOOKER
                #ifdef LED_ENABLE
                if(!lastSeen && OpenFIREdiamond.seen())
                    OF_RGB::LedOff();
                else if(lastSeen && !OpenFIREdiamond.seen())
                    OF_RGB::SetLedPackedColor(OF_RGB::IRSeen0Color);
                #endif // LED_ENABLE

            lastSeen = OpenFIREdiamond.seen();
        }
    } else {
        if(lastSeen != OpenFIREsquare.seen()) {
            #ifdef MAMEHOOKER
            if(!OF_Serial::serialMode) {
            #endif // MAMEHOOKER
                #ifdef LED_ENABLE
                if(!lastSeen && OpenFIREsquare.seen())
                    OF_RGB::LedOff();
                else if(lastSeen && !OpenFIREsquare.seen())
                    OF_RGB::SetLedPackedColor(OF_RGB::IRSeen0Color);
                #endif // LED_ENABLE
            #ifdef MAMEHOOKER
            }
            #endif // MAMEHOOKER
            lastSeen = OpenFIREsquare.seen();
        }
    }
}

bool FW_Common::SelectCalProfile(const int &profile)
{
    if(profile >= PROFILE_COUNT)
        return false;

    if(OF_Prefs::currentProfile != profile) {
        stateFlags |= FW_Const::StateFlag_PrintSelectedProfile;
        OF_Prefs::currentProfile = profile;
    }

    OpenFIREper.source(OF_Prefs::profiles[profile].adjX, OF_Prefs::profiles[profile].adjY);                                                          
    OpenFIREper.deinit(0);

    // set IR sensitivity
    if(OF_Prefs::profiles[profile].irSens <= DFRobotIRPositionEx::Sensitivity_Max)
        SetIrSensitivity((DFRobotIRPositionEx::Sensitivity_e)OF_Prefs::profiles[profile].irSens);

    // set run mode
    if(OF_Prefs::profiles[profile].runMode < FW_Const::RunMode_Count)
        SetRunMode((FW_Const::RunMode_e)OF_Prefs::profiles[profile].runMode);

    #ifdef USES_DISPLAY
        if(gunMode != FW_Const::GunMode_Docked)
            OLED.TopPanelUpdate("Using ", OF_Prefs::profiles[profile].name);
    #endif // USES_DISPLAY
 
    #ifdef LED_ENABLE
        SetLedColorFromMode();
    #endif // LED_ENABLE

    // enable save to allow setting new default profile
    stateFlags |= FW_Const::StateFlag_SavePreferencesEn;
    return true;
}

#ifdef LED_ENABLE
void FW_Common::SetLedColorFromMode()
{
    switch(gunMode) {
    case FW_Const::GunMode_Calibration:
        OF_RGB::SetLedPackedColor(OF_RGB::CalModeColor);
        break;
    case FW_Const::GunMode_Pause:
        OF_RGB::SetLedPackedColor(OF_Prefs::profiles[OF_Prefs::currentProfile].color);
        break;
    case FW_Const::GunMode_Run:
        if(lastSeen)
             OF_RGB::LedOff();
        else OF_RGB::SetLedPackedColor(OF_RGB::IRSeen0Color);
        break;
    default:
        break;
    }
}
#endif // LED_ENABLE

#ifdef USES_DISPLAY
void FW_Common::RedrawDisplay()
{
    if(gunMode == FW_Const::GunMode_Docked)
        OLED.ScreenModeChange(ExtDisplay::Screen_Docked);
    else if(gunMode == FW_Const::GunMode_Pause) {
        OLED.ScreenModeChange(ExtDisplay::Screen_Pause);
        if(OF_Prefs::toggles[OF_Const::simplePause])
            OLED.PauseListUpdate(ExtDisplay::ScreenPause_Save);
        else OLED.PauseScreenShow(OF_Prefs::currentProfile,
                                  OF_Prefs::profiles[0].name,
                                  OF_Prefs::profiles[1].name,
                                  OF_Prefs::profiles[2].name,
                                  OF_Prefs::profiles[3].name);
    }
}
#endif // USES_DISPLAY

void FW_Common::SetIrSensitivity(const int &sensitivity)
{
    if(sensitivity > DFRobotIRPositionEx::Sensitivity_Max)
        return;

    if(OF_Prefs::profiles[OF_Prefs::currentProfile].irSens != sensitivity) {
        OF_Prefs::profiles[OF_Prefs::currentProfile].irSens = sensitivity;
        stateFlags |= FW_Const::StateFlag_SavePreferencesEn;
    }

    if(dfrIRPos != nullptr) dfrIRPos->sensitivityLevel((DFRobotIRPositionEx::Sensitivity_e)sensitivity);
    //if(!(stateFlags & FW_Const::StateFlag_PrintSelectedProfile))
        //PrintIrSensitivity();
}

void FW_Common::SetIrLayout(const int &layout)
{
    // TODO: we need an enum for layout types available
    if(layout >= OF_Const::layoutTypes)
        return;

    if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout != layout) {
        OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout = layout;
        stateFlags |= FW_Const::StateFlag_SavePreferencesEn;
    }
}

int FW_Common::SavePreferences()
{
    // Unless the user's Docked,
    // Only allow one write per pause state until something changes.
    // Extra protection to ensure the same data can't write a bunch of times.
    if(gunMode != FW_Const::GunMode_Docked) {
        if(!(stateFlags & FW_Const::StateFlag_SavePreferencesEn))
            return OF_Prefs::Error_Success;

        stateFlags &= ~FW_Const::StateFlag_SavePreferencesEn;

        #ifdef USES_DISPLAY
            if(OLED.display != nullptr)
                OLED.ScreenModeChange(ExtDisplay::Screen_Saving);
        #endif // USES_DISPLAY
    }

    if(OF_Prefs::SaveProfiles() == OF_Prefs::Error_Success) {
        #ifdef USES_DISPLAY
            OLED.ScreenModeChange(ExtDisplay::Screen_SaveSuccess);
        #endif // USES_DISPLAY

        if(gunMode != FW_Const::GunMode_Docked) Serial.println("Settings saved to Flash"), Serial.flush();
        else Serial.printf("%c%c (Successfully saved to LittleFS Storage)", OF_Const::sSave, true), Serial.flush();
        OF_Prefs::SaveToggles();

        if(OF_Prefs::toggles[OF_Const::customPins])
            OF_Prefs::SavePins();

        OF_Prefs::SaveSettings();
        OF_Prefs::SaveButtons();
        OF_Prefs::SaveUSBID();

        #ifdef LED_ENABLE
            for(uint i = 0; i < 3; ++i) {
                OF_RGB::LedUpdate(25,25,255);
                delay(55);
                OF_RGB::LedOff();
                delay(40);
            }
        #endif // LED_ENABLE

        #ifdef USES_DISPLAY
            RedrawDisplay();
        #endif // USES_DISPLAY

        return OF_Prefs::Error_Success;
    } else {
        #ifdef USES_DISPLAY
            OLED.ScreenModeChange(ExtDisplay::Screen_SaveError);
        #endif // USES_DISPLAY

        // TODO: reimpl a detailed error string
        if(gunMode != FW_Const::GunMode_Docked) Serial.println("Error saving Preferences to Flash.");
        else Serial.printf("%c%c (Failed to save to LittleFS Storage)", OF_Const::sSave, false), Serial.flush();

        /*
        if(nvPrefsError != OF_Prefs::Error_Success) {
            Serial.print(NVRAMlabel);
            Serial.print(" error: ");
        #ifdef SAMCO_FLASH_ENABLE
            Serial.println(OF_Prefs::ErrorCodeToString(nvPrefsError));
        #else
            Serial.println(nvPrefsError);
        #endif // SAMCO_FLASH_ENABLE
        }*/

        #ifdef LED_ENABLE
            for(uint i = 0; i < 2; ++i) {
                OF_RGB::LedUpdate(255,10,5);
                delay(145);
                OF_RGB::LedOff();
                delay(60);
            }
        #endif // LED_ENABLE

        #ifdef USES_DISPLAY
            RedrawDisplay();
        #endif // USES_DISPLAY

        return OF_Prefs::Error_Write;
    }
}

void FW_Common::UpdateBindings(const bool &rebindStrSel)
{
    switch(gunMode) {
    case FW_Const::GunMode_Run:
        buttons.ReleaseAll();
        break;
    case FW_Const::GunMode_Docked:
    case FW_Const::GunMode_Init:
        // Updates pins
        for(int i = 0; i < ButtonCount; ++i)
            LightgunButtons::ButtonDesc[i].pin = OF_Prefs::pins[i];
        break;
    default:
        break;
    }

    if(rebindStrSel) {
        #if defined(PLAYER_START) && defined(PLAYER_SELECT)
        playerStartBtn = PLAYER_START;
        playerSelectBtn = PLAYER_SELECT;
        #else
        if(OF_Prefs::usb.devicePID > 0 && OF_Prefs::usb.devicePID < 5) {
            playerStartBtn = OF_Prefs::usb.devicePID + '0';
            playerSelectBtn = OF_Prefs::usb.devicePID + '4';
        } else {
            playerStartBtn = '1';
            playerSelectBtn = '5';
        }
        #endif // PLAYER_NUMBER
    }

    for(int i = 0; i < ButtonCount; ++i)
        memcpy(&LightgunButtons::ButtonDesc[i].reportType,
               OF_Prefs::backupButtonDesc[i],
               sizeof(OF_Prefs::backupButtonDesc[0]));

    // Updates button functions for low-button mode
    if(OF_Prefs::toggles[OF_Const::lowButtonsMode]) {
        memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_A].reportType2,
               OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_Start],
               2);
        memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_B].reportType2,
               OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_Select],
               2);
    }

    #ifdef MAMEHOOKER
    if(OF_Serial::serialMappingsOffscreenShot) {
        memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Trigger].reportType2,
               OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_A],
               sizeof(LightgunButtons::Desc_s::reportType)*2);

        // remap bindings for low button users to make e.g. VCop 3 playable with 1 btn + pedal
        if(OF_Prefs::toggles[OF_Const::lowButtonsMode])
            memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_A].reportType,
                   OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_Reload],
                   sizeof(LightgunButtons::Desc_s::reportType)*2);
    }

    if(OF_Serial::serialMappingsPedalMode) switch(OF_Serial::serialMappingsPedalMode) {
    case 1:
        memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Pedal].reportType,
               OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_A],
               sizeof(OF_Prefs::backupButtonDesc[0]));
        break;
    case 2:
        memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Pedal].reportType,
               OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_B],
               sizeof(OF_Prefs::backupButtonDesc[0]));
        break;
    default:
        break;
    }
    #endif // MAMEHOOKER

    UpdateStartSelect();
}

void FW_Common::UpdateStartSelect()
{
    uint8_t *btnMatchedPtr;
    for(int i = 0; i < ButtonCount; ++i) {
        do {
            btnMatchedPtr = (uint8_t*)memchr(&LightgunButtons::ButtonDesc[i].reportCode, 0xFF, sizeof(OF_Prefs::backupButtonDesc[i])-1);
            if(btnMatchedPtr != nullptr)
                *btnMatchedPtr = playerStartBtn;
        } while(btnMatchedPtr != nullptr);

        do {
            btnMatchedPtr = (uint8_t*)memchr(&LightgunButtons::ButtonDesc[i].reportCode, 0xFE, sizeof(OF_Prefs::backupButtonDesc[i])-1);
            if(btnMatchedPtr != nullptr)
                *btnMatchedPtr = playerSelectBtn;
        } while(btnMatchedPtr != nullptr);
    }
}



// ============ 696969 ========== ripristino di Serial dopo definizione per connessione seriali ==============
#ifdef OPENFIRE_WIRELESS_ENABLE
    #undef Serial
    #ifdef AUX_SERIAL
        #define Serial AUX_SERIAL
        #undef AuxSerial
    #endif
#endif // OPENFIRE_WIRELESS_ENABLE
// ============ 696969 ===== fine ripristino di Serial dopo definizione per connessione seriali ==============