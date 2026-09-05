/*!
 * @file OpenFIREcommon.h
 * @brief Shared methods used throughout the OpenFIRE project.
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
 * @copyright That One Seong, 2025
 * @copyright GNU Lesser General Public License
 */ 

#include <Arduino.h>
#include <Wire.h>



#include "OpenFIREcommon.h"
#include "OpenFIREFeedback.h"
#include "OpenFIRElights.h"
#include "OpenFIREserial.h"

#ifdef ARDUINO_ARCH_ESP32
    #include "esp32-hal-tinyusb.h"
#endif

// ============ [ESP32_PORT] ============
// Definition of Serial for managing wireless serial connections / redifinizione di Serial per gestire le connessione wireless seriali
#ifdef OPENFIRE_WIRELESS_ENABLE
    extern Stream* Serial_OpenFIRE_Stream;
    #ifdef Serial
        #define AUX_SERIAL Serial
        #undef Serial
    #endif
    #define Serial (*Serial_OpenFIRE_Stream)
#endif // OPENFIRE_WIRELESS_ENABLE
// ============ [ESP32_PORT] ============
// End definition of Serial for managing wireless serial connections / fine redifinizione di Serial per gestire le connessione wireless seriali

#ifdef ARDUINO_ARCH_ESP32  // [ESP32_PORT]
    #define delay(ms) vTaskDelay(pdMS_TO_TICKS(ms))                    
#endif //ARDUINO_ARCH_ESP32

#ifdef ARDUINO_ARCH_ESP32
    ESP32FIFO esp32_fifo(8);
#endif // ARDUINO_ARCH_ESP32



// button object instance (defined in OpenFIREcommon.h/OpenFIREprefs.h)
LightgunButtons FW_Common::buttons(lgbData, ButtonCount);

void FW_Common::RebootToBootloader()
{
    #ifdef USES_DISPLAY
    if(OLED.display != nullptr) {
        OLED.display->clearDisplay();
        OLED.display->setTextColor(WHITE, BLACK);
        OLED.display->setTextSize(1);

        FW_Common::OLED.display->setCursor(37, 16);
        FW_Common::OLED.display->print("Ready for");

        FW_Common::OLED.display->setCursor(40, 30);
        FW_Common::OLED.display->print("firmware");

        FW_Common::OLED.display->setCursor(46, 44);
        FW_Common::OLED.display->print("update");

        // Transfer the framebuffer to the OLED before rebooting.
        OLED.display->display();
    }
    #endif

    #ifdef ARDUINO_ARCH_ESP32
        usb_persist_restart(RESTART_BOOTLOADER);
    #elif defined(ARDUINO_ARCH_RP2040)
        rp2040.rebootToBootloader();
    #endif
}

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
            if(!OLED.Begin()) { 
                if(OLED.display != nullptr) {
                    delete OLED.display; 
                    OLED.display = nullptr; 
                }
            }
        }
    #endif // USES_DISPLAY
    }
}

void FW_Common::PinsReset(const int8_t *pinMap)
{
    if(pinMap == nullptr)
        pinMap = OF_Prefs::pins;

    OpenFIRECamera::End();

    #ifdef USES_RUMBLE
        if(pinMap[OF_Const::rumblePin] >= 0)
            pinMode(pinMap[OF_Const::rumblePin], INPUT);
    #endif

    #ifdef USES_SOLENOID
        if(pinMap[OF_Const::solenoidPin] >= 0)
            pinMode(pinMap[OF_Const::solenoidPin], INPUT);
    #endif

    #ifdef USES_SWITCHES
        #ifdef USES_RUMBLE
            if(pinMap[OF_Const::rumbleSwitch] >= 0)
                pinMode(pinMap[OF_Const::rumbleSwitch], INPUT);
        #endif

        #ifdef USES_SOLENOID
            if(pinMap[OF_Const::solenoidSwitch] >= 0)
                pinMode(pinMap[OF_Const::solenoidSwitch], INPUT);
        #endif

        if(pinMap[OF_Const::autofireSwitch] >= 0)
            pinMode(pinMap[OF_Const::autofireSwitch], INPUT);
    #endif

    #ifdef LED_ENABLE
        // LedOff() uses OF_Prefs::pins. Call it only when that array still
        // describes the hardware currently active.
        if(pinMap == OF_Prefs::pins)
            OF_RGB::LedOff();

        #ifdef FOURPIN_LED
            if(ledIsValid) {
                if(pinMap[OF_Const::ledR] >= 0)
                    pinMode(pinMap[OF_Const::ledR], INPUT);
                if(pinMap[OF_Const::ledG] >= 0)
                    pinMode(pinMap[OF_Const::ledG], INPUT);
                if(pinMap[OF_Const::ledB] >= 0)
                    pinMode(pinMap[OF_Const::ledB], INPUT);
            }
            ledIsValid = false;
        #endif

        #ifdef CUSTOM_NEOPIXEL
            if(OF_RGB::externPixel != nullptr) {
                OF_RGB::externPixel->clear();
                OF_RGB::externPixel->show();
                delete OF_RGB::externPixel;
                OF_RGB::externPixel = nullptr;
            }
        #endif
    #endif

    #ifdef USES_DISPLAY
        if(OLED.display != nullptr)
            OLED.Stop();
    #endif
}

void FW_Common::CameraSet()
{
    if (!OpenFIRECamera::Begin()) {
        PrintIrError();
        return;
    }

    const CameraProfile& profile = OpenFIRECamera::Profile();

    OF_Prefs::InitProfileDefaults(profile);

    OpenFIREsquare.configure(profile);
    OpenFIREdiamond.configure(profile);
    OpenFIREper.configure(profile);

    #ifdef USE_MULTI_ONE_EURO_FILTER
        oef_multi.configure(profile);
    #endif

    OpenFIREper.source(
        OF_Prefs::profiles[OF_Prefs::currentProfile].adjX,
        OF_Prefs::profiles[OF_Prefs::currentProfile].adjY);
    OpenFIREper.deinit(0);
    
    camNotAvailable = false;
}



void FW_Common::SetMode(const FW_Const::GunMode_e &newMode)
{
    #ifdef USE_SQUARE_ADVANCED
        OpenFIREsquare.setCalibrationMode(newMode == FW_Const::GunMode_Calibration);

        if(newMode != FW_Const::GunMode_Calibration && OF_Prefs::currentProfile < PROFILE_COUNT && !OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout)
            OpenFIREsquare.setWideLayout((OF_Prefs::profiles[OF_Prefs::currentProfile].TRled - OF_Prefs::profiles[OF_Prefs::currentProfile].TLled) > res_y);
    #endif

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
        // Docked mode has no hardware resources that need explicit release.
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
    }
}

void FW_Common::ExecCalMode(const bool &fromDesktop)
{
    buttons.ReportDisable();

    uint8_t calStage = 0;
    bool communicationFailed = false;
    uint8_t caliPayload[5];

    // Queste sostituiscono i vecchi valori hardcoded (512 e 384) e 
    // garantiscono una calibrazione perfetta sia per la DFRobot (4:3) che per la PixArt (1:1).
    // Center of the unified Mouse coordinate space.
    const CameraProfile& cameraProfile = OpenFIRECamera::Profile();
    const int CENTER_X = cameraProfile.mouseResX / 2;
    const int CENTER_Y = cameraProfile.mouseResY / 2;

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

    // Force center mouse to center (Absolute HID range 0-32767)
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
    if(fromDesktop) {
        const uint8_t stage = FW_Const::Cali_Init;
        if(!OF_Serial::AppSerialSendResponse(OF_Const::sCaliStageUpd, &stage, 1))
            goto calibration_failed;
    }

    while(gunMode == FW_Const::GunMode_Calibration) {
        buttons.Poll(1);

        if(fromDesktop)
            OF_Serial::SerialProcessingDocked();
        const bool desktopCancel = fromDesktop && OF_Serial::AppSerialTakeCalibrationCancel();

        if(irPosUpdateTick) {
            irPosUpdateTick = 0;
            GetPosition();
        }

        if(fromDesktop && camNotAvailable)
            goto calibration_failed;

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
        if(((buttons.pressedReleased & (FW_Const::ExitPauseModeBtnMask | FW_Const::ExitPauseModeHoldBtnMask)) && !justBooted) ||
           desktopCancel) {
            goto calibration_cancelled;
        } else if(buttons.pressed == FW_Const::BtnMask_Trigger && !mouseMoving) {
            ++calStage;
            if(fromDesktop && !OF_Serial::AppSerialSendResponse(OF_Const::sCaliStageUpd, &calStage, 1))
                goto calibration_failed;

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

                    // Set Cam center offsets using dynamic CENTER_X and CENTER_Y
                    if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) {
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = (OpenFIREdiamond.testMedianX() - CENTER_X) * cos(OpenFIREdiamond.Ang()) -
                                                                             (OpenFIREdiamond.testMedianY() - CENTER_Y) * sin(OpenFIREdiamond.Ang()) + CENTER_X;
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = (OpenFIREdiamond.testMedianX() - CENTER_X) * sin(OpenFIREdiamond.Ang()) +
                                                                             (OpenFIREdiamond.testMedianY() - CENTER_Y) * cos(OpenFIREdiamond.Ang()) + CENTER_Y;
                    } else {
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = (OpenFIREsquare.testMedianX() - CENTER_X) * cos(OpenFIREsquare.Ang()) -
                                                                             (OpenFIREsquare.testMedianY() - CENTER_Y) * sin(OpenFIREsquare.Ang()) + CENTER_X;
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = (OpenFIREsquare.testMedianX() - CENTER_X) * sin(OpenFIREsquare.Ang()) +
                                                                             (OpenFIREsquare.testMedianY() - CENTER_Y) * cos(OpenFIREsquare.Ang()) + CENTER_Y;
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

                    caliPayload[0] = 1;
                    memcpy(&caliPayload[1], &topOffset, sizeof(int));
                    if(fromDesktop && !OF_Serial::AppSerialSendResponse(OF_Const::sCaliInfoUpd, caliPayload, sizeof(caliPayload)))
                        goto calibration_failed;

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

                    caliPayload[0] = 2;
                    memcpy(&caliPayload[1], &bottomOffset, sizeof(int));
                    if(fromDesktop && !OF_Serial::AppSerialSendResponse(OF_Const::sCaliInfoUpd, caliPayload, sizeof(caliPayload)))
                        goto calibration_failed;

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

                    caliPayload[0] = 3;
                    memcpy(&caliPayload[1], &leftOffset, sizeof(int));
                    if(fromDesktop && !OF_Serial::AppSerialSendResponse(OF_Const::sCaliInfoUpd, caliPayload, sizeof(caliPayload)))
                        goto calibration_failed;

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

                    caliPayload[0] = 4;
                    memcpy(&caliPayload[1], &rightOffset, sizeof(int));
                    if(fromDesktop && !OF_Serial::AppSerialSendResponse(OF_Const::sCaliInfoUpd, caliPayload, sizeof(caliPayload)))
                        goto calibration_failed;

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
                    // Apply new Cam center offsets with Offsets applied using dynamic CENTER_X and CENTER_Y
                    if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) {
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = (OpenFIREdiamond.testMedianX() - CENTER_X) * cos(OpenFIREdiamond.Ang()) -
                                                                             (OpenFIREdiamond.testMedianY() - CENTER_Y) * sin(OpenFIREdiamond.Ang()) + CENTER_X;
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = (OpenFIREdiamond.testMedianX() - CENTER_X) * sin(OpenFIREdiamond.Ang()) +
                                                                             (OpenFIREdiamond.testMedianY() - CENTER_Y) * cos(OpenFIREdiamond.Ang()) + CENTER_Y;
                    } else {
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = (OpenFIREsquare.testMedianX() - CENTER_X) * cos(OpenFIREsquare.Ang()) -
                                                                             (OpenFIREsquare.testMedianY() - CENTER_Y) * sin(OpenFIREsquare.Ang()) + CENTER_X;
                        OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = (OpenFIREsquare.testMedianX() - CENTER_X) * sin(OpenFIREsquare.Ang()) +
                                                                             (OpenFIREsquare.testMedianY() - CENTER_Y) * cos(OpenFIREsquare.Ang()) + CENTER_Y;
                    }

                    caliPayload[0] = 5;
                    memcpy(&caliPayload[1], &OF_Prefs::profiles[OF_Prefs::currentProfile].TLled, sizeof(float));
                    if(fromDesktop && !OF_Serial::AppSerialSendResponse(OF_Const::sCaliInfoUpd, caliPayload, sizeof(caliPayload)))
                        goto calibration_failed;

                    caliPayload[0] = 6;
                    memcpy(&caliPayload[1], &OF_Prefs::profiles[OF_Prefs::currentProfile].TRled, sizeof(float));
                    if(fromDesktop && !OF_Serial::AppSerialSendResponse(OF_Const::sCaliInfoUpd, caliPayload, sizeof(caliPayload)))
                        goto calibration_failed;

                    // Update Cam centre in perspective library
                    OpenFIREper.source(OF_Prefs::profiles[OF_Prefs::currentProfile].adjX, OF_Prefs::profiles[OF_Prefs::currentProfile].adjY);
                    OpenFIREper.deinit(0);

                    // Let the user test.
                    SetMode(FW_Const::GunMode_Verification);
                    while(gunMode == FW_Const::GunMode_Verification) {
                        buttons.Poll();

                        if(fromDesktop)
                            OF_Serial::SerialProcessingDocked();
                        const bool verificationCancel = fromDesktop && OF_Serial::AppSerialTakeCalibrationCancel();

                        if(irPosUpdateTick) {
                            irPosUpdateTick = 0;
                            GetPosition();
                        }

                        if(fromDesktop && camNotAvailable)
                            goto calibration_failed;

                        // Cancellation wins over a simultaneous trigger/restart.
                        if(verificationCancel ||
                           ((buttons.pressedReleased & FW_Const::ExitPauseModeBtnMask) && !justBooted))
                            goto calibration_cancelled;

                        // If it's good, move onto calibration finish.
                        if(buttons.pressed == FW_Const::BtnMask_Trigger) {
                            calStage++;
                            // Stay in Verification Mode; the code outside of the calibration loop will catch us.
                            break;
                        // Press A/B to restart calibration for current profile
                        } else if(buttons.pressedReleased & FW_Const::ExitPauseModeHoldBtnMask) {
                            calStage = 0;
                            const uint8_t stage = FW_Const::Cali_Init;
                            if(fromDesktop && !OF_Serial::AppSerialSendResponse(OF_Const::sCaliStageUpd, &stage, 1))
                                goto calibration_failed;

                            // (Re)set current values to factory defaults
                            OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset = 0;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset = 0;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset = 0;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset = 0;
                            
                            // RESET DEL CENTRO OTTICO USANDO I VALORI DINAMICI
                            OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = CENTER_X;
                            OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = CENTER_Y;
                            
                            SetMode(FW_Const::GunMode_Calibration);
                            AbsMouse5.move(32768/2, 32768/2);
                            AbsMouse5.report();

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
                analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0);  // [ESP32_PORT] per EPS32
            #else // rp2040
            digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
            #endif
            delay(50);
            analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength]);
            delay(125);
            #ifdef ARDUINO_ARCH_ESP32
                analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0);  // [ESP32_PORT] per ESP32
            #else // rp2040            
            digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
            #endif
        }
    #endif // USES_RUMBLE

    // The trigger has confirmed the new calibration. A lost final ACK must
    // not undo values that the App may already have accepted.
    calStage = FW_Const::Cali_Verify + 1;
    if(fromDesktop && !OF_Serial::AppSerialSendResponse(OF_Const::sCaliStageUpd, &calStage, 1))
        OF_Serial::AppSerialSendError();
    return;

calibration_failed:
    communicationFailed = true;

calibration_cancelled:
    // One restore path for user cancellation and failed intermediate updates.
    OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset = _topOffset;
    OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset = _bottomOffset;
    OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset = _leftOffset;
    OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset = _rightOffset;
    OF_Prefs::profiles[OF_Prefs::currentProfile].TLled = _TLled;
    OF_Prefs::profiles[OF_Prefs::currentProfile].TRled = _TRled;
    OF_Prefs::profiles[OF_Prefs::currentProfile].adjX = _adjX;
    OF_Prefs::profiles[OF_Prefs::currentProfile].adjY = _adjY;
    OpenFIREper.source(_adjX, _adjY);
    OpenFIREper.deinit(0);
    stateFlags |= FW_Const::StateFlag_PrintSelectedProfile;

    SetMode(fromDesktop ? FW_Const::GunMode_Docked : FW_Const::GunMode_Run);

    if(fromDesktop) {
        if(!communicationFailed) {
            // Reuse existing messages: reset the App's provisional values
            // before End, so it cannot report a cancelled calibration as saved.
            calStage = FW_Const::Cali_Init;
            if(!OF_Serial::AppSerialSendResponse(OF_Const::sCaliStageUpd, &calStage, 1))
                communicationFailed = true;
            else {
                calStage = FW_Const::Cali_Verify + 1;
                if(!OF_Serial::AppSerialSendResponse(OF_Const::sCaliStageUpd, &calStage, 1))
                    communicationFailed = true;
            }
        }

        if(communicationFailed)
            OF_Serial::AppSerialSendError();
    }
}

void FW_Common::GetPosition()
{
    const CameraProfile& cameraProfile = OpenFIRECamera::Profile();
    const int CAM_COORD_RES_X = cameraProfile.mouseResX;
    const int CAM_COORD_RES_Y = cameraProfile.mouseResY;

    // Perspective code uses 2 extra precision bits (res_x/res_y are << 2).
    constexpr int PERSPECTIVE_SCALE = 1 << 2;
    constexpr int SCREEN_RES_X = res_x / PERSPECTIVE_SCALE;
    constexpr int SCREEN_RES_Y = res_y / PERSPECTIVE_SCALE;
    const int CAM_TEST_WIDTH = cameraProfile.testWidth;
    const int CAM_TEST_HEIGHT = cameraProfile.testHeight;
    const int CAM_TEST_OFFSET_X = cameraProfile.testOffsetX;
    const int CAM_TEST_OFFSET_Y = cameraProfile.testOffsetY;

    // Target aspect ratio used by Serial AR correction: 4:3 content.
    constexpr int AR_CORRECTION_W = 4;
    constexpr int AR_CORRECTION_H = 3;

    if(OpenFIRECamera::IsReady()) {
        int error = OpenFIRECamera::Read();
        if(error == OpenFIRECamera::Error_Success) {
           
            // if diamond layout, or square
            if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) { // layoutDiamond = 1
                OpenFIREdiamond.begin(OpenFIRECamera::XPositions(), OpenFIRECamera::YPositions(), OpenFIRECamera::Seen());

                OpenFIREper.warp(OpenFIREdiamond.X(0), OpenFIREdiamond.Y(0),
                                 OpenFIREdiamond.X(1), OpenFIREdiamond.Y(1),
                                 OpenFIREdiamond.X(2), OpenFIREdiamond.Y(2),
                                 OpenFIREdiamond.X(3), OpenFIREdiamond.Y(3),
                                 res_x / 2, 0, 0,
                                 res_y / 2, res_x / 2,
                                 res_y, res_x, res_y / 2);
            } else { // layoutSquare = 0
                OpenFIREsquare.begin(OpenFIRECamera::XPositions(), OpenFIRECamera::YPositions(), OpenFIRECamera::Seen());               
               
                #ifdef USE_MULTI_ONE_EURO_FILTER
                    X_in[0] = OpenFIREsquare.X(0);
                    Y_in[0] = OpenFIREsquare.Y(0);
                    X_in[1] = OpenFIREsquare.X(1);
                    Y_in[1] = OpenFIREsquare.Y(1);
                    X_in[2] = OpenFIREsquare.X(2);
                    Y_in[2] = OpenFIREsquare.Y(2);
                    X_in[3] = OpenFIREsquare.X(3);
                    Y_in[3] = OpenFIREsquare.Y(3);          

                    oef_multi.process(X_in, Y_in, X_out, Y_out);

                    OpenFIREper.warp(X_out[0], Y_out[0],
                                 X_out[1], Y_out[1],
                                 X_out[2], Y_out[2],
                                 X_out[3], Y_out[3],
                                 OF_Prefs::profiles[OF_Prefs::currentProfile].TLled, 0,
                                 OF_Prefs::profiles[OF_Prefs::currentProfile].TRled, 0,
                                 OF_Prefs::profiles[OF_Prefs::currentProfile].TLled, res_y,
                                 OF_Prefs::profiles[OF_Prefs::currentProfile].TRled, res_y);
   
                #else
                    OpenFIREper.warp(OpenFIREsquare.X(0), OpenFIREsquare.Y(0),
                                 OpenFIREsquare.X(1), OpenFIREsquare.Y(1),
                                 OpenFIREsquare.X(2), OpenFIREsquare.Y(2),
                                 OpenFIREsquare.X(3), OpenFIREsquare.Y(3),
                                 OF_Prefs::profiles[OF_Prefs::currentProfile].TLled, 0,
                                 OF_Prefs::profiles[OF_Prefs::currentProfile].TRled, 0,
                                 OF_Prefs::profiles[OF_Prefs::currentProfile].TLled, res_y,
                                 OF_Prefs::profiles[OF_Prefs::currentProfile].TRled, res_y);
                #endif // USE_MULTI_ONE_EURO_FILTER

            }

            // Output mapped to screen resolution because offsets are measured in pixels
            mouseX = map(OpenFIREper.getX(), 0, res_x, (0 - OF_Prefs::profiles[OF_Prefs::currentProfile].leftOffset), (res_x + OF_Prefs::profiles[OF_Prefs::currentProfile].rightOffset));                 
            mouseY = map(OpenFIREper.getY(), 0, res_y, (0 - OF_Prefs::profiles[OF_Prefs::currentProfile].topOffset), (res_y + OF_Prefs::profiles[OF_Prefs::currentProfile].bottomOffset));         

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
                    case OF_Const::ar16_9: {
                        // Fit centered 4:3 content inside a 16:9 display.
                        constexpr int correctedWidth =
                            (res_x * AR_CORRECTION_W * 9) /
                            (AR_CORRECTION_H * 16);

                        constexpr int correctionX =
                            (res_x - correctedWidth) / 2;

                        conMoveX = map(conMoveX, correctionX, res_x - correctionX, 0, 32767);
                        conMoveX = constrain(conMoveX, 0, 32767);
                        conMoveY = map(conMoveY, 0, res_y, 0, 32767);
                        break;
                    }

                    case OF_Const::ar16_10: {
                        // Fit centered 4:3 content inside a 16:10 display.
                        constexpr int correctedWidth =
                            (res_x * AR_CORRECTION_W * 10) /
                            (AR_CORRECTION_H * 16);

                        constexpr int correctionX =
                            (res_x - correctedWidth) / 2;

                        conMoveX = map(conMoveX, correctionX, res_x - correctionX, 0, 32767);
                        conMoveX = constrain(conMoveX, 0, 32767);
                        conMoveY = map(conMoveY, 0, res_y, 0, 32767);
                        break;
                    }

                    case OF_Const::ar3_2: {
                        // Fit centered 4:3 content inside a 3:2 display.
                        constexpr int correctedWidth =
                            (res_x * AR_CORRECTION_W * 2) /
                            (AR_CORRECTION_H * 3);

                        constexpr int correctionX =
                            (res_x - correctedWidth) / 2;

                        conMoveX = map(conMoveX, correctionX, res_x - correctionX, 0, 32767);
                        conMoveX = constrain(conMoveX, 0, 32767);
                        conMoveY = map(conMoveY, 0, res_y, 0, 32767);
                        break;
                    }

                    case OF_Const::ar5_4: {
                        // A 5:4 display is narrower than 4:3, so the correction
                        // is applied vertically instead of horizontally.
                        constexpr int correctedHeight =
                            (res_y * AR_CORRECTION_H * 5) /
                            (AR_CORRECTION_W * 4);

                        constexpr int correctionY =
                            (res_y - correctedHeight) / 2;

                        conMoveX = map(conMoveX, 0, res_x, 0, 32767);
                        conMoveY = map(conMoveY, correctionY, res_y - correctionY, 0, 32767);
                        conMoveY = constrain(conMoveY, 0, 32767);
                        break;
                    }

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

                if(buttons.analogOutput)
                     Gamepad16.moveCam(conMoveX, conMoveY);
                else AbsMouse5.move(conMoveX, conMoveY);

            } else {
                if(gunMode == FW_Const::GunMode_Verification) {
                    // Output mapped to Mouse resolution
                    conMoveX = map(conMoveX, 0, res_x, 0, 32767);
                    conMoveY = map(conMoveY, 0, res_y, 0, 32767);

                    AbsMouse5.move(conMoveX, conMoveY);
                    AbsMouse5.report();
                }

                if(millis() - testLastStamp > 50) {
                    testLastStamp = millis();
                    // RAW Camera Output mapped to screen res (1920x1080)
                    // Screen resolution is now dynamically derived from res_x/res_y.
                    int rawX[4];
                    int rawY[4];
                    bool outsideFov[4];

                    // RAW Output for viewing in processing sketch mapped to 1920x1080 screen resolution
                    for (int i = 0; i < 4; ++i) {
                        int pointX;
                        int pointY;

                        if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) {
                            pointX = OpenFIREdiamond.X(i);
                            pointY = OpenFIREdiamond.Y(i);

                            rawX[i] = map(pointX, 0, CAM_COORD_RES_X,
                                          CAM_TEST_OFFSET_X + CAM_TEST_WIDTH, CAM_TEST_OFFSET_X);
                            rawY[i] = map(pointY, 0, CAM_COORD_RES_Y,
                                          CAM_TEST_OFFSET_Y, CAM_TEST_OFFSET_Y + CAM_TEST_HEIGHT);
                        } else {
                            pointX = OpenFIREsquare.X(i);
                            pointY = OpenFIREsquare.Y(i);

                            rawX[i] = map(pointX, 0, CAM_COORD_RES_X,
                                          CAM_TEST_OFFSET_X, CAM_TEST_OFFSET_X + CAM_TEST_WIDTH);
                            rawY[i] = map(pointY, 0, CAM_COORD_RES_Y,
                                          CAM_TEST_OFFSET_Y, CAM_TEST_OFFSET_Y + CAM_TEST_HEIGHT);
                        }

                        outsideFov[i] =
                            pointX < 0 || pointX > cameraProfile.mouseMaxX ||
                            pointY < 0 || pointY > cameraProfile.mouseMaxY;
                    }

                    if(runMode == FW_Const::RunMode_Processing) {
                        int mouseXscaled = mouseX / PERSPECTIVE_SCALE;
                        int mouseYscaled = mouseY / PERSPECTIVE_SCALE;

                        // Encode outside-FOV state in bit 0 of each transmitted X coordinate.
                        // Even X = inside FOV, odd X = outside FOV. rawX/rawY stay unchanged.
                        int serialX[4];
                        for(int i = 0; i < 4; ++i)
                            serialX[i] = rawX[i] * 2 + (outsideFov[i] ? 1 : 0);

                        if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout) {
                            int testMedianX = map(OpenFIREdiamond.testMedianX(), 0, CAM_COORD_RES_X,
                                                  CAM_TEST_OFFSET_X + CAM_TEST_WIDTH, CAM_TEST_OFFSET_X);
                            int testMedianY = map(OpenFIREdiamond.testMedianY(), 0, CAM_COORD_RES_Y,
                                                  CAM_TEST_OFFSET_Y, CAM_TEST_OFFSET_Y + CAM_TEST_HEIGHT);
                            uint8_t payload[sizeof(int) * 12];
                            memcpy(&payload[0],  &serialX[0],   sizeof(int));
                            memcpy(&payload[4],  &rawY[0],      sizeof(int));
                            memcpy(&payload[8],  &serialX[1],   sizeof(int));
                            memcpy(&payload[12], &rawY[1],      sizeof(int));
                            memcpy(&payload[16], &serialX[2],   sizeof(int));
                            memcpy(&payload[20], &rawY[2],      sizeof(int));
                            memcpy(&payload[24], &serialX[3],   sizeof(int));
                            memcpy(&payload[28], &rawY[3],      sizeof(int));
                            memcpy(&payload[32], &mouseXscaled, sizeof(int));
                            memcpy(&payload[36], &mouseYscaled, sizeof(int));
                            memcpy(&payload[40], &testMedianX,  sizeof(int));
                            memcpy(&payload[44], &testMedianY,  sizeof(int));
                            OF_Serial::AppSerialSendEvent(OF_Const::sTestCoords, payload, sizeof(payload));
                        } else {
                            int testMedianX = map(OpenFIREsquare.testMedianX(), 0, CAM_COORD_RES_X,
                                                  CAM_TEST_OFFSET_X, CAM_TEST_OFFSET_X + CAM_TEST_WIDTH);
                            int testMedianY = map(OpenFIREsquare.testMedianY(), 0, CAM_COORD_RES_Y,
                                                  CAM_TEST_OFFSET_Y, CAM_TEST_OFFSET_Y + CAM_TEST_HEIGHT);
                            uint8_t payload[sizeof(int) * 12];
                            memcpy(&payload[0],  &serialX[0],   sizeof(int));
                            memcpy(&payload[4],  &rawY[0],      sizeof(int));
                            memcpy(&payload[8],  &serialX[1],   sizeof(int));
                            memcpy(&payload[12], &rawY[1],      sizeof(int));
                            memcpy(&payload[16], &serialX[2],   sizeof(int));
                            memcpy(&payload[20], &rawY[2],      sizeof(int));
                            memcpy(&payload[24], &serialX[3],   sizeof(int));
                            memcpy(&payload[28], &rawY[3],      sizeof(int));
                            memcpy(&payload[32], &mouseXscaled, sizeof(int));
                            memcpy(&payload[36], &mouseYscaled, sizeof(int));
                            memcpy(&payload[40], &testMedianX,  sizeof(int));
                            memcpy(&payload[44], &testMedianY,  sizeof(int));
                            OF_Serial::AppSerialSendEvent(OF_Const::sTestCoords, payload, sizeof(payload));
                        }
                    }

                    #ifdef USES_DISPLAY
                        /*
                        for(int i = 0; i < 4; ++i) {
                            rawX[i] = map(rawX[i],
                                CAM_TEST_OFFSET_X,
                                CAM_TEST_OFFSET_X + CAM_TEST_WIDTH,
                                0, SCREEN_RES_X);

                            rawY[i] = map(rawY[i],
                                CAM_TEST_OFFSET_Y,
                                CAM_TEST_OFFSET_Y + CAM_TEST_HEIGHT,
                                0, SCREEN_RES_Y);
                        }
                        */
                        OLED.DrawVisibleIR(rawX, rawY);
                    #endif // USES_DISPLAY
                }
            }
        } else if(error != OpenFIRECamera::Error_DataMismatch)
            PrintIrError();
    } else PrintIrError();
}

void FW_Common::PrintIrError()
{
    // set flag to warn desktop app when docking
    const bool firstError = !camNotAvailable;
    if(firstError)
        camNotAvailable = true;

    const uint8_t error = OF_Const::sErrCam;
    if(firstError && OF_Serial::AppSerialSendEvent(OF_Const::sError, &error, 1)) {
        camWarningTimestamp = millis();
        return;
    }

    if(millis() - camWarningTimestamp > CAM_WARNING_INTERVAL) {
        if(!OF_Serial::AppSerialSendEvent(OF_Const::sError, &error, 1))
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

    SetMode(gunMode);

    OpenFIREper.source(OF_Prefs::profiles[profile].adjX, OF_Prefs::profiles[profile].adjY);                                                          
    OpenFIREper.deinit(0);

    // set IR sensitivity
    if(OF_Prefs::profiles[profile].irSens <= OpenFIRECamera::Sensitivity_Max)
        SetIrSensitivity((OpenFIRECamera::Sensitivity_e)OF_Prefs::profiles[profile].irSens);

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
    if(sensitivity > OpenFIRECamera::Sensitivity_Max)
        return;

    if(OF_Prefs::profiles[OF_Prefs::currentProfile].irSens != sensitivity) {
        OF_Prefs::profiles[OF_Prefs::currentProfile].irSens = sensitivity;
        stateFlags |= FW_Const::StateFlag_SavePreferencesEn;
    }

    OpenFIRECamera::SetSensitivity((OpenFIRECamera::Sensitivity_e)sensitivity);

}

void FW_Common::SetIrLayout(const int &layout)
{
    // TODO: we need an enum for layout types available
    if(layout >= OF_Const::layoutTypes)
        return;

    if(OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout != layout) {
        OF_Prefs::profiles[OF_Prefs::currentProfile].irLayout = layout;

        SetMode(gunMode);

        OpenFIREper.deinit(0);

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

    int saveResult = OF_Prefs::SaveProfiles();

    if(saveResult == OF_Prefs::Error_Success) {
        int result = OF_Prefs::SaveToggles();
        if(saveResult == OF_Prefs::Error_Success &&
           result != OF_Prefs::Error_Success)
             saveResult = result;

        if(OF_Prefs::toggles[OF_Const::customPins]) {
            result = OF_Prefs::SavePins();
            if(saveResult == OF_Prefs::Error_Success &&
               result != OF_Prefs::Error_Success)
                saveResult = result;
        }

        result = OF_Prefs::SaveSettings();
        if(saveResult == OF_Prefs::Error_Success &&
           result != OF_Prefs::Error_Success)
            saveResult = result;

        result = OF_Prefs::SaveButtons();
        if(saveResult == OF_Prefs::Error_Success &&
           result != OF_Prefs::Error_Success)
            saveResult = result;

        result = OF_Prefs::SaveUSBID();
        if(saveResult == OF_Prefs::Error_Success &&
           result != OF_Prefs::Error_Success)
            saveResult = result;
    }

    // During an App commit the pin map in RAM may not be active yet.
    // Report the result through the protocol; do not drive LEDs/OLED here.
    if(dockedSaving)
        return saveResult;

    if(saveResult == OF_Prefs::Error_Success) {

        #ifdef USES_DISPLAY
            OLED.ScreenModeChange(ExtDisplay::Screen_SaveSuccess);
        #endif // USES_DISPLAY

        if(gunMode != FW_Const::GunMode_Docked)
            Serial.println("Settings saved to Flash"), Serial.flush();
        
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

        SetMode(gunMode);

        return OF_Prefs::Error_Success;
    } else {
        #ifdef USES_DISPLAY
            OLED.ScreenModeChange(ExtDisplay::Screen_SaveError);
        #endif // USES_DISPLAY

        // TODO: reimpl a detailed error string
        if(gunMode != FW_Const::GunMode_Docked)
            Serial.println("Error saving Preferences to Flash.");

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

        return saveResult;

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


// ============ [ESP32_PORT] ============
// Restoration of Serial after definition for serial connections / ripristino di Serial dopo definizione per connessione seriali
#ifdef OPENFIRE_WIRELESS_ENABLE
    #undef Serial
    #ifdef AUX_SERIAL
        #define Serial AUX_SERIAL
        #undef AuxSerial
    #endif
#endif // OPENFIRE_WIRELESS_ENABLE
// ============ [ESP32_PORT] ============
// End restoration of Serial after definition for serial connections / fine ripristino di Serial dopo definizione per connessione seriali
