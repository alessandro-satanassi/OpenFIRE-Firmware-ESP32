 /*!
 * @file OpenFIREserial.cpp
 * @brief Serial RX buffer reading routines.
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

#include "OpenFIREserial.h"
#include "OpenFIREprefs.h"
#include "OpenFIREFeedback.h"
#include "OpenFIRElights.h"
#include "boards/OpenFIREshared.h"
#include "OpenFIREcommon.h"

// ============ [ESP32_PORT] ============
// redefinition of Serial to handle wireless serial connections / redifinizione di Serial per gestire le connessione wireless seriali
#ifdef OPENFIRE_WIRELESS_ENABLE
    extern Stream* Serial_OpenFIRE_Stream;
    #ifdef Serial
        #define AUX_SERIAL Serial
        #undef Serial
    #endif
    #define Serial (*Serial_OpenFIRE_Stream)
#endif // OPENFIRE_WIRELESS_ENABLE
// ============ [ESP32_PORT] ============
// END redefinition of Serial to handle wireless serial connections / fine redifinizione di Serial per gestire le connessione wireless seriali ========

#ifdef ARDUINO_ARCH_ESP32  // [ESP32_PORT]
    #define delay(ms) vTaskDelay(pdMS_TO_TICKS(ms))                    
#endif //ARDUINO_ARCH_ESP32


#ifdef MAMEHOOKER
void OF_Serial::SerialProcessing()
{
    // For more info about Serial commands, see the OpenFIRE repo wiki.

    switch(Serial.read()) {
        // Start Signal
        case 'S':
          if(serialMode)
              Serial.println("SERIALREAD: Detected Serial Start command while already in Serial handoff mode!");
          // TODO: handle the other Start line bits? this assume equiv of `S6`
          else {
              serialMode = true;
              OF_FFB::FFBShutdown();

              #ifdef USES_SOLENOID
              serialSolCustomHoldLength = 0;
              serialSolCustomPauseLength = 0;
              #endif // USES_SOLENOID

              #ifdef USES_RUMBLE
              serialRumbCustomHoldLength = 0;
              serialRumbCustomPauseLength = 0;
              #endif // USES_RUMBLE

              #ifdef LED_ENABLE
                  // Set the LEDs to a mid-intense white.
                  OF_RGB::LedUpdate(127, 127, 127);
              #endif // LED_ENABLE

              #ifdef USES_DISPLAY
                  // init basic display to show mamehook icon
                  if(FW_Common::gunMode == FW_Const::GunMode_Run)
                      FW_Common::OLED.ScreenModeChange(ExtDisplay::Screen_Mamehook_Single, FW_Common::buttons.analogOutput);
              #endif // USES_DISPLAY
          }
          break;
        // Modesetting Signal
        case 'M':
          switch(Serial.read()) {
              // input mode
              case '0':
                Serial.read(); // nomf
                switch(Serial.read()) {
                    case '2': // "hybrid" - just use the default m&kb mode
                    case '0': // mouse & kb 
                      FW_Common::buttons.analogOutput = false;
                      break;
                    // gamepad
                    case '1':
                      FW_Common::buttons.analogOutput = true;
                      Gamepad16.stickRight = (Serial.peek() == 'L') ? true: false;
                      break;
                    // official "MiSTer optimized" mode
                    case '9':
                      FW_Common::buttons.analogOutput = true;
                      Gamepad16.stickRight = true;
                      // HACK SHACK - testing MiSTer-friendly default gamepad maps
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Trigger].reportCode3 = PAD_A,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_A].reportCode3       = PAD_B,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_B].reportCode3       = PAD_X,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Reload].reportCode3  = PAD_Y,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Start].reportCode3   = PAD_START,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Select].reportCode3  = PAD_SELECT,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Up].reportCode3      = PAD_UP,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Down].reportCode3    = PAD_DOWN,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Left].reportCode3    = PAD_LEFT,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Right].reportCode3   = PAD_RIGHT,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Pedal].reportCode3   = PAD_LB,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Pedal2].reportCode3  = PAD_RB,
                      LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Pump].reportCode3    = PAD_C;
                      #ifdef USES_DISPLAY
                          FW_Common::OLED.mister = true;
                      #endif // USES_DISPLAY
                      break;
                }
                FW_Common::buttons.ReleaseAll();
                #ifdef USES_DISPLAY
                    if(!serialMode && FW_Common::gunMode == FW_Const::GunMode_Run)
                        FW_Common::OLED.ScreenModeChange(ExtDisplay::Screen_Normal, FW_Common::buttons.analogOutput);
                    else if(serialMode && FW_Common::gunMode == FW_Const::GunMode_Run &&
                            FW_Common::OLED.serialDisplayType > ExtDisplay::ScreenSerial_None &&
                            FW_Common::OLED.serialDisplayType < ExtDisplay::ScreenSerial_Both) {
                        FW_Common::OLED.ScreenModeChange(ExtDisplay::Screen_Mamehook_Single, FW_Common::buttons.analogOutput);
                    }
                #endif // USES_DISPLAY
                break;
              // offscreen button mode
              case '1':
                Serial.read(); // nomf
                switch(Serial.read()) {
                    // cursor in bottom left - just use disabled
                    case '1':
                    // "true offscreen shot" mode - just use disabled for now
                    case '3':
                    // disabled
                    case '0':
                      memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Trigger].reportType2,
                             &OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_Trigger][2],
                             sizeof(LightgunButtons::Desc_s::reportType)*2);
                      // reset remapping for low button users if M1x2 was previously called
                      if(OF_Prefs::toggles[OF_Const::lowButtonsMode])
                          memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_A].reportType,
                                 OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_A],
                                 sizeof(LightgunButtons::Desc_s::reportType)*2);
                      serialMappingsOffscreenShot = false;
                      break;
                    // offscreen button
                    case '2':
                      memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Trigger].reportType2,
                             OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_A],
                             sizeof(LightgunButtons::Desc_s::reportType)*2);
                      // remap bindings for low button users to make e.g. VCop 3 playable with 1 btn + pedal
                      if(OF_Prefs::toggles[OF_Const::lowButtonsMode])
                          memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_A].reportType,
                                 OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_Reload],
                                 sizeof(LightgunButtons::Desc_s::reportType)*2);
                      serialMappingsOffscreenShot = true;
                      break;
                }
                FW_Common::UpdateBindings(false);
                break;
              // pedal functionality
              case '2':
                Serial.read();                                         // nomf
                switch(Serial.read()) {
                    // separate button (default to original binds)
                    case '0':
                      memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Pedal].reportType,
                             OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_Pedal],
                             sizeof(OF_Prefs::backupButtonDesc[0]));
                      serialMappingsPedalMode = 0;
                      break;
                    // make reload button (mapping of Button A)
                    case '1':
                      memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Pedal].reportType,
                             OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_A],
                             sizeof(OF_Prefs::backupButtonDesc[0]));
                      serialMappingsPedalMode = 1;
                      break;
                    // make middle mouse button (mapping of Button B, useful for low buttons mode & e.g. using VCop3 ES mode)
                    case '2':
                      memcpy(&LightgunButtons::ButtonDesc[FW_Const::BtnIdx_Pedal].reportType,
                             OF_Prefs::backupButtonDesc[FW_Const::BtnIdx_B],
                             sizeof(OF_Prefs::backupButtonDesc[0]));
                      serialMappingsPedalMode = 2;
                      break;
                }
                FW_Common::UpdateBindings(false);
                break;
              // aspect ratio correction
              case '3':
                Serial.read(); // nomf
                serialARcorrection = Serial.read() - '0';
                if(!serialMode) {
                    if(serialARcorrection) { Serial.println("Setting 4:3 correction on!"); }
                    else { Serial.println("Setting 4:3 correction off!"); }
                }
                break;
              #ifdef USES_TEMP
              // temp sensor disabling (why?)
              case '4':
                Serial.read(); // nomf
                Serial.read();
                break;
              #endif // USES_TEMP
              // autoreload (TODO: maybe?)
              case '5':
                Serial.read(); // nomf
                Serial.read();
                break;
              // rumble only mode (enable Rumble FF)
              case '6':
                Serial.read(); // nomf
                switch(Serial.read()) {
                    // disable
                    case '0':
                      if(OF_Prefs::pins[OF_Const::solenoidSwitch] == -1 && OF_Prefs::pins[OF_Const::solenoidPin] >= 0)
                          OF_Prefs::toggles[OF_Const::solenoid] = true;
                      if(OF_Prefs::pins[OF_Const::rumblePin] >= 0) 
                          OF_Prefs::toggles[OF_Const::rumbleFF] = false;
                      break;
                    // enable
                    case '1':
                      if(OF_Prefs::pins[OF_Const::rumbleSwitch] == -1 && OF_Prefs::pins[OF_Const::rumblePin] >= 0) { OF_Prefs::toggles[OF_Const::rumble] = true; }
                      if(OF_Prefs::pins[OF_Const::solenoidSwitch] == -1 && OF_Prefs::pins[OF_Const::solenoidPin] >= 0) { OF_Prefs::toggles[OF_Const::solenoid] = false; }
                      if(OF_Prefs::pins[OF_Const::rumblePin] >= 0) { OF_Prefs::toggles[OF_Const::rumbleFF] = true; }
                      break;
                }
                OF_FFB::FFBShutdown();
                break;
              #ifdef USES_SOLENOID
              // solenoid automatic mode
              case '8':
                Serial.read(); // Nomf the padding bit.
                switch(Serial.read()) {
                // "auto"
                case '1':
                    OF_FFB::burstFireActive = true;
                    OF_Prefs::toggles[OF_Const::autofire] = false;
                    break;
                // "always on"
                case '2':
                    OF_Prefs::toggles[OF_Const::autofire] = true;
                    OF_FFB::burstFireActive = false;
                    break;
                // disabled
                case '0':
                    OF_Prefs::toggles[OF_Const::autofire] = false;
                    OF_FFB::burstFireActive = false;
                    break;
                }
                break;
              #endif // USES_SOLENOID
              #ifdef USES_DISPLAY
              case 'D':
                Serial.read(); // Nomf padding byte
                switch(Serial.read()) {
                    case '0':
                      FW_Common::OLED.serialDisplayType = ExtDisplay::ScreenSerial_None;
                      break;
                    case '1':
                      FW_Common::OLED.serialDisplayType = ExtDisplay::ScreenSerial_Life;
                      break;
                    case '2':
                      FW_Common::OLED.serialDisplayType = ExtDisplay::ScreenSerial_Ammo;
                      break;
                    case '3':
                      FW_Common::OLED.serialDisplayType = ExtDisplay::ScreenSerial_Both;
                      break;
                }
                
                if(Serial.read() == 'B') {
                    FW_Common::OLED.lifeBar = true;
		            FW_Common::dispMaxLife = 0;
                } else FW_Common::OLED.lifeBar = false;

                // prevent glitching if currently in pause mode
                if(FW_Common::gunMode == FW_Const::GunMode_Run) {
                    if(FW_Common::OLED.serialDisplayType == ExtDisplay::ScreenSerial_Both)
                        FW_Common::OLED.ScreenModeChange(ExtDisplay::Screen_Mamehook_Dual);
                    else if(FW_Common::OLED.serialDisplayType > ExtDisplay::ScreenSerial_None)
                        FW_Common::OLED.ScreenModeChange(ExtDisplay::Screen_Mamehook_Single, FW_Common::buttons.analogOutput);
                }
                break;
              #endif // USES_DISPLAY
              default:
                if(!serialMode) Serial.println("SERIALREAD: Serial modesetting command found, but no valid set bit found!");
                break;
          }
          break;
        // End Signal
        // Check to make sure that 'E' is not actually a glitched command bit
        // by ensuring that there's no adjacent bit.
        case 'E':
          // Either "no buffer" (-1) or any of the non-visual control bits, mainly Carriage Returns/Newlines
          // (as Windows sends these implicitly in `echo` commands by default)
          if(Serial.peek() <= 32) {
              if(!serialMode) Serial.println("SERIALREAD: Detected Serial End command while Serial Handoff mode is already off!");
              else {
                  serialMode = false;
                  memset(serialQueue, false, sizeof(serialQueue));
                  serialARcorrection = false;
                  serialMappingsOffscreenShot = false;
                  serialMappingsPedalMode = 0;
                  #ifdef USES_DISPLAY
                      FW_Common::OLED.serialDisplayType = ExtDisplay::ScreenSerial_None;
                      if(FW_Common::gunMode == FW_Const::GunMode_Run) FW_Common::OLED.ScreenModeChange(ExtDisplay::Screen_Normal, FW_Common::buttons.analogOutput);
                  #endif // USES_DISPLAY
                  #ifdef LED_ENABLE
                      // Clear any stale serial LED pulses
                      serialLEDPulseColorMap = 0b00000000;
                      serialLEDPulses = 0;
                      serialLEDPulsesLast = 0;
                      serialLEDPulseRising = true;
                      serialLEDR = 0;
                      serialLEDG = 0;
                      serialLEDB = 0;
                      serialLEDChange = false;
                      if(FW_Common::gunMode == FW_Const::GunMode_Run)
                        OF_RGB::LedOff(); // Turn it off, and let lastSeen handle it from here.
                  #endif // LED_ENABLE
                  #ifdef USES_RUMBLE
                      #ifdef ARDUINO_ARCH_ESP32
                        analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0); // [ESP32_PORT] per ESP32
                      #else // rp2040
                      digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
                      #endif
                      serialRumbPulseStage = 0;
                      serialRumbPulses = 0;
                      serialRumbPulsesLast = 0;
                      serialRumbCustomHoldLength = 0;
                      serialRumbCustomPauseLength = 0;
                  #endif // USES_RUMBLE
                  #ifdef USES_SOLENOID
                      OF_FFB::SetSolenoid(LOW);
                      serialSolPulses = 0;
                      serialSolPulsesLast = 0;
                      serialSolCustomHoldLength = 0;
                      serialSolCustomPauseLength = 0;
                  #endif // USES_SOLENOID
                  FW_Common::buttons.ReleaseAll();
                  // remap back to defaults, in case they were changed
                  FW_Common::UpdateBindings(true);
                  Serial.println("Received end serial pulse, releasing FF override.");
              }
              break;
          }
          break;
        // owo SPECIAL SETUP EH?
        case 'X':
          switch(Serial.read()) {
              // Set Autofire Interval Length
              case 'I':
                OF_FFB::autofireDoubleLengthWait = Serial.read() - '0';
                break;
              // Remap player numbers
              case 'R':
              {
                char serialInput = Serial.read();
                if(serialInput >= '1' && serialInput <= '4') {
                    FW_Common::playerStartBtn = serialInput;
                    FW_Common::playerSelectBtn = serialInput + 4;
                    FW_Common::UpdateBindings(false);
                } else Serial.println("SERIALREAD: Player remap command called, but an invalid or no slot number was declared!");
                break;
              }
              default:
                Serial.println("SERIALREAD: Internal setting cmd detected, but not valid!");
                Serial.println("Internally recognized commands are:");
                Serial.println("I(nterval Autofire)2/3/4 / R(emap)1/2/3/4");
                break;
          }
          // End of 'X'
          break;
        // Enter Docked Mode
        case OF_Const::sDock1:
          if(Serial_available(1) && Serial.read() == OF_Const::sDock2) {
            AppSerialSessionBegin();
            #if /*defined(ARDUINO_ARCH_RP2040) &&*/ defined(DUAL_CORE) // This may be being run from Core 1, so signal if running in main Run Mode.
            if(FW_Common::gunMode == FW_Const::GunMode_Run) {
                #ifdef ARDUINO_ARCH_ESP32 
                esp32_fifo.push(FW_Const::GunMode_Docked);
                esp32_fifo.pop();
                #else //rp2040
                rp2040.fifo.push(FW_Const::GunMode_Docked);
                rp2040.fifo.pop();
                #endif
            }
            else FW_Common::SetMode(FW_Const::GunMode_Docked);
            #else
            FW_Common::SetMode(FW_Const::GunMode_Docked);
            #endif // DUAL_CORE
          }
          break;
        // Force Feedback
        case 'F':
          switch(Serial.read()) {
              #ifdef USES_SOLENOID
              // Solenoid bits
              case '0':
                Serial.read(); // nomf the padding
                switch(Serial.read()) {
                // Solenoid "on" command
                case '1':
                    #ifdef USES_TEMP
                    // block every other signal if temp is at warning threshold
                    if(OF_Prefs::pins[OF_Const::tempPin] > -1) {
                        switch(OF_FFB::tempStatus) {
                        case OF_FFB::Temp_Safe:
                            serialQueue[SerialQueue_Solenoid] = true;
                            break;
                        case OF_FFB::Temp_Warning:
                            if(serialSolTempBuffer) serialQueue[SerialQueue_Solenoid] = true;
                            serialSolTempBuffer = !serialSolTempBuffer;
                            break;
                        case OF_FFB::Temp_Fatal:
                        default:
                            break;
                        }
                    } else 
                    #endif // USES_TEMP
                    serialQueue[SerialQueue_Solenoid] = true;
                    break;
                // Solenoid "pulse" command (only if not already pulsing)
                case '2':
                    if(!serialQueue[SerialQueue_SolPulse]) {
                        Serial.read(); // nomf the padding bit.
                        if(Serial.peek() >= '0' & Serial.peek() <= '9') {
                            serialQueue[SerialQueue_SolPulse] = true;
                            char serialInputS[4] = {0,0,0,0};
                            for(uint n = 0; n < 3; ++n) {
                                serialInputS[n] = Serial.read();
                                if(Serial.peek() < '0' || Serial.peek() > '9')
                                    break;
                            }
                            serialSolPulses = atoi(serialInputS);
                            if(!serialSolPulses) serialSolPulses++;
                            serialSolPulsesLast = 0;
                        }
                    }
                    break;
                // Solenoid "off" command
                case '0':
                    serialQueue[SerialQueue_Solenoid] = false, serialQueue[SerialQueue_SolPulse] = false;
                    break;
                }
                break;
              #endif // USES_SOLENOID
              #ifdef USES_RUMBLE
              // Rumble bits
              case '1':
                Serial.read(); // nomf the padding
                switch(Serial.read()) {
                // Rumble "on" command
                case '1':
                    serialQueue[SerialQueue_Rumble] = true;
                    break;
                // Rumble "pulse" command (only if not already pulsing)
                case '2':
                    if(!serialQueue[SerialQueue_RumbPulse]) {
                        Serial.read(); // nomf the padding
                        if(Serial.peek() >= '0' && Serial.peek() <= '9') {
                            serialQueue[SerialQueue_RumbPulse] = true;
                            char serialInputS[4] = {0,0,0,0};
                            for(uint n = 0; n < 3; ++n) {
                                serialInputS[n] = Serial.read();
                                if(Serial.peek() < '0' || Serial.peek() > '9')
                                    break;
                            }
                            serialRumbPulses = atoi(serialInputS);
                            if(!serialRumbPulses) serialRumbPulses++;
                            serialRumbPulsesLast = 0;
                        }
                    }
                    break;
                // Rumble "off" command
                case '0':
                    serialQueue[SerialQueue_Rumble] = false, serialQueue[SerialQueue_RumbPulse] = false;
                    break;
                }
                break;
              #endif // USES_RUMBLE
              #ifdef LED_ENABLE
              // LED Red bits
              case '2':
                Serial.read(); // nomf the padding
                switch(Serial.read()) {
                // LED Red "static on" command
                case '1':
                    Serial.read(); // nomf
                    if(Serial.peek() >= '0' && Serial.peek() <= '9') {
                        serialLEDChange = true;
                        serialQueue[SerialQueue_Red] = true;
                        char serialInputS[4] = {0,0,0,0};
                        for(uint n = 0; n < 3; ++n) {
                            serialInputS[n] = Serial.read();
                            if(Serial.peek() < '0' || Serial.peek() > '9')
                                break;
                        }
                        // Static emitting overrides pulse bits
                        serialLEDR = atoi(serialInputS);
                        serialQueue[SerialQueue_LEDPulse] = false;
                        serialLEDPulseColorMap = 0;
                    }
                    break;
                // LED Red "pulse" command (only if not already pulsing)
                case '2':
                    if(!serialQueue[SerialQueue_LEDPulse]) {
                        Serial.read(); // nomf
                        if(Serial.peek() >= '0' && Serial.peek() <= '9') {
                            serialLEDChange = true, serialQueue[SerialQueue_LEDPulse] = true,
                            serialLEDPulseColorMap = 0b00000001; // Set the R LED as the one pulsing only (overwrites the others).
                            char serialInputS[4] = {0,0,0,0};
                            for(uint n = 0; n < 3; ++n) {
                                serialInputS[n] = Serial.read();
                                if(Serial.peek() < '0' || Serial.peek() > '9')
                                    break;
                            }
                            serialLEDPulses = atoi(serialInputS);
                            serialLEDPulsesLast = 0;
                        }
                    }
                    break;
                // LED Red "off" command
                case '0':
                    serialLEDChange = true, serialQueue[SerialQueue_Red] = false, serialQueue[SerialQueue_LEDPulse] = false,
                    serialLEDR = 0, serialLEDPulseColorMap = 0;
                    break;
                }
                break;
              // LED Green bits
              case '3':
                Serial.read(); // nomf the padding
                switch(Serial.read()) {
                // LED Green "static on" command
                case '1':
                    Serial.read(); // nomf
                    if(Serial.peek() >= '0' && Serial.peek() <= '9') {
                        serialLEDChange = true, serialQueue[SerialQueue_Green] = true;
                        char serialInputS[4] = {0,0,0,0};
                        for(uint n = 0; n < 3; ++n) {
                            serialInputS[n] = Serial.read();
                            if(Serial.peek() < '0' || Serial.peek() > '9')
                                break;
                        }
                        serialLEDG = atoi(serialInputS);
                        serialQueue[SerialQueue_LEDPulse] = false, serialLEDPulseColorMap = 0;
                    }
                    break;
                // LED Green "pulse" command (only if not already pulsing)
                case '2':
                    if(!serialQueue[SerialQueue_LEDPulse]) {
                        Serial.read(); // nomf
                        if(Serial.peek() >= '0' && Serial.peek() <= '9') {
                            serialLEDChange = true, serialQueue[SerialQueue_LEDPulse] = true,
                            serialLEDPulseColorMap = 0b00000010; // Set the G LED as the one pulsing only (overwrites the others).
                            char serialInputS[4] = {0,0,0,0};
                            for(uint n = 0; n < 3; ++n) {
                                serialInputS[n] = Serial.read();
                                if(Serial.peek() < '0' || Serial.peek() > '9')
                                    break;
                            }
                            serialLEDPulses = atoi(serialInputS);
                            serialLEDPulsesLast = 0;
                        }
                    }
                    break;
                // LED Green "off" command
                case '0':
                    serialLEDChange = true,
                    serialQueue[SerialQueue_Green] = false, serialQueue[SerialQueue_LEDPulse] = false,
                    serialLEDG = 0, serialLEDPulseColorMap = 0;
                    break;
                }
                break;
              // LED Blue bits
              case '4':
                Serial.read(); // nomf the padding
                switch(Serial.read()) {
                // LED Blue "static on" command
                case '1':
                    Serial.read(); // nomf
                    if(Serial.peek() >= '0' && Serial.peek() <= '9') {
                        serialLEDChange = true, serialQueue[SerialQueue_Blue] = true;
                        char serialInputS[4] = {0,0,0,0};
                        for(uint n = 0; n < 3; ++n) {
                            serialInputS[n] = Serial.read();
                            if(Serial.peek() < '0' || Serial.peek() > '9')
                                break;
                        }
                        serialLEDB = atoi(serialInputS);
                        serialQueue[SerialQueue_LEDPulse] = false;
                        serialLEDPulseColorMap = 0;
                    }
                    break;
                // LED Blue "pulse" command (only if not already pulsing)
                case '2':
                    if(!serialQueue[SerialQueue_LEDPulse]) {
                        Serial.read(); // nomf
                        if(Serial.peek() >= '0' && Serial.peek() <= '9') {
                            serialLEDChange = true, serialQueue[SerialQueue_LEDPulse] = true,
                            serialLEDPulseColorMap = 0b00000100; // Set the B LED as the one pulsing only (overwrites the others).
                            char serialInputS[4] = {0,0,0,0};
                            for(uint n = 0; n < 3; ++n) {
                                serialInputS[n] = Serial.read();
                                if(Serial.peek() < '0' || Serial.peek() > '9')
                                    break;
                            }
                            serialLEDPulses = atoi(serialInputS);
                            serialLEDPulsesLast = 0;
                        }
                    }
                    break;
                // LED Blue "off" command
                case '0':
                    serialLEDChange = true,
                    serialQueue[SerialQueue_Blue] = false, serialQueue[SerialQueue_LEDPulse] = false,
                    serialLEDB = 0, serialLEDPulseColorMap = 0;
                    break;
                }
                break;
              #endif // LED_ENABLE
              #ifdef USES_DISPLAY
              case 'D':
                switch(Serial.read()) {
                case 'A':
                    Serial.read(); // nomf the padding
                    if(Serial.peek() >= '0' && Serial.peek() <= '9') {
                        char serialInputS[4] = {0,0,0,0};
                        for(uint n = 0; n < 3; ++n) {
                            serialInputS[n] = Serial.read();
                            if(Serial.peek() < '0' || Serial.peek() > '9')
                                break;
                        }
                        serialAmmoCount = atoi(serialInputS);
                        serialAmmoCount = constrain(serialAmmoCount, 0, 99);
                        serialDisplayChange = true;
                    }
                    break;
                case 'L':
                    Serial.read(); // nomf the padding
                    if(Serial.peek() >= '0' && Serial.peek() <= '9') {
                        char serialInputS[4] = {0,0,0,0};
                        for(uint n = 0; n < 3; ++n) {
                            serialInputS[n] = Serial.read();
                            if(Serial.peek() < '0' || Serial.peek() > '9')
                                break;
                        }

                        serialLifeCount = atoi(serialInputS);
                        if(FW_Common::OLED.lifeBar) {
                            if(serialLifeCount > FW_Common::dispMaxLife)
                                FW_Common::dispMaxLife = serialLifeCount;
                            FW_Common::dispLifePercentage = (100 * serialLifeCount) / FW_Common::dispMaxLife; // Calculate the Life % to show 
                        }
                        serialDisplayChange = true;
                    }
                    break;
                }
                break;
              #endif // USES_DISPLAY
              #if !defined(USES_SOLENOID) && !defined(USES_RUMBLE) && !defined(LED_ENABLE)
              default:
                //Serial.println("SERIALREAD: Feedback command detected, but no feedback devices are built into this firmware!");
                break;
              #endif
          }
          // End of 'F'
          break;
        // Custom Pulse Overrides
        case 'R':
          switch(Serial.read()) {
              // Solenoid
              case '0':
                Serial.read(); // nomf
                if(Serial.peek() >= '0' && Serial.peek() <= '2') {
                    char serialInput = Serial.read();
                    Serial.read(); // nomf
                    if(Serial.peek() >= '0' && Serial.peek() <='9') {
                        char serialInputS[4] = {0,0,0,0};
                        for(uint n = 0; n < 3; ++n) {
                            serialInputS[n] = Serial.read();
                            if(Serial.peek() < '0' || Serial.peek() > '9')
                                break;
                        }

                        switch(serialInput) {
                        // hold length
                        case '0':
                            serialSolCustomHoldLength = atoi(serialInputS);
                            break;
                        // pause length
                        case '1':
                            serialSolCustomPauseLength = atoi(serialInputS);
                            break;
                        // analog(?)
                        case '2':
                        default:
                            break;
                        }
                    }
                }
                break;
              // Rumble
              case '1':
                Serial.read(); // nomf
                if(Serial.peek() >= '0' && Serial.peek() <= '2') {
                    char serialInput = Serial.read();
                    Serial.read(); // nomf
                    if(Serial.peek() >= '0' && Serial.peek() <= '9') {
                        char serialInputS[4] = {0,0,0,0};
                        for(uint n = 0; n < 3; ++n) {
                            serialInputS[n] = Serial.read();
                            if(Serial.peek() < '0' || Serial.peek() > '9')
                                break;
                        }

                        switch(serialInput) {
                        // hold length
                        case '0':
                            serialRumbCustomHoldLength = atoi(serialInputS);
                            break;
                        // pause length
                        case '1':
                            serialRumbCustomPauseLength = atoi(serialInputS);
                            break;
                        // analog(?)
                        case '2':
                        default:
                            break;
                        }
                    }
                }
                break;
              // LED Red
              case '2':
              // LED Green
              case '3':
              // LED Blue
              case '4':
              default:
                break;
          }
          // End of 'R'
          break;
    }
}

void OF_Serial::SerialHandling()
{
    // The Mamehook feedback system handles most of the timing for us.
    // For the most part, all we have to do is just read and process what it sends us at face value.
    // Solenoid "normal" enable bits need to be monitored to ensure it isn't on for too long, and force shutdown if it is.
    // Solenoid "pulse" bits will borrow from current force feedback settings.
    // Rumble pulse bits are also something we do need to calculate ourselves.
    // The display (if enabled) is handled in the normal Core 0 gunmode run method.

    #ifdef USES_SOLENOID
      if(OF_Prefs::toggles[OF_Const::solenoid]) {
          // Solenoid "on" command
          if(serialQueue[SerialQueue_Solenoid]) {
              if(OF_FFB::GetSolenoid()) {
                  if(millis() - serialSolTimestamp > SERIAL_SOLENOID_MAXSHUTOFF) {
                      OF_FFB::SetSolenoid(LOW);
                      serialQueue[SerialQueue_Solenoid] = false;
                  }
              } else {
                  OF_FFB::SetSolenoid(HIGH);
                  serialSolTimestamp = millis();
              }
          // Solenoid "pulse" command
          } else if(serialQueue[SerialQueue_SolPulse]) {
              if(!serialSolPulsesLast) {                            // Have we started pulsing?
                  OF_FFB::SetSolenoid(HIGH);  // Start pulsing it on!
                  serialSolPulsesLast++;                                 // Start the sequence.
                  serialSolPulsesLastUpdate = millis();                  // timestamp
              } else if(serialSolPulsesLast <= serialSolPulses) {   // Have we met the pulses quota?
                  if(OF_FFB::GetSolenoid()) {
                      // custom hold length
                      if(serialSolCustomHoldLength) {
                          if(millis() - serialSolPulsesLastUpdate >= serialSolCustomHoldLength) {
                              OF_FFB::SetSolenoid(LOW);  // Start pulsing it off.
                              if(serialSolPulsesLast >= serialSolPulses)
                                  serialQueue[SerialQueue_SolPulse] = false;
                              else serialSolPulsesLast++, serialSolPulsesLastUpdate = millis();  // Timestamp our last pulse event.
                          }
                      // current settings hold length
                      } else if(millis() - serialSolPulsesLastUpdate >= OF_Prefs::settings[OF_Const::solenoidOnLength]) {
                          OF_FFB::SetSolenoid(LOW);  // Start pulsing it off.
                          if(serialSolPulsesLast >= serialSolPulses)
                              serialQueue[SerialQueue_SolPulse] = false;
                          else serialSolPulsesLast++, serialSolPulsesLastUpdate = millis();  // Timestamp our last pulse event.
                      }
                  } else {
                      // custom pause length
                      if(serialSolCustomPauseLength) {
                          if(millis() - serialSolPulsesLastUpdate >= serialSolCustomPauseLength) {
                              OF_FFB::SetSolenoid(HIGH); // Start pulsing it on.
                              serialSolPulsesLastUpdate = millis();          // Timestamp our last pulse event.
                          }
                      // current settings pause length
                      } else if(millis() - serialSolPulsesLastUpdate >=
                                OF_Prefs::settings[OF_Const::solenoidOffLength] << OF_FFB::autofireDoubleLengthWait ? 1 : 0) {
                          OF_FFB::SetSolenoid(HIGH); // Start pulsing it on.
                          serialSolPulsesLastUpdate = millis();          // Timestamp our last pulse event.
                      }
                  }
              }
          // Solenoid "off" command
          } else OF_FFB::SetSolenoid(LOW);
      // solenoid toggle not allowed, just force it off.
      } else OF_FFB::SetSolenoid(LOW);
  #endif // USES_SOLENOID

  #ifdef USES_RUMBLE
      if(OF_Prefs::toggles[OF_Const::rumble]) {
          // Rumble "on" command
          if(serialQueue[SerialQueue_Rumble]) {
              analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength]); // turn/keep it on.
          // Rumble "pulse" command
          } else if(serialQueue[SerialQueue_RumbPulse]) {
              // Pulses start
              if(!serialRumbPulsesLast) {
                  if(serialRumbCustomHoldLength && serialRumbCustomPauseLength)
                       analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength]);
                  else analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength] / 3);
                  serialRumbPulseStage = 0;                              // Set that we're at stage 0.
                  serialRumbPulsesLast++;
                  serialRumbPulsesLastUpdate = millis();
              // Pulses processing
              } else if(serialRumbPulsesLast <= serialRumbPulses) {
                  // G4IR-style on/off style ramping
                  if(serialRumbCustomHoldLength && serialRumbCustomPauseLength) {
                      if(!serialRumbPulseStage) {
                          if(millis() - serialRumbPulsesLastUpdate > serialRumbCustomHoldLength) {
                              serialRumbPulseStage = 0;
                              #ifdef ARDUINO_ARCH_ESP32
                                analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0); // [ESP32_PORT] per ESP32
                              #else //rp2040
                              digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
                              #endif
                              if(serialRumbPulsesLast >= serialRumbPulses)
                                  serialQueue[SerialQueue_RumbPulse] = false;
                          }
                      } else if(millis() - serialRumbPulsesLastUpdate > serialRumbCustomPauseLength) {
                          serialRumbPulseStage++;
                          analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength]);
                      }
                  // OF-style analog ramping
                  } else if(millis() - serialRumbPulsesLastUpdate > serialRumbPulsesLength) { // have we waited enough time between pulse stages?
                      switch(serialRumbPulseStage) {
                          // Rising to Sustain
                          case 0:
                              analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength]);
                              serialRumbPulseStage++;                    // Increments the stage of the pulse.
                              serialRumbPulsesLastUpdate = millis();     // and timestamps when we've had updated this last.
                              break;                                     // Then quits until next pulse stage
                          // Sustain to Falling
                          case 1:
                              analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength] / 2);
                              serialRumbPulseStage++;
                              serialRumbPulsesLastUpdate = millis();
                              break;
                          // Falloff
                          case 2:
                              analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength] / 3);
                              serialRumbPulseStage++;
                              serialRumbPulsesLastUpdate = millis();
                              break;
                          // Check
                          case 3:
                              if(serialRumbPulsesLast >= serialRumbPulses) {
                                  #ifdef ARDUINO_ARCH_ESP32
                                    analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0); // [ESP32_PORT] per ESP32
                                  #else //rp2040
                                  digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
                                  #endif
                                  serialQueue[SerialQueue_RumbPulse] = false;
                              } else serialRumbPulsesLast++, serialRumbPulseStage = 0;
                              break;
                      }
                  }
              }
          // Rumble "off"
          } else 
                #ifdef ARDUINO_ARCH_ESP32
                    analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0); // [ESP32_PORT] per ESP32
                #else //rp2040
                    digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
                #endif        
      // Rumble disabled, not allowed to be on
      } else 
                #ifdef ARDUINO_ARCH_ESP32
                    analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0); // [ESP32_PORT] per ESP32
                #else //rp2040
                    digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
                #endif  
  #endif // USES_RUMBLE

  #ifdef LED_ENABLE
    if(serialLEDChange) {                                     // Has the LED command state changed?
        // LED "pulse" command
        if(serialQueue[SerialQueue_LEDPulse]) {
            // LED pulsing start
            if(!serialLEDPulsesLast) {                        // Are we just starting?
                serialQueue[SerialQueue_Red] = false, serialQueue[SerialQueue_Green] = false, serialQueue[SerialQueue_Blue] = false;
                serialLEDPulsesLast = 1;                           // Set that we have started.
                serialLEDPulseRising = true;                       // Set the LED cycle to rising.
                // Reset all the LEDs to zero, the color map will tell us which one to focus on.
                serialLEDR = 0, serialLEDG = 0, serialLEDB = 0;
            // LED pulsing processing
            } else if(serialLEDPulsesLast <= serialLEDPulses) {
                if(millis() - serialLEDPulsesLastUpdate > serialLEDPulsesLength) { // have we waited enough time between pulse stages?
                    switch(serialLEDPulseColorMap) {           // Check the color map
                        case 0b00000001:                       // Basically for R, G, or B,
                            if(serialLEDPulseRising) {
                                serialLEDR += 3;                   // Set the LED value up by three (it's easiest to do blindly like this without over/underflowing tbh)
                                if(serialLEDR == 255)       // If we've reached the max value,
                                    serialLEDPulseRising = false;  // Set that we're in the falling state now.
                            } else {
                                serialLEDR -= 3;                   // Decrement the value.
                                if(serialLEDR == 0) {         // If the LED value has reached the lowest point,
                                    serialLEDPulseRising = true;   // Set that we should be in the rising part of a new cycle.
                                    serialLEDPulsesLast++;         // This was a pulse cycle, so increment that.
                                }
                            }
                            serialLEDPulsesLastUpdate = millis(); // Timestamp this event.
                            break;                             // And get out.
                        case 0b00000010:
                            if(serialLEDPulseRising) {
                                serialLEDG += 3;
                                if(serialLEDG == 255)
                                    serialLEDPulseRising = false;
                            } else {
                                serialLEDG -= 3;
                                if(serialLEDG == 0) {
                                    serialLEDPulseRising = true;
                                    serialLEDPulsesLast++;
                                }
                            }
                            serialLEDPulsesLastUpdate = millis();
                            break;
                        case 0b00000100:
                            if(serialLEDPulseRising) {
                                serialLEDB += 3;
                                if(serialLEDB == 255)
                                    serialLEDPulseRising = false;
                            } else {
                                serialLEDB -= 3;
                                if(serialLEDB == 0) {
                                    serialLEDPulseRising = true;
                                    serialLEDPulsesLast++;
                                }
                            }
                            serialLEDPulsesLastUpdate = millis();
                            break;
                    }
                    // Then, commit the changed value.
                    OF_RGB::LedUpdate(serialLEDR, serialLEDG, serialLEDB);
                }
            // LED pulsing finishing
            } else serialLEDPulseColorMap = 0b00000000, serialQueue[SerialQueue_LEDPulse] = false;
        // Any LED static bits
        } else if(serialQueue[SerialQueue_Red] ||           // Are either the R,
                  serialQueue[SerialQueue_Green] ||         // G,
                  serialQueue[SerialQueue_Blue]) {          // OR B digital bits set to on?
            // Command the LED to change/turn on with the values serialProcessing set for us.
            OF_RGB::LedUpdate(serialLEDR, serialLEDG, serialLEDB);
            serialLEDChange = false;                               // Set the bit to off.
        // LEDs off
        } else OF_RGB::LedOff(), serialLEDChange = false;     // We've done the change, so set it off to reduce redundant LED updates.
    }
    #endif // LED_ENABLE
}
#endif // MAMEHOOKER

// ==================== FINE CODICE PER MAMEHOOKER =======================================

// ===== Desktop App framed serial protocol ===================================

// Pin map that was active before an App commit started. It is needed because
// OF_Prefs::pins is overwritten by incoming records before PinsReset() runs.
static int8_t appSerialPinsBeforeCommit[OF_Const::boardInputsCount] = {};
static bool appSerialPinsBeforeCommitValid = false;

uint8_t OF_Serial::AppSerialCRC8(const uint8_t *data, uint16_t length)
{
    uint8_t crc = 0;

    while(length--) {
        crc ^= *data++;
        for(uint8_t bit = 0; bit < 8; ++bit)
            crc = (crc & 0x80) ? (uint8_t)((crc << 1) ^ 0x9B) : (uint8_t)(crc << 1);
    }

    return crc;
}

uint8_t OF_Serial::AppSerialNextSequence()
{
    ++appSerialTxSequence;
    if(appSerialTxSequence == 0)
        appSerialTxSequence = 1;

    return appSerialTxSequence;
}

void OF_Serial::AppSerialSessionBegin()
{
    appSerialSessionActive = true;
    appSerialRxLength = 0;
    appSerialRxTimestamp = 0;
    appSerialTxSequence = 0;
    appSerialRawDockState = 0;
    appSerialLastRxValid = false;
    appSerialCommitActive = false;
    appSerialCommitFailed = false;
    appSerialCalibrationCancel = false;
    // A new serial session must not discard RAM awaiting a successful save.
    FW_Common::dockedSaving = appSerialPinsBeforeCommitValid;

    appSerialWaitingForAck = false;
    appSerialDispatching = false;
    appSerialProcessingDeferred = false;
    appSerialDeferredValid = false;
}

void OF_Serial::AppSerialSessionEnd()
{
    appSerialSessionActive = false;
    appSerialRxLength = 0;
    appSerialRxTimestamp = 0;
    appSerialRawDockState = 0;
    appSerialLastRxValid = false;
    appSerialCommitActive = false;
    appSerialCommitFailed = false;
    appSerialCalibrationCancel = false;
    FW_Common::dockedSaving = appSerialPinsBeforeCommitValid;

    appSerialWaitingForAck = false;
    appSerialDeferredValid = false;
}

bool OF_Serial::AppSerialWriteFrame(uint8_t typeFlags, uint8_t command, uint8_t sequence, const void *payload, uint8_t length)
{
    if(length > APP_SERIAL_MAX_PAYLOAD || (length > 0 && payload == nullptr))
        return false;

    appSerialTxBuffer[0] = APP_SERIAL_START_1;
    appSerialTxBuffer[1] = APP_SERIAL_START_2;
    appSerialTxBuffer[2] = typeFlags;
    appSerialTxBuffer[3] = command;
    appSerialTxBuffer[4] = sequence;
    appSerialTxBuffer[5] = length;

    if(length > 0)
        memcpy(&appSerialTxBuffer[6], payload, length);

    appSerialTxBuffer[6 + length] = AppSerialCRC8(&appSerialTxBuffer[2], (uint16_t)(4 + length));
    const uint16_t frameLength = (uint16_t)(APP_SERIAL_OVERHEAD + length);

    const bool written = Serial.write(appSerialTxBuffer, frameLength) == frameLength;

    // Reliable frames are followed immediately by a wait for the peer's ACK.
    // Push them out now instead of relying on the USB/serial buffer latency.
    // Best-effort real-time events deliberately remain non-blocking.
    if(written && (typeFlags & APP_SERIAL_TYPE_MASK) != APP_SERIAL_TYPE_EVENT)
        Serial.flush();

    return written;
}

void OF_Serial::AppSerialSendAck(uint8_t command, uint8_t sequence)
{
    AppSerialWriteFrame(APP_SERIAL_TYPE_ACK, command, sequence, nullptr, 0);
}

bool OF_Serial::AppSerialReadFrame(AppSerialFrame_s &frame)
{
    if(appSerialRxLength > 0 && millis() - appSerialRxTimestamp > APP_SERIAL_FRAME_TIMEOUT)
        appSerialRxLength = 0;

    for(;;) {
        while(appSerialRxLength >= 2) {
            uint16_t start = 0;
            while(start + 1 < appSerialRxLength &&
                  (appSerialRxBuffer[start] != APP_SERIAL_START_1 || appSerialRxBuffer[start + 1] != APP_SERIAL_START_2))
                ++start;

            if(start + 1 >= appSerialRxLength) {
                if(appSerialRxBuffer[appSerialRxLength - 1] == APP_SERIAL_START_1) {
                    appSerialRxBuffer[0] = APP_SERIAL_START_1;
                    appSerialRxLength = 1;
                } else appSerialRxLength = 0;
                break;
            }

            if(start > 0) {
                memmove(appSerialRxBuffer, &appSerialRxBuffer[start], appSerialRxLength - start);
                appSerialRxLength -= start;
            }

            if(appSerialRxLength < 6)
                break;

            const uint8_t typeFlags = appSerialRxBuffer[2];
            const uint8_t type = typeFlags & APP_SERIAL_TYPE_MASK;
            const uint8_t sequence = appSerialRxBuffer[4];
            const uint8_t length = appSerialRxBuffer[5];
            const bool invalidFlags = (typeFlags & (uint8_t)~(APP_SERIAL_TYPE_MASK | APP_SERIAL_FLAG_FINAL)) != 0;
            const bool invalidFinal = (typeFlags & APP_SERIAL_FLAG_FINAL) &&
                                      (type != APP_SERIAL_TYPE_RESPONSE || length != 0);
            const bool invalidSequence = (type == APP_SERIAL_TYPE_EVENT) ? (sequence != 0) : (sequence == 0);
            const bool invalidAck = type == APP_SERIAL_TYPE_ACK && length != 0;

            if(invalidFlags || invalidFinal || invalidSequence || invalidAck || length > APP_SERIAL_MAX_PAYLOAD) {
                memmove(appSerialRxBuffer, &appSerialRxBuffer[1], --appSerialRxLength);
                continue;
            }

            const uint16_t frameLength = (uint16_t)(APP_SERIAL_OVERHEAD + length);
            if(appSerialRxLength < frameLength)
                break;

            const uint8_t receivedCRC = appSerialRxBuffer[6 + length];
            const uint8_t calculatedCRC = AppSerialCRC8(&appSerialRxBuffer[2], (uint16_t)(4 + length));
            if(receivedCRC != calculatedCRC) {
                memmove(appSerialRxBuffer, &appSerialRxBuffer[1], --appSerialRxLength);
                continue;
            }

            frame.typeFlags = typeFlags;
            frame.command = appSerialRxBuffer[3];
            frame.sequence = sequence;
            frame.length = length;
            frame.crc = receivedCRC;
            if(length > 0)
                memcpy(frame.payload, &appSerialRxBuffer[6], length);

            appSerialRxLength -= frameLength;
            if(appSerialRxLength > 0)
                memmove(appSerialRxBuffer, &appSerialRxBuffer[frameLength], appSerialRxLength);

            return true;
        }

        if(!Serial.available())
            return false;

        const int incoming = Serial.read();
        if(incoming < 0)
            return false;

        if(appSerialRxLength >= APP_SERIAL_MAX_FRAME) {
            memmove(appSerialRxBuffer, &appSerialRxBuffer[1], APP_SERIAL_MAX_FRAME - 1);
            appSerialRxLength = APP_SERIAL_MAX_FRAME - 1;
        }

        appSerialRxBuffer[appSerialRxLength++] = (uint8_t)incoming;
        appSerialRxTimestamp = millis();
    }
}

bool OF_Serial::AppSerialRequestMatchesLast(const AppSerialFrame_s &frame)
{
    return appSerialLastRxValid &&
           appSerialLastRxCommand == frame.command &&
           appSerialLastRxSequence == frame.sequence &&
           appSerialLastRxLength == frame.length &&
           appSerialLastRxCRC == frame.crc;
}

void OF_Serial::AppSerialProcessDeferredRequest()
{
    if(appSerialWaitingForAck ||
       appSerialDispatching ||
       appSerialProcessingDeferred)
        return;

    appSerialProcessingDeferred = true;

    while(appSerialSessionActive &&
          appSerialDeferredValid &&
          !appSerialWaitingForAck) {
        const AppSerialFrame_s frame = appSerialDeferredFrame;
        appSerialDeferredValid = false;
        AppSerialHandleFrame(frame);
    }

    appSerialProcessingDeferred = false;
}

void OF_Serial::AppSerialHandleFrame(const AppSerialFrame_s &frame)
{
    const uint8_t type = frame.typeFlags & APP_SERIAL_TYPE_MASK;

    if(type != APP_SERIAL_TYPE_REQUEST)
        return;

    const bool duplicate = AppSerialRequestMatchesLast(frame);

    if(!duplicate) {
        appSerialLastRxValid = true;
        appSerialLastRxCommand = frame.command;
        appSerialLastRxSequence = frame.sequence;
        appSerialLastRxLength = frame.length;
        appSerialLastRxCRC = frame.crc;
    }

    // The complete request frame has been received and validated.
    AppSerialSendAck(frame.command, frame.sequence);

    if(!duplicate) {
        const bool alreadyDispatching = appSerialDispatching;
        appSerialDispatching = true;
        AppSerialDispatchRequest(frame);
        appSerialDispatching = alreadyDispatching;
    }

    AppSerialProcessDeferredRequest();
}

bool OF_Serial::AppSerialSendReliable(uint8_t typeFlags,
                                      uint8_t command,
                                      const void *payload,
                                      uint8_t length,
                                      uint8_t sequence)
{
    if(!appSerialSessionActive || length > APP_SERIAL_MAX_PAYLOAD)
        return false;

    // Single-result commit replies echo the request sequence. Other response
    // streams retain their independent per-frame sequence, as before.
    if(sequence == 0)
        sequence = AppSerialNextSequence();

    bool success = false;

    appSerialWaitingForAck = true;

    for(uint8_t attempt = 0;
        attempt <= APP_SERIAL_MAX_RETRIES &&
        appSerialSessionActive &&
        !success;
        ++attempt) {
        if(!AppSerialWriteFrame(typeFlags,
                                command,
                                sequence,
                                payload,
                                length))
            break;

        const unsigned long started = millis();

        while(appSerialSessionActive &&
              millis() - started < APP_SERIAL_ACK_TIMEOUT &&
              !success) {
            AppSerialFrame_s frame;

            while(AppSerialReadFrame(frame)) {
                const uint8_t type =
                    frame.typeFlags & APP_SERIAL_TYPE_MASK;

                if(type == APP_SERIAL_TYPE_ACK) {
                    if(frame.command == command &&
                       frame.sequence == sequence) {
                        success = true;
                        break;
                    }

                    continue;
                }

                if(type == APP_SERIAL_TYPE_REQUEST) {
                    if(AppSerialRequestMatchesLast(frame)) {
                        // The App is retrying the request whose response is
                        // currently waiting for confirmation.
                        AppSerialSendAck(frame.command, frame.sequence);
                    } else {
                        // A following valid request proves that the stop-and-
                        // wait peer has advanced beyond this response. Retain
                        // the request, complete this send and dispatch it after
                        // unwinding the current command.
                        if(!appSerialDeferredValid) {
                            appSerialDeferredFrame = frame;
                            appSerialDeferredValid = true;
                        }

                        success = true;
                        break;
                    }
                }
            }

            if(!success)
                yield();
        }
    }

    appSerialWaitingForAck = false;

    if(!appSerialDispatching)
        AppSerialProcessDeferredRequest();

    return success;
}

bool OF_Serial::AppSerialSendResponse(uint8_t command,
                                      const void *payload,
                                      uint8_t length,
                                      bool final)
{
    if(final && length != 0)
        return false;

    const uint8_t typeFlags =
        APP_SERIAL_TYPE_RESPONSE |
        (final ? APP_SERIAL_FLAG_FINAL : 0);

    return AppSerialSendReliable(typeFlags, command, payload, length);
}

bool OF_Serial::AppSerialSendEvent(uint8_t command, const void *payload, uint8_t length)
{
    if(!appSerialSessionActive)
        return false;

    return AppSerialWriteFrame(APP_SERIAL_TYPE_EVENT, command, 0, payload, length);
}

void OF_Serial::AppSerialSendError(uint8_t error)
{
    if(error)
        AppSerialSendResponse(OF_Const::sError, &error, 1);
    else AppSerialSendResponse(OF_Const::sError);
}

void OF_Serial::AppSerialSendCommitError(const AppSerialFrame_s &frame)
{
    const uint8_t payload[] = {APP_SERIAL_ERR_COMMIT_RETRY, frame.command};
    AppSerialSendReliable(APP_SERIAL_TYPE_RESPONSE,
                          OF_Const::sError,
                          payload,
                          sizeof(payload),
                          frame.sequence);
}

bool OF_Serial::AppSerialTakeCalibrationCancel()
{
    const bool requested = appSerialCalibrationCancel;
    appSerialCalibrationCancel = false;
    return requested;
}

// Serial Buffer in Docked Mode should always be read by the main core on multicore systems.
void OF_Serial::SerialProcessingDocked()
{

    // ExecCalMode() is still inside the dispatch of the calibration command.
    // A request received while a calibration response is waiting for its ACK
    // is deferred, so process it explicitly at the next calibration poll.
    if(appSerialSessionActive &&
       appSerialDeferredValid &&
       !appSerialWaitingForAck &&
       (FW_Common::gunMode == FW_Const::GunMode_Calibration ||
        FW_Common::gunMode == FW_Const::GunMode_Verification)) {
        const AppSerialFrame_s deferredFrame = appSerialDeferredFrame;

        appSerialDeferredValid = false;
        AppSerialHandleFrame(deferredFrame);

        if(!appSerialSessionActive)
            return;
    }

    // The two-byte docking handshake intentionally remains unframed. Once it
    // succeeds, every App byte in both directions uses the framed protocol.
    if(!appSerialSessionActive) {
        if(appSerialRawDockState != 0 &&
           millis() - appSerialRxTimestamp > APP_SERIAL_FRAME_TIMEOUT)
            appSerialRawDockState = 0;

        while(Serial.available()) {
            const int incoming = Serial.read();
            if(incoming < 0)
                return;

            if(appSerialRawDockState != 0 &&
               incoming == OF_Const::sDock2) {
                AppSerialSessionBegin();
                FW_Common::SetMode(FW_Const::GunMode_Docked);
                break;
            }

            appSerialRawDockState =
                incoming == OF_Const::sDock1 ? 1 : 0;

            if(appSerialRawDockState != 0)
                appSerialRxTimestamp = millis();
        }
    }

    if(appSerialSessionActive) {
        AppSerialFrame_s frame;
        while(appSerialSessionActive && AppSerialReadFrame(frame))
            AppSerialHandleFrame(frame);
    }
}

bool OF_Serial::AppSerialDispatchCommit(const AppSerialFrame_s &frame,
                                        bool fullDisconnect)
{
    // A new commit also restarts an interrupted attempt.
    // Preserve the calibration data and the original hardware pin map
    // until the complete configuration has been saved successfully.
    if(frame.command == OF_Const::sCommitStart) {
        if(!appSerialPinsBeforeCommitValid) {
            memcpy(appSerialPinsBeforeCommit,
                   OF_Prefs::pins,
                   sizeof(appSerialPinsBeforeCommit));

            appSerialPinsBeforeCommitValid = true;
            FW_Common::buttons.Unset();
        }

        // Start the new transfer from the board defaults.
        // The App subsequently overwrites the custom pins, when enabled,
        // and always transmits the complete button mapping.
        OF_Prefs::LoadPresets();

        FW_Common::dockedSaving = true;
        appSerialCommitActive = true;
        appSerialCommitFailed = false;

        if(!AppSerialSendReliable(APP_SERIAL_TYPE_RESPONSE,
                                  frame.command,
                                  nullptr,
                                  0,
                                  frame.sequence)) {
            // Keep the received data and original hardware pin map so that
            // the App can start another complete transfer.
            appSerialCommitActive = false;
        }

        return true;
    }

    // Explicit reset commands remain available even if a transfer stopped.
    if(!appSerialCommitActive ||
       frame.command == OF_Const::sClearFlash ||
       frame.command == OF_Const::sRebootToBootloader)
        return false;

    switch(frame.command) {
        case OF_Const::sSave: {
            const bool complete = !appSerialCommitFailed;

            const bool saved =
                complete &&
                FW_Common::SavePreferences() == OF_Prefs::Error_Success;

            if(saved) {
                // Deinitialize the hardware using the pin map that was active
                // before the configuration transfer.
                FW_Common::PinsReset(appSerialPinsBeforeCommit);

                // Apply all newly saved runtime settings.
                FW_Common::CameraSet();
                FW_Common::SetMode(FW_Common::gunMode); // Riallineamento Square Advanced dopo il salvataggio. Con USE_SQUARE_ADVANCED disabilitato è praticamente ininfluente.
                FW_Common::FeedbackSet();
                FW_Common::UpdateBindings(true);

                #ifdef LED_ENABLE
                if(FW_Common::gunMode == FW_Const::GunMode_Docked) {
                    OF_RGB::LedUpdate(127, 127, 255);
                } else if(FW_Common::gunMode == FW_Const::GunMode_Pause) {
                    OF_RGB::SetLedPackedColor(
                        OF_Prefs::profiles[
                            OF_Prefs::currentProfile
                        ].color);
                }
                #endif // LED_ENABLE

                FW_Common::buttons.Begin();

                appSerialPinsBeforeCommitValid = false;
                FW_Common::dockedSaving = false;
            }

            // On failure, retain the configuration in RAM, the original
            // hardware pin map and the recovery guard for another attempt.
            appSerialCommitActive = false;

            static const char successText[] =
                " (Successfully saved to LittleFS Storage)";
            static const char failureText[] =
                " (Failed to save to LittleFS Storage)";

            const char *text = saved ? successText : failureText;
            const uint8_t textLength = (uint8_t)strlen(text);

            uint8_t payload[1 + sizeof(successText)];
            payload[0] = saved;
            memcpy(&payload[1], text, textLength);

            AppSerialSendReliable(APP_SERIAL_TYPE_RESPONSE,
                                  OF_Const::sSave,
                                  payload,
                                  (uint8_t)(1 + textLength),
                                  frame.sequence);

            #ifdef USES_DISPLAY
            // SavePreferences() skips its visual effects during an App commit.
            // Send the result before holding the OLED message so that the
            // display delay does not extend the App response time.
            if(FW_Common::OLED.display != nullptr) {
                FW_Common::OLED.ScreenModeChange(
                    saved ? ExtDisplay::Screen_SaveSuccess :
                            ExtDisplay::Screen_SaveError);

                // Same duration as the original save feedback patterns.
                delay(saved ? 285 : 410);
                FW_Common::RedrawDisplay();
            }
            #endif // USES_DISPLAY

            break;
        }

        case OF_Const::serialTerminator: {
            // End only this attempt: never reload flash over unsaved calibration.
            appSerialCommitActive = false;
            appSerialCommitFailed = false;

            AppSerialSendResponse(OF_Const::serialTerminator, nullptr, 0, true);
            if(fullDisconnect)
                AppSerialSessionEnd(); // Stay Docked with pending RAM.
            break;
        }

        case OF_Const::sCommitID:
            if(frame.length == sizeof(OF_Prefs::USBMap_t))
                memcpy(&OF_Prefs::usb, frame.payload, sizeof(OF_Prefs::USBMap_t));
            else appSerialCommitFailed = true;
            break;

        case OF_Const::sCommitToggles:
            if(!AppSerialReceiveRecord(frame.payload, frame.length,
                                       OF_Prefs::toggles,
                                       OF_Prefs::OFPresets.boolTypes_Strings,
                                       sizeof(OF_Prefs::toggles) / OF_Const::boolTypesCount))
                appSerialCommitFailed = true;
            break;

        case OF_Const::sCommitPins:
            if(!AppSerialReceiveRecord(frame.payload, frame.length,
                                       OF_Prefs::pins,
                                       OF_Prefs::OFPresets.boardInputs_Strings,
                                       sizeof(OF_Prefs::pins) / OF_Const::boardInputsCount))
                appSerialCommitFailed = true;
            break;

        case OF_Const::sCommitSettings:
            if(!AppSerialReceiveRecord(frame.payload, frame.length,
                                       OF_Prefs::settings,
                                       OF_Prefs::OFPresets.settingsTypes_Strings,
                                       sizeof(OF_Prefs::settings) / OF_Const::settingsTypesCount))
                appSerialCommitFailed = true;
            break;

        case OF_Const::sCommitBtns:
            if(!AppSerialReceiveRecord(frame.payload, frame.length,
                                       OF_Prefs::backupButtonDesc,
                                       OF_Prefs::OFPresets.boardInputs_Strings,
                                       sizeof(OF_Prefs::backupButtonDesc) / ButtonCount))
                appSerialCommitFailed = true;
            break;

        case OF_Const::sCommitProfile:
            if(!AppSerialReceiveRecord(frame.payload, frame.length,
                                       OF_Prefs::profiles,
                                       OF_Prefs::OFPresets.profSettingTypes_Strings,
                                       sizeof(uint32_t),
                                       true))
                appSerialCommitFailed = true;
            break;

        default:
            appSerialCommitFailed = true;
            appSerialCommitActive = false;
            AppSerialSendCommitError(frame);
            break;
    }

    return true;
}

void OF_Serial::AppSerialRestoreRunState()
{
    if(!FW_Common::justBooted)
        FW_Common::SetMode(FW_Const::GunMode_Run);
    else
        FW_Common::SetMode(FW_Const::GunMode_Init);

    FW_Common::SetRunMode(
        (FW_Const::RunMode_e)
        OF_Prefs::profiles[OF_Prefs::currentProfile].runMode);
}

void OF_Serial::AppSerialDispatchRequest(const AppSerialFrame_s &frame)
{
    const bool fullDisconnect =
        frame.command == OF_Const::serialTerminator &&
        frame.length == 2 &&
        frame.payload[0] == OF_Const::serialTerminator &&
        frame.payload[1] == OF_Const::serialTerminator;

    // Keep the profile stable until ExecCalMode() has restored or accepted it.
    if((FW_Common::gunMode == FW_Const::GunMode_Calibration ||
        FW_Common::gunMode == FW_Const::GunMode_Verification) &&
       frame.command != OF_Const::serialTerminator) {
        appSerialCalibrationCancel = true;
        AppSerialSendError();
        return;
    }

    if(AppSerialDispatchCommit(frame, fullDisconnect))
        return;

    // Do not use hardware with a partly updated configuration.
    if(appSerialPinsBeforeCommitValid) {
        switch(frame.command) {
        case OF_Const::sGetToggles:
        case OF_Const::sGetPins:
        case OF_Const::sGetSettings:
        case OF_Const::sGetBtns:
        case OF_Const::sGetProfile:
        case OF_Const::serialTerminator:
        case OF_Const::sClearFlash:
        case OF_Const::sRebootToBootloader:
            break;

        default:
            AppSerialSendCommitError(frame);
            return;
        }
    }

    switch(frame.command) {  

    case OF_Const::serialTerminator:
        if(FW_Common::gunMode == FW_Const::GunMode_Calibration ||
           FW_Common::gunMode == FW_Const::GunMode_Verification) {
            appSerialCalibrationCancel = true;

            // ExecCalMode() must consume the cancellation flag before
            // AppSerialSessionEnd() clears it.
            if(fullDisconnect) {
                AppSerialSendResponse(OF_Const::serialTerminator,
                                      nullptr,
                                      0,
                                      true);

                // Do not call AppSerialSessionEnd() here: it would clear
                // appSerialCalibrationCancel before ExecCalMode() consumes it.
                appSerialSessionActive = false;
            }

            break;
        }

        // A cancel can cross a calibration End already confirmed by trigger.
        // Only the explicit FE FE FE request disconnects an idle session.
        if(!fullDisconnect)
            break;

        AppSerialSendResponse(OF_Const::serialTerminator,
                              nullptr,
                              0,
                              true);

        AppSerialSessionEnd();

        if(appSerialPinsBeforeCommitValid)
            break; // Do not run with pending, partly updated configuration.

        AppSerialRestoreRunState();

        break;

    case OF_Const::sGetToggles:
        AppSerialSendRecords(frame.command, OF_Prefs::toggles,
                             OF_Prefs::OFPresets.boolTypes_Strings,
                             sizeof(OF_Prefs::toggles) / OF_Const::boolTypesCount);
        break;

    case OF_Const::sGetPins:
        AppSerialSendRecords(frame.command, OF_Prefs::pins,
                             OF_Prefs::OFPresets.boardInputs_Strings,
                             sizeof(OF_Prefs::pins) / OF_Const::boardInputsCount);
        break;

    case OF_Const::sGetSettings:
        AppSerialSendRecords(frame.command, OF_Prefs::settings,
                             OF_Prefs::OFPresets.settingsTypes_Strings,
                             sizeof(OF_Prefs::settings) / OF_Const::settingsTypesCount);
        break;

    case OF_Const::sGetBtns:
        AppSerialSendRecords(frame.command, OF_Prefs::backupButtonDesc,
                             OF_Prefs::OFPresets.boardInputs_Strings,
                             sizeof(OF_Prefs::backupButtonDesc) / ButtonCount);
        break;

    case OF_Const::sGetProfile:
    {
        bool sent = true;
        for(int prof = 0; prof < PROFILE_COUNT && sent; ++prof)
            sent = AppSerialSendRecords(frame.command,
                                        &OF_Prefs::profiles[prof],
                                        OF_Prefs::OFPresets.profSettingTypes_Strings,
                                        sizeof(uint32_t),
                                        prof,
                                        false);
        if(sent && AppSerialSendResponse(frame.command, nullptr, 0, true) &&
           appSerialPinsBeforeCommitValid)
            AppSerialSendCommitError(frame); // A reconnected App must save again.
        break;
    }

    case OF_Const::sIRTest:
        if(FW_Common::camNotAvailable) {
            AppSerialSendError(OF_Const::sErrCam);
        } else if(frame.payload[0]) {
            FW_Common::SetRunMode(FW_Const::RunMode_Processing);
        } else if(FW_Common::runMode == FW_Const::RunMode_Processing) {
            switch(OF_Prefs::profiles[OF_Prefs::currentProfile].runMode) {
            case FW_Const::RunMode_Normal:
                FW_Common::SetRunMode(FW_Const::RunMode_Normal);
                break;

            case FW_Const::RunMode_Average:
                FW_Common::SetRunMode(FW_Const::RunMode_Average);
                break;

            case FW_Const::RunMode_Average2:
                FW_Common::SetRunMode(FW_Const::RunMode_Average2);
                break;
            }
        }
        break;

    case OF_Const::sCaliProfile: {
        const uint8_t operation = frame.payload[0];
        const uint8_t profile = frame.payload[1];

        FW_Common::SelectCalProfile(profile);

        const uint8_t currentProfile = OF_Prefs::currentProfile;

        if(!AppSerialSendResponse(
                OF_Const::sCurrentProf,
                &currentProfile,
                sizeof(currentProfile))) {
            AppSerialSendError();
            break;
        }

        switch(operation) {
        case OF_Const::sCaliProfile:
            // When used inside the payload, sCaliProfile means:
            // select the requested profile without starting calibration.
            break;

        case OF_Const::sCaliStart:
            if(FW_Common::camNotAvailable) {
                AppSerialSendError(OF_Const::sErrCam);
            } else {
                const uint8_t caliSettings = frame.payload[2];

                FW_Common::SetIrSensitivity(caliSettings & 0x0F);
                FW_Common::SetIrLayout(caliSettings >> 4);
                FW_Common::SetMode(FW_Const::GunMode_Calibration);
                FW_Common::ExecCalMode(true);

                // A full App disconnect received inside ExecCalMode()
                // first cancels calibration so that the backed-up profile
                // is restored, then completes the normal undock here.
                if(!appSerialSessionActive) {
                    AppSerialSessionEnd();
                    AppSerialRestoreRunState();
                }
            }
            break;
        }
        break;
    }

    #ifdef USES_SOLENOID
    case OF_Const::sTestSolenoid:
        OF_FFB::SetSolenoid(HIGH);
        delay(OF_Prefs::settings[OF_Const::solenoidOnLength]);
        OF_FFB::SetSolenoid(LOW);
        break;
    #endif

    #ifdef USES_RUMBLE
    case OF_Const::sTestRumble:
        analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength]);
        delay(OF_Prefs::settings[OF_Const::rumbleInterval]);
        #ifdef ARDUINO_ARCH_ESP32
        analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0);
        #else
        digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);
        #endif
        break;
    #endif

    #ifdef LED_ENABLE
    case OF_Const::sTestLEDR:
        OF_RGB::LedUpdate(255, 0, 0);
        break;
    case OF_Const::sTestLEDG:
        OF_RGB::LedUpdate(0, 255, 0);
        break;
    case OF_Const::sTestLEDB:
        OF_RGB::LedUpdate(0, 0, 255);
        break;
    #endif

    case OF_Const::sClearFlash:
        OF_Prefs::ResetPreferences();
        Serial.flush();
        #ifdef ARDUINO_ARCH_ESP32
            #ifdef OPENFIRE_WIRELESS_ENABLE
            if(TinyUSBDevices.onBattery) {
                // The dongle reset notification can be added separately.
            }
            #endif
            //esp_restart();
            esp_rom_software_reset_system();
        #else
            rp2040.reboot();
        #endif
        break;

    case OF_Const::sRebootToBootloader:
        Serial.flush();
        FW_Common::RebootToBootloader();
        break;

    default:
        break;
    }
}

bool OF_Serial::AppSerialSendRecords(uint8_t command,
                                     const void *data,
                                     const std::unordered_map<std::string_view, int> &fields,
                                     size_t dataSize,
                                     int profile,
                                     bool sendFinal)
{
    uint8_t payload[APP_SERIAL_MAX_PAYLOAD];

    const bool profileData = profile >= 0;
    const bool buttonData = data == OF_Prefs::backupButtonDesc;

    for(const auto &pair : fields) {
        if(pair.second < 0)
            continue;

        // CurrentProf belongs to the complete profile set and must therefore
        // be transmitted only once, together with the first profile.
        if(profileData &&
           pair.second == OF_Const::profCurrent &&
           profile != 0)
            continue;

        const size_t nameLength = pair.first.length();
        size_t pos = nameLength + 1;

        // Space for the field name, its terminator and the record metadata.
        if(pos + 2 > sizeof(payload)) {
            AppSerialSendError();
            return false;
        }

        memcpy(payload, pair.first.data(), nameLength);
        payload[nameLength] = '\0';

        if(profileData) {
            if(pair.second == OF_Const::profCurrent) {
                payload[pos++] = (uint8_t)sizeof(uint8_t);
                payload[pos++] = (uint8_t)OF_Prefs::currentProfile;
            } else {
                const size_t valueSize =
                    pair.second == OF_Const::profName ?
                    sizeof(OF_Prefs::ProfileData_s::name) :
                    dataSize;

                if(pos + 2 + valueSize > sizeof(payload)) {
                    AppSerialSendError();
                    return false;
                }

                payload[pos++] = (uint8_t)valueSize;
                payload[pos++] = (uint8_t)profile;

                memcpy(&payload[pos],
                       (const uint8_t*)data + (dataSize * pair.second),
                       valueSize);

                pos += valueSize;
            }
        } else {
            // The last button entry is intentionally not managed by the App,
            // preserving the behaviour of the original implementation.
            if(buttonData &&
               pair.second >= ButtonCount - 1)
                continue;

            if(pos + 1 + dataSize > sizeof(payload)) {
                AppSerialSendError();
                return false;
            }

            payload[pos++] = (uint8_t)dataSize;

            memcpy(&payload[pos],
                   (const uint8_t*)data + (dataSize * pair.second),
                   dataSize);

            pos += dataSize;
        }

        if(!AppSerialSendResponse(command,
                                  payload,
                                  (uint8_t)pos))
            return false;
    }

    if(!sendFinal)
        return true;

    return AppSerialSendResponse(command, nullptr, 0, true);
}

bool OF_Serial::AppSerialReceiveRecord(
    const uint8_t *payload,
    uint8_t length,
    void *data,
    const std::unordered_map<std::string_view, int> &fields,
    size_t dataSize,
    bool profileData)
{
    if(payload == nullptr || length < 2)
        return false;

    const bool buttonData = data == OF_Prefs::backupButtonDesc;

    uint8_t nameLength = 0;

    while(nameLength < length &&
          payload[nameLength] != '\0')
        ++nameLength;

    // At least the string terminator and value-size byte must be present.
    if(nameLength >= length ||
       nameLength + 1 >= length)
        return false;

    const std::string_view fieldName(
        (const char*)payload,
        nameLength);

    size_t pos = nameLength + 1;
    const uint8_t valueSize = payload[pos++];

    // One lookup is sufficient and avoids scanning the field name again.
    const auto field = fields.find(fieldName);

    // Every frame contains exactly one record, therefore an unknown field
    // can be ignored safely for forward compatibility.
    if(field == fields.end())
        return true;

    const int index = field->second;

    if(index < 0)
        return true;

    if(profileData) {
        if(index == OF_Const::profCurrent) {
            // CurrentProf is transported as one byte even though the runtime
            // variable uses the native unsigned-integer type.
            if(valueSize != sizeof(uint8_t) ||
               pos + valueSize != length)
                return false;

            const uint8_t currentProfile = payload[pos];

            if(currentProfile >= PROFILE_COUNT)
                return false;

            OF_Prefs::currentProfile = currentProfile;
            return true;
        }

        // Every ordinary profile record also contains the profile number.
        if(pos >= length)
            return false;

        const uint8_t profNum = payload[pos++];

        if(profNum >= PROFILE_COUNT)
            return false;

        const size_t expectedSize =
            index == OF_Const::profName ?
            sizeof(OF_Prefs::ProfileData_s::name) :
            dataSize;

        if(valueSize != expectedSize ||
           pos + valueSize != length)
            return false;

        memcpy((uint8_t*)&OF_Prefs::profiles[profNum] +
                   (dataSize * index),
               &payload[pos],
               valueSize);

        return true;
    }

    if(valueSize != dataSize ||
       pos + valueSize != length)
        return false;

    // The last button entry is intentionally not managed by the App.
    if(buttonData &&
       index >= ButtonCount - 1)
        return true;

    memcpy((uint8_t*)data + (dataSize * index),
           &payload[pos],
           valueSize);

    return true;
}

void OF_Serial::PrintResults()
{
    if(millis() - lastPrintMillis < 100)
        return;

    #ifdef OPENFIRE_WIRELESS_ENABLE        
        if (!(TinyUSBDevices.onBattery ? SerialWireless : TinyUSBDevice.mounted())) { // [ESP32_PORT] poi decide come sistemare per bene ma così dovrebbe andare bene
    #else
    if(!Serial) {
    #endif
        FW_Common::stateFlags |= FW_Const::StateFlagsDtrReset;
        return;
    }

    if(FW_Common::stateFlags & FW_Const::StateFlag_PrintPreferences) {
        FW_Common::stateFlags &= ~FW_Const::StateFlag_PrintPreferences;

        // Prints basic storage device information
        // an estimation of storage used, though doesn't take extended prefs into account.
        if(FW_Common::stateFlags & FW_Const::StateFlag_PrintPreferencesStorage) {
            FW_Common::stateFlags &= ~FW_Const::StateFlag_PrintPreferencesStorage;
            
            #ifdef SAMCO_FLASH_ENABLE
                unsigned int required = OF_Prefs::Size();

            #ifndef PRINT_VERBOSE
                if(required < flash.size())
                    return;
            #endif

                Serial.print("NV Storage capacity: ");
                Serial.print(flash.size());
                Serial.print(", required size: ");
                Serial.println(required);

            #ifdef PRINT_VERBOSE
                Serial.print("Profile struct size: ");
                Serial.print((unsigned int)sizeof(OF_Prefs::profileData_t));
                Serial.print(", Profile data array size: ");
                Serial.println((unsigned int)sizeof(OF_Prefs::profiles));
            #endif

            #endif // SAMCO_FLASH_ENABLE
        }

        // prints all stored preferences information in a table
        Serial.print("Default Profile: ");
        Serial.println(OF_Prefs::profiles[OF_Prefs::currentProfile].name);
        
        Serial.println("Profiles:");
        for(uint i = 0; i < PROFILE_COUNT; ++i) {
            // report if a profile has been cal'd
            if(OF_Prefs::profiles[i].topOffset && OF_Prefs::profiles[i].bottomOffset &&
              OF_Prefs::profiles[i].leftOffset && OF_Prefs::profiles[i].rightOffset) {
                size_t len = strlen(OF_Prefs::profiles[i].name);
                Serial.print(OF_Prefs::profiles[i].name);
                while(len < 18) {
                    Serial.print(' ');
                    ++len;
                }
                Serial.print("Top: ");
                Serial.print(OF_Prefs::profiles[i].topOffset);
                Serial.print(", Bottom: ");
                Serial.print(OF_Prefs::profiles[i].bottomOffset);
                Serial.print(", Left: ");
                Serial.print(OF_Prefs::profiles[i].leftOffset);
                Serial.print(", Right: ");
                Serial.print(OF_Prefs::profiles[i].rightOffset);
                Serial.print(", TLled: ");
                Serial.print(OF_Prefs::profiles[i].TLled);
                Serial.print(", TRled: ");
                Serial.print(OF_Prefs::profiles[i].TRled);
                //Serial.print(", AdjX: ");
                //Serial.print(OF_Prefs::profiles[i].adjX);
                //Serial.print(", AdjY: ");
                //Serial.print(OF_Prefs::profiles[i].adjY);
                Serial.print(" IR: ");
                Serial.print((unsigned int)OF_Prefs::profiles[i].irSens);
                Serial.print(" Mode: ");
                Serial.print((unsigned int)OF_Prefs::profiles[i].runMode);
                Serial.print(" Layout: ");
                if(OF_Prefs::profiles[i].irLayout)
                    Serial.println("Diamond");
                else Serial.println("Square");
            }
        }
        /*
        Serial.print(finalX);
        Serial.print(" (");
        Serial.print(MoveXAxis);
        Serial.print("), ");
        Serial.print(finalY);
        Serial.print(" (");
        Serial.print(MoveYAxis);
        Serial.print("), H ");
        Serial.println(mySamco.H());*/

        //Serial.print("conMove ");
        //Serial.print(conMoveXAxis);
        //Serial.println(conMoveYAxis);
        
        if(FW_Common::stateFlags & FW_Const::StateFlag_PrintSelectedProfile) {
            FW_Common::stateFlags &= ~FW_Const::StateFlag_PrintSelectedProfile;

            // Print selected profile
            Serial.print("Profile: ");
            Serial.println(OF_Prefs::profiles[OF_Prefs::currentProfile].name);

            // Print current sensitivity
            Serial.print("IR Camera Sensitivity: ");
            Serial.println((int)OF_Prefs::profiles[OF_Prefs::currentProfile].irSens);

            // Subroutine that prints current runmode
            if(FW_Common::runMode < FW_Const::RunMode_Count) {
                Serial.print("Mode: ");
                Serial.println(FW_Const::RunModeLabels[FW_Common::runMode]);
            }

            #ifdef USES_RUMBLE
                Serial.print("Rumble enabled: ");
                if(OF_Prefs::toggles[OF_Const::rumble])
                    Serial.println("True");
                else Serial.println("False");
            #endif // USES_RUMBLE

            #ifdef USES_SOLENOID
                Serial.print("Solenoid enabled: ");
                if(OF_Prefs::toggles[OF_Const::solenoid]) {
                    Serial.println("True");
                    Serial.print("Rapid fire enabled: ");
                    if(OF_Prefs::toggles[OF_Const::autofire])
                        Serial.println("True");
                    else Serial.println("False");

                    Serial.print("Burst fire enabled: ");
                    if(OF_FFB::burstFireActive)
                        Serial.println("True");
                    else Serial.println("False");
                }
                else Serial.println("False");
            #endif // USES_SOLENOID

            // #ifdef ARDUINO_ARCH_RP2040 // [ESP32_PORT] per ESP32
            #ifdef DUAL_CORE
                Serial.println("Running on dual cores.");
            #else
                Serial.println("Running on one core.");
            #endif // DUAL_CORE
            // #endif // ARDUINO_ARCH_RP2040 // [ESP32_PORT] per ESP32

            Serial.printf("Firmware version: v%.1f"
                          #ifdef GIT_HASH
                          "-%s"
                          #endif // GIT_ HASH);
                          , OPENFIRE_VERSION
                          #ifdef GIT_HASH
                          , GIT_HASH
                          #endif // GIT_HASH
                          );
        }
                    
        lastPrintMillis = millis();
    }
}

#ifdef DEBUG_SERIAL
void OF_Serial::PrintDebugSerial()
{
    // only print every second
    if(millis() - serialDbMs >= 1000 && Serial) {
        Serial.print("mode ");
        Serial.print(FW_Common::gunMode);
        Serial.print(", IR pos fps ");
        Serial.print(irPosCount);
        Serial.print(", loop/sec ");
        Serial.print(frameCount);

        /*
        Serial.print(", Mouse X,Y ");
        Serial.print(FW_Common::conMoveXAxis);
        Serial.print(",");
        Serial.println(FW_Common::conMoveYAxis);
        */
        
        frameCount = 0;
        irPosCount = 0;
        serialDbMs = millis();
    }
}
#endif // DEBUG_SERIAL

bool OF_Serial::Serial_available(uint8_t min) 
{
    // in futuro valutare di togliere questa funzione
    if ((Serial.available() >= min)) return true;
    else {
        unsigned long timer_out = millis();
        while ((Serial.available() < min) && (millis() - timer_out < 1000)) yield();
        return Serial.available() >= min ? true : false;
    }
}

// ============ [ESP32_PORT] ============
// restore Serial after it was redefined for serial connections / ripristino di Serial dopo definizione per connessione seriali ==============
#ifdef OPENFIRE_WIRELESS_ENABLE
    #undef Serial
    #ifdef AUX_SERIAL
        #define Serial AUX_SERIAL
        #undef AuxSerial
    #endif
#endif // OPENFIRE_WIRELESS_ENABLE
// ============ [ESP32_PORT] ============
// restore Serial after it was redefined for serial connections / fine ripristino di Serial dopo definizione per connessione seriali ==============
