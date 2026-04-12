/*!
* @file  OpenFIREshared.h
* @brief Shared board assets for use between OpenFIRE microcontroller clients
*        and configuration apps for the OpenFIRE platform.
*
* @copyright That One Seong, 2025
*
*  OpenFIREshared is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef _OPENFIRESHARED_H_
#define _OPENFIRESHARED_H_

#include <string>
#include <map>
#include <unordered_map>
#include <vector>

//// BOARD IDENTIFIERS (for Desktop App identification and determining presets)

#ifdef ARDUINO_ADAFRUIT_ITSYBITSY_RP2040
    #define OPENFIRE_BOARD "adafruitItsyRP2040"
#elifdef ARDUINO_ADAFRUIT_KB2040_RP2040
    #define OPENFIRE_BOARD "adafruitKB2040"
#elifdef ARDUINO_NANO_RP2040_CONNECT
    #define OPENFIRE_BOARD "arduinoNanoRP2040"
#elifdef ARDUINO_WAVESHARE_RP2040_ZERO
    #define OPENFIRE_BOARD "waveshareZero"
#elifdef ARDUINO_YD_RP2040
    #define OPENFIRE_BOARD "vccgndYD"
#elifdef ARDUINO_RASPBERRY_PI_PICO
    #define OPENFIRE_BOARD "rpipico"
#elifdef ARDUINO_RASPBERRY_PI_PICO_W
    #define OPENFIRE_BOARD "rpipicow"
#elifdef ARDUINO_RASPBERRY_PI_PICO_2
    #define OPENFIRE_BOARD "rpipico2"
#elifdef ARDUINO_RASPBERRY_PI_PICO_2W
    #define OPENFIRE_BOARD "rpipico2w"
#elifdef ARDUINO_ESP32_S3_WROOM1_DevKitC_1_N16R8
    #define OPENFIRE_BOARD "esp32-s3-devkitc-1"
#elifdef ARDUINO_WAVESHARE_ESP32_S3_PICO
    #define OPENFIRE_BOARD "waveshare-esp32-s3-pico"
#elifdef ARDUINO_GENERIC_RP2350
    #define OPENFIRE_BOARD "generic-rp2350"
#else
    #define OPENFIRE_BOARD "generic"
#endif // board

class OF_Const
{
public:
    // Any new non-btn slots should ideally be added at the bottom, above the "count" line
    // Inputs Map indices
    enum {
        unavailable = -2,
        btnUnmapped = -1,
        btnTrigger = 0,
        btnGunA,
        btnGunB,
        btnStart,
        btnSelect,
        btnGunC,
        btnGunUp,
        btnGunDown,
        btnGunLeft,
        btnGunRight,
        btnPedal,
        btnPedal2,
        btnPump,
        btnHome,
        // ^ btn inputs
        rumblePin,
        solenoidPin,
        rumbleSwitch,
        solenoidSwitch,
        autofireSwitch,
        neoPixel,
        ledR,
        ledG,
        ledB,
        wiiClockGen,
        camSDA,
        camSCL,
        periphSDA,
        periphSCL,
        analogX,
        analogY,
        tempPin,
        // Add non-button inputs here
        boardInputsCount
    } boardInputs_e;

    // Note: names should be written in plain english,
    // as these name strings are shared with Apps
    const std::unordered_map<std::string_view, int> boardInputs_Strings = {
        {"Unmapped",            btnUnmapped     },
        {"Trigger",             btnTrigger      },
        {"Button A",            btnGunA         },
        {"Button B",            btnGunB         },
        {"Button C",            btnGunC         },
        {"Start",               btnStart        },
        {"Select",              btnSelect       },
        {"D-Pad Up",            btnGunUp        },
        {"D-Pad Down",          btnGunDown      },
        {"D-Pad Left",          btnGunLeft      },
        {"D-Pad Right",         btnGunRight     },
        {"Pedal",               btnPedal        },
        {"Alt Pedal",           btnPedal2       },
        {"Pump Action",         btnPump         },
        {"Home Button",         btnHome         },
        {"Rumble Signal",       rumblePin       },
        {"Solenoid Signal",     solenoidPin     },
        {"Rumble Switch",       rumbleSwitch    },
        {"Solenoid Switch",     solenoidSwitch  },
        {"Autofire Switch",     autofireSwitch  },
        {"External NeoPixel",   neoPixel        },
        {"RGB LED Red",         ledR            },
        {"RGB LED Green",       ledG            },
        {"RGB LED Blue",        ledB            },
        {"Wii Cam Clock",       wiiClockGen     },
        {"Camera SDA",          camSDA          },
        {"Camera SCL",          camSCL          },
        {"Peripherals SDA",     periphSDA       },
        {"Peripherals SCL",     periphSCL       },
        {"Analog Stick X",      analogX         },
        {"Analog Stick Y",      analogY         },
        {"Temperature Sensor",  tempPin         },
    };

    // For Apps to use for lists of pin functions
    char* boardInputs_sortedStr[boardInputsCount+1];

    // Constructor
    OF_Const() {
#ifdef OF_APP // generate strings list for the available board inputs
        for(auto &func : boardInputs_Strings)
            boardInputs_sortedStr[func.second+1] = (char*)func.first.data();
#endif // OF_APP
    }

    // Boolean/toggle settings indices
    enum {
        customPins = 0,
        rumble,
        solenoid,
        autofire,
        simplePause,
        holdToPause,
        commonAnode,
        lowButtonsMode,
        rumbleFF,
        invertStaticPixels,
        i2cOLED,
        i2cOLEDaltAddr,
        // Add here
        boolTypesCount
    } boolTypes_e;

    const std::unordered_map<std::string_view, int> boolTypes_Strings = {
        {"CustomPins",          customPins          },
        {"Rumble",              rumble              },
        {"Solenoid",            solenoid            },
        {"Autofire",            autofire            },
        {"SimplePause",         simplePause         },
        {"HoldToPause",         holdToPause         },
        {"LEDAnode",            commonAnode         },
        {"LowButtons",          lowButtonsMode      },
        {"RumbFFB",             rumbleFF            },
        {"InvertStaticPixels",  invertStaticPixels  },
        {"I2COLEDEnabled",      i2cOLED             },
        {"I2COLEDAltAddr",      i2cOLEDaltAddr      },
    };

    // Variable settings indices
    enum {
        rumbleStrength = 0,
        rumbleInterval,
        solenoidOnLength,
        solenoidOffLength,
        solenoidHoldLength,
        holdToPauseLength,
        customLEDcount,
        customLEDstatic,
        customLEDcolor1,
        customLEDcolor2,
        customLEDcolor3,
        tempWarning,
        tempShutdown,
        analogMode,
        // Add here
        settingsTypesCount
    } settingsTypes_e;

    const std::unordered_map<std::string_view, int> settingsTypes_Strings = {
        {"RumbPwr",             rumbleStrength      },
        {"RumbTime",            rumbleInterval      },
        {"SolOn",               solenoidOnLength    },
        {"SolOff",              solenoidOffLength   },
        {"SolHold",             solenoidHoldLength  },
        {"HoldToPauseLength",   holdToPauseLength   },
        {"CtmPixelsCount",      customLEDcount      },
        {"StaticPixels",        customLEDstatic     },
        {"StaticColor1",        customLEDcolor1     },
        {"StaticColor2",        customLEDcolor2     },
        {"StaticColor3",        customLEDcolor3     },
        {"TempWarning",         tempWarning         },
        {"TempDanger",          tempShutdown        },
        {"AnalogMode",          analogMode          },
    };

    enum {
        analogModeStick = 0,
        analogModeDpad,
        analogModeKeys
    } analogModeSettings_e;

    // Profile data type indices
    // this MUST match the order of ProfileData_s in (FW)OpenFIREprefs
    // as ProfData is accessed by struct offset.
    enum {
        profTopOffset = 0,
        profBottomOffset,
        profLeftOffset,
        profRightOffset,
        profTLled,
        profTRled,
        profAdjX,
        profAdjY,
        profIrSens,
        profRunMode,
        profIrLayout,
        profAR,
        profColor,
        profName,
        profDataTypes,
        profCurrent = 0xFD,
    } profSyncTypes_e;

    const std::unordered_map<std::string_view, int> profSettingTypes_Strings = {
        {"TopOffset",   profTopOffset       },
        {"BtmOffset",   profBottomOffset    },
        {"LftOffset",   profLeftOffset      },
        {"RhtOffset",   profRightOffset     },
        {"TLLed",       profTLled           },
        {"TRLed",       profTRled           },
        {"AdjX",        profAdjX            },
        {"AdjY",        profAdjY            },
        {"IrSens",      profIrSens          },
        {"IrRunMode",   profRunMode         },
        {"IrLayout",    profIrLayout        },
        {"AspectRatio", profAR              },
        {"Color",       profColor           },
        {"Name",        profName            },
        {"CurrentProf", profCurrent         },
    };

    // Layout types indices
    enum {
        layoutSquare = 0,
        layoutDiamond,
        // Add here
        layoutTypes
    } layoutTypes_e;

    // Aspect ratio indices
    enum {
        ar16_9 = 0,
        ar16_10,
        ar3_2,
        ar5_4,
        // Add here
        ar4_3,
        aspectRatiosCount
    } aspectRatios_e;

    /* ////
     * Shared serial control/signal codes for both boards and app.
     * For purposes of app-side debugability: ASCII 128+ should be for the board to send,
     * and invisible ASCII characters/control codes 0-32 should be for the app to send.
     *
     * ASCII 33-127 should be avoided whenever possible.
     */////
    enum {
        // Docking commands
        sDock1 = 1,
        sDock2,

        // Mode toggles from app
        sIRTest = 5,
        sCaliProfile,
        sCaliStart,
        sCaliSens,
        sCaliLayout,

        // Test signals from app
        sTestSolenoid = 15,
        sTestRumble,
        sTestLEDR,
        sTestLEDG,
        sTestLEDB,

        // Error types from board (with sError)
        sErrCam = 0x80, // 128
        sErrPeriphGeneric,

        // Status updates from board
        sBtnPressed = 0x90, // 144
        sBtnReleased,
        sAnalogPosUpd,
        sTemperatureUpd,
        sCaliStageUpd,
        sCaliInfoUpd,
        sTestCoords,
        sCurrentProf,

        // Push settings to board
        sCommitStart = 0xAA, // 170
        sCommitToggles,
        sCommitPins,
        sCommitSettings,
        sCommitProfile,
        sCommitBtns,
        sCommitID,

        // Grab settings from board
        sGetPins = 0xC8, // 200
        sGetToggles,
        sGetSettings,
        sGetProfile,
        sGetBtns,

        // for non-RP2040 boards that don't have a magic number-type reset
        sRebootToBootloader = 0xF0, // 245

        sError = 0xFA, // 250
        sSave = 0xFC, // 252
        sClearFlash = 0xFD, // 253
        // Terminates out of any current mode, or undocks
        serialTerminator = 0xFE // 254
    } serialCmdTypes_e;

    enum {
        usbPID = 0,
        usbName,
    } usbIdSyncTypes_e;

    /// @brief      Map of default pin mappings for each supported board
    /// @details    Key = board, int array maps to RP2040 GPIO where each value is a FW function (or unmapped).
    /// @note       /*xx*/ indicates the number of the GPIO, e.g. /*02*/ for GPIO-02
    const std::unordered_map<std::string_view, std::vector<int>> boardsPresetsMap = {
        //=====================================================================================================================
        // Raspberry Pi Pico
        // Board Type: RP2040
        // Notes: Raspberry Pi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipico",                 {/*00*/ btnGunA,        btnGunB,        btnGunC,        btnStart,       btnSelect,
                                     /*05*/ btnHome,        btnGunUp,       btnGunDown,     btnGunLeft,     btnGunRight,
                                     /*10*/ ledR,           ledG,           ledB,           btnPump,        btnPedal,
                                     /*15*/ btnTrigger,     solenoidPin,    rumblePin,      periphSDA,      periphSCL,
                                     /*20*/ camSDA,         camSCL,         btnUnmapped,    unavailable,    unavailable,
                                     /*25*/ unavailable,    btnUnmapped,    btnUnmapped,    tempPin,        unavailable     }},
        //=====================================================================================================================
        // Raspberry Pi Pico W
        // Board Type: RP2040
        // Notes: Raspberry Pi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipicow",                {/*00*/ btnGunA,        btnGunB,        btnGunC,        btnStart,       btnSelect,
                                     /*05*/ btnHome,        btnGunUp,       btnGunDown,     btnGunLeft,     btnGunRight,
                                     /*10*/ periphSDA,      periphSCL,      btnUnmapped,    btnPump,        btnPedal,
                                     /*15*/ btnTrigger,     solenoidPin,    rumblePin,      btnUnmapped,    btnUnmapped,
                                     /*20*/ camSDA,         camSCL,         btnUnmapped,    unavailable,    unavailable,
                                     /*25*/ unavailable,    analogY,        analogX,        tempPin,        unavailable     }},
        //=====================================================================================================================
        // Raspberry Pi Pico 2
        // Board Type: RP2350A
        // Notes: Pico 2 boards are pin-identical to Pico 1 boards
        {"rpipico2",                {/*00*/ btnGunA,        btnGunB,        btnGunC,        btnStart,       btnSelect,
                                     /*05*/ btnHome,        btnGunUp,       btnGunDown,     btnGunLeft,     btnGunRight,
                                     /*10*/ ledR,           ledG,           ledB,           btnPump,        btnPedal,
                                     /*15*/ btnTrigger,     solenoidPin,    rumblePin,      periphSDA,      periphSCL,
                                     /*20*/ camSDA,         camSCL,         btnUnmapped,    unavailable,    unavailable,
                                     /*25*/ unavailable,    btnUnmapped,    btnUnmapped,    tempPin,        unavailable     }},
        //=====================================================================================================================
        // Raspberry Pi Pico 2W
        // Board Type: RP2350A
        // Notes: Pico 2 boards are pin-identical to Pico 1 boards
        {"rpipico2w",               {/*00*/ btnGunA,        btnGunB,        btnGunC,        btnStart,       btnSelect,
                                     /*05*/ btnHome,        btnGunUp,       btnGunDown,     btnGunLeft,     btnGunRight,
                                     /*10*/ ledR,           ledG,           ledB,           btnPump,        btnPedal,
                                     /*15*/ btnTrigger,     solenoidPin,    rumblePin,      periphSDA,      periphSCL,
                                     /*20*/ camSDA,         camSCL,         btnUnmapped,    unavailable,    unavailable,
                                     /*25*/ unavailable,    btnUnmapped,    btnUnmapped,    tempPin,        unavailable     }},
        //=====================================================================================================================
        // Adafruit ItsyBitsy RP2040
        // Board Type: RP2040
        // Notes: pins 13-17 & 21-23 are unexposed
        {"adafruitItsyRP2040",      {/*00*/ btnUnmapped,    btnUnmapped,    camSDA,         camSCL,         btnPedal,
                                     /*05*/ btnUnmapped,    btnTrigger,     btnGunDown,     btnGunLeft,     btnGunUp,
                                     /*10*/ btnGunRight,    btnGunC,        btnUnmapped,    unavailable,    unavailable,
                                     /*15*/ unavailable,    unavailable,    unavailable,    btnUnmapped,    btnUnmapped,
                                     /*20*/ btnUnmapped,    unavailable,    unavailable,    unavailable,    rumblePin,
                                     /*25*/ solenoidPin,    btnGunB,        btnGunA,        btnStart,       btnSelect       }},
        //=====================================================================================================================
        // Adafruit "Kee Boar" KB2040
        // Board Type: RP2040
        // Notes: pins 11-17 & 21-25 are unexposed
        {"adafruitKB2040",          {/*00*/ btnUnmapped,    btnUnmapped,    camSDA,         camSCL,         btnGunB,
                                     /*05*/ rumblePin,      btnGunC,        solenoidPin,    btnSelect,      btnStart,
                                     /*10*/ btnGunRight,    unavailable,    unavailable,    unavailable,    unavailable,
                                     /*15*/ unavailable,    unavailable,    unavailable,    btnGunUp,       btnGunLeft,
                                     /*20*/ btnGunDown,     unavailable,    unavailable,    unavailable,    unavailable,
                                     /*25*/ unavailable,    tempPin,        btnHome,        btnTrigger,     btnGunA         }},
        //=====================================================================================================================
        // Arduino Nano RP2040 Connect
        // Notes: pins 2-3, 8-11, 14, & 22-24 are unexposed;
        //        some other pins are analog, but controlled by NINA and thus are unavailable in OpenFIRE for the moment
        {"arduinoNanoRP2040",       {/*00*/ btnTrigger,     btnPedal,       unavailable,    unavailable,    btnGunA,
                                     /*05*/ btnGunC,        btnUnmapped,    btnGunB,        unavailable,    unavailable,
                                     /*10*/ unavailable,    unavailable,    camSDA,         camSCL,         unavailable,
                                     /*15*/ btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                     /*20*/ btnUnmapped,    btnUnmapped,    unavailable,    unavailable,    unavailable,
                                     /*25*/ btnUnmapped,    btnUnmapped,    btnUnmapped,    tempPin,        btnUnmapped     }},
        //=====================================================================================================================
        // Waveshare Zero RP2040
        // Note: pin 16 is reserved for the board's builtin NeoPixel (not currently used?);
        //       pins 17-25 are underside pads which are not exposed in the app for layout reasons;
        {"waveshareZero",           {/*00*/ btnTrigger,     btnGunA,        btnGunB,        btnGunC,        btnStart,
                                     /*05*/ btnSelect,      btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                     /*10*/ btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,    camSDA,
                                     /*15*/ camSCL,         unavailable,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                     /*20*/ btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                     /*25*/ btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,    tempPin         }},
        //=====================================================================================================================
        /* more RP-series boards should be added here */
        //=====================================================================================================================
        // Espressif ESP32 S3 WROOM-1 DevkitC-1 N16R8
        // Board Type: ESP32
        {"esp32-s3-devkitc-1",      {/*00*/ unavailable,    btnTrigger,     btnGunRight,    btnUnmapped,    analogX,
                                     /*05*/ analogY,        tempPin,        btnUnmapped,    camSDA,         camSCL,
                                     /*10*/ btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,    neoPixel,
                                     /*15*/ periphSCL,      rumblePin,      solenoidPin,    periphSDA,      btnUnmapped,
                                     /*20*/ btnUnmapped,    btnGunC,        unavailable,    unavailable,    unavailable,
                                     /*25*/ unavailable,    unavailable,    unavailable,    unavailable,    unavailable,
                                     /*30*/ unavailable,    unavailable,    unavailable,    unavailable,    unavailable,
                                     /*35*/ btnHome,        btnGunA,        btnGunB,        btnSelect,      btnStart,
                                     /*40*/ btnGunUp,       btnGunDown,     btnGunLeft,     unavailable,    unavailable,
                                     /*45*/ btnPump,        btnPedal2,      btnPedal,       unavailable                        }},
        //=====================================================================================================================
        // Waveshare ESP32 S3 Pico
        // Board Type: ESP32
        {"waveshare-esp32-s3-pico", {/*00*/ unavailable,    btnUnmapped,    btnUnmapped,    unavailable,    camSDA,
                                     /*05*/ camSCL,         btnPedal2,      analogY,        analogX,        tempPin,
                                     /*10*/ btnUnmapped,    btnGunA,        btnGunB,        btnGunC,        btnStart,
                                     /*15*/ btnSelect,      btnHome,        btnGunUp,       btnGunDown,     unavailable,
                                     /*20*/ unavailable,    unavailable,    unavailable,    unavailable,    unavailable,
                                     /*25*/ unavailable,    unavailable,    unavailable,    unavailable,    unavailable,
                                     /*30*/ unavailable,    unavailable,    unavailable,    btnGunLeft,     btnGunRight,
                                     /*35*/ periphSDA,      periphSCL,      neoPixel,       btnPump,        btnPedal,
                                     /*40*/ btnTrigger,     rumblePin,      solenoidPin,    unavailable,    unavailable,
                                     /*45*/ unavailable,    unavailable,    unavailable,    unavailable                     }},
        //=====================================================================================================================
        /* more ESP boards should be added here */
    };

    static const unsigned int TEMPERATURE_SENSOR_ERROR_VALUE = 125; // ADC reading indicating sensor fault / disconnection.

// Only needed for the Desktop App, don't build for microcontroller firmware!
#ifdef OF_APP

    const std::map<std::string_view, const char *> boardNames = {
        {"rpipico",                 "Raspberry Pi Pico (RP2040)"},
        {"rpipicow",                "Raspberry Pi Pico W (RP2040)"},
        {"rpipico2",                "Raspberry Pi Pico 2 (RP2350)"},
        {"rpipico2w",               "Raspberry Pi Pico 2W (RP2350)"},
        {"adafruitItsyRP2040",      "Adafruit ItsyBitsy RP2040"},
        {"adafruitKB2040",          "Adafruit Keeboar KB2040"},
        {"arduinoNanoRP2040",       "Arduino Nano Connect RP2040"},
        {"waveshareZero",           "Waveshare Zero RP2040"},
        {"esp32-s3-devkitc-1",      "ESP32-S3 WROOM-1 DevkitC-1 (N16R8)"},
        {"waveshare-esp32-s3-pico", "Waveshare ESP32-S3-Pico"},
        // Add more here!
        {"generic-rp2350",          "Unknown RP2350 Board"},
        {"generic",                 "Unknown RP2040 Board"}
    };

    /// @brief      Types of board architectures
    /// @details    Board archs are to be dictated by the application on a per-board basis,
    ///             to be used for defining which pins are capable of what.
    const char* boardArchs[2] = {
        "rp2040_235X",
        "esp32-s3"
    };

    /// @brief      Indices for the boardArchs strings above
    enum {
        boardRP,
        boardESP32_S3
    } boardArchs_e;

    /// @brief      Indices to be used as a bitmap for defining pin capabilities
    /// @note       Default (0) assumes the pin is digital only with no I2C or SPI capability
    enum {
        pinDigital = 0,          // macro for no special function at all
        pinSystem  = 0,          // macro for pin that shouldn't be exposed to the user
        pinHasADC  = 0b00000001, // whether pin is connected to an ADC or not
        pinI2C0SDA = 0b00000010, // pin has I2C0 SDA
        pinI2C0SCL = 0b00000110, // pin has I2C0 SCL
        pinI2C1SDA = 0b00001010, // pin has I2C1 SDA
        pinI2C1SCL = 0b00001110, // pin has I2C1 SCL
        pinSPI0RX  = 0b00100000, // pin has SPI0 RX
        pinSPI0TX  = 0b01000000, // pin has SPI0 TX
        pinSPI0SCK = 0b01100000, // pin has SPI0 SCK
        pinSPI0CSn = 0b10000000, // pin has SPI0 CSn
        pinSPI1RX  = 0b00110000, // pin has SPI1 RX
        pinSPI1TX  = 0b01010000, // pin has SPI1 TX
        pinSPI1SCK = 0b01110000, // pin has SPI1 SCK
        pinSPI1CSn = 0b10010000, // pin has SPI1 CSn
        // ESP can software define I2C & SPI pins and their channels
        // So just let them map any function from the app.
        pinAnyI2C  = 0b00000001 << 8,
        pinAnySPI  = 0b00000010 << 8,

        // below are bitmasks to be used in Application
        pinCanSPI   = 0b11110000,
        pinCanI2C   = 1 << 1,
        // check if SCL, else SDA
        pinIsI2CSCL = 1 << 2,
        // check if I2C1, else I2C0
        pinIsI2C1   = 1 << 3,
        // check if SPI1, else SPI0
        pinIsSPI1   = 1 << 4
    } pinCapabilities_e;

    /// @brief      Map of capabilities of each pin for a board type
    /// @details    Dictates what types of functions a pin can be mapped to, based on its capabilities
    ///             This applies to ALL boards using a specific architecture.
    const std::unordered_map<std::string_view, std::vector<int>> mcuCapableMaps = {
        //====================================================
        // Base Microcontroller: RP2040 & RP235X(A|B)
        // GPIO: 30(RP2040/RP2350A) / 48(RP2350B)
        // Notes: Because RP235X-series shares the lower 30 GPIO, both RPI MCUs can share the same map
        {boardArchs[boardRP],   {/*00*/ pinI2C0SDA | pinSPI0RX,                     pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI0RX,
                                 /*05*/ pinI2C0SCL | pinSPI0CSn,                    pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI1RX,                 pinI2C0SCL | pinSPI1CSn,
                                 /*10*/ pinI2C1SDA | pinSPI1SCK,                    pinI2C1SCL | pinSPI1TX,                 pinI2C0SDA | pinSPI1RX,                 pinI2C0SCL | pinSPI1CSn,                pinI2C1SDA | pinSPI1SCK,
                                 /*15*/ pinI2C1SCL | pinSPI1TX,                     pinI2C0SDA | pinSPI0RX,                 pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,
                                 /*20*/ pinI2C0SDA | pinSPI0RX,                     pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI1RX,
                                 /*25*/ pinI2C0SCL | pinSPI1CSn,                    pinI2C1SDA | pinSPI1SCK | pinHasADC,    pinI2C1SCL | pinSPI1TX  | pinHasADC,    pinI2C0SDA | pinSPI1RX | pinHasADC,     pinI2C0SCL | pinSPI1CSn | pinHasADC,
                                 /*30*/ pinI2C1SDA | pinSPI1SCK,                    pinI2C1SCL | pinSPI1TX,                 pinI2C0SDA | pinSPI0RX,                 pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,
                                 /*35*/ pinI2C1SCL | pinSPI0TX,                     pinI2C0SDA | pinSPI0RX,                 pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,
                                 /*40*/ pinI2C0SDA | pinSPI1RX  | pinHasADC,        pinI2C0SCL | pinSPI1CSn | pinHasADC,    pinI2C1SDA | pinSPI1SCK | pinHasADC,    pinI2C1SCL | pinSPI1TX | pinHasADC,     pinI2C0SDA | pinSPI1RX  | pinHasADC,
                                 /*45*/ pinI2C0SCL | pinSPI1CSn | pinHasADC,        pinI2C1SDA | pinSPI1SCK | pinHasADC,    pinI2C1SCL | pinSPI1TX  | pinHasADC                                                                                     }},
        //====================================================
        // Base Microcontroller: ESP32-S3
        // GPIO: 49
        // Notes: Only ADC1 (GP1-10) can reliably be used for user-mappable analog reading.
        //        SPI0/1 (GP26-32) are reserved for the system, so only SPI2/3 are available for user mapping.
        //        GP0, 3, & 45-46 are also system-reserved and shouldn't be exposed.
        {boardArchs[boardESP32_S3],{/*00*/ pinSystem,                               pinHasADC | pinAnyI2C | pinAnySPI,      pinHasADC | pinAnyI2C | pinAnySPI,      pinSystem,                              pinHasADC | pinAnyI2C | pinAnySPI,
                                    /*05*/ pinHasADC | pinAnyI2C | pinAnySPI,       pinHasADC | pinAnyI2C | pinAnySPI,      pinHasADC | pinAnyI2C | pinAnySPI,      pinHasADC | pinAnyI2C | pinAnySPI,      pinHasADC | pinAnyI2C | pinAnySPI,
                                    /*10*/ pinHasADC | pinAnyI2C | pinAnySPI,       pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,
                                    /*15*/ pinAnyI2C | pinAnySPI,                   pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,
                                    /*20*/ pinAnyI2C | pinAnySPI,                   pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,
                                    /*25*/ pinAnyI2C | pinAnySPI,                   pinSystem,                              pinSystem,                              pinSystem,                              pinSystem,
                                    /*30*/ pinSystem,                               pinSystem,                              pinSystem,                              pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,
                                    /*35*/ pinAnyI2C | pinAnySPI,                   pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,
                                    /*40*/ pinAnyI2C | pinAnySPI,                   pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI,
                                    /*45*/ pinSystem,                               pinSystem,                              pinAnyI2C | pinAnySPI,                  pinAnyI2C | pinAnySPI                                                           }},
        //====================================================
        // Board Overrides: Raspberry Pi Pico (Non-/W, 1&2)
        // Some pins that should have I2C or SPI functions apparently aren't allowed on rpipico[2](w)?
        {"rpipico",             {/*00*/ pinI2C0SDA | pinSPI0RX,     pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI0RX,
                                 /*05*/ pinI2C0SCL | pinSPI0CSn,    pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI1RX,                 pinI2C0SCL | pinSPI1CSn,
                                 /*10*/ pinI2C1SDA | pinSPI1SCK,    pinI2C1SCL | pinSPI1TX,                 pinI2C0SDA | pinSPI1RX,                 pinI2C0SCL | pinSPI1CSn,                pinI2C1SDA | pinSPI1SCK,
                                 /*15*/ pinI2C1SCL | pinSPI1TX,     pinI2C0SDA | pinSPI0RX,                 pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,
                                 /*20*/ pinI2C0SDA | pinSPI0RX,     pinI2C0SCL | pinSPI0CSn,                pinDigital,                             pinDigital,                             pinDigital,
                                 /*25*/ pinDigital,                 pinI2C1SDA | pinSPI1SCK | pinHasADC,    pinI2C1SCL | pinSPI1TX | pinHasADC,     pinSPI1RX  | pinHasADC,                 pinDigital                  }},
        {"rpipicow",            {/*00*/ pinI2C0SDA | pinSPI0RX,     pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI0RX,
                                 /*05*/ pinI2C0SCL | pinSPI0CSn,    pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI1RX,                 pinI2C0SCL | pinSPI1CSn,
                                 /*10*/ pinI2C1SDA | pinSPI1SCK,    pinI2C1SCL | pinSPI1TX,                 pinI2C0SDA | pinSPI1RX,                 pinI2C0SCL | pinSPI1CSn,                pinI2C1SDA | pinSPI1SCK,
                                 /*15*/ pinI2C1SCL | pinSPI1TX,     pinI2C0SDA | pinSPI0RX,                 pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,
                                 /*20*/ pinI2C0SDA | pinSPI0RX,     pinI2C0SCL | pinSPI0CSn,                pinDigital,                             pinDigital,                             pinDigital,
                                 /*25*/ pinDigital,                 pinI2C1SDA | pinSPI1SCK | pinHasADC,    pinI2C1SCL | pinSPI1TX | pinHasADC,     pinSPI1RX  | pinHasADC,                 pinDigital                  }},
        {"rpipico2",            {/*00*/ pinI2C0SDA | pinSPI0RX,     pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI0RX,
                                 /*05*/ pinI2C0SCL | pinSPI0CSn,    pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI1RX,                 pinI2C0SCL | pinSPI1CSn,
                                 /*10*/ pinI2C1SDA | pinSPI1SCK,    pinI2C1SCL | pinSPI1TX,                 pinI2C0SDA | pinSPI1RX,                 pinI2C0SCL | pinSPI1CSn,                pinI2C1SDA | pinSPI1SCK,
                                 /*15*/ pinI2C1SCL | pinSPI1TX,     pinI2C0SDA | pinSPI0RX,                 pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,
                                 /*20*/ pinI2C0SDA | pinSPI0RX,     pinI2C0SCL | pinSPI0CSn,                pinDigital,                             pinDigital,                             pinDigital,
                                 /*25*/ pinDigital,                 pinI2C1SDA | pinSPI1SCK | pinHasADC,    pinI2C1SCL | pinSPI1TX | pinHasADC,     pinSPI1RX  | pinHasADC,                 pinDigital                  }},
        {"rpipico2w",           {/*00*/ pinI2C0SDA | pinSPI0RX,     pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI0RX,
                                 /*05*/ pinI2C0SCL | pinSPI0CSn,    pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,                 pinI2C0SDA | pinSPI1RX,                 pinI2C0SCL | pinSPI1CSn,
                                 /*10*/ pinI2C1SDA | pinSPI1SCK,    pinI2C1SCL | pinSPI1TX,                 pinI2C0SDA | pinSPI1RX,                 pinI2C0SCL | pinSPI1CSn,                pinI2C1SDA | pinSPI1SCK,
                                 /*15*/ pinI2C1SCL | pinSPI1TX,     pinI2C0SDA | pinSPI0RX,                 pinI2C0SCL | pinSPI0CSn,                pinI2C1SDA | pinSPI0SCK,                pinI2C1SCL | pinSPI0TX,
                                 /*20*/ pinI2C0SDA | pinSPI0RX,     pinI2C0SCL | pinSPI0CSn,                pinDigital,                             pinDigital,                             pinDigital,
                                 /*25*/ pinDigital,                 pinI2C1SDA | pinSPI1SCK | pinHasADC,    pinI2C1SCL | pinSPI1TX | pinHasADC,     pinSPI1RX  | pinHasADC,                 pinDigital                  }},
        };

    enum {
        posNothing  = 0,
        posLeft     = 0b00000001 << 8,
        posRight    = 0b00000010 << 8,
        posMiddle   = 0b00000100 << 8,
        posCheck    = posLeft | posRight | posMiddle
    } boardBoxPositions_e;

    /// @brief      Map of graphical placement for each pin in the application
    /// @details    Key = board, int vector maps to microcontroller GPIO.
    ///             Each pin should be a combination of grid layout slot it should be in,
    ///             OR'd by the grid it should belong to.
    ///             Unexposed pins should use only posNothing (0).
    ///             (Yep, bitpacking! Three least significant bits of the second byte determine left/right/under position)
    ///             (If anyone is aware of a better way of doing this, please let me know/send a PR!)
    const std::unordered_map<std::string_view, std::vector<unsigned int>> boardsBoxPositions = {
        //=====================================================================================================================
        // Raspberry Pi Pico: 15 pins left, rest of the pins right. Mostly linear order save for the reserved pins.
        // Board Type: RP2040
        // Notes: Raspberry Pi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipico",                 {/*00*/ 1  | posLeft,   2  | posLeft,   4  | posLeft,   5  | posLeft,   6  | posLeft,
                                     /*05*/ 7  | posLeft,   9  | posLeft,   10 | posLeft,   11 | posLeft,   12 | posLeft,
                                     /*10*/ 14 | posLeft,   15 | posLeft,   16 | posLeft,   17 | posLeft,   19 | posLeft,
                                     /*15*/ 20 | posLeft,   20 | posRight,  19 | posRight,  17 | posRight,  16 | posRight,
                                     /*20*/ 15 | posRight,  14 | posRight,  12 | posRight,    posNothing,     posNothing,
                                     /*25*/   posNothing,   10 | posRight,  9  | posRight,  7  | posRight,    posNothing    }},
        //=====================================================================================================================
        // Raspberry Pi Pico W: same as non-W Pico.
        // Board Type: RP2040
        // Notes: Raspberry Pi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipicow",                {/*00*/ 1  | posLeft,   2  | posLeft,   4  | posLeft,   5  | posLeft,   6  | posLeft,
                                     /*05*/ 7  | posLeft,   9  | posLeft,   10 | posLeft,   11 | posLeft,   12 | posLeft,
                                     /*10*/ 14 | posLeft,   15 | posLeft,   16 | posLeft,   17 | posLeft,   19 | posLeft,
                                     /*15*/ 20 | posLeft,   20 | posRight,  19 | posRight,  17 | posRight,  16 | posRight,
                                     /*20*/ 15 | posRight,  14 | posRight,  12 | posRight,    posNothing,     posNothing,
                                     /*25*/   posNothing,   10 | posRight,  9  | posRight,  7  | posRight,    posNothing    }},
        //=====================================================================================================================
        // Raspberry Pi Pico 2: same as Pico 1
        // Board Type: RP2350A
        // Notes: Raspberry Pi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipico2",                {/*00*/ 1  | posLeft,   2  | posLeft,   4  | posLeft,   5  | posLeft,   6  | posLeft,
                                     /*05*/ 7  | posLeft,   9  | posLeft,   10 | posLeft,   11 | posLeft,   12 | posLeft,
                                     /*10*/ 14 | posLeft,   15 | posLeft,   16 | posLeft,   17 | posLeft,   19 | posLeft,
                                     /*15*/ 20 | posLeft,   20 | posRight,  19 | posRight,  17 | posRight,  16 | posRight,
                                     /*20*/ 15 | posRight,  14 | posRight,  12 | posRight,    posNothing,     posNothing,
                                     /*25*/   posNothing,   10 | posRight,  9  | posRight,  7  | posRight,    posNothing    }},
        //=====================================================================================================================
        // Raspberry Pi Pico 2W: same as non-W Pico 1&2.
        // Board Type: RP2350A
        // Notes: Raspberry Pi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipico2w",               {/*00*/ 1  | posLeft,   2  | posLeft,   4  | posLeft,   5  | posLeft,   6  | posLeft,
                                     /*05*/ 7  | posLeft,   9  | posLeft,   10 | posLeft,   11 | posLeft,   12 | posLeft,
                                     /*10*/ 14 | posLeft,   15 | posLeft,   16 | posLeft,   17 | posLeft,   19 | posLeft,
                                     /*15*/ 20 | posLeft,   20 | posRight,  19 | posRight,  17 | posRight,  16 | posRight,
                                     /*20*/ 15 | posRight,  14 | posRight,  12 | posRight,    posNothing,     posNothing,
                                     /*25*/   posNothing,   10 | posRight,  9  | posRight,  7  | posRight,    posNothing    }},
        //=====================================================================================================================
        // Adafruit ItsyBitsy RP2040: A very cluttered and kind of unfriendly layout tbh :(
        // Board Type: RP2040
        // Notes: pins 13-17 & 21-23 are unexposed
        {"adafruitItsyRP2040",      {/*00*/ 13 | posRight,  14 | posRight,  12 | posRight,  11 | posRight,  2  | posMiddle,
                                     /*05*/ 1  | posMiddle, 9  | posRight,  8  | posRight,  7  | posRight,  6  | posRight,
                                     /*10*/ 5  | posRight,  4  | posRight,  14 | posLeft,     posNothing,     posNothing,
                                     /*15*/   posNothing,     posNothing,     posNothing,   11 | posLeft,   12 | posLeft,
                                     /*20*/ 13 | posLeft,     posNothing,     posNothing,     posNothing,   9  | posLeft,
                                     /*25*/ 10 | posLeft,   5  | posLeft,   6  | posLeft,   7  | posLeft,   8  | posLeft    }},
        //=====================================================================================================================
        // Adafruit KB2040: Like the Itsy with more padding.
        // Board Type: RP2040
        // Notes: pins 11-17 & 21-25 are unexposed
        {"adafruitKB2040",          {/*00*/ 3  | posLeft,   4  | posLeft,   7  | posLeft,   8  | posLeft,   9  | posLeft,
                                     /*05*/ 10 | posLeft,   11 | posLeft,   12 | posLeft,   13 | posLeft,   14 | posLeft,
                                     /*10*/ 14 | posRight,    posNothing,     posNothing,     posNothing,     posNothing,
                                     /*15*/   posNothing,     posNothing,     posNothing,   11 | posRight,  13 | posRight,
                                     /*20*/ 12 | posRight,    posNothing,     posNothing,     posNothing,     posNothing,
                                     /*25*/   posNothing,   10 | posRight,  9  | posRight,  8  | posRight,  7  | posRight   }},
        //=====================================================================================================================
        // Arduino Nano RP2040
        // Board Type: RP2040
        // Notes: pins 2-3, 8-11, 14, & 22-24 are unexposed;
        //        some other pins are analog, but controlled by NINA and thus are unavailable in OpenFIRE for the moment
        {"arduinoNanoRP2040",       {/*00*/ 18 | posRight,  17 | posRight,    posNothing,     posNothing,   4  | posRight,
                                     /*05*/ 6  | posRight,  4  | posLeft,   5  | posRight,    posNothing,     posNothing,
                                     /*10*/   posNothing,     posNothing,   11 | posLeft,   12 | posLeft,     posNothing,
                                     /*15*/ 13 | posRight,  12 | posRight,  11 | posRight,  10 | posRight,  9  | posRight,
                                     /*20*/ 8  | posRight,  7  | posRight,    posNothing,     posNothing,     posNothing,
                                     /*25*/ 14 | posRight,  7  | posLeft,   8  | posLeft,   9  | posLeft,   10 | posLeft    }},
        //=====================================================================================================================
        // Waveshare Zero RP2040: Clockwise layout
        // Board Type: RP2040
        // Note: pin 16 is reserved for the board's builtin NeoPixel (not currently used?);
        //       pins 17-25 are underside pads which are not exposed in the app for layout reasons;
        {"waveshareZero",           {/*00*/ 2  | posRight,  3  | posRight,  4  | posRight,  5  | posRight,  6  | posRight,
                                     /*00*/ 7  | posRight,  8  | posRight,  9  | posRight,  10 | posRight,  11 | posRight,
                                     /*00*/ 3  | posMiddle, 2  | posMiddle, 11 | posLeft,   10 | posLeft,   9  | posLeft,
                                     /*00*/ 8  | posLeft,     posNothing,     posNothing,     posNothing,     posNothing,
                                     /*00*/   posNothing,     posNothing,     posNothing,     posNothing,     posNothing,
                                     /*00*/   posNothing,   7  | posLeft,   6  | posLeft,   5  | posLeft,   4  | posLeft    }},
        //=====================================================================================================================
        // Insert new ESP layouts below this one!
        // Feel free to use any of the above as a template.
        // ***


        //=====================================================================================================================
        // Espressif ESP32 S3 WROOM-1 DevkitC-1 N16R8
        // Board Type: ESP32-S3
        {"esp32-s3-devkitc-1",      {/*00*/   posNothing,   19 | posLeft,   18 | posLeft,  10  | posRight,   19 | posRight,
                                     /*05*/ 18 | posRight,  17 | posRight, 16  | posRight,   11 | posRight,  8  | posRight,
                                     /*10*/ 7  | posRight,  6  | posRight, 5  | posRight,   4  | posRight,  3  | posRight,
                                     /*15*/ 15 | posRight,  14 | posRight,  13 | posRight,  12 | posRight,    3  | posLeft,
                                     /*20*/  4  | posLeft,   5  | posLeft,     posNothing,     posNothing,     posNothing,
                                     /*25*/   posNothing,     posNothing,     posNothing,     posNothing,     posNothing,
                                     /*30*/   posNothing,     posNothing,     posNothing,     posNothing,     posNothing,
                                     /*35*/ 10 | posLeft,   11 | posLeft,   12 | posLeft,   13 | posLeft,   14 | posLeft,
                                     /*40*/ 15 | posLeft,   16 | posLeft,   17 | posLeft,     posNothing,     posNothing,
                                     /*45*/ 8  | posLeft,  9  | posRight,   6  | posLeft,    posNothing                    }},
        //=====================================================================================================================
        // Waveshare ESP32 S3 Pico
        // Board Type: ESP32-S3
        {"waveshare-esp32-s3-pico", {/*00*/   posNothing,   17 | posRight,  16 | posRight,    posNothing,   15 | posRight,
                                     /*05*/ 14 | posRight,  12 | posRight,  10 | posRight,  9  | posRight,  7  | posRight,
                                     /*10*/ 6  | posRight,  1  | posLeft,   2  | posLeft,   4  | posLeft,   5  | posLeft,
                                     /*15*/ 6  | posLeft,   7  | posLeft,   9  | posLeft,   10 | posLeft,     posNothing,
                                     /*20*/   posNothing,     posNothing,     posNothing,     posNothing,     posNothing,
                                     /*25*/   posNothing,     posNothing,     posNothing,     posNothing,     posNothing,
                                     /*30*/   posNothing,     posNothing,     posNothing,   11 | posLeft,   12 | posLeft,
                                     /*35*/ 14 | posLeft,   15 | posLeft,   16 | posLeft,   17 | posLeft,   19 | posLeft,
                                     /*40*/ 20 | posLeft,   19 | posRight,  20 | posRight,    posNothing,     posNothing,
                                     /*45*/   posNothing,     posNothing,     posNothing,     posNothing                    }},
        //=====================================================================================================================
        // Insert new ESP layouts below this one!
        // Feel free to use any of the above as a template.
        // ***


        //=====================================================================================================================
        // Generic (RP2040/2350) layout
        // Just reveal all pins; user assumes full responsibility if something goes wrong here
        {"generic",                 {/*00*/ 1  | posLeft,   2  | posLeft,   3  | posLeft,   4  | posLeft,   5  | posLeft,
                                     /*05*/ 6  | posLeft,   7  | posLeft,   8  | posLeft,   9  | posLeft,   10 | posLeft,
                                     /*10*/ 11 | posLeft,   12 | posLeft,   13 | posLeft,   14 | posLeft,   15 | posLeft,
                                     /*15*/ 16 | posLeft,   16 | posRight,  15 | posRight,  14 | posRight,  13 | posRight,
                                     /*20*/ 12 | posRight,  11 | posRight,  10 | posRight,  9  | posRight,  8  | posRight,
                                     /*25*/ 7  | posRight,  6  | posRight,  5  | posRight,  4  | posRight,  3  | posRight   }}
    };

    typedef struct {
        const char * name;
        std::vector<int> pin;
    } boardAltPresetsMap_t;

    /// @brief      MultiMap of alternative pin mappings for supported boards to show in the application.
    /// @details    Key = board (one board can be multiple), string literal label, int array maps to RP2040 GPIO where each value is a FW function (or unmapped).
    const std::unordered_multimap<std::string_view, boardAltPresetsMap_t> boardsAltPresets = {
        //=====================================================================================================
        // Raspberry Pi Pico Presets (currently a test)
        // Notes: rpi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        /*
        {"rpipico",                 {"Test",
                                    {   btnPump,       btnPedal,       btnUnmapped,    btnUnmapped,    btnUnmapped,
                                        btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                        btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                        btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                        btnUnmapped,   btnUnmapped,    btnUnmapped,    unavailable,    unavailable,
                                        unavailable,   btnUnmapped,    btnUnmapped,    btnUnmapped,    unavailable}}},

        {"rpipico",                 {"Test 2",
                                    {   btnGunA,       btnTrigger,     btnUnmapped,    btnUnmapped,    btnUnmapped,
                                        btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                        btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                        btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                        btnUnmapped,   btnUnmapped,    btnUnmapped,    unavailable,    unavailable,
                                        unavailable,   btnUnmapped,    btnUnmapped,    btnUnmapped,    unavailable}}},
        */

        //=====================================================================================================
        // Adafruit ItsyBitsy RP2040 Presets
        // Notes: pins 13-17 & 21-23 are unexposed
        {"adafruitItsyRP2040",      {"SAMCO 2.0 (Btn C as Home)",
                                    {   btnUnmapped,   btnUnmapped,    camSDA,         camSCL,         btnPedal,
                                        btnUnmapped,   btnTrigger,     btnGunDown,     btnGunLeft,     btnGunUp,
                                        btnGunRight,   btnHome,        btnUnmapped,    unavailable,    unavailable,
                                        unavailable,   unavailable,    unavailable,    btnUnmapped,    btnUnmapped,
                                        btnUnmapped,   unavailable,    unavailable,    unavailable,    rumblePin,
                                        solenoidPin,   btnGunB,        btnGunA,        btnStart,       btnSelect}}},

        {"adafruitItsyRP2040",      {"SAMCO 1.1",
                                    {   btnUnmapped,   btnUnmapped,    camSDA,         camSCL,         btnUnmapped,
                                        btnUnmapped,   btnGunA,        btnGunB,        rumblePin,      btnHome,
                                        btnTrigger,    btnUnmapped,    btnUnmapped,    unavailable,    unavailable,
                                        unavailable,   unavailable,    unavailable,    btnUnmapped,    btnUnmapped,
                                        btnUnmapped,   unavailable,    unavailable,    unavailable,    btnUnmapped,
                                        btnUnmapped,   btnUnmapped,    btnPedal,       btnUnmapped,    btnUnmapped}}},
    };

#endif
};


#endif // _OPENFIRESHARED_H_
