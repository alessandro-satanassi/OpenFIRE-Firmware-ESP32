/*!
* @file OpenFIREshared.h
* @brief Shared board assets.
*
* @copyright That One Seong, 2024
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

#include <cstdint>
#include <string>
#include <unordered_map>

// just to detect we're using Qt, and thus building Desktop App and not FW
#ifdef QT_GUI_LIB
#include <QMultiMap>
#endif

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
#elif defined(ARDUINO_ESP32_S3_WROOM1_DevKitC_1_N16R8)  // 696969 aggiunto da me
    #define OPENFIRE_BOARD "esp32-s3-devkitc-1"
#else
    #define OPENFIRE_BOARD "generic"
#endif // board

class OF_Const
{
public:
    // !!! These orders should remain the same to maintain backwards compatibility!!!
    // Any new slots should explicitly be added at the bottom above the "count" line

    // Inputs
    enum {
        unavailable = -2,
        btnUnmapped = -1,
        btnTrigger = 0,
        btnGunA,
        btnGunB,
        btnGunC,
        btnStart,
        btnSelect,
        btnGunUp,
        btnGunDown,
        btnGunLeft,
        btnGunRight,
        btnPedal,
        btnPedal2,
        btnHome,
        btnPump,
        rumblePin,
        solenoidPin,
        rumbleSwitch,
        solenoidSwitch,
        autofireSwitch,
        neoPixel,
        ledR,
        ledG,
        ledB,
        camSDA,
        camSCL,
        periphSDA,
        periphSCL,
        battery,
        analogX,
        analogY,
        tempPin,
        boardInputsCount
    } boardInputs_e;

    // Boolean/toggle settings
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
        boolTypesCount
    } boolTypes_e;

    // Variable settings
    enum {
        rumbleStrength = 0,
        rumbleInterval,
        solenoidNormalInterval,
        solenoidFastInterval,
        solenoidHoldLength,
        autofireWaitFactor,
        holdToPauseLength,
        customLEDcount,
        customLEDstatic,
        customLEDcolor1,
        customLEDcolor2,
        customLEDcolor3,
        settingsTypesCount
    } settingsTypes_e;

    // Layout types
    enum {
        layoutSquare = 0,
        layoutDiamond
    } layoutTypes_e;

    typedef struct {
        int8_t pin[50]; // 30]; modificato da me 696969 per gestire ESP32S3 necessitano 48 pin per rp2040 ne bastavano 30 .. l'esp32 espone 48 pin da 1 in avanti .. arrotondato a 50
    } boardMap_t;

    /// @brief      Map of default pin mappings for each supported board
    /// @details    Key = board, int array maps to RP2040 GPIO where each value is a FW function (or unmapped).
    inline static const std::unordered_map<std::string, boardMap_t> boardsPresetsMap = {
        //=====================================================================================================
        // Notes: rpi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipico",             {btnGunA,       btnGunB,        btnGunC,        btnStart,       btnSelect,
                                 btnHome,       btnGunUp,       btnGunDown,     btnGunLeft,     btnGunRight,
                                 ledR,          ledG,           ledB,           btnPump,        btnPedal,
                                 btnTrigger,    solenoidPin,    rumblePin,      btnUnmapped,    btnUnmapped,
                                 camSDA,        camSCL,         btnUnmapped,    unavailable,    unavailable,
                                 unavailable,   btnUnmapped,    btnUnmapped,    tempPin,        unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable}},
        //=====================================================================================================
        // Notes: rpi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipicow",            {btnGunA,       btnGunB,        btnGunC,        btnStart,       btnSelect,
                                 btnHome,       btnGunUp,       btnGunDown,     btnGunLeft,     btnGunRight,
                                 ledR,          ledG,           ledB,           btnPump,        btnPedal,
                                 btnTrigger,    solenoidPin,    rumblePin,      periphSDA,      periphSCL,
                                 camSDA,        camSCL,         btnUnmapped,    unavailable,    unavailable,
                                 unavailable,   analogY,        analogX,        tempPin,        unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable}},
        //=====================================================================================================
        // Notes: pins 13-17 & 21-23 are unexposed
        {"adafruitItsyRP2040",  {btnUnmapped,   btnUnmapped,    camSDA,         camSCL,         btnPedal,
                                 btnUnmapped,   btnTrigger,     btnGunDown,     btnGunLeft,     btnGunUp,
                                 btnGunRight,   btnGunC,        btnUnmapped,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   unavailable,    unavailable,    unavailable,    rumblePin,
                                 solenoidPin,   btnGunB,        btnGunA,        btnStart,       btnSelect,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable}},
        //=====================================================================================================
        // Notes: pins 11-17 & 21-25 are unexposed
        {"adafruitKB2040",      {btnUnmapped,   btnUnmapped,    camSDA,         camSCL,         btnGunB,
                                 rumblePin,     btnGunC,        solenoidPin,    btnSelect,      btnStart,
                                 btnGunRight,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    btnGunUp,       btnGunLeft,
                                 btnGunDown,    unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   tempPin,        btnHome,        btnTrigger,     btnGunA,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable}},
        //=====================================================================================================
        // Notes: pins 2-3, 8-11, 14, & 22-24 are unexposed;
        //        some other pins are analog, but controlled by NINA and thus are unavailable in OpenFIRE for the moment
        {"arduinoNanoRP2040",   {btnTrigger,    btnPedal,       unavailable,    unavailable,    btnGunA,
                                 btnGunC,       btnUnmapped,    btnGunB,        unavailable,    unavailable,
                                 unavailable,   unavailable,    camSDA,         camSCL,         unavailable,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    unavailable,    unavailable,    unavailable,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    tempPin,        btnUnmapped,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable}},
        //=====================================================================================================
        // Note: pin 16 is reserved for the board's builtin NeoPixel (not currently used?);
        //       pins 17-25 are underside pads which are not exposed in the app for layout reasons;
        {"waveshareZero",       {btnTrigger,    btnGunA,        btnGunB,        btnGunC,        btnStart,
                                 btnSelect,     btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    camSDA,
                                 camSCL,        unavailable,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    tempPin,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable}},
        //=====================================================================================================
        // Notes: ESP32 // 696969 aggiunto da me  /*xx*/ indica il numero del GPIO es. /*02*/ per GPIO2
        {"esp32-s3-devkitc-1",  {/*00*/unavailable,   /*01*/btnTrigger,     /*02*/btnGunRight,    /*03*/btnUnmapped,    /*04*/analogX,
                                 /*05*/analogY,       /*06*/tempPin,        /*07*/btnUnmapped,    /*08*/camSDA,         /*09*/camSCL,
                                 /*10*/btnUnmapped,   /*11*/btnUnmapped,    /*12*/btnUnmapped,    /*13*/btnUnmapped,    /*14*/btnUnmapped,
                                 /*15*/periphSCL,     /*16*/rumblePin,      /*17*/solenoidPin,    /*18*/periphSDA,      /*19*/unavailable,
                                 /*20*/unavailable,   /*21*/btnGunC,        /*22*/unavailable,    /*23*/unavailable,    /*24*/unavailable,
                                 /*25*/unavailable,   /*26*/unavailable,    /*27*/unavailable,    /*28*/unavailable,    /*29*/unavailable,
                                 /*30*/unavailable,   /*31*/unavailable,    /*32*/unavailable,    /*33*/unavailable,    /*34*/unavailable,
                                 /*35*/btnHome,       /*36*/btnGunA,        /*37*/btnGunB,        /*38*/btnSelect,      /*39*/btnStart,
                                 /*40*/btnGunUp,      /*41*/btnGunDown,     /*42*/btnGunLeft,     /*43*/unavailable,    /*44*/unavailable,
                                 /*45*/btnPump,       /*46*/unavailable,    /*47*/btnPedal,       /*48*/neoPixel,       /*49*/unavailable}},
        //=====================================================================================================
    };

// Only needed for the Desktop App, don't build for microcontroller firmware!
#ifdef QT_VERSION

    // Used by pinBoxes, matching boardInputs_e (except "unavailable")
    inline static const QStringList valuesNameList = {
        "Unmapped",
        "Trigger",
        "Button A",
        "Button B",
        "Button C",
        "Start",
        "Select",
        "D-Pad Up",
        "D-Pad Down",
        "D-Pad Left",
        "D-Pad Right",
        "Pedal",
        "Alt Pedal",
        "Home Button",
        "Pump Action",
        "Rumble Signal",
        "Solenoid Signal",
        "Rumble Switch",
        "Solenoid Switch",
        "Autofire Switch",
        "External NeoPixel",
        "RGB LED Red",
        "RGB LED Green",
        "RGB LED Blue",
        "Camera SDA",
        "Camera SCL",
        "Peripherals SDA",
        "Peripherals SCL",
        "Battery Sensor",
        "Analog Pin X",
        "Analog Pin Y",
        "Temp Sensor"
    };

    inline static const QMap<std::string, const char *> boardNames = {
        {"rpipico",             "Raspberry Pi Pico (RP2040)"},
        {"rpipicow",            "Raspberry Pi Pico W (RP2040)"},
        {"adafruitItsyRP2040",  "Adafruit ItsyBitsy RP2040"},
        {"adafruitKB2040",      "Adafruit Keeboar KB2040"},
        {"arduinoNanoRP2040",   "Arduino Nano Connect RP2040"},
        {"waveshareZero",       "Waveshare Zero RP2040"},
        {"esp32-s3-devkitc-1",  "Esp32-S3 Devkitc-1"},                   // 696969 aggiunto da me
        // Add more here!
        {"generic",             "Unknown Board"}
    };

    enum {
        posNothing = 0,
        posLeft = 32,
        posRight = 64,
        posMiddle = 128
    } boardBoxPositions_e;

    typedef struct {
        uint8_t pin[50]; // 30]; modificato da me 696969 per gestire ESP32S3 necessitano 48 pin per rp2040 ne bastavano 30 .. l'esp32 espone 48 pin da 1 in avanti .. arrotondato a 50
    } boardBoxPosMap_t;

    /// @brief      Map of graphical placement for each pin in the application
    /// @details    Key = board, int array maps to RP2040 GPIO.
    ///             Each pin should be a combination of grid layout slot it should be in,
    ///             added by the grid it should belong to.
    ///             Unexposed pins should use only posNothing (0).
    ///             (Yep, bitpacking! Three most significant bits determine left/right/under position)
    ///             (If anyone is aware of a better way of doing this, please let me know/send a PR!)
    inline static const QMap<std::string, boardBoxPosMap_t> boardsBoxPositions = {
        //=====================================================================================================
        // Raspberry Pi Pico: 15 pins left, rest of the pins right. Mostly linear order save for the reserved pins.
        // Notes: rpi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipico",             {1+posLeft,     2+posLeft,      4+posLeft,      5+posLeft,      6+posLeft,
                                 7+posLeft,     9+posLeft,      10+posLeft,     11+posLeft,     12+posLeft,
                                 14+posLeft,    15+posLeft,     16+posLeft,     17+posLeft,     19+posLeft,
                                 20+posLeft,    20+posRight,    19+posRight,    17+posRight,    16+posRight,
                                 15+posRight,   14+posRight,    12+posRight,    posNothing,     posNothing,
                                 posNothing,    10+posRight,    9+posRight,     7+posRight,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing}},
        //=====================================================================================================
        // Raspberry Pi Pico W: same as non-W Pico.
        // Notes: rpi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipicow",            {1+posLeft,     2+posLeft,      4+posLeft,      5+posLeft,      6+posLeft,
                                 7+posLeft,     9+posLeft,      10+posLeft,     11+posLeft,     12+posLeft,
                                 14+posLeft,    15+posLeft,     16+posLeft,     17+posLeft,     19+posLeft,
                                 20+posLeft,    20+posRight,    19+posRight,    17+posRight,    16+posRight,
                                 15+posRight,   14+posRight,    12+posRight,    posNothing,     posNothing,
                                 posNothing,    10+posRight,    9+posRight,     7+posRight,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing}},
        //=====================================================================================================
        // Adafruit ItsyBitsy RP2040: A very cluttered and kind of unfriendly layout tbh :(
        // Notes: pins 13-17 & 21-23 are unexposed
        {"adafruitItsyRP2040",  {13+posRight,   14+posRight,    12+posRight,    11+posRight,    2+posMiddle,
                                 1+posMiddle,   9+posRight,     8+posRight,     7+posRight,     6+posRight,
                                 5+posRight,    4+posRight,     14+posLeft,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     11+posLeft,     12+posLeft,
                                 13+posLeft,    posNothing,     posNothing,     posNothing,     9+posLeft,
                                 10+posLeft,    5+posLeft,      6+posLeft,      7+posLeft,      8+posLeft,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing}},
        //=====================================================================================================
        // Adafruit KB2040: Like the Itsy with more padding.
        // Notes: pins 11-17 & 21-25 are unexposed
        {"adafruitKB2040",      {3+posLeft,     4+posLeft,      7+posLeft,      8+posLeft,      9+posLeft,
                                 10+posLeft,    11+posLeft,     12+posLeft,     13+posLeft,     14+posLeft,
                                 14+posRight,   posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     11+posRight,    13+posRight,
                                 12+posRight,   posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    10+posRight,    9+posRight,     8+posRight,     7+posRight,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing}},
        //=====================================================================================================
        // Arduino Nano RP2040: gweh
        // Notes: pins 2-3, 8-11, 14, & 22-24 are unexposed;
        //        some other pins are analog, but controlled by NINA and thus are unavailable in OpenFIRE for the moment
        {"arduinoNanoRP2040",   {18+posRight,   17+posRight,    posNothing,     posNothing,     4+posRight,
                                 6+posRight,    4+posLeft,      5+posRight,     posNothing,     posNothing,
                                 posNothing,    posNothing,     11+posLeft,     12+posLeft,     posNothing,
                                 13+posRight,   12+posRight,    11+posRight,    10+posRight,    9+posRight,
                                 8+posRight,    7+posRight,     posNothing,     posNothing,     posNothing,
                                 14+posRight,   7+posLeft,      8+posLeft,      9+posLeft,      10+posLeft,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing}},
        //=====================================================================================================
        // Waveshare Zero RP2040: Clockwise layout
        // Note: pin 16 is reserved for the board's builtin NeoPixel (not currently used?);
        //       pins 17-25 are underside pads which are not exposed in the app for layout reasons;
        {"waveshareZero",       {2+posRight,    3+posRight,     4+posRight,     5+posRight,     6+posRight,
                                 7+posRight,    8+posRight,     9+posRight,     10+posRight,    11+posRight,
                                 3+posMiddle,   2+posMiddle,    11+posLeft,     10+posLeft,     9+posLeft,
                                 8+posLeft,     posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    7+posLeft,      6+posLeft,      5+posLeft,      4+posLeft,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing}},
        //=====================================================================================================
        // Insert new layouts below this one!
        // Feel free to use any of the above as a template.
        // ***

         //=====================================================================================================
        // Esp32 S3 Devkitc 1
        // Notes: 696969 Esp32   /*xx*/ indica il numero del GPIO es. /*02*/ per GPIO2
        {"esp32-s3-devkitc-1",  {/*00*/posNothing,    /*01*/19+posLeft,     /*02*/18+posLeft,     /*03*/posNothing,     /*04*/19+posRight,
                                 /*05*/18+posRight,   /*06*/17+posRight,    /*07*/posNothing,     /*08*/11+posRight,    /*09*/8+posRight,
                                 /*10*/posNothing,    /*11*/posNothing,     /*12*/posNothing,     /*13*/posNothing,     /*14*/posNothing,
                                 /*15*/15+posRight,   /*16*/14+posRight,    /*17*/13+posRight,    /*18*/12+posRight,    /*19*/posNothing,
                                 /*20*/posNothing,    /*21*/5+posLeft,      /*22*/posNothing,     /*23*/posNothing,     /*24*/posNothing,
                                 /*25*/posNothing,    /*26*/posNothing,     /*27*/posNothing,     /*28*/posNothing,     /*29*/posNothing,
                                 /*30*/posNothing,    /*31*/posNothing,     /*32*/posNothing,     /*33*/posNothing,     /*34*/posNothing,
                                 /*35*/10+posLeft,    /*36*/11+posLeft,     /*37*/12+posLeft,     /*38*/13+posLeft,     /*39*/14+posLeft,
                                 /*40*/15+posLeft,    /*41*/16+posLeft,     /*42*/17+posLeft,     /*43*/posNothing,     /*44*/posNothing,
                                 /*45*/8+posLeft,     /*46*/posNothing,     /*47*/6+posLeft,      /*48*/7+posLeft,      /*49*/posNothing}},
        //=====================================================================================================

        //=====================================================================================================
        // Generic layout
        // Just reveal all pins; user assumes full responsibility if something goes wrong here
        {"generic",             {1+posLeft,     2+posLeft,      3+posLeft,      4+posLeft,      5+posLeft,
                                 6+posLeft,     7+posLeft,      8+posLeft,      9+posLeft,      10+posLeft,
                                 11+posLeft,    12+posLeft,     13+posLeft,     14+posLeft,     15+posLeft,
                                 16+posLeft,    16+posRight,    15+posRight,    14+posRight,    13+posRight,
                                 12+posRight,   11+posRight,    10+posRight,    9+posRight,     8+posRight,
                                 7+posRight,    6+posRight,     5+posRight,     4+posRight,     3+posRight,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing,
                                 posNothing,    posNothing,     posNothing,     posNothing,     posNothing}}
    };

    typedef struct {
        const char * name;
        int8_t pin[50];  // 30]; modificato da me 696969 per gestire ESP32S3 necessitano 48 pin per rp2040 ne bastavano 30 .. l'esp32 espone 48 pin da 1 in avanti .. arrotondato a 50
    } boardAltPresetsMap_t;

    /// @brief      MultiMap of alternative pin mappings for supported boards to show in the application.
    /// @details    Key = board (one board can be multiple), string literal label, int array maps to RP2040 GPIO where each value is a FW function (or unmapped).
    inline static const QMultiMap<std::string, boardAltPresetsMap_t> boardsAltPresets = {
        //=====================================================================================================
        // Raspberry Pi Pico Presets (currently a test)
        // Notes: rpi boards do not expose pins 23-25; pin 29/A3 is used for builtin chipset temp monitor
        {"rpipico",             {"Test",
                                 btnPump,       btnPedal,       btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    unavailable,    unavailable,
                                 unavailable,   btnUnmapped,    btnUnmapped,    btnUnmapped,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable}},

        {"rpipico",             {"Test 2",
                                 btnGunA,       btnTrigger,     btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnUnmapped,    unavailable,    unavailable,
                                 unavailable,   btnUnmapped,    btnUnmapped,    btnUnmapped,    unavailable,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable}},


        //=====================================================================================================
        // Adafruit ItsyBitsy RP2040 Presets
        // Notes: pins 13-17 & 21-23 are unexposed
        {"adafruitItsyRP2040",  {"SAMCO 2.0",
                                 btnUnmapped,   btnUnmapped,    camSDA,         camSCL,         btnPedal,
                                 btnUnmapped,   btnTrigger,     btnGunDown,     btnGunLeft,     btnGunUp,
                                 btnGunRight,   btnHome,        btnUnmapped,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   unavailable,    unavailable,    unavailable,    rumblePin,
                                 solenoidPin,   btnGunB,        btnGunA,        btnStart,       btnSelect,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable}},


        {"adafruitItsyRP2040",  {"SAMCO 1.1",
                                 btnUnmapped,   btnUnmapped,    camSDA,         camSCL,         btnUnmapped,
                                 btnUnmapped,   btnGunA,        btnGunB,        rumblePin,      btnHome,
                                 btnTrigger,    btnUnmapped,    btnUnmapped,    unavailable,    unavailable,
                                 unavailable,   unavailable,    unavailable,    btnUnmapped,    btnUnmapped,
                                 btnUnmapped,   unavailable,    unavailable,    unavailable,    btnUnmapped,
                                 btnUnmapped,   btnUnmapped,    btnPedal,       btnUnmapped,    btnUnmapped,
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable, 
                                 unavailable,   unavailable,    unavailable,    unavailable,    unavailable}},

    };

#endif
};


#endif // _OPENFIRESHARED_H_
