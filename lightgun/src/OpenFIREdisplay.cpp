/*!
 * @file OpenFIREdisplay.cpp
 * @brief Macros for lightgun HUD display (primarily for SSD1306 OLED modules).
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
 * @copyright That One Seong, 2024
 * @copyright GNU Lesser General Public License
 */

// we're using our own splash screen kthx ada
#define SSD1306_NO_SPLASH

#include <Arduino.h>
#include <string.h>

#ifdef USE_LOVYAN_GFX
  // nulla 
#else
  // #include <Adafruit_GFX.h>
#endif

#include <Wire.h>
#include <TinyUSB_Devices.h>

#include "OpenFIREdisplay.h"
#include "OpenFIREconstant.h"
#include "OpenFIREprefs.h"
#include "OpenFIREFeedback.h"
#include "OpenFIREDefines.h"

#include "OpenFIREweb.h"

bool ExtDisplay::Begin()
{
    if(display != nullptr) {
        delete display;
        display = nullptr;
    }

    // TODO: for some reason, doing this AFTER saving updated pins settings (even when doing it from defaults and there's no default mappings for peripheral pins)
    // causes the board to hang. Even though this is all correct (and any display objects should get deleted from the above, so don't think it can be a new object thing)...
    if(OF_Prefs::pins[OF_Const::periphSCL] >= 0 && OF_Prefs::pins[OF_Const::periphSDA] >= 0) {
      #ifdef ARDUINO_ARCH_ESP32
        Wire1.setPins(OF_Prefs::pins[OF_Const::periphSDA], OF_Prefs::pins[OF_Const::periphSCL]);  // [ESP32_PORT] per esp32
        #ifdef USE_LOVYAN_GFX
        display = new LGFX_SSD1306(1/*i2c_port wire_1*/,OF_Prefs::pins[OF_Const::periphSDA], OF_Prefs::pins[OF_Const::periphSCL], OF_Prefs::toggles[OF_Const::i2cOLEDaltAddr] ? 0x3D : 0x3C, SCREEN_WIDTH, SCREEN_HEIGHT);
        #else
        display = new Adafruit_SSD1306(SCREEN_WIDTH, SCREEN_HEIGHT, &Wire1, -1);
        #endif
      #else // rp2040
        if(bitRead(OF_Prefs::pins[OF_Const::periphSCL], 1) && bitRead(OF_Prefs::pins[OF_Const::periphSDA], 1)) {
            // I2C1
            if(bitRead(OF_Prefs::pins[OF_Const::periphSCL], 0) && !bitRead(OF_Prefs::pins[OF_Const::periphSDA], 0)) {
                Wire1.end();
                // SDA/SCL are indeed on verified correct pins
                Wire1.setSDA(OF_Prefs::pins[OF_Const::periphSDA]);
                Wire1.setSCL(OF_Prefs::pins[OF_Const::periphSCL]);
                #ifdef USE_LOVYAN_GFX
                display = new LGFX_SSD1306(1/*i2c_port wire_1*/,OF_Prefs::pins[OF_Const::periphSDA], OF_Prefs::pins[OF_Const::periphSCL], OF_Prefs::toggles[OF_Const::i2cOLEDaltAddr] ? 0x3D : 0x3C, SCREEN_WIDTH, SCREEN_HEIGHT);
                #else
                //Wire1.end();
                // SDA/SCL are indeed on verified correct pins
                //Wire1.setSDA(OF_Prefs::pins[OF_Const::periphSDA]);
                //Wire1.setSCL(OF_Prefs::pins[OF_Const::periphSCL]);
                display = new Adafruit_SSD1306(SCREEN_WIDTH, SCREEN_HEIGHT, &Wire1, -1);
                #endif
            } else return false;
        } else if(!bitRead(OF_Prefs::pins[OF_Const::periphSCL], 1) && !bitRead(OF_Prefs::pins[OF_Const::periphSDA], 1)) {
            // I2C0
            if(bitRead(OF_Prefs::pins[OF_Const::periphSCL], 0) && !bitRead(OF_Prefs::pins[OF_Const::periphSDA], 0)) {
                Wire.end();
                // SDA/SCL are indeed on verified correct pins
                Wire.setSDA(OF_Prefs::pins[OF_Const::periphSDA]);
                Wire.setSCL(OF_Prefs::pins[OF_Const::periphSCL]);
                #ifdef USE_LOVYAN_GFX
                display = new LGFX_SSD1306(0/*i2c_port wire*/,OF_Prefs::pins[OF_Const::periphSDA], OF_Prefs::pins[OF_Const::periphSCL], OF_Prefs::toggles[OF_Const::i2cOLEDaltAddr] ? 0x3D : 0x3C, SCREEN_WIDTH, SCREEN_HEIGHT);
                #else
                //Wire.end();
                // SDA/SCL are indeed on verified correct pins
                //Wire.setSDA(OF_Prefs::pins[OF_Const::periphSDA]);
                //Wire.setSCL(OF_Prefs::pins[OF_Const::periphSCL]);
                display = new Adafruit_SSD1306(SCREEN_WIDTH, SCREEN_HEIGHT, &Wire, -1);
                #endif
            } else return false;
        } else return false;
      #endif
    } else return false;

    #ifdef USE_LOVYAN_GFX
    if(display->init()) {
    #else
    if(display->begin(SSD1306_SWITCHCAPVCC, OF_Prefs::toggles[OF_Const::i2cOLEDaltAddr] ? 0x3D : 0x3C)) {
    #endif    
        display->clearDisplay();
        ScreenModeChange(Screen_None);
        return true;
    } else return false;
}

void ExtDisplay::Stop()
{
    irOverlayValid = 0;
    if(display != nullptr) {
        delete display;
        display = nullptr;
    }
}

// Icona 16x16 "Impostazioni Windows" (Buco centrale 4x4 smussato)
static const unsigned char PROGMEM webGearIco[] = {
  0x01, 0x80, 0x01, 0x80, 0x33, 0xcc, 0x3f, 0xfc, 
  0x1f, 0xf8, 0x1f, 0xf8, 0x7e, 0x7e, 0x7c, 0x3e, 
  0x7c, 0x3e, 0x7e, 0x7e, 0x1f, 0xf8, 0x1f, 0xf8, 
  0x3f, 0xfc, 0x33, 0xcc, 0x01, 0x80, 0x01, 0x80
};

void ExtDisplay::TopPanelUpdate(const char *textPrefix, const char *profText)
{
    if(display != nullptr) {
        bool isWeb = false;
#ifdef ARDUINO_ARCH_ESP32
        if (OF_WebConfigModeActive) isWeb = true;
#endif

        // 1. Sfondo e linea della barra (bianca se WebApp, altrimenti nera)
        display->fillRect(0, 0, 128, 16, isWeb ? WHITE : BLACK);
        display->drawFastHLine(0, 15, 128, isWeb ? BLACK : WHITE);
        
        display->setCursor(2, 2);
        display->setTextSize(1);
        
        // 2. Colore del testo (nero se WebApp, altrimenti bianco)
        display->setTextColor(isWeb ? BLACK : WHITE, isWeb ? WHITE : BLACK);
        
        display->print(textPrefix);
        if(profText != nullptr)
            display->print(profText); 
            
        // 3. Il "Badge" WEB a destra con l'icona
        if (isWeb) {
            // Disegna il quadratino nero all'estrema destra (largo 20 pixel, alto 16)
            display->fillRect(106, 0, 22, 16, BLACK); 
            
            // Disegna l'icona del mappamondo in BIANCO sopra allo sfondo nero
            // (La posizioniamo a X=110, Y=0, grandezza 16x16)
            display->drawBitmap(110, 0, webGearIco, 16, 16, WHITE);
        }
        
        display->display();
    }
}

void ExtDisplay::ScreenModeChange(const int &screenMode, const bool &isAnalog)
{
    // The new screen replaces the background; do not restore pixels from the old one.
    irOverlayValid = 0;
    if(display != nullptr) {
        idleTimeStamp = millis();
        display->fillRect(0, 16, 128, 48, BLACK);
        if(screenState >= Screen_Mamehook_Single &&
           screenMode == Screen_Normal) {
            currentAmmo = 0, currentLife = 0;
        }
        screenState = screenMode;
        display->setTextColor(WHITE, BLACK);
        switch(screenMode) {
          case Screen_Normal:
            if(mister) display->drawBitmap(48, 23, misterIco, MISTERKUN_WIDTH, MISTERKUN_HEIGHT, WHITE);
            else {
                #ifdef ARDUINO_ARCH_ESP32
                if(TinyUSBDevices.onBattery) { display->drawBitmap(2, 46, wifiConnectIco, CONNECTION_WIDTH, CONNECTION_HEIGHT, WHITE); }
                #else //rp2040
                if(TinyUSBDevices.onBattery) { display->drawBitmap(2, 46, btConnectIco, CONNECTION_WIDTH, CONNECTION_HEIGHT, WHITE); }
                #endif
                else { display->drawBitmap(2, 46, usbConnectIco, CONNECTION_WIDTH, CONNECTION_HEIGHT, WHITE); }
                if(isAnalog) { display->drawBitmap(108, 49, gamepadIco, GAMEPAD_WIDTH, GAMEPAD_HEIGHT, WHITE); }
                else { display->drawBitmap(109, 48, mouseIco, MOUSE_WIDTH, MOUSE_HEIGHT, WHITE); }
            }
            break;
          case Screen_None:
          case Screen_Docked:
            display->fillRect(0, 0, 128, 16, BLACK);
            display->drawBitmap(24, 0, customSplashBanner, CUSTSPLASHBANN_WIDTH, CUSTSPLASHBANN_HEIGHT, WHITE);
            display->drawBitmap(40, 16, customSplash, CUSTSPLASH_WIDTH, CUSTSPLASH_HEIGHT, WHITE);
            display->display();
            break;
          case Screen_Init:
            display->setTextSize(2);
            display->setCursor(20, 18);
            display->println("Welcome!");
            display->setTextSize(1);
            display->setCursor(12, 40);
            display->println(" Pull trigger to");
            display->setCursor(12, 52);
            display->println("start calibration!");
            break;
          case Screen_IRTest:
            TopPanelUpdate("", "IR Test");
            break;
          case Screen_Saving:
            TopPanelUpdate("", "Saving Profiles");
            display->setTextSize(2);
            display->setCursor(16, 18);
            display->println("Saving...");
            break;
          case Screen_SaveSuccess:
            display->setTextSize(2);
            display->setCursor(30, 18);
            display->println("Save");
            display->setCursor(4, 40);
            display->println("successful");
            break;
          case Screen_SaveError:
            display->setTextSize(2);
            display->setCursor(30, 18);
            display->setTextColor(BLACK, WHITE);
            display->println("Save");
            display->setCursor(22, 40);
            display->println("failed");
            break;
          case Screen_Mamehook_Single:
            #ifdef ARDUINO_ARCH_ESP32  
            if(TinyUSBDevices.onBattery) { display->drawBitmap(2, 46, wifiConnectIco, CONNECTION_WIDTH, CONNECTION_HEIGHT, WHITE); }
            #else //rp2040
            if(TinyUSBDevices.onBattery) { display->drawBitmap(2, 46, btConnectIco, CONNECTION_WIDTH, CONNECTION_HEIGHT, WHITE); }
            #endif
            else { display->drawBitmap(2, 46, usbConnectIco, CONNECTION_WIDTH, CONNECTION_HEIGHT, WHITE); }
            if(isAnalog) { display->drawBitmap(108, 49, gamepadIco, GAMEPAD_WIDTH, GAMEPAD_HEIGHT, WHITE); }
            else { display->drawBitmap(109, 48, mouseIco, MOUSE_WIDTH, MOUSE_HEIGHT, WHITE); }
            if(serialDisplayType == ScreenSerial_Life && lifeBar) {
              display->drawBitmap(52, 23, lifeBarBanner, LIFEBAR_BANNER_WIDTH, LIFEBAR_BANNER_HEIGHT, WHITE);
              display->drawBitmap(11, 35, lifeBarLarge, LIFEBAR_LARGE_WIDTH, LIFEBAR_LARGE_HEIGHT, WHITE);
              PrintLife(currentLife);
            } else if(serialDisplayType == ScreenSerial_Ammo) {
              PrintAmmo(currentAmmo);
            } else {
              display->drawBitmap(0, 17, mamehookIco, MAMEHOOK_WIDTH, MAMEHOOK_HEIGHT, WHITE);
            }
            break;
          case Screen_Mamehook_Dual:
            display->drawBitmap(63, 16, dividerLine, DIVIDER_WIDTH, DIVIDER_HEIGHT, WHITE);
            if(lifeBar) {
              display->drawBitmap(20, 23, lifeBarBanner, LIFEBAR_BANNER_WIDTH, LIFEBAR_BANNER_HEIGHT, WHITE);
              display->drawBitmap(2, 37, lifeBarSmall, LIFEBAR_SMALL_WIDTH, LIFEBAR_SMALL_HEIGHT, WHITE);
            }
            PrintAmmo(currentAmmo);
            PrintLife(currentLife);
            break;
        }
        display->display();
    }
}

void ExtDisplay::IdleOps()
{
    if(display != nullptr) {
        switch(screenState) {
        case Screen_Normal:
        case Screen_Mamehook_Single:
        case Screen_Mamehook_Dual:
          #ifdef USES_TEMP
          if(OF_Prefs::pins[OF_Const::tempPin] > -1) {
              if(millis() - idleTimeStamp > OLED_IDLEUPD_INTERVAL) {
                  idleTimeStamp = millis();
                  if(showingTemp) {
                      TopPanelUpdate("Prof: ", OF_Prefs::profiles[OF_Prefs::currentProfile].name);
                      showingTemp = false;
                  } else {
                      idleTempStamp = idleTimeStamp;
                      ShowTemp();
                      showingTemp = true;
                  }
              } else if(showingTemp) {
                  if(millis() - idleTempStamp > OLED_TEMPUPD_INTERVAL) {
                      idleTempStamp = millis();
                      if(currentTemp != OF_FFB::temperatureCurrent) {
                          ShowTemp();
                      }
                  }
              }
          }
          #endif // USES_TEMP
          break;
        case Screen_Pause:
        case Screen_Profile:
        case Screen_Saving:
        case Screen_Calibrating:
        default:
          break;
        }
    }
}

#ifdef USES_TEMP
void ExtDisplay::ShowTemp()
{
    if (OF_FFB::temperatureCurrent == OF_Const::TEMPERATURE_SENSOR_ERROR_VALUE) {
      // Analog read maxed out, likely disconnected/temp sensor fault
      TopPanelUpdate("Temp: ", "Fault!");
      return;
    } else if(OF_FFB::temperatureCurrent < 10) {
        tempString[0] = OF_FFB::temperatureCurrent + '0';
        tempString[1] = ' ';
        tempString[2] = 'C';
        tempString[3] = '\0';
    } else {
        int tempLeft = OF_FFB::temperatureCurrent / 10;
        int tempRight = OF_FFB::temperatureCurrent - (tempLeft * 10);
        tempString[0] = tempLeft + '0';
        tempString[1] = tempRight + '0';
        tempString[2] = ' ';
        tempString[3] = 'C';
        tempString[4] = '\0';
    }

    TopPanelUpdate("Temp: ", tempString);
    currentTemp = OF_FFB::temperatureCurrent;
}
#endif // USES_TEMP

void ExtDisplay::CalibrationStageUpdate(uint8_t stage)
{
    if(display == nullptr || stage > FW_Const::Cali_Verify)
        return;

    // Static drawing on stage changes only, never at the camera frame rate.
    // Leave the profile header (rows 0..15) untouched.
    irOverlayValid = 0;
    screenState = Screen_Calibrating;
    display->fillRect(0, 16, SCREEN_WIDTH, SCREEN_HEIGHT - 16, BLACK);
    display->setTextSize(1);
    display->setTextColor(WHITE, BLACK);

    const char *lines[3];
    char progress[] = {char('1' + stage), '/', '6', '\0'};
    int textX, textWidth, textY, lineSpacing;
    if(stage == FW_Const::Cali_Verify) {
        lines[0] = "CHECK AIM";
        lines[1] = "TRIGGER: CONFIRM";
        lines[2] = "A/B: REPEAT";
        textX = 0;
        textWidth = SCREEN_WIDTH;
        textY = 21;
        lineSpacing = 15;
    } else {
        static const char *const directions[] = {
            "CENTER", "TOP", "BOTTOM", "LEFT", "RIGHT", "CENTER"
        };
        int x = 33, y = 40;
        switch(stage) {
            case FW_Const::Cali_Top:    y = 24; break;
            case FW_Const::Cali_Bottom: y = 56; break;
            case FW_Const::Cali_Left:   x = 5;  break;
            case FW_Const::Cali_Right:  x = 61; break;
            default: break;
        }

        display->drawRect(5, 24, 57, 33, WHITE);
        // A small stationary reticle, including at the edges of the TV.
        display->fillRect(x - 5, y - 5, 11, 11, BLACK);
        display->drawCircle(x, y, 4, WHITE);
        display->drawFastHLine(x - 5, y, 3, WHITE);
        display->drawFastHLine(x + 3, y, 3, WHITE);
        display->drawFastVLine(x, y - 5, 3, WHITE);
        display->drawFastVLine(x, y + 3, 3, WHITE);
        display->drawFastHLine(x - 1, y, 3, WHITE);
        display->drawFastVLine(x, y - 1, 3, WHITE);

        lines[0] = progress;
        lines[1] = directions[stage];
        lines[2] = "SHOOT";
        textX = 74;
        textWidth = 48;
        textY = 20;
        lineSpacing = 13;
    }

    for(uint8_t i = 0; i < 3; ++i) {
        // Both supported backends use the default 6x8 font at size 1.
        const int width = int(strlen(lines[i])) * 6;
        display->setCursor(textX + (textWidth - width) / 2, textY + i * lineSpacing);
        display->print(lines[i]);
    }
    display->display();
}

// Warning: SLOOOOW, should only be used where the mouse isn't being updated.
// Use at your own discression.
void ExtDisplay::DrawVisibleIR(int32_t *pointX, int32_t *pointY, IRDrawMode_e mode)
{
    if(mode == IRDraw_None || display == nullptr)
        return;

    if(mode == IRDraw_Overlay) {
        DrawVisibleIROverlay(pointX, pointY);
        return;
    }

    // Keep the original default mode, including in-place mapping and clipping.
    if(mode == IRDraw_Clear && screenState != Screen_Calibrating) {
        irOverlayValid = 0;
        display->fillRect(0, 16, 128, 48, BLACK);
        for(uint8_t i = 0; i < 4; ++i) {
          if(pointX[i] < 0 || pointX[i] > 1920 || pointY[i] < 0 || pointY[i] > 1080)
            continue;
          pointX[i] = map(pointX[i], 0, 1920, 0, 128);
          pointY[i] = map(pointY[i], 0, 1080, 16, 64);
          //pointY[i] = constrain(pointY[i], 16, 64);
          display->fillCircle(pointX[i], pointY[i], 1, WHITE);
        }
        display->display();
    }
}

void ExtDisplay::DrawVisibleIROverlay(const int32_t *pointX, const int32_t *pointY)
{
    static const int8_t offsetX[5] = {0, -1, 1, 0, 0};
    static const int8_t offsetY[5] = {0, 0, 0, -1, 1};
    uint8_t nextX[4] = {};
    uint8_t nextY[4] = {};
    uint8_t nextValid = 0;
    uint8_t moved = 0;

    for(uint8_t i = 0; i < 4; ++i) {
        if(pointX[i] < 0 || pointX[i] > 1920 || pointY[i] < 0 || pointY[i] > 1080)
            continue;
        nextX[i] = (uint8_t)map(pointX[i], 0, 1920, 0, 128);
        nextY[i] = (uint8_t)map(pointY[i], 0, 1080, 16, 64);
        nextValid |= (uint8_t)(1U << i);
        if(!(irOverlayValid & (1U << i)) ||
           nextX[i] != irOverlay[i].x || nextY[i] != irOverlay[i].y)
            moved = 1;
    }

    // Camera movement below one OLED pixel needs no pixel work or I2C refresh.
    if(!moved && nextValid == irOverlayValid)
        return;

    // Batch pixel writes: LovyanGFX OLED panels can auto-refresh at endWrite().
    // Adafruit's startWrite/endWrite are compatible no-ops for this buffered OLED.
    display->startWrite();

    // Restore ALL old points before sampling ANY new ones: footprints can overlap.
    for(uint8_t i = 0; i < 4; ++i) {
        if(!(irOverlayValid & (1U << i)))
            continue;
        for(uint8_t bit = 0; bit < 5; ++bit) {
            const int16_t x = (int16_t)irOverlay[i].x + offsetX[bit];
            const int16_t y = (int16_t)irOverlay[i].y + offsetY[bit];
            if(x >= 0 && x < SCREEN_WIDTH && y >= 16 && y < SCREEN_HEIGHT)
                display->drawPixel(x, y, (irOverlay[i].pixels & (1U << bit)) ? WHITE : BLACK);
        }
    }

    // Read only the RAM framebuffer, through the selected graphics backend.
    for(uint8_t i = 0; i < 4; ++i) {
        if(!(nextValid & (1U << i)))
            continue;
        irOverlay[i].x = nextX[i];
        irOverlay[i].y = nextY[i];
        irOverlay[i].pixels = 0;
        for(uint8_t bit = 0; bit < 5; ++bit) {
            const int16_t x = (int16_t)nextX[i] + offsetX[bit];
            const int16_t y = (int16_t)nextY[i] + offsetY[bit];
            if(x < 0 || x >= SCREEN_WIDTH || y < 16 || y >= SCREEN_HEIGHT)
                continue;
            #ifdef USE_LOVYAN_GFX
                const uint8_t lit = display->readPixel(x, y) != 0;
            #else
                const uint8_t lit = display->getPixel(x, y);
            #endif
            if(lit)
                irOverlay[i].pixels |= (uint8_t)(1U << bit);
        }
    }

    // Sample all backgrounds first, then draw the same radius-1 white points.
    for(uint8_t i = 0; i < 4; ++i) {
        if(!(nextValid & (1U << i)))
            continue;
        for(uint8_t bit = 0; bit < 5; ++bit) {
            const int16_t x = (int16_t)nextX[i] + offsetX[bit];
            const int16_t y = (int16_t)nextY[i] + offsetY[bit];
            if(x >= 0 && x < SCREEN_WIDTH && y >= 16 && y < SCREEN_HEIGHT)
                display->drawPixel(x, y, WHITE);
        }
    }
    irOverlayValid = nextValid;
    display->endWrite();
    display->display();
}

void ExtDisplay::PauseScreenShow(const int &currentProf, const char* name1, const char* name2, const char* name3, const char* name4)
{
    if(display != nullptr) {
        const char* namesList[] = { name1, name2, name3, name4 };
        TopPanelUpdate("Using ", namesList[currentProf]); // names are placeholder
        display->fillRect(0, 16, 128, 48, BLACK);
        display->setTextSize(1);
        display->setCursor(0, 17);
        display->print(" A > ");
        display->println(name1);
        display->setCursor(0, 17+11);
        display->print(" B > ");
        display->println(name2);
        display->setCursor(0, 17+(11*2));
        display->print("Str> ");
        display->println(name3);
        display->setCursor(0, 17+(11*3));
        display->print("Sel> ");
        display->println(name4);
        display->display();
    }
}

void ExtDisplay::PauseListUpdate(const int &selection)
{
    if(display != nullptr) {
        display->fillRect(0, 16, 128, 48, BLACK);
        display->drawBitmap(60, 18, upArrowGlyph, ARROW_WIDTH, ARROW_HEIGHT, WHITE);
        display->drawBitmap(60, 59, downArrowGlyph, ARROW_WIDTH, ARROW_HEIGHT, WHITE);
        display->setTextSize(1);
        // Seong Note: Yeah, some of these are pretty out-of-bounds-esque behavior,
        // but pause mode selection in actual use would prevent some of these extremes from happening.
        // Just covering our asses.
        switch(selection) {
          case ScreenPause_Calibrate:
            display->setTextColor(WHITE, BLACK);
            display->setCursor(0, 25);
            display->println(" Send Escape Keypress");
            display->setTextColor(BLACK, WHITE);
            display->setCursor(0, 36);
            display->println(" Calibrate ");
            display->setTextColor(WHITE, BLACK);
            display->setCursor(0, 47);
            display->println(" Profile Select ");
            break;
          case ScreenPause_ProfileSelect:
            display->setTextColor(WHITE, BLACK);
            display->setCursor(0, 25);
            display->println(" Calibrate ");
            display->setTextColor(BLACK, WHITE);
            display->setCursor(0, 36);
            display->println(" Profile Select ");
            display->setTextColor(WHITE, BLACK);
            display->setCursor(0, 47);
            display->println(" Save Gun Settings ");
            break;
          case ScreenPause_Save:
            display->setTextColor(WHITE, BLACK);
            display->setCursor(0, 25);
            display->println(" Profile Select ");
            display->setTextColor(BLACK, WHITE);
            display->setCursor(0, 36);
            display->println(" Save Gun Settings ");
            display->setTextColor(WHITE, BLACK);
            display->setCursor(0, 47);
            if(OF_Prefs::pins[OF_Const::rumblePin] >= 0 && OF_Prefs::pins[OF_Const::rumbleSwitch] == -1) {
              display->println(" Rumble Toggle ");
            } else if(OF_Prefs::pins[OF_Const::solenoidPin] >= 0 && OF_Prefs::pins[OF_Const::solenoidSwitch] == -1) {
              display->println(" Solenoid Toggle ");
            } else {
              display->println(" Send Escape Keypress");
            }
            break;
          case ScreenPause_Rumble:
            display->setTextColor(WHITE, BLACK);
            display->setCursor(0, 25);
            display->println(" Save Gun Settings ");
            display->setTextColor(BLACK, WHITE);
            display->setCursor(0, 36);
            if(OF_Prefs::pins[OF_Const::rumblePin] >= 0 && OF_Prefs::pins[OF_Const::rumbleSwitch] == -1) {
              display->println(" Rumble Toggle ");
              display->setTextColor(WHITE, BLACK);
              display->setCursor(0, 47);
              if(OF_Prefs::pins[OF_Const::solenoidPin] >= 0 && OF_Prefs::pins[OF_Const::solenoidSwitch] == -1) {
                display->println(" Solenoid Toggle ");
              } else {
                display->println(" Send Escape Keypress");
              }
            } else if(OF_Prefs::pins[OF_Const::solenoidPin] >= 0 && OF_Prefs::pins[OF_Const::solenoidSwitch] == -1) {
              display->println(" Solenoid Toggle ");
              display->setTextColor(WHITE, BLACK);
              display->setCursor(0, 47);
              display->println(" Send Escape Keypress");
            } else {
              display->println(" Send Escape Keypress");
              display->setTextColor(WHITE, BLACK);
              display->setCursor(0, 47);
              display->println(" Calibrate ");
            }
            break;
          case ScreenPause_Solenoid:
            display->setTextColor(WHITE, BLACK);
            display->setCursor(0, 25);
            if(OF_Prefs::pins[OF_Const::rumblePin] >= 0 && OF_Prefs::pins[OF_Const::rumbleSwitch] == -1) {
              display->println(" Rumble Toggle ");
              display->setTextColor(BLACK, WHITE);
              display->setCursor(0, 36);
              if(OF_Prefs::pins[OF_Const::solenoidPin] >= 0 && OF_Prefs::pins[OF_Const::solenoidSwitch] == -1) {
                display->println(" Solenoid Toggle ");
                display->setTextColor(WHITE, BLACK);
                display->setCursor(0, 47);
                display->println(" Send Escape Keypress");
              } else {
                display->println(" Send Escape Keypress");
                display->setTextColor(WHITE, BLACK);
                display->setCursor(0, 47);
                display->println("Calibrate");
              }
            } else if(OF_Prefs::pins[OF_Const::solenoidPin] >= 0 && OF_Prefs::pins[OF_Const::solenoidSwitch] == -1) {
              display->println(" Save Gun Settings");
              display->setTextColor(BLACK, WHITE);
              display->setCursor(0, 36);
              display->println(" Solenoid Toggle ");
              display->setTextColor(WHITE, BLACK);
              display->setCursor(0, 47);
              display->println(" Send Escape Keypress");
            } else {
              display->println(" Send Escape Keypress");
              display->setTextColor(BLACK, WHITE);
              display->setCursor(0, 36);
              display->println(" Calibrate ");
              display->setTextColor(WHITE, BLACK);
              display->setCursor(0, 47);
              display->println(" Profile Select ");
            }
            break;
          case ScreenPause_EscapeKey:
            display->setTextColor(WHITE, BLACK);
            display->setCursor(0, 25);
            if(OF_Prefs::pins[OF_Const::solenoidPin] >= 0 && OF_Prefs::pins[OF_Const::solenoidSwitch] == -1) {
              display->println(" Solenoid Toggle ");
              display->setTextColor(BLACK, WHITE);
              display->setCursor(0, 36);
              display->println(" Send Escape Keypress");
              display->setTextColor(WHITE, BLACK);
              display->setCursor(0, 47);
              display->println(" Calibrate ");
            } else if(OF_Prefs::pins[OF_Const::rumblePin] >= 0 && OF_Prefs::pins[OF_Const::rumbleSwitch] == -1) {
              display->println(" Rumble Toggle ");
              display->setTextColor(BLACK, WHITE);
              display->setCursor(0, 36);
              display->println(" Send Escape Keypress");
              display->setTextColor(WHITE, BLACK);
              display->setCursor(0, 47);
              display->println(" Calibrate ");
            } else {
              display->println(" Save Gun Settings ");
              display->setTextColor(BLACK, WHITE);
              display->setCursor(0, 36);
              display->println(" Send Escape Keypress");
              display->setTextColor(WHITE, BLACK);
              display->setCursor(0, 47);
              display->println(" Calibrate ");
            }
            break;
        }
        display->display();
    }
}

void ExtDisplay::PauseProfileUpdate(const int &selection, const char* name1, const char* name2, const char* name3, const char* name4)
{
    if(display != nullptr) {
        display->fillRect(0, 16, 128, 48, BLACK);
        display->drawBitmap(60, 18, upArrowGlyph, ARROW_WIDTH, ARROW_HEIGHT, WHITE);
        display->drawBitmap(60, 59, downArrowGlyph, ARROW_WIDTH, ARROW_HEIGHT, WHITE);
        display->setTextSize(1);
        switch(selection) {
          case 0: // Profile #0, etc.
            display->setTextColor(WHITE, BLACK);
            display->setCursor(4, 25);
            display->println(name4);
            display->setTextColor(BLACK, WHITE);
            display->setCursor(4, 36);
            display->println(name1);
            display->setTextColor(WHITE, BLACK);
            display->setCursor(4, 47);
            display->println(name2);
            break;
          case 1:
            display->setTextColor(WHITE, BLACK);
            display->setCursor(4, 25);
            display->println(name1);
            display->setTextColor(BLACK, WHITE);
            display->setCursor(4, 36);
            display->println(name2);
            display->setTextColor(WHITE, BLACK);
            display->setCursor(4, 47);
            display->println(name3);
            break;
          case 2:
            display->setTextColor(WHITE, BLACK);
            display->setCursor(4, 25);
            display->println(name2);
            display->setTextColor(BLACK, WHITE);
            display->setCursor(4, 36);
            display->println(name3);
            display->setTextColor(WHITE, BLACK);
            display->setCursor(4, 47);
            display->println(name4);
            break;
          case 3:
            display->setTextColor(WHITE, BLACK);
            display->setCursor(4, 25);
            display->println(name3);
            display->setTextColor(BLACK, WHITE);
            display->setCursor(4, 36);
            display->println(name4);
            display->setTextColor(WHITE, BLACK);
            display->setCursor(4, 47);
            display->println(name1);
            break;
        }
        display->display();
    }
}

void ExtDisplay::SaveScreen(const int &status)
{
    if(display != nullptr) {
        display->fillRect(0, 16, 128, 48, BLACK);
        display->setTextColor(WHITE, BLACK);
        display->setTextSize(2);
        display->setCursor(24, 24);
        display->println("Saving...");
        display->display();
    }
}

void ExtDisplay::PrintAmmo(const uint &ammo)
{
    if(display != nullptr) {
        currentAmmo = ammo;

        // use the rounding error to get the left & right digits
        uint ammoLeft = ammo / 10;
        uint ammoRight = ammo - (ammoLeft * 10);

        ammoEmpty = ammo ? false : true;

        if(screenState == Screen_Mamehook_Single) {
            display->fillRect(40, 22, (NUMBER_GLYPH_WIDTH*2)+6, NUMBER_GLYPH_HEIGHT, BLACK);
            display->drawBitmap(40,                      22, numbers[ammoLeft],  NUMBER_GLYPH_WIDTH, NUMBER_GLYPH_HEIGHT, WHITE);
            display->drawBitmap(40+6+NUMBER_GLYPH_WIDTH, 22, numbers[ammoRight], NUMBER_GLYPH_WIDTH, NUMBER_GLYPH_HEIGHT, WHITE);
        } else if(screenState == Screen_Mamehook_Dual) {
            display->fillRect(72, 22, (NUMBER_GLYPH_WIDTH*2)+6, NUMBER_GLYPH_HEIGHT, BLACK);
            display->drawBitmap(72,                      22, numbers[ammoLeft],  NUMBER_GLYPH_WIDTH, NUMBER_GLYPH_HEIGHT, WHITE);
            display->drawBitmap(72+6+NUMBER_GLYPH_WIDTH, 22, numbers[ammoRight], NUMBER_GLYPH_WIDTH, NUMBER_GLYPH_HEIGHT, WHITE);
        }

        display->display();
    }
}

void ExtDisplay::PrintLife(const uint &life)
{
    if(display != nullptr) {
        currentLife = life;
        lifeEmpty = life ? false : true;
        if(screenState == Screen_Mamehook_Single) {
            if(lifeBar) {
                display->fillRect(14, 37, 100, 9, BLACK);
                display->fillRect(52, 51, 30, 8, BLACK);
                display->fillRect(14, 37, life, 9, WHITE);

                if(life) {
                  display->setTextSize(1);
                  display->setCursor(52, 51);
                  display->setTextColor(WHITE, BLACK);
                  display->print(life);
                  display->println(" %");
                }

                display->display();
            } else {
                display->fillRect(22, 19, HEART_LARGE_WIDTH*5+4, HEART_LARGE_HEIGHT+22+HEART_LARGE_HEIGHT, BLACK);
                switch(life) {
                  case 9:
                    display->drawBitmap(22, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+1+HEART_LARGE_WIDTH, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+2+HEART_LARGE_WIDTH*2, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+3+HEART_LARGE_WIDTH*3, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+4+HEART_LARGE_WIDTH*4, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);

                    display->drawBitmap(30, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(30+1+HEART_LARGE_WIDTH, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(30+2+HEART_LARGE_WIDTH*2, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(30+3+HEART_LARGE_WIDTH*3, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    break;
                  case 8:
                    display->drawBitmap(22, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+1+HEART_LARGE_WIDTH, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+2+HEART_LARGE_WIDTH*2, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+3+HEART_LARGE_WIDTH*3, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+4+HEART_LARGE_WIDTH*4, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);

                    display->drawBitmap(39, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(39+1+HEART_LARGE_WIDTH, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(39+2+HEART_LARGE_WIDTH*2, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    break;
                  case 7:
                    display->drawBitmap(22, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+1+HEART_LARGE_WIDTH, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+2+HEART_LARGE_WIDTH*2, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+3+HEART_LARGE_WIDTH*3, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+4+HEART_LARGE_WIDTH*4, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);

                    display->drawBitmap(48, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(48+1+HEART_LARGE_WIDTH, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    break;
                  case 6:
                    display->drawBitmap(22, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+1+HEART_LARGE_WIDTH, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+2+HEART_LARGE_WIDTH*2, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+3+HEART_LARGE_WIDTH*3, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+4+HEART_LARGE_WIDTH*4, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);

                    display->drawBitmap(56, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    break;
                  case 5:
                    display->drawBitmap(22, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+HEART_LARGE_WIDTH, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+HEART_LARGE_WIDTH*2, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+HEART_LARGE_WIDTH*3, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+HEART_LARGE_WIDTH*4, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    break;
                  case 4:
                    display->drawBitmap(30, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(30+1+HEART_LARGE_WIDTH, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(30+2+HEART_LARGE_WIDTH*2, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(30+3+HEART_LARGE_WIDTH*3, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    break;
                  case 3:
                    display->drawBitmap(39, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(39+1+HEART_LARGE_WIDTH, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(39+2+HEART_LARGE_WIDTH*2, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    break;
                  case 2:
                    display->drawBitmap(48, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(48+1+HEART_LARGE_WIDTH, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    break;
                  case 1:
                    display->drawBitmap(56, 30, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    break;
                  case 0:
                    break;
                  default: // 10+
                    display->drawBitmap(22, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+1+HEART_LARGE_WIDTH, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+2+HEART_LARGE_WIDTH*2, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+3+HEART_LARGE_WIDTH*3, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+4+HEART_LARGE_WIDTH*4, 19, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);

                    display->drawBitmap(22, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+1+HEART_LARGE_WIDTH, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+2+HEART_LARGE_WIDTH*2, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+3+HEART_LARGE_WIDTH*3, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    display->drawBitmap(22+4+HEART_LARGE_WIDTH*4, 41, lifeIcoLarge, HEART_LARGE_WIDTH, HEART_LARGE_HEIGHT, WHITE);
                    break;
                }
                display->display();
            }
        } else if(screenState == Screen_Mamehook_Dual) {
            if(lifeBar) {
                display->fillRect(4, 39, 55, 5, BLACK);
                display->fillRect(20, 51, 30, 8, BLACK);
                display->fillRect(4, 39, map(life, 0, 100, 0, 55), 5, WHITE);

                if(life) {
                  display->setTextSize(1);
                  display->setCursor(20, 51);
                  display->setTextColor(WHITE, BLACK);
                  display->print(life);
                  display->println(" %");
                }

                display->display();
            } else {
                display->fillRect(1, 22, HEART_SMALL_WIDTH*5, HEART_SMALL_HEIGHT+20+HEART_SMALL_HEIGHT, BLACK);
                switch(life) {
                  case 9:
                    display->drawBitmap(1, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*2, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*3, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*4, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);

                    display->drawBitmap(7, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(7+HEART_SMALL_WIDTH, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(7+HEART_SMALL_WIDTH*2, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(7+HEART_SMALL_WIDTH*3, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    break;
                  case 8:
                    display->drawBitmap(1, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*2, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*3, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*4, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);

                    display->drawBitmap(13, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(13+HEART_SMALL_WIDTH, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(13+HEART_SMALL_WIDTH*2, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    break;
                  case 7:
                    display->drawBitmap(1, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*2, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*3, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*4, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);

                    display->drawBitmap(19, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(19+HEART_SMALL_WIDTH, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    break;
                  case 6:
                    display->drawBitmap(1, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*2, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*3, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*4, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);

                    display->drawBitmap(25, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    break;
                  case 5:
                    display->drawBitmap(1, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*2, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*3, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*4, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    break;
                  case 4:
                    display->drawBitmap(7, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(7+HEART_SMALL_WIDTH, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(7+HEART_SMALL_WIDTH*2, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(7+HEART_SMALL_WIDTH*3, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    break;
                  case 3:
                    display->drawBitmap(13, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(13+HEART_SMALL_WIDTH, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(13+HEART_SMALL_WIDTH*2, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    break;
                  case 2:
                    display->drawBitmap(19, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(19+HEART_SMALL_WIDTH, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    break;
                  case 1:
                    display->drawBitmap(25, 32, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    break;
                  case 0:
                    break;
                  default: // 10+
                    display->drawBitmap(1, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*2, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*3, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*4, 22, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);

                    display->drawBitmap(1, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*2, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*3, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    display->drawBitmap(1+HEART_SMALL_WIDTH*4, 42, lifeIcoSmall, HEART_SMALL_WIDTH, HEART_SMALL_HEIGHT, WHITE);
                    break;
                }
                display->display();
            }
        }
    }
}
