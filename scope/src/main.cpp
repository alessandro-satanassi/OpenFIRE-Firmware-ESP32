/*!
 * @file main.cpp
 * @brief main SCOPE
 * @n CPP main SCOPE
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#include <Arduino.h>
#include <SPI.h>

// ===================================================================================
// HARDWARE ABSTRACTION: OPTIONAL GRAPHICS MANAGEMENT / GESTIONE GRAFICA OPZIONALE
// ===================================================================================
// The scope uses a GC9A01 240x240 round display driven by LovyanGFX.
// Preprocessor directives (#ifdef) keep in the final binary only the code
// needed by the target board.
// /
// Lo scope usa un display rotondo GC9A01 240x240 pilotato da LovyanGFX.
// Le direttive preprocessore (#ifdef) lasciano nel binario finale solo il codice
// necessario alla board target.

#ifdef USES_DISPLAY
  #include "OpenFIRE_logo.h"
  #ifdef USE_LOVYAN_GFX
    #define LGFX_USE_V1
    #include <LovyanGFX.hpp>
    #include "LGFX_128_GC9A01_240x240.hpp"

    LGFX tft;

    #define WHITE            TFT_WHITE
    #define BLACK            TFT_BLACK
    #define BLUE             TFT_BLUE
    #define RED              TFT_RED
    #define MAGENTA          TFT_MAGENTA
    #define GREEN            TFT_GREEN
    #define CYAN             TFT_CYAN
    #define YELLOW           TFT_YELLOW
    #define BROWN            TFT_BROWN
    #define GRAY             TFT_LIGHTGRAY
    #define ORANGE           TFT_ORANGE
  #endif // USE_LOVYAN_GFX

#endif // USES_DISPLAY

//#include "OpenFIRE-SCOPE-version.h"

/*
CPU frequency: uncomment to force it (80, 160 or 240 MHz). Commented out = 240 MHz.
/
Frequenza della CPU: decommentare per forzarla (80, 160 o 240 MHz). Commentato = 240 MHz.
*/
//#define SCOPE_CPU_FREQ_MHZ 240

#ifdef USES_DISPLAY
  bool display_init = false;

  // Rotation of the display: 0 -> none, 1 -> 90°, 2 -> 180°, 3 -> 270° / Rotazione del display: 0 -> nessuna, 1 -> 90°, 2 -> 180°, 3 -> 270°
  #define SCOPE_DISPLAY_ROTATION 0

  // How long the OpenFIRE logo stays on screen at power on / Quanto resta a schermo il logo OpenFIRE all'accensione
  #define SCOPE_LOGO_TIME_MS 3000
#endif // USES_DISPLAY

// ===================================================================================
// MAIN SETUP
// ===================================================================================

// The main show!
void setup() { 

  #ifdef SCOPE_CPU_FREQ_MHZ
    setCpuFrequencyMhz(SCOPE_CPU_FREQ_MHZ);
  #endif // SCOPE_CPU_FREQ_MHZ

  // 1. HARDWARE INITIALIZATION / INIZIALIZZAZIONE HARDWARE (Video Feedback)
  #ifdef USES_DISPLAY
    #ifdef USE_LOVYAN_GFX
      display_init = tft.init();
      if(display_init) {
        tft.setSwapBytes(true);
        tft.setColorDepth(16);
        tft.setRotation(SCOPE_DISPLAY_ROTATION);
        tft.fillScreen(BLACK);

        // OpenFIRE logo in the centre of the round display / logo OpenFIRE al centro del display rotondo
        tft.pushImage((tft.width()  - LOGO_RGB_ALPHA_OPEN_WIDTH)  / 2,
                      (tft.height() - LOGO_RGB_ALPHA_OPEN_HEIGHT) / 2,
                      LOGO_RGB_ALPHA_OPEN_WIDTH, LOGO_RGB_ALPHA_OPEN_HEIGHT, (uint16_t *)logo_rgb_alpha_open);

        vTaskDelay(pdMS_TO_TICKS(SCOPE_LOGO_TIME_MS));
      }
    #endif // USE_LOVYAN_GFX
  #endif // USES_DISPLAY

  // ===================================================================================
  // WIRELESS / WIRELESS
  // ===================================================================================
  // Next step / prossimo passo

}

// ===================================================================================
// MAIN LOOP
// ===================================================================================

void loop()
{
  vTaskDelay(pdMS_TO_TICKS(1)); // Let's give the microcontroller some breathing room to avoid resets triggered by the Task Watchdog. / diamo un po' di respiro al micro per non incorrere in reset dal Task Watchdog
}
