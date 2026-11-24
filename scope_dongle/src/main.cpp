/*!
 * @file main.cpp
 * @brief main SCOPE_DONGLE
 * @n CPP main SCOPE_DONGLE
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
// The scope dongle uses the same ST7735 80x160 display of the lightgun dongle, driven by LovyanGFX.
// Preprocessor directives (#ifdef) keep in the final binary only the code
// needed by the target board.
// /
// Il dongle del mirino usa lo stesso display ST7735 80x160 del dongle della lightgun, pilotato da LovyanGFX.
// Le direttive preprocessore (#ifdef) lasciano nel binario finale solo il codice
// necessario alla board target.

#ifdef USES_DISPLAY
  #include "OpenFIRE_logo.h"
  #ifdef USE_LOVYAN_GFX
    #define LGFX_USE_V1
    #include <LovyanGFX.hpp>
    #include "LGFX_096_ST7735S_80x160.hpp"

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

//#include "OpenFIRE-SCOPE-DONGLE-version.h"

#include "OpenFIRE_Wireless_Scope.h"

/*
CPU frequency: uncomment to force it (80, 160 or 240 MHz). Commented out = 240 MHz.
/
Frequenza della CPU: decommentare per forzarla (80, 160 o 240 MHz). Commentato = 240 MHz.
*/
//#define SCOPE_DONGLE_CPU_FREQ_MHZ 240

#ifdef USES_DISPLAY
  bool display_init = false;

  // Rotation of the display, as in the lightgun dongle / Rotazione del display, come nel dongle della lightgun
  #define SCOPE_DONGLE_DISPLAY_ROTATION 3

  // How long the OpenFIRE logo stays on screen at power on / Quanto resta a schermo il logo OpenFIRE all'accensione
  #define SCOPE_DONGLE_LOGO_TIME_MS 3000
#endif // USES_DISPLAY

// EN: Channel and state of the dongle on the display, to the right of the small logo: "FREE" or the paired player.
// IT: Canale e stato del dongle sul display, a destra del logo piccolo: "FREE" oppure il giocatore abbinato.
void showScopeDongleState() {
  #ifdef USES_DISPLAY
    #ifdef USE_LOVYAN_GFX
      if (!display_init) return;
      tft.fillRect(64, 20, 96, 60, BLACK);
      tft.setCursor(70, 28);
      tft.printf("CH %u", WirelessScope.channel);
      tft.setCursor(70, 52);
      if (WirelessScope.link_state >= SCOPE_LINK_STATE::SCOPE_LINK_PAIRED) tft.printf("P%u", WirelessScope.player);
      else tft.print("FREE");
    #endif // USE_LOVYAN_GFX
  #endif // USES_DISPLAY
}

// EN: Messages to the PC daemon over the USB (channel chosen, paired with player and scope address, free again):
//     see the "Specifica del demone PC OpenFIRE".
// IT: Messaggi al demone del PC sulla USB (canale scelto, abbinato con giocatore e indirizzo del mirino, tornato libero):
//     vedi la "Specifica del demone PC OpenFIRE".
void notifyScopeDongleState() {
  // TO DO / DA FARE: OF_SCOPE:... with WirelessScope.channel, WirelessScope.player, WirelessScope.ip_peer

}

// ===================================================================================
// MAIN SETUP
// ===================================================================================

// The main show!
void setup() { 

  #ifdef SCOPE_DONGLE_CPU_FREQ_MHZ
    setCpuFrequencyMhz(SCOPE_DONGLE_CPU_FREQ_MHZ);
  #endif // SCOPE_DONGLE_CPU_FREQ_MHZ

  // 1. HARDWARE INITIALIZATION / INIZIALIZZAZIONE HARDWARE (Video Feedback)
  #ifdef USES_DISPLAY
    #ifdef USE_LOVYAN_GFX
      display_init = tft.init();
      if(display_init) {
        tft.setSwapBytes(true);
        tft.setColorDepth(16);
        tft.setRotation(SCOPE_DONGLE_DISPLAY_ROTATION);

        // OpenFIRE logo, same position as the lightgun dongle / logo OpenFIRE, stessa posizione del dongle della lightgun
        tft.fillScreen(BLACK);
        tft.pushImage(10, 1, LOGO_RGB_ALPHA_OPEN_WIDTH, LOGO_RGB_ALPHA_OPEN_HEIGHT, (uint16_t *)logo_rgb_alpha_open);
        vTaskDelay(pdMS_TO_TICKS(SCOPE_DONGLE_LOGO_TIME_MS));

        // Banner and small logo, shown during the channel search / banner e logo piccolo, mostrati durante la ricerca del canale
        tft.fillScreen(BLACK);
        tft.drawBitmap(40, 0, customSplashBanner, CUSTSPLASHBANN_WIDTH, CUSTSPLASHBANN_HEIGHT, BLUE);
        tft.pushImage(9, 24, LOGO_RGB_ALPHA_WIDTH, LOGO_RGB_ALPHA_HEIGHT, (uint16_t *)logo_rgb_alpha);
        tft.setTextSize(2);
        tft.setTextColor(WHITE);
      }
    #endif // USE_LOVYAN_GFX
  #endif // USES_DISPLAY

  // ===================================================================================
  // WIRELESS / WIRELESS
  // ===================================================================================
  // The dongle chooses the freest channel, starts its access point and waits for a scope.
  // /
  // Il dongle sceglie il canale piu' libero, avvia il suo access point e aspetta un mirino.
  WirelessScope.begin();                    // channel search and access point / ricerca del canale e access point
  showScopeDongleState();                   // channel, free / canale, libero
  notifyScopeDongleState();

  WirelessScope.connection_scope_dongle();  // waits for a scope and pairs (blocking) / aspetta un mirino e si abbina (bloccante)
  showScopeDongleState();                   // paired player / giocatore abbinato
  notifyScopeDongleState();

}

// ===================================================================================
// MAIN LOOP
// ===================================================================================

void loop()
{
  // EN: The scope left the network: the dongle is free and visible again and waits for a new pairing.
  // IT: Il mirino e' uscito dalla rete: il dongle torna libero e visibile e aspetta un nuovo abbinamento.
  if (WirelessScope.link_state == SCOPE_LINK_STATE::SCOPE_LINK_SEARCH) {
    showScopeDongleState();
    notifyScopeDongleState();
    WirelessScope.connection_scope_dongle();
    showScopeDongleState();
    notifyScopeDongleState();
  }

  vTaskDelay(pdMS_TO_TICKS(1)); // Let's give the microcontroller some breathing room to avoid resets triggered by the Task Watchdog. / diamo un po' di respiro al micro per non incorrere in reset dal Task Watchdog
}
