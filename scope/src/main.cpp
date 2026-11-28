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

#include "OpenFIRE-SCOPE-version.h"

#include "OpenFIRE_Wireless_Scope.h"
#include "OpenFIRE_Scope_Video.h"   // reception of the frames, rotation and display / ricezione dei frame, rotazione e display

// ===================================================================================
// PLAYER / GIOCATORE
// ===================================================================================
// EN: Player number (1-4) sent to the dongle in the Confirm of the pairing. 0 = not received yet.
//     It arrives from the S3 of the lightgun in the first valid Player over the UART and is fixed at the
//     pairing: a hot change of player does not exist.
// IT: Numero del giocatore (1-4) mandato al dongle nella Conferma dell'abbinamento. 0 = non ancora ricevuto.
//     Arriva dall'S3 della lightgun nel primo Giocatore valido sulla UART e si fissa all'abbinamento:
//     il cambio di giocatore a caldo non esiste.
volatile uint8_t scope_player = 0;
volatile bool scope_display_on = true;   // from the Display of the S3 / dal Display dell'S3

// ===================================================================================
// ROLL OF THE LIGHTGUN / ROLL DELLA LIGHTGUN
// ===================================================================================
// EN: Rotation of the lightgun in degrees, from -90 to +90 (0 = no rotation), positive = turned clockwise as seen
//     by the shooter. The video turns every frame the opposite way, with the last value at the start of the rotation.
//     Values outside -90..+90 are taken as -90 or +90. It arrives from the S3 of the lightgun over the UART, only when it changes.
// IT: Rotazione della lightgun in gradi, da -90 a +90 (0 = nessuna rotazione), positivo = girata in senso orario
//     vista da chi spara. Il video ruota ogni frame al contrario, con l'ultimo valore all'inizio della rotazione.
//     I valori fuori da -90..+90 valgono -90 o +90. Arriva dall'S3 della lightgun sulla UART, solo quando cambia.
volatile int8_t scope_roll = 0;

// EN: true while the sprite of the turning logo (51 KB) is in memory: the video takes its memory only after it.
// IT: true finche' lo sprite del logo che gira (51 KB) e' in memoria: il video prende la sua memoria solo dopo.
volatile bool scope_logo_spin_on = false;

// EN: Bench tests without the S3: fixed player, the pairing starts at once. Comment it out when the UART is ready.
// IT: Prove al banco senza l'S3: giocatore fisso, l'abbinamento parte subito. Commentarlo quando la UART e' pronta.
#define SCOPE_TEST_PLAYER 1

// ===================================================================================
// UART WITH THE ESP32-S3 OF THE LIGHTGUN / UART CON L'ESP32-S3 DELLA LIGHTGUN
// ===================================================================================
// EN: 115200 baud 8N1, same frame in both directions: 0xA5, type, length, data, CRC-8
//     (polynomial 0x9B, initial value 0, no reflection, no final XOR, on type, length and data, without 0xA5).
//     Types 0x01-0x7F from the S3 to the scope, 0x81-0xFF from the scope to the S3. Unknown types are ignored.
//     The S3 starts: after its USB or dongle connection it sends Hello (at most 3 times), the scope answers Present,
//     then the S3 sends the Player once and the Roll only when it changes. No keepalive, no answer to the Roll.
// IT: 115200 baud 8N1, stessa trama nei due sensi: 0xA5, tipo, lunghezza, dati, CRC-8
//     (polinomio 0x9B, valore iniziale 0, senza riflessione, senza XOR finale, su tipo, lunghezza e dati, senza lo 0xA5).
//     Tipi da 0x01 a 0x7F dall'S3 al mirino, da 0x81 a 0xFF dal mirino all'S3. I tipi sconosciuti si ignorano.
//     Comincia l'S3: dopo il suo collegamento USB o dongle manda Ciao (al massimo 3 volte), il mirino risponde Presente,
//     poi l'S3 manda una volta il Giocatore e il Roll solo quando cambia. Niente keepalive, nessuna risposta al Roll.
#if defined(SCOPE_UART_TX_GPIO) && defined(SCOPE_UART_RX_GPIO)
  #define SCOPE_UART               Serial1
  #define SCOPE_UART_BAUD          115200
  #define SCOPE_UART_START         0xA5
  #define SCOPE_UART_PROTOCOL      1      // protocol version / versione del protocollo
  #define SCOPE_UART_MAX_DATA      16

  // S3 -> scope / S3 -> mirino
  #define SCOPE_UART_CMD_HELLO     0x01   // protocol version / versione protocollo
  #define SCOPE_UART_CMD_PLAYER    0x02   // player (1-4) / giocatore (1-4)
  #define SCOPE_UART_CMD_ROLL      0x03   // int8 degrees -90..+90 (see scope_roll) / gradi int8 -90..+90 (vedi scope_roll)
  #define SCOPE_UART_CMD_DISPLAY   0x04   // 0 off, 1 on / 0 spento, 1 acceso
  // scope -> S3 / mirino -> S3
  #define SCOPE_UART_MSG_PRESENT   0x81   // protocol version, firmware major, minor, patch / versione protocollo, firmware major, minor, patch
  #define SCOPE_UART_MSG_INPUTS    0x83   // reserved: buttons and wheel / riservato: tasti e rotellina

  static uint8_t scopeUartCrc8Step(uint8_t crc, uint8_t value) {
    crc ^= value;
    for (uint8_t bit = 0; bit < 8; ++bit)
      crc = (crc & 0x80) ? (uint8_t)((crc << 1) ^ 0x9B) : (uint8_t)(crc << 1);
    return crc;
  }

  static void scopeUartSend(uint8_t type, const uint8_t *data, uint8_t length) {
    uint8_t frame[4 + SCOPE_UART_MAX_DATA];
    frame[0] = SCOPE_UART_START;
    frame[1] = type;
    frame[2] = length;
    memcpy(&frame[3], data, length);
    uint8_t crc = 0;
    for (uint8_t i = 1; i < 3 + length; ++i) crc = scopeUartCrc8Step(crc, frame[i]);
    frame[3 + length] = crc;
    SCOPE_UART.write(frame, 4 + length);
  }

  static void scopeUartHandle(uint8_t type, const uint8_t *data, uint8_t length) {
    switch (type) {
      case SCOPE_UART_CMD_HELLO: {
        const uint8_t present[4] = { SCOPE_UART_PROTOCOL, OPENFIRE_SCOPE_VERSION_MAJOR,
                                     OPENFIRE_SCOPE_VERSION_MINOR, OPENFIRE_SCOPE_VERSION_PATCH };
        scopeUartSend(SCOPE_UART_MSG_PRESENT, present, sizeof(present));
        break;
      }
      case SCOPE_UART_CMD_PLAYER:
        // EN: only the first one: the player is fixed at the pairing / IT: solo il primo: il giocatore si fissa all'abbinamento
        if (length >= 1 && scope_player == 0 && data[0] >= 1 && data[0] <= 4) scope_player = data[0];
        break;
      case SCOPE_UART_CMD_ROLL:
        if (length >= 1) {
          const int8_t roll = (int8_t)data[0];
          scope_roll = (roll > 90) ? 90 : ((roll < -90) ? -90 : roll);
        }
        break;
      case SCOPE_UART_CMD_DISPLAY:
        if (length >= 1) scope_display_on = (data[0] != 0);
        break;
      default:   // unknown type: ignored / tipo sconosciuto: ignorato
        break;
    }
  }

  // EN: Own task: the UART keeps working also while the pairing (blocking) is running.
  //     Wait for 0xA5, then type, length, data and CRC-8; wrong CRC -> frame discarded.
  // IT: Task proprio: la UART continua a lavorare anche mentre l'abbinamento (bloccante) e' in corso.
  //     Aspetta 0xA5, poi tipo, lunghezza, dati e CRC-8; CRC sbagliato -> trama scartata.
  static void scopeUartTask(void *param) {
    SCOPE_UART.begin(SCOPE_UART_BAUD, SERIAL_8N1, SCOPE_UART_RX_GPIO, SCOPE_UART_TX_GPIO);
    uint8_t state = 0, type = 0, length = 0, count = 0, crc = 0;
    uint8_t data[SCOPE_UART_MAX_DATA];

    for (;;) {
      while (SCOPE_UART.available() > 0) {
        const uint8_t b = (uint8_t)SCOPE_UART.read();
        switch (state) {
          case 0:   // start / inizio
            if (b == SCOPE_UART_START) { crc = 0; state = 1; }
            break;
          case 1:   // type / tipo
            type = b; crc = scopeUartCrc8Step(crc, b); state = 2;
            break;
          case 2:   // length / lunghezza
            if (b > SCOPE_UART_MAX_DATA) { state = 0; break; }
            length = b; count = 0; crc = scopeUartCrc8Step(crc, b); state = length ? 3 : 4;
            break;
          case 3:   // data / dati
            data[count++] = b; crc = scopeUartCrc8Step(crc, b);
            if (count == length) state = 4;
            break;
          default:  // CRC
            if (b == crc) scopeUartHandle(type, data, length);
            state = 0;
            break;
        }
      }
      vTaskDelay(pdMS_TO_TICKS(1));
    }
  }
#endif // SCOPE_UART_TX_GPIO && SCOPE_UART_RX_GPIO

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
// PAIRING STATE ON THE DISPLAY / STATO DELL'ABBINAMENTO SUL DISPLAY
// ===================================================================================
// EN: After the logo the display is drawn only by this task. Under the logo a state line:
//       "WAIT"      yellow  waiting for the player from the S3
//       "SEARCH P1" cyan    searching the dongle: meanwhile the logo turns slowly
//       "P1 CH 36"  green   paired: the logo stops straight
//     Paired and without video for 500 ms (no valid pause) a red "NO SIGNAL" line appears above the logo.
//     Own task: the animation does not slow down the UART and the pairing (blocking) does not stop the animation.
//     The panel is shared with the video: while the video is on screen this task draws nothing, then redraws everything.
// IT: Dopo il logo il display lo disegna solo questo task. Sotto il logo una riga di stato:
//       "WAIT"      gialla   attesa del giocatore dall'S3
//       "SEARCH P1" azzurra  ricerca del dongle: intanto il logo gira lentamente
//       "P1 CH 36"  verde    abbinato: il logo si ferma dritto
//     Abbinato e senza video da 500 ms (nessuna pausa valida) compare sopra il logo una riga rossa "NO SIGNAL".
//     Task proprio: l'animazione non rallenta la UART e l'abbinamento (bloccante) non ferma l'animazione.
//     Il pannello e' condiviso con il video: mentre il video e' sullo schermo questo task non disegna, poi ridisegna tutto.
#if defined(USES_DISPLAY) && defined(USE_LOVYAN_GFX)
  #define SCOPE_STATE_Y            213     // centre of the state line, below the turning logo / centro della riga di stato, sotto il logo che gira
  #define SCOPE_NOSIGNAL_Y         27      // centre of the "NO SIGNAL" line, above the logo / centro della riga "NO SIGNAL", sopra il logo
  #define SCOPE_LOGO_SPIN_TURN_MS  6000    // one turn of the logo / un giro del logo
  #define SCOPE_LOGO_SPIN_FRAME_MS 40      // 25 frames per second / 25 fotogrammi al secondo
  #define SCOPE_LOGO_SPIN_SIZE     160     // square that holds the logo at any angle (diagonal 139x78 = 159) / quadrato che contiene il logo a ogni angolo (diagonale 139x78 = 159)

  // EN: Sprite of the turning logo (160x160x2 = 51 KB), only during the search, in the array of the video frame.
  // IT: Sprite del logo che gira (160x160x2 = 51 KB), solo durante la ricerca, nell'array del frame del video.
  LGFX_Sprite logo_spin(&tft);
  static_assert(SCOPE_LOGO_SPIN_SIZE * SCOPE_LOGO_SPIN_SIZE * 2 <= SCOPE_VIDEO_SHARED_BYTES, "logo sprite / sprite del logo");

  static void drawScopeLogoStraight() {
    tft.fillRect((tft.width() - SCOPE_LOGO_SPIN_SIZE) / 2, (tft.height() - SCOPE_LOGO_SPIN_SIZE) / 2,
                 SCOPE_LOGO_SPIN_SIZE, SCOPE_LOGO_SPIN_SIZE, BLACK);
    tft.pushImage((tft.width()  - LOGO_RGB_ALPHA_OPEN_WIDTH)  / 2,
                  (tft.height() - LOGO_RGB_ALPHA_OPEN_HEIGHT) / 2,
                  LOGO_RGB_ALPHA_OPEN_WIDTH, LOGO_RGB_ALPHA_OPEN_HEIGHT, (uint16_t *)logo_rgb_alpha_open);
  }

  static void drawScopeLogoTurned(float angle) {
    logo_spin.fillSprite(BLACK);
    logo_spin.pushImageRotateZoom(SCOPE_LOGO_SPIN_SIZE / 2.0f, SCOPE_LOGO_SPIN_SIZE / 2.0f,
                                  LOGO_RGB_ALPHA_OPEN_WIDTH / 2.0f, LOGO_RGB_ALPHA_OPEN_HEIGHT / 2.0f,
                                  angle, 1.0f, 1.0f,
                                  LOGO_RGB_ALPHA_OPEN_WIDTH, LOGO_RGB_ALPHA_OPEN_HEIGHT, (const uint16_t *)logo_rgb_alpha_open);
    logo_spin.pushSprite((tft.width() - SCOPE_LOGO_SPIN_SIZE) / 2, (tft.height() - SCOPE_LOGO_SPIN_SIZE) / 2);
  }

  static void drawScopeStateLine(uint8_t state, uint8_t player_number) {
    char text[16];
    uint16_t color;
    if (state == SCOPE_LINK_STATE::SCOPE_LINK_SEARCH) {
      snprintf(text, sizeof(text), "SEARCH P%u", player_number);
      color = CYAN;
    } else if (state >= SCOPE_LINK_STATE::SCOPE_LINK_PAIRED) {
      snprintf(text, sizeof(text), "P%u CH %u", WirelessScope.player, WirelessScope.channel);
      color = GREEN;
    } else {
      snprintf(text, sizeof(text), "WAIT");
      color = YELLOW;
    }
    tft.fillRect(0, SCOPE_STATE_Y - 10, tft.width(), 20, BLACK);
    tft.setTextSize(2);
    tft.setTextDatum(textdatum_t::middle_center);
    tft.setTextColor(color, BLACK);
    tft.drawString(text, tft.width() / 2, SCOPE_STATE_Y);
  }

  static void drawScopeNoSignalLine(bool no_signal) {
    tft.fillRect(0, SCOPE_NOSIGNAL_Y - 10, tft.width(), 20, BLACK);
    if (!no_signal) return;
    tft.setTextSize(2);
    tft.setTextDatum(textdatum_t::middle_center);
    tft.setTextColor(RED, BLACK);
    tft.drawString("NO SIGNAL", tft.width() / 2, SCOPE_NOSIGNAL_Y);
  }

  static void scopeDisplayTask(void *param) {
    uint8_t last_state = 0xFF;
    uint8_t last_player = 0xFF;
    bool last_no_signal = false;
    bool had_panel = true;   // the logo of the setup is already on screen / il logo del setup e' gia' sullo schermo
    bool spinning = false;
    float angle = 0.0f;

    for (;;) {
      // EN: Video on screen: nothing to draw; when it ends everything is redrawn. / IT: video sullo schermo: niente da disegnare; quando finisce si ridisegna tutto.
      if (!scopePanelTryLock()) {
        had_panel = false;
        vTaskDelay(pdMS_TO_TICKS(SCOPE_LOGO_SPIN_FRAME_MS));
        continue;
      }
      scopePanelPowerApply();   // display on/off as the Display of the S3 / display acceso/spento come il Display dell'S3

      const uint8_t state = WirelessScope.link_state;
      const uint8_t player_number = scope_player;
      const bool no_signal = scopeVideoNoSignal();

      if (!had_panel) {
        had_panel = true;
        tft.fillScreen(BLACK);
        if (!spinning) drawScopeLogoStraight();
        drawScopeStateLine(state, player_number);
        drawScopeNoSignalLine(no_signal);
        last_no_signal = no_signal;
      }

      if (no_signal != last_no_signal) {
        drawScopeNoSignalLine(no_signal);
        last_no_signal = no_signal;
      }

      const bool search = (state == SCOPE_LINK_STATE::SCOPE_LINK_SEARCH);
      if (search && !spinning) {
        // EN: Search: the sprite goes into the array of the video frame as soon as the video has given it back (tried
        //     at every turn); until then the logo stays still.
        // IT: Ricerca: lo sprite va nell'array del frame del video appena il video lo ha restituito (si prova a ogni
        //     giro); fino ad allora il logo resta fermo.
        uint8_t *mem = scopeLogoMemoryTake();
        if (mem != nullptr) {
          logo_spin.setColorDepth(16);
          logo_spin.setSwapBytes(true);
          logo_spin.setBuffer(mem, SCOPE_LOGO_SPIN_SIZE, SCOPE_LOGO_SPIN_SIZE);
          spinning = true;
          angle = 0.0f;
        }
      } else if (!search && spinning) {
        // EN: End of the search: the array goes back to the video and the logo stops straight.
        // IT: Fine della ricerca: l'array torna al video e il logo si ferma dritto.
        logo_spin.deleteSprite();   // external buffer: nothing freed / buffer esterno: niente liberato
        spinning = false;
        scopeLogoMemoryGive();
        drawScopeLogoStraight();
      }

      if (state != last_state || player_number != last_player) {
        drawScopeStateLine(state, player_number);
        last_state = state;
        last_player = player_number;
      }

      if (spinning) {
        angle += 360.0f * SCOPE_LOGO_SPIN_FRAME_MS / SCOPE_LOGO_SPIN_TURN_MS;
        if (angle >= 360.0f) angle -= 360.0f;
        drawScopeLogoTurned(angle);
      }

      scopePanelUnlock();
      vTaskDelay(pdMS_TO_TICKS(SCOPE_LOGO_SPIN_FRAME_MS));
    }
  }
#endif // USES_DISPLAY && USE_LOVYAN_GFX

// ===================================================================================
// MAIN SETUP
// ===================================================================================

// The main show!
void setup() { 

  #ifdef SCOPE_CPU_FREQ_MHZ
    setCpuFrequencyMhz(SCOPE_CPU_FREQ_MHZ);
  #endif // SCOPE_CPU_FREQ_MHZ

  // EN: UART with the S3 first of all: the Hello of the S3 is answered also during the logo.
  // IT: UART con l'S3 prima di tutto: il Ciao dell'S3 ha risposta anche durante il logo.
  #if defined(SCOPE_UART_TX_GPIO) && defined(SCOPE_UART_RX_GPIO)
    xTaskCreate(scopeUartTask, "scopeUart", 4096, NULL, 2, NULL);
  #endif // SCOPE_UART_TX_GPIO && SCOPE_UART_RX_GPIO

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

        // EN: From here the display is drawn only by the video and by the task of the state screens, one at a time.
        // IT: Da qui il display lo disegnano solo il video e il task delle schermate di stato, uno alla volta.
        scopeVideoBegin(&tft);
        xTaskCreate(scopeDisplayTask, "scopeDisplay", 4096, NULL, 1, NULL);
      }
    #endif // USE_LOVYAN_GFX
  #endif // USES_DISPLAY

  // ===================================================================================
  // UART WITH THE S3 / UART CON L'S3
  // ===================================================================================
  #ifdef SCOPE_TEST_PLAYER
    scope_player = SCOPE_TEST_PLAYER;
  #endif // SCOPE_TEST_PLAYER

  // ===================================================================================
  // WIRELESS / WIRELESS
  // ===================================================================================
  WirelessScope.begin();                     // Wi-Fi 6, 5 GHz, 20 MHz / Wi-Fi 6, 5 GHz, 20 MHz

  // EN: The pairing starts at the first valid Player from the S3, also with the display off.
  // IT: L'abbinamento parte al primo Giocatore valido dall'S3, anche a display spento.
  while (scope_player == 0) vTaskDelay(pdMS_TO_TICKS(10));

  WirelessScope.connection_scope(scope_player); // searches a free dongle and pairs (blocking) / cerca un dongle libero e si abbina (bloccante)

}

// ===================================================================================
// MAIN LOOP
// ===================================================================================

void loop()
{
  // EN: The scope left the dongle's network (link lost): same pairing again, with the same player.
  // IT: Il mirino e' uscito dalla rete del dongle (collegamento perso): stesso abbinamento di nuovo, con lo stesso giocatore.
  if (WirelessScope.link_state == SCOPE_LINK_STATE::SCOPE_LINK_SEARCH) WirelessScope.connection_scope(scope_player);

  vTaskDelay(pdMS_TO_TICKS(1)); // Let's give the microcontroller some breathing room to avoid resets triggered by the Task Watchdog. / diamo un po' di respiro al micro per non incorrere in reset dal Task Watchdog
}
