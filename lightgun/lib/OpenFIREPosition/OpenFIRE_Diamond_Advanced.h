#ifdef USE_DIAMOND_ADVANCED
/*!
 * @file OpenFIRE_Diamond_Advanced.h
 * @brief Light Gun library for 4 LED setup
 * @n CPP file for Samco Light Gun 4 LED setup (Diamond layout)
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
 * @copyright Samco, https://github.com/samuelballantyne, 2024
 * @copyright GNU Lesser General Public License
 *
 * @author [Sam Ballantyne](samuelballantyne@hotmail.com)
 * @version V1.0
 * @date 2024
 */

#ifndef _OpenFIRE_Diamond_Advanced_h_
#define _OpenFIRE_Diamond_Advanced_h_

#include <stdint.h>
#include "OpenFIREConst.h"
#include <OpenFIRECameraProfile.h>
#include "OpenFIRE_Anchor.h"

// Stessa interfaccia pubblica di OpenFIRE_Diamond.h, cosi' OpenFIREcommon non cambia.
// Indici dei LED (spazio camera grezzo, X non specchiata):
// 0 = alto, 1 = X massima, 2 = basso, 3 = X minima.
class OpenFIRE_Diamond {

    unsigned int see[4] = {0, 0, 0, 0};

    int medianY = 0;
    int medianX = 0;

    // Geometria pubblicata (con la molla applicata)
    int FinalX[4] = {0, 0, 0, 0};
    int FinalY[4] = {0, 0, 0, 0};

    float angle = 0;
    float height2 = 0;
    float width2 = 0;

    unsigned int start = 0;

    unsigned int seenFlags = 0;

    uint8_t camToMouseShift = 0;
    int camToMouseMult = 1;

    // --- STATO DEL TRACKING ---
    // Geometria del frame precedente senza molla (misurata o ricostruita)
    float prevGeomX[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float prevGeomY[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    // Ultima geometria vista con tutti e 4 i LED (ancora esatta per ricostruzioni a 2-3 LED)
    float anchorX[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float anchorY[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    bool anchor_ok = false;
    uint8_t prevNumSeen = 0;        // LED visti nel frame precedente accettato (0 = buio)
    uint8_t darkFrames = 255;       // fotogrammi di buio dopo l'ultimo frame accettato, fermo a 255
    uint8_t shortDarkFrames = 8;    // durata massima del buio breve, in fotogrammi (dagli fps del profilo)
    float minLedDist2 = 64.0f;      // distanza minima al quadrato tra 2 LED (un pixel del sensore)
    uint8_t prevSeenMask = 0;       // quali indici erano misurati nel frame precedente

    // --- VARIABILI MOTORE KINEMATICO (MOLLA) ---
    // "Debito" tra la geometria calcolata e quella mostrata, creato quando cambia
    // l'insieme dei LED misurati e riassorbito gradualmente (come in Square_Advanced).
    float offset_X[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float offset_Y[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    // true se il debito attuale e' nato da un aumento dei LED misurati (ritorno proporzionale)
    bool spring_fast = false;
    // Movimento medio filtrato di ogni LED, usato per la velocita' di ritorno della molla
    float motion_X[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float motion_Y[4] = {0.0f, 0.0f, 0.0f, 0.0f};

    float FPS_NORMALIZATION = 1.0f;
    float COSTANTE_MOLLA = 0.25f;

    // Rapporto atteso diagonale orizzontale / verticale (dal profilo: segue il monitor).
    // Serve all'aggancio a freddo; non si aggiorna dal rapporto apparente deformato dalla prospettiva.
    float diagRatioRef = 1.0f;

    // Focale nominale e centro dell'immagine (unita' mouse): ricostruzione dei LED nascosti
    AnchorReconstruction anchorReconstruction;
    float focal = 0.0f;
    float imageCX = 0.0f;
    float imageCY = 0.0f;

    static bool ConvexCycle(const float* x, const float* y);

public:

    void configure(const CameraProfile& profile);

    /// @brief Main function to calculate X, Y, and H
    void begin(const int* px, const int* py, unsigned int seen);

    int X(int index) const { return FinalX[index]; }
    int Y(int index) const { return FinalY[index]; }
    unsigned int testSee(int index) const { return see[index]; }
    int testMedianX() const { return medianX; }
    int testMedianY() const { return medianY; }

    /// @brief Height
    float H() const { return height2; }

    /// @brief Height
    float W() const { return width2; }

    /// @brief Angle
    float Ang() const { return angle; }

    /// @brief Bit mask of positions the camera saw
    unsigned int seen() const { return seenFlags; }
};

#endif // _OpenFIRE_Diamond_Advanced_h_
#endif // USE_DIAMOND_ADVANCED
