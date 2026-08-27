#ifdef USE_SQUARE_ADVANCED
/*!
 * @file OpenFIRE_Square_Advanced.h
 * @brief Light Gun library for 4 LED setup
 * @n CPP file for Samco Light Gun 4 LED setup
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

#ifndef _OpenFIRE_Square_Advanced_h_
#define _OpenFIRE_Square_Advanced_h_

#include <stdint.h>
#include "OpenFIREConst.h"
#include <OpenFIRECameraProfile.h>


class OpenFIRE_Square {
public:
    //================================================================
    // ENUM PUBBLICO E API
    //================================================================

    void configure(const CameraProfile& profile);

    // Metodo principale di elaborazione
    void begin(const int* px, const int* py, unsigned int seen);
    
    // Metodi 'Getter' per l'accesso allo stato finale
    int X(int index) const { return FinalX[index]; }
    int Y(int index) const { return FinalY[index]; }
    unsigned int testSee(int index) const { return see[index]; }
    int testMedianX() const { return medianX; }
    int testMedianY() const { return medianY; }
    float H() const { return height; }
    float W() const { return width; }
    float Ang() const { return angle; }
    unsigned int seen() const { return seenFlags; }

    void setCalibrationMode(bool enabled) { calibrationMode = enabled; }
    //void setWideLayout(bool enabled) { wideLayout = enabled; }

    void setWideLayout(bool enabled) {
        if (wideLayout != enabled) {
            wideLayout = enabled;

            start = 0;
            is_tracking_stable = false;
            prev_num_points_seen = 0;
            prev_point_seen_mask = 0;
            current_point_seen_mask = 0;

            for (uint8_t i = 0; i < 4; i++) {
                offset_X[i] = 0.0f;
                offset_Y[i] = 0.0f;
            }
        }
    }


private:
    //================================================================
    // STATO DEL TRACKER
    //================================================================
    
    // Le coordinate finali sono conservate in interi perché rappresentano i pixel virtuali
    // del mouse. Vengono traslate usando il bit-shift (CamToMouseShift) per garantire
    // la massima velocità della ALU rispetto all'uso dei float.
    /*
    int  FinalX[4] = {400 * CamToMouseMult, 623 * CamToMouseMult, 400 * CamToMouseMult, 623 * CamToMouseMult};
    int  FinalY[4] = {200 * CamToMouseMult, 200 * CamToMouseMult, 568 * CamToMouseMult, 568 * CamToMouseMult};
    */
    int FinalX[4] = {0, 0, 0, 0};
    int FinalY[4] = {0, 0, 0, 0};
   
    int  medianX = 0;
    int  medianY = 0;
    
    // Registro a scorrimento (shift register) ereditato dalle API originali di Samco.
    // Viene mantenuto per non rompere la retrocompatibilità con funzioni esterne
    // che si aspettano di leggere uno storico di stabilità a 5 frame per singolo sensore.
    unsigned int see[4] = {0, 0, 0, 0};
    
    /*
    float height = (568 * CamToMouseMult) - (200 * CamToMouseMult);
    float width = (623 * CamToMouseMult) - (400 * CamToMouseMult);
    */
    float height = (float)(FinalY[2] - FinalY[0]);
    float width  = (float)(FinalX[1] - FinalX[0]);  
     
    float angle = 0;

    unsigned int start = 0;
    unsigned int seenFlags = 0;

    // is_tracking_stable agisce da "frizione": impedisce al motore predittivo (1-2 punti)
    // di inventare geometrie se nel frame precedente eravamo ciechi. Richiede un aggancio
    // solido (almeno 3 punti) prima di consentire le predizioni parziali.
    bool is_tracking_stable = false; 
    
    // ideal_aspect_ratio conserva l'ultima proporzione reale nota del rettangolo IR.
    // È vitale per la gestione a 2 punti, dove dobbiamo decidere se stiamo guardando
    // il lato lungo o corto del sensore senza avere tutti i riferimenti ottici.
    /*
    float ideal_aspect_ratio = 16.0f / 9.0f; 
    */
    float ideal_aspect_ratio = IR_SQUARE_BASE_HEIGHT_RATIO;
    
    uint8_t prev_num_points_seen = 0; 
    


    // Usiamo maschere di bit (bitmask) per lo storico visivo perché permettono di
    // verificare stati complessi e fare confronti incrociati in un singolo ciclo di clock
    // usando operatori bit a bit (&, |), risultando enormemente più veloci degli array.
    uint8_t prev_point_seen_mask = 0b00000000; // A->00001000  B->00000100  C->00000010  D->00000001
    uint8_t current_point_seen_mask = 0;

    // A---B  //TxIR superiori orizzontali 
    // |   |  // la base AB-CD è sempre più corta dell'altezza AC-BD
    // |   |  // le diagonali del rettangolo sono AD-BC
    // C---D  //TxIR inferiori orizzontali

    int prev2_medianX;
    int prev2_medianY;
 
    /*
    float height_left = (568 * CamToMouseMult) - (200 * CamToMouseMult);
    float height_right = (568 * CamToMouseMult) - (200 * CamToMouseMult);
    float width_top = (623 * CamToMouseMult) - (400 * CamToMouseMult);
    float width_bottom = (623 * CamToMouseMult) - (400 * CamToMouseMult);
    */
    float height_left  = height;
    float height_right = height;
    float width_top    = width;
    float width_bottom = width;


    /////////////////////////////////////////////////////////////////////

    // --- VARIABILI MOTORE KINEMATICO (MOLLA) ---
    // Questi array conservano il "debito" di scostamento tra la geometria pura calcolata
    // e la posizione reale mostrata all'utente. Devono essere float per consentire
    // assorbimenti frazionari (sub-pixel) senza incappare nell'arrotondamento prematuro a zero.
    float offset_X[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float offset_Y[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    // Definisce la reattività del filtro anti-glitch. Un valore di 0.25 significa che
    // ogni errore ottico istantaneo viene assorbito al 25% per frame, distribuendo 
    // l'anomalia su ~4 frame. Mantiene il feeling "snappy" ma uccide il tremolio.
    float FPS_NORMALIZATION = 1.0f;
    float COSTANTE_MOLLA = 0.25f; 
    
    // Questo è un muro matematico. Quando mancano dei punti ottici, l'algoritmo valuta 6 combinazioni.
    // PENALITA_STORICA garantisce che il codice preferisca mantenere l'associazione logica 
    // del frame precedente piuttosto che inventare rotazioni spazialmente assurde.
    const int32_t PENALITA_STORICA = 1000000; 

    int prev_GeomX[4] = {0, 0, 0, 0};
    int prev_GeomY[4] = {0, 0, 0, 0};

    /////////////////////////////////////////////////////////////////////

    int camMaxX = 0;
    uint8_t camToMouseShift = 0;
    int camToMouseMult = 1;
    int mouseResX = 0;
    int mouseResY = 0;

    bool calibrationMode = false;    // true = è in corso la calibrazione
    bool wideLayout = false;         // true = configurazione con base maggiore di altezza ad esempio per sensori sugli angoli del monitor

};

#endif // _OpenFIRE_Square_Advanced_h_
#endif //USE_SQUARE_ADVANCED