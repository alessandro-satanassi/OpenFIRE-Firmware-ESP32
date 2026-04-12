#ifdef USE_SQUARE_ADVANCED
/*!
 * @file OpenFIRE_Square_Advanced.h
 * @brief Light Gun library for 4 LED setup
 * @n CPP file for Samco Light Gun 4 LED setup
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2025
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2025
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


class OpenFIRE_Square {
public:
    //================================================================
    // ENUM PUBBLICO E API
    //================================================================

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

private:
    //================================================================
    // STATO DEL TRACKER
    //================================================================
    int  FinalX[4] = {400 * CamToMouseMult, 623 * CamToMouseMult, 400 * CamToMouseMult, 623 * CamToMouseMult};
    int  FinalY[4] = {200 * CamToMouseMult, 200 * CamToMouseMult, 568 * CamToMouseMult, 568 * CamToMouseMult};
    int  medianX = MouseMaxX / 2;
    int  medianY = MouseMaxY / 2;
    unsigned int see[4];
    float height = (568 * CamToMouseMult) - (200 * CamToMouseMult);
    float width = (623 * CamToMouseMult) - (400 * CamToMouseMult);
    float angle = 0;

    unsigned int start = 0;
    unsigned int seenFlags = 0;

    bool is_tracking_stable = false; // s efalse vul dire che al frame precedente non sono stati visti punti o solo 1
    float ideal_aspect_ratio = 16.0f / 9.0f; // rapporto tra larghezza e altezza del rettangolo formato dai TxIR
    uint8_t prev_num_points_seen = 0; // numero di TxIR visti
    uint8_t prev_point_seen_mask = 0b00000000; // A->00001000  B->00000100  C->00000010  D->00000001
    uint8_t last_identified_idx; //Per ricordare quale punto era visibile con solo un sensore visto
    uint8_t current_point_seen_mask = 0;
    // A---B  //TxIR superiori orrizzontali 
    // |   |  // la base AB-CD è sempre più corta dell'altezza AC-BD
    // |   |  // le diagonali del rettangolo sono AD-BC
    // C---D  //TxIR inferiori orrizzontali

    int prev2_medianX;
    int prev2_medianY;
    float prev2_height;
    float prev2_width;
    float prev2_angle;

    float shape_vec_AB_x, shape_vec_AB_y; // Vettore dal vertice A al vertice B
    float shape_vec_AD_x, shape_vec_AD_y; // Vettore dal vertice A al vertice D
    int point4_X, point4_Y;
    uint8_t point4_idx;

    //int shape_correction_x;
    //int shape_correction_y;
    float height_left = (568 * CamToMouseMult) - (200 * CamToMouseMult);
    float height_right = (568 * CamToMouseMult) - (200 * CamToMouseMult);
    float width_top = (623 * CamToMouseMult) - (400 * CamToMouseMult);
    float width_bottom = (623 * CamToMouseMult) - (400 * CamToMouseMult);

};

#endif // _OpenFIRE_Square_Advanced_h_
#endif //USE_SQUARE_ADVANCED