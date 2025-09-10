#ifdef USE_POS_ONE_EURO_FILTER
/*!
 * @file OpenFIRE_One_Euro_Filter.h
 * @brief Filtro OneEuro per singolo punto – dichiarazione classe
 * @n CPP file for Filtro OneEuro
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2025
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2025
 */

// OpenFIRE_One_Euro_Filter.h
// Definizione della classe filtro Kalman per singolo punto

#ifndef OPENFIRE_ONE_EURO_FILTER_H
#define OPENFIRE_ONE_EURO_FILTER_H

#include <stdint.h>
#include "OpenFIREConst.h"


class OpenFIRE_One_Euro_Filter {
public:
    // Costruttore: inizializza lo stato interno
    OpenFIRE_One_Euro_Filter();

    // Applica il filtro Kalman al punto (in/out)
    void One_Euro_Filter(int &outX, int &outY);

private:
    //================================================================
    // COSTANTI DI TUNING FILTRO ONE EURO
    //================================================================
    static constexpr float OE_FREQ           = 200.0f;           // Frequenza campionamento (Hz)
    static constexpr float OE_MIN_CUTOFF     = 1.0f;             // Cutoff minimo (Hz)
    static constexpr float OE_BETA           = 0.007f;           // Guadagno velocità
    static constexpr float OE_D_CUTOFF       = 1.0f;             // Cutoff derivata (Hz)
    static constexpr float OE_DT             = 1.0f / OE_FREQ;   // Intervallo campionamento (s)

    //================================================================
    // STATO INTERNO FILTRO ONE EURO (Pos/Vel) – singolo punto
    //================================================================
    // Stato X
    float oe_x_prev;        // Ultima misura raw X
    float oe_x_hat;         // Ultimo valore filtrato X
    float oe_vel_hat_x;     // Velocità stimata filtrata X

    // Stato Y
    float oe_y_prev;        // Ultima misura raw Y
    float oe_y_hat;         // Ultimo valore filtrato Y
    float oe_vel_hat_y;     // Velocità stimata filtrata Y    

};

#endif // OPENFIRE_ONE_EURO_FILTER_H
#endif //USE_POS_ONE_EURO_FILTER