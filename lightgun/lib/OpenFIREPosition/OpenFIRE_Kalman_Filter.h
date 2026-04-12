#ifdef USE_POS_KALMAN_FILTER
/*!
 * @file OpenFIRE_Kalman_Filter.h
 * @brief Filtro Kalman per singolo punto – dichiarazione classe
 * @n CPP file for Filtro Kalman
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2025
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2025
 */

// OpenFIRE_Kalman_Filter.h
// Definizione della classe filtro Kalman per singolo punto

#ifndef OPENFIRE_KALMAN_FILTER_H
#define OPENFIRE_KALMAN_FILTER_H

#include <stdint.h>
#include "OpenFIREConst.h"


class OpenFIRE_Kalman_Filter {
public:
    // Costruttore: inizializza lo stato interno
    OpenFIRE_Kalman_Filter();

    // Applica il filtro Kalman al punto (in/out)
    void Kalman_Filter(int &outX, int &outY);

private:
    //================================================================
    // COSTANTI DI TUNING FILTRO KALMAN (Pos/Vel) – versione singolo punto
    //================================================================
    static constexpr float kf_X_MIN              = 0.0f;
    static constexpr float kf_X_MAX              = static_cast<float>(MouseMaxX * 3);
    static constexpr float kf_Y_MIN              = 0.0f;
    static constexpr float kf_Y_MAX              = static_cast<float>(MouseMaxY * 3);

    static constexpr float kf_X_CENTER           = (kf_X_MAX + kf_X_MIN) / 2.0f;
    static constexpr float kf_Y_CENTER           = (kf_Y_MAX + kf_Y_MIN) / 2.0f;
    static constexpr float kf_HALF_WIDTH         = (kf_X_MAX - kf_X_MIN) / 2.0f;
    static constexpr float kf_HALF_HEIGHT        = (kf_Y_MAX - kf_Y_MIN) / 2.0f;

    static constexpr float kf_INITIAL_P_POS_VALUE = 100.0f;
    static constexpr float kf_INITIAL_P_VEL_VALUE = 10.0f;

    static constexpr float kf_MIN_COVARIANCE_VALUE = 1e-6f;
    static constexpr float kf_MAX_P_VALUE          = 1e6f;

    static constexpr float kf_Q_MIN_PROCESS        = 0.2f;
    static constexpr float kf_Q_MAX_PROCESS        = 8.0f;
    static constexpr float kf_ACCEL_Q_THRESH_START = 5.0f;
    static constexpr float kf_ACCEL_Q_THRESH_END   = 50.0f;

    static constexpr float kf_R_MIN                = 0.05f;
    static constexpr float kf_R_MAX                = 1000.0f;
    static constexpr float kf_ACCEL_R_THRESH_START = 10.0f;
    static constexpr float kf_ACCEL_R_THRESH_END   = 100.0f;

    static constexpr float kf_R_X_EDGE_SMOOTH_START = 0.4f;
    static constexpr float kf_R_X_EDGE_SMOOTH_END   = 1.0f;
    static constexpr float kf_R_Y_EDGE_SMOOTH_START = 0.4f;
    static constexpr float kf_R_Y_EDGE_SMOOTH_END   = 1.0f;

    static constexpr float kf_R_AT_X_EDGE          = 2000.0f;
    static constexpr float kf_R_AT_Y_EDGE          = 2000.0f;
    static constexpr float kf_R_AT_X_EDGE_FOR_Y    = 500.0f;
    static constexpr float kf_R_AT_Y_EDGE_FOR_X    = 500.0f;

    static constexpr float kf_R_CROSS_AXIS_INFLUENCE_X = 0.2f;
    static constexpr float kf_R_CROSS_AXIS_INFLUENCE_Y = 0.2f;

    //================================================================
    // STATO INTERNO FILTRO KALMAN (Pos/Vel) – singolo punto
    //================================================================
    float kf_x_state[2];       // [posizione, velocità]
    float kf_y_state[2];       // [posizione, velocità]

    float kf_p_x_00;
    float kf_p_x_01;
    float kf_p_x_10;
    float kf_p_x_11;

    float kf_p_y_00;
    float kf_p_y_01;
    float kf_p_y_10;
    float kf_p_y_11;

    float kf_last_measured_x;
    float kf_last_measured_y;
    float kf_last_vx_raw;
    float kf_last_vy_raw;

    float kf_inv_halfW;
    float kf_inv_halfH;
    float kf_inv_dQ;
    float kf_inv_dR;

    bool  kf_is_initialized;
};

#endif // OPENFIRE_KALMAN_FILTER_H
#endif //USE_POS_KALMAN_FILTER