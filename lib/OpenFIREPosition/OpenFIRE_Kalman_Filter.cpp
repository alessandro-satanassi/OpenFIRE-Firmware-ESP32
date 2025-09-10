#ifdef USE_POS_KALMAN_FILTER
/*!
 * @file OpenFIRE_Kalman_Filter.cpp
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


// versione ottimizzata
#include <Arduino.h>
#include "OpenFIRE_Kalman_Filter.h"

// Helper inlined
static inline float cnsf(float v, float lo, float hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

static inline float ssf(float e0, float e1, float x) {
    x = cnsf((x - e0) / (e1 - e0), 0.0f, 1.0f);
    return x * x * (3.0f - 2.0f * x);
}

static inline float lerpf(float a, float b, float t) {
    return a + (b - a) * t;
}

OpenFIRE_Kalman_Filter::OpenFIRE_Kalman_Filter()
{
    // Stato iniziale
    kf_x_state[0] = kf_x_state[1] = 0.0f;
    kf_y_state[0] = kf_y_state[1] = 0.0f;

    kf_p_x_00 = kf_INITIAL_P_POS_VALUE;
    kf_p_x_01 = 0.0f;
    kf_p_x_10 = 0.0f;
    kf_p_x_11 = kf_INITIAL_P_VEL_VALUE;

    kf_p_y_00 = kf_INITIAL_P_POS_VALUE;
    kf_p_y_01 = 0.0f;
    kf_p_y_10 = 0.0f;
    kf_p_y_11 = kf_INITIAL_P_VEL_VALUE;

    kf_last_measured_x = kf_last_measured_y = 0.0f;
    kf_last_vx_raw = kf_last_vy_raw = 0.0f;

    kf_inv_halfW = 1.0f / kf_HALF_WIDTH;
    kf_inv_halfH = 1.0f / kf_HALF_HEIGHT;
    kf_inv_dQ    = 1.0f / (kf_ACCEL_Q_THRESH_END - kf_ACCEL_Q_THRESH_START);
    kf_inv_dR    = 1.0f / (kf_ACCEL_R_THRESH_END - kf_ACCEL_R_THRESH_START);

    kf_is_initialized = false;
}

void OpenFIRE_Kalman_Filter::Kalman_Filter(int &outX, int &outY)
{
    const float mx = static_cast<float>(outX + MouseMaxX);
    const float my = static_cast<float>(outY + MouseMaxY);

    if (!kf_is_initialized) {
        kf_x_state[0] = mx;  kf_x_state[1] = 0.0f;
        kf_y_state[0] = my;  kf_y_state[1] = 0.0f;
        kf_last_measured_x = mx;
        kf_last_measured_y = my;
        kf_last_vx_raw = kf_last_vy_raw = 0.0f;
        kf_is_initialized = true;
        return;
    }

    // Velocità e accelerazione grezze
    const float vx = mx - kf_last_measured_x;
    const float vy = my - kf_last_measured_y;
    const float ax = vx - kf_last_vx_raw;
    const float ay = vy - kf_last_vy_raw;
    const float mag = fabsf(vx) + fabsf(vy) + fabsf(ax) + fabsf(ay);

    // Modulazione Q
    const float tQ = cnsf((mag - kf_ACCEL_Q_THRESH_START) * kf_inv_dQ, 0.0f, 1.0f);
    const float Q  = lerpf(kf_Q_MIN_PROCESS, kf_Q_MAX_PROCESS, tQ * tQ);

    // Modulazione R base
    const float tR = cnsf((mag - kf_ACCEL_R_THRESH_START) * kf_inv_dR, 0.0f, 1.0f);
    const float Rb = lerpf(kf_R_MAX, kf_R_MIN, tR * tR * tR);

    // Influenza bordi
    const float dx = (kf_x_state[0] - kf_X_CENTER) * kf_inv_halfW;
    const float dy = (kf_y_state[0] - kf_Y_CENTER) * kf_inv_halfH;
    const float fx = ssf(kf_R_X_EDGE_SMOOTH_START, kf_R_X_EDGE_SMOOTH_END, fabsf(dx));
    const float fy = ssf(kf_R_Y_EDGE_SMOOTH_START, kf_R_Y_EDGE_SMOOTH_END, fabsf(dy));

    const float Rx = cnsf(
        lerpf(Rb, kf_R_AT_X_EDGE, fx) +
        lerpf(0.0f, kf_R_AT_Y_EDGE_FOR_X, fy * kf_R_CROSS_AXIS_INFLUENCE_X),
        kf_MIN_COVARIANCE_VALUE, kf_R_MAX * 5.0f
    );
    const float Ry = cnsf(
        lerpf(Rb, kf_R_AT_Y_EDGE, fy) +
        lerpf(0.0f, kf_R_AT_X_EDGE_FOR_Y, fx * kf_R_CROSS_AXIS_INFLUENCE_Y),
        kf_MIN_COVARIANCE_VALUE, kf_R_MAX * 5.0f
    );

    // Precalcolo Q matrix terms
    const float q00 = Q * 0.25f;
    const float q01 = Q * 0.5f;
    const float q11 = Q;

    // --- X ---
    kf_x_state[0] += kf_x_state[1];
    float p00 = kf_p_x_00 + kf_p_x_01 + kf_p_x_10 + kf_p_x_11 + q00;
    float p01 = kf_p_x_01 + kf_p_x_11 + q01;
    float p11 = kf_p_x_11 + q11;
    p00 = cnsf(p00, kf_MIN_COVARIANCE_VALUE, kf_MAX_P_VALUE);
    p11 = cnsf(p11, kf_MIN_COVARIANCE_VALUE, kf_MAX_P_VALUE);
    const float denomX = p00 + Rx;
    const float K0 = p00 / denomX;
    const float K1 = p01 / denomX;
    const float inx = mx - kf_x_state[0];
    kf_x_state[0] += K0 * inx;
    kf_x_state[1] += K1 * inx;
    kf_p_x_00 = p00 - K0 * p00;
    kf_p_x_01 = p01 - K0 * p01;
    kf_p_x_10 = kf_p_x_01;
    kf_p_x_11 = p11 - K1 * p01;

    // --- Y ---
    kf_y_state[0] += kf_y_state[1];
    float py00 = kf_p_y_00 + kf_p_y_01 + kf_p_y_10 + kf_p_y_11 + q00;
    float py01 = kf_p_y_01 + kf_p_y_11 + q01;
    float py11 = kf_p_y_11 + q11;
    py00 = cnsf(py00, kf_MIN_COVARIANCE_VALUE, kf_MAX_P_VALUE);
    py11 = cnsf(py11, kf_MIN_COVARIANCE_VALUE, kf_MAX_P_VALUE);
    const float denomY = py00 + Ry;
    const float L0 = py00 / denomY;
    const float L1 = py01 / denomY;
    const float iny = my - kf_y_state[0];
    kf_y_state[0] += L0 * iny;
    kf_y_state[1] += L1 * iny;
    kf_p_y_00 = py00 - L0 * py00;
    kf_p_y_01 = py01 - L0 * py01;
    kf_p_y_10 = kf_p_y_01;
    kf_p_y_11 = py11 - L1 * py01;

    // Salva stato
    kf_last_measured_x = mx;
    kf_last_vx_raw     = vx;
    kf_last_measured_y = my;
    kf_last_vy_raw     = vy;

    // Output
    outX = static_cast<int>(kf_x_state[0]) - MouseMaxX;
    outY = static_cast<int>(kf_y_state[0]) - MouseMaxY;
}



 // versione originale
#ifdef COMMENTO
#include <Arduino.h>
#include "OpenFIRE_Kalman_Filter.h"


// Costruttore: inizializza lo stato interno
OpenFIRE_Kalman_Filter::OpenFIRE_Kalman_Filter()
{
    // Stato
    kf_x_state[0] = 0.0f;
    kf_x_state[1] = 0.0f;
    kf_y_state[0] = 0.0f;
    kf_y_state[1] = 0.0f;

    kf_p_x_00 = kf_INITIAL_P_POS_VALUE;
    kf_p_x_01 = 0.0f;
    kf_p_x_10 = 0.0f;
    kf_p_x_11 = kf_INITIAL_P_VEL_VALUE;

    kf_p_y_00 = kf_INITIAL_P_POS_VALUE;
    kf_p_y_01 = 0.0f;
    kf_p_y_10 = 0.0f;
    kf_p_y_11 = kf_INITIAL_P_VEL_VALUE;

    kf_last_measured_x = 0.0f;
    kf_last_measured_y = 0.0f;
    kf_last_vx_raw     = 0.0f;
    kf_last_vy_raw     = 0.0f;

    kf_is_initialized  = false;
}

// Implementazione filtro Kalman per un singolo punto
void OpenFIRE_Kalman_Filter::Kalman_Filter(int &outX, int &outY)
{
    float mx = static_cast<float>(outX+MouseMaxX);
    float my = static_cast<float>(outY+MouseMaxY);

    // Parametri pre-calcolati
    const float inv_halfW = 1.0f / kf_HALF_WIDTH;
    const float inv_halfH = 1.0f / kf_HALF_HEIGHT;
    const float inv_dQ    = 1.0f / (kf_ACCEL_Q_THRESH_END - kf_ACCEL_Q_THRESH_START);
    const float inv_dR    = 1.0f / (kf_ACCEL_R_THRESH_END - kf_ACCEL_R_THRESH_START);

    // Inizializzazione al primo frame
    if (!kf_is_initialized) {
        kf_x_state[0] = mx;
        kf_x_state[1] = 0.0f;
        kf_y_state[0] = my;
        kf_y_state[1] = 0.0f;
        //outX = static_cast<int>(mx-MouseMaxX);
        //outY = static_cast<int>(my-MouseMaxY);
        kf_last_measured_x = mx;
        kf_last_measured_y = my;
        kf_last_vx_raw     = 0.0f;
        kf_last_vy_raw     = 0.0f;
        kf_is_initialized  = true;
        return;
    }

    // Funzioni di supporto
    auto cns = [](float v, float lo, float hi) {
        return v < lo ? lo : (v > hi ? hi : v);
    };
    auto ss = [&](float e0, float e1, float x) {
        x = cns((x - e0) / (e1 - e0), 0.0f, 1.0f);
        return x * x * (3.0f - 2.0f * x);
    };
    auto lerp = [](float a, float b, float t) {
        return a + (b - a) * t;
    };

    // 1) Velocità e accelerazione grezze
    float vx = mx - kf_last_measured_x;
    float vy = my - kf_last_measured_y;
    float ax = vx - kf_last_vx_raw;
    float ay = vy - kf_last_vy_raw;
    float mag = fabsf(vx) + fabsf(vy) + fabsf(ax) + fabsf(ay);

    // 2) Modulazione Q
    float tQ = cns((mag - kf_ACCEL_Q_THRESH_START) * inv_dQ, 0.0f, 1.0f);
    float Q  = lerp(kf_Q_MIN_PROCESS, kf_Q_MAX_PROCESS, tQ * tQ);

    // 3) Modulazione R
    float tR = cns((mag - kf_ACCEL_R_THRESH_START) * inv_dR, 0.0f, 1.0f);
    float Rb = lerp(kf_R_MAX, kf_R_MIN, tR * tR * tR);

    // 4) Influenza bordi
    float dx = (kf_x_state[0] - kf_X_CENTER) * inv_halfW;
    float dy = (kf_y_state[0] - kf_Y_CENTER) * inv_halfH;
    float fx = ss(kf_R_X_EDGE_SMOOTH_START, kf_R_X_EDGE_SMOOTH_END, fabsf(dx));
    float fy = ss(kf_R_Y_EDGE_SMOOTH_START, kf_R_Y_EDGE_SMOOTH_END, fabsf(dy));
    float Rx = cns(lerp(Rb, kf_R_AT_X_EDGE, fx)
                 + lerp(0.0f, kf_R_AT_Y_EDGE_FOR_X, fy * kf_R_CROSS_AXIS_INFLUENCE_X),
                 kf_MIN_COVARIANCE_VALUE, kf_R_MAX * 5.0f);
    float Ry = cns(lerp(Rb, kf_R_AT_Y_EDGE, fy)
                 + lerp(0.0f, kf_R_AT_X_EDGE_FOR_Y, fx * kf_R_CROSS_AXIS_INFLUENCE_Y),
                 kf_MIN_COVARIANCE_VALUE, kf_R_MAX * 5.0f);

    // 5) Predizione & aggiornamento X
    kf_x_state[0] += kf_x_state[1];
    float q00 = Q * 0.25f, q01 = Q * 0.5f, q11 = Q;
    float p00 = kf_p_x_00 + kf_p_x_01 + kf_p_x_10 + kf_p_x_11 + q00;
    float p01 = kf_p_x_01 + kf_p_x_11 + q01;
    float p11 = kf_p_x_11 + q11;
    p00 = cns(p00, kf_MIN_COVARIANCE_VALUE, kf_MAX_P_VALUE);
    p11 = cns(p11, kf_MIN_COVARIANCE_VALUE, kf_MAX_P_VALUE);
    float denomX = p00 + Rx;
    float K0 = p00 / denomX, K1 = p01 / denomX;
    float inx = mx - kf_x_state[0];
    kf_x_state[0] += K0 * inx;
    kf_x_state[1] += K1 * inx;
    kf_p_x_00 = p00 - K0 * p00;
    kf_p_x_01 = p01 - K0 * p01;
    kf_p_x_10 = kf_p_x_01;
    kf_p_x_11 = p11 - K1 * p01;

    // 6) Predizione & aggiornamento Y
    kf_y_state[0] += kf_y_state[1];
    float py00 = kf_p_y_00 + kf_p_y_01 + kf_p_y_10 + kf_p_y_11 + q00;
    float py01 = kf_p_y_01 + kf_p_y_11 + q01;
    float py11 = kf_p_y_11 + q11;
    py00 = cns(py00, kf_MIN_COVARIANCE_VALUE, kf_MAX_P_VALUE);
    py11 = cns(py11, kf_MIN_COVARIANCE_VALUE, kf_MAX_P_VALUE);
    float denomY = py00 + Ry;
    float L0 = py00 / denomY, L1 = py01 / denomY;
    float iny = my - kf_y_state[0];
    kf_y_state[0] += L0 * iny;
    kf_y_state[1] += L1 * iny;
    kf_p_y_00 = py00 - L0 * py00;
    kf_p_y_01 = py01 - L0 * py01;
    kf_p_y_10 = kf_p_y_01;
    kf_p_y_11 = py11 - L1 * py01;

    // 7) Salvataggio stato precedente
    kf_last_measured_x = mx;
    kf_last_vx_raw     = vx;
    kf_last_measured_y = my;
    kf_last_vy_raw     = vy;

    // 8) Output filtrato
    outX = static_cast<int>(kf_x_state[0]) - MouseMaxX;
    outY = static_cast<int>(kf_y_state[0]) - MouseMaxY;
}
#endif COMMENTO

#endif //USE_POS_KALMAN_FILTER