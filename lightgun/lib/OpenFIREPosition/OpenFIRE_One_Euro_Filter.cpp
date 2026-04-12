#ifdef USE_POS_ONE_EURO_FILTER
/*!
 * @file OpenFIRE_One_Euro_Filter.cpp
 * @brief Filtro OneEuro per singolo punto – implementazione
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2025
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2025
 */

//versione ottimizzata al massimo, senza funzioni
#include <cmath>
#include "OpenFIRE_One_Euro_Filter.h"

OpenFIRE_One_Euro_Filter::OpenFIRE_One_Euro_Filter()
    : oe_x_prev(0.0f), oe_x_hat(0.0f), oe_vel_hat_x(0.0f),
      oe_y_prev(0.0f), oe_y_hat(0.0f), oe_vel_hat_y(0.0f)
{
}

void OpenFIRE_One_Euro_Filter::One_Euro_Filter(int &outX, int &outY)
{
    // Misure grezze
    const float rawX = static_cast<float>(outX);
    const float rawY = static_cast<float>(outY);

    // Precalcolo 2π * dt (costante comune)
    const float two_pi_dt = 2.0f * 3.14159265358979323846f * OE_DT;

    // Alpha per cutoff di derivata
    float r_d = OE_D_CUTOFF * two_pi_dt;
    const float alpha_d = r_d / (r_d + 1.0f);

    // ---------- X ----------
    float dx = (rawX - oe_x_prev) / OE_DT;
    oe_vel_hat_x = alpha_d * dx + (1.0f - alpha_d) * oe_vel_hat_x;

    float r_x = (OE_MIN_CUTOFF + OE_BETA * fabsf(oe_vel_hat_x)) * two_pi_dt;
    float alpha_x = r_x / (r_x + 1.0f);

    oe_x_hat = alpha_x * rawX + (1.0f - alpha_x) * oe_x_hat;
    oe_x_prev = rawX;

    // ---------- Y ----------
    float dy = (rawY - oe_y_prev) / OE_DT;
    oe_vel_hat_y = alpha_d * dy + (1.0f - alpha_d) * oe_vel_hat_y;

    float r_y = (OE_MIN_CUTOFF + OE_BETA * fabsf(oe_vel_hat_y)) * two_pi_dt;
    float alpha_y = r_y / (r_y + 1.0f);

    oe_y_hat = alpha_y * rawY + (1.0f - alpha_y) * oe_y_hat;
    oe_y_prev = rawY;

    // Output finali
    outX = static_cast<int>(oe_x_hat);
    outY = static_cast<int>(oe_y_hat);
}


 // versione ottimizzata
#ifdef COMMENTO
#include <cmath>
#include "OpenFIRE_One_Euro_Filter.h"

// Coefficiente alpha ottimizzato (precalcolo di 2π * dt)
static inline float oe_alpha_fast(float cutoff, float two_pi_dt) {
    float r = cutoff * two_pi_dt;
    return r / (r + 1.0f);
}

OpenFIRE_One_Euro_Filter::OpenFIRE_One_Euro_Filter()
    : oe_x_prev(0.0f), oe_x_hat(0.0f), oe_vel_hat_x(0.0f),
      oe_y_prev(0.0f), oe_y_hat(0.0f), oe_vel_hat_y(0.0f)
{
}

void OpenFIRE_One_Euro_Filter::One_Euro_Filter(int &outX, int &outY)
{
    const float rawX = static_cast<float>(outX);
    const float rawY = static_cast<float>(outY);

    // Precalcolo costanti comuni
    const float two_pi_dt = 2.0f * 3.14159265358979323846f * OE_DT;
    const float alpha_d = oe_alpha_fast(OE_D_CUTOFF, two_pi_dt);

    // ----- Asse X -----
    float dx = (rawX - oe_x_prev) / OE_DT;
    oe_vel_hat_x = alpha_d * dx + (1.0f - alpha_d) * oe_vel_hat_x;
    float cutoffX = OE_MIN_CUTOFF + OE_BETA * fabsf(oe_vel_hat_x);
    float ax = oe_alpha_fast(cutoffX, two_pi_dt);
    oe_x_hat = ax * rawX + (1.0f - ax) * oe_x_hat;
    oe_x_prev = rawX;

    // ----- Asse Y -----
    float dy = (rawY - oe_y_prev) / OE_DT;
    oe_vel_hat_y = alpha_d * dy + (1.0f - alpha_d) * oe_vel_hat_y;
    float cutoffY = OE_MIN_CUTOFF + OE_BETA * fabsf(oe_vel_hat_y);
    float ay = oe_alpha_fast(cutoffY, two_pi_dt);
    oe_y_hat = ay * rawY + (1.0f - ay) * oe_y_hat;
    oe_y_prev = rawY;

    // Output filtrato
    outX = static_cast<int>(oe_x_hat);
    outY = static_cast<int>(oe_y_hat);
}
#endif // commento

// versione originale
#ifdef COMMENTO
#include <cmath>
#include "OpenFIRE_One_Euro_Filter.h"

// Calcola coefficiente alpha da cutoff
static inline float oe_alpha(float cutoff, float dt) {
    const float PI = 3.14159265358979323846f;
    float r = 2.0f * PI * cutoff * dt;
    return r / (r + 1.0f);
}

// Costruttore: inizializza lo stato interno
OpenFIRE_One_Euro_Filter::OpenFIRE_One_Euro_Filter()
    : oe_x_prev(0.0f), oe_x_hat(0.0f), oe_vel_hat_x(0.0f),
      oe_y_prev(0.0f), oe_y_hat(0.0f), oe_vel_hat_y(0.0f)
{
}

// Implementazione filtro One Euro per un singolo punto
void OpenFIRE_One_Euro_Filter::One_Euro_Filter(int &outX, int &outY)
{
    // Misura grezza
    float rawX = static_cast<float>(outX);
    float rawY = static_cast<float>(outY);

    // Filtro asse X
    float dx = (rawX - oe_x_prev) / OE_DT;
    float ad = oe_alpha(OE_D_CUTOFF, OE_DT);
    oe_vel_hat_x = ad * dx + (1.0f - ad) * oe_vel_hat_x;
    float cutoffX = OE_MIN_CUTOFF + OE_BETA * std::fabs(oe_vel_hat_x);
    float ax = oe_alpha(cutoffX, OE_DT);
    oe_x_hat = ax * rawX + (1.0f - ax) * oe_x_hat;
    oe_x_prev = rawX;

    // Filtro asse Y
    float dy = (rawY - oe_y_prev) / OE_DT;
    float ady = oe_alpha(OE_D_CUTOFF, OE_DT);
    oe_vel_hat_y = ady * dy + (1.0f - ady) * oe_vel_hat_y;
    float cutoffY = OE_MIN_CUTOFF + OE_BETA * std::fabs(oe_vel_hat_y);
    float ay = oe_alpha(cutoffY, OE_DT);
    oe_y_hat = ay * rawY + (1.0f - ay) * oe_y_hat;
    oe_y_prev = rawY;

    // Output filtrato
    outX = static_cast<int>(oe_x_hat);
    outY = static_cast<int>(oe_y_hat);
}
#endif // commento


#endif //USE_POS_ONE_EURO_FILTER