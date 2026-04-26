
#ifdef USE_MULTI_ONE_EURO_FILTER

#ifndef OpenFIRE_Multi_One_Euro_Filter_h
#define OpenFIRE_Multi_One_Euro_Filter_h

#include <Arduino.h>
#include "OpenFIREConst.h"

class OpenFIRE_One_Euro_Multi {
private:
    struct FilterState {
        float x_prev, x_hat;
        float y_prev, y_hat;
        float vel_x_hat, vel_y_hat;
    };

    FilterState states[4];
    unsigned long lastMicros;
    bool initialized = false;

    // Variabili per il centro dello schermo (calcolate una sola volta al boot)
    float inv_center_x;
    float inv_center_y;

    // --- PARAMETRI DI TUNING ---
    const float min_cutoff = 1.0f; 
    const float d_cutoff = 10.0f;   
    const float max_cutoff = 30.0f; 
    const float beta_base = (0.011f * (float)CamResX) / (float)MouseResX;

    const float OEF_TWO_PI = 6.28318530718f;

    inline float fast_alpha(float cutoff, float dt_two_pi) {
        float te = cutoff * dt_two_pi;
        return te / (te + 1.0f);
    }

public:
    OpenFIRE_One_Euro_Multi();
    void process(int* x, int* y);
};

#endif // OpenFIRE_Multi_One_Euro_Filter_h

#endif // USE_MULTI_ONE_EURO_FILTER

