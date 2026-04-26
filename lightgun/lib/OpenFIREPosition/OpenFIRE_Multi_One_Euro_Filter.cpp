#ifdef USE_MULTI_ONE_EURO_FILTER

#include "OpenFIRE_Multi_One_Euro_Filter.h"

OpenFIRE_One_Euro_Multi::OpenFIRE_One_Euro_Multi() {
    lastMicros = 0;
    initialized = false;
    for(int i = 0; i < 4; i++) {
        states[i] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    }
    
    // PRE-CALCOLO ASSOLUTO: Eseguito una sola volta all'avvio.
    // Troviamo i reciproci del centro schermo per usare moltiplicazioni iper-veloci.
    inv_center_x = 1.0f / ((float)MouseMaxX * 0.5f);
    inv_center_y = 1.0f / ((float)MouseMaxY * 0.5f);
}

void OpenFIRE_One_Euro_Multi::process(int* x, int* y) {
    unsigned long currentMicros = micros();
    
    float dt = ((float)(currentMicros - lastMicros)) * 0.000001f; 
    if (dt <= 0.0f || dt > 0.1f) dt = 0.005f; 
    lastMicros = currentMicros;

    if (!initialized) {
        for (int i = 0; i < 4; i++) {
            states[i].x_prev = states[i].x_hat = (float)x[i];
            states[i].y_prev = states[i].y_hat = (float)y[i];
            states[i].vel_x_hat = states[i].vel_y_hat = 0.0f;
        }
        initialized = true;
        return;
    }

    const float dt_two_pi = dt * OEF_TWO_PI;
    const float a_d = fast_alpha(d_cutoff, dt_two_pi);
    const float inv_dt = 1.0f / dt; 

    for (int i = 0; i < 4; i++) {
        // --- 1. LOGICA SPAZIALE ---
        // Sfruttiamo i centri calcolati nel costruttore
        float distH = fabsf((float)x[i] - ((float)MouseMaxX * 0.5f)) * inv_center_x;
        float distV = fabsf((float)y[i] - ((float)MouseMaxY * 0.5f)) * inv_center_y;
        float maxDist = (distH > distV) ? distH : distV;
        
        if (maxDist > 1.0f) maxDist = 1.0f; 
        float adaptiveBeta = beta_base * (1.0f - (maxDist * 0.10f));

        // --- 2. FILTRAGGIO ASSE X ---
        float dx = ((float)x[i] - states[i].x_prev) * inv_dt;
        
        // OTTIMIZZAZIONE EMA: Forma matematica ridotta (risparmia CPU)
        states[i].vel_x_hat += a_d * (dx - states[i].vel_x_hat);
        
        float cutoff_x = min_cutoff + adaptiveBeta * fabsf(states[i].vel_x_hat);
        if (cutoff_x > max_cutoff) cutoff_x = max_cutoff; 
        
        float a_x = fast_alpha(cutoff_x, dt_two_pi);
        
        // OTTIMIZZAZIONE EMA: Forma matematica ridotta
        states[i].x_hat += a_x * ((float)x[i] - states[i].x_hat);
        states[i].x_prev = (float)x[i];

        // --- 3. FILTRAGGIO ASSE Y ---
        float dy = ((float)y[i] - states[i].y_prev) * inv_dt;
        
        states[i].vel_y_hat += a_d * (dy - states[i].vel_y_hat);
        
        float cutoff_y = min_cutoff + adaptiveBeta * fabsf(states[i].vel_y_hat);
        if (cutoff_y > max_cutoff) cutoff_y = max_cutoff;
        
        float a_y = fast_alpha(cutoff_y, dt_two_pi);
        
        states[i].y_hat += a_y * ((float)y[i] - states[i].y_hat);
        states[i].y_prev = (float)y[i];

        // --- 4. SAFETY NET HARDWARE ---
        // Se un glitch cosmico genera un Not-A-Number, lo resettiamo istantaneamente
        if (isnan(states[i].x_hat) || isnan(states[i].y_hat)) {
            states[i].x_hat = (float)x[i];
            states[i].y_hat = (float)y[i];
            states[i].vel_x_hat = 0.0f;
            states[i].vel_y_hat = 0.0f;
        }

        x[i] = (int)roundf(states[i].x_hat);
        y[i] = (int)roundf(states[i].y_hat);
    }
}

#endif // USE_MULTI_ONE_EURO_FILTER

