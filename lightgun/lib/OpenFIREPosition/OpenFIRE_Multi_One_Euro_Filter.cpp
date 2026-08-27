#ifdef USE_MULTI_ONE_EURO_FILTER
/*!
 * @file OpenFIRE_Multi_One_Euro_Filter.cpp
 * @brief Library for One Euro Filter for 4 LED
 * @n CPP Library for One Euro Filter for 4 LED
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V2.0
 * @date 2026
 */

#include "OpenFIRE_Multi_One_Euro_Filter.h"

OpenFIRE_One_Euro_Multi::OpenFIRE_One_Euro_Multi() {
    lastMicros = 0;
    initialized = false;
    for(int i = 0; i < 4; i++) {
        states[i] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    }
}

void OpenFIRE_One_Euro_Multi::configure(const CameraProfile& profile) {
    sensor_noise_unit_x = ((float)profile.mouseResX / (float)profile.sensorResX) * profile.noiseFactor;
    sensor_noise_unit_y = ((float)profile.mouseResY / (float)profile.sensorResY) * profile.noiseFactor;

    snap_base_x = snap_sensor_gain * sensor_noise_unit_x;
    snap_base_y = snap_sensor_gain * sensor_noise_unit_y;
    snap_edge_multiplier_x = snap_edge_sensor_gain * sensor_noise_unit_x;
    snap_edge_multiplier_y = snap_edge_sensor_gain * sensor_noise_unit_y;
    beta_base_x = beta_sensor_gain / sensor_noise_unit_x;
    beta_base_y = beta_sensor_gain / sensor_noise_unit_y;

    mouseMaxX = profile.mouseMaxX;
    mouseMaxY = profile.mouseMaxY;
    inv_center_x = 1.0f / ((float)mouseMaxX * 0.5f);
    inv_center_y = 1.0f / ((float)mouseMaxY * 0.5f);
    fallbackDt = 1.0f / (float)profile.fps;

    initialized = false;
    lastMicros = 0;
}

void OpenFIRE_One_Euro_Multi::process(int* x_in, int* y_in, float* x_out, float* y_out) {
    // Calcolo del Delta-Time (dt) tramite il micro-clock dell'ESP32.
    // L'aritmetica unsigned assorbe nativamente l'overflow del timer.
    unsigned long currentMicros = micros();
    float dt = ((float)(currentMicros - lastMicros)) * 0.000001f; 
    lastMicros = currentMicros;

    // Protezione base per divisioni per zero o blocchi anomali
    // if (dt <= 0.0f) dt = 0.005f; 
    if (dt <= 0.0f) dt = fallbackDt;

    // Cold start dinamico: se passa troppo tempo (es. perdita del tracking),
    // forziamo il reset per evitare di calcolare velocità inerziali spaziali.
    if (dt > 0.1f) {
        initialized = false;
    }

    // Inizializzazione al primissimo frame valido
    if (!initialized) {
        for (int i = 0; i < 4; i++) {
            states[i].x_prev = states[i].x_hat = (float)x_in[i];
            states[i].y_prev = states[i].y_hat = (float)y_in[i];
            states[i].vel_x_hat = states[i].vel_y_hat = 0.0f;
            x_out[i] = states[i].x_hat;
            y_out[i] = states[i].y_hat;
        }
        initialized = true;
        return;
    }

    const float dt_two_pi = dt * OEF_TWO_PI;
    
    // --- GESTIONE ASIMMETRICA DELLA DERIVATA ---
    const float a_d_base = fast_alpha(d_cutoff_base, dt_two_pi);  
    const float a_d_snap = fast_alpha(d_cutoff_snap, dt_two_pi); 
    
    const float inv_dt = 1.0f / dt; 

    // --- Trova l'Anello Forte (Il LED più pulito guida la reattività) ---
    float best_edge_attenuation = 0.20f;

    for (int i = 0; i < 4; i++) {
        // Valutiamo la penalita ottica solo per coordinate interne al FOV.
        // Le coordinate stimate fuori limite non subiscono penalizzazione da bordo/lente.
        if (x_in[i] >= 0 && x_in[i] <= mouseMaxX && y_in[i] >= 0 && y_in[i] <= mouseMaxY) {
            float distH = fabsf((float)x_in[i] - ((float)mouseMaxX * 0.5f)) * inv_center_x;
            float distV = fabsf((float)y_in[i] - ((float)mouseMaxY * 0.5f)) * inv_center_y;
            
            float maxDistSq = (distH * distH) + (distV * distV);
            float edge_attenuation = 1.0f - (maxDistSq * 0.8f); // 0.8f lascia un margine vitale agli angoli estremi
            
            if (edge_attenuation < 0.20f) edge_attenuation = 0.20f; 
            
            if (edge_attenuation > best_edge_attenuation) {
                best_edge_attenuation = edge_attenuation;
            }
        }
    }

    // Beta adattivo e soglie dinamiche, scalati separatamente sui due assi.
    // DFRobot mantiene esattamente il tuning originale; le altre CAM vengono
    // adattate in base alla granularita fisica del sensore e a CamNoiseFactor.
    float adaptiveBetaX = beta_base_x * best_edge_attenuation;
    float adaptiveBetaY = beta_base_y * best_edge_attenuation;

    float edge_factor = 1.0f - best_edge_attenuation;
    float dynamic_snap_x = snap_base_x + (snap_edge_multiplier_x * edge_factor);
    float dynamic_snap_y = snap_base_y + (snap_edge_multiplier_y * edge_factor);

    float lower_bound_x = dynamic_snap_x * 0.8f;
    float upper_bound_x = dynamic_snap_x * 1.2f;
    float inv_lerp_range_x = 1.0f / (dynamic_snap_x * 0.4f);

    float lower_bound_y = dynamic_snap_y * 0.8f;
    float upper_bound_y = dynamic_snap_y * 1.2f;
    float inv_lerp_range_y = 1.0f / (dynamic_snap_y * 0.4f);

    // --- Trova la Massima Energia Cinetica (Reattività Globale) ---
    float max_abs_vel_x = 0.0f;
    float max_abs_vel_y = 0.0f;

    for (int i = 0; i < 4; i++) {
        // Aggiorniamo le derivate per l'Asse X
        float dx = ((float)x_in[i] - states[i].x_prev) * inv_dt;
        float diff_x = fabsf(dx - states[i].vel_x_hat);
        
        float a_d_current_x;
        if (diff_x <= lower_bound_x) {
            a_d_current_x = a_d_base;
        } else if (diff_x >= upper_bound_x) {
            a_d_current_x = a_d_snap;
        } else {
            a_d_current_x = a_d_base + ((diff_x - lower_bound_x) * inv_lerp_range_x) * (a_d_snap - a_d_base);
        }
        
        states[i].vel_x_hat += a_d_current_x * (dx - states[i].vel_x_hat);
        
        float abs_vel_x = fabsf(states[i].vel_x_hat);
        if (abs_vel_x > max_abs_vel_x) max_abs_vel_x = abs_vel_x;

        // Aggiorniamo le derivate per l'Asse Y
        float dy = ((float)y_in[i] - states[i].y_prev) * inv_dt;
        float diff_y = fabsf(dy - states[i].vel_y_hat);
        
        float a_d_current_y;
        if (diff_y <= lower_bound_y) {
            a_d_current_y = a_d_base;
        } else if (diff_y >= upper_bound_y) {
            a_d_current_y = a_d_snap;
        } else {
            a_d_current_y = a_d_base + ((diff_y - lower_bound_y) * inv_lerp_range_y) * (a_d_snap - a_d_base);
        }
        
        states[i].vel_y_hat += a_d_current_y * (dy - states[i].vel_y_hat);
        
        float abs_vel_y = fabsf(states[i].vel_y_hat);
        if (abs_vel_y > max_abs_vel_y) max_abs_vel_y = abs_vel_y;
    }

    // --- Calcolo dell'Alpha Unificato per il Corpo Rigido ---
    float cutoff_x = min_cutoff + adaptiveBetaX * max_abs_vel_x;
    if (cutoff_x > max_cutoff) cutoff_x = max_cutoff; 
    float a_x = fast_alpha(cutoff_x, dt_two_pi);

    float cutoff_y = min_cutoff + adaptiveBetaY * max_abs_vel_y;
    if (cutoff_y > max_cutoff) cutoff_y = max_cutoff;
    float a_y = fast_alpha(cutoff_y, dt_two_pi);

    // --- Filtraggio all'Unisono ---
    for (int i = 0; i < 4; i++) {
        
        // Il Micro-Snap deve sapere solo se il raw e rimasto identico al frame precedente.
        // Evitiamo di ricalcolare due derivate gia non necessarie in questa fase.
        bool x_still = ((float)x_in[i] == states[i].x_prev);
        bool y_still = ((float)y_in[i] == states[i].y_prev);

        states[i].x_hat += a_x * ((float)x_in[i] - states[i].x_hat);
        states[i].x_prev = (float)x_in[i];

        states[i].y_hat += a_y * ((float)y_in[i] - states[i].y_hat);
        states[i].y_prev = (float)y_in[i];

        // --- MICRO-SNAP (TAGLIO DELLA CODA ASINTOTICA BLINDATO) ---
        if (x_still && fabsf((float)x_in[i] - states[i].x_hat) < micro_snap_threshold) {
            states[i].x_hat = (float)x_in[i];
        }
        if (y_still && fabsf((float)y_in[i] - states[i].y_hat) < micro_snap_threshold) {
            states[i].y_hat = (float)y_in[i];
        }

        // --- SAFETY NET HARDWARE ---
        if (isnan(states[i].x_hat) || isnan(states[i].y_hat)) {
            states[i].x_hat = (float)x_in[i];
            states[i].y_hat = (float)y_in[i];
            states[i].vel_x_hat = 0.0f;
            states[i].vel_y_hat = 0.0f;
        }

        // Output sub-pixel diretto nell'array float del main
        x_out[i] = states[i].x_hat;
        y_out[i] = states[i].y_hat;
    }
}

#endif // USE_MULTI_ONE_EURO_FILTER
