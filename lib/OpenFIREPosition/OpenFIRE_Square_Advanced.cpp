#ifdef USE_SQUARE_ADVANCED
/*!
 * @file OpenFIRE_Square_Advanced.cpp
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

#include <Arduino.h>
#include "OpenFIRE_Square_Advanced.h"

#include <cfloat>


/*
 * ==============================================================================
 * Sezione Helper: Funzioni di utilità per l'ottimizzazione
 *
 * Definire queste funzioni all'inizio del file o in un file header condiviso.
 * L'uso di 'inline' suggerisce al compilatore di sostituire la chiamata alla
 * funzione con il suo codice, eliminando l'overhead e massimizzando la velocità.
 * ==============================================================================
 */

/**
 * @brief Calcola la radice quadrata inversa (1/sqrt(x)) in modo approssimato ma estremamente veloce.
 * @param number Il numero di cui calcolare la radice quadrata inversa.
 * @return Un valore approssimato di 1 / sqrt(number).
 */
inline float fast_inv_sqrtf(float number) {
    union {
        float f;
        uint32_t i;
    } conv = { .f = number };
    conv.i = 0x5f3759df - (conv.i >> 1); // Stima iniziale con "magic number"
    // Iterazione di Newton-Raphson per raffinare la stima e migliorare la precisione
    conv.f *= 1.5f - (number * 0.5f * conv.f * conv.f);
    return conv.f;
}

/**
 * @brief Arrotonda un valore float all'intero più vicino in modo veloce.
 * Funziona correttamente sia con valori positivi che negativi.
 * @param val Il valore float da arrotondare.
 * @return L'intero arrotondato.
 */
inline int fast_roundf(float val) {
    // Aggiunge o sottrae 0.5 e poi tronca (casting a int).
    // Esempio: fast_roundf(-3.7) -> (int)(-3.7 - 0.5) -> (int)(-4.2) -> -4
    // Esempio: fast_roundf(3.7)  -> (int)(3.7 + 0.5)  -> (int)(4.2)  -> 4
    return (int)(val + (val > 0.0f ? 0.5f : -0.5f));
}



// Definizione degli indici per chiarezza
#define A 0
#define B 1
#define C 2
#define D 3

void OpenFIRE_Square::begin(const int* px, const int* py, unsigned int seen) {

    // --- FASE 0: Inizializzazione e Controllo di Visibilità ---
    
    seenFlags = seen;

    // Il sistema richiede di vedere tutti e 4 i sensori almeno una volta per inizializzarsi.
    if (seenFlags == 0x0F) { // 0x0F in binario è 1111, tutti i sensori visti.
        start = 0xFF;
    } else if (!start) {
        return; // Esci se non siamo ancora stati inizializzati.
    }

    // --- FASE 1: Estrazione e Pre-elaborazione dei Punti Visti ---

    uint8_t num_points_seen = 0;
    int positionXX[4];
    int positionYY[4];

    // Estrae i punti visibili, li mette negli array di lavoro e applica la trasformazione.
    for (uint8_t i = 0; i < 4; ++i) {
        if ((seenFlags >> i) & 0x01) {
            positionXX[num_points_seen] = (CamMaxX - px[i]) << CamToMouseShift;
            positionYY[num_points_seen] = py[i] << CamToMouseShift;
            num_points_seen++;
        }
    }

    // Procede solo se abbiamo abbastanza dati per lavorare (almeno 2 punti).
    //////if (num_points_seen > 1) {
    if (num_points_seen >= 1) {

        //just_reacquired_tracking = true; // Abbiamo appena ritrovato i punti
    
        // Dichiarazione di tutte le variabili locali usate nei vari blocchi logici.
        // Questo risolve gli errori di "scope" (variabile non dichiarata).
        int32_t dx, dy;
        uint8_t a, b, c, d;
        uint8_t P1, P2;


        if (num_points_seen == 1)
        {
            // L'algoritmo procede solo se ha una posizione di riferimento dal frame precedente.
            if (prev_num_points_seen >= 1) 
            {
                // --- FASE 1: IDENTIFICAZIONE DEL VERTICE PIÙ VICINO ---

                const int x1 = positionXX[0]; 
                const int y1 = positionYY[0];

                uint8_t identified_idx = 0;
                int32_t min_dist_sq = INT32_MAX; 

                for (uint8_t i = 0; i < 4; i++) 
                {
                    const int dx = x1 - FinalX[i];
                    const int dy = y1 - FinalY[i];
            
                    // MODIFICA: La distanza al quadrato è esplicitamente int32_t
                    const int32_t dist_sq = (int32_t)dx * dx + (int32_t)dy * dy;

                    if (dist_sq < min_dist_sq) 
                    {
                        min_dist_sq = dist_sq;
                        identified_idx = i;
                    }
                }

                // --- CONTROLLO DI SICUREZZA (ANTI-JUMP) ---
        
                const int max_jump_distance = ((int)width * 3) / 4; 
            
                const int32_t MAX_ALLOWED_DISTANCE_SQ = (int32_t)max_jump_distance * max_jump_distance;

                if (min_dist_sq > MAX_ALLOWED_DISTANCE_SQ) {
                    prev_num_points_seen = 0;
                    return; 
                }

                // --- FASE 2: TRASCINAMENTO RIGIDO (Rigid Drag) ---
        
                const int delta_x = x1 - FinalX[identified_idx];
                const int delta_y = y1 - FinalY[identified_idx];

                for (uint8_t i = 0; i < 4; i++) {
                    //FinalX[i] += delta_x;
                    //FinalY[i] += delta_y;
                    //positionXX[i] = FinalX[i];
                    //positionYY[i] = FinalY[i];
                    positionXX[i] = FinalX[i] + delta_x;
                    positionYY[i] = FinalY[i] + delta_y;


                }
            }
            else // Se non c'era nessun punto nel frame precedente, non possiamo fare nulla.
            {
                prev_num_points_seen = 0;
                return;
            }
        }


        #ifdef COMMENTO
        // il migliore fino al 22/06/2025 ore 19:44
        // =================================================================================================
        // GESTIONE 2 PUNTI VISTI (VERSIONE FINALE UNIFICATA CON PESI A 3 LIVELLI)
        //
        // Strategia:
        // La struttura di calcolo è unica. L'intelligenza è demandata alla selezione
        // dei pesi (W_POS, W_SHAPE, W_ANGLE) che vengono impostati a zero per i costi
        // che si vogliono escludere in un determinato contesto.
        // =================================================================================================

        if (num_points_seen == 2)
        {
            // --- FASE 1: IDENTIFICAZIONE ---

            const int x1 = positionXX[0], y1 = positionYY[0];
            const int x2 = positionXX[1], y2 = positionYY[1];

            uint8_t best_idx1 = 0, best_idx2 = 1;
            float min_total_cost = -1.0f;

            // --- Calcolo di TUTTI i possibili componenti di costo ---

            // 1. Costo di Posizione
            const uint8_t pairs[6][2] = { {0,1}, {2,3}, {0,2}, {1,3}, {0,3}, {1,2} };
            float pos_costs[12];
            for (int i = 0; i < 6; i++) {
                const uint8_t v1 = pairs[i][0], v2 = pairs[i][1];
                float dx1a=(float)x1-FinalX[v1], dy1a=(float)y1-FinalY[v1];
                float dx2a=(float)x2-FinalX[v2], dy2a=(float)y2-FinalY[v2];
                pos_costs[i*2] = dx1a*dx1a + dy1a*dy1a + dx2a*dx2a + dy2a*dy2a;
                float dx1b=(float)x1-FinalX[v2], dy1b=(float)y1-FinalY[v2];
                float dx2b=(float)x2-FinalX[v1], dy2b=(float)y2-FinalY[v1];
                pos_costs[i*2 + 1] = dx1b*dx1b + dy1b*dy1b + dx2b*dx2b + dy2b*dy2b;
            }
            const float inv_pos_norm = 1.0f / (width*width + height*height + 1e-6f);

            // 2. Costo di Forma e Angolo (calcolati anche se poi non usati)
            const float d2H = height * height;
            const float estW = height * ideal_aspect_ratio;
            const float d2W = estW * estW;
            const float d2D = d2W + d2H;
            const float expected_d2s[6] = { d2W, d2W, d2H, d2H, d2D, d2D };
            const float dxm = (float)x1 - (float)x2;
            const float dym = (float)y1 - (float)y2;
            const float d2m = dxm*dxm + dym*dym;
            const float r = ideal_aspect_ratio;
            const float r2 = r * r;
            const float cos_theta_W = (1.0f - r2) / (1.0f + r2 + 1e-6f);
            const float cos_theta_H = (r2 - 1.0f) / (1.0f + r2 + 1e-6f);
            const float cdx1 = (float)x1 - medianX, cdy1 = (float)y1 - medianY;
            const float cdx2 = (float)x2 - medianX, cdy2 = (float)y2 - medianY;
            const float dot = cdx1 * cdx2 + cdy1 * cdy2;
            float mag = sqrtf((cdx1*cdx1 + cdy1*cdy1) * (cdx2*cdx2 + cdy2*cdy2));
            mag = fmaxf(mag, 1e-6f);
            float cos_theta = dot / mag;
            cos_theta = fminf(1.0f, fmaxf(-1.0f, cos_theta));

            // =========================================================================================
            // MODELLO DEFINITIVO v2 - Bilanciamento Finale
            // Sintesi tra stabilità temporale e robustezza geometrica.
            // =========================================================================================

            float W_POS, W_SHAPE, W_ANGLE;

            if (prev_num_points_seen == 4) {
                // LIVELLO 1 (4/4): GOLD. Posizione e Forma sono ugualmente importanti.
                W_POS   = 1.2f;
                W_SHAPE = 1.2f;  // <- Alzato per eguagliare la posizione.
                W_ANGLE = 1.0f;

            } else if (prev_num_points_seen == 3) {
                // LIVELLO 2 (3/4): SILVER. Leggermente meno fiducia nella posizione,
                // ma fiducia massima nella forma "ideale" della memoria.
                W_POS   = 1.0f;
                W_SHAPE = 1.4f;  // <- Leggermente più alto di W_POS per la perfezione della forma.
                W_ANGLE = 1.1f;

            } else if (prev_num_points_seen == 2) {
                // LIVELLO 3 (2/4): BRONZE. Meno fiducia nella posizione (basata su 2pt),
                // più fiducia nella geometria misurata ORA.
                W_POS   = 0.8f;
                W_SHAPE = 1.3f;  // <- Alto per compensare la fiducia ridotta in W_POS.
                W_ANGLE = 1.2f;

            } else { // prev_num_points_seen <= 1
                // LIVELLO 4 (1/4 o 0/4): INAFFIDABILE. La posizione non conta quasi nulla.
                // Ci si affida solo alla geometria dell'input corrente.
                W_POS   = 0.1f;
                W_SHAPE = 1.8f;
                W_ANGLE = 0.1f; //05f; //0.1f; //1.5f;
            }


            // --- Calcolo del costo totale con un'unica logica ---
            for (int i = 0; i < 6; i++) {
                //const float shape_cost = fabsf(d2m - expected_d2s[i]) / fmaxf(expected_d2s[i], 25.0f);
                const float shape_cost = fabsf(d2m - expected_d2s[i]) / fmaxf(d2m, fmaxf(expected_d2s[i], 1.0f));

        
                const float expected_cos = (i < 2) ? cos_theta_W : (i < 4) ? cos_theta_H : -1.0f;
                const float diff_cos = cos_theta - expected_cos;
                const float angle_cost = (diff_cos * diff_cos) / 4.0f;
                const float normalized_pos_cost1 = pos_costs[i*2] * inv_pos_norm;
                const float normalized_pos_cost2 = pos_costs[i*2 + 1] * inv_pos_norm;

                const float cost1 = normalized_pos_cost1 * W_POS + shape_cost * W_SHAPE + angle_cost * W_ANGLE;
                const float cost2 = normalized_pos_cost2 * W_POS + shape_cost * W_SHAPE + angle_cost * W_ANGLE;

                if (min_total_cost < 0 || cost1 < min_total_cost) {
                    min_total_cost = cost1;
                    best_idx1 = pairs[i][0]; best_idx2 = pairs[i][1];
                }
                if (min_total_cost < 0 || cost2 < min_total_cost) {
                    min_total_cost = cost2;
                    best_idx1 = pairs[i][1]; best_idx2 = pairs[i][0];
                }
            }

            // --- FASE 2: RICOSTRUZIONE GEOMETRICA ---

            float Vx[4], Vy[4];
            for (int i = 0; i < 4; i++) { Vx[i] = FinalX[i]; Vy[i] = FinalY[i]; }
            Vx[best_idx1] = (float)x1; Vy[best_idx1] = (float)y1;
            Vx[best_idx2] = (float)x2; Vy[best_idx2] = (float)y2;
    
            const uint8_t idx1 = best_idx1, idx2 = best_idx2;
            const bool top    = (idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0);
            const bool bottom = (idx1 == 2 && idx2 == 3) || (idx1 == 3 && idx2 == 2);
            const bool left   = (idx1 == 0 && idx2 == 2) || (idx1 == 2 && idx2 == 0);
            const bool right  = (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 1);

            if (top || bottom) {
                uint8_t p1 = top ? 0 : 2, p2 = top ? 1 : 3;
                uint8_t p3 = top ? 2 : 0, p4 = top ? 3 : 1;
                if (Vx[p1] > Vx[p2]) { uint8_t t=p1; p1=p2; p2=t; }
                float dx = Vx[p2]-Vx[p1], dy = Vy[p2]-Vy[p1];
                float base = sqrtf(dx*dx + dy*dy);
                if (base > 2.0f) {
                    float cos_angle = fabsf(dx) / base;
                    base = base / fmaxf(cos_angle, 0.3f);
                    float newH = fmaxf(base / ideal_aspect_ratio, height * 0.5f);
                    float inv_base = 1.0f / base;
                    float ndx = -dy*inv_base, ndy = dx*inv_base;
                    float sx = ndx*newH, sy = ndy*newH;
                    if (top) {
                        Vx[p3]=Vx[p1]+sx; Vy[p3]=Vy[p1]+sy; Vx[p4]=Vx[p2]+sx; Vy[p4]=Vy[p2]+sy;
                    } else {
                        Vx[p3]=Vx[p1]-sx; Vy[p3]=Vy[p1]-sy; Vx[p4]=Vx[p2]-sx; Vy[p4]=Vy[p2]-sy;
                    }
                }
            } else if (left || right) {
                uint8_t p1 = left ? 0 : 1, p2 = left ? 2 : 3;
                uint8_t p3 = left ? 1 : 0, p4 = left ? 3 : 2;
                if (Vy[p1] > Vy[p2]) { uint8_t t=p1; p1=p2; p2=t; }
                float dx = Vx[p2]-Vx[p1], dy = Vy[p2]-Vy[p1];
                float hgt = sqrtf(dx*dx+dy*dy);
                if (hgt > 2.0f) {
                    float newW = hgt * ideal_aspect_ratio;
                    float inv_hgt = 1.0f/hgt;
                    float ndx = -dy*inv_hgt, ndy = dx*inv_hgt;
                    float sx = ndx*newW, sy = ndy*newW;
                    if (left) {
                        Vx[p3]=Vx[p1]-sx; Vy[p3]=Vy[p1]-sy; Vx[p4]=Vx[p2]-sx; Vy[p4]=Vy[p2]-sy;
                    } else {
                        Vx[p3]=Vx[p1]+sx; Vy[p3]=Vy[p1]+sy; Vx[p4]=Vx[p2]+sx; Vy[p4]=Vy[p2]+sy;
                    }
                }
            } else { // Diagonale
                bool diagAD = (idx1 == 0 && idx2 == 3) || (idx1 == 3 && idx2 == 0);
                float px1 = diagAD ? Vx[0] : Vx[1], py1 = diagAD ? Vy[0] : Vy[1];
                float px2 = diagAD ? Vx[3] : Vx[2], py2 = diagAD ? Vy[3] : Vy[2];
                float dx = px2 - px1, dy = py2 - py1;
                float L = sqrtf(dx*dx + dy*dy);
                if (L > 2.0f) {
                    float H = L / sqrtf(ideal_aspect_ratio*ideal_aspect_ratio + 1.0f);
                    float W = H * ideal_aspect_ratio;
                    float phi_m = atan2f(dy, dx);
                    float phi_i = diagAD ? atan2f(H, W) : atan2f(H, -W);
                    float theta = phi_m - phi_i;
                    float c=cosf(theta), s=sinf(theta);
                    float wx=c*W, wy=s*W;
                    float hx=-s*H, hy=c*H;
                    if (diagAD) {
                        Vx[1]=Vx[0]+wx; Vy[1]=Vy[0]+wy; Vx[2]=Vx[0]+hx; Vy[2]=Vy[0]+hy;
                    } else {
                        Vx[0]=Vx[1]-wx; Vy[0]=Vy[1]-wy; Vx[3]=Vx[1]+hx; Vy[3]=Vy[1]+hy;
                    }
                }
            }

            // --- FASE 3: ASSEGNAZIONE FINALE ---

            for (int i = 0; i < 4; i++) {
                FinalX[i] = roundf(Vx[i]);
                FinalY[i] = roundf(Vy[i]);
                positionXX[i] = FinalX[i];
                positionYY[i] = FinalY[i];
            }
            current_point_seen_mask = (1 << (3 - best_idx1)) | (1 << (3 - best_idx2));
        }
        #endif



        #ifdef COMMENTO
        /*
         * BLOCCO DI CODICE COMPLETO PER GESTIONE DI 2 PUNTI RILEVATI
         * Ottimizzato per microcontrollori (ESP32, RP2040) e senza dipendenze da struct o funzioni esterne.
         *
         * PRESUPPONE CHE LE SEGUENTI VARIABILI ESISTANO E SIANO POPOLATE NEL SUO SCOPE:
         * - Input: positionXX[0], positionXX[1], positionYY[0], positionYY[1]
         * - Stato t-1: FinalX[4], FinalY[4], medianX, medianY, height, width, angle, prev_num_points_seen
         * - Stato t-2: prev2_medianX, prev2_medianY
         * - Config: ideal_aspect_ratio
         * - Output: Il codice aggiorna FinalX, FinalY, positionXX, positionYY e calcola current_point_seen_mask
        */
        if (num_points_seen == 2)
        {
            
            /////////////// NUOVO PESO VELOCITA ////////////////////////////////////////////////
            // --- FASE 0: PREDIZIONE CINEMATICA (IL NUOVO "ARBITRO") ---
            float velocita_x = 0.0f;
            float velocita_y = 0.0f;
    
            // Calcola la velocità solo se abbiamo una storia valida (evita salti al primo avvio)
            if (prev2_medianX != 0 || prev2_medianY != 0) {
                velocita_x = (float)medianX - (float)prev2_medianX;
                velocita_y = (float)medianY - (float)prev2_medianY;
            }
            const float centroide_predetto_x = (float)medianX + velocita_x;
            const float centroide_predetto_y = (float)medianY + velocita_y;
            /////////////////// FINE PER NUOVO PESO VELOCITA /////////////////////////////////////
            
            
            // --- FASE 1: IDENTIFICAZIONE E CALCOLO COSTI ---
            const int x1 = positionXX[0], y1 = positionYY[0];
            const int x2 = positionXX[1], y2 = positionYY[1];

            uint8_t best_idx1 = 0, best_idx2 = 1;
            float min_total_cost = -1.0f;

            // Calcolo di TUTTI i possibili componenti di costo
            // 1. Costo di Posizione
            const uint8_t pairs[6][2] = { {0,1}, {2,3}, {0,2}, {1,3}, {0,3}, {1,2} };
            int32_t pos_costs[12];
            for (uint8_t i = 0; i < 6; i++) {
                const uint8_t v1 = pairs[i][0], v2 = pairs[i][1];
                int32_t dx1a=(float)x1-FinalX[v1], dy1a=(float)y1-FinalY[v1];
                int32_t dx2a=(float)x2-FinalX[v2], dy2a=(float)y2-FinalY[v2];
                pos_costs[i*2] = dx1a*dx1a + dy1a*dy1a + dx2a*dx2a + dy2a*dy2a;
                int32_t dx1b=(float)x1-FinalX[v2], dy1b=(float)y1-FinalY[v2];
                int32_t dx2b=(float)x2-FinalX[v1], dy2b=(float)y2-FinalY[v1];
                pos_costs[i*2 + 1] = dx1b*dx1b + dy1b*dy1b + dx2b*dx2b + dy2b*dy2b;
            }
            const float inv_pos_norm = 1.0f / (width*width + height*height + 1e-6f);

            // 2. Costo di Forma e Angolo
            const float d2H = height * height;
            const float estW = height * ideal_aspect_ratio;
            const float d2W = estW * estW;
            const float d2D = d2W + d2H;
            const float expected_d2s[6] = { d2W, d2W, d2H, d2H, d2D, d2D };
            const float dxm = (float)x1 - (float)x2;
            const float dym = (float)y1 - (float)y2;
            const float d2m = dxm*dxm + dym*dym;
            const float r = ideal_aspect_ratio;
            const float r2 = r * r;
            const float cos_theta_W = (1.0f - r2) / (1.0f + r2 + 1e-6f);
            const float cos_theta_H = (r2 - 1.0f) / (1.0f + r2 + 1e-6f);
            const float cdx1 = (float)x1 - medianX, cdy1 = (float)y1 - medianY;
            const float cdx2 = (float)x2 - medianX, cdy2 = (float)y2 - medianY;
            const float dot = cdx1 * cdx2 + cdy1 * cdy2;
            float mag = sqrtf((cdx1*cdx1 + cdy1*cdy1) * (cdx2*cdx2 + cdy2*cdy2));
            mag = fmaxf(mag, 1e-6f);
            float cos_theta = dot / mag;
            cos_theta = fminf(1.0f, fmaxf(-1.0f, cos_theta));

            // 3. Pesi dei Costi basati sulla fiducia dello stato precedente
            float W_POS, W_SHAPE, W_ANGLE;
            if (prev_num_points_seen == 4) {
                W_POS = 1.2f; W_SHAPE = 1.2f; W_ANGLE = 1.0f;
            } else if (prev_num_points_seen == 3) {
                W_POS = 1.0f; W_SHAPE = 1.4f; W_ANGLE = 1.1f;
            } else if (prev_num_points_seen == 2) {
                W_POS = 0.8f; W_SHAPE = 1.3f; W_ANGLE = 1.2f;
            } else {
                W_POS = 0.1f; W_SHAPE = 1.8f; W_ANGLE = 0.1f;
            }
    
            //////////////////////////////////////////////////////////////////////
            // NUOVO: Peso per il costo cinematico, da usare come "tie-breaker"
            const float W_KINEMATIC = 0.7f;
            //////////////////////////////////////////////////////////////////////

            // --- Ciclo Principale: Valuta tutte le 12 ipotesi ---
            for (int i = 0; i < 6; i++) {
                ////////////////////////////////////////////////////////////////////////////////////////
                // NUOVO: Calcolo del costo di coerenza cinematica
                const float ipotetico_centroide_x = ((float)x1 + (float)x2) / 2.0f;
                const float ipotetico_centroide_y = ((float)y1 + (float)y2) / 2.0f;
                const float dx_kin = ipotetico_centroide_x - centroide_predetto_x;
                const float dy_kin = ipotetico_centroide_y - centroide_predetto_y;
                const float kinematic_cost = (dx_kin * dx_kin + dy_kin * dy_kin) * inv_pos_norm;
                ////////////////////////////////////////////////////////////////////////////////////////


                // Calcolo costi esistenti
                const float shape_cost = fabsf(d2m - expected_d2s[i]) / fmaxf(d2m, fmaxf(expected_d2s[i], 1.0f));
                const float expected_cos = (i < 2) ? cos_theta_W : (i < 4) ? cos_theta_H : -1.0f;
                const float diff_cos = cos_theta - expected_cos;
                const float angle_cost = (diff_cos * diff_cos) / 4.0f;
                const float normalized_pos_cost1 = pos_costs[i*2] * inv_pos_norm;
                const float normalized_pos_cost2 = pos_costs[i*2 + 1] * inv_pos_norm;

                // MODIFICATO: Costo totale con il nuovo termine cinematico
                const float cost1 = normalized_pos_cost1 * W_POS + shape_cost * W_SHAPE + angle_cost * W_ANGLE + kinematic_cost * W_KINEMATIC;
                const float cost2 = normalized_pos_cost2 * W_POS + shape_cost * W_SHAPE + angle_cost * W_ANGLE + kinematic_cost * W_KINEMATIC;

                if (min_total_cost < 0 || cost1 < min_total_cost) {
                    min_total_cost = cost1;
                    best_idx1 = pairs[i][0]; best_idx2 = pairs[i][1];
                }
                if (min_total_cost < 0 || cost2 < min_total_cost) {
                    min_total_cost = cost2;
                    best_idx1 = pairs[i][1]; best_idx2 = pairs[i][0];
                }
            }

            // --- FASE 2: RICOSTRUZIONE GEOMETRICA ---
            float Vx[4], Vy[4];
            for (int i = 0; i < 4; i++) { Vx[i] = FinalX[i]; Vy[i] = FinalY[i]; }
            Vx[best_idx1] = (float)x1; Vy[best_idx1] = (float)y1;
            Vx[best_idx2] = (float)x2; Vy[best_idx2] = (float)y2;
    
            const uint8_t idx1 = best_idx1, idx2 = best_idx2;
            const bool top      = (idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0);
            const bool bottom   = (idx1 == 2 && idx2 == 3) || (idx1 == 3 && idx2 == 2);
            const bool left     = (idx1 == 0 && idx2 == 2) || (idx1 == 2 && idx2 == 0);
            const bool right    = (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 1);

            if (top || bottom) {
                uint8_t p1 = top ? 0 : 2, p2 = top ? 1 : 3;
                uint8_t p3 = top ? 2 : 0, p4 = top ? 3 : 1;
                if (Vx[p1] > Vx[p2]) { uint8_t t=p1; p1=p2; p2=t; }
                float dx = Vx[p2]-Vx[p1], dy = Vy[p2]-Vy[p1];
                float base = sqrtf(dx*dx + dy*dy);
                if (base > 2.0f) {
                    float cos_angle = fabsf(dx) / base;
                    base = base / fmaxf(cos_angle, 0.3f);
                    float newH = fmaxf(base / ideal_aspect_ratio, height * 0.5f);
                    float inv_base = 1.0f / base;
                    float ndx = -dy*inv_base, ndy = dx*inv_base;
                    float sx = ndx*newH, sy = ndy*newH;
                    if (top) {
                        Vx[p3]=Vx[p1]+sx; Vy[p3]=Vy[p1]+sy; Vx[p4]=Vx[p2]+sx; Vy[p4]=Vy[p2]+sy;
                    } else {
                        Vx[p3]=Vx[p1]-sx; Vy[p3]=Vy[p1]-sy; Vx[p4]=Vx[p2]-sx; Vy[p4]=Vy[p2]-sy;
                    }
                }
            } else if (left || right) {
                uint8_t p1 = left ? 0 : 1, p2 = left ? 2 : 3;
                uint8_t p3 = left ? 1 : 0, p4 = left ? 3 : 2;
                if (Vy[p1] > Vy[p2]) { uint8_t t=p1; p1=p2; p2=t; }
                float dx = Vx[p2]-Vx[p1], dy = Vy[p2]-Vy[p1];
                float hgt = sqrtf(dx*dx+dy*dy);
                if (hgt > 2.0f) {
                    float newW = hgt * ideal_aspect_ratio;
                    float inv_hgt = 1.0f/hgt;
                    float ndx = -dy*inv_hgt, ndy = dx*inv_hgt;
                    float sx = ndx*newW, sy = ndy*newW;
                    if (left) {
                        Vx[p3]=Vx[p1]-sx; Vy[p3]=Vy[p1]-sy; Vx[p4]=Vx[p2]-sx; Vy[p4]=Vy[p2]-sy;
                    } else {
                        Vx[p3]=Vx[p1]+sx; Vy[p3]=Vy[p1]+sy; Vx[p4]=Vx[p2]+sx; Vy[p4]=Vy[p2]+sy;
                    }
                }
            } else { // Diagonale
                bool diagAD = (idx1 == 0 && idx2 == 3) || (idx1 == 3 && idx2 == 0);
                float px1 = diagAD ? Vx[0] : Vx[1], py1 = diagAD ? Vy[0] : Vy[1];
                float px2 = diagAD ? Vx[3] : Vx[2], py2 = diagAD ? Vy[3] : Vy[2];
                float dx = px2 - px1, dy = py2 - py1;
                float L = sqrtf(dx*dx + dy*dy);
                if (L > 2.0f) {
                    float H = L / sqrtf(ideal_aspect_ratio*ideal_aspect_ratio + 1.0f);
                    float W = H * ideal_aspect_ratio;
                    float phi_m = atan2f(dy, dx);
                    float phi_i = diagAD ? atan2f(H, W) : atan2f(H, -W);
                    float theta = phi_m - phi_i;
                    float c=cosf(theta), s=sinf(theta);
                    float wx=c*W, wy=s*W;
                    float hx=-s*H, hy=c*H;
                    if (diagAD) {
                        Vx[1]=Vx[0]+wx; Vy[1]=Vy[0]+wy; Vx[2]=Vx[0]+hx; Vy[2]=Vy[0]+hy;
                    } else {
                        Vx[0]=Vx[1]-wx; Vy[0]=Vy[1]-wy; Vx[3]=Vx[1]+hx; Vy[3]=Vy[1]+hy;
                    }
                }
            }

            // --- FASE 3: ASSEGNAZIONE FINALE ---
            for (int i = 0; i < 4; i++) {
                FinalX[i] = roundf(Vx[i]);
                FinalY[i] = roundf(Vy[i]);
                positionXX[i] = FinalX[i];
                positionYY[i] = FinalY[i];
            }
            current_point_seen_mask = (1 << (3 - best_idx1)) | (1 << (3 - best_idx2));
        }
        #endif

        #ifdef COMMENTO
        /*
         * BLOCCO DI CODICE COMPLETO v7 - RICOSTRUZIONE CON MEMORIA DI FORMA
         * Mantiene la forma non-rettangolare appresa a 4 punti per eliminare
         * gli "scatti" del centroide durante la transizione 4->2 punti.
         * Utilizza la logica di costo della v5 per l'identificazione e una nuova FASE 2 per la ricostruzione.
         *
         * PRESUPPONE CHE LE NUOVE VARIABILI shape_vec... ESISTANO E SIANO STATE POPOLATE
        */
        if (num_points_seen == 2)
        {
            // --- FASE 0 & 0.5: PREDIZIONE CINEMATICA E PESI (Dalla v5) ---
            float velocita_x = 0.0f;
            float velocita_y = 0.0f;
            if (prev2_medianX != 0 || prev2_medianY != 0) {
                velocita_x = (float)medianX - (float)prev2_medianX;
                velocita_y = (float)medianY - (float)prev2_medianY;
            }
            const float centroide_predetto_x = (float)medianX + velocita_x;
            const float centroide_predetto_y = (float)medianY + velocita_y;

            const float speed = hypotf(velocita_x, velocita_y);
            const float object_diagonal = hypotf(width, height);
            const float motion_intensity = speed / (object_diagonal + 1e-6f);
    
            const float W_KINEMATIC_base = 0.7f;
            const float W_VERTICAL_CONSISTENCY_base = 1.5f;
            const float W_KINEMATIC_dynamic = W_KINEMATIC_base * (1.0f + motion_intensity);
            const float W_VERTICAL_CONSISTENCY_dynamic = W_VERTICAL_CONSISTENCY_base * (1.0f + motion_intensity * 2.0f);
            const float W_SCALE_CONSISTENCY = 3.0f;

            // --- FASE 1: IDENTIFICAZIONE DEI VERTICI (Dalla v5) ---
            const int x1 = positionXX[0], y1 = positionYY[0];
            const int x2 = positionXX[1], y2 = positionYY[1];

            uint8_t best_idx1 = 0, best_idx2 = 1;
            float min_total_cost = -1.0f;

            const uint8_t pairs[6][2] = { {0,1}, {2,3}, {0,2}, {1,3}, {0,3}, {1,2} };
            float pos_costs[12];
            for (int i = 0; i < 6; i++) {
                const uint8_t v1 = pairs[i][0], v2 = pairs[i][1];
                float dx1a=(float)x1-FinalX[v1], dy1a=(float)y1-FinalY[v1];
                float dx2a=(float)x2-FinalX[v2], dy2a=(float)y2-FinalY[v2];
                pos_costs[i*2] = dx1a*dx1a + dy1a*dy1a + dx2a*dx2a + dy2a*dy2a;
                float dx1b=(float)x1-FinalX[v2], dy1b=(float)y1-FinalY[v2];
                float dx2b=(float)x2-FinalX[v1], dy2b=(float)y2-FinalY[v1];
                pos_costs[i*2 + 1] = dx1b*dx1b + dy1b*dy1b + dx2b*dx2b + dy2b*dy2b;
            }
            const float inv_pos_norm = 1.0f / (width*width + height*height + 1e-6f);

            const float d2H = height * height;
            const float estW = height * ideal_aspect_ratio;
            const float d2W = estW * estW;
            const float d2D = d2W + d2H;
            const float expected_d2s[6] = { d2W, d2W, d2H, d2H, d2D, d2D };
            const float dxm = (float)x1 - (float)x2;
            const float dym = (float)y1 - (float)y2;
            const float d2m = dxm*dxm + dym*dym;
    
            float W_POS, W_SHAPE, W_ANGLE;
            if (prev_num_points_seen == 4) {
                W_POS = 1.2f; W_SHAPE = 1.2f; W_ANGLE = 1.0f;
            } else if (prev_num_points_seen == 3) {
                W_POS = 1.0f; W_SHAPE = 1.4f; W_ANGLE = 1.1f;
            } else if (prev_num_points_seen == 2) {
                W_POS = 0.8f; W_SHAPE = 1.3f; W_ANGLE = 1.2f;
            } else {
                W_POS = 0.1f; W_SHAPE = 1.8f; W_ANGLE = 0.1f;
            }

            const float r = ideal_aspect_ratio;
            const float r2 = r * r;
            const float cos_theta_W = (1.0f - r2) / (1.0f + r2 + 1e-6f);
            const float cos_theta_H = (r2 - 1.0f) / (1.0f + r2 + 1e-6f);
            const float cdx1 = (float)x1 - medianX, cdy1 = (float)y1 - medianY;
            const float cdx2 = (float)x2 - medianX, cdy2 = (float)y2 - medianY;
            const float dot = cdx1 * cdx2 + cdy1 * cdy2;
            float mag = sqrtf((cdx1*cdx1 + cdy1*cdy1) * (cdx2*cdx2 + cdy2*cdy2));
            mag = fmaxf(mag, 1e-6f);
            float cos_theta = dot / mag;
            cos_theta = fminf(1.0f, fmaxf(-1.0f, cos_theta));

            for (int i = 0; i < 6; i++) {
                const float ipotetico_centroide_x = ((float)x1 + (float)x2) / 2.0f;
                const float ipotetico_centroide_y = ((float)y1 + (float)y2) / 2.0f;
                const float dx_kin = ipotetico_centroide_x - centroide_predetto_x;
                const float dy_kin = ipotetico_centroide_y - centroide_predetto_y;
                const float kinematic_cost = (dx_kin * dx_kin + dy_kin * dy_kin) * inv_pos_norm;

                float vertical_consistency_cost = 0.0f;
                if (i < 2) {
                    const float old_avg_y = (i == 0) ? (FinalY[0] + FinalY[1]) / 2.0f : (FinalY[2] + FinalY[3]) / 2.0f;
                    const float new_avg_y = ((float)y1 + (float)y2) / 2.0f;
                    const float measured_dy = new_avg_y - old_avg_y;
                    vertical_consistency_cost = fabsf(measured_dy - velocita_y) / (height + 10.0f);
                }

                const float expected_length = sqrtf(expected_d2s[i]);
                const float measured_length = sqrtf(d2m);
                const float scaling_consistency_cost = fabsf(measured_length - expected_length) / (object_diagonal + 10.0f);

                const float shape_cost = fabsf(d2m - expected_d2s[i]) / fmaxf(d2m, fmaxf(expected_d2s[i], 1.0f));
                const float expected_cos = (i < 2) ? cos_theta_W : (i < 4) ? cos_theta_H : -1.0f;
                const float diff_cos = cos_theta - expected_cos;
                const float angle_cost = (diff_cos * diff_cos) / 4.0f;
                const float normalized_pos_cost1 = pos_costs[i*2] * inv_pos_norm;
                const float normalized_pos_cost2 = pos_costs[i*2 + 1] * inv_pos_norm;

                const float cost1 = normalized_pos_cost1 * W_POS + shape_cost * W_SHAPE + angle_cost * W_ANGLE + kinematic_cost * W_KINEMATIC_dynamic + vertical_consistency_cost * W_VERTICAL_CONSISTENCY_dynamic + scaling_consistency_cost * W_SCALE_CONSISTENCY;
                const float cost2 = normalized_pos_cost2 * W_POS + shape_cost * W_SHAPE + angle_cost * W_ANGLE + kinematic_cost * W_KINEMATIC_dynamic + vertical_consistency_cost * W_VERTICAL_CONSISTENCY_dynamic + scaling_consistency_cost * W_SCALE_CONSISTENCY;

                if (min_total_cost < 0 || cost1 < min_total_cost) {
                    min_total_cost = cost1;
                    best_idx1 = pairs[i][0]; best_idx2 = pairs[i][1];
                }
                if (min_total_cost < 0 || cost2 < min_total_cost) {
                    min_total_cost = cost2;
                    best_idx1 = pairs[i][1]; best_idx2 = pairs[i][0];
                }
            }


            // --- FASE 2: RICOSTRUZIONE CON MEMORIA DI FORMA (Logica nuova) ---
            float Vx[4], Vy[4];

            // Controlliamo se la coppia migliore identificata è il lato superiore (A-B)
            bool is_top_edge_best_match = (best_idx1 == 0 && best_idx2 == 1) || (best_idx1 == 1 && best_idx2 == 0);

            if (is_top_edge_best_match && (shape_vec_AB_x != 0 || shape_vec_AB_y != 0)) {
                // --- NUOVA LOGICA DI RICOSTRUZIONE TRAMITE TRASFORMAZIONE ---
        
                // Identifichiamo quale punto misurato (p1, p2) è il nuovo A e quale il nuovo B
                float pA_x, pA_y, pB_x, pB_y;
                float dist1_to_A = powf((float)positionXX[0] - FinalX[0], 2) + powf((float)positionYY[0] - FinalY[0], 2);
                float dist2_to_A = powf((float)positionXX[1] - FinalX[0], 2) + powf((float)positionYY[1] - FinalY[0], 2);
        
                if (dist1_to_A < dist2_to_A) {
                    pA_x = (float)positionXX[0]; pA_y = (float)positionYY[0]; // p1 è il nuovo A
                    pB_x = (float)positionXX[1]; pB_y = (float)positionYY[1]; // p2 è il nuovo B
                } else {
                    pA_x = (float)positionXX[1]; pA_y = (float)positionYY[1]; // p2 è il nuovo A
                    pB_x = (float)positionXX[0]; pB_y = (float)positionYY[0]; // p1 è il nuovo B
                }

                // Assegnamo i punti noti
                Vx[0] = pA_x; Vy[0] = pA_y;
                Vx[1] = pB_x; Vy[1] = pB_y;
        
                // Vettore misurato A->B
                const float measured_vec_x = pB_x - pA_x;
                const float measured_vec_y = pB_y - pA_y;
        
                // Vettore memorizzato A->B (dalla forma appresa)
                const float old_vec_x = shape_vec_AB_x;
                const float old_vec_y = shape_vec_AB_y;

                // Calcoliamo la trasformazione (scala + rotazione) per mappare il vecchio vettore sul nuovo
                const float old_len = hypotf(old_vec_x, old_vec_y);
                const float new_len = hypotf(measured_vec_x, measured_vec_y);
                const float scale = (old_len > 1e-6f) ? new_len / old_len : 1.0f;
        
                const float old_ang = atan2f(old_vec_y, old_vec_x);
                const float new_ang = atan2f(measured_vec_y, measured_vec_x);
                const float delta_ang = new_ang - old_ang;
                const float cos_delta = cosf(delta_ang);
                const float sin_delta = sinf(delta_ang);

                // Applichiamo la stessa trasformazione all'altro vettore memorizzato (A->D) per trovarne la nuova versione
                // 1. Applica la scala
                const float scaled_AD_x = shape_vec_AD_x * scale;
                const float scaled_AD_y = shape_vec_AD_y * scale;
                // 2. Applica la rotazione
                const float new_AD_x = scaled_AD_x * cos_delta - scaled_AD_y * sin_delta;
                const float new_AD_y = scaled_AD_x * sin_delta + scaled_AD_y * cos_delta;
        
                // Ricostruisci i punti mancanti basandoti sulla forma trasformata
                Vx[3] = pA_x + new_AD_x; // D_new = A_new + transformed_vector_AD
                Vy[3] = pA_y + new_AD_y;
                Vx[2] = Vx[3] - measured_vec_x; // C_new = D_new - measured_vector_AB
                Vy[2] = Vy[3] - measured_vec_y;

            } else {
                // --- FALLBACK ALLA VECCHIA RICOSTRUZIONE GEOMETRICA per tutti gli altri casi ---
                for (int i = 0; i < 4; i++) { Vx[i] = FinalX[i]; Vy[i] = FinalY[i]; }
                Vx[best_idx1] = (float)positionXX[0]; Vy[best_idx1] = (float)positionYY[0];
                Vx[best_idx2] = (float)positionXX[1]; Vy[best_idx2] = (float)positionYY[1];

                const uint8_t idx1 = best_idx1, idx2 = best_idx2;
                const bool top      = (idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0);
                const bool bottom   = (idx1 == 2 && idx2 == 3) || (idx1 == 3 && idx2 == 2);
                const bool left     = (idx1 == 0 && idx2 == 2) || (idx1 == 2 && idx2 == 0);
                const bool right    = (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 1);

                if (top || bottom) {
                    uint8_t p1 = top ? 0 : 2, p2 = top ? 1 : 3;
                    uint8_t p3 = top ? 2 : 0, p4 = top ? 3 : 1;
                    if (Vx[p1] > Vx[p2]) { uint8_t t=p1; p1=p2; p2=t; }
                    float dx = Vx[p2]-Vx[p1], dy = Vy[p2]-Vy[p1];
                    float base = sqrtf(dx*dx + dy*dy);
                    if (base > 2.0f) {
                        float cos_angle = fabsf(dx) / base;
                        base = base / fmaxf(cos_angle, 0.3f);
                        float newH = fmaxf(base / ideal_aspect_ratio, height * 0.5f);
                        float inv_base = 1.0f / base;
                        float ndx = -dy*inv_base, ndy = dx*inv_base;
                        float sx = ndx*newH, sy = ndy*newH;
                        if (top) { Vx[p3]=Vx[p1]+sx; Vy[p3]=Vy[p1]+sy; Vx[p4]=Vx[p2]+sx; Vy[p4]=Vy[p2]+sy; }
                        else { Vx[p3]=Vx[p1]-sx; Vy[p3]=Vy[p1]-sy; Vx[p4]=Vx[p2]-sx; Vy[p4]=Vy[p2]-sy; }
                    }
                } else if (left || right) {
                    uint8_t p1 = left ? 0 : 1, p2 = left ? 2 : 3;
                    uint8_t p3 = left ? 1 : 0, p4 = left ? 3 : 2;
                    if (Vy[p1] > Vy[p2]) { uint8_t t=p1; p1=p2; p2=t; }
                    float dx = Vx[p2]-Vx[p1], dy = Vy[p2]-Vy[p1];
                    float hgt = sqrtf(dx*dx+dy*dy);
                    if (hgt > 2.0f) {
                        float newW = hgt * ideal_aspect_ratio;
                        float inv_hgt = 1.0f/hgt;
                        float ndx = -dy*inv_hgt, ndy = dx*inv_hgt;
                        float sx = ndx*newW, sy = ndy*newW;
                        if (left) { Vx[p3]=Vx[p1]-sx; Vy[p3]=Vy[p1]-sy; Vx[p4]=Vx[p2]-sx; Vy[p4]=Vy[p2]-sy; }
                        else { Vx[p3]=Vx[p1]+sx; Vy[p3]=Vy[p1]+sy; Vx[p4]=Vx[p2]+sx; Vy[p4]=Vy[p2]+sy; }
                    }
                } else { // Diagonale
                    bool diagAD = (idx1 == 0 && idx2 == 3) || (idx1 == 3 && idx2 == 0);
                    float px1 = diagAD ? Vx[0] : Vx[1], py1 = diagAD ? Vy[0] : Vy[1];
                    float px2 = diagAD ? Vx[3] : Vx[2], py2 = diagAD ? Vy[3] : Vy[2];
                    float dx = px2 - px1, dy = py2 - py1;
                    float L = sqrtf(dx*dx + dy*dy);
                    if (L > 2.0f) {
                        float H = L / sqrtf(ideal_aspect_ratio*ideal_aspect_ratio + 1.0f);
                        float W = H * ideal_aspect_ratio;
                        float phi_m = atan2f(dy, dx);
                        float phi_i = diagAD ? atan2f(H, W) : atan2f(H, -W);
                        float theta = phi_m - phi_i;
                        float c=cosf(theta), s=sinf(theta);
                        float wx=c*W, wy=s*W;
                        float hx=-s*H, hy=c*H;
                        if (diagAD) { Vx[1]=Vx[0]+wx; Vy[1]=Vy[0]+wy; Vx[2]=Vx[0]+hx; Vy[2]=Vy[0]+hy; }
                        else { Vx[0]=Vx[1]-wx; Vy[0]=Vy[1]-wy; Vx[3]=Vx[1]+hx; Vy[3]=Vy[1]+hy; }
                    }
                }
            }


            // --- FASE 3: ASSEGNAZIONE FINALE ---
            for (int i = 0; i < 4; i++) {
                FinalX[i] = roundf(Vx[i]);
                FinalY[i] = roundf(Vy[i]);
                positionXX[i] = FinalX[i];
                positionYY[i] = FinalY[i];
            }
            current_point_seen_mask = (1 << (3 - best_idx1)) | (1 << (3 - best_idx2));
        }
        #endif //COMMENTO

        #ifndef COMMENTO
        if (num_points_seen == 2)
        {
            // --- FASE 1: IDENTIFICAZIONE ---

            const int x1 = positionXX[0], y1 = positionYY[0];
            const int x2 = positionXX[1], y2 = positionYY[1];

            uint8_t best_idx1 = 0, best_idx2 = 1;
            int32_t min_total_cost = INT32_MAX;

            // --- Calcolo di TUTTI i possibili componenti di costo ---

            // 1. Costo di Posizione
            //                            (   BASI   )  ( ALTEZZE  )  (DIAGONALI )
            //                             A,B    C,D    A,C    B,D    A,D    B,C
            const uint8_t pairs[6][2] = { {0,1}, {2,3}, {0,2}, {1,3}, {0,3}, {1,2} };         
            int32_t pos_costs[12];
            // pos_costs
            //  0 -> AB (base) TOP
            //  1 -> BA (base) TOP
            //  2 -> CD (base) BOTTOM
            //  3 -> DC (base) BOTTOM
            //  4 -> AC (altezza) LEFT
            //  5 -> CA (altezza) LEFT
            //  6 -> BD (altezza) RIGHT
            //  7 -> DB (altezza) RIGHT
            //  8 -> AD (diagonale) left -> right
            //  9 -> DA (diagonale) left -> right
            // 10 -> BC (diagonale) right -> left
            // 11 -> CB (diagonale) right -> left

            for (uint8_t i = 0; i < 6; i++) {
                const uint8_t v1 = pairs[i][0], v2 = pairs[i][1];
                int32_t dx1a=x1-FinalX[v1], dy1a=y1-FinalY[v1];
                int32_t dx2a=x2-FinalX[v2], dy2a=y2-FinalY[v2];
                pos_costs[i*2] = dx1a*dx1a + dy1a*dy1a + dx2a*dx2a + dy2a*dy2a;
                if (pos_costs[i*2] < min_total_cost) {
                    min_total_cost = pos_costs[i*2];
                    best_idx1 = pairs[i][0]; best_idx2 = pairs[i][1];
                }
                int32_t dx1b=x1-FinalX[v2], dy1b=y1-FinalY[v2];
                int32_t dx2b=x2-FinalX[v1], dy2b=y2-FinalY[v1];
                pos_costs[i*2 + 1] = dx1b*dx1b + dy1b*dy1b + dx2b*dx2b + dy2b*dy2b;
                if (pos_costs[i*2 + 1] < min_total_cost) {
                    min_total_cost = pos_costs[i*2 + 1];
                    best_idx1 = pairs[i][1]; best_idx2 = pairs[i][0];
                }
            }

            // --- FASE 2: RICOSTRUZIONE GEOMETRICA ---

            int Vx[4], Vy[4];
            for (uint8_t i = 0; i < 4; i++) { Vx[i] = FinalX[i]; Vy[i] = FinalY[i]; }
            Vx[best_idx1] = x1; Vy[best_idx1] = y1;
            Vx[best_idx2] = x2; Vy[best_idx2] = y2;
    
            const uint8_t idx1 = best_idx1, idx2 = best_idx2;
            bool top    = (idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0);
            bool bottom = (idx1 == 2 && idx2 == 3) || (idx1 == 3 && idx2 == 2);
            bool left   = (idx1 == 0 && idx2 == 2) || (idx1 == 2 && idx2 == 0);
            bool right  = (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 1);

            // --- INTEGRAZIONE FRAME-2 PER VERIFICA DIREZIONE ---
            int dx_pred = medianX - prev2_medianX; // movimento previsto X
            int dy_pred = medianY - prev2_medianY; // movimento previsto Y
            // Soglia minima per considerare il movimento come significativo
            const int soglia_mov = 100; // puoi tararlo in base alla sensibilità
            // Correzione in base alla direzione prevista
            if (top && dy_pred > soglia_mov) { top = false; bottom = true; }
            else if (bottom && dy_pred < -soglia_mov) { bottom = false; top = true; }
            else if (left && dx_pred > soglia_mov) { left = false; right = true; }
            else if (right && dx_pred < -soglia_mov) { right = false; left = true; }
            
            if (top || bottom) {
                //                 A   C             B   D
                uint8_t p1 = top ? 0 : 2, p2 = top ? 1 : 3;
                //                 C   A             D   B
                uint8_t p3 = top ? 2 : 0, p4 = top ? 3 : 1;
                
                if (Vx[p1] > Vx[p2]) { uint8_t t=p1; p1=p2; p2=t; }
                float dx = Vx[p2]-Vx[p1], dy = Vy[p2]-Vy[p1];
                float base = sqrtf(dx*dx + dy*dy);
                if (base > 2.0f) {
                    float cos_angle = fabsf(dx) / base;
                    base = base / fmaxf(cos_angle, 0.3f);
                    float newH = fmaxf(base / ideal_aspect_ratio, height * 0.5f);
                    float inv_base = 1.0f / base;
                    float ndx = -dy*inv_base, ndy = dx*inv_base;
                    float sx = ndx*newH, sy = ndy*newH;
                    if (top) {
                        Vx[p3]=Vx[p1]+sx; Vy[p3]=Vy[p1]+sy; Vx[p4]=Vx[p2]+sx; Vy[p4]=Vy[p2]+sy;
                    } else {
                        Vx[p3]=Vx[p1]-sx; Vy[p3]=Vy[p1]-sy; Vx[p4]=Vx[p2]-sx; Vy[p4]=Vy[p2]-sy;
                    }
                }
            } else if (left || right) {
                //                  A   B              C   D
                uint8_t p1 = left ? 0 : 1, p2 = left ? 2 : 3;
                //                  B   A              D   C
                uint8_t p3 = left ? 1 : 0, p4 = left ? 3 : 2;
                
                if (Vy[p1] > Vy[p2]) { uint8_t t=p1; p1=p2; p2=t; }
                float dx = Vx[p2]-Vx[p1], dy = Vy[p2]-Vy[p1];
                float hgt = sqrtf(dx*dx+dy*dy);
                if (hgt > 2.0f) {
                    float newW = hgt * ideal_aspect_ratio;
                    float inv_hgt = 1.0f/hgt;
                    float ndx = -dy*inv_hgt, ndy = dx*inv_hgt;
                    float sx = ndx*newW, sy = ndy*newW;
                    if (left) {
                        Vx[p3]=Vx[p1]-sx; Vy[p3]=Vy[p1]-sy; Vx[p4]=Vx[p2]-sx; Vy[p4]=Vy[p2]-sy;
                    } else {
                        Vx[p3]=Vx[p1]+sx; Vy[p3]=Vy[p1]+sy; Vx[p4]=Vx[p2]+sx; Vy[p4]=Vy[p2]+sy;
                    }
                }
            } else { // Diagonale
                bool diagAD = (idx1 == 0 && idx2 == 3) || (idx1 == 3 && idx2 == 0);
                float px1 = diagAD ? Vx[0] : Vx[1], py1 = diagAD ? Vy[0] : Vy[1];
                float px2 = diagAD ? Vx[3] : Vx[2], py2 = diagAD ? Vy[3] : Vy[2];
                float dx = px2 - px1, dy = py2 - py1;
                float L = sqrtf(dx*dx + dy*dy);
                if (L > 2.0f) {
                    float H = L / sqrtf(ideal_aspect_ratio*ideal_aspect_ratio + 1.0f);
                    float W = H * ideal_aspect_ratio;
                    float phi_m = atan2f(dy, dx);
                    float phi_i = diagAD ? atan2f(H, W) : atan2f(H, -W);
                    float theta = phi_m - phi_i;
                    float c=cosf(theta), s=sinf(theta);
                    float wx=c*W, wy=s*W;
                    float hx=-s*H, hy=c*H;
                    if (diagAD) {
                        Vx[1]=Vx[0]+wx; Vy[1]=Vy[0]+wy; Vx[2]=Vx[0]+hx; Vy[2]=Vy[0]+hy;
                    } else {
                        Vx[0]=Vx[1]-wx; Vy[0]=Vy[1]-wy; Vx[3]=Vx[1]+hx; Vy[3]=Vy[1]+hy;
                    }
                }
            }

            // --- FASE 3: ASSEGNAZIONE FINALE ---

            for (uint8_t i = 0; i < 4; i++) {
                FinalX[i] = roundf(Vx[i]);
                FinalY[i] = roundf(Vy[i]);
                positionXX[i] = FinalX[i];
                positionYY[i] = FinalY[i];
            }
            current_point_seen_mask = (1 << (3 - best_idx1)) | (1 << (3 - best_idx2));
            
            
        }
        #endif //COMMENTO
        else
        if (num_points_seen == 3)
        {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////
            // ==================== logica per 3 punti disordinati ============================================== //
            ////////////////////////////////////////////////////////////////////////////////////////////////////////            
            
            // Se vediamo 3 punti, eseguiamo la stima semplice e veloce del quarto punto.
            int32_t d01_sq, d12_sq, d02_sq;
            
            dx = positionXX[0] - positionXX[1]; dy = positionYY[0] - positionYY[1];
            d01_sq = dx * dx + dy * dy;

            dx = positionXX[1] - positionXX[2]; dy = positionYY[1] - positionYY[2];
            d12_sq = dx * dx + dy * dy;

            dx = positionXX[0] - positionXX[2]; dy = positionYY[0] - positionYY[2];
            d02_sq = dx * dx + dy * dy;
            
            uint8_t a_idx;  // uno dei due punti del lato più lungo del triangolo che corrisponde ad una diagonale del rettangolo
            uint8_t b_idx;  // vertice tra lato lungo e lato corto del triangolo 
            uint8_t c_idx;  // uno dei due punti del lato più lungo del triangolo che corrisponde ad una diagonale del rettangolo

            if (d01_sq >= d12_sq && d01_sq >= d02_sq) {
                // d01 è la diagonale del rettangolo
                // b, il punto tra lato lungo e lato corto è per forza il punto 2 
                a_idx = 0; c_idx = 1; b_idx = 2;
                //current_point_seen_mask = 0b00001111;
            }
            else if (d12_sq >= d02_sq) {
                //d12 è la diagonale  del rettangolo
                // b, il punto tra lato lungo e lato corto è per forza il punto 0
                a_idx = 1; c_idx = 2; b_idx = 0;
                //current_point_seen_mask = 0b00001111;
            }
            else {
                //d02 è la diagonale  del rettangolo
                // b, il punto tra lato lungo e lato corto è per forza il punto 2 
                a_idx = 0; c_idx = 2; b_idx = 1;
                //current_point_seen_mask = 0b00001111;
            }
            
            // utilizzando le proprietà del Parallelogramma calcolo il 4 vertice mancante, 
            // senza comunque conoscerne l'esatta posizione logica, che verrà individuata successiva dalla gestione dei 4 punti
            positionXX[3] = positionXX[a_idx] + positionXX[c_idx] - positionXX[b_idx];
            positionYY[3] = positionYY[a_idx] + positionYY[c_idx] - positionYY[b_idx];

            point4_X = positionXX[3];
            point4_Y = positionYY[3];

        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // ==================== logica per 4 punti disordinati ============================================== //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // --- FASE 3: Ordinamento Finale dei 4 Punti ---
        // Questa sezione viene eseguita sempre, per garantire che FinalX/Y
        // abbiano sempre un ordine coerente A,B,C,D.
        
        uint8_t orderX[4] = {0, 1, 2, 3};
        uint8_t orderY[4] = {0, 1, 2, 3};

        // ordinanmento più verboso, ma più performante
        // Ordinamento degli indici in base alla coordinata X (sorting network ottimizzato)
        {
            uint8_t tmp;
            if (positionXX[orderX[0]] > positionXX[orderX[1]]) { tmp = orderX[0]; orderX[0] = orderX[1]; orderX[1] = tmp; }
            if (positionXX[orderX[2]] > positionXX[orderX[3]]) { tmp = orderX[2]; orderX[2] = orderX[3]; orderX[3] = tmp; }
            if (positionXX[orderX[0]] > positionXX[orderX[2]]) { tmp = orderX[0]; orderX[0] = orderX[2]; orderX[2] = tmp; }
            if (positionXX[orderX[1]] > positionXX[orderX[3]]) { tmp = orderX[1]; orderX[1] = orderX[3]; orderX[3] = tmp; }
            if (positionXX[orderX[1]] > positionXX[orderX[2]]) { tmp = orderX[1]; orderX[1] = orderX[2]; orderX[2] = tmp; }
        }

        // ordinanmento più verboso, ma più performante
        // Ordinamento degli indici in base alla coordinata Y (sorting network ottimizzato)
        {
            uint8_t tmp;
            if (positionYY[orderY[0]] > positionYY[orderY[1]]) { tmp = orderY[0]; orderY[0] = orderY[1]; orderY[1] = tmp; }
            if (positionYY[orderY[2]] > positionYY[orderY[3]]) { tmp = orderY[2]; orderY[2] = orderY[3]; orderY[3] = tmp; }
            if (positionYY[orderY[0]] > positionYY[orderY[2]]) { tmp = orderY[0]; orderY[0] = orderY[2]; orderY[2] = tmp; }
            if (positionYY[orderY[1]] > positionYY[orderY[3]]) { tmp = orderY[1]; orderY[1] = orderY[3]; orderY[3] = tmp; }
            if (positionYY[orderY[1]] > positionYY[orderY[2]]) { tmp = orderY[1]; orderY[1] = orderY[2]; orderY[2] = tmp; }
        }

        // Assegnazione dei vertici A,B,C,D tramite euristiche
        int32_t dist_sq1, dist_sq2;

        dx = positionXX[orderY[0]] - positionXX[orderX[0]];
        dy = positionYY[orderY[0]] - positionYY[orderX[0]];
        dist_sq1 = (dx * dx) + (dy * dy);

        dx = positionXX[orderY[3]] - positionXX[orderX[0]];
        dy = positionYY[orderY[3]] - positionYY[orderX[0]];
        dist_sq2 = (dx * dx) + (dy * dy);
        
        const int CRITICAL_ZONE = (30 * CamToMouseMult); // 30 è un valore provato in tante situazioni e pare andare bene

        if ((positionYY[orderY[1]] - positionYY[orderY[0]]) > CRITICAL_ZONE)
        {
            // Caso Normale
            if (dist_sq1 < dist_sq2) {
                a = orderX[0]; d = orderX[3];
                if (orderX[1] == orderY[3]) { c = orderX[1]; b = orderX[2]; }
                else { b = orderX[1]; c = orderX[2]; }
            } else {
                c = orderX[0]; b = orderX[3];
                if (orderX[1] == orderY[3]) { d = orderX[1]; a = orderX[2]; }
                else { a = orderX[1]; d = orderX[2]; }
            }
        }
        else
        {
            // Caso Zona Critica (rettangolo quasi verticale)
            a = orderY[0]; b = orderY[1];
            c = orderY[2]; d = orderY[3];
        }

        // sarebbe sufficiente per la zona critica questo ulteriore controllo, ma poichè ha
        // un piccolo peso computazionale, lo facciamo sempre per casi limiti e come maggiore robustezza 
        // Correzione finale per garantire la convenzione (A=TL, B=TR, C=BL, D=BR)
        {
            uint8_t aux_swap;
            if (positionYY[a] > positionYY[c]) { aux_swap = a; a = c; c = aux_swap; }
            if (positionYY[b] > positionYY[d]) { aux_swap = b; b = d; d = aux_swap; }
            if (positionXX[a] > positionXX[b]) { aux_swap = a; a = b; b = aux_swap; }
            if (positionXX[c] > positionXX[d]) { aux_swap = c; c = d; d = aux_swap; }
        }
              
        // =========================================================================//
        // ==========   GESTIONE 3 PUNTI VISTI ED UNO DA STIMARE   =================//
        // ========== OTTIMIZZAZIONE COORDINATE DEL PUNTO MANCANTE =================//
        // =========================================================================//
        
        if (num_points_seen == 3) {
            uint8_t idx1, idx2, idx3;
            uint8_t point_idx_aux;
            if ((positionXX[a] == point4_X) && (positionYY[a] == point4_Y)) {
                point4_idx = A; point_idx_aux = a;
                idx1 = B; idx2 = C; idx3 = D;
            } else if ((positionXX[b] == point4_X) && (positionYY[b] == point4_Y)) {
                point4_idx = B; point_idx_aux = b;
                idx1 = A; idx2 = D; idx3 = C;
            } else if ((positionXX[c] == point4_X) && (positionYY[c] == point4_Y)) {
                point4_idx = C; point_idx_aux = c;
                idx1 = A; idx2 = D; idx3 = B;
            } else {
                point4_idx = D; point_idx_aux = d;
                idx1 = B; idx2 = C; idx3 = A;
            }  

            // Calcolo della Base nel frame precedente
            const int prev_base_x = FinalX[idx1] + FinalX[idx2] - FinalX[idx3];
            const int prev_base_y = FinalY[idx1] + FinalY[idx2] - FinalY[idx3];

            // Calcolo del vettore di Errore/Correzione
            const int error_x = FinalX[point4_idx] - prev_base_x;
            const int error_y = FinalY[point4_idx] - prev_base_y;

            // risultato finale 
            positionXX[point_idx_aux] += error_x;
            positionYY[point_idx_aux] += error_y;  
        }

        // =========================================================================//
        // ======== FINE OTTIMIZZAZIONE COORDINATE DEL PUNTO MANCANTE ==============//
        // =========================================================================//
 

        // --- FASE 4: Assegnazione Finale e Calcoli Derivati ---
        FinalX[A] = positionXX[a]; FinalY[A] = positionYY[a];
        FinalX[B] = positionXX[b]; FinalY[B] = positionYY[b];
        FinalX[C] = positionXX[c]; FinalY[C] = positionYY[c];
        FinalX[D] = positionXX[d]; FinalY[D] = positionYY[d];


        //////////////////////////////////////////////////////////////////////////////////////////
        // ==== calcolo medianX, medianY, height, width, angle ... per compatibilità ========== //
        //////////////////////////////////////////////////////////////////////////////////////////
        
        prev2_medianX = medianX;
        prev2_medianY = medianY;
        //prev2_height = height;
        //prev2_width = width;
        //prev2_angle = angle;

        // Calcolo del centroide con arrotondamento (il +2 è un trucco per l'arrotondamento all'intero più vicino)
        medianX = (FinalX[A] + FinalX[B] + FinalX[C] + FinalX[D] + 2) / 4;
        medianY = (FinalY[A] + FinalY[B] + FinalY[C] + FinalY[D] + 2) / 4;

        // Calcoli finali mantenuti per compatibilità
        height_left = hypotf((float)FinalY[A] - FinalY[C], (float)FinalX[A] - FinalX[C]); // lunghezza lato AC
        height_right = hypotf((float)FinalY[B] - FinalY[D], (float)FinalX[B] - FinalX[D]); // lunghezza lato BD
        width_top = hypotf((float)FinalY[A] - FinalY[B], (float)FinalX[A] - FinalX[B]); // lunghezza lato AB
        width_bottom = hypotf((float)FinalY[C] - FinalY[D], (float)FinalX[C] - FinalX[D]); // lunghezza lato CD
        height = (height_left + height_right) / 2.0f;
        width = (width_top + width_bottom) / 2.0f;
        
        if (num_points_seen == 4) {
            // NUOVO: Aggiorna la nostra "verità assoluta" sulla forma
            if (height > 1e-6f) { // Evita divisione per zero
                ideal_aspect_ratio = width / height;
            }

            // Vettore A->B
            shape_vec_AB_x = FinalX[1] - FinalX[0];
            shape_vec_AB_y = FinalY[1] - FinalY[0];

            // Vettore A->D
            shape_vec_AD_x = FinalX[3] - FinalX[0];
            shape_vec_AD_y = FinalY[3] - FinalY[0];

            // DA INSERIRE ALLA FINE DEL BLOCCO if (num_points_seen == 4)

            // Calcola dove sarebbe il punto C se fosse un parallelogramma perfetto
            //int parallelogram_C_x = FinalX[0] + FinalX[2] - FinalX[1]; // In realtà è A+C-B per trovare D, quindi C=B+D-A
            int parallelogram_C_x = FinalX[1] + FinalX[3] - FinalX[0]; // C = B + D - A
            int parallelogram_C_y = FinalY[1] + FinalY[3] - FinalY[0];
        }

        angle = (atan2f((float)FinalY[A] - FinalY[B], (float)FinalX[B] - FinalX[A]) + atan2f((float)FinalY[C] - FinalY[D], (float)FinalX[D] - FinalX[C])) / 2.0f;      

        is_tracking_stable = true;
        prev_point_seen_mask = current_point_seen_mask; 
        if (num_points_seen >= 3) current_point_seen_mask = 0b00001111;
    }
    else {
    is_tracking_stable = false;
    prev_point_seen_mask = 0;
}
    prev_num_points_seen = num_points_seen;
}

#endif //USE_SQUARE_ADVANCED