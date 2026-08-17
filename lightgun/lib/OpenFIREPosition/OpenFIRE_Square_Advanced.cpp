#ifdef USE_SQUARE_ADVANCED
/*!
 * @file OpenFIRE_Square_Advanced.cpp
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

#include <Arduino.h>
#include "OpenFIRE_Square_Advanced.h"
#include <cfloat>

/**
 * @brief Calcola la radice quadrata inversa (1/sqrt(x)) in modo approssimato ma estremamente veloce.
 * Nota Architetturale: Attualmente disattivata per non corrompere la precisione sub-pixel (errore 0.17%). 
 * Mantenuta a scopo di documentazione per eventuali downgrade prestazionali critici.
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
 * Sostituisce la funzione standard round() per risparmiare overhead chiamate e sfruttare
 * la cast-truncation nativa dell'ALU dell'ESP32. L'addizione/sottrazione di 0.5f sposta
 * la soglia geometrica senza richiedere rami condizionali pesanti.
 * @param val Il valore float da arrotondare.
 * @return L'intero arrotondato.
 */
inline int fast_roundf(float val) {
    return (int)(val + (val > 0.0f ? 0.5f : -0.5f));
}

// Definizione degli indici per chiarezza
#define A 0
#define B 1
#define C 2
#define D 3


void OpenFIRE_Square::begin(const int* px, const int* py, unsigned int seen) {

    // =========================================================================//
    // ================= FASE 0: INIT E ACQUISIZIONE DATI REALI ================//
    // =========================================================================//
    
    // --- FASE 0: Inizializzazione e Controllo di Visibilità ---
    
    seenFlags = seen;

    // Il sistema richiede un lock visivo completo su tutti i LED almeno una volta (0x0F = 1111).
    // Questo garantisce che la metrica geometrica iniziale (ideal_aspect_ratio) sia calcolata 
    // su dati reali e non su approssimazioni di partenza.

    // Il sistema richiede di vedere tutti e 4 i sensori almeno una volta per inizializzarsi.
    if (seenFlags == 0x0F) { // 0x0F in binario è 1111, tutti i sensori visti.
        start = 0xFF;
    } else if (!start) {
        return; // Esci se non siamo ancora stati inizializzati.
    }

    uint8_t num_points_seen = 0;
    int positionXX[4];
    int positionYY[4];
    
    // real_x/y conservano i dati "crudi" letti dalla telecamera prima di qualsiasi predizione.
    // Ci serviranno alla fine nel Setaccio Universale per capire quali punti della matrice 
    // finale sono ottici e quali sono stati inventati dalla matematica.
    int real_x[4]; 
    int real_y[4]; 

    // Fondiamo estrazione, shift e backup in un unico ciclo per saturare la cache L1 
    // e azzerare cicli e salti logici inutili della CPU.

    // Estrae i punti visibili, li mette negli array di lavoro e applica la trasformazione.
    for (uint8_t i = 0; i < 4; ++i) {
        if ((seenFlags >> i) & 0x01) {
            int calc_x = (CamMaxX - px[i]) << CamToMouseShift;
            int calc_y = py[i] << CamToMouseShift;
            
            positionXX[num_points_seen] = calc_x;
            positionYY[num_points_seen] = calc_y;
            real_x[num_points_seen] = calc_x;
            real_y[num_points_seen] = calc_y;
            
            num_points_seen++;
        }
    }
    
    
    if (num_points_seen >= 1) {

        // =========================================================================//
        // ======== CANCELLO DI SICUREZZA (COLD START DOPO IL BUIO) ================//
        // =========================================================================//
        // Se usciamo da un momento di cecità parziale o totale, rifiutiamo di indovinare 
        // la posizione con meno di 3 punti per evitare salti impazziti del mirino.
        if (!is_tracking_stable && num_points_seen < 3) {
            prev_num_points_seen = 0; 
            return; 
        }

        int32_t dx, dy;

        // =========================================================================//
        // ========================= GESTIONE 1 PUNTO ==============================//
        // =========================================================================//
        // Quando vediamo 1 solo sensore, abbandoniamo l'analisi angolare. Assumiamo per 
        // inerzia pura che l'arma abbia traslato senza ruotare sull'asse Z.
        if (num_points_seen == 1) {
            if (prev_num_points_seen >= 1) {
                const int x1 = positionXX[0]; 
                const int y1 = positionYY[0];

                uint8_t identified_idx = 0;
                int32_t min_dist_sq = INT32_MAX; 

                // Cerchiamo a quale dei 4 vecchi punti appartiene quello nuovo valutando
                // la distanza al quadrato (risparmia la FPU evitando sqrtf).
                for (uint8_t i = 0; i < 4; i++) {
                    const int dx_dist = x1 - FinalX[i];
                    const int dy_dist = y1 - FinalY[i];
                    const int32_t dist_sq = (int32_t)dx_dist * dx_dist + (int32_t)dy_dist * dy_dist;

                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        identified_idx = i;
                    }
                }

                // Meccanismo anti-teletrasporto: se la distanza del punto trovato dal frame 
                // precedente è fisicamente irrealistica, invalidiamo il tracking.
                const int max_jump_distance = ((int)width * 3) / 4; 
                const int32_t MAX_ALLOWED_DISTANCE_SQ = (int32_t)max_jump_distance * max_jump_distance;

                if (min_dist_sq > MAX_ALLOWED_DISTANCE_SQ) {
                    prev_num_points_seen = 0;
                    return; 
                }

                const int delta_x = x1 - FinalX[identified_idx];
                const int delta_y = y1 - FinalY[identified_idx];

                for (uint8_t i = 0; i < 4; i++) {
                    positionXX[i] = FinalX[i] + delta_x;
                    positionYY[i] = FinalY[i] + delta_y;
                }
            } else {
                prev_num_points_seen = 0;
                return;
            }
        }
        
        // =========================================================================//
        // ========================= GESTIONE 2 PUNTI ==============================//
        // =========================================================================//
        // Risolviamo l'ambiguità ottica (Quale lato del rettangolo sto guardando?).
        else if (num_points_seen == 2) {
            const int x1 = positionXX[0], y1 = positionYY[0];
            const int x2 = positionXX[1], y2 = positionYY[1];

            uint8_t best_idx1 = 0, best_idx2 = 1;
            int32_t min_total_cost = INT32_MAX;
            
            // Le 6 possibili combinazioni geometriche dei lati del rettangolo.
            //                            (   BASI   )  ( ALTEZZE  )  (DIAGONALI )
            //                             A,B    C,D    A,C    B,D    A,D    B,C            
            const uint8_t pairs[6][2] = { {0,1}, {2,3}, {0,2}, {1,3}, {0,3}, {1,2} };         

            int32_t read_dx = x2 - x1;
            int32_t read_dy = y2 - y1;
            int32_t read_dist_sq = read_dx * read_dx + read_dy * read_dy;

            // OTTIMIZZAZIONE 2: Pre-calcolo fuori dal ciclo per alleggerire la ALU.
            const int32_t w_sq = (int32_t)(width * width);
            const int32_t h_sq = (int32_t)(height * height);
            const int32_t expected_dist_arr[6] = {w_sq, w_sq, h_sq, h_sq, w_sq + h_sq, w_sq + h_sq};

            // Il calcolo dei costi è un sistema a punteggio che confronta l'ipotesi corrente
            // contro lo stato geometrico e storico noto, premiando le deduzioni conservative.
            for (uint8_t i = 0; i < 6; i++) {
                const uint8_t v1 = pairs[i][0], v2 = pairs[i][1];
                int32_t penalty = 0;
                
                if (prev_num_points_seen >= 1 && prev_point_seen_mask != 0) {
                    if (!(prev_point_seen_mask & (1 << (3 - v1)))) penalty += PENALITA_STORICA;
                    if (!(prev_point_seen_mask & (1 << (3 - v2)))) penalty += PENALITA_STORICA;
                }

                penalty += abs(read_dist_sq - expected_dist_arr[i]);

                int32_t dx1a=x1-FinalX[v1], dy1a=y1-FinalY[v1];
                int32_t dx2a=x2-FinalX[v2], dy2a=y2-FinalY[v2];
                int32_t costA = dx1a*dx1a + dy1a*dy1a + dx2a*dx2a + dy2a*dy2a + penalty;
                if (costA < min_total_cost) { min_total_cost = costA; best_idx1 = v1; best_idx2 = v2; }

                int32_t dx1b=x1-FinalX[v2], dy1b=y1-FinalY[v2];
                int32_t dx2b=x2-FinalX[v1], dy2b=y2-FinalY[v1];
                int32_t costB = dx1b*dx1b + dy1b*dy1b + dx2b*dx2b + dy2b*dy2b + penalty;
                if (costB < min_total_cost) { min_total_cost = costB; best_idx1 = v2; best_idx2 = v1; }
            }

            float Vx[4], Vy[4];
            for (uint8_t i = 0; i < 4; i++) { Vx[i] = (float)FinalX[i]; Vy[i] = (float)FinalY[i]; }
            Vx[best_idx1] = (float)x1; Vy[best_idx1] = (float)y1;
            Vx[best_idx2] = (float)x2; Vy[best_idx2] = (float)y2;
    
            const uint8_t idx1 = best_idx1, idx2 = best_idx2;
            bool top    = (idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0);
            bool bottom = (idx1 == 2 && idx2 == 3) || (idx1 == 3 && idx2 == 2);
            bool left   = (idx1 == 0 && idx2 == 2) || (idx1 == 2 && idx2 == 0);
            bool right  = (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 1);

            // Per prevenire il "Bug Clessidra" (incrocio errato dei vettori derivati),
            // viene riallineato dinamicamente l'indice dei punti e scalata la proiezione
            // imponendo con la forza matematica il rispetto dell'aspect_ratio originale.
            if (top || bottom) {
                //                 A   C             B   D
                uint8_t p1 = top ? 0 : 2, p2 = top ? 1 : 3;
                //                 C   A             D   B
                uint8_t p3 = top ? 2 : 0, p4 = top ? 3 : 1;
                if (Vx[p1] > Vx[p2]) { 
                    uint8_t t=p1; p1=p2; p2=t; t=p3; p3=p4; p4=t; 
                }
                float dx_v = Vx[p2]-Vx[p1], dy_v = Vy[p2]-Vy[p1];
                float base = sqrtf(dx_v*dx_v + dy_v*dy_v);
                
                if (base > 2.0f) {
                    float newH = fmaxf(base / ideal_aspect_ratio, height * 0.5f);
                    float inv_base = 1.0f / base;
                    float sx = -dy_v*inv_base*newH, sy = dx_v*inv_base*newH;
                    if (top) { Vx[p3]=Vx[p1]+sx; Vy[p3]=Vy[p1]+sy; Vx[p4]=Vx[p2]+sx; Vy[p4]=Vy[p2]+sy; } 
                    else { Vx[p3]=Vx[p1]-sx; Vy[p3]=Vy[p1]-sy; Vx[p4]=Vx[p2]-sx; Vy[p4]=Vy[p2]-sy; }
                }
            } else if (left || right) {
                //                  A   B              C   D
                uint8_t p1 = left ? 0 : 1, p2 = left ? 2 : 3;
                //                  B   A              D   C
                uint8_t p3 = left ? 1 : 0, p4 = left ? 3 : 2;
                if (Vy[p1] > Vy[p2]) { 
                    uint8_t t=p1; p1=p2; p2=t; t=p3; p3=p4; p4=t; 
                }
                float dx_v = Vx[p2]-Vx[p1], dy_v = Vy[p2]-Vy[p1];
                float hgt = sqrtf(dx_v*dx_v+dy_v*dy_v);
                if (hgt > 2.0f) {
                    float newW = hgt * ideal_aspect_ratio;
                    float inv_hgt = 1.0f/hgt;
                    float sx = -dy_v*inv_hgt*newW, sy = dx_v*inv_hgt*newW;
                    if (left) { Vx[p3]=Vx[p1]-sx; Vy[p3]=Vy[p1]-sy; Vx[p4]=Vx[p2]-sx; Vy[p4]=Vy[p2]-sy; } 
                    else { Vx[p3]=Vx[p1]+sx; Vy[p3]=Vy[p1]+sy; Vx[p4]=Vx[p2]+sx; Vy[p4]=Vy[p2]+sy; }
                }
            } else { // Gestione Diagonale
                bool diagAD = (idx1 == 0 && idx2 == 3) || (idx1 == 3 && idx2 == 0);
                float px1 = diagAD ? Vx[0] : Vx[1], py1 = diagAD ? Vy[0] : Vy[1];
                float px2 = diagAD ? Vx[3] : Vx[2], py2 = diagAD ? Vy[3] : Vy[2];
                float dx_v = px2 - px1, dy_v = py2 - py1;
                float L = sqrtf(dx_v*dx_v + dy_v*dy_v);
                if (L > 2.0f) {
                    float H = L / sqrtf(ideal_aspect_ratio*ideal_aspect_ratio + 1.0f);
                    float W = H * ideal_aspect_ratio;
                    float phi_m = atan2f(dy_v, dx_v);
                    float phi_i = diagAD ? atan2f(H, W) : atan2f(H, -W);
                    float theta = phi_m - phi_i;
                    float c_a=cosf(theta), s_a=sinf(theta);
                    float wx=c_a*W, wy=s_a*W;
                    float hx=-s_a*H, hy=c_a*H;
                    if (diagAD) { Vx[1]=Vx[0]+wx; Vy[1]=Vy[0]+wy; Vx[2]=Vx[0]+hx; Vy[2]=Vy[0]+hy; } 
                    else { Vx[0]=Vx[1]-wx; Vy[0]=Vy[1]-wy; Vx[3]=Vx[1]+hx; Vy[3]=Vy[1]+hy; }
                }
            }

            for (uint8_t i = 0; i < 4; i++) { positionXX[i] = fast_roundf(Vx[i]); positionYY[i] = fast_roundf(Vy[i]); }
        }
        
        // =========================================================================//
        // ========================= GESTIONE 3 PUNTI ==============================//
        // =========================================================================//
        // Sfrutta il teorema del parallelogramma. Individuato il vertice opposto (che
        // sottende la diagonale più lunga), genera il quarto vertice mancante mediante
        // pura addizione vettoriale. Molto economico per la CPU.
        else if (num_points_seen == 3) {
            // Se vediamo 3 punti, eseguiamo la stima semplice e veloce del quarto punto.
            
            dx = positionXX[0] - positionXX[1]; dy = positionYY[0] - positionYY[1];
            int32_t d01_sq = dx * dx + dy * dy;

            dx = positionXX[1] - positionXX[2]; dy = positionYY[1] - positionYY[2];
            int32_t d12_sq = dx * dx + dy * dy;

            dx = positionXX[0] - positionXX[2]; dy = positionYY[0] - positionYY[2];
            int32_t d02_sq = dx * dx + dy * dy;
            
            uint8_t a_idx;  // uno dei due punti del lato più lungo del triangolo che corrisponde ad una diagonale del rettangolo
            uint8_t b_idx;  // vertice tra lato lungo e lato corto del triangolo 
            uint8_t c_idx;  // uno dei due punti del lato più lungo del triangolo che corrisponde ad una diagonale del rettangolo

            if (d01_sq >= d12_sq && d01_sq >= d02_sq) { 
                // d01 è la diagonale del rettangolo
                // b, il punto tra lato lungo e lato corto è per forza il punto 2 
                a_idx = 0; c_idx = 1; b_idx = 2; 
            } 
            else if (d12_sq >= d02_sq) { 
                //d12 è la diagonale  del rettangolo
                // b, il punto tra lato lungo e lato corto è per forza il punto 0
                a_idx = 1; c_idx = 2; b_idx = 0; 
            } 
            else { 
                //d02 è la diagonale  del rettangolo
                // b, il punto tra lato lungo e lato corto è per forza il punto 2 
                a_idx = 0; c_idx = 2; b_idx = 1; 
            }
                      
            // utilizzando le proprietà del Parallelogramma calcolo il 4 vertice mancante, 
            // senza comunque conoscerne l'esatta posizione logica, che verrà individuata successiva dalla gestione dei 4 punti            
            positionXX[3] = positionXX[a_idx] + positionXX[c_idx] - positionXX[b_idx];
            positionYY[3] = positionYY[a_idx] + positionYY[c_idx] - positionYY[b_idx];
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // =============== logica per l'ordinamento dei 4 punti disordinati ================================= //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // =========================================================================//
        // ======== ORDINAMENTO PER TUTTI I PUNTI E IDENTIFICAZIONE A,B,C,D ========//
        // =========================================================================//
        
        // Questa sezione viene eseguita sempre, per garantire che FinalX/Y
        // abbiano sempre un ordine coerente A,B,C,D.

        uint8_t orderX[4] = {0, 1, 2, 3}, orderY[4] = {0, 1, 2, 3}, a, b, c, d;

        // Implementazione di un "Sorting Network" statico a 5 confronti. 
        // È matematicamente provato essere la strada più rapida per ordinare 4 elementi 
        // senza incorrere nei pesanti overhead logici tipici dei loop di ordinamento classici (es. Bubble Sort).
        
        // ordinanmento più verboso, ma più performante
        // Ordinamento degli indici in base alla coordinata X (sorting network ottimizzato)
        { uint8_t tmp;
            if (positionXX[orderX[0]] > positionXX[orderX[1]]) { tmp = orderX[0]; orderX[0] = orderX[1]; orderX[1] = tmp; }
            if (positionXX[orderX[2]] > positionXX[orderX[3]]) { tmp = orderX[2]; orderX[2] = orderX[3]; orderX[3] = tmp; }
            if (positionXX[orderX[0]] > positionXX[orderX[2]]) { tmp = orderX[0]; orderX[0] = orderX[2]; orderX[2] = tmp; }
            if (positionXX[orderX[1]] > positionXX[orderX[3]]) { tmp = orderX[1]; orderX[1] = orderX[3]; orderX[3] = tmp; }
            if (positionXX[orderX[1]] > positionXX[orderX[2]]) { tmp = orderX[1]; orderX[1] = orderX[2]; orderX[2] = tmp; }
        }

        // ordinanmento più verboso, ma più performante
        // Ordinamento degli indici in base alla coordinata Y (sorting network ottimizzato)        
        { uint8_t tmp;
            if (positionYY[orderY[0]] > positionYY[orderY[1]]) { tmp = orderY[0]; orderY[0] = orderY[1]; orderY[1] = tmp; }
            if (positionYY[orderY[2]] > positionYY[orderY[3]]) { tmp = orderY[2]; orderY[2] = orderY[3]; orderY[3] = tmp; }
            if (positionYY[orderY[0]] > positionYY[orderY[2]]) { tmp = orderY[0]; orderY[0] = orderY[2]; orderY[2] = tmp; }
            if (positionYY[orderY[1]] > positionYY[orderY[3]]) { tmp = orderY[1]; orderY[1] = orderY[3]; orderY[3] = tmp; }
            if (positionYY[orderY[1]] > positionYY[orderY[2]]) { tmp = orderY[1]; orderY[1] = orderY[2]; orderY[2] = tmp; }
        }

        // Assegnazione dei vertici A,B,C,D tramite euristiche
        dx = positionXX[orderY[0]] - positionXX[orderX[0]]; dy = positionYY[orderY[0]] - positionYY[orderX[0]];
        int32_t dist_sq1 = (dx * dx) + (dy * dy);

        dx = positionXX[orderY[3]] - positionXX[orderX[0]]; dy = positionYY[orderY[3]] - positionYY[orderX[0]];
        int32_t dist_sq2 = (dx * dx) + (dy * dy);
        
        /*
        const int CRITICAL_ZONE = (30 * CamToMouseMult); // 30 è un valore provato in tante situazioni e pare andare bene
        */

        // Percentuale storica testata sulla DFRobot (30 pixel su una altezza di 768)
        constexpr float CRITICAL_ZONE_PERCENT = 30.0f / (float)768;
        // Calcoliamo il valore dinamico sulla risoluzione della CAM e lo portiamo nello spazio unificato
        const int CRITICAL_ZONE = (int)((float)CamResY * CRITICAL_ZONE_PERCENT * (float)CamToMouseMult);

        // Se l'arma è inclinata severamente (oltre ~45 gradi, superando la CRITICAL_ZONE), 
        // l'ordinamento ingenuo Y non basta. Valutiamo le distanze relative per mantenere coerenti
        // le etichette A, B, C, D e impedire l'incrocio letale degli assi.
        if ((positionYY[orderY[1]] - positionYY[orderY[0]]) > CRITICAL_ZONE) {
            // Caso Normale
            if (dist_sq1 < dist_sq2) {
                a = orderX[0]; d = orderX[3];
                if (orderX[1] == orderY[3]) { c = orderX[1]; b = orderX[2]; } else { b = orderX[1]; c = orderX[2]; }
            } else {
                c = orderX[0]; b = orderX[3];
                if (orderX[1] == orderY[3]) { d = orderX[1]; a = orderX[2]; } else { a = orderX[1]; d = orderX[2]; }
            }
        } else { 
            // Caso Zona Critica (rettangolo quasi verticale)
            a = orderY[0]; b = orderY[1]; c = orderY[2]; d = orderY[3]; 
        }

        // sarebbe sufficiente per la zona critica questo ulteriore controllo, ma poichè ha
        // un piccolo peso computazionale, lo facciamo sempre per casi limiti e come maggiore robustezza 
        // Correzione finale per garantire la convenzione (A=TL, B=TR, C=BL, D=BR)        
        { uint8_t aux_swap;
            if (positionYY[a] > positionYY[c]) { aux_swap = a; a = c; c = aux_swap; }
            if (positionYY[b] > positionYY[d]) { aux_swap = b; b = d; d = aux_swap; }
            if (positionXX[a] > positionXX[b]) { aux_swap = a; a = b; b = aux_swap; }
            if (positionXX[c] > positionXX[d]) { aux_swap = c; c = d; d = aux_swap; }
        }

        // =========================================================================//
        // ======== COSTRUZIONE MASCHERA: QUALI PUNTI ORDINATI SONO VERI? ==========//
        // =========================================================================//
        // confronta i punti ordinati finali
        // con il backup "real_x/y". Segnala alla logica quali di questi 
        // sono veri punti fisici e quali sono il risultato delle predizioni matematiche 
        // calcolate sopra. Essenziale per pilotare la molla cinematica.

        if (num_points_seen == 4) current_point_seen_mask = 0b1111;
        else {
            current_point_seen_mask = 0;

            for (uint8_t i = 0; i < num_points_seen; i++) {
                const int rx = real_x[i];
                const int ry = real_y[i];

                if (positionXX[a] == rx && positionYY[a] == ry) current_point_seen_mask |= 0b1000;
                else if (positionXX[b] == rx && positionYY[b] == ry) current_point_seen_mask |= 0b0100;
                else if (positionXX[c] == rx && positionYY[c] == ry) current_point_seen_mask |= 0b0010;
                else if (positionXX[d] == rx && positionYY[d] == ry) current_point_seen_mask |= 0b0001;
            }
        }

        // OTTIMIZZAZIONE: Unrolling Booleano puro per l'emulatore Shift-Register.
        // Mantiene intatta l'interfaccia verso le API originali Samco che si aspettavano 
        // di leggere cicli progressivi, ma qui lo eseguiamo senza i salti condizionali di un ciclo for.
        see[0] = (see[0] << 1) | ((current_point_seen_mask >> 3) & 1);
        see[1] = (see[1] << 1) | ((current_point_seen_mask >> 2) & 1);
        see[2] = (see[2] << 1) | ((current_point_seen_mask >> 1) & 1);
        see[3] = (see[3] << 1) | (current_point_seen_mask & 1);

        // =========================================================================//
        // ======== MOTORE UNIVERSALE DELLA MOLLA KINEMATICA =======================//
        // =========================================================================//
        // Qui gestiamo l'assorbimento dell'errore visivo. Quando un LED viene perso e 
        // ricalcolato, la nuova coordinata ideale viene caricata con l'errore generato (offset).
        // L'offset viene poi fatto tendere allo zero dolcemente (smorzamento) consumandolo nel tempo.
        int GeomX[4], GeomY[4];
        GeomX[A] = positionXX[a]; GeomY[A] = positionYY[a];
        GeomX[B] = positionXX[b]; GeomY[B] = positionYY[b];
        GeomX[C] = positionXX[c]; GeomY[C] = positionYY[c];
        GeomX[D] = positionXX[d]; GeomY[D] = positionYY[d];

        if (num_points_seen >= 2 && prev_num_points_seen >= 1) {     
            int move_x = 0, move_y = 0;
            uint8_t stable_count = 0;
            for (uint8_t i = 0; i < 4; i++) {
                if ((prev_point_seen_mask & (1 << (3 - i))) && (current_point_seen_mask & (1 << (3 - i)))) {
                    move_x += (GeomX[i] - prev_GeomX[i]); 
                    move_y += (GeomY[i] - prev_GeomY[i]);
                    stable_count++;
                }
            }
            
            int avg_move_x = 0, avg_move_y = 0;
            if (stable_count > 0) {
                avg_move_x = move_x / stable_count; 
                avg_move_y = move_y / stable_count;
            } else if (prev_num_points_seen > 0) {
                avg_move_x = medianX - prev2_medianX;
                avg_move_y = medianY - prev2_medianY;
            }

            bool model_changed = (num_points_seen != prev_num_points_seen);
            for (uint8_t i = 0; i < 4; i++) {
                bool was_seen = (prev_point_seen_mask & (1 << (3 - i))) != 0;
                bool is_seen = (current_point_seen_mask & (1 << (3 - i))) != 0;

                if (is_seen != was_seen || (!is_seen && model_changed)) {
                    int predicted_x = FinalX[i] + avg_move_x;
                    int predicted_y = FinalY[i] + avg_move_y;
                    
                    offset_X[i] = (float)(predicted_x - GeomX[i]);
                    offset_Y[i] = (float)(predicted_y - GeomY[i]);
                }
            }

            int spostamento = abs(avg_move_x) + abs(avg_move_y);
            // FORZIAMO UN DEBITO MINIMO DI 1.0f per evitare il congelamento dell'offset
            //float consumo = (float)spostamento * COSTANTE_MOLLA; 
            float consumo = fmaxf(1.0f, (float)spostamento * COSTANTE_MOLLA);
            
            // OTTIMIZZAZIONE: Branchless Math per la Molla. 
            // Usa le istruzioni hardware min/max per non far spezzare la pipeline della CPU 
            // all'indovinare i salti logici negativi/positivi del debito.
            for (uint8_t i = 0; i < 4; i++) {
                if (offset_X[i] > 0.0f) offset_X[i] = fmaxf(0.0f, offset_X[i] - consumo);
                else if (offset_X[i] < 0.0f) offset_X[i] = fminf(0.0f, offset_X[i] + consumo);
                
                if (offset_Y[i] > 0.0f) offset_Y[i] = fmaxf(0.0f, offset_Y[i] - consumo);
                else if (offset_Y[i] < 0.0f) offset_Y[i] = fminf(0.0f, offset_Y[i] + consumo);
            }
        } else {
            for (uint8_t i = 0; i < 4; i++) { offset_X[i] = 0.0f; offset_Y[i] = 0.0f; }
        }

        // =========================================================================//
        // ======== APPLICAZIONE FINALE E SALVATAGGIO STATO ========================//
        // =========================================================================//
        
        // OTTIMIZZAZIONE: Lookup array per evitare lunghe clausole if/else nidificate.
        // Assegna l'offset alla geometria, spingendo la correzione smorzata sui finali di frame.
        
        // --- Assegnazione Finale e Calcoli Derivati ---
        
        uint8_t assign_map[4] = {a, b, c, d};
        for (uint8_t i = 0; i < 4; i++) {
            FinalX[i] = GeomX[i] + fast_roundf(offset_X[i]);
            FinalY[i] = GeomY[i] + fast_roundf(offset_Y[i]);
            
            positionXX[assign_map[i]] = FinalX[i];
            positionYY[assign_map[i]] = FinalY[i];
            
            prev_GeomX[i] = GeomX[i];
            prev_GeomY[i] = GeomY[i];
        }

        prev2_medianX = medianX;
        prev2_medianY = medianY;

        // Calcolo del centroide con arrotondamento (il +2 è un trucco per l'arrotondamento all'intero più vicino)
        medianX = (FinalX[A] + FinalX[B] + FinalX[C] + FinalX[D] + 2) / 4;
        medianY = (FinalY[A] + FinalY[B] + FinalY[C] + FinalY[D] + 2) / 4;
        
        // Calcoli finali mantenuti per compatibilità
        // OTTIMIZZAZIONE: Sottrazione eseguita sui registri interi *prima* della
        // conversione per non caricare la FPU (Floating Point Unit hardware) di stress doppio.
        float dX_AC = (float)(FinalX[A] - FinalX[C]), dY_AC = (float)(FinalY[A] - FinalY[C]);
        height_left = sqrtf(dX_AC*dX_AC + dY_AC*dY_AC); // lunghezza lato AC
        
        float dX_BD = (float)(FinalX[B] - FinalX[D]), dY_BD = (float)(FinalY[B] - FinalY[D]);
        height_right = sqrtf(dX_BD*dX_BD + dY_BD*dY_BD); // lunghezza lato BD
        
        float dX_AB = (float)(FinalX[A] - FinalX[B]), dY_AB = (float)(FinalY[A] - FinalY[B]);
        width_top = sqrtf(dX_AB*dX_AB + dY_AB*dY_AB); // lunghezza lato AB
        
        float dX_CD = (float)(FinalX[C] - FinalX[D]), dY_CD = (float)(FinalY[C] - FinalY[D]);
        width_bottom = sqrtf(dX_CD*dX_CD + dY_CD*dY_CD); // lunghezza lato CD
        
        // OTTIMIZZAZIONE: In algebra computazionale IEEE-754 la moltiplicazione per 0.5f 
        // è radicalmente più veloce rispetto alla divisione in virgola mobile e produce esatto output.
        height = (height_left + height_right) * 0.5f;
        width = (width_top + width_bottom) * 0.5f;
        
        if (num_points_seen == 4) {
            // Aggiorna la nostra "verità assoluta" sulla forma
            if (height > 1e-6f) { ideal_aspect_ratio = width / height; }
        }

        angle = (atan2f((float)(FinalY[A] - FinalY[B]), (float)(FinalX[B] - FinalX[A])) + 
                 atan2f((float)(FinalY[C] - FinalY[D]), (float)(FinalX[D] - FinalX[C]))) * 0.5f;      

        is_tracking_stable = true;
        prev_point_seen_mask = current_point_seen_mask; 
    }
    else {
        is_tracking_stable = false;
        prev_point_seen_mask = 0;
    }
    prev_num_points_seen = num_points_seen;
}

#endif //USE_SQUARE_ADVANCED