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

// Numerical limit, not the camera field of view: predicted LEDs may be off-screen.
// With |coordinate| <= 16383, every squared 2-D distance used by the ordering
// fits int32_t: 2 * (2 * 16383)^2 < INT32_MAX. Nonnegative costs use uint32_t
// with bounded additions; no 64-bit arithmetic is needed.
static constexpr int SQUARE_COORD_LIMIT = 16383;

static bool SquareCoordinatesValid(const int *x, const int *y) {
    for (uint8_t i = 0; i < 4; ++i) {
        if (x[i] < -SQUARE_COORD_LIMIT || x[i] > SQUARE_COORD_LIMIT ||
            y[i] < -SQUARE_COORD_LIMIT || y[i] > SQUARE_COORD_LIMIT)
            return false;
    }
    return true;
}

// A projected rectangle must remain strictly convex. Also rejects duplicate
// vertices and any three collinear vertices, without imposing a minimum angle
// or aspect ratio that would unnecessarily limit oblique views.
static bool SquareGeometryValid(const int *x, const int *y) {
    if (!SquareCoordinatesValid(x, y))
        return false;

    const uint8_t perimeter[4] = {A, B, D, C};
    int32_t first_turn = 0;
    for (uint8_t i = 0; i < 4; ++i) {
        const uint8_t p = perimeter[i];
        const uint8_t q = perimeter[(i + 1) & 3];
        const uint8_t r = perimeter[(i + 2) & 3];
        const int32_t turn = (x[q] - x[p]) * (y[r] - y[q]) -
                             (y[q] - y[p]) * (x[r] - x[q]);
        if (turn == 0 || (i != 0 && (turn > 0) != (first_turn > 0)))
            return false;
        first_turn = turn;
    }
    return true;
}

// Ricostruzione dei LED nascosti dall'ultima geometria vista con 4 LED, che ne conserva la forma
// reale (prospettiva compresa): una similitudine (n = 2 LED visti) o un'affinita' (n = 3) porta
// i LED visti allora (px, py) su quelli visti ora (qx, qy) e si applica agli m nascosti (hx, hy:
// in ingresso nell'ultima geometria, in uscita ora; validi solo se torna true). Per non deformare
// la forma con la prospettiva, ciascuna delle due viste si guarda da una camera virtuale ruotata
// (rotazione minima) verso la direzione media dei suoi LED visti: se l'arma ha solo ruotato le
// due viste differiscono allora per una rotazione attorno all'asse e la ricostruzione e' esatta.
// Basta una focale approssimata; con focal = 0 si lavora direttamente sull'immagine.
// false anche se un nascosto esce dal limite numerico delle coordinate.
namespace {
struct SquareView {
    float f, ox, oy;                                 // focale e centro dell'immagine
    float r[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};        // rotazione della camera virtuale (per righe)

    SquareView(float focal, float cx, float cy) : f(focal), ox(cx), oy(cy) {}

    bool aim(const float* x, const float* y, uint8_t n) {
        if (f <= 0.0f) return true;
        float mx = 0.0f, my = 0.0f, mz = 0.0f;       // somma delle direzioni (versori)
        for (uint8_t i = 0; i < n; ++i) {
            const float dx = x[i] - ox, dy = y[i] - oy;
            const float k = 1.0f / sqrtf(dx * dx + dy * dy + f * f);
            mx += dx * k; my += dy * k; mz += f * k;
        }
        const float k = 1.0f / sqrtf(mx * mx + my * my + mz * mz);
        mx *= k; my *= k; mz *= k;
        if (!(mz > 0.0f)) return false;
        // rotazione minima che porta (mx, my, mz) sull'asse ottico: asse (my, -mx, 0)
        const float vx = my, vy = -mx, c = 1.0f / (1.0f + mz);
        r[0] = 1.0f - vy * vy * c; r[1] = vx * vy * c;         r[2] = vy;
        r[3] = vx * vy * c;         r[4] = 1.0f - vx * vx * c; r[5] = -vx;
        r[6] = -vy;                 r[7] = vx;                 r[8] = 1.0f - (vx * vx + vy * vy) * c;
        return true;
    }
    // dall'immagine al piano della camera virtuale (unita' mouse) e ritorno
    bool to(float x, float y, float& u, float& v) const {
        const float dx = x - ox, dy = y - oy;
        if (f <= 0.0f) { u = dx; v = dy; return true; }
        const float wz = r[6] * dx + r[7] * dy + r[8] * f;
        if (!(wz > 1e-3f * f)) return false;
        const float k = f / wz;
        u = (r[0] * dx + r[1] * dy + r[2] * f) * k;
        v = (r[3] * dx + r[4] * dy + r[5] * f) * k;
        return true;
    }
    bool from(float u, float v, float& x, float& y) const {
        if (f <= 0.0f) { x = u + ox; y = v + oy; return true; }
        const float wz = r[2] * u + r[5] * v + r[8] * f;
        if (!(wz > 1e-3f * f)) return false;
        const float k = f / wz;
        x = ox + (r[0] * u + r[3] * v + r[6] * f) * k;
        y = oy + (r[1] * u + r[4] * v + r[7] * f) * k;
        return true;
    }
};
}

static bool SquareFromAnchor(float focal, float ox, float oy, uint8_t n,
                             const float* px, const float* py, const float* qx, const float* qy,
                             uint8_t m, float* hx, float* hy) {
    SquareView P(focal, ox, oy), Q(focal, ox, oy);
    if (!P.aim(px, py, n) || !Q.aim(qx, qy, n)) return false;
    float pu[3], pv[3], qu[3], qv[3];
    for (uint8_t i = 0; i < n; ++i)
        if (!P.to(px[i], py[i], pu[i], pv[i]) || !Q.to(qx[i], qy[i], qu[i], qv[i])) return false;
    // Base dalla prima posizione: (P1 - P0, P2 - P0) allora, (Q1 - Q0, Q2 - Q0) ora. Con 2 LED il
    // secondo vettore e' il primo ruotato di 90 gradi (similitudine: scala e rotazione).
    const float a1x = pu[1] - pu[0], a1y = pv[1] - pv[0], b1x = qu[1] - qu[0], b1y = qv[1] - qv[0];
    const float a2x = n == 2 ? -a1y : pu[2] - pu[0], a2y = n == 2 ? a1x : pv[2] - pv[0];
    const float b2x = n == 2 ? -b1y : qu[2] - qu[0], b2y = n == 2 ? b1x : qv[2] - qv[0];
    const float det = a1x * a2y - a1y * a2x;
    if (!(fabsf(det) >= 4.0f)) return false;
    const float inv = 1.0f / det;
    if (n == 3) {
        // Mappa lineare (allora -> ora) p q / r s: un movimento reale dell'arma la lascia vicina a
        // scala e rotazione; un taglio forte vuol dire LED visti con etichette sbagliate.
        const float p = (b1x * a2y - b2x * a1y) * inv, q = (b2x * a1x - b1x * a2x) * inv;
        const float r = (b1y * a2y - b2y * a1y) * inv, t = (b2y * a1x - b1y * a2x) * inv;
        const float cs = (p + t) * (p + t) + (r - q) * (r - q), ac = (p - t) * (p - t) + (q + r) * (q + r);
        if (!(ac <= 0.0625f * cs)) return false;   // parte di taglio fino a 1/4 di quella di scala e rotazione
    }
    for (uint8_t j = 0; j < m; ++j) {
        float u, v;
        if (!P.to(hx[j], hy[j], u, v)) return false;
        const float ex = u - pu[0], ey = v - pv[0];
        const float s = (ex * a2y - ey * a2x) * inv, t = (a1x * ey - a1y * ex) * inv;
        if (!Q.from(qu[0] + s * b1x + t * b2x, qv[0] + s * b1y + t * b2y, hx[j], hy[j])) return false;
        if (!(fabsf(hx[j]) <= SQUARE_COORD_LIMIT && fabsf(hy[j]) <= SQUARE_COORD_LIMIT)) return false;
    }
    return true;
}

void OpenFIRE_Square::configure(const CameraProfile& profile) {
    // Ripartenza pulita (es. dopo il salvataggio, che riavvia la camera): lo storico del
    // tracking precedente non deve sopravvivere alla nuova geometria predefinita.
    start = 0;
    is_tracking_stable = false;
    prev_num_points_seen = 0;
    prev_point_seen_mask = 0;
    current_point_seen_mask = 0;
    for (uint8_t i = 0; i < 4; ++i) {
        motion_X[i] = motion_Y[i] = 0;
        offset_X[i] = offset_Y[i] = 0.0f;
        see[i] = 0;
    }
    spring_fast = false;
    dark_frames = 255;
    anchor_ok = false;

    camMaxX = profile.camMaxX;
    camToMouseShift = profile.camToMouseShift;
    camToMouseMult = profile.camToMouseMult;
    mouseResX = profile.mouseResX;
    mouseResY = profile.mouseResY;
    // Focale nominale in unita' mouse (prospettiva: forma del rettangolo, vedi trueAspectOf, e
    // LED nascosti con 2-3 LED, vedi SquareFromAnchor), dal rapporto focale / semilarghezza del
    // sensore indicato nel profilo. Basta un valore approssimato; 0 se il profilo non lo indica.
    focal = mouseResX * 0.5f * profile.lensFocalRatio;
    FPS_NORMALIZATION = 209.0f / (float)profile.fps;

    FinalX[0] = profile.squareTLX * camToMouseMult;
    FinalX[1] = profile.squareTRX * camToMouseMult;
    FinalX[2] = profile.squareBLX * camToMouseMult;
    FinalX[3] = profile.squareBRX * camToMouseMult;
    FinalY[0] = profile.squareTLY * camToMouseMult;
    FinalY[1] = profile.squareTRY * camToMouseMult;
    FinalY[2] = profile.squareBLY * camToMouseMult;
    FinalY[3] = profile.squareBRY * camToMouseMult;

    medianX = profile.mouseMaxX / 2;
    medianY = profile.mouseMaxY / 2;
    height = (float)(FinalY[2] - FinalY[0]);
    width = (float)(FinalX[1] - FinalX[0]);
    height_left = height;
    height_right = height;
    width_top = width;
    width_bottom = width;
}

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
    if (seenFlags != 0x0F && !start) {
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
            // Validate visible camera coordinates before the signed left shifts.
            if (px[i] < 0 || px[i] > camMaxX || py[i] < 0 ||
                py[i] > ((mouseResY - 1) >> camToMouseShift)) {
                return;
            }
            int calc_x = (camMaxX - px[i]) << camToMouseShift;
            int calc_y = py[i] << camToMouseShift;
            if (calc_x > SQUARE_COORD_LIMIT || calc_y > SQUARE_COORD_LIMIT) {
                return;
            }
            
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
            if (dark_frames < 255) dark_frames++;
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
                const float jump_distance = ((((int)(wideLayout ? height : width) * 3) / 4) * FPS_NORMALIZATION);
                const int max_jump_distance = jump_distance < 46340.0f ? (int)jump_distance : 46340;
                const int32_t MAX_ALLOWED_DISTANCE_SQ = jump_distance < 46340.0f ?
                    (int32_t)max_jump_distance * max_jump_distance : INT32_MAX;

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
            if (!(ideal_aspect_ratio > 0.0f && ideal_aspect_ratio < FLT_MAX)) {
                return;
            }
            const int x1 = positionXX[0], y1 = positionYY[0];
            const int x2 = positionXX[1], y2 = positionYY[1];

            uint8_t best_idx1 = 0, best_idx2 = 1;
            uint32_t min_total_cost = UINT32_MAX;

            // Le 6 possibili combinazioni geometriche dei lati del rettangolo.
            //                             A,B    C,D    A,C    B,D    A,D    B,C
            const uint8_t pairs[6][2] = { {0,1}, {2,3}, {0,2}, {1,3}, {0,3}, {1,2} };

            const int32_t read_dx = x2 - x1;
            const int32_t read_dy = y2 - y1;
            const uint32_t read_dist_sq = read_dx * read_dx + read_dy * read_dy;

            // Keep the original float-to-integer rounding of the side lengths.
            const uint32_t w_sq = (int32_t)(width * width);
            const uint32_t h_sq = (int32_t)(height * height);
            const uint32_t diagonal_sq = w_sq + h_sq;
            const uint32_t expected_dist_arr[6] = {w_sq, w_sq, h_sq, h_sq, diagonal_sq, diagonal_sq};

            // Each observed/previous vertex distance is reused by three pairs.
            // A single squared 2-D distance fits int32_t; two fit uint32_t.
            uint32_t distance1[4], distance2[4], history_penalty[4];
            const bool use_history = prev_num_points_seen >= 1 && prev_point_seen_mask != 0;
            for (uint8_t i = 0; i < 4; ++i) {
                const int32_t dx1 = x1 - FinalX[i], dy1 = y1 - FinalY[i];
                const int32_t dx2 = x2 - FinalX[i], dy2 = y2 - FinalY[i];
                distance1[i] = dx1 * dx1 + dy1 * dy1;
                distance2[i] = dx2 * dx2 + dy2 * dy2;
                history_penalty[i] = use_history && !(prev_point_seen_mask & (1 << (3 - i))) ?
                                     PENALITA_STORICA : 0;
            }

            for (uint8_t i = 0; i < 6; ++i) {
                const uint8_t v1 = pairs[i][0], v2 = pairs[i][1];
                uint32_t penalty = history_penalty[v1] + history_penalty[v2];
                const uint32_t expected = expected_dist_arr[i];
                const uint32_t distance_error = read_dist_sq >= expected ?
                                               read_dist_sq - expected : expected - read_dist_sq;

                // All terms are nonnegative. Test the remaining budget BEFORE
                // adding, so an overflowing cost can never become a false winner.
                if (distance_error >= min_total_cost || penalty >= min_total_cost - distance_error)
                    continue;
                penalty += distance_error;

                const uint32_t costA = distance1[v1] + distance2[v2];
                if (costA < min_total_cost - penalty) {
                    min_total_cost = costA + penalty;
                    best_idx1 = v1;
                    best_idx2 = v2;
                }

                const uint32_t costB = distance1[v2] + distance2[v1];
                if (costB < min_total_cost - penalty) {
                    min_total_cost = costB + penalty;
                    best_idx1 = v2;
                    best_idx2 = v1;
                }
            }

            // Current camera profiles guarantee at least one representable cost.
            // Keep the last valid geometry if a future profile exceeds that range.
            if (min_total_cost == UINT32_MAX)
                return;

            float Vx[4], Vy[4];
            for (uint8_t i = 0; i < 4; i++) { Vx[i] = (float)FinalX[i]; Vy[i] = (float)FinalY[i]; }
            Vx[best_idx1] = (float)x1; Vy[best_idx1] = (float)y1;
            Vx[best_idx2] = (float)x2; Vy[best_idx2] = (float)y2;
    
            const uint8_t idx1 = best_idx1, idx2 = best_idx2;
            bool top    = (idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0);
            bool bottom = (idx1 == 2 && idx2 == 3) || (idx1 == 3 && idx2 == 2);
            bool left   = (idx1 == 0 && idx2 == 2) || (idx1 == 2 && idx2 == 0);
            bool right  = (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 1);

            // Con l'ultima geometria vista con 4 LED i due nascosti ne conservano la forma reale
            // (SquareFromAnchor); altrimenti si costruisce un rettangolo piatto.
            bool from_anchor = false;
            if (anchor_ok) {
                uint8_t hid[2], nh = 0;
                for (uint8_t i = 0; i < 4; i++) if (i != idx1 && i != idx2) hid[nh++] = i;
                const float ax[2] = {(float)anchorX[idx1], (float)anchorX[idx2]};
                const float ay[2] = {(float)anchorY[idx1], (float)anchorY[idx2]};
                const float qx[2] = {(float)x1, (float)x2}, qy[2] = {(float)y1, (float)y2};
                float hx[2] = {(float)anchorX[hid[0]], (float)anchorX[hid[1]]};
                float hy[2] = {(float)anchorY[hid[0]], (float)anchorY[hid[1]]};
                from_anchor = SquareFromAnchor(focal, mouseResX * 0.5f, mouseResY * 0.5f, 2, ax, ay, qx, qy, 2, hx, hy);
                if (from_anchor)
                    for (uint8_t j = 0; j < 2; j++) { Vx[hid[j]] = hx[j]; Vy[hid[j]] = hy[j]; }
            }

            // Per prevenire il "Bug Clessidra" (incrocio errato dei vettori derivati),
            // viene riallineato dinamicamente l'indice dei punti e scalata la proiezione
            // imponendo con la forza matematica il rispetto dell'aspect_ratio originale.
            if (from_anchor) {
                // I due LED nascosti sono gia' ricostruiti.
            } else if (top || bottom) {
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

            for (uint8_t i = 0; i < 4; i++) {
                // Ordered comparisons reject NaN/Inf too, before conversion to int.
                if (!(Vx[i] >= -SQUARE_COORD_LIMIT && Vx[i] <= SQUARE_COORD_LIMIT &&
                      Vy[i] >= -SQUARE_COORD_LIMIT && Vy[i] <= SQUARE_COORD_LIMIT)) {
                    return;
                }
                positionXX[i] = fast_roundf(Vx[i]);
                positionYY[i] = fast_roundf(Vy[i]);
            }
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

        if (!SquareCoordinatesValid(positionXX, positionYY)) {
            return;
        }

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

        // ===========================================================================================
        // ===========================================================================================
        // ===========================================================================================

        // =========================================================================//
        // ============================ CALIBRAZIONE ===============================//
        // =========================================================================//

        if (calibrationMode) {
            // In calibrazione la lightgun è frontale e con rotazione prossima a zero.
            a = orderY[0]; b = orderY[1]; c = orderY[2]; d = orderY[3];
        }

        // =========================================================================//
        // ============================ LAYOUT WIDE ================================//
        // =========================================================================//

        else if (wideLayout) {
            // Keep the existing upper threshold, but scale it down for small targets.
            const int camera_zone = (mouseResX * 5) / 128;
            const int geometry_zone = (positionXX[orderX[3]] - positionXX[orderX[0]]) / 8;
            const int CRITICAL_ZONE_WIDE = geometry_zone < camera_zone ? geometry_zone : camera_zone;

            // distanza del punto più alto dall'estremo destro
            dx = positionXX[orderX[3]] - positionXX[orderY[0]];
            dy = positionYY[orderX[3]] - positionYY[orderY[0]];
            int32_t dist_sq1_wide = (dx * dx) + (dy * dy);

            // distanza del punto più alto dall'estremo sinistro
            dx = positionXX[orderX[0]] - positionXX[orderY[0]];
            dy = positionYY[orderX[0]] - positionYY[orderY[0]];
            int32_t dist_sq2_wide = (dx * dx) + (dy * dy);

            if ((positionXX[orderX[3]] - positionXX[orderX[2]]) > CRITICAL_ZONE_WIDE) {

                if (dist_sq1_wide < dist_sq2_wide) {
                    b = orderY[0];
                    c = orderY[3];

                    if (orderY[1] == orderX[0]) {
                        a = orderY[1];
                        d = orderY[2];
                    } else {
                        d = orderY[1];
                        a = orderY[2];
                    }

                } else {
                    a = orderY[0];
                    d = orderY[3];

                    if (orderY[1] == orderX[0]) {
                        c = orderY[1];
                        b = orderY[2];
                    } else {
                        b = orderY[1];
                        c = orderY[2];
                    }
                }

            } else {
                // Caso Zona Critica
                a = orderX[1];
                b = orderX[3];
                c = orderX[0];
                d = orderX[2];
            }
        }

        // =========================================================================//
        // ============================ LAYOUT NORMALE =============================//
        // =========================================================================//

        else {

            // Assegnazione dei vertici A,B,C,D tramite euristiche
            dx = positionXX[orderY[0]] - positionXX[orderX[0]];
            dy = positionYY[orderY[0]] - positionYY[orderX[0]];
            int32_t dist_sq1 = (dx * dx) + (dy * dy);

            dx = positionXX[orderY[3]] - positionXX[orderX[0]];
            dy = positionYY[orderY[3]] - positionYY[orderX[0]];
            int32_t dist_sq2 = (dx * dx) + (dy * dy);

            /*
            const int CRITICAL_ZONE = (30 * camToMouseMult); // 30 è un valore provato in tante situazioni e pare andare bene
            */

            // Percentuale storica testata sulla DFRobot (30 pixel su una altezza di 768)
            // Calcoliamo il valore dinamico sulla risoluzione della CAM e lo portiamo nello spazio unificato
            const int camera_zone = (mouseResY * 5) / 128;
            const int geometry_zone = (positionYY[orderY[3]] - positionYY[orderY[0]]) / 8;
            const int CRITICAL_ZONE = geometry_zone < camera_zone ? geometry_zone : camera_zone;

            if ((positionYY[orderY[1]] - positionYY[orderY[0]]) > CRITICAL_ZONE) {

                if (dist_sq1 < dist_sq2) {
                    a = orderX[0];
                    d = orderX[3];

                    if (orderX[1] == orderY[3]) {
                        c = orderX[1];
                        b = orderX[2];
                    } else {
                        b = orderX[1];
                        c = orderX[2];
                    }

                } else {
                    c = orderX[0];
                    b = orderX[3];

                    if (orderX[1] == orderY[3]) {
                        d = orderX[1];
                        a = orderX[2];
                    } else {
                        a = orderX[1];
                        d = orderX[2];
                    }
                }

            } else {
                // Caso Zona Critica (rettangolo quasi verticale)
                a = orderY[0];
                b = orderY[1];
                c = orderY[2];
                d = orderY[3];
            }
        }

        // ===========================================================================================
        // ===========================================================================================
        // ===========================================================================================

        // sarebbe sufficiente per la zona critica questo ulteriore controllo, ma poichè ha
        // un piccolo peso computazionale, lo facciamo sempre per casi limiti e come maggiore robustezza
        // Correzione finale per garantire la convenzione (A=TL, B=TR, C=BL, D=BR)
        { uint8_t aux_swap;
            if (positionYY[a] > positionYY[c]) { aux_swap = a; a = c; c = aux_swap; }
            if (positionYY[b] > positionYY[d]) { aux_swap = b; b = d; d = aux_swap; }
            if (positionXX[a] > positionXX[b]) { aux_swap = a; a = b; b = aux_swap; }
            if (positionXX[c] > positionXX[d]) { aux_swap = c; c = d; d = aux_swap; }
        }     

        // Con 1-2 LED la traslazione o la ricostruzione ha gia' assegnato le etichette (A, B, C, D
        // in ordine): con la forma prospettica il riordino dei punti potrebbe scambiarle.
        if (num_points_seen <= 2 && !calibrationMode) { a = A; b = B; c = C; d = D; }

        // =========================================================================//
        // ======== CONTINUITA' DELL'ORDINAMENTO (3-4 PUNTI) =======================//
        // =========================================================================//
        // Con prospettiva forte e rotazione l'euristica puo' alternare le etichette da un
        // frame all'altro, o darne una specchiata (impossibile per la camera), soprattutto
        // con 3 LED dove il quarto punto e' ricostruito. Il perimetro A, B, D, C si ricava
        // dall'angolo attorno al centroide, che e' sempre orientato bene; resta da scegliere
        // quale punto e' A:
        // - in tracking: la rotazione ciclica piu' vicina alla geometria precedente (solo
        //   sui LED visti). Si usa se l'euristica e' specchiata, oppure se e' nettamente
        //   migliore (costo < 1/4) e la forma (base/altezza come nel layout) la conferma o,
        //   se la forma non decide, non implica una rotazione dell'arma di oltre 80 gradi
        //   superiore a quella dell'euristica: cosi' un aggancio sbagliato (es. partenza
        //   con l'arma molto ruotata) si corregge. La forma puo' anche dare ragione
        //   all'euristica: con 4 LED veri se lo conferma la forma vera del rettangolo, corretta
        //   dalla prospettiva (con l'arma molto di lato, sopra o sotto lo schermo la forma vista
        //   si puo' invertire); con 3 LED solo nel primo secondo circa dopo una partenza o
        //   ripresa (aggancio ancora incerto) e se la continuita' implica l'arma piu' ruotata;
        // - dopo un buio breve, come in tracking se i 4 LED sono vicini all'ultima geometria;
        // - alla partenza, solo se l'euristica e' specchiata: la rotazione con l'arma piu' dritta.
        // In calibrazione non si applica.
        if (num_points_seen >= 3 && !calibrationMode) {
            const uint8_t heur[4] = {a, b, d, c};
            float cxm = 0.0f, cym = 0.0f;
            for (uint8_t i = 0; i < 4; i++) { cxm += positionXX[i]; cym += positionYY[i]; }
            cxm *= 0.25f; cym *= 0.25f;
            float an[4]; uint8_t ord[4] = {0, 1, 2, 3};
            for (uint8_t i = 0; i < 4; i++) an[i] = atan2f(positionYY[i] - cym, positionXX[i] - cxm);
            for (uint8_t i = 1; i < 4; i++)
                for (uint8_t j = i; j > 0 && an[ord[j]] < an[ord[j - 1]]; j--) { const uint8_t t = ord[j]; ord[j] = ord[j - 1]; ord[j - 1] = t; }
            // Area con segno dell'etichettatura dell'euristica: negativa = specchiata
            const float area_h = (float)(positionXX[d] - positionXX[a]) * (float)(positionYY[c] - positionYY[b]) -
                                 (float)(positionYY[d] - positionYY[a]) * (float)(positionXX[c] - positionXX[b]);
            // Rotazione apparente di un'etichettatura (p = A, B, D, C), dai lati verticali A-C e B-D:
            // con la prospettiva convergono verso il loro punto di fuga, e la direzione dal centro
            // dell'immagine al punto di fuga e' la rotazione dell'arma attorno alla canna qualunque
            // sia l'angolo di mira (con i lati paralleli e' la loro direzione). E' esatta solo per
            // una lente ideale: con la distorsione e il punto ricostruito (3 LED) e' approssimata,
            // ma il margine di 80 gradi del controllo la tollera.
            auto rollOf = [&](const uint8_t* p) {
                const float ox = mouseResX * 0.5f, oy = mouseResY * 0.5f;   // centro dell'immagine
                const float ax = positionXX[p[0]] - ox, ay = positionYY[p[0]] - oy;
                const float bx = positionXX[p[1]] - ox, by = positionYY[p[1]] - oy;
                const float dx = positionXX[p[2]] - ox, dy = positionYY[p[2]] - oy;
                const float cx = positionXX[p[3]] - ox, cy = positionYY[p[3]] - oy;
                // rette A-C e B-D in coordinate omogenee e loro intersezione (punto di fuga)
                const float l1x = ay - cy, l1y = cx - ax, l1z = ax * cy - ay * cx;
                const float l2x = by - dy, l2y = dx - bx, l2z = bx * dy - by * dx;
                float vx = l1y * l2z - l1z * l2y, vy = l1z * l2x - l1x * l2z;
                const float sx = (cx - ax) + (dx - bx), sy = (cy - ay) + (dy - by);   // da A, B verso C, D
                if (vx * sx + vy * sy < 0.0f) { vx = -vx; vy = -vy; }
                if (fabsf(vx) + fabsf(vy) < 1.0f) { vx = sx; vy = sy; }   // caso degenere
                return fabsf(atan2f(-vx, vy));
            };
            const uint8_t labs[4] = {A, B, D, C};
            auto costOf = [&](const uint8_t* p) {
                uint32_t cst = 0;
                for (uint8_t j = 0; j < 4; ++j) {
                    if (p[j] >= num_points_seen) continue;  // punto ricostruito (solo con 3 LED)
                    const int32_t ddx = positionXX[p[j]] - prev_GeomX[labs[j]];
                    const int32_t ddy = positionYY[p[j]] - prev_GeomY[labs[j]];
                    const uint32_t q = (uint32_t)(ddx * ddx + ddy * ddy);
                    cst = (cst > UINT32_MAX - q) ? UINT32_MAX : cst + q;
                }
                return cst;
            };
            // Buio breve (al massimo 8 fotogrammi, es. un LED che sfarfalla o una mano che passa):
            // con i 4 LED di nuovo visibili si resta in continuita' se la rotazione ciclica piu'
            // vicina all'ultima geometria dista poco (somma dei quadrati < 1/16 del lato piu' corto
            // al quadrato); altrimenti si riparte come alla partenza.
            bool resume = false;
            if (!is_tracking_stable && num_points_seen == 4 && dark_frames <= 8) {
                uint32_t best = UINT32_MAX; uint8_t bk = 0;
                for (uint8_t k = 0; k < 4; ++k) {
                    const uint8_t cand[4] = {ord[k], ord[(k + 1) & 3], ord[(k + 2) & 3], ord[(k + 3) & 3]};
                    const uint32_t cst = costOf(cand);
                    if (cst < best) { best = cst; bk = k; }
                }
                float minSide2 = FLT_MAX;
                for (uint8_t j = 0; j < 4; ++j) {
                    const float dx = (float)(positionXX[ord[(bk + j) & 3]] - positionXX[ord[(bk + j + 1) & 3]]);
                    const float dy = (float)(positionYY[ord[(bk + j) & 3]] - positionYY[ord[(bk + j + 1) & 3]]);
                    minSide2 = fminf(minSide2, dx * dx + dy * dy);
                }
                resume = (float)best < minSide2 * 0.0625f;
            }
            int chosen = -1;
            const bool continuity = (is_tracking_stable && prev_num_points_seen >= 3) || resume;
            if (!continuity) frames_since_cold = 0;
            else if (frames_since_cold < 255) frames_since_cold++;
            if (continuity) {
                uint32_t best = UINT32_MAX; uint8_t bk = 0;
                for (uint8_t k = 0; k < 4; ++k) {
                    const uint8_t cand[4] = {ord[k], ord[(k + 1) & 3], ord[(k + 2) & 3], ord[(k + 3) & 3]};
                    const uint32_t cst = costOf(cand);
                    if (cst < best) { best = cst; bk = k; }
                }
                const uint8_t cand[4] = {ord[bk], ord[(bk + 1) & 3], ord[(bk + 2) & 3], ord[(bk + 3) & 3]};
                const bool same = cand[0] == heur[0] && cand[1] == heur[1] && cand[2] == heur[2] && cand[3] == heur[3];
                if (!same) {
                    if (area_h <= 0.0f) chosen = bk;
                    else if (best < costOf(heur) / 4) {
                        // Prima decide la forma, se e' chiara: con le etichette giuste la base e l'altezza
                        // stanno come nel layout (wideLayout), con un aggancio ruotato di 90 gradi sono
                        // invertite. Cosi' l'arma inclinata di 85-90 gradi non scambia i LED, e un aggancio
                        // sbagliato si corregge anche con l'arma inclinata di 20-45 gradi.
                        // Forma: (AB^2 + CD^2 - AC^2 - BD^2) / somma, positiva = base piu' lunga dell'altezza.
                        auto aspectOf = [&](const uint8_t* p) {
                            auto d2 = [&](uint8_t i, uint8_t j) {
                                const float dx = (float)(positionXX[i] - positionXX[j]), dy = (float)(positionYY[i] - positionYY[j]);
                                return dx * dx + dy * dy;
                            };
                            const float h = d2(p[1], p[0]) + d2(p[2], p[3]), v = d2(p[3], p[0]) + d2(p[2], p[1]);
                            const float q = (h - v) / fmaxf(h + v, 1.0f);
                            return wideLayout ? q : -q;
                        };
                        // Con l'arma molto di lato, o molto sopra o sotto lo schermo, la prospettiva accorcia
                        // un lato fino a invertire la forma vista. Con 4 LED veri si ricava allora la forma
                        // vera del rettangolo, sul piano dello schermo, dalla prospettiva dei 4 punti e dalla
                        // focale: nello spazio D - A = (B - A) + (C - A), da cui le profondita' di B e di C
                        // rispetto ad A (rapporti di aree nell'immagine) e quindi i lati veri A-B e A-C.
                        // Stessa scala e stesso segno di aspectOf; 0 se il quadrilatero e' degenere. Basta
                        // valutarla per la continuita': l'ordine dell'euristica puo' essere incrociato.
                        // Se il profilo non indica la focale si usa la regola dei 3 LED.
                        auto trueAspectOf = [&](const uint8_t* p) {
                            const float ox = mouseResX * 0.5f, oy = mouseResY * 0.5f;   // centro dell'immagine
                            const float xa = positionXX[p[0]] - ox, ya = positionYY[p[0]] - oy;
                            const float xb = positionXX[p[1]] - ox, yb = positionYY[p[1]] - oy;
                            const float xd = positionXX[p[2]] - ox, yd = positionYY[p[2]] - oy;
                            const float xc = positionXX[p[3]] - ox, yc = positionYY[p[3]] - oy;
                            // doppia area con segno del triangolo (0, 1, 2)
                            auto area2 = [](float x0, float y0, float x1, float y1, float x2, float y2) {
                                return (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
                            };
                            const float sB = area2(xb, yb, xd, yd, xc, yc), sC = area2(xc, yc, xd, yd, xb, yb);
                            if (fabsf(sB) < 1.0f || fabsf(sC) < 1.0f) return 0.0f;
                            const float kB = area2(xa, ya, xd, yd, xc, yc) / sB;   // profondita' di B / profondita' di A
                            const float kC = area2(xa, ya, xd, yd, xb, yb) / sC;   // profondita' di C / profondita' di A
                            if (!(kB > 0.0f && kC > 0.0f)) return 0.0f;
                            const float hx = kB * xb - xa, hy = kB * yb - ya, hz = focal * (kB - 1.0f);   // lato A-B
                            const float vx = kC * xc - xa, vy = kC * yc - ya, vz = focal * (kC - 1.0f);   // lato A-C
                            const float h = hx * hx + hy * hy + hz * hz, v = vx * vx + vy * vy + vz * vz;
                            const float q = (h - v) / fmaxf(h + v, 1.0f);
                            return wideLayout ? q : -q;
                        };
                        const float q_cand = aspectOf(cand), q_heur = aspectOf(heur);
                        bool ok;
                        if (q_cand > 0.1f && q_heur < -0.1f) ok = true;         // la forma conferma la continuita'
                        else if (q_cand < -0.1f && q_heur > 0.1f &&
                                 (num_points_seen == 4 && focal > 0.0f ? trueAspectOf(cand) < -0.1f
                                                       : (frames_since_cold < 200 && rollOf(cand) > rollOf(heur)))) ok = false;   // la forma conferma l'euristica
                        else ok = rollOf(cand) < rollOf(heur) + 80.0f * 0.0174533f;  // forma incerta: rotazione
                        if (ok) chosen = bk;
                    }
                }
            } else if (area_h < 0.0f) {
                float bestRoll = 10.0f;
                for (uint8_t k = 0; k < 4; ++k) {
                    const uint8_t cand[4] = {ord[k], ord[(k + 1) & 3], ord[(k + 2) & 3], ord[(k + 3) & 3]};
                    const float r = rollOf(cand);
                    if (r < bestRoll) { bestRoll = r; chosen = k; }
                }
            }
            if (chosen >= 0) { a = ord[chosen]; b = ord[(chosen + 1) & 3]; d = ord[(chosen + 2) & 3]; c = ord[(chosen + 3) & 3]; }
        }

        // =========================================================================//
        // ======== COSTRUZIONE MASCHERA: QUALI PUNTI ORDINATI SONO VERI? ==========//
        // =========================================================================//
        // confronta i punti ordinati finali
        // con il backup "real_x/y". Segnala alla logica quali di questi 
        // sono veri punti fisici e quali sono il risultato delle predizioni matematiche 
        // calcolate sopra. Essenziale per pilotare la molla cinematica.

        uint8_t next_point_seen_mask;
        if (num_points_seen == 4) next_point_seen_mask = 0b1111;
        else {
            next_point_seen_mask = 0;

            for (uint8_t i = 0; i < num_points_seen; i++) {
                const int rx = real_x[i];
                const int ry = real_y[i];

                if (positionXX[a] == rx && positionYY[a] == ry) next_point_seen_mask |= 0b1000;
                else if (positionXX[b] == rx && positionYY[b] == ry) next_point_seen_mask |= 0b0100;
                else if (positionXX[c] == rx && positionYY[c] == ry) next_point_seen_mask |= 0b0010;
                else if (positionXX[d] == rx && positionYY[d] == ry) next_point_seen_mask |= 0b0001;
            }
        }

        // MIGLIORAMENTO DELL POSIZIONE DEL QUARTO PUNTO QUANTO SI VEDONO SOLO 3 PUNTI
        // Con l'ultima geometria vista con 4 LED il quarto punto ne conserva la forma reale
        // (SquareFromAnchor); altrimenti l'iniezione prospettica storica qui sotto.
        bool hidden_from_anchor = false;
        if (num_points_seen == 3 && anchor_ok) {
            const uint8_t missing_mask = (~next_point_seen_mask) & 0x0F;
            const uint8_t lab[4] = {a, b, c, d};   // indice in positionXX di A, B, C, D
            uint8_t hl = 0, vl[3], nv = 0, nh = 0;
            for (uint8_t L = 0; L < 4; L++) {
                if (missing_mask & (1 << (3 - L))) { hl = L; nh++; }
                else if (nv < 3) vl[nv++] = L;
            }
            if (nh == 1) {
                float ax[3], ay[3], qx[3], qy[3];
                for (uint8_t k = 0; k < 3; k++) {
                    ax[k] = (float)anchorX[vl[k]]; ay[k] = (float)anchorY[vl[k]];
                    qx[k] = (float)positionXX[lab[vl[k]]]; qy[k] = (float)positionYY[lab[vl[k]]];
                }
                float hx = (float)anchorX[hl], hy = (float)anchorY[hl];
                if (SquareFromAnchor(focal, mouseResX * 0.5f, mouseResY * 0.5f, 3, ax, ay, qx, qy, 1, &hx, &hy)) {
                    positionXX[lab[hl]] = fast_roundf(hx);
                    positionYY[lab[hl]] = fast_roundf(hy);
                    hidden_from_anchor = true;
                }
            }
        }

        // =========================================================================//
        // === (3 PUNTI): INIEZIONE PROSPETTICA STORICA (Zero-Cost Math) ====//
        // =========================================================================//
        if (!hidden_from_anchor && num_points_seen == 3 && prev_num_points_seen >= 3) {
            uint8_t missing_mask = (~next_point_seen_mask) & 0x0F;

            // Manca D (BR). I punti certi sono A, B, C.
            if (missing_mask == 0b0001) { 
                int affine_dx = positionXX[b] + positionXX[c] - positionXX[a];
                int affine_dy = positionYY[b] + positionYY[c] - positionYY[a];
                int prev_affine_dx = prev_GeomX[B] + prev_GeomX[C] - prev_GeomX[A];
                int prev_affine_dy = prev_GeomY[B] + prev_GeomY[C] - prev_GeomY[A];
                positionXX[d] = affine_dx + (prev_GeomX[D] - prev_affine_dx);
                positionYY[d] = affine_dy + (prev_GeomY[D] - prev_affine_dy);
            }
            // Manca C (BL). I punti certi sono A, B, D.
            else if (missing_mask == 0b0010) {
                int affine_cx = positionXX[a] + positionXX[d] - positionXX[b];
                int affine_cy = positionYY[a] + positionYY[d] - positionYY[b];
                int prev_affine_cx = prev_GeomX[A] + prev_GeomX[D] - prev_GeomX[B];
                int prev_affine_cy = prev_GeomY[A] + prev_GeomY[D] - prev_GeomY[B];
                positionXX[c] = affine_cx + (prev_GeomX[C] - prev_affine_cx);
                positionYY[c] = affine_cy + (prev_GeomY[C] - prev_affine_cy);
            }
            // Manca B (TR). I punti certi sono A, C, D.
            else if (missing_mask == 0b0100) {
                int affine_bx = positionXX[a] + positionXX[d] - positionXX[c];
                int affine_by = positionYY[a] + positionYY[d] - positionYY[c]; 
                int prev_affine_bx = prev_GeomX[A] + prev_GeomX[D] - prev_GeomX[C];
                int prev_affine_by = prev_GeomY[A] + prev_GeomY[D] - prev_GeomY[C];
                positionXX[b] = affine_bx + (prev_GeomX[B] - prev_affine_bx);
                positionYY[b] = affine_by + (prev_GeomY[B] - prev_affine_by);
            }
            // Manca A (TL). I punti certi sono B, C, D.
            else if (missing_mask == 0b1000) {
                int affine_ax = positionXX[b] + positionXX[c] - positionXX[d];
                int affine_ay = positionYY[b] + positionYY[c] - positionYY[d];
                int prev_affine_ax = prev_GeomX[B] + prev_GeomX[C] - prev_GeomX[D];
                int prev_affine_ay = prev_GeomY[B] + prev_GeomY[C] - prev_GeomY[D];
                positionXX[a] = affine_ax + (prev_GeomX[A] - prev_affine_ax);
                positionYY[a] = affine_ay + (prev_GeomY[A] - prev_affine_ay);
            }
        }
        
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

        if (!SquareGeometryValid(GeomX, GeomY)) {
            // Retry only four measured vertices. Partial tracking and calibration
            // keep their existing rejection behavior.
            if (num_points_seen != 4 || calibrationMode)
                return;

            // Reuse the existing Y sort, ordering left/right inside each pair.
            // Accept this fallback only when its axes agree with the layout.
            a = orderY[0]; b = orderY[1]; c = orderY[2]; d = orderY[3];
            if (positionXX[a] > positionXX[b]) { uint8_t t = a; a = b; b = t; }
            if (positionXX[c] > positionXX[d]) { uint8_t t = c; c = d; d = t; }

            // Sign of (AB² + CD²) - (AC² + BD²), without computing
            // four lengths: half of that difference is (B-C) dot (D-A).
            // The coordinate limit keeps this signed 32-bit dot product safe.
            const int32_t dx_bc = positionXX[b] - positionXX[c];
            const int32_t dy_bc = positionYY[b] - positionYY[c];
            const int32_t dx_da = positionXX[d] - positionXX[a];
            const int32_t dy_da = positionYY[d] - positionYY[a];
            const int32_t axis_score = dx_bc * dx_da + dy_bc * dy_da;
            if (axis_score == 0 || (axis_score > 0) != wideLayout)
                return;

            GeomX[A] = positionXX[a]; GeomY[A] = positionYY[a];
            GeomX[B] = positionXX[b]; GeomY[B] = positionYY[b];
            GeomX[C] = positionXX[c]; GeomY[C] = positionYY[c];
            GeomX[D] = positionXX[d]; GeomY[D] = positionYY[d];
            if (!SquareGeometryValid(GeomX, GeomY))
                return;
        }

        // From here every path commits either smoothed or pure valid geometry.
        // No later rejection may leave partly updated smoothing state behind.

        if (num_points_seen >= 2 && prev_num_points_seen >= 1) {     
            int move_x = 0, move_y = 0;
            uint8_t stable_count = 0;
            for (uint8_t i = 0; i < 4; i++) {
                if ((prev_point_seen_mask & (1 << (3 - i))) && (next_point_seen_mask & (1 << (3 - i)))) {
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
                bool is_seen = (next_point_seen_mask & (1 << (3 - i))) != 0;

                if (is_seen != was_seen || (!is_seen && model_changed)) {
                    int predicted_x = FinalX[i] + avg_move_x;
                    int predicted_y = FinalY[i] + avg_move_y;
                    
                    offset_X[i] = (float)(predicted_x - GeomX[i]);
                    offset_Y[i] = (float)(predicted_y - GeomY[i]);
                    spring_fast = num_points_seen > prev_num_points_seen;
                }
            }

            int spostamento = abs(avg_move_x) + abs(avg_move_y);
            // Default consumption for partial visibility; four-LED recovery
            // below uses its filtered motion estimate instead.
            float consumo = fmaxf(FPS_NORMALIZATION, (float)spostamento * COSTANTE_MOLLA);

            // With four measured LEDs, blend smoothly at low motion and
            // converge quickly at high motion. With 1-3 LEDs the return stays
            // linear (slow for large offsets), except after LEDs were added:
            // then it is proportional at the slowest rate (MIN_RETURN), because
            // the old geometry was less reliable than the new one.
            if (num_points_seen == 4 || spring_fast) {
                float largest_offset = 0.0f;
                for (uint8_t i = 0; i < 4; ++i) {
                    if (offset_X[i] != 0.0f)
                        largest_offset = fmaxf(largest_offset, fabsf(offset_X[i]));
                    if (offset_Y[i] != 0.0f)
                        largest_offset = fmaxf(largest_offset, fabsf(offset_Y[i]));
                }

                if (largest_offset > 0.0f) {
                    // Estimate motion only while recovering four measured LEDs.
                    // Reset when the number of seen LEDs changes: old velocity must
                    // not speed up a new, slow return. Newly visible LEDs start from zero.
                    int32_t motion_sum = 0;
                    for (uint8_t i = 0; i < 4; ++i) {
                        if (prev_num_points_seen != num_points_seen)
                            motion_X[i] = motion_Y[i] = 0;
                        if (prev_point_seen_mask & (1 << (3 - i))) {
                            const int32_t dx = GeomX[i] - prev_GeomX[i];
                            const int32_t dy = GeomY[i] - prev_GeomY[i];
                            motion_X[i] = (motion_X[i] + dx) / 2;
                            motion_Y[i] = (motion_Y[i] + dy) / 2;
                            motion_sum += abs(motion_X[i]) + abs(motion_Y[i]);
                        }
                    }

                    // Average each LED's signed motion over time BEFORE taking magnitudes:
                    // suppress jitter without cancelling rotation or zoom.
                    // Integer rounding here affects only the motion estimate,
                    // not the measured geometry or the fractional spring offsets.
                    const float movement = stable_count ?
                        (float)motion_sum / (float)stable_count : (float)spostamento;
                    const float MIN_RETURN = 0.03125f;
                    const float MAX_RETURN = 0.5f;
                    const float MOTION_GAIN = 32.0f;
                    const float motion = MOTION_GAIN * movement;
                    const float scale = fmaxf(height, 1.0f) * FPS_NORMALIZATION;

                    float return_rate;
                    if (num_points_seen != 4 || motion <= scale * MIN_RETURN)
                        return_rate = MIN_RETURN;
                    else if (motion >= scale * MAX_RETURN)
                        return_rate = MAX_RETURN;
                    else
                        return_rate = motion / scale;

                    consumo = fmaxf(0.5f * FPS_NORMALIZATION,
                        fmaxf(movement * COSTANTE_MOLLA, largest_offset * return_rate));
                }
            }

            // Molla: tutti gli offset vengono ridotti dello stesso fattore, cosi' la
            // componente piu' grande cala di "consumo" e la forma della correzione resta
            // la stessa. Sottrarre "consumo" a ogni componente azzerava prima quelle
            // piccole, deformando il quadrilatero e causando scatti del cursore.
            {
                float largest = 0.0f;
                for (uint8_t i = 0; i < 4; i++) {
                    largest = fmaxf(largest, fabsf(offset_X[i]));
                    largest = fmaxf(largest, fabsf(offset_Y[i]));
                }
                const float factor = largest > consumo ? 1.0f - consumo / largest : 0.0f;
                for (uint8_t i = 0; i < 4; i++) {
                    offset_X[i] *= factor;
                    offset_Y[i] *= factor;
                }
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
        
        int nextFinalX[4], nextFinalY[4];
        bool geometry_changed = false;
        for (uint8_t i = 0; i < 4; i++) {
            // Offsets are bounded combinations of validated integer coordinates.
            const int offset_x = fast_roundf(offset_X[i]);
            const int offset_y = fast_roundf(offset_Y[i]);
            nextFinalX[i] = GeomX[i] + offset_x;
            nextFinalY[i] = GeomY[i] + offset_y;
            geometry_changed |= offset_x != 0 || offset_y != 0;
        }
        // GeomX/Y already passed validation. Recheck only if smoothing changed them.
        if (geometry_changed && !SquareGeometryValid(nextFinalX, nextFinalY)) {
            // The measured/reconstructed geometry is valid; only the smoothing
            // offsets made it unsafe. Use that geometry without smoothing for
            // this update. Holding the old offsets here could reject the same
            // valid four-LED frame forever after reacquisition.
            for (uint8_t i = 0; i < 4; ++i) {
                offset_X[i] = 0.0f;
                offset_Y[i] = 0.0f;
                nextFinalX[i] = GeomX[i];
                nextFinalY[i] = GeomY[i];
            }
        }

        // Publish only after the entire candidate geometry passed validation.
        for (uint8_t i = 0; i < 4; i++) {
            FinalX[i] = nextFinalX[i];
            FinalY[i] = nextFinalY[i];
            
            prev_GeomX[i] = GeomX[i];
            prev_GeomY[i] = GeomY[i];
        }

        current_point_seen_mask = next_point_seen_mask;
        // Update visibility history only when accepting the new geometry.
        see[0] = (see[0] << 1) | ((next_point_seen_mask >> 3) & 1);
        see[1] = (see[1] << 1) | ((next_point_seen_mask >> 2) & 1);
        see[2] = (see[2] << 1) | ((next_point_seen_mask >> 1) & 1);
        see[3] = (see[3] << 1) | (next_point_seen_mask & 1);

        if (num_points_seen == 4) {
            start = 0xFF;
            for (uint8_t i = 0; i < 4; i++) { anchorX[i] = GeomX[i]; anchorY[i] = GeomY[i]; }
            anchor_ok = true;
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
        prev_point_seen_mask = next_point_seen_mask; 
        dark_frames = num_points_seen >= 3 ? 0 : 255;
    }
    else {
        is_tracking_stable = false;
        prev_point_seen_mask = 0;
        if (dark_frames < 255) dark_frames++;
        anchor_ok = false;   // dopo il buio la posizione dell'arma rispetto allo schermo puo' essere cambiata
    }
    prev_num_points_seen = num_points_seen;
}

#endif //USE_SQUARE_ADVANCED
