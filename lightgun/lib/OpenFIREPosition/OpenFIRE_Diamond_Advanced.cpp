#ifdef USE_DIAMOND_ADVANCED
/*!
 * @file OpenFIRE_Diamond_Advanced.cpp
 * @brief Light Gun library for 4 LED setup
 * @n CPP file for Samco Light Gun 4 LED setup (Diamond layout)
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
#include "OpenFIRE_Diamond_Advanced.h"
#include <cfloat>

void OpenFIRE_Diamond::configure(const CameraProfile& profile)
{
    camToMouseShift = profile.camToMouseShift;
    camToMouseMult = profile.camToMouseMult;
    FPS_NORMALIZATION = profile.fps ? 209.0f / (float)profile.fps : 1.0f;

    FinalX[0] = profile.diamondTCX * camToMouseMult;
    FinalX[1] = profile.diamondRCX * camToMouseMult;
    FinalX[2] = profile.diamondBCX * camToMouseMult;
    FinalX[3] = profile.diamondLCX * camToMouseMult;
    FinalY[0] = profile.diamondTCY * camToMouseMult;
    FinalY[1] = profile.diamondRCY * camToMouseMult;
    FinalY[2] = profile.diamondBCY * camToMouseMult;
    FinalY[3] = profile.diamondLCY * camToMouseMult;

    medianX = profile.mouseMaxX / 2;
    medianY = profile.mouseMaxY / 2;

    {
        const float hd = (float)(FinalX[1] - FinalX[3]);
        const float vd = (float)(FinalY[2] - FinalY[0]);
        diagRatioRef = (hd > 0.0f && vd > 0.0f) ? hd / vd : 1.0f;
    }

    start = 0;
    prevNumSeen = 0;
    darkFrames = 255;
    prevSeenMask = 0;
    for (uint8_t i = 0; i < 4; i++) {
        prevGeomX[i] = (float)FinalX[i];
        prevGeomY[i] = (float)FinalY[i];
        offset_X[i] = offset_Y[i] = 0.0f;
        motion_X[i] = motion_Y[i] = 0.0f;
        see[i] = 0;
    }
    spring_fast = false;
    angle = 0;
    height2 = width2 = 0;
}

// Vero se il ciclo 0-1-2-3 e' un quadrilatero convesso (tutte le svolte nello stesso verso).
bool OpenFIRE_Diamond::ConvexCycle(const float* x, const float* y)
{
    float first = 0.0f;
    for (uint8_t i = 0; i < 4; i++) {
        const uint8_t p = i, q = (i + 1) & 3, r = (i + 2) & 3;
        const float turn = (x[q] - x[p]) * (y[r] - y[q]) - (y[q] - y[p]) * (x[r] - x[q]);
        if (turn == 0.0f || (i && ((turn > 0.0f) != (first > 0.0f)))) return false;
        first = turn;
    }
    return true;
}

void OpenFIRE_Diamond::begin(const int* px, const int* py, unsigned int seen)
{
    seenFlags = seen;

    // Punti visti, nello spazio del mouse
    float sx[4], sy[4];
    uint8_t n = 0;
    for (uint8_t i = 0; i < 4; i++) {
        if (seen & (1u << i)) {
            sx[n] = (float)(px[i] << camToMouseShift);
            sy[n] = (float)(py[i] << camToMouseShift);
            n++;
        }
    }

    // Si parte solo dopo aver visto tutti e 4 i LED
    if (!start && n != 4) return;
    // Buio: si mantiene l'ultima posizione e si ricorda che la geometria non e' piu' aggiornata
    if (n == 0) { prevNumSeen = 0; if (darkFrames < 255) darkFrames++; return; }
    // Dopo il buio servono almeno 3 LED: con 1-2 l'identificazione non e' affidabile
    if (prevNumSeen == 0 && n < 3) { if (darkFrames < 255) darkFrames++; return; }

    float Gx[4], Gy[4];
    uint8_t mask = 0;

    if (n == 4) {
        // =========================================================================//
        // ======== 4 LED: ORDINAMENTO CICLICO ATTORNO AL CENTROIDE ================//
        // =========================================================================//
        // L'angolo cresce alto -> X max -> basso -> X min (Y verso il basso), qualunque
        // sia la distanza o la prospettiva. Resta solo da scegliere quale punto e' "alto".
        const float cx = (sx[0] + sx[1] + sx[2] + sx[3]) * 0.25f;
        const float cy = (sy[0] + sy[1] + sy[2] + sy[3]) * 0.25f;
        uint8_t ord[4] = {0, 1, 2, 3};
        float an[4];
        for (uint8_t i = 0; i < 4; i++) an[i] = atan2f(sy[i] - cy, sx[i] - cx);
        for (uint8_t i = 1; i < 4; i++)
            for (uint8_t j = i; j > 0 && an[ord[j]] < an[ord[j - 1]]; j--) {
                const uint8_t t = ord[j]; ord[j] = ord[j - 1]; ord[j - 1] = t;
            }

        // Regola di partenza (o ripresa dopo il buio). Il rombo e' piu' largo che alto (segue
        // il monitor), quindi la diagonale piu' lunga e' quella orizzontale anche con l'arma
        // molto ruotata; "alto" e' l'estremo della diagonale verticale con Y minore.
        // Se le due diagonali sono simili (rombo quasi quadrato o prospettiva forte)
        // si usa il punto con Y minima, come prima.
        uint8_t cold = 0;
        bool coldShape = false;  // true se la regola di partenza ha deciso con le diagonali
        {
            const float ax = sx[ord[0]] - sx[ord[2]], ay = sy[ord[0]] - sy[ord[2]];
            const float bx = sx[ord[1]] - sx[ord[3]], by = sy[ord[1]] - sy[ord[3]];
            const float da = ax * ax + ay * ay, db = bx * bx + by * by;
            const bool wide = diagRatioRef > 1.15f, tall = diagRatioRef < 0.87f;
            if ((wide || tall) && (da > 1.3225f * db || db > 1.3225f * da)) {
                // verticale = la diagonale piu' corta (rombo largo) o piu' lunga (rombo alto)
                const bool vertIsA = wide ? (da < db) : (da > db);
                const uint8_t v0 = vertIsA ? 0 : 1;
                cold = (sy[ord[v0]] <= sy[ord[v0 + 2]]) ? v0 : (uint8_t)(v0 + 2);
                coldShape = true;
            } else {
                for (uint8_t r = 1; r < 4; r++) if (sy[ord[r]] < sy[ord[cold]]) cold = r;
            }
        }

        uint8_t best = cold;
        // Buio breve (al massimo 8 fotogrammi): si resta in continuita' se la rotazione ciclica
        // piu' vicina all'ultima posizione mostrata dista poco (somma dei quadrati < 1/16 del
        // lato piu' corto al quadrato); altrimenti vale la regola di partenza.
        bool resume = false;
        if (prevNumSeen == 0 && start && darkFrames <= 8) {
            float bestCost = FLT_MAX;
            uint8_t r0 = 0;
            for (uint8_t r = 0; r < 4; r++) {
                float cost = 0.0f;
                for (uint8_t k = 0; k < 4; k++) {
                    const float dx = sx[ord[(r + k) & 3]] - (float)FinalX[k];
                    const float dy = sy[ord[(r + k) & 3]] - (float)FinalY[k];
                    cost += dx * dx + dy * dy;
                }
                if (cost < bestCost) { bestCost = cost; r0 = r; }
            }
            float minSide2 = FLT_MAX;
            for (uint8_t k = 0; k < 4; k++) {
                const float dx = sx[ord[(r0 + k) & 3]] - sx[ord[(r0 + k + 1) & 3]];
                const float dy = sy[ord[(r0 + k) & 3]] - sy[ord[(r0 + k + 1) & 3]];
                minSide2 = fminf(minSide2, dx * dx + dy * dy);
            }
            resume = bestCost < minSide2 * 0.0625f;
        }
        if (prevNumSeen != 0 || resume) {
            // Continuita': la rotazione ciclica piu' vicina all'ultima posizione mostrata,
            // cosi' l'identita' dei LED non salta anche con l'arma molto ruotata.
            float bestCost = FLT_MAX;
            for (uint8_t r = 0; r < 4; r++) {
                float cost = 0.0f;
                for (uint8_t k = 0; k < 4; k++) {
                    const float dx = sx[ord[(r + k) & 3]] - (float)FinalX[k];
                    const float dy = sy[ord[(r + k) & 3]] - (float)FinalY[k];
                    cost += dx * dx + dy * dy;
                }
                if (cost < bestCost) { bestCost = cost; best = r; }
            }
            // Se la continuita' implica una rotazione dell'arma di oltre 80 gradi superiore a
            // quella della regola di partenza, si torna alla regola di partenza: cosi' un
            // aggancio sbagliato (partenza con l'arma molto ruotata) si corregge tornando dritti.
            if (best != cold) {
                // rotazione apparente: diagonale 0->2 ruotata di 90 gradi. Con l'arma molto di lato la
                // diagonale orizzontale si inclina per la prospettiva, quella verticale molto meno
                // (solo con la mira laterale e l'arma insieme molto sopra o sotto lo schermo).
                auto rollOf = [&](uint8_t r) {
                    const float vx = sx[ord[(r + 2) & 3]] - sx[ord[r]], vy = sy[ord[(r + 2) & 3]] - sy[ord[r]];
                    return fabsf(atan2f(-vx, vy));
                };
                if (rollOf(best) >= rollOf(cold) + 80.0f * 0.0174533f) best = cold;
            }
            // Se le diagonali sono chiare e la continuita' scambia l'orizzontale con la verticale
            // (rotazione dispari rispetto alla regola di partenza), l'aggancio e' sbagliato di 90
            // gradi qualunque sia l'inclinazione dell'arma: si torna alla regola di partenza.
            if (coldShape && ((best ^ cold) & 1)) best = cold;
        }
        for (uint8_t k = 0; k < 4; k++) { Gx[k] = sx[ord[(best + k) & 3]]; Gy[k] = sy[ord[(best + k) & 3]]; }
        if (!ConvexCycle(Gx, Gy)) return;
        mask = 0x0F;
    } else {
        // =========================================================================//
        // ======== 1-3 LED: IDENTIFICAZIONE E RICOSTRUZIONE =======================//
        // =========================================================================//
        // Ogni punto visto riceve un indice diverso, scegliendo l'assegnazione piu'
        // vicina all'ultima posizione mostrata. Con 3 LED:
        // - le assegnazioni specchiate (impossibili per una camera) sono scartate;
        // - alla ripresa dopo il buio si confronta solo la forma (posizioni relative al
        //   baricentro), perche' nel frattempo la mira puo' essersi spostata molto.
        uint8_t bl[3] = {0, 1, 2};
        float bestCost = FLT_MAX;
        const bool centered = (n == 3 && prevNumSeen == 0);
        float scx = 0.0f, scy = 0.0f;
        if (centered) { scx = (sx[0] + sx[1] + sx[2]) / 3.0f; scy = (sy[0] + sy[1] + sy[2]) / 3.0f; }
        for (uint8_t a = 0; a < 4; a++) {
            for (uint8_t b = 0; b < 4; b++) {
                if (n >= 2 && b == a) continue;
                for (uint8_t c = 0; c < 4; c++) {
                    if (n >= 3 && (c == a || c == b)) continue;
                    const uint8_t lab[3] = {a, b, c};
                    if (n == 3) {
                        // I 3 punti presi nell'ordine degli indici devono girare nello stesso verso
                        // del rombo (lo stesso dell'ordinamento ciclico usato con 4 LED).
                        uint8_t o[3] = {0, 1, 2};
                        for (uint8_t i = 1; i < 3; i++)
                            for (uint8_t j = i; j > 0 && lab[o[j]] < lab[o[j - 1]]; j--) { const uint8_t t = o[j]; o[j] = o[j - 1]; o[j - 1] = t; }
                        const float turn = (sx[o[1]] - sx[o[0]]) * (sy[o[2]] - sy[o[1]]) - (sy[o[1]] - sy[o[0]]) * (sx[o[2]] - sx[o[1]]);
                        if (turn <= 0.0f) continue;
                    }
                    float fcx = 0.0f, fcy = 0.0f;
                    if (centered) {
                        fcx = ((float)FinalX[a] + (float)FinalX[b] + (float)FinalX[c]) / 3.0f;
                        fcy = ((float)FinalY[a] + (float)FinalY[b] + (float)FinalY[c]) / 3.0f;
                    }
                    float cost = 0.0f;
                    for (uint8_t i = 0; i < n; i++) {
                        const float dx = (sx[i] - scx) - ((float)FinalX[lab[i]] - fcx);
                        const float dy = (sy[i] - scy) - ((float)FinalY[lab[i]] - fcy);
                        cost += dx * dx + dy * dy;
                    }
                    if (cost < bestCost) { bestCost = cost; bl[0] = a; bl[1] = b; bl[2] = c; }
                    if (n < 3) break;
                }
                if (n < 2) break;
            }
        }
        if (bestCost == FLT_MAX) return;  // 3 punti allineati: nessuna assegnazione valida
        for (uint8_t k = 0; k < 4; k++) { Gx[k] = prevGeomX[k]; Gy[k] = prevGeomY[k]; }
        for (uint8_t i = 0; i < n; i++) { Gx[bl[i]] = sx[i]; Gy[bl[i]] = sy[i]; mask |= 1u << bl[i]; }

        if (n == 3) {
            // Trasformazione affine esatta: il punto mancante mantiene, rispetto ai 3 visti,
            // le stesse coordinate relative che aveva nella geometria precedente.
            uint8_t m = 0;
            while (mask & (1u << m)) m++;
            const uint8_t k0 = (m + 1) & 3, k1 = (m + 2) & 3, k2 = (m + 3) & 3;
            const float ux = prevGeomX[k1] - prevGeomX[k0], uy = prevGeomY[k1] - prevGeomY[k0];
            const float vx = prevGeomX[k2] - prevGeomX[k0], vy = prevGeomY[k2] - prevGeomY[k0];
            const float det = ux * vy - uy * vx;
            if (fabsf(det) > 1.0f) {
                const float qx = prevGeomX[m] - prevGeomX[k0], qy = prevGeomY[m] - prevGeomY[k0];
                const float cu = (qx * vy - qy * vx) / det, cv = (ux * qy - uy * qx) / det;
                Gx[m] = Gx[k0] + cu * (Gx[k1] - Gx[k0]) + cv * (Gx[k2] - Gx[k0]);
                Gy[m] = Gy[k0] + cu * (Gy[k1] - Gy[k0]) + cv * (Gy[k2] - Gy[k0]);
            } else {
                // Geometria precedente degenere: parallelogramma
                Gx[m] = Gx[k0] + Gx[k2] - Gx[k1];
                Gy[m] = Gy[k0] + Gy[k2] - Gy[k1];
            }
        } else if (n == 2) {
            // Similitudine (rotazione + scala uniforme + traslazione) che porta i due punti
            // della geometria precedente sui due punti visti, applicata agli altri due.
            const uint8_t i1 = bl[0], i2 = bl[1];
            const float vx = prevGeomX[i2] - prevGeomX[i1], vy = prevGeomY[i2] - prevGeomY[i1];
            const float wx = Gx[i2] - Gx[i1], wy = Gy[i2] - Gy[i1];
            const float vv = vx * vx + vy * vy;
            if (vv < 64.0f) return;
            const float sa = (wx * vx + wy * vy) / vv, sb = (vx * wy - vy * wx) / vv;
            const float s2 = sa * sa + sb * sb;
            if (s2 < 0.25f || s2 > 4.0f) return;  // scala non plausibile in un frame
            for (uint8_t k = 0; k < 4; k++) {
                if (k == i1 || k == i2) continue;
                const float qx = prevGeomX[k] - prevGeomX[i1], qy = prevGeomY[k] - prevGeomY[i1];
                Gx[k] = Gx[i1] + sa * qx - sb * qy;
                Gy[k] = Gy[i1] + sb * qx + sa * qy;
            }
        } else {
            // Un solo LED: traslazione della geometria precedente
            const uint8_t i1 = bl[0];
            const float dx = Gx[i1] - prevGeomX[i1], dy = Gy[i1] - prevGeomY[i1];
            for (uint8_t k = 0; k < 4; k++) {
                if (k == i1) continue;
                Gx[k] = prevGeomX[k] + dx;
                Gy[k] = prevGeomY[k] + dy;
            }
        }
        if (!ConvexCycle(Gx, Gy)) return;
    }

    // =========================================================================//
    // ======== MOLLA: TRANSIZIONE MORBIDA AL CAMBIO DEI LED MISURATI ==========//
    // =========================================================================//
    if (start && prevNumSeen > 0) {
        // Movimento dei LED misurati sia ora sia nel frame precedente
        const uint8_t common = mask & prevSeenMask;
        float mx = 0.0f, my = 0.0f;
        uint8_t cnt = 0;
        for (uint8_t k = 0; k < 4; k++) {
            if (common & (1u << k)) { mx += Gx[k] - prevGeomX[k]; my += Gy[k] - prevGeomY[k]; cnt++; }
        }
        if (cnt) { mx /= cnt; my /= cnt; }

        // Cambia l'insieme dei LED misurati: si riparte dalla posizione mostrata
        // (piu' il movimento attuale), e la differenza diventa il debito della molla.
        if (mask != prevSeenMask) {
            spring_fast = n > prevNumSeen;
            for (uint8_t k = 0; k < 4; k++) {
                offset_X[k] = (float)FinalX[k] + mx - Gx[k];
                offset_Y[k] = (float)FinalY[k] + my - Gy[k];
            }
        }

        float largest_offset = 0.0f;
        for (uint8_t k = 0; k < 4; k++) {
            largest_offset = fmaxf(largest_offset, fabsf(offset_X[k]));
            largest_offset = fmaxf(largest_offset, fabsf(offset_Y[k]));
        }

        if (largest_offset > 0.0f) {
            // Movimento filtrato di ogni LED (media del segno prima del modulo: riduce il
            // jitter senza annullare rotazione o zoom). Si azzera quando cambia il numero
            // di LED visti, cosi' una velocita' vecchia non accelera un nuovo ritorno lento.
            float motion_sum = 0.0f;
            for (uint8_t k = 0; k < 4; k++) {
                if (prevNumSeen != n) motion_X[k] = motion_Y[k] = 0.0f;
                if (common & (1u << k)) {
                    motion_X[k] = (motion_X[k] + Gx[k] - prevGeomX[k]) * 0.5f;
                    motion_Y[k] = (motion_Y[k] + Gy[k] - prevGeomY[k]) * 0.5f;
                    motion_sum += fabsf(motion_X[k]) + fabsf(motion_Y[k]);
                }
            }
            const float movement = cnt ? motion_sum / (float)cnt : 0.0f;

            // Stessa legge di Square_Advanced: con 4 LED il ritorno e' lento da fermi e
            // veloce in movimento (dove lo scostamento non si nota); con 1-3 LED e' lineare
            // (lento per gli scostamenti grandi), salvo dopo un aumento dei LED misurati:
            // allora e' proporzionale alla velocita' minima.
            const float MIN_RETURN = 0.03125f;
            const float MAX_RETURN = 0.5f;
            const float MOTION_GAIN = 32.0f;
            const float hx = Gx[2] - Gx[0], hy = Gy[2] - Gy[0];
            const float scale = fmaxf(sqrtf(hx * hx + hy * hy), 1.0f) * FPS_NORMALIZATION;
            const float motion = MOTION_GAIN * movement;

            float return_rate;
            if (n != 4 || motion <= scale * MIN_RETURN)
                return_rate = MIN_RETURN;
            else if (motion >= scale * MAX_RETURN)
                return_rate = MAX_RETURN;
            else
                return_rate = motion / scale;

            const float consumo = fmaxf(0.5f * FPS_NORMALIZATION,
                fmaxf(movement * COSTANTE_MOLLA, (n == 4 || spring_fast) ? largest_offset * return_rate : 0.0f));

            // Tutti gli offset si riducono dello stesso fattore: la forma della correzione
            // resta la stessa e non ci sono scatti.
            const float factor = largest_offset > consumo ? 1.0f - consumo / largest_offset : 0.0f;
            for (uint8_t k = 0; k < 4; k++) {
                offset_X[k] *= factor;
                offset_Y[k] *= factor;
            }
        }
    } else {
        for (uint8_t k = 0; k < 4; k++) offset_X[k] = offset_Y[k] = 0.0f;
    }

    // =========================================================================//
    // ======== APPLICAZIONE FINALE E SALVATAGGIO STATO ========================//
    // =========================================================================//
    for (uint8_t k = 0; k < 4; k++) {
        FinalX[k] = (int)lroundf(Gx[k] + offset_X[k]);
        FinalY[k] = (int)lroundf(Gy[k] + offset_Y[k]);
        prevGeomX[k] = Gx[k];
        prevGeomY[k] = Gy[k];
        see[k] = (see[k] << 1) | ((mask >> k) & 1);
    }
    prevSeenMask = mask;
    prevNumSeen = n;
    darkFrames = 0;
    if (n == 4) start = 0xFF;

    medianX = (FinalX[0] + FinalX[1] + FinalX[2] + FinalX[3]) / 4;
    medianY = (FinalY[0] + FinalY[1] + FinalY[2] + FinalY[3]) / 4;
    height2 = (float)(FinalY[2] - FinalY[0]);
    width2 = (float)(FinalX[1] - FinalX[3]);
    // Rotazione della diagonale orizzontale, come nel codice originale (usata in calibrazione)
    angle = atan2f((float)(FinalY[1] - FinalY[3]), (float)(FinalX[1] - FinalX[3]));
}

#endif // USE_DIAMOND_ADVANCED
