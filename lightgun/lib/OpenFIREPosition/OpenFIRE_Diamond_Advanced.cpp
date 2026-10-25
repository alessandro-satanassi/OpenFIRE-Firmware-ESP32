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
#include "OpenFIRE_Anchor.h"
#include <cfloat>

void OpenFIRE_Diamond::configure(const CameraProfile& profile)
{
    camToMouseShift = profile.camToMouseShift;
    camToMouseMult = profile.camToMouseMult;
    FPS_NORMALIZATION = profile.fps ? 209.0f / (float)profile.fps : 1.0f;
    shortDarkFrames = profile.fps / 26;   // buio breve: circa 40 ms (8 fotogrammi a 209 fps)
    {
        // Distanza minima tra 2 LED: un pixel reale del sensore, in unita' mouse
        const float sensorPixel = profile.sensorResX ? (float)profile.mouseResX / (float)profile.sensorResX : 1.0f;
        minLedDist2 = sensorPixel * sensorPixel;
    }
    // Focale nominale e centro dell'immagine in unita' mouse (ricostruzione dei LED nascosti)
    focal = profile.mouseResX * 0.5f * profile.lensFocalRatio;
    imageCX = profile.mouseResX * 0.5f;
    imageCY = profile.mouseResY * 0.5f;
    anchorReconstruction.configure(focal, imageCX, imageCY);

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
    anchor_ok = false;
    for (uint8_t i = 0; i < 4; i++) {
        prevGeomX[i] = (float)FinalX[i];
        prevGeomY[i] = (float)FinalY[i];
        anchorX[i] = (float)FinalX[i];
        anchorY[i] = (float)FinalY[i];
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
    if (n == 0) { prevNumSeen = 0; anchor_ok = false; if (darkFrames < 255) darkFrames++; return; }
    // Dopo il buio servono almeno 3 LED: con 1-2 l'identificazione non e' affidabile
    if (prevNumSeen == 0 && n < 3) { anchor_ok = false; if (darkFrames < 255) darkFrames++; return; }

    float Gx[4], Gy[4];
    uint8_t mask = 0;

    if (n == 4) {
        // =========================================================================//
        // ======== 4 LED: ORDINAMENTO CICLICO ATTORNO AL CENTROIDE ================//
        // =========================================================================//
        // L'angolo cresce alto -> X max -> basso -> X min (Y verso il basso), qualunque
        // sia la distanza o la prospettiva. Resta solo da scegliere quale punto e' "alto".
        uint8_t ord[4];
        if (!OrderFourPoints(sx,sy,ord)) return;

        // Preserve a recent, close correspondence. Shape must not turn a valid lock
        // by 90 degrees merely because perspective shortened one diagonal.
        uint8_t best = 0;
        float bestCost = FLT_MAX, minSide2 = FLT_MAX;
        for (uint8_t k = 0; k < 4; ++k) {
            const float dx = sx[ord[k]] - sx[ord[(k+1)&3]], dy = sy[ord[k]] - sy[ord[(k+1)&3]];
            minSide2 = fminf(minSide2, dx*dx + dy*dy);
        }
        if (start) {
            for (uint8_t r = 0; r < 4; ++r) {
                float cost = 0.0f;
                for (uint8_t k = 0; k < 4; ++k) {
                    const float dx = sx[ord[(r+k)&3]] - prevGeomX[k];
                    const float dy = sy[ord[(r+k)&3]] - prevGeomY[k];
                    cost += dx*dx + dy*dy;
                }
                if (cost < bestCost) { bestCost=cost; best=r; }
            }
        }
        const float step = COSTANTE_MOLLA * FPS_NORMALIZATION;
        const bool recent = start && (prevNumSeen == 4 || (prevNumSeen == 0 && darkFrames <= shortDarkFrames));
        if (!(recent && bestCost < minSide2 * step * step)) {
            // The two screen diagonals share a physical midpoint. Their intersection
            // gives relative endpoint depths, hence their lengths in camera space.
            // No square root: only compare squared lengths.
            const float ax=sx[ord[0]], ay=sy[ord[0]], bx=sx[ord[1]], by=sy[ord[1]];
            const float ux=sx[ord[2]]-ax, uy=sy[ord[2]]-ay;
            const float vx=sx[ord[3]]-bx, vy=sy[ord[3]]-by;
            float da=ux*ux+uy*uy, db=vx*vx+vy*vy;
            float verticalA=sy[ord[2]]-ay, verticalB=sy[ord[3]]-by;
            CameraQuad quad;
            if (quad.project(focal,imageCX,imageCY,sx,sy,ord)) {
                da=quad.distance2(0,2);db=quad.distance2(1,3);
                verticalA=quad.y[2]-quad.y[0];verticalB=quad.y[3]-quad.y[1];
            }
            const float ref=diagRatioRef>=1?diagRatioRef:1/diagRatioRef;
            const float threshold=1+(ref-1)*COSTANTE_MOLLA;
            if (da>threshold*threshold*db || db>threshold*threshold*da) {
                const bool vertA = diagRatioRef>=1 ? da<db : da>db;
                best=vertA?0:1;
                if ((vertA?verticalA:verticalB)<0) best+=2;
            } else {
                best=0;for(uint8_t k=1;k<4;k++)if(sy[ord[k]]<sy[ord[best]])best=k;
            }
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
            const float* refX = anchor_ok ? anchorX : prevGeomX;
            const float* refY = anchor_ok ? anchorY : prevGeomY;
            uint8_t m = 0;
            while (mask & (1u << m)) ++m;
            const uint8_t k0 = (m+1)&3, k1 = (m+2)&3, k2 = (m+3)&3;
            const float ax[3] = {refX[k0],refX[k1],refX[k2]}, ay[3] = {refY[k0],refY[k1],refY[k2]};
            const float qx[3] = {Gx[k0],Gx[k1],Gx[k2]}, qy[3] = {Gy[k0],Gy[k1],Gy[k2]};
            float hx = refX[m], hy = refY[m];
            if (anchorReconstruction.reconstruct(3,ax,ay,qx,qy,1,&hx,&hy)) {
                Gx[m] = hx; Gy[m] = hy;
            } else {
                // Affine fallback: calculate it only when the virtual view is unsuitable.
                const float ux = refX[k1]-refX[k0], uy = refY[k1]-refY[k0];
                const float vx = refX[k2]-refX[k0], vy = refY[k2]-refY[k0];
                const float det = ux*vy-uy*vx;
                if (fabsf(det)>1.0f) {
                    const float dx = refX[m]-refX[k0], dy = refY[m]-refY[k0];
                    const float cu = (dx*vy-dy*vx)/det, cv = (ux*dy-uy*dx)/det;
                    Gx[m] = Gx[k0]+cu*(Gx[k1]-Gx[k0])+cv*(Gx[k2]-Gx[k0]);
                    Gy[m] = Gy[k0]+cu*(Gy[k1]-Gy[k0])+cv*(Gy[k2]-Gy[k0]);
                } else {
                    Gx[m] = Gx[k0]+Gx[k2]-Gx[k1]; Gy[m] = Gy[k0]+Gy[k2]-Gy[k1];
                }
            }
        } else if (n == 2) {
            const float* refX = anchor_ok ? anchorX : prevGeomX;
            const float* refY = anchor_ok ? anchorY : prevGeomY;
            const uint8_t i1 = bl[0], i2 = bl[1];
            const float vx = refX[i2]-refX[i1], vy = refY[i2]-refY[i1];
            const float wx = Gx[i2]-Gx[i1], wy = Gy[i2]-Gy[i1];
            const float vv = vx*vx+vy*vy;
            if (vv<minLedDist2) return;
            // The reference preserves shape; plausibility compares consecutive frames.
            const float pvx = prevGeomX[i2]-prevGeomX[i1], pvy = prevGeomY[i2]-prevGeomY[i1];
            const float pv2 = pvx*pvx+pvy*pvy, w2 = wx*wx+wy*wy;
            if (pv2<minLedDist2 || w2<pv2*COSTANTE_MOLLA || w2*COSTANTE_MOLLA>pv2) return;
            uint8_t hidden[2], nh = 0;
            for (uint8_t k=0;k<4;++k) if (k!=i1 && k!=i2) hidden[nh++]=k;
            const float ax[2] = {refX[i1],refX[i2]}, ay[2] = {refY[i1],refY[i2]};
            const float qx[2] = {Gx[i1],Gx[i2]}, qy[2] = {Gy[i1],Gy[i2]};
            float hx[2] = {refX[hidden[0]],refX[hidden[1]]}, hy[2] = {refY[hidden[0]],refY[hidden[1]]};
            if (anchorReconstruction.reconstruct(2,ax,ay,qx,qy,2,hx,hy)) {
                for (uint8_t j=0;j<2;++j) { Gx[hidden[j]]=hx[j]; Gy[hidden[j]]=hy[j]; }
            } else {
                // Similarity fallback. Avoid its two divisions on the normal virtual path.
                const float sa = (wx*vx+wy*vy)/vv, sb = (vx*wy-vy*wx)/vv;
                for (uint8_t j=0;j<2;++j) {
                    const uint8_t k=hidden[j];const float dx=refX[k]-refX[i1],dy=refY[k]-refY[i1];
                    Gx[k] = Gx[i1]+sa*dx-sb*dy; Gy[k] = Gy[i1]+sb*dx+sa*dy;
                }
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
        // Limite numerico, non il campo visivo: un LED ricostruito puo' stare fuori schermo, ma un
        // valore enorme viene solo da una geometria precedente quasi degenere.
        for (uint8_t k = 0; k < 4; k++)
            if (!(fabsf(Gx[k]) <= 16383.0f && fabsf(Gy[k]) <= 16383.0f)) return;
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
    if (n == 4) {
        start = 0xFF;
        for (uint8_t k = 0; k < 4; k++) { anchorX[k] = Gx[k]; anchorY[k] = Gy[k]; }
        anchor_ok = true;
    }

    medianX = (FinalX[0] + FinalX[1] + FinalX[2] + FinalX[3]) / 4;
    medianY = (FinalY[0] + FinalY[1] + FinalY[2] + FinalY[3]) / 4;
    height2 = (float)(FinalY[2] - FinalY[0]);
    width2 = (float)(FinalX[1] - FinalX[3]);
    // Rotazione della diagonale orizzontale, come nel codice originale (usata in calibrazione)
    angle = atan2f((float)(FinalY[1] - FinalY[3]), (float)(FinalX[1] - FinalX[3]));
}

#endif // USE_DIAMOND_ADVANCED
