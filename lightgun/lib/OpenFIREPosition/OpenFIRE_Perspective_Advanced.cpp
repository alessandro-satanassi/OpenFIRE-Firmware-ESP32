#ifdef USE_PERSPECTIVE_ADVANCED
/*!
 * @file OpenFIRE_Perspective_Advanced.cpp
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
 * Derived from Wiimote Whiteboard library:
 * Copyright 2021 88hcsif
 * Copyright (c) 2008 Stephane Duchesneau
 * by Stephane Duchesneau <stephane.duchesneau@gmail.com>
 * Ported from Johnny Lee's C# WiiWhiteboard project (Warper.cs file)
 */

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>
// EN: Rescale the former 10000-based determinant tolerance to the binary normalization.
// IT: Riportiamo la tolleranza del determinante dalla precedente scala 10000 a quella binaria.
static constexpr float NORMALIZED_DET_EPS = 1e-8f * (10000.0f / NORM_SCALE) * (10000.0f / NORM_SCALE);

// EN: Cached unit-square to destination transform; vertices follow the perimeter.
// IT: Trasformazione memorizzata dal quadrato unitario alla destinazione, con vertici ciclici.
inline bool computeSquareToQuad(float* mat, float x0, float y0, float x1, float y1, float x2, float y2,
                                float x3, float y3) {
    float dx1 = x1 - x2;
    float dy1 = y1 - y2;
    float dx2 = x3 - x2;
    float dy2 = y3 - y2;
    float sx = x0 - x1 + x2 - x3;
    float sy = y0 - y1 + y2 - y3;

    float det = (dx1 * dy2 - dx2 * dy1);
    float g = 0.0f, h = 0.0f;

    if (!isfinite(det) || det == 0.0f || (fabsf(det) <= NORMALIZED_DET_EPS && (sx != 0.0f || sy != 0.0f)))
        return false;

    if (fabsf(det) > NORMALIZED_DET_EPS) {
        float invDet = 1.0f / det;
        g = (sx * dy2 - dx2 * sy) * invDet;
        h = (dx1 * sy - sx * dy1) * invDet;
    }

    mat[0] = x1 - x0 + g * x1;
    mat[1] = y1 - y0 + g * y1;
    mat[2] = g;
    mat[3] = x3 - x0 + h * x3;
    mat[4] = y3 - y0 + h * y3;
    mat[5] = h;
    mat[6] = x0;
    mat[7] = y0;
    mat[8] = 1.0f;

    for (unsigned i = 0; i < 9; ++i) {
        if (!isfinite(mat[i]))
            return false;
    }
    return true;
}

// EN: Evaluate the inverse at (0,0). Only the three homogeneous coordinates are needed.
// IT: Valutiamo l'inversa in (0,0): servono soltanto le tre coordinate omogenee.
// EN: Their common determinant scale cancels in the final perspective division.
// IT: La scala comune del determinante si elimina nella divisione prospettica finale.
inline bool computeOriginToSquare(float* query, float x0, float y0, float x1, float y1, float x2, float y2,
                                  float x3, float y3, float& sourceDet) {
    float dx1 = x1 - x2, dy1 = y1 - y2, dx2 = x3 - x2, dy2 = y3 - y2;
    float sx = x0 - x1 + x2 - x3, sy = y0 - y1 + y2 - y3;
    float det = dx1 * dy2 - dx2 * dy1, g = 0, h = 0;
    if (fabsf(det) <= NORMALIZED_DET_EPS && (sx != 0.0f || sy != 0.0f))
        return false;
    if (fabsf(det) > NORMALIZED_DET_EPS) {
        float inv = 1.0f / det;
        g = (sx * dy2 - dx2 * sy) * inv;
        h = (dx1 * sy - sx * dy1) * inv;
    }
    float a = x1 - x0 + g * x1, b = x3 - x0 + h * x3, c = x0;
    float d = y1 - y0 + g * y1, e = y3 - y0 + h * y3, f = y0;
    float C = b * f - c * e, F = c * d - a * f, I = a * e - b * d;
    sourceDet = I + g * C + h * F;
    if (fabsf(sourceDet) <= NORMALIZED_DET_EPS)
        return false;
    query[0] = C;
    query[1] = F;
    query[2] = I;
    return true;
}

// EN/IT: Camera geometry and lens coefficients / geometria camera e coefficienti della lente.
void OpenFIRE_Perspective::configure(const CameraProfile& profile) {
    cx = (float)profile.mouseResX * 0.5f;
    cy = (float)profile.mouseResY * 0.5f;
    invCx = 1.0f / cx;
    invCy = 1.0f / cy;
    minQuadArea = (float)profile.mouseResX * (float)profile.mouseResY * MIN_QUAD_AREA_RATIO;
    srcX = cx;
    srcY = cy;
    k1 = profile.lensRadialK1;
    k2 = profile.lensRadialK2;
    lensCorrectionMax = profile.lensCorrectionMax;
    lensCorrectionMin = profile.lensCorrectionMin;
}

// EN: Use the same profile-defined radial map for the LEDs and the aim point.
// IT: Usiamo la stessa correzione radiale del profilo per LED e punto di mira.
inline void OpenFIRE_Perspective::applyLensCorrection(float& x, float& y) {
    if (k1 == 0.0f && k2 == 0.0f)
        return;

    float dx = x - cx;
    float dy = y - cy;

    float nx = dx * invCx;
    float ny = dy * invCy;

    float r2 = (nx * nx) + (ny * ny);
    float distortion = 1.0f + (k1 * r2);

    if (k2 != 0.0f)
        distortion += k2 * r2 * r2;

    if (distortion > lensCorrectionMax)
        distortion = lensCorrectionMax;
    if (distortion < lensCorrectionMin)
        distortion = lensCorrectionMin;

    x = cx + (dx * distortion);
    y = cy + (dy * distortion);
}

// EN/IT: Area of a cyclic quadrilateral / area di un quadrilatero con vertici ciclici.
inline float OpenFIRE_Perspective::calculateQuadArea(float x0, float y0, float x1, float y1, float x2,
                                                     float y2, float x3, float y3) {
    return 0.5f * fabsf((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0));
}

void OpenFIRE_Perspective::warp(float x0, float y0, float x1, float y1, float x2, float y2, float x3,
                                float y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2,
                                float dx3, float dy3) {

    // EN: Convert caller order 0,1,2,3 to perimeter order 0,1,3,2.
    // IT: Convertiamo l'ordine del chiamante 0,1,2,3 nell'ordine ciclico 0,1,3,2.
    // EN/IT: Diamond callers adapt their vertex order / i chiamanti Diamond adattano i propri vertici.
    float pt_TL_x = x0, pt_TL_y = y0;
    float pt_TR_x = x1, pt_TR_y = y1;
    float pt_BR_x = x3, pt_BR_y = y3;
    float pt_BL_x = x2, pt_BL_y = y2;

    // EN: Reject a collapsed raw shape before lens correction can make it appear two-dimensional.
    // IT: Scartiamo una forma raw collassata prima che la lente la faccia apparire bidimensionale.
    float rawArea = calculateQuadArea(pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);
    if (rawArea <= minQuadArea)
        return;

    applyLensCorrection(pt_TL_x, pt_TL_y);
    applyLensCorrection(pt_TR_x, pt_TR_y);
    applyLensCorrection(pt_BR_x, pt_BR_y);
    applyLensCorrection(pt_BL_x, pt_BL_y);

    float currentArea =
        calculateQuadArea(pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);

    // EN: Cache only the destination transform and the parallax area reference.
    // IT: Memorizziamo soltanto la destinazione e l'area di riferimento per la parallasse.
    // EN: The current camera view is solved again each frame, including oblique views.
    // IT: La vista corrente della camera viene ricalcolata a ogni frame, anche se obliqua.
    if (!init) {
        if (currentArea > minQuadArea) {

            if (!computeSquareToQuad(dstmatrix, dx0 * INV_NORM_SCALE, dy0 * INV_NORM_SCALE,
                                     dx1 * INV_NORM_SCALE, dy1 * INV_NORM_SCALE, dx3 * INV_NORM_SCALE,
                                     dy3 * INV_NORM_SCALE, dx2 * INV_NORM_SCALE, dy2 * INV_NORM_SCALE)) {
                return;
            }
            baseArea = currentArea;
            smoothedArea = currentArea;
            init = true;
        }
        if (!init)
            return;
    }

    // EN/IT: Keep the existing area filter / conserviamo il filtro dell'area esistente.
    if (currentArea > minQuadArea) {
        smoothedArea = (0.1f * currentArea) + (0.9f * smoothedArea);
    }

    float dynamicSrcY = srcY;

    if (parallaxFactor != 0.0f && baseArea > minQuadArea && smoothedArea > minQuadArea) {
        float distanceRatio = sqrtf(baseArea / smoothedArea);
        dynamicSrcY += parallaxFactor * (distanceRatio - 1.0f);
    }

    float correctedSrcX = srcX;

    applyLensCorrection(correctedSrcX, dynamicSrcY);

    // EN: Center on the actual lens-corrected aim, including calibration and parallax.
    // IT: Centriamo sulla mira effettiva corretta per la lente, incluse calibrazione e parallasse.
    // EN: The query is now (0,0); binary normalization avoids additional decimal scaling error.
    // IT: Il punto da trasformare e ora (0,0); la scala binaria evita ulteriori errori decimali.
    pt_TL_x = (pt_TL_x - correctedSrcX) * INV_NORM_SCALE;
    pt_TL_y = (pt_TL_y - dynamicSrcY) * INV_NORM_SCALE;
    pt_TR_x = (pt_TR_x - correctedSrcX) * INV_NORM_SCALE;
    pt_TR_y = (pt_TR_y - dynamicSrcY) * INV_NORM_SCALE;
    pt_BR_x = (pt_BR_x - correctedSrcX) * INV_NORM_SCALE;
    pt_BR_y = (pt_BR_y - dynamicSrcY) * INV_NORM_SCALE;
    pt_BL_x = (pt_BL_x - correctedSrcX) * INV_NORM_SCALE;
    pt_BL_y = (pt_BL_y - dynamicSrcY) * INV_NORM_SCALE;

    float sourceQuery[3], sourceDet;
    if (!computeOriginToSquare(sourceQuery, pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x,
                               pt_BL_y, sourceDet)) {
        return;
    }

    // EN: Apply the destination directly to the homogeneous point, without composing matrices.
    // IT: Applichiamo direttamente la destinazione al punto omogeneo, senza comporre matrici.
    const float u = sourceQuery[0];
    const float v = sourceQuery[1];
    const float w = sourceQuery[2];
    const float r0 = u * dstmatrix[0] + v * dstmatrix[3] + w * dstmatrix[6];
    const float r1 = u * dstmatrix[1] + v * dstmatrix[4] + w * dstmatrix[7];
    const float r3 = u * dstmatrix[2] + v * dstmatrix[5] + w * dstmatrix[8];

    // EN: Preserve the original pole tolerance after removing the inverse determinant division.
    // IT: Conserviamo la tolleranza sul polo dopo aver eliminato la divisione per il determinante.
    if (fabsf(r3) > 1e-8f * fabsf(sourceDet)) {
        float invR3 = 1.0f / r3;

        float dstX_float = roundf((r0 * invR3) * NORM_SCALE);
        float dstY_float = roundf((r1 * invR3) * NORM_SCALE);

        // EN/IT: Keep the last valid output on non-finite results / manteniamo l'ultima uscita valida.
        if (!isfinite(dstX_float) || !isfinite(dstY_float))
            return;

        // EN/IT: Bound before converting to int32 / limitiamo i valori prima della conversione a int32.
        if (dstX_float > 2000000000.0f)
            dstX_float = 2000000000.0f;
        if (dstX_float < -2000000000.0f)
            dstX_float = -2000000000.0f;
        if (dstY_float > 2000000000.0f)
            dstY_float = 2000000000.0f;
        if (dstY_float < -2000000000.0f)
            dstY_float = -2000000000.0f;

        dstX = (int)dstX_float;
        dstY = (int)dstY_float;
    }
}

void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setLensCorrection(float coefficientK1, float coefficientK2) {
    k1 = coefficientK1;
    k2 = coefficientK2;
}
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) {
    srcX = adjustedX;
    srcY = adjustedY;
}
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif
