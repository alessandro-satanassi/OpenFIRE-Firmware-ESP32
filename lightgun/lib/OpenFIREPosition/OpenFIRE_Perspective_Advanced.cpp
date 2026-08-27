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

// ========================================================================================
// MOTORE ALGEBRICO E REGISTRI CPU (Zero Overhead)
// ========================================================================================

// Calcola la matrice di trasformazione prospettica proiettando un quadrato ideale
// (la calibrazione TV) sul quadrilatero irregolare letto dalla telecamera IR.
inline bool computeSquareToQuad(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float dx1 = x1 - x2; float dy1 = y1 - y2;
  float dx2 = x3 - x2; float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3; float sy = y0 - y1 + y2 - y3;

  float det = (dx1 * dy2 - dx2 * dy1);
  float g = 0.0f, h = 0.0f;
  
  // Limite di degenerazione matematica (1e-8 è adatto alla scala normalizzata).
  // Un quadrilatero realmente degenere non deve produrre una matrice parziale.
  if (fabsf(det) <= 1e-8f && (sx != 0.0f || sy != 0.0f)) return false;

  if (fabsf(det) > 1e-8f) {
    float invDet = 1.0f / det;
    g = (sx * dy2 - dx2 * sy) * invDet;
    h = (dx1 * sy - sx * dy1) * invDet;
  }

  mat[0] = x1 - x0 + g * x1; mat[1] = y1 - y0 + g * y1; mat[2] = g;
  mat[3] = x3 - x0 + h * x3; mat[4] = y3 - y0 + h * y3; mat[5] = h;
  mat[6] = x0;               mat[7] = y0;               mat[8] = 1.0f;
  return true;
}

// Inversa della funzione precedente. Genera la matrice dei co-fattori.
// Indispensabile per mappare il singolo pixel letto sul sensore (il "proiettile")
// sulla superficie fisica e non distorta dello schermo della TV.
inline bool computeQuadToSquare(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float dx1 = x1 - x2; float dy1 = y1 - y2;
  float dx2 = x3 - x2; float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3; float sy = y0 - y1 + y2 - y3;

  float det = (dx1 * dy2 - dx2 * dy1);
  float g = 0.0f, h = 0.0f;

  // Se il termine prospettico è singolare ma il quadrilatero non è affine,
  // non è possibile ricavare una trasformazione valida.
  if (fabsf(det) <= 1e-8f && (sx != 0.0f || sy != 0.0f)) return false;
  
  if (fabsf(det) > 1e-8f) {
    float invDet = 1.0f / det;
    g = (sx * dy2 - dx2 * sy) * invDet;
    h = (dx1 * sy - sx * dy1) * invDet;
  }

  float a = x1 - x0 + g * x1; float b = x3 - x0 + h * x3; float c = x0;
  float d = y1 - y0 + g * y1; float e = y3 - y0 + h * y3; float f = y0;

  float A = e - f * h; float B = c * h - b; float C = b * f - c * e;
  float D = f * g - d; float E = a - c * g; float F = c * d - a * f;
  float G = d * h - e * g; float H = b * g - a * h; float I = a * e - b * d;
  
  float invDetQuad = (a * A + b * D + c * G);
  
  // Se la camera restituisce punti anomali/allineati (determinante nullo),
  // segnaliamo esplicitamente il fallimento. warp() manterrà l'ultimo output valido.
  if (fabsf(invDetQuad) <= 1e-8f) return false;

  float idet = 1.0f / invDetQuad;
  mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
  mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
  mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
  return true;
}

// Unrolling manuale del calcolo della matrice. Non usare cicli FOR qui:
// l'ESP32 esegue questa istruzione flat in una manciata di cicli di clock.
inline void multMats(const float* a, const float* b, float* res) {
  res[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
  res[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
  res[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
  res[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
  res[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  res[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
  res[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
  res[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
  res[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

// ========================================================================================
// FISICA E CALCOLO AREA
// ========================================================================================

// Algoritmo di de-warping spaziale di Brown-Conrady semplificato.
// Necessario perché ai bordi estremi del campo visivo della telecamera IR
// i LED "spanciano" a causa della lente sferica, falsando il calcolo prospettico.
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
}

inline void OpenFIRE_Perspective::applyLensCorrection(float &x, float &y) {
  if (k1 == 0.0f && k2 == 0.0f) return;
  
  float dx = x - cx; 
  float dy = y - cy;
  
  float nx = dx * invCx;
  float ny = dy * invCy;
  
  float r2 = (nx * nx) + (ny * ny);
  float distortion = 1.0f + (k1 * r2);

  // Secondo termine radiale: con k2 = 0.0f il percorso DFRobot resta
  // matematicamente identico alla versione precedente.
  if (k2 != 0.0f) distortion += k2 * r2 * r2;
  
  // Limiti di sicurezza intrinseci per evitare deformazioni distruttive del quadrilatero
  // nel caso in cui i parametri di distorsione sfuggano al controllo dell'utente.
  if (distortion > 1.2f) distortion = 1.2f;
  if (distortion < 0.8f) distortion = 0.8f;
  
  x = cx + (dx * distortion); 
  y = cy + (dy * distortion);
}

// Calcolo dell'area del poligono tramite prodotto vettoriale incrociato (Shoelace Formula).
// Viene usata come indicatore Z-Depth per l'effetto parallasse.
inline float OpenFIRE_Perspective::calculateQuadArea(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  return 0.5f * fabsf((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0));
}

// ========================================================================================
// ESECUZIONE WARP (CUORE DEL SISTEMA)
// ========================================================================================

void OpenFIRE_Perspective::warp(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  // Riorganizzazione Ciclica Pura: TL, TR, BR, BL.
  // Essenziale. L'algoritmo matematico esplode se i vertici del quadrilatero 
  // si incrociano ad X (Bug Clessidra). Questa mappatura fissa garantisce 
  // un perimetro sempre chiuso e ordinato per il proiettore matriciale.
  float pt_TL_x = x0,  pt_TL_y = y0;
  float pt_TR_x = x1,  pt_TR_y = y1;
  float pt_BR_x = x3,  pt_BR_y = y3; 
  float pt_BL_x = x2,  pt_BL_y = y2; 

  // Prima della correzione ottica verifichiamo che i quattro punti descrivano
  // ancora una superficie 2D significativa. Una linea o un punto non possono
  // definire un'omografia; la correzione radiale non deve trasformarli
  // artificialmente in un quadrilatero apparentemente valido.
  float rawArea = calculateQuadArea(pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);
  if (rawArea <= minQuadArea) return;

  applyLensCorrection(pt_TL_x, pt_TL_y);
  applyLensCorrection(pt_TR_x, pt_TR_y);
  applyLensCorrection(pt_BR_x, pt_BR_y);
  applyLensCorrection(pt_BL_x, pt_BL_y);

  float currentArea = calculateQuadArea(pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);

  // Snapshot della calibrazione iniziale.
  // Fissa il "punto zero" della distanza e della rotazione. Tutto il gioco
  // viene calcolato come differenza proporzionale da questo esatto fotogramma.
  if (!init) {
    if (currentArea > minQuadArea) {
      // Inizializzazione della base dello schermo con lo stesso ordine ciclico (TL, TR, BR, BL)
      if (!computeSquareToQuad(dstmatrix, dx0 * INV_NORM_SCALE, dy0 * INV_NORM_SCALE,
                                      dx1 * INV_NORM_SCALE, dy1 * INV_NORM_SCALE,
                                      dx3 * INV_NORM_SCALE, dy3 * INV_NORM_SCALE,
                                      dx2 * INV_NORM_SCALE, dy2 * INV_NORM_SCALE)) {
        return;
      }
      baseArea = currentArea;
      smoothedArea = currentArea;
      init = true;
    }
    if (!init) return; 
  }

  // Freeze di sicurezza: se la telecamera viene coperta accidentalmente (area vicina a 0), 
  // impediamo al calcolo di assorbire l'errore, salvando la stabilità della z-depth.
  if (currentArea > minQuadArea) {
    smoothedArea = (0.1f * currentArea) + (0.9f * smoothedArea);
  }

  float dynamicSrcY = srcY;
  
  // Applicazione Parallasse: se la distanza del giocatore diminuisce (area sale rispetto alla base),
  // il centro ottico (dynamicSrcY) viene spinto verso l'alto/basso per compensare 
  // fisicamente l'altezza della canna rispetto al sensore montato.
  if (parallaxFactor != 0.0f && baseArea > minQuadArea && smoothedArea > minQuadArea) {
    float distanceRatio = sqrtf(baseArea / smoothedArea);
    dynamicSrcY += parallaxFactor * (distanceRatio - 1.0f);
  }

  pt_TL_x *= INV_NORM_SCALE; pt_TL_y *= INV_NORM_SCALE;
  pt_TR_x *= INV_NORM_SCALE; pt_TR_y *= INV_NORM_SCALE;
  pt_BR_x *= INV_NORM_SCALE; pt_BR_y *= INV_NORM_SCALE;
  pt_BL_x *= INV_NORM_SCALE; pt_BL_y *= INV_NORM_SCALE;

  if (!computeQuadToSquare(srcmatrix, pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y)) {
    return;
  }
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  float normSrcX = srcX * INV_NORM_SCALE;
  float normSrcY = dynamicSrcY * INV_NORM_SCALE;

  // Moltiplicazione del vettore proiettile [X, Y, 1] contro la matrice di trasformazione fusa.
  float r0 = (normSrcX * warpmatrix[0] + normSrcY * warpmatrix[3] + warpmatrix[6]);
  float r1 = (normSrcX * warpmatrix[1] + normSrcY * warpmatrix[4] + warpmatrix[7]);
  float r3 = (normSrcX * warpmatrix[2] + normSrcY * warpmatrix[5] + warpmatrix[8]);
  
  if (fabsf(r3) > 1e-8f) {
    float invR3 = 1.0f / r3;
    // Conversione divisione prospettica (w) in spazio bidimensionale (x,y)
    float dstX_float = roundf((r0 * invR3) * NORM_SCALE);
    float dstY_float = roundf((r1 * invR3) * NORM_SCALE);

    // Clamp fisico Assoluto: blocca le coordinate finali prima del casting ad intero.
    // Se spariamo parallelamente allo schermo, la proiezione prospettica tende all'infinito.
    // Questo clamp previene l'overflow int32 che causerebbe il crash UB (Undefined Behavior) del core ESP32.
    if (dstX_float > 2000000000.0f) dstX_float = 2000000000.0f;
    if (dstX_float < -2000000000.0f) dstX_float = -2000000000.0f;
    if (dstY_float > 2000000000.0f) dstY_float = 2000000000.0f;
    if (dstY_float < -2000000000.0f) dstY_float = -2000000000.0f;

    dstX = (int)dstX_float;
    dstY = (int)dstY_float;
  }
}

void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setLensCorrection(float coefficientK1, float coefficientK2) { k1 = coefficientK1; k2 = coefficientK2; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = adjustedX; srcY = adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED
