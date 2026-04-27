
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// ========================================================================================
// 1. MOTORE ALGEBRICO E REGISTRI CPU (Zero Overhead)
// ========================================================================================

inline void computeSquareToQuad(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float dx1 = x1 - x2; float dy1 = y1 - y2;
  float dx2 = x3 - x2; float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3; float sy = y0 - y1 + y2 - y3;

  float det = (dx1 * dy2 - dx2 * dy1);
  float g = 0.0f, h = 0.0f;
  
  // Limite di degenerazione matematica (1e-8 è perfetto per la scala normalizzata a 1.0)
  if (fabsf(det) > 1e-8f) {
    float invDet = 1.0f / det;
    g = (sx * dy2 - dx2 * sy) * invDet;
    h = (dx1 * sy - sx * dy1) * invDet;
  }

  mat[0] = x1 - x0 + g * x1; mat[1] = y1 - y0 + g * y1; mat[2] = g;
  mat[3] = x3 - x0 + h * x3; mat[4] = y3 - y0 + h * y3; mat[5] = h;
  mat[6] = x0;               mat[7] = y0;               mat[8] = 1.0f;
}

inline void computeQuadToSquare(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float dx1 = x1 - x2; float dy1 = y1 - y2;
  float dx2 = x3 - x2; float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3; float sy = y0 - y1 + y2 - y3;

  float det = (dx1 * dy2 - dx2 * dy1);
  float g = 0.0f, h = 0.0f;
  
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
  
  // Se la fotocamera restituisce punti anomali/allineati, questa condizione fallisce 
  // e la matrice NON viene scritta. Il mirino si congela sull'ultimo frame valido.
  if (fabsf(invDetQuad) > 1e-8f) {
    float idet = 1.0f / invDetQuad;
    mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
    mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
    mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
  }
}

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
// 2. FISICA E CALCOLO AREA
// ========================================================================================

inline void OpenFIRE_Perspective::applyLensCorrection(float &x, float &y) {
  if (k1 == 0.0f) return; 
  
  float dx = x - CX; 
  float dy = y - CY;
  
  float nx = dx * INV_CX;
  float ny = dy * INV_CY;
  
  float r2 = (nx * nx) + (ny * ny);
  float distortion = 1.0f + (k1 * r2);
  
  if (distortion > 1.2f) distortion = 1.2f;
  if (distortion < 0.8f) distortion = 0.8f;
  
  x = CX + (dx * distortion); 
  y = CY + (dy * distortion);
}

inline float OpenFIRE_Perspective::calculateQuadArea(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  return 0.5f * fabsf((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0));
}

// ========================================================================================
// 3. ESECUZIONE WARP (CUORE DEL SISTEMA)
// ========================================================================================

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  // Riorganizzazione Ciclica Pura: TL, TR, BR, BL 
  // Essenziale per creare un perimetro geometricamente inattaccabile
  float pt_TL_x = (float)x0,  pt_TL_y = (float)y0;
  float pt_TR_x = (float)x1,  pt_TR_y = (float)y1;
  float pt_BR_x = (float)x3,  pt_BR_y = (float)y3; 
  float pt_BL_x = (float)x2,  pt_BL_y = (float)y2; 

  applyLensCorrection(pt_TL_x, pt_TL_y);
  applyLensCorrection(pt_TR_x, pt_TR_y);
  applyLensCorrection(pt_BR_x, pt_BR_y);
  applyLensCorrection(pt_BL_x, pt_BL_y);

  float currentArea = calculateQuadArea(pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);

  if (!init) {
    if (currentArea > 10.0f) {
      // Inizializzazione della base dello schermo con lo stesso ordine ciclico (TL, TR, BR, BL)
      computeSquareToQuad(dstmatrix, dx0 * INV_NORM_SCALE, dy0 * INV_NORM_SCALE, 
                                     dx1 * INV_NORM_SCALE, dy1 * INV_NORM_SCALE, 
                                     dx3 * INV_NORM_SCALE, dy3 * INV_NORM_SCALE, 
                                     dx2 * INV_NORM_SCALE, dy2 * INV_NORM_SCALE);
      baseArea = currentArea;
      smoothedArea = currentArea;
      init = true;
    }
    if (!init) return; 
  }

  // Freeze di sicurezza: se la telecamera viene coperta (area vicina a 0), 
  // le variabili di parallasse non si degradano.
  if (currentArea > 10.0f) {
    smoothedArea = (0.1f * currentArea) + (0.9f * smoothedArea);
  }

  float dynamicSrcY = srcY;
  // Parallasse ottimizzato: hardware FPU assoluto
  if (parallaxFactor != 0.0f && baseArea > 10.0f && smoothedArea > 10.0f) {
    float distanceRatio = sqrtf(baseArea / smoothedArea);
    dynamicSrcY += parallaxFactor * (distanceRatio - 1.0f);
  }

  pt_TL_x *= INV_NORM_SCALE; pt_TL_y *= INV_NORM_SCALE;
  pt_TR_x *= INV_NORM_SCALE; pt_TR_y *= INV_NORM_SCALE;
  pt_BR_x *= INV_NORM_SCALE; pt_BR_y *= INV_NORM_SCALE;
  pt_BL_x *= INV_NORM_SCALE; pt_BL_y *= INV_NORM_SCALE;

  computeQuadToSquare(srcmatrix, pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  float normSrcX = srcX * INV_NORM_SCALE;
  float normSrcY = dynamicSrcY * INV_NORM_SCALE;

  float r0 = (normSrcX * warpmatrix[0] + normSrcY * warpmatrix[3] + warpmatrix[6]);
  float r1 = (normSrcX * warpmatrix[1] + normSrcY * warpmatrix[4] + warpmatrix[7]);
  float r3 = (normSrcX * warpmatrix[2] + normSrcY * warpmatrix[5] + warpmatrix[8]);
  
  if (fabsf(r3) > 1e-8f) {
    float invR3 = 1.0f / r3;
    float dstX_float = roundf((r0 * invR3) * NORM_SCALE);
    float dstY_float = roundf((r1 * invR3) * NORM_SCALE);

    // Clamp fisico: blocca le coordinate oltre lo schermo per prevenire crash UB 
    if (dstX_float > 2000000000.0f) dstX_float = 2000000000.0f;
    if (dstX_float < -2000000000.0f) dstX_float = -2000000000.0f;
    if (dstY_float > 2000000000.0f) dstY_float = 2000000000.0f;
    if (dstY_float < -2000000000.0f) dstY_float = -2000000000.0f;

    dstX = (int)dstX_float;
    dstY = (int)dstY_float;
  }
}

void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = adjustedX; srcY = adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED


#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// ========================================================================================
// 1. MOTORE ALGEBRICO (ESTREMA OTTIMIZZAZIONE DEI REGISTRI CPU)
// ========================================================================================

// Usato solo all'inizializzazione (dstmatrix)
inline void computeSquareToQuad(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float dx1 = x1 - x2; float dy1 = y1 - y2;
  float dx2 = x3 - x2; float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3; float sy = y0 - y1 + y2 - y3;

  float det = (dx1 * dy2 - dx2 * dy1);
  float g = 0.0f, h = 0.0f;
  
  if (fabsf(det) > 1e-8f) {
    float invDet = 1.0f / det;
    g = (sx * dy2 - dx2 * sy) * invDet;
    h = (dx1 * sy - sx * dy1) * invDet;
  }

  mat[0] = x1 - x0 + g * x1; mat[1] = y1 - y0 + g * y1; mat[2] = g;
  mat[3] = x3 - x0 + h * x3; mat[4] = y3 - y0 + h * y3; mat[5] = h;
  mat[6] = x0;               mat[7] = y0;               mat[8] = 1.0f;
}

// Usato ogni frame (srcmatrix). Lavora nei registri della CPU ed evita accessi inutili in RAM.
inline void computeQuadToSquare(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float dx1 = x1 - x2; float dy1 = y1 - y2;
  float dx2 = x3 - x2; float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3; float sy = y0 - y1 + y2 - y3;

  float det = (dx1 * dy2 - dx2 * dy1);
  float g = 0.0f, h = 0.0f;
  
  if (fabsf(det) > 1e-8f) {
    float invDet = 1.0f / det;
    g = (sx * dy2 - dx2 * sy) * invDet;
    h = (dx1 * sy - sx * dy1) * invDet;
  }

  float a = x1 - x0 + g * x1; float b = x3 - x0 + h * x3; float c = x0;
  float d = y1 - y0 + g * y1; float e = y3 - y0 + h * y3; float f = y0;

  // Inversione geometrica calcolata direttamente in memoria volatile
  float A = e - f * h; float B = c * h - b; float C = b * f - c * e;
  float D = f * g - d; float E = a - c * g; float F = c * d - a * f;
  float G = d * h - e * g; float H = b * g - a * h; float I = a * e - b * d;
  
  float invDetQuad = (a * A + b * D + c * G);
  if (fabsf(invDetQuad) > 1e-8f) {
    float idet = 1.0f / invDetQuad;
    // Unica scrittura nella memoria RAM
    mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
    mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
    mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
  }
}

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
// 2. FISICA: LENTE E AREA
// ========================================================================================

inline void OpenFIRE_Perspective::applyLensCorrection(float &x, float &y) {
  if (k1 == 0.0f) return; 
  
  float dx = x - CX; 
  float dy = y - CY;
  
  float nx = dx * INV_CX;
  float ny = dy * INV_CY;
  
  float r2 = (nx * nx) + (ny * ny);
  float distortion = 1.0f + (k1 * r2);
  
  if (distortion > 1.2f) distortion = 1.2f;
  if (distortion < 0.8f) distortion = 0.8f;
  
  x = CX + (dx * distortion); 
  y = CY + (dy * distortion);
}

inline float OpenFIRE_Perspective::calculateQuadArea(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  return 0.5f * fabsf((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0));
}

// ========================================================================================
// 3. ESECUZIONE WARP (CUORE DEL SISTEMA)
// ========================================================================================

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  float pt_TL_x = (float)x0,  pt_TL_y = (float)y0;
  float pt_TR_x = (float)x1,  pt_TR_y = (float)y1;
  float pt_BR_x = (float)x3,  pt_BR_y = (float)y3; 
  float pt_BL_x = (float)x2,  pt_BL_y = (float)y2; 

  applyLensCorrection(pt_TL_x, pt_TL_y);
  applyLensCorrection(pt_TR_x, pt_TR_y);
  applyLensCorrection(pt_BR_x, pt_BR_y);
  applyLensCorrection(pt_BL_x, pt_BL_y);

  float currentArea = calculateQuadArea(pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);

  if (!init) {
    if (currentArea > 10.0f) {
      computeSquareToQuad(dstmatrix, dx0 * INV_NORM_SCALE, dy0 * INV_NORM_SCALE, 
                                     dx1 * INV_NORM_SCALE, dy1 * INV_NORM_SCALE, 
                                     dx3 * INV_NORM_SCALE, dy3 * INV_NORM_SCALE, 
                                     dx2 * INV_NORM_SCALE, dy2 * INV_NORM_SCALE);
      baseArea = currentArea;
      smoothedArea = currentArea;
      init = true;
    }
    if (!init) return; 
  }

  if (currentArea > 10.0f) {
    smoothedArea = (0.1f * currentArea) + (0.9f * smoothedArea);
  }

  float dynamicSrcY = srcY;
  // Blocco Short-Circuit: se parallaxFactor è 0, tutto questo viene saltato.
  if (parallaxFactor != 0.0f && baseArea > 10.0f && smoothedArea > 10.0f) {
    // Una singola divisione diretta passata alla FPU hardware, seguita dalla radice veloce
    float distanceRatio = sqrtf(baseArea / smoothedArea);
    dynamicSrcY += parallaxFactor * (distanceRatio - 1.0f);
  }

  pt_TL_x *= INV_NORM_SCALE; pt_TL_y *= INV_NORM_SCALE;
  pt_TR_x *= INV_NORM_SCALE; pt_TR_y *= INV_NORM_SCALE;
  pt_BR_x *= INV_NORM_SCALE; pt_BR_y *= INV_NORM_SCALE;
  pt_BL_x *= INV_NORM_SCALE; pt_BL_y *= INV_NORM_SCALE;

  computeQuadToSquare(srcmatrix, pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  float normSrcX = srcX * INV_NORM_SCALE;
  float normSrcY = dynamicSrcY * INV_NORM_SCALE;

  float r0 = (normSrcX * warpmatrix[0] + normSrcY * warpmatrix[3] + warpmatrix[6]);
  float r1 = (normSrcX * warpmatrix[1] + normSrcY * warpmatrix[4] + warpmatrix[7]);
  float r3 = (normSrcX * warpmatrix[2] + normSrcY * warpmatrix[5] + warpmatrix[8]);
  
  if (fabsf(r3) > 1e-8f) {
    float invR3 = 1.0f / r3;
    float dstX_float = roundf((r0 * invR3) * NORM_SCALE);
    float dstY_float = roundf((r1 * invR3) * NORM_SCALE);

    if (dstX_float > 2000000000.0f) dstX_float = 2000000000.0f;
    if (dstX_float < -2000000000.0f) dstX_float = -2000000000.0f;
    if (dstY_float > 2000000000.0f) dstY_float = 2000000000.0f;
    if (dstY_float < -2000000000.0f) dstY_float = -2000000000.0f;

    dstX = (int)dstX_float;
    dstY = (int)dstY_float;
  }
}

void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = adjustedX; srcY = adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED
#endif


#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// ========================================================================================
// 1. MOTORE ALGEBRICO (OTTIMIZZATO PER FPU HARDWARE ESP32)
// ========================================================================================

inline void computeSquareToQuad(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float dx1 = x1 - x2; float dy1 = y1 - y2;
  float dx2 = x3 - x2; float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3; float sy = y0 - y1 + y2 - y3;

  float det = (dx1 * dy2 - dx2 * dy1);
  float g = 0.0f, h = 0.0f;
  
  if (fabsf(det) > 1e-8f) {
    // Sostituita la doppia divisione lenta con una moltiplicazione per il reciproco
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
}

inline void computeQuadToSquare(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  computeSquareToQuad(mat, x0, y0, x1, y1, x2, y2, x3, y3);

  float a = mat[0], d = mat[1], g = mat[2];
  float b = mat[3], e = mat[4], h = mat[5];
  float c = mat[6], f = mat[7];

  float A = e - f * h; float B = c * h - b; float C = b * f - c * e;
  float D = f * g - d; float E = a - c * g; float F = c * d - a * f;
  float G = d * h - e * g; float H = b * g - a * h; float I = a * e - b * d;
  
  float det = (a * A + b * D + c * G);
  if (fabsf(det) > 1e-8f) {
    float idet = 1.0f / det;
    mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
    mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
    mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
  }
}

// LOOP UNROLLING: Rimosso il ciclo "for" annidato. Pipeline CPU al 100%.
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
// 2. FISICA: LENTE E AREA
// ========================================================================================

inline void OpenFIRE_Perspective::applyLensCorrection(float &x, float &y) {
  if (k1 == 0.0f) return; 
  
  float dx = x - CX; 
  float dy = y - CY;
  
  // Moltiplicazioni invece di divisioni per normalizzare
  float nx = dx * INV_CX;
  float ny = dy * INV_CY;
  
  float r2 = (nx * nx) + (ny * ny);
  float distortion = 1.0f + (k1 * r2);
  
  if (distortion > 1.2f) distortion = 1.2f;
  if (distortion < 0.8f) distortion = 0.8f;
  
  x = CX + (dx * distortion); 
  y = CY + (dy * distortion);
}

inline float OpenFIRE_Perspective::calculateQuadArea(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  return 0.5f * fabsf((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0));
}

// ========================================================================================
// 3. ESECUZIONE WARP
// ========================================================================================

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  float pt_TL_x = (float)x0,  pt_TL_y = (float)y0;
  float pt_TR_x = (float)x1,  pt_TR_y = (float)y1;
  float pt_BR_x = (float)x3,  pt_BR_y = (float)y3; 
  float pt_BL_x = (float)x2,  pt_BL_y = (float)y2; 

  float scr_TL_x = dx0, scr_TL_y = dy0;
  float scr_TR_x = dx1, scr_TR_y = dy1;
  float scr_BR_x = dx3, scr_BR_y = dy3; 
  float scr_BL_x = dx2, scr_BL_y = dy2; 

  applyLensCorrection(pt_TL_x, pt_TL_y);
  applyLensCorrection(pt_TR_x, pt_TR_y);
  applyLensCorrection(pt_BR_x, pt_BR_y);
  applyLensCorrection(pt_BL_x, pt_BL_y);

  float currentArea = calculateQuadArea(pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);

  if (!init) {
    if (currentArea > 10.0f) {
      // Normalizzazione via moltiplicazione (velocissima)
      computeSquareToQuad(dstmatrix, scr_TL_x * INV_NORM_SCALE, scr_TL_y * INV_NORM_SCALE, 
                                     scr_TR_x * INV_NORM_SCALE, scr_TR_y * INV_NORM_SCALE, 
                                     scr_BR_x * INV_NORM_SCALE, scr_BR_y * INV_NORM_SCALE, 
                                     scr_BL_x * INV_NORM_SCALE, scr_BL_y * INV_NORM_SCALE);
      baseArea = currentArea;
      smoothedArea = currentArea;
      init = true;
    }
    if (!init) return; 
  }

  if (currentArea > 10.0f) {
    smoothedArea = (0.1f * currentArea) + (0.9f * smoothedArea);
  }

  float dynamicSrcY = srcY;
  if (parallaxFactor != 0.0f && baseArea > 10.0f && smoothedArea > 10.0f) {
    // Evitata divisione per area smussata
    float invSmoothed = 1.0f / smoothedArea;
    float distanceRatio = sqrtf(baseArea * invSmoothed);
    dynamicSrcY += parallaxFactor * (distanceRatio - 1.0f);
  }

  pt_TL_x *= INV_NORM_SCALE; pt_TL_y *= INV_NORM_SCALE;
  pt_TR_x *= INV_NORM_SCALE; pt_TR_y *= INV_NORM_SCALE;
  pt_BR_x *= INV_NORM_SCALE; pt_BR_y *= INV_NORM_SCALE;
  pt_BL_x *= INV_NORM_SCALE; pt_BL_y *= INV_NORM_SCALE;

  computeQuadToSquare(srcmatrix, pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  float normSrcX = srcX * INV_NORM_SCALE;
  float normSrcY = dynamicSrcY * INV_NORM_SCALE;

  float r0 = (normSrcX * warpmatrix[0] + normSrcY * warpmatrix[3] + warpmatrix[6]);
  float r1 = (normSrcX * warpmatrix[1] + normSrcY * warpmatrix[4] + warpmatrix[7]);
  float r3 = (normSrcX * warpmatrix[2] + normSrcY * warpmatrix[5] + warpmatrix[8]);
  
  if (fabsf(r3) > 1e-8f) {
    // Convertita la divisione finale
    float invR3 = 1.0f / r3;
    float dstX_float = roundf((r0 * invR3) * NORM_SCALE);
    float dstY_float = roundf((r1 * invR3) * NORM_SCALE);

    if (dstX_float > 2000000000.0f) dstX_float = 2000000000.0f;
    if (dstX_float < -2000000000.0f) dstX_float = -2000000000.0f;
    if (dstY_float > 2000000000.0f) dstY_float = 2000000000.0f;
    if (dstY_float < -2000000000.0f) dstY_float = -2000000000.0f;

    dstX = (int)dstX_float;
    dstY = (int)dstY_float;
  }
}

void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = adjustedX; srcY = adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// ========================================================================================
// 1. MOTORE ALGEBRICO (FLOAT 32-BIT HARDWARE-ACCELERATO)
// ========================================================================================

void computeSquareToQuad(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float dx1 = x1 - x2; float dy1 = y1 - y2;
  float dx2 = x3 - x2; float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3; float sy = y0 - y1 + y2 - y3;

  float det = (dx1 * dy2 - dx2 * dy1);
  float g = 0.0f, h = 0.0f;
  
  if (fabsf(det) > 1e-8f) {
    g = (sx * dy2 - dx2 * sy) / det;
    h = (dx1 * sy - sx * dy1) / det;
  }

  float a = x1 - x0 + g * x1; float b = x3 - x0 + h * x3; float c = x0;
  float d = y1 - y0 + g * y1; float e = y3 - y0 + h * y3; float f = y0;

  mat[0] = a; mat[1] = d; mat[2] = g;
  mat[3] = b; mat[4] = e; mat[5] = h;
  mat[6] = c; mat[7] = f; mat[8] = 1.0f;
}

void computeQuadToSquare(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  computeSquareToQuad(mat, x0, y0, x1, y1, x2, y2, x3, y3);

  float a = mat[0], d = mat[1], g = mat[2];
  float b = mat[3], e = mat[4], h = mat[5];
  float c = mat[6], f = mat[7];

  float A = e - f * h; float B = c * h - b; float C = b * f - c * e;
  float D = f * g - d; float E = a - c * g; float F = c * d - a * f;
  float G = d * h - e * g; float H = b * g - a * h; float I = a * e - b * d;
  
  float det = (a * A + b * D + c * G);
  float idet = 0.0f;
  if (fabsf(det) > 1e-8f) idet = 1.0f / det;

  mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
  mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
  mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
}

void multMats(float* a, float* b, float* res) {
  for (int r = 0; r < 3; r++) {
    int ri = r * 3;
    for (int c = 0; c < 3; c++) {
      res[ri + c] = (a[ri + 0] * b[0 + c] + a[ri + 1] * b[3 + c] + a[ri + 2] * b[6 + c]);
    }
  }
}

// ========================================================================================
// 2. FISICA: LENTE E AREA
// ========================================================================================

void OpenFIRE_Perspective::applyLensCorrection(float &x, float &y) {
  if (k1 == 0.0f) return; 
  
  float dx = x - centerX; 
  float dy = y - centerY;
  
  float nx = dx / centerX;
  float ny = dy / centerY;
  
  float r2 = (nx * nx) + (ny * ny);
  float distortion = 1.0f + (k1 * r2);
  
  // Clamp di sicurezza
  if (distortion > 1.2f) distortion = 1.2f;
  if (distortion < 0.8f) distortion = 0.8f;
  
  x = centerX + (dx * distortion); 
  y = centerY + (dy * distortion);
}

float OpenFIRE_Perspective::calculateQuadArea(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  return 0.5f * fabsf((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0));
}

// ========================================================================================
// 3. ESECUZIONE WARP (CON NORMALIZZAZIONE SCALATA)
// ========================================================================================

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  // A. Estrazione e Riordino Convesso (TL, TR, BR, BL)
  float pt_TL_x = (float)x0,  pt_TL_y = (float)y0;
  float pt_TR_x = (float)x1,  pt_TR_y = (float)y1;
  float pt_BR_x = (float)x3,  pt_BR_y = (float)y3; 
  float pt_BL_x = (float)x2,  pt_BL_y = (float)y2; 

  float scr_TL_x = dx0, scr_TL_y = dy0;
  float scr_TR_x = dx1, scr_TR_y = dy1;
  float scr_BR_x = dx3, scr_BR_y = dy3; 
  float scr_BL_x = dx2, scr_BL_y = dy2; 

  // B. Correzione Lente (scala reale)
  applyLensCorrection(pt_TL_x, pt_TL_y);
  applyLensCorrection(pt_TR_x, pt_TR_y);
  applyLensCorrection(pt_BR_x, pt_BR_y);
  applyLensCorrection(pt_BL_x, pt_BL_y);

  // C. Area e Parallasse (scala reale)
  float currentArea = calculateQuadArea(pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);

  if (!init) {
    if (currentArea > 10.0f) {
      // Inizializza la matrice dello schermo normalizzandola!
      computeSquareToQuad(dstmatrix, scr_TL_x / NORM_SCALE, scr_TL_y / NORM_SCALE, 
                                     scr_TR_x / NORM_SCALE, scr_TR_y / NORM_SCALE, 
                                     scr_BR_x / NORM_SCALE, scr_BR_y / NORM_SCALE, 
                                     scr_BL_x / NORM_SCALE, scr_BL_y / NORM_SCALE);
      baseArea = currentArea;
      smoothedArea = currentArea;
      init = true;
    }
    if (!init) return; 
  }

  if (currentArea > 10.0f) {
    smoothedArea = (0.1f * currentArea) + (0.9f * smoothedArea);
  }

  float dynamicSrcY = srcY;
  if (parallaxFactor != 0.0f && baseArea > 10.0f && smoothedArea > 10.0f) {
    float distanceRatio = sqrtf(baseArea / smoothedArea);
    dynamicSrcY += parallaxFactor * (distanceRatio - 1.0f);
  }

  // D. NORMALIZZAZIONE DEI DATI CAMERA
  pt_TL_x /= NORM_SCALE; pt_TL_y /= NORM_SCALE;
  pt_TR_x /= NORM_SCALE; pt_TR_y /= NORM_SCALE;
  pt_BR_x /= NORM_SCALE; pt_BR_y /= NORM_SCALE;
  pt_BL_x /= NORM_SCALE; pt_BL_y /= NORM_SCALE;

  // E. MOTORE PROSPETTICO SUI DATI NORMALIZZATI
  computeQuadToSquare(srcmatrix, pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  // Il punto del mirino deve essere anch'esso normalizzato!
  float normSrcX = srcX / NORM_SCALE;
  float normSrcY = dynamicSrcY / NORM_SCALE;

  float r0 = (normSrcX * warpmatrix[0] + normSrcY * warpmatrix[3] + warpmatrix[6]);
  float r1 = (normSrcX * warpmatrix[1] + normSrcY * warpmatrix[4] + warpmatrix[7]);
  float r3 = (normSrcX * warpmatrix[2] + normSrcY * warpmatrix[5] + warpmatrix[8]);
  
  // F. DE-NORMALIZZAZIONE E SICUREZZA
  if (fabsf(r3) > 1e-8f) {
    // Ricalcoliamo la coordinata finale e la rimoltiplichiamo per 10000.0f
    float dstX_float = roundf((r0 / r3) * NORM_SCALE);
    float dstY_float = roundf((r1 / r3) * NORM_SCALE);

    // Blocco anti-overflow per C++ (Clamp fisico a 2 miliardi)
    if (dstX_float > 2000000000.0f) dstX_float = 2000000000.0f;
    if (dstX_float < -2000000000.0f) dstX_float = -2000000000.0f;
    if (dstY_float > 2000000000.0f) dstY_float = 2000000000.0f;
    if (dstY_float < -2000000000.0f) dstY_float = -2000000000.0f;

    dstX = (int)dstX_float;
    dstY = (int)dstY_float;
  }
}

void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = adjustedX; srcY = adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// ========================================================================================
// 1. MOTORE ALGEBRICO E MATRICI (64-BIT)
// ========================================================================================

void computeSquareToQuad(double* mat, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  double dx1 = x1 - x2; double dy1 = y1 - y2;
  double dx2 = x3 - x2; double dy2 = y3 - y2;
  double sx = x0 - x1 + x2 - x3; double sy = y0 - y1 + y2 - y3;

  double det = (dx1 * dy2 - dx2 * dy1);
  double g = 0.0, h = 0.0;
  
  if (fabs(det) > 1e-12) {
    g = (sx * dy2 - dx2 * sy) / det;
    h = (dx1 * sy - sx * dy1) / det;
  }

  double a = x1 - x0 + g * x1; double b = x3 - x0 + h * x3; double c = x0;
  double d = y1 - y0 + g * y1; double e = y3 - y0 + h * y3; double f = y0;

  mat[0] = a; mat[1] = d; mat[2] = g;
  mat[3] = b; mat[4] = e; mat[5] = h;
  mat[6] = c; mat[7] = f; mat[8] = 1.0;
}

void computeQuadToSquare(double* mat, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  computeSquareToQuad(mat, x0, y0, x1, y1, x2, y2, x3, y3);

  double a = mat[0], d = mat[1], g = mat[2];
  double b = mat[3], e = mat[4], h = mat[5];
  double c = mat[6], f = mat[7];

  double A = e - f * h; double B = c * h - b; double C = b * f - c * e;
  double D = f * g - d; double E = a - c * g; double F = c * d - a * f;
  double G = d * h - e * g; double H = b * g - a * h; double I = a * e - b * d;
  
  double det = (a * A + b * D + c * G);
  double idet = 0.0;
  if (fabs(det) > 1e-12) idet = 1.0 / det;

  mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
  mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
  mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
}

void multMats(double* a, double* b, double* res) {
  for (int r = 0; r < 3; r++) {
    int ri = r * 3;
    for (int c = 0; c < 3; c++) {
      res[ri + c] = (a[ri + 0] * b[0 + c] + a[ri + 1] * b[3 + c] + a[ri + 2] * b[6 + c]);
    }
  }
}

// ========================================================================================
// 2. FISICA: LENTE E AREA
// ========================================================================================

void OpenFIRE_Perspective::applyLensCorrection(double &x, double &y) {
  if (k1 == 0.0f) return; 
  
  double dx = x - centerX; 
  double dy = y - centerY;
  
  double nx = dx / centerX;
  double ny = dy / centerY;
  
  double r2 = (nx * nx) + (ny * ny);
  double distortion = 1.0 + ((double)k1 * r2);
  
  if (distortion > 1.2) distortion = 1.2;
  if (distortion < 0.8) distortion = 0.8;
  
  x = centerX + (dx * distortion); 
  y = centerY + (dy * distortion);
}

// Calcolo Area per poligono circolare convesso (TL, TR, BR, BL)
double OpenFIRE_Perspective::calculateQuadArea(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  return 0.5 * fabs((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0));
}

// ========================================================================================
// 3. ESECUZIONE WARP (LOGICA PRINCIPALE)
// ========================================================================================

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  // A. ISOLAMENTO E RIORDINO CICLICO (TL, TR, BR, BL)
  // Trasformiamo i parametri (TL, TR, BL, BR) in un quadrato perfetto convesso.
  double pt_TL_x = (double)x0,  pt_TL_y = (double)y0;
  double pt_TR_x = (double)x1,  pt_TR_y = (double)y1;
  double pt_BR_x = (double)x3,  pt_BR_y = (double)y3; // Scambiato
  double pt_BL_x = (double)x2,  pt_BL_y = (double)y2; // Scambiato

  double scr_TL_x = (double)dx0, scr_TL_y = (double)dy0;
  double scr_TR_x = (double)dx1, scr_TR_y = (double)dy1;
  double scr_BR_x = (double)dx3, scr_BR_y = (double)dy3; // Scambiato
  double scr_BL_x = (double)dx2, scr_BL_y = (double)dy2; // Scambiato

  // B. CORREZIONE LENTE SUI PUNTI CAMERA
  applyLensCorrection(pt_TL_x, pt_TL_y);
  applyLensCorrection(pt_TR_x, pt_TR_y);
  applyLensCorrection(pt_BR_x, pt_BR_y);
  applyLensCorrection(pt_BL_x, pt_BL_y);

  // C. CALCOLO AREA (Geometria protetta)
  double currentArea = calculateQuadArea(pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);

  // D. INIZIALIZZAZIONE SICURA
  if (!init) {
    if (currentArea > 10.0) {
      computeSquareToQuad(dstmatrix, scr_TL_x, scr_TL_y, scr_TR_x, scr_TR_y, scr_BR_x, scr_BR_y, scr_BL_x, scr_BL_y);
      baseArea = currentArea;
      smoothedArea = currentArea;
      init = true;
    }
    if (!init) return; 
  }

  // E. FILTRO AREA E PARALLASSE
  if (currentArea > 10.0) {
    smoothedArea = (0.1 * currentArea) + (0.9 * smoothedArea);
  }

  double dynamicSrcY = srcY;
  if (parallaxFactor != 0.0f && baseArea > 10.0 && smoothedArea > 10.0) {
    double distanceRatio = sqrt(baseArea / smoothedArea);
    dynamicSrcY += (double)parallaxFactor * (distanceRatio - 1.0);
  }

  // F. MOTORE PROSPETTICO
  computeQuadToSquare(srcmatrix, pt_TL_x, pt_TL_y, pt_TR_x, pt_TR_y, pt_BR_x, pt_BR_y, pt_BL_x, pt_BL_y);
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  double r0 = (srcX * warpmatrix[0] + dynamicSrcY * warpmatrix[3] + warpmatrix[6]);
  double r1 = (srcX * warpmatrix[1] + dynamicSrcY * warpmatrix[4] + warpmatrix[7]);
  double r3 = (srcX * warpmatrix[2] + dynamicSrcY * warpmatrix[5] + warpmatrix[8]);
  
  // G. PROTEZIONE OUTPUT ESTREMI (Evita crash C++ Undefined Behavior per divisione int/float off-screen)
  if (fabs(r3) > 1e-12) {
    double dstX_double = round(r0 / r3);
    double dstY_double = round(r1 / r3);

    // Blocco fisico a 2 miliardi (INT_MAX è 2.14 miliardi)
    if (dstX_double > 2000000000.0) dstX_double = 2000000000.0;
    if (dstX_double < -2000000000.0) dstX_double = -2000000000.0;
    if (dstY_double > 2000000000.0) dstY_double = 2000000000.0;
    if (dstY_double < -2000000000.0) dstY_double = -2000000000.0;

    dstX = (int)dstX_double;
    dstY = (int)dstY_double;
  }
}

// Metodi Accessori
void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = (double)adjustedX; srcY = (double)adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// ========================================================================================
// 1. MOTORE ALGEBRICO (ROW-MAJOR @ 64-BIT)
// ========================================================================================

void computeSquareToQuad(double* mat, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  double dx1 = x1 - x2; double dy1 = y1 - y2;
  double dx2 = x3 - x2; double dy2 = y3 - y2;
  double sx = x0 - x1 + x2 - x3; double sy = y0 - y1 + y2 - y3;

  double det = (dx1 * dy2 - dx2 * dy1);
  double g = 0.0, h = 0.0;
  
  if (fabs(det) > 1e-12) {
    g = (sx * dy2 - dx2 * sy) / det;
    h = (dx1 * sy - sx * dy1) / det;
  }

  double a = x1 - x0 + g * x1; double b = x3 - x0 + h * x3; double c = x0;
  double d = y1 - y0 + g * y1; double e = y3 - y0 + h * y3; double f = y0;

  mat[0] = a; mat[1] = d; mat[2] = g;
  mat[3] = b; mat[4] = e; mat[5] = h;
  mat[6] = c; mat[7] = f; mat[8] = 1.0;
}

void computeQuadToSquare(double* mat, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  computeSquareToQuad(mat, x0, y0, x1, y1, x2, y2, x3, y3);

  double a = mat[0], d = mat[1], g = mat[2];
  double b = mat[3], e = mat[4], h = mat[5];
  double c = mat[6], f = mat[7];

  double A = e - f * h; double B = c * h - b; double C = b * f - c * e;
  double D = f * g - d; double E = a - c * g; double F = c * d - a * f;
  double G = d * h - e * g; double H = b * g - a * h; double I = a * e - b * d;
  
  double det = (a * A + b * D + c * G);
  double idet = 0.0;
  if (fabs(det) > 1e-12) idet = 1.0 / det;

  mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
  mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
  mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
}

void multMats(double* a, double* b, double* res) {
  for (int r = 0; r < 3; r++) {
    int ri = r * 3;
    for (int c = 0; c < 3; c++) {
      res[ri + c] = (a[ri + 0] * b[0 + c] + a[ri + 1] * b[3 + c] + a[ri + 2] * b[6 + c]);
    }
  }
}

// ========================================================================================
// 2. FISICA, LENTE (RISOLUZIONE-INDIPENDENTE) E CALCOLO AREA
// ========================================================================================

void OpenFIRE_Perspective::applyLensCorrection(double &x, double &y) {
  if (k1 == 0.0f) return; 
  
  // Distanza dal centro
  double dx = x - centerX; 
  double dy = y - centerY;
  
  // Normalizzazione: trasformiamo la distanza in un rapporto (0.0 al centro, ~1.0 ai bordi)
  double nx = dx / centerX;
  double ny = dy / centerY;
  
  double r2 = (nx * nx) + (ny * ny);
  double distortion = 1.0 + ((double)k1 * r2);
  
  // Freno fisico (Clamp): previene implosioni matematiche per punti stimati molto fuori schermo
  if (distortion > 1.2) distortion = 1.2;
  if (distortion < 0.8) distortion = 0.8;
  
  // Riapplichiamo la distorsione alla scala originale
  x = centerX + (dx * distortion); 
  y = centerY + (dy * distortion);
}

double OpenFIRE_Perspective::calculateQuadArea(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  return 0.5 * fabs((x3 - x0) * (y2 - y1) - (x2 - x1) * (y3 - y0));
}

// ========================================================================================
// 3. ESECUZIONE WARP
// ========================================================================================

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  double fx0 = (double)x0, fy0 = (double)y0;
  double fx1 = (double)x1, fy1 = (double)y1;
  double fx2 = (double)x2, fy2 = (double)y2;
  double fx3 = (double)x3, fy3 = (double)y3;

  applyLensCorrection(fx0, fy0);
  applyLensCorrection(fx1, fy1);
  applyLensCorrection(fx2, fy2);
  applyLensCorrection(fx3, fy3);

  double currentArea = calculateQuadArea(fx0, fy0, fx1, fy1, fx2, fy2, fx3, fy3);

  if (!init) {
    if (currentArea > 10.0) {
      computeSquareToQuad(dstmatrix, (double)dx0, (double)dy0, (double)dx1, (double)dy1, (double)dx2, (double)dy2, (double)dx3, (double)dy3);
      baseArea = currentArea;
      smoothedArea = currentArea;
      init = true;
    }
    if (!init) return; 
  }

  if (currentArea > 10.0) {
    smoothedArea = (0.1 * currentArea) + (0.9 * smoothedArea);
  }

  double dynamicSrcY = srcY;
  if (parallaxFactor != 0.0f && baseArea > 10.0 && smoothedArea > 10.0) {
    double distanceRatio = sqrt(baseArea / smoothedArea);
    dynamicSrcY += (double)parallaxFactor * (distanceRatio - 1.0);
  }

  computeQuadToSquare(srcmatrix, fx0, fy0, fx1, fy1, fx2, fy2, fx3, fy3);
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  double r0 = (srcX * warpmatrix[0] + dynamicSrcY * warpmatrix[3] + warpmatrix[6]);
  double r1 = (srcX * warpmatrix[1] + dynamicSrcY * warpmatrix[4] + warpmatrix[7]);
  double r3 = (srcX * warpmatrix[2] + dynamicSrcY * warpmatrix[5] + warpmatrix[8]);
  
  if (fabs(r3) > 1e-12) {
    dstX = (int)round(r0 / r3);
    dstY = (int)round(r1 / r3);
  }
}

void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = (double)adjustedX; srcY = (double)adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// ========================================================================================
// 1. MOTORE ALGEBRICO ORIGINALE (ROW-MAJOR @ 64-BIT)
// ========================================================================================

void computeSquareToQuad(double* mat, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  double dx1 = x1 - x2; double dy1 = y1 - y2;
  double dx2 = x3 - x2; double dy2 = y3 - y2;
  double sx = x0 - x1 + x2 - x3; double sy = y0 - y1 + y2 - y3;

  double det = (dx1 * dy2 - dx2 * dy1);
  double g = 0.0, h = 0.0;
  
  if (fabs(det) > 1e-12) {
    g = (sx * dy2 - dx2 * sy) / det;
    h = (dx1 * sy - sx * dy1) / det;
  }

  double a = x1 - x0 + g * x1; double b = x3 - x0 + h * x3; double c = x0;
  double d = y1 - y0 + g * y1; double e = y3 - y0 + h * y3; double f = y0;

  mat[0] = a; mat[1] = d; mat[2] = g;
  mat[3] = b; mat[4] = e; mat[5] = h;
  mat[6] = c; mat[7] = f; mat[8] = 1.0;
}

void computeQuadToSquare(double* mat, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  computeSquareToQuad(mat, x0, y0, x1, y1, x2, y2, x3, y3);

  double a = mat[0], d = mat[1], g = mat[2];
  double b = mat[3], e = mat[4], h = mat[5];
  double c = mat[6], f = mat[7];

  double A = e - f * h; double B = c * h - b; double C = b * f - c * e;
  double D = f * g - d; double E = a - c * g; double F = c * d - a * f;
  double G = d * h - e * g; double H = b * g - a * h; double I = a * e - b * d;
  
  double det = (a * A + b * D + c * G);
  double idet = 0.0;
  if (fabs(det) > 1e-12) idet = 1.0 / det;

  mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
  mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
  mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
}

void multMats(double* a, double* b, double* res) {
  for (int r = 0; r < 3; r++) {
    int ri = r * 3;
    for (int c = 0; c < 3; c++) {
      res[ri + c] = (a[ri + 0] * b[0 + c] + a[ri + 1] * b[3 + c] + a[ri + 2] * b[6 + c]);
    }
  }
}

// ========================================================================================
// 2. FISICA, LENTE E CALCOLO AREA
// ========================================================================================

void OpenFIRE_Perspective::applyLensCorrection(double &x, double &y) {
  if (k1 == 0.0f) return; 
  double dx = x - centerX; 
  double dy = y - centerY;
  double r2 = (dx * dx) + (dy * dy);
  double distortion = 1.0 + ((double)k1 * r2);
  
  // FIX CRITICO: Freno fisico per i punti stimati estremamente fuori schermo.
  // Evita che la formula della lente disintegri la matrice se il punto è a x=15000.
  if (distortion > 1.2) distortion = 1.2;
  if (distortion < 0.8) distortion = 0.8;
  
  x = centerX + (dx * distortion); 
  y = centerY + (dy * distortion);
}

// Formula robusta per l'area calibrata sull'ordine TL(0), TR(1), BL(2), BR(3)
double OpenFIRE_Perspective::calculateQuadArea(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  return 0.5 * fabs((x3 - x0) * (y2 - y1) - (x2 - x1) * (y3 - y0));
}

// ========================================================================================
// 3. ESECUZIONE WARP
// ========================================================================================

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  // Iniezione diretta a 64-bit
  double fx0 = (double)x0, fy0 = (double)y0;
  double fx1 = (double)x1, fy1 = (double)y1;
  double fx2 = (double)x2, fy2 = (double)y2;
  double fx3 = (double)x3, fy3 = (double)y3;

  // Correzione Lente a monte (protetta dal Fix)
  applyLensCorrection(fx0, fy0);
  applyLensCorrection(fx1, fy1);
  applyLensCorrection(fx2, fy2);
  applyLensCorrection(fx3, fy3);

  // Calcolo Area geometrica per Parallasse
  double currentArea = calculateQuadArea(fx0, fy0, fx1, fy1, fx2, fy2, fx3, fy3);

  // Inizializzazione protetta: la camera deve inquadrare chiaramente prima di salvare la base
  if (!init) {
    if (currentArea > 10.0) {
      computeSquareToQuad(dstmatrix, (double)dx0, (double)dy0, (double)dx1, (double)dy1, (double)dx2, (double)dy2, (double)dx3, (double)dy3);
      baseArea = currentArea;
      smoothedArea = currentArea;
      init = true;
    }
    if (!init) return; 
  }

  // Filtro EMA (Ammortizzatore asse Z)
  if (currentArea > 10.0) {
    smoothedArea = (0.1 * currentArea) + (0.9 * smoothedArea);
  }

  // Calcolo Parallasse (Alzo mirino)
  double dynamicSrcY = srcY;
  if (parallaxFactor != 0.0f && baseArea > 10.0 && smoothedArea > 10.0) {
    double distanceRatio = sqrt(baseArea / smoothedArea);
    dynamicSrcY += (double)parallaxFactor * (distanceRatio - 1.0);
  }

  // Motore Matematico Principale
  computeQuadToSquare(srcmatrix, fx0, fy0, fx1, fy1, fx2, fy2, fx3, fy3);
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  // Vettore Post-Moltiplicato
  double r0 = (srcX * warpmatrix[0] + dynamicSrcY * warpmatrix[3] + warpmatrix[6]);
  double r1 = (srcX * warpmatrix[1] + dynamicSrcY * warpmatrix[4] + warpmatrix[7]);
  double r3 = (srcX * warpmatrix[2] + dynamicSrcY * warpmatrix[5] + warpmatrix[8]);
  
  // Sicurezza anti-crash e output
  if (fabs(r3) > 1e-12) {
    dstX = (int)round(r0 / r3);
    dstY = (int)round(r1 / r3);
  }
}

// Metodi Accessori
void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = (double)adjustedX; srcY = (double)adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// ========================================================================================
// IL TUO MOTORE ALGEBRICO ORIGINALE (Aggiornato a 64-bit)
// ========================================================================================

void computeSquareToQuad(double* mat, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  double dx1 = x1 - x2; double dy1 = y1 - y2;
  double dx2 = x3 - x2; double dy2 = y3 - y2;
  double sx = x0 - x1 + x2 - x3; double sy = y0 - y1 + y2 - y3;

  double det = (dx1 * dy2 - dx2 * dy1);
  double g = 0.0, h = 0.0;
  
  if (fabs(det) > 1e-12) {
    g = (sx * dy2 - dx2 * sy) / det;
    h = (dx1 * sy - sx * dy1) / det;
  }

  double a = x1 - x0 + g * x1; double b = x3 - x0 + h * x3; double c = x0;
  double d = y1 - y0 + g * y1; double e = y3 - y0 + h * y3; double f = y0;

  mat[0] = a; mat[1] = d; mat[2] = g;
  mat[3] = b; mat[4] = e; mat[5] = h;
  mat[6] = c; mat[7] = f; mat[8] = 1.0;
}

void computeQuadToSquare(double* mat, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  computeSquareToQuad(mat, x0, y0, x1, y1, x2, y2, x3, y3);

  double a = mat[0], d = mat[1], g = mat[2];
  double b = mat[3], e = mat[4], h = mat[5];
  double c = mat[6], f = mat[7];

  double A = e - f * h; double B = c * h - b; double C = b * f - c * e;
  double D = f * g - d; double E = a - c * g; double F = c * d - a * f;
  double G = d * h - e * g; double H = b * g - a * h; double I = a * e - b * d;
  
  double det = (a * A + b * D + c * G);
  double idet = 0.0;
  if (fabs(det) > 1e-12) idet = 1.0 / det;

  mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
  mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
  mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
}

void multMats(double* a, double* b, double* res) {
  for (int r = 0; r < 3; r++) {
    int ri = r * 3;
    for (int c = 0; c < 3; c++) {
      res[ri + c] = (a[ri + 0] * b[0 + c] + a[ri + 1] * b[3 + c] + a[ri + 2] * b[6 + c]);
    }
  }
}

// ========================================================================================
// FISICA E CALCOLO AREA (64-bit)
// ========================================================================================

void OpenFIRE_Perspective::applyLensCorrection(double &x, double &y) {
  if (k1 == 0.0f) return; 
  double dx = x - centerX; double dy = y - centerY;
  double r2 = (dx * dx) + (dy * dy);
  double distortion = 1.0 + ((double)k1 * r2);
  x = centerX + (dx * distortion); 
  y = centerY + (dy * distortion);
}

double OpenFIRE_Perspective::calculateQuadArea(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  return 0.5 * fabs((x3 - x0) * (y2 - y1) - (x2 - x1) * (y3 - y0));
}

// ========================================================================================
// ESECUZIONE WARP
// ========================================================================================

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  double fx0 = (double)x0, fy0 = (double)y0;
  double fx1 = (double)x1, fy1 = (double)y1;
  double fx2 = (double)x2, fy2 = (double)y2;
  double fx3 = (double)x3, fy3 = (double)y3;

  applyLensCorrection(fx0, fy0);
  applyLensCorrection(fx1, fy1);
  applyLensCorrection(fx2, fy2);
  applyLensCorrection(fx3, fy3);

  double currentArea = calculateQuadArea(fx0, fy0, fx1, fy1, fx2, fy2, fx3, fy3);

  if (!init) {
    if (currentArea > 10.0) {
      computeSquareToQuad(dstmatrix, (double)dx0, (double)dy0, (double)dx1, (double)dy1, (double)dx2, (double)dy2, (double)dx3, (double)dy3);
      baseArea = currentArea;
      smoothedArea = currentArea;
      init = true;
    }
    if (!init) return; 
  }

  if (currentArea > 10.0) {
    smoothedArea = (0.1 * currentArea) + (0.9 * smoothedArea);
  }

  double dynamicSrcY = srcY;
  if (parallaxFactor != 0.0f && baseArea > 10.0 && smoothedArea > 10.0) {
    double distanceRatio = sqrt(baseArea / smoothedArea);
    dynamicSrcY += (double)parallaxFactor * (distanceRatio - 1.0);
  }

  computeQuadToSquare(srcmatrix, fx0, fy0, fx1, fy1, fx2, fy2, fx3, fy3);
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  double r0 = (srcX * warpmatrix[0] + dynamicSrcY * warpmatrix[3] + warpmatrix[6]);
  double r1 = (srcX * warpmatrix[1] + dynamicSrcY * warpmatrix[4] + warpmatrix[7]);
  double r3 = (srcX * warpmatrix[2] + dynamicSrcY * warpmatrix[5] + warpmatrix[8]);
  
  if (fabs(r3) > 1e-12) {
    dstX = (int)round(r0 / r3);
    dstY = (int)round(r1 / r3);
  }
}

void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = (double)adjustedX; srcY = (double)adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// ========================================================================================
// 1. IL MOTORE ORIGINALE (ROW-MAJOR) - INTATTO E INVIOLATO
// ========================================================================================

void computeSquareToQuad(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float dx1 = x1 - x2; float dy1 = y1 - y2;
  float dx2 = x3 - x2; float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3; float sy = y0 - y1 + y2 - y3;

  float det = (dx1 * dy2 - dx2 * dy1);
  float g = 0.0f, h = 0.0f;
  
  if (fabsf(det) > 1e-6f) {
    g = (sx * dy2 - dx2 * sy) / det;
    h = (dx1 * sy - sx * dy1) / det;
  }

  float a = x1 - x0 + g * x1; float b = x3 - x0 + h * x3; float c = x0;
  float d = y1 - y0 + g * y1; float e = y3 - y0 + h * y3; float f = y0;

  mat[0] = a; mat[1] = d; mat[2] = g;
  mat[3] = b; mat[4] = e; mat[5] = h;
  mat[6] = c; mat[7] = f; mat[8] = 1.0f;
}

void computeQuadToSquare(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  computeSquareToQuad(mat, x0, y0, x1, y1, x2, y2, x3, y3);

  float a = mat[0], d = mat[1], g = mat[2];
  float b = mat[3], e = mat[4], h = mat[5];
  float c = mat[6], f = mat[7];

  float A = e - f * h; float B = c * h - b; float C = b * f - c * e;
  float D = f * g - d; float E = a - c * g; float F = c * d - a * f;
  float G = d * h - e * g; float H = b * g - a * h; float I = a * e - b * d;
  
  float det = (a * A + b * D + c * G);
  float idet = 0.0f;
  if (fabsf(det) > 1e-6f) idet = 1.0f / det;

  mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
  mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
  mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
}

void multMats(float* a, float* b, float* res) {
  for (int r = 0; r < 3; r++) {
    int ri = r * 3;
    for (int c = 0; c < 3; c++) {
      res[ri + c] = (a[ri + 0] * b[0 + c] + a[ri + 1] * b[3 + c] + a[ri + 2] * b[6 + c]);
    }
  }
}

// ========================================================================================
// 2. FUNZIONI FISICHE
// ========================================================================================

void OpenFIRE_Perspective::applyLensCorrection(float &x, float &y) {
  if (k1 == 0.0f) return; 
  float dx = x - centerX; float dy = y - centerY;
  float r2 = (dx * dx) + (dy * dy);
  float distortion = 1.0f + (k1 * r2);
  x = centerX + (dx * distortion); 
  y = centerY + (dy * distortion);
}

// Calcolo dell'area specifico per il tuo ordine: 0(TL), 1(TR), 2(BL), 3(BR)
// Usa il prodotto incrociato delle vere diagonali (0->3 e 1->2) per essere invincibile
float OpenFIRE_Perspective::calculateQuadArea(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  return 0.5f * fabsf((x3 - x0) * (y2 - y1) - (x2 - x1) * (y3 - y0));
}

// ========================================================================================
// 3. ESECUZIONE (WARP)
// ========================================================================================

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  float fx0 = (float)x0, fy0 = (float)y0;
  float fx1 = (float)x1, fy1 = (float)y1;
  float fx2 = (float)x2, fy2 = (float)y2;
  float fx3 = (float)x3, fy3 = (float)y3;

  // 1. CORREZIONE LENTE A MONTE
  applyLensCorrection(fx0, fy0);
  applyLensCorrection(fx1, fy1);
  applyLensCorrection(fx2, fy2);
  applyLensCorrection(fx3, fy3);

  // 2. CALCOLO AREA (Per Parallasse)
  float currentArea = calculateQuadArea(fx0, fy0, fx1, fy1, fx2, fy2, fx3, fy3);

  // 3. INIZIALIZZAZIONE DELLA MATRICE SCHERMO E DELL'AREA BASE
  if (!init) {
    if (currentArea > 10.0f) {
      computeSquareToQuad(dstmatrix, dx0, dy0, dx1, dy1, dx2, dy2, dx3, dy3);
      baseArea = currentArea;
      smoothedArea = currentArea;
      init = true;
    }
    if (!init) return; // Aspetta finché non vede bene i LED
  }

  // 4. FILTRO EMA (Ammortizzatore per il Parallasse)
  if (currentArea > 10.0f) {
    smoothedArea = (0.1f * currentArea) + (0.9f * smoothedArea);
  }

  // 5. CALCOLO PARALLASSE (Spostamento dinamico del mirino sull'asse Y)
  float dynamicSrcY = srcY;
  if (parallaxFactor != 0.0f && baseArea > 10.0f && smoothedArea > 10.0f) {
    float distanceRatio = sqrtf(baseArea / smoothedArea);
    dynamicSrcY += parallaxFactor * (distanceRatio - 1.0f);
  }

  // 6. MOTORE MATEMATICO ORIGINALE
  computeQuadToSquare(srcmatrix, fx0, fy0, fx1, fy1, fx2, fy2, fx3, fy3);
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  // Vettore Post-Moltiplicato usando il dynamicSrcY
  float r0 = (srcX * warpmatrix[0] + dynamicSrcY * warpmatrix[3] + warpmatrix[6]);
  float r1 = (srcX * warpmatrix[1] + dynamicSrcY * warpmatrix[4] + warpmatrix[7]);
  float r3 = (srcX * warpmatrix[2] + dynamicSrcY * warpmatrix[5] + warpmatrix[8]);
  
  // 7. RESTITUZIONE DATI
  if (fabsf(r3) > 1e-6f) {
    // Usiamo roundf per centrare perfettamente il pixel
    dstX = (int)roundf(r0 / r3);
    dstY = (int)roundf(r1 / r3);
  }
}

// Setter e Getter
void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = adjustedX; srcY = adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// Algoritmo Heckbert puro: Mappatura dal Quadrato al Quadrilatero
void computeSquareToQuad(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float dx1 = x1 - x2;
  float dy1 = y1 - y2;
  float dx2 = x3 - x2;
  float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3;
  float sy = y0 - y1 + y2 - y3;

  float det = (dx1 * dy2 - dx2 * dy1);
  float g = 0.0f;
  float h = 0.0f;
  
  // Protezione base per non dividere per zero
  if (fabsf(det) > 1e-6f) {
    g = (sx * dy2 - dx2 * sy) / det;
    h = (dx1 * sy - sx * dy1) / det;
  }

  float a = x1 - x0 + g * x1;
  float b = x3 - x0 + h * x3;
  float c = x0;
  float d = y1 - y0 + g * y1;
  float e = y3 - y0 + h * y3;
  float f = y0;

  // Compressione 3x3 Row-Major esatta dell'array 4x4 originale
  mat[0] = a; mat[1] = d; mat[2] = g;
  mat[3] = b; mat[4] = e; mat[5] = h;
  mat[6] = c; mat[7] = f; mat[8] = 1.0f;
}

// Inversione geometrica
void computeQuadToSquare(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  
  computeSquareToQuad(mat, x0, y0, x1, y1, x2, y2, x3, y3);

  // Mappatura fedele agli indici della vecchia matrice
  float a = mat[0]; float d = mat[1]; float g = mat[2];
  float b = mat[3]; float e = mat[4]; float h = mat[5];
  float c = mat[6]; float f = mat[7];

  // Matrice dei cofattori originale (Traspone e inverte in un sol colpo)
  float A =     e - f * h;
  float B = c * h - b;
  float C = b * f - c * e;
  float D = f * g - d;
  float E =     a - c * g;
  float F = c * d - a * f;
  float G = d * h - e * g;
  float H = b * g - a * h;
  float I = a * e - b * d;
  
  float det = (a * A + b * D + c * G);
  float idet = 0.0f;
  if (fabsf(det) > 1e-6f) {
    idet = 1.0f / det;
  }

  mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
  mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
  mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
}

// Moltiplicazione Row-Major
void multMats(float* a, float* b, float* res) {
  for (int r = 0; r < 3; r++) {
    int ri = r * 3;
    for (int c = 0; c < 3; c++) {
      res[ri + c] = (
        a[ri + 0] * b[0 + c] +
        a[ri + 1] * b[3 + c] +
        a[ri + 2] * b[6 + c]
      );
    }
  }
}

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  if (!init) {
    computeSquareToQuad(dstmatrix, dx0, dy0, dx1, dy1, dx2, dy2, dx3, dy3);
    init = true;
  }
  
  // Cast identico all'originale
  computeQuadToSquare(srcmatrix, (float)x0, (float)y0, (float)x1, (float)y1, (float)x2, (float)y2, (float)x3, (float)y3);
  
  multMats(srcmatrix, dstmatrix, warpmatrix);
  
  // Vettore Post-Moltiplicato [srcX, srcY, 1] * Matrice
  // Questa è la chiave di volta che impedisce i salti col tuo ordine di punti
  float r0 = (srcX * warpmatrix[0] + srcY * warpmatrix[3] + warpmatrix[6]);
  float r1 = (srcX * warpmatrix[1] + srcY * warpmatrix[4] + warpmatrix[7]);
  float r3 = (srcX * warpmatrix[2] + srcY * warpmatrix[5] + warpmatrix[8]);
  
  if (fabsf(r3) > 1e-6f) {
    dstX = (int)(r0 / r3);
    dstY = (int)floorf(r1 / r3);
  }
}

void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) {
  srcX = adjustedX;
  srcY = adjustedY;
}

void OpenFIRE_Perspective::deinit (bool set) {
  init = set;
}

int OpenFIRE_Perspective::getX() {
  return dstX;
}

int OpenFIRE_Perspective::getY() {
  return dstY;
}

#endif // USE_PERSPECTIVE_ADVANCED
#endif


#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

// ========================================================================================
// RISOLUTORE 64-BIT: DIRECT LINEAR TRANSFORM (DLT)
// ========================================================================================
bool OpenFIRE_Perspective::computeHomography(double src[4][2], double dst[4][2], double H[9]) {
    double P[8][9] = {0.0};

    // Costruzione del sistema lineare
    for (int i = 0; i < 4; i++) {
        double x = src[i][0];
        double y = src[i][1];
        double u = dst[i][0];
        double v = dst[i][1];

        P[i*2][0] = x; P[i*2][1] = y; P[i*2][2] = 1.0;
        P[i*2][6] = -u * x; P[i*2][7] = -u * y; P[i*2][8] = u;

        P[i*2+1][3] = x; P[i*2+1][4] = y; P[i*2+1][5] = 1.0;
        P[i*2+1][6] = -v * x; P[i*2+1][7] = -v * y; P[i*2+1][8] = v;
    }

    // Risoluzione via Gauss-Jordan
    for (int i = 0; i < 8; i++) {
        int pivot = i;
        double maxVal = fabs(P[i][i]);
        
        for (int j = i + 1; j < 8; j++) {
            if (fabs(P[j][i]) > maxVal) {
                maxVal = fabs(P[j][i]);
                pivot = j;
            }
        }

        if (maxVal < 1e-12) return false; 

        if (pivot != i) {
            for (int j = i; j < 9; j++) {
                double temp = P[i][j];
                P[i][j] = P[pivot][j];
                P[pivot][j] = temp;
            }
        }

        double div = P[i][i];
        for (int j = i; j < 9; j++) P[i][j] /= div;

        for (int j = 0; j < 8; j++) {
            if (j != i) {
                double mul = P[j][i];
                for (int k = i; k < 9; k++) P[j][k] -= mul * P[i][k];
            }
        }
    }

    H[0] = P[0][8]; H[1] = P[1][8]; H[2] = P[2][8];
    H[3] = P[3][8]; H[4] = P[4][8]; H[5] = P[5][8];
    H[6] = P[6][8]; H[7] = P[7][8]; H[8] = 1.0;

    return true;
}

void OpenFIRE_Perspective::applyLensCorrection(double &x, double &y) {
  if (k1 == 0.0f) return; 
  double dx = x - (double)centerX; 
  double dy = y - (double)centerY;
  double r2 = (dx * dx) + (dy * dy);
  double distortion = 1.0 + ((double)k1 * r2);
  x = (double)centerX + (dx * distortion); 
  y = (double)centerY + (dy * distortion);
}

// ========================================================================================
// WARP PRINCIPALE
// ========================================================================================
void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  
  // 1. Calcolo Area sui dati GREZZI ma a 64-bit per evitare troncamenti.
  // Ordine perimetro: 0(TL) -> 1(TR) -> 3(BR) -> 2(BL)
  double d_x0 = x0, d_y0 = y0;
  double d_x1 = x1, d_y1 = y1;
  double d_x2 = x2, d_y2 = y2;
  double d_x3 = x3, d_y3 = y3;

  double currentArea = 0.5 * fabs(
      (d_x0*d_y1 - d_x1*d_y0) + 
      (d_x1*d_y3 - d_x3*d_y1) + 
      (d_x3*d_y2 - d_x2*d_y3) + 
      (d_x2*d_y0 - d_x0*d_y2)
  );

  // 2. NORMALIZZAZIONE: Riduciamo la scala per proteggere la stabilità di Gauss-Jordan
  double src[4][2] = {
      { d_x0 / NORM_SCALE, d_y0 / NORM_SCALE },
      { d_x1 / NORM_SCALE, d_y1 / NORM_SCALE },
      { d_x2 / NORM_SCALE, d_y2 / NORM_SCALE },
      { d_x3 / NORM_SCALE, d_y3 / NORM_SCALE }
  };

  double dst[4][2] = {
      { (double)dx0 / NORM_SCALE, (double)dy0 / NORM_SCALE },
      { (double)dx1 / NORM_SCALE, (double)dy1 / NORM_SCALE },
      { (double)dx2 / NORM_SCALE, (double)dy2 / NORM_SCALE },
      { (double)dx3 / NORM_SCALE, (double)dy3 / NORM_SCALE }
  };

  for(int i = 0; i < 4; i++) {
      applyLensCorrection(src[i][0], src[i][1]);
  }

  // 3. Gestione Inizializzazione Parallasse
  if (!init) {
      if (currentArea > 1.0) {
          baseArea = currentArea;
          smoothedArea = currentArea;
          init = true;
      }
      if (!init) return;
  }

  if (currentArea > 1.0) {
      smoothedArea = (0.1 * currentArea) + (0.9 * smoothedArea);
  }

  double dynamicSrcY = (double)srcY;
  if (parallaxFactor != 0.0f && baseArea > 1.0 && smoothedArea > 1.0) {
      double ratio = sqrt(baseArea / smoothedArea);
      dynamicSrcY += (double)parallaxFactor * (ratio - 1.0);
  }

  // 4. Motore DLT
  if (computeHomography(src, dst, warpmatrix)) {
      
      // Proiezione del punto scalato
      double nx = (double)srcX / NORM_SCALE;
      double ny = dynamicSrcY / NORM_SCALE;

      double w = (warpmatrix[6] * nx) + (warpmatrix[7] * ny) + warpmatrix[8];
      
      if (fabs(w) > 1e-12) {
          double r0 = (warpmatrix[0] * nx) + (warpmatrix[1] * ny) + warpmatrix[2];
          double r1 = (warpmatrix[3] * nx) + (warpmatrix[4] * ny) + warpmatrix[5];
          
          // De-Normalizzazione: riportiamo i dati alla scala originale
          dstX = (int)round((r0 / w) * NORM_SCALE);
          dstY = (int)round((r1 / w) * NORM_SCALE);
      }
  }
}

void OpenFIRE_Perspective::setLensCorrection(float coefficientK1) { k1 = coefficientK1; }
void OpenFIRE_Perspective::setDynamicParallax(float hardwareOffset) { parallaxFactor = hardwareOffset; }
void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) { srcX = adjustedX; srcY = adjustedY; }
void OpenFIRE_Perspective::deinit(bool set) { init = set; }
int OpenFIRE_Perspective::getX() { return dstX; }
int OpenFIRE_Perspective::getY() { return dstY; }

#endif // USE_PERSPECTIVE_ADVANCED
#endif


#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include <math.h>

void computeSquareToQuad(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {

  float dx1 = x1 - x2;
  float dy1 = y1 - y2;
  float dx2 = x3 - x2;
  float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3;
  float sy = y0 - y1 + y2 - y3;
  
  float det = (dx1 * dy2 - dx2 * dy1);
  if (fabsf(det) < 1e-6f) return; // Evita crash per divisione per zero

  float g = (sx * dy2 - dx2 * sy) / det;
  float h = (dx1 * sy - sx * dy1) / det;
  float a = x1 - x0 + g * x1;
  float b = x3 - x0 + h * x3;
  float c = x0;
  float d = y1 - y0 + g * y1;
  float e = y3 - y0 + h * y3;
  float f = y0;

  // Mappatura corretta per emulare il layout originale [X, Y, 1]
  mat[0] = a; mat[1] = d; mat[2] = g;
  mat[3] = b; mat[4] = e; mat[5] = h;
  mat[6] = c; mat[7] = f; mat[8] = 1.0f;
}

void computeQuadToSquare(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  
  computeSquareToQuad(mat, x0, y0, x1, y1, x2, y2, x3, y3);

  float a = mat[0], d = mat[1], g = mat[2];
  float b = mat[3], e = mat[4], h = mat[5];
  float c = mat[6], f = mat[7];

  float A =     e - f * h;
  float B = c * h - b;
  float C = b * f - c * e;
  float D = f * g - d;
  float E =     a - c * g;
  float F = c * d - a * f;
  float G = d * h - e * g;
  float H = b * g - a * h;
  float I = a * e - b * d;
  
  float det = a * A + b * D + c * G;
  if (fabsf(det) < 1e-6f) return;
  float idet = 1.0f / det;

  mat[0] = A * idet; mat[1] = D * idet; mat[2] = G * idet;
  mat[3] = B * idet; mat[4] = E * idet; mat[5] = H * idet;
  mat[6] = C * idet; mat[7] = F * idet; mat[8] = I * idet;
}

void multMats(const float* m1, const float* m2, float* res) {
  res[0] = m1[0]*m2[0] + m1[1]*m2[3] + m1[2]*m2[6];
  res[1] = m1[0]*m2[1] + m1[1]*m2[4] + m1[2]*m2[7];
  res[2] = m1[0]*m2[2] + m1[1]*m2[5] + m1[2]*m2[8];
  
  res[3] = m1[3]*m2[0] + m1[4]*m2[3] + m1[5]*m2[6];
  res[4] = m1[3]*m2[1] + m1[4]*m2[4] + m1[5]*m2[7];
  res[5] = m1[3]*m2[2] + m1[4]*m2[5] + m1[5]*m2[8];
  
  res[6] = m1[6]*m2[0] + m1[7]*m2[3] + m1[8]*m2[6];
  res[7] = m1[6]*m2[1] + m1[7]*m2[4] + m1[8]*m2[7];
  res[8] = m1[6]*m2[2] + m1[7]*m2[5] + m1[8]*m2[8];
}

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  if (!init) {
    computeSquareToQuad(dstmatrix, dx0, dy0, dx1, dy1, dx2, dy2, dx3, dy3);
    init = true;
  }
  computeQuadToSquare(srcmatrix, (float)x0, (float)y0, (float)x1, (float)y1, (float)x2, (float)y2, (float)x3, (float)y3);
  multMats(srcmatrix, dstmatrix, warpmatrix);

  // Calcolo finale che emula perfettamente [srcX, srcY, 1] * warpmatrix
  float r0 = (srcX * warpmatrix[0] + srcY * warpmatrix[3] + warpmatrix[6]);
  float r1 = (srcX * warpmatrix[1] + srcY * warpmatrix[4] + warpmatrix[7]);
  float r3 = (srcX * warpmatrix[2] + srcY * warpmatrix[5] + warpmatrix[8]);
  
  if (fabsf(r3) > 1e-6f) {
    dstX = (int)(r0 / r3);           // Mantenuto il cast a int originario
    dstY = (int)floorf(r1 / r3);     // Mantenuto il floorf originario
  }
}

void OpenFIRE_Perspective::source(float adjustedX, float adjustedY) {
  srcX = adjustedX;
  srcY = adjustedY;
}

void OpenFIRE_Perspective::deinit(bool set) {
  init = set;
}

int OpenFIRE_Perspective::getX() {
  return dstX;
}

int OpenFIRE_Perspective::getY() {
  return dstY;
}

#endif // USE_PERSPECTIVE_ADVANCED
#endif
///////////////////////////////////


#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#include "OpenFIRE_Perspective_Advanced.h"
#include "math.h"

void computeSquareToQuad(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {

  float dx1 = x1 - x2;
  float dy1 = y1 - y2;
  float dx2 = x3 - x2;
  float dy2 = y3 - y2;
  float sx = x0 - x1 + x2 - x3;
  float sy = y0 - y1 + y2 - y3;
  float g = (sx * dy2 - dx2 * sy) / (dx1 * dy2 - dx2 * dy1);
  float h = (dx1 * sy - sx * dy1) / (dx1 * dy2 - dx2 * dy1);
  float a = x1 - x0 + g * x1;
  float b = x3 - x0 + h * x3;
  float c = x0;
  float d = y1 - y0 + g * y1;
  float e = y3 - y0 + h * y3;
  float f = y0;

  mat[ 0] = a;
  mat[ 1] = d;
  mat[ 2] = 0.0f;
  mat[ 3] = g;
  mat[ 4] = b;
  mat[ 5] = e;
  mat[ 6] = 0.0f;
  mat[ 7] = h;
  mat[ 8] = 0.0f;
  mat[ 9] = 0.0f;
  mat[10] = 1.0f;
  mat[11] = 0.0f;
  mat[12] = c;
  mat[13] = f;
  mat[14] = 0.0f;
  mat[15] = 1.0f;
}

void computeQuadToSquare(float* mat, float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  
  computeSquareToQuad(mat, x0, y0, x1, y1, x2, y2, x3, y3);

  float a = mat[ 0];
  float d = mat[ 1];
  float g = mat[ 3];
  float b = mat[ 4];
  float e = mat[ 5];
  float h = mat[ 7];
  float c = mat[12];
  float f = mat[13];

  float A =     e - f * h;
  float B = c * h - b;
  float C = b * f - c * e;
  float D = f * g - d;
  float E =     a - c * g;
  float F = c * d - a * f;
  float G = d * h - e * g;
  float H = b * g - a * h;
  float I = a * e - b * d;
  float idet = 1.0 / (a * A + b * D + c * G);

  mat[ 0] = A * idet;
  mat[ 1] = D * idet;
  mat[ 2] = 0.0f;
  mat[ 3] = G * idet;

  mat[ 4] = B * idet;
  mat[ 5] = E * idet;
  mat[ 6] = 0.0f;
  mat[ 7] = H * idet;

  mat[ 8] = 0.0f;
  mat[ 9] = 0.0f;
  mat[10] = 1.0f;
  mat[11] = 0.0f;

  mat[12] = C * idet;
  mat[13] = F * idet;
  mat[14] = 0.0f;
  mat[15] = I * idet;
}

void multMats(float* a, float* b, float* res) {

  for (int r = 0; r < 4; r++) {
    int ri = r * 4;
    for (int c = 0; c < 4; c++) {
      res[ri + c] = (
        a[ri + 0] * b[c +  0] +
        a[ri + 1] * b[c +  4] +
        a[ri + 2] * b[c +  8] +
        a[ri + 3] * b[c + 12]);
    }
  }
}

void OpenFIRE_Perspective::warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3) {
  if (!init) {
   // float dx0 = 567;
   // float dx1 = 1360;
   // float dy2 = 1080;
    computeSquareToQuad(dstmatrix, dx0, dy0, dx1, dy1, dx2, dy2, dx3, dy3);
    init = true;
  }
  computeQuadToSquare(srcmatrix, float(x0), float(y0), float(x1), float(y1), float(x2), float(y2), float(x3), float(y3));
  multMats(srcmatrix, dstmatrix, warpmatrix);
  float r0 = (srcX * warpmatrix[0] + srcY * warpmatrix[4] + warpmatrix[12]);
  float r1 = (srcX * warpmatrix[1] + srcY * warpmatrix[5] + warpmatrix[13]);
  float r3 = (srcX * warpmatrix[3] + srcY * warpmatrix[7] + warpmatrix[15]);
  dstX = (r0/r3);
  dstY = floorf(r1/r3);
}

void OpenFIRE_Perspective:: source(float adjustedX, float adjustedY) {
  srcX = adjustedX;
  srcY = adjustedY;
}

void OpenFIRE_Perspective::deinit (bool set) {
  init = set;
}

int OpenFIRE_Perspective::getX() {
  return dstX;
}

int OpenFIRE_Perspective::getY() {
  return dstY;
}


#endif // USE_PERSPECTIVE_ADVANCED

#endif // COMMENTO