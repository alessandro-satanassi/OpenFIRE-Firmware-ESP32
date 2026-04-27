
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

