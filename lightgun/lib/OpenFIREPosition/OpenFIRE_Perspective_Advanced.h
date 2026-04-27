
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

#include "OpenFIREConst.h" 

// ==============================================================================
// PARAMETRI FISICI DELLA LIGHTGUN
// ==============================================================================
#define DEFAULT_LENS_K1 0.006f       
#define DEFAULT_PARALLAX_FACTOR 0.0f  

// Costanti matematiche assolute (Calcolate dal compilatore, costo zero su ESP32)
static constexpr float NORM_SCALE = 10000.0f;  
static constexpr float INV_NORM_SCALE = 1.0f / 10000.0f; 

static constexpr float CX = (float)MouseResX * 0.5f;
static constexpr float CY = (float)MouseResY * 0.5f;
static constexpr float INV_CX = 1.0f / CX;
static constexpr float INV_CY = 1.0f / CY;
// ==============================================================================

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  float srcmatrix[9];
  float dstmatrix[9];
  float warpmatrix[9];

  float srcX = CX;
  float srcY = CY;

  float k1 = DEFAULT_LENS_K1;              
  float parallaxFactor = DEFAULT_PARALLAX_FACTOR;  
  float baseArea = 0.0f;        
  float smoothedArea = 0.0f; 

  int dstX = 0;
  int dstY = 0;

  inline void applyLensCorrection(float &x, float &y);
  inline float calculateQuadArea(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3);

public:
  void warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3);
  void source(float adjustedX, float adjustedY);
  
  void setLensCorrection(float coefficientK1);
  void setDynamicParallax(float hardwareOffset); 
  
  void deinit(bool set);
  int getX();
  int getY();
};

#endif

#endif // USE_PERSPECTIVE_ADVANCED
