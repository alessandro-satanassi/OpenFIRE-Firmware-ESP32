#ifdef USE_PERSPECTIVE_ADVANCED
/*!
 * @file OpenFIRE_Perspective_Advanced.h
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


#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

#include "OpenFIREConst.h" 

// ==============================================================================
// PARAMETRI MATEMATICI E FISICI
// ==============================================================================
// Regola l'altezza della mira (asse Y) in base alla distanza del giocatore dalla TV.
// Serve a compensare la differenza fisica di altezza tra la canna della pistola
// e il sensore ottico posizionato al suo interno (o sopra).
static constexpr float DEFAULT_PARALLAX_FACTOR = 0.0f;

// Scala numerica interna: non rappresenta né la risoluzione della camera né quella
// del monitor. Mantenerla fissa preserva esattamente il comportamento numerico
// già validato con DFRobot e limita le perdite di precisione nelle omografie.
static constexpr float NORM_SCALE = 10000.0f;
static constexpr float INV_NORM_SCALE = 1.0f / NORM_SCALE;

// Centro geometrico del sensore nello spazio Mouse.
static constexpr float CX = (float)MouseResX * 0.5f;
static constexpr float CY = (float)MouseResY * 0.5f;
static constexpr float INV_CX = 1.0f / CX;
static constexpr float INV_CY = 1.0f / CY;

// Soglia minima dell'area espressa come frazione dell'intero spazio Mouse.
// Il coefficiente è scelto in modo che con la DFRobot storica (4096x3072 Mouse)
// MIN_QUAD_AREA sia esattamente 10.0f, mantenendo il comportamento originale.
static constexpr float MIN_QUAD_AREA_RATIO = 7.947286e-7f;
static constexpr float MIN_QUAD_AREA = (float)MouseResX * (float)MouseResY * MIN_QUAD_AREA_RATIO;

static_assert(MouseResX > 0 && MouseResY > 0, "MouseResX/MouseResY must be greater than zero");
// ==============================================================================

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  // Matrici lineari (1D) da 9 elementi invece che [3][3]. 
  // Migliora drasticamente il caching della CPU e permette cicli di accesso 
  // diretti in memoria contigua durante le pesanti moltiplicazioni matriciali.
  float srcmatrix[9];
  float dstmatrix[9];
  float warpmatrix[9];

  float srcX = CX;
  float srcY = CY;

  // Coefficienti di distorsione radiale definiti dalla camera in OpenFIREConst.h.
  float k1 = CamLensRadialK1;
  float k2 = CamLensRadialK2;
  float parallaxFactor = DEFAULT_PARALLAX_FACTOR;
  
  // Tracciamento dell'area per il parallasse. L'area calcolata viene usata come proxy 
  // affidabile e leggero (al posto della trigonometria) per stimare i cambiamenti 
  // di distanza (asse Z) del giocatore rispetto allo schermo.
  float baseArea = 0.0f;        
  float smoothedArea = 0.0f; 

  int dstX = 0;
  int dstY = 0;

  inline void applyLensCorrection(float &x, float &y);
  inline float calculateQuadArea(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3);

public:
  void warp(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3);
  void source(float adjustedX, float adjustedY);
  
  // Mantiene la vecchia API K1-only e aggiunge l'overload K1+K2.
  void setLensCorrection(float coefficientK1);
  void setLensCorrection(float coefficientK1, float coefficientK2);
  void setDynamicParallax(float hardwareOffset); 
  
  void deinit(bool set);
  int getX();
  int getY();
};

#endif

#endif // USE_PERSPECTIVE_ADVANCED
