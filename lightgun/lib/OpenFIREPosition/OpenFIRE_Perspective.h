#ifndef USE_PERSPECTIVE_ADVANCED
/*
 * @file OpenFIRE_Perspective.h
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

#ifndef OpenFIRE_Perspective_h
#define OpenFIRE_Perspective_h

#include <OpenFIRECameraProfile.h>

class OpenFIRE_Perspective {

private:


  bool init = false;
  float srcmatrix[16];
  float dstmatrix[16];
  float warpmatrix[16];

  float dx0;
  float dy0;
  float dx1;
  float dy1; 
  float dx2; 
  float dy2; 
  float dx3; 
  float dy3;

  /*
  float srcX = 512.0f;
  float srcY = 384.0f;
  */
  // Centro dinamico della telecamera fisica (FOV Center)
  float srcX = 0.0f;
  float srcY = 0.0f;
  /*
  // Centro dinamico basato sullo spazio virtuale unificato (già normalizzato)
  float srcX = (float)MouseResX * 0.5f;
  float srcY = (float)MouseResY * 0.5f;
  */

  int dstX;
  int dstY;

public:
  void configure(const CameraProfile& profile) { srcX = (float)profile.camResX * 0.5f; srcY = (float)profile.camResY * 0.5f; }
  void warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3);
  void source(float adjustedX, float adjustedY);
  void deinit (bool set);
  int getX();
  int getY();
};

#endif

#endif // USE_PERSPECTIVE_ADVANCED