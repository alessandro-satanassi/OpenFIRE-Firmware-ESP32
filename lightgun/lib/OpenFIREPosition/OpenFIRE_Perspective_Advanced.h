
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


#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

#include "OpenFIREConst.h" 

// ==============================================================================
// PARAMETRI FISICI DELLA LIGHTGUN
// ==============================================================================
#define DEFAULT_LENS_K1 0.006f       
#define DEFAULT_PARALLAX_FACTOR 150.0f //0.0f // 150.0f per canna 3 cm sotto congne di mira  

// Costanti pre-calcolate (Zero costo per la CPU in tempo reale)
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
#endif

#endif // USE_PERSPECTIVE_ADVANCED

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

#include "OpenFIREConst.h" 

// ==============================================================================
// PARAMETRI FISICI DELLA LIGHTGUN
// ==============================================================================
#define DEFAULT_LENS_K1 0.006f       
#define DEFAULT_PARALLAX_FACTOR 0.0f  

// Costanti pre-calcolate a tempo di compilazione (Zero costo per la CPU)
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

  // Funzioni inline per eliminare il costo di chiamata (Call Overhead)
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
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

#include "OpenFIREConst.h" 

// ==============================================================================
// PARAMETRI FISICI DELLA LIGHTGUN
// ==============================================================================
#define DEFAULT_LENS_K1 0.006f       
#define DEFAULT_PARALLAX_FACTOR 0.0f  

// Il vero segreto della stabilità sui float a 32-bit
#define NORM_SCALE 10000.0f  
// ==============================================================================

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  // Tornati a float (32-bit) per sfruttare la FPU Hardware dell'ESP32
  float srcmatrix[9];
  float dstmatrix[9];
  float warpmatrix[9];

  float srcX = (float)MouseResX / 2.0f;
  float srcY = (float)MouseResY / 2.0f;

  float k1 = DEFAULT_LENS_K1;              
  float parallaxFactor = DEFAULT_PARALLAX_FACTOR;  
  float baseArea = 0.0f;        
  float smoothedArea = 0.0f; 

  const float centerX = (float)MouseResX / 2.0f;
  const float centerY = (float)MouseResY / 2.0f;

  int dstX = 0;
  int dstY = 0;

  void applyLensCorrection(float &x, float &y);
  float calculateQuadArea(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3);

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
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

#include "OpenFIREConst.h" 

// ==============================================================================
// PARAMETRI FISICI DELLA LIGHTGUN
// ==============================================================================
// K1 è una percentuale relativa. 0.006 = espansione dello 0.6% ai bordi
#define DEFAULT_LENS_K1 0.006f       
#define DEFAULT_PARALLAX_FACTOR 0.0f  
// ==============================================================================

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  // Matrici a 64-bit (double) per precisione pura
  double srcmatrix[9];
  double dstmatrix[9];
  double warpmatrix[9];

  double srcX = (double)MouseResX / 2.0;
  double srcY = (double)MouseResY / 2.0;

  float k1 = DEFAULT_LENS_K1;              
  float parallaxFactor = DEFAULT_PARALLAX_FACTOR;  
  double baseArea = 0.0;        
  double smoothedArea = 0.0; 

  // Centro ottico dinamico 
  const double centerX = (double)MouseResX / 2.0;
  const double centerY = (double)MouseResY / 2.0;

  int dstX = 0;
  int dstY = 0;

  void applyLensCorrection(double &x, double &y);
  double calculateQuadArea(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);

public:
  // API Esterna: x0(TL), x1(TR), x2(BL), x3(BR)
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
#endif


#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

// Assicurati di includere il tuo file con le costanti
#include "OpenFIREConst.h" 

// ==============================================================================
// PARAMETRI FISICI DELLA LIGHTGUN
// ==============================================================================
// Ora è una percentuale! 0.006 = 0.6% di spinta verso l'esterno ai bordi
#define DEFAULT_LENS_K1 0.006f       
#define DEFAULT_PARALLAX_FACTOR 0.0f  
// ==============================================================================

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  // Matrici a 64-bit (double) per precisione infinita
  double srcmatrix[9];
  double dstmatrix[9];
  double warpmatrix[9];

  // Il centro del mirino parte dal centro della telecamera "scalata"
  double srcX = (double)MouseResX / 2.0;
  double srcY = (double)MouseResY / 2.0;

  float k1 = DEFAULT_LENS_K1;              
  float parallaxFactor = DEFAULT_PARALLAX_FACTOR;  
  double baseArea = 0.0;        
  double smoothedArea = 0.0; 

  // Centro ottico dinamico basato sulle tue costanti
  const double centerX = (double)MouseResX / 2.0;
  const double centerY = (double)MouseResY / 2.0;

  int dstX = 0;
  int dstY = 0;

  void applyLensCorrection(double &x, double &y);
  double calculateQuadArea(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);

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
#endif

#endif // USE_PERSPECTIVE_ADVANCED

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

// ==============================================================================
// PARAMETRI FISICI DELLA LIGHTGUN
// ==============================================================================
// Valore ideale per la DFRobot SEN0158 (PixArt). Corregge il barilotto.
#define DEFAULT_LENS_K1 1.5e-8f       
#define DEFAULT_PARALLAX_FACTOR 0.0f  
// ==============================================================================

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  // Matrici elevate a 64-bit (double) per precisione infinita e stabilità totale
  double srcmatrix[9];
  double dstmatrix[9];
  double warpmatrix[9];

  double srcX = 512.0;
  double srcY = 384.0;

  float k1 = DEFAULT_LENS_K1;              
  float parallaxFactor = DEFAULT_PARALLAX_FACTOR;  
  double baseArea = 0.0;        
  double smoothedArea = 0.0; 

  // Centro geometrico del sensore IR (tipicamente 1024x768)
  const double centerX = 512.0;
  const double centerY = 384.0;

  int dstX = 0;
  int dstY = 0;

  void applyLensCorrection(double &x, double &y);
  double calculateQuadArea(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);

public:
  // L'API nativa (x0=TL, x1=TR, x2=BL, x3=BR) 
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
#endif


#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

// ==============================================================================
// PARAMETRI FISICI DELLA LIGHTGUN
// ==============================================================================
#define DEFAULT_LENS_K1 0.0f          
#define DEFAULT_PARALLAX_FACTOR 0.0f  
// ==============================================================================

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  // Matrici elevate a 64-bit (double) per precisione infinita
  double srcmatrix[9];
  double dstmatrix[9];
  double warpmatrix[9];

  double srcX = 512.0;
  double srcY = 384.0;

  float k1 = DEFAULT_LENS_K1;              
  float parallaxFactor = DEFAULT_PARALLAX_FACTOR;  
  double baseArea = 0.0;        
  double smoothedArea = 0.0; 

  const double centerX = 512.0;
  const double centerY = 384.0;

  int dstX = 0;
  int dstY = 0;

  void applyLensCorrection(double &x, double &y);
  double calculateQuadArea(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);

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
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

// ==============================================================================
// PARAMETRI FISICI DELLA LIGHTGUN
// ==============================================================================
#define DEFAULT_LENS_K1 0.0f          // Es: 5.0e-8f per raddrizzare il barilotto
#define DEFAULT_PARALLAX_FACTOR 0.0f  // Es: 15.0f per correggere l'offset cam/canna
// ==============================================================================

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  // Matrici Row-Major (Identiche all'algoritmo originale)
  float srcmatrix[9];
  float dstmatrix[9];
  float warpmatrix[9];

  float srcX = 512.0f;
  float srcY = 384.0f;

  // Variabili per le correzioni fisiche
  float k1 = DEFAULT_LENS_K1;              
  float parallaxFactor = DEFAULT_PARALLAX_FACTOR;  
  float baseArea = 0.0f;        
  float smoothedArea = 0.0f; 

  const float centerX = 512.0f;
  const float centerY = 384.0f;

  int dstX = 0;
  int dstY = 0;

  // Helpers per la fisica
  void applyLensCorrection(float &x, float &y);
  float calculateQuadArea(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3);

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
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  // Array a 9 elementi ottimizzato, layout Row-Major fedele al Wiimote Warper originale
  float srcmatrix[9];
  float dstmatrix[9];
  float warpmatrix[9];

  float srcX = 512.0f;
  float srcY = 384.0f;

  int dstX = 0;
  int dstY = 0;

public:
  // L'API accetta i punti nativamente nell'ordine TL, TR, BL, BR
  void warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3);
  void source(float adjustedX, float adjustedY);
  void deinit(bool set);
  int getX();
  int getY();
};

#endif

#endif // USE_PERSPECTIVE_ADVANCED
#endif


#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

// ==============================================================================
// PARAMETRI FISICI DELLA LIGHTGUN
// ==============================================================================
#define DEFAULT_LENS_K1 0.0f          
#define DEFAULT_PARALLAX_FACTOR 0.0f  
// Scala di normalizzazione per proteggere la FPU (copre coordinate fino a +-100000)
#define NORM_SCALE 10000.0  
// ==============================================================================

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  // La matrice ora è a 64-bit (double) per non perdere MAI precisione
  double warpmatrix[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  float srcX = 512.0f;
  float srcY = 384.0f;
  
  float k1 = DEFAULT_LENS_K1;              
  float parallaxFactor = DEFAULT_PARALLAX_FACTOR;  
  
  // Anche le aree usano i 64-bit per gestire calcoli oltre i 100 milioni
  double baseArea = 0.0;        
  double smoothedArea = 0.0;  
  
  const float centerX = 512.0f;
  const float centerY = 384.0f;

  int dstX = 0;
  int dstY = 0;

  void applyLensCorrection(double &x, double &y);
  bool computeHomography(double src[4][2], double dst[4][2], double H[9]);

public:
  void warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3);
  void source(float adjustedX, float adjustedY);
  
  void setLensCorrection(float coefficientK1);
  void setDynamicParallax(float hardwareOffset); 
  
  void deinit (bool set);
  int getX();
  int getY();
};

#endif

#endif // USE_PERSPECTIVE_ADVANCED
#endif

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

class OpenFIRE_Perspective {

private:
  bool init = false;
  
  // Matrici compatte 3x3
  float srcmatrix[9];
  float dstmatrix[9];
  float warpmatrix[9];

  float srcX = 512.0f;
  float srcY = 384.0f;

  int dstX = 0;
  int dstY = 0;

public:
  void warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3);
  void source(float adjustedX, float adjustedY);
  void deinit (bool set);
  int getX();
  int getY();
};

#endif

#endif // USE_PERSPECTIVE_ADVANCED
#endif

/////////////////////////////////////////////

#ifdef COMMENTO
#ifdef USE_PERSPECTIVE_ADVANCED

#ifndef OpenFIRE_Perspective_Advanced_h
#define OpenFIRE_Perspective_Advanced_h

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

  float srcX = 512.0f;
  float srcY = 384.0f;

  int dstX;
  int dstY;

public:
  void warp(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, float dx0, float dy0, float dx1, float dy1, float dx2, float dy2, float dx3, float dy3);
  void source(float adjustedX, float adjustedY);
  void deinit (bool set);
  int getX();
  int getY();
};

#endif


#endif // USE_PERSPECTIVE_ADVANCED
# endif // COMMENTO