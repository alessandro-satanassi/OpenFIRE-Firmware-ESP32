/*!
 * @file OpenFIRE_Scope_Decoder.cpp
 * @brief Decoder of the video slices and roll rotation of the scope
 * @n CPP Decoder of the video slices and roll rotation of the scope
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

// EN: Techniques of the prototype dec5.c: 32-bit bit reader refilled 16 bits at a time, fast 9-bit AC table,
//     dequantization multiplied by the AAN factors, integer AAN IDCT with reduced variants, colour through tables
//     with clamp and byte swap included. The hot code is in IRAM, the tables in internal RAM.
// IT: Tecniche del prototipo dec5.c: lettore di bit a 32 bit ricaricato 16 bit alla volta, tabella AC veloce a 9 bit,
//     dequantizzazione moltiplicata per i fattori AAN, IDCT AAN intera con varianti ridotte, colore con tabelle
//     che includono clamp e scambio dei byte. Il codice caldo sta in IRAM, le tabelle in RAM interna.

#include "OpenFIRE_Scope_Decoder.h"
#include <string.h>
#include <math.h>

#if defined(ARDUINO_ARCH_ESP32)
  #include "esp_attr.h"
#endif
#ifndef IRAM_ATTR
  #define IRAM_ATTR
#endif
#ifndef DRAM_ATTR
  #define DRAM_ATTR
#endif
// EN: forced inline of the hot helpers; MSVC (scope daemon on Windows) has no __attribute__.
// IT: inline forzato delle funzioni calde; MSVC (demone dei mirini su Windows) non ha __attribute__.
#if defined(_MSC_VER)
  #define SCOPE_INLINE static __forceinline
#else
  #define SCOPE_INLINE static inline __attribute__((always_inline))
#endif

#define W  SCOPE_DEC_W
#define H  SCOPE_DEC_H

// ===================================================================================
// STANDARD TABLES / TABELLE STANDARD (ITU-T T.81)
// ===================================================================================
static DRAM_ATTR const uint8_t ZZ[64] = {
   0,  1,  8, 16,  9,  2,  3, 10, 17, 24, 32, 25, 18, 11,  4,  5, 12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13,  6,  7, 14, 21, 28,
  35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51, 58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63 };

// EN: K.1 luminance and K.2 chrominance, natural order (row by row). / IT: K.1 luminanza e K.2 crominanza, ordine naturale (riga per riga).
static const uint8_t Q_LUM[64] = {
  16, 11, 10, 16, 24, 40, 51, 61,  12, 12, 14, 19, 26, 58, 60, 55,  14, 13, 16, 24, 40, 57, 69, 56,  14, 17, 22, 29, 51, 87, 80, 62,
  18, 22, 37, 56, 68,109,103, 77,  24, 35, 55, 64, 81,104,113, 92,  49, 64, 78, 87,103,121,120,101,  72, 92, 95, 98,112,100,103, 99 };
static const uint8_t Q_CHR[64] = {
  17, 18, 24, 47, 99, 99, 99, 99,  18, 21, 26, 66, 99, 99, 99, 99,  24, 26, 56, 99, 99, 99, 99, 99,  47, 66, 99, 99, 99, 99, 99, 99,
  99, 99, 99, 99, 99, 99, 99, 99,  99, 99, 99, 99, 99, 99, 99, 99,  99, 99, 99, 99, 99, 99, 99, 99,  99, 99, 99, 99, 99, 99, 99, 99 };

// EN: K.3-K.6: number of codes per length (1-16) and values. / IT: K.3-K.6: numero di codici per lunghezza (1-16) e valori.
static const uint8_t DC_L_BITS[16] = { 0, 1, 5, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 };
static const uint8_t DC_C_BITS[16] = { 0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 };
static DRAM_ATTR const uint8_t DC_VALS[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
static const uint8_t AC_L_BITS[16] = { 0, 2, 1, 3, 3, 2, 4, 3, 5, 5, 4, 4, 0, 0, 1, 0x7d };
static DRAM_ATTR const uint8_t AC_L_VALS[162] = {
  0x01,0x02,0x03,0x00,0x04,0x11,0x05,0x12,0x21,0x31,0x41,0x06,0x13,0x51,0x61,0x07,0x22,0x71,0x14,0x32,0x81,0x91,0xa1,0x08,
  0x23,0x42,0xb1,0xc1,0x15,0x52,0xd1,0xf0,0x24,0x33,0x62,0x72,0x82,0x09,0x0a,0x16,0x17,0x18,0x19,0x1a,0x25,0x26,0x27,0x28,
  0x29,0x2a,0x34,0x35,0x36,0x37,0x38,0x39,0x3a,0x43,0x44,0x45,0x46,0x47,0x48,0x49,0x4a,0x53,0x54,0x55,0x56,0x57,0x58,0x59,
  0x5a,0x63,0x64,0x65,0x66,0x67,0x68,0x69,0x6a,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7a,0x83,0x84,0x85,0x86,0x87,0x88,0x89,
  0x8a,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9a,0xa2,0xa3,0xa4,0xa5,0xa6,0xa7,0xa8,0xa9,0xaa,0xb2,0xb3,0xb4,0xb5,0xb6,
  0xb7,0xb8,0xb9,0xba,0xc2,0xc3,0xc4,0xc5,0xc6,0xc7,0xc8,0xc9,0xca,0xd2,0xd3,0xd4,0xd5,0xd6,0xd7,0xd8,0xd9,0xda,0xe1,0xe2,
  0xe3,0xe4,0xe5,0xe6,0xe7,0xe8,0xe9,0xea,0xf1,0xf2,0xf3,0xf4,0xf5,0xf6,0xf7,0xf8,0xf9,0xfa };
static const uint8_t AC_C_BITS[16] = { 0, 2, 1, 2, 4, 4, 3, 4, 7, 5, 4, 4, 0, 1, 2, 0x77 };
static DRAM_ATTR const uint8_t AC_C_VALS[162] = {
  0x00,0x01,0x02,0x03,0x11,0x04,0x05,0x21,0x31,0x06,0x12,0x41,0x51,0x07,0x61,0x71,0x13,0x22,0x32,0x81,0x08,0x14,0x42,0x91,
  0xa1,0xb1,0xc1,0x09,0x23,0x33,0x52,0xf0,0x15,0x62,0x72,0xd1,0x0a,0x16,0x24,0x34,0xe1,0x25,0xf1,0x17,0x18,0x19,0x1a,0x26,
  0x27,0x28,0x29,0x2a,0x35,0x36,0x37,0x38,0x39,0x3a,0x43,0x44,0x45,0x46,0x47,0x48,0x49,0x4a,0x53,0x54,0x55,0x56,0x57,0x58,
  0x59,0x5a,0x63,0x64,0x65,0x66,0x67,0x68,0x69,0x6a,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7a,0x82,0x83,0x84,0x85,0x86,0x87,
  0x88,0x89,0x8a,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9a,0xa2,0xa3,0xa4,0xa5,0xa6,0xa7,0xa8,0xa9,0xaa,0xb2,0xb3,0xb4,
  0xb5,0xb6,0xb7,0xb8,0xb9,0xba,0xc2,0xc3,0xc4,0xc5,0xc6,0xc7,0xc8,0xc9,0xca,0xd2,0xd3,0xd4,0xd5,0xd6,0xd7,0xd8,0xd9,0xda,
  0xe2,0xe3,0xe4,0xe5,0xe6,0xe7,0xe8,0xe9,0xea,0xf2,0xf3,0xf4,0xf5,0xf6,0xf7,0xf8,0xf9,0xfa };

// EN: AAN factors (IJG jidctfst) scaled by 2^14, natural order. / IT: fattori AAN (IJG jidctfst) scalati di 2^14, ordine naturale.
static const int16_t AAN[64] = {
  16384, 22725, 21407, 19266, 16384, 12873,  8867,  4520,  22725, 31521, 29692, 26722, 22725, 17855, 12299,  6270,
  21407, 29692, 27969, 25172, 21407, 16819, 11585,  5906,  19266, 26722, 25172, 22654, 19266, 15137, 10426,  5315,
  16384, 22725, 21407, 19266, 16384, 12873,  8867,  4520,  12873, 17855, 16819, 15137, 12873, 10114,  6967,  3552,
   8867, 12299, 11585, 10426,  8867,  6967,  4799,  2446,   4520,  6270,  5906,  5315,  4520,  3552,  2446,  1247 };

// ===================================================================================
// TABLES BUILT AT POWER ON / TABELLE COSTRUITE ALL'ACCENSIONE
// ===================================================================================
#define LUTB 9   // fast tables: 9 bits / tabelle veloci: 9 bit

typedef struct {
  uint16_t lut[1 << LUTB];   // (length << 8) | symbol, 0 = code longer than 9 bits / 0 = codice piu' lungo di 9 bit
  int32_t  maxcode[18];
  int32_t  valoff[18];
  const uint8_t *vals;
} ScopeHuff;

static ScopeHuff HDL, HDC, HAL, HAC;
static int32_t FACL[1 << LUTB], FACC[1 << LUTB];   // fast AC: value << 8 | run << 4 | bits / AC veloce: valore << 8 | run << 4 | bit
static int32_t QSL[64], QSC[64];                   // quantization x AAN / quantizzazione x AAN
static int     g_q = -1;                           // q of the tables above / q delle tabelle qui sopra
static uint8_t  CLAMP_TAB[1024];
static uint16_t R5[1024], G6[1024], B5[1024];      // with clamp and byte swap / con clamp e scambio dei byte
static int16_t  CRR16[256], CBB16[256];
static int32_t  CRG[256], CBG[256];

// EN: Circle: visible pixels (radius 120, for the rotation) and pixels written by the decoder (radius 121 inside the
//     coded MCU: the nearest-neighbour rounding of the rotation can read up to 0.71 pixels beyond the visible circle).
// IT: Cerchio: pixel visibili (raggio 120, per la rotazione) e pixel scritti dal decoder (raggio 121 dentro gli MCU
//     codificati: l'arrotondamento nearest-neighbour della rotazione puo' leggere fino a 0,71 pixel oltre il cerchio visibile).
static int16_t VXS[H], VXE[H];   // visible span of row y / tratto visibile della riga y
static int16_t WXS[H], WXE[H];   // written span of row y / tratto scritto della riga y
static uint8_t MCU_OK[SCOPE_DEC_SLICES][15];
static uint8_t MCU_FULL[SCOPE_DEC_SLICES][15];
static uint16_t OUT_X[SCOPE_DEC_SLICES], OUT_W[SCOPE_DEC_SLICES];
// EN: YCbCr frame: number of coded MCU, MCU for the MCU not coded (black), window of every output band of 8 rows.
// IT: Frame YCbCr: numero di MCU codificati, MCU per gli MCU non codificati (nero), finestra di ogni fascia di uscita da 8 righe.
static uint16_t MCU_N = 0;
static uint8_t  BLACK_MCU[SCOPE_DEC_MCU_BYTES];
static uint16_t BAND_X[SCOPE_ROT_BANDS], BAND_W[SCOPE_ROT_BANDS];

static void huffBuild(ScopeHuff *h, const uint8_t *bits, const uint8_t *vals) {
  int code = 0, k = 0;
  memset(h->lut, 0, sizeof(h->lut));
  h->vals = vals;
  for (int l = 1; l <= 16; l++) {
    h->valoff[l] = k - code;
    for (int i = 0; i < bits[l - 1]; i++, code++, k++) {
      if (l <= LUTB) {
        const int sh = LUTB - l;
        for (int j = 0; j < (1 << sh); j++) h->lut[(code << sh) | j] = (uint16_t)((l << 8) | vals[k]);
      }
    }
    h->maxcode[l] = bits[l - 1] ? code - 1 : -1;
    code <<= 1;
  }
  h->valoff[17] = 0;
  h->maxcode[17] = 0x7fffffff;   // stops the search: 17 = invalid code / ferma la ricerca: 17 = codice non valido
}

static void facBuild(int32_t *fac, const ScopeHuff *h) {
  for (int i = 0; i < (1 << LUTB); i++) {
    fac[i] = 0;
    const uint32_t e = h->lut[i];
    const int len = e >> 8, sym = e & 0xFF, run = sym >> 4, size = sym & 15;
    if (!len || !size || len + size > LUTB) continue;
    const int bits = (i >> (LUTB - len - size)) & ((1 << size) - 1);
    const int v = bits < (1 << (size - 1)) ? bits - (1 << size) + 1 : bits;
    fac[i] = v * 256 + run * 16 + (len + size);
  }
}

static inline bool inCircle(int x, int y, double r) {
  const double dx = x - 119.5, dy = y - 119.5;
  return dx * dx + dy * dy <= r * r;
}

void scopeDecInit() {
  for (int i = 0; i < 1024; i++) { const int v = i - 384; CLAMP_TAB[i] = v < 0 ? 0 : v > 255 ? 255 : v; }
  for (int i = 0; i < 1024; i++) {
    const int v = CLAMP_TAB[i];
    const uint16_t r = (uint16_t)((v >> 3) << 11), g = (uint16_t)((v >> 2) << 5), b = (uint16_t)(v >> 3);
    R5[i] = (uint16_t)((r >> 8) | (r << 8)); G6[i] = (uint16_t)((g >> 8) | (g << 8)); B5[i] = (uint16_t)((b >> 8) | (b << 8));
  }
  for (int i = 0; i < 256; i++) {
    const int x = i - 128;
    CRR16[i] = (int16_t)((int16_t)(1.402 * x + (x >= 0 ? 0.5 : -0.5)) + 384);   // +384: clamp offset / offset del clamp
    CBB16[i] = (int16_t)((int16_t)(1.772 * x + (x >= 0 ? 0.5 : -0.5)) + 384);
    CRG[i] = (int32_t)(-0.714136 * 65536 * x);
    CBG[i] = (int32_t)(-0.344136 * 65536 * x) + 32768 + (384 << 16);
  }
  huffBuild(&HDL, DC_L_BITS, DC_VALS);   huffBuild(&HDC, DC_C_BITS, DC_VALS);
  huffBuild(&HAL, AC_L_BITS, AC_L_VALS); huffBuild(&HAC, AC_C_BITS, AC_C_VALS);
  facBuild(FACL, &HAL); facBuild(FACC, &HAC);
  g_q = -1;

  // EN: Coded MCU: at least one pixel within radius 120.5 (same table as the encoder).
  // IT: MCU codificato: almeno un pixel entro il raggio 120,5 (stessa tabella dell'encoder).
  for (int my = 0; my < SCOPE_DEC_SLICES; my++) for (int mx = 0; mx < 15; mx++) {
    uint8_t ok = 0;
    for (int y = my * 16; y < my * 16 + 16 && !ok; y++) for (int x = mx * 16; x < mx * 16 + 16; x++) if (inCircle(x, y, 120.5)) { ok = 1; break; }
    MCU_OK[my][mx] = ok;
  }
  for (int y = 0; y < H; y++) {
    int vs = W, ve = -1, ws = W, we = -1;
    for (int x = 0; x < W; x++) {
      if (inCircle(x, y, 120.0)) { if (x < vs) vs = x; ve = x; }
      if (inCircle(x, y, 121.0) && MCU_OK[y >> 4][x >> 4]) { if (x < ws) ws = x; we = x; }
    }
    VXS[y] = (int16_t)vs; VXE[y] = (int16_t)ve; WXS[y] = (int16_t)ws; WXE[y] = (int16_t)we;
  }
  for (int my = 0; my < SCOPE_DEC_SLICES; my++) {
    for (int mx = 0; mx < 15; mx++) {
      uint8_t full = MCU_OK[my][mx];
      for (int y = my * 16; y < my * 16 + 16; y++) if (WXS[y] > mx * 16 || WXE[y] < mx * 16 + 15) full = 0;
      MCU_FULL[my][mx] = full;
    }
    int x0 = W, x1 = -1;
    for (int y = my * 16; y < my * 16 + 16; y++) { if (VXS[y] < x0) x0 = VXS[y]; if (VXE[y] > x1) x1 = VXE[y]; }
    OUT_X[my] = (uint16_t)x0; OUT_W[my] = (uint16_t)(x1 - x0 + 1);
  }
  MCU_N = 0;
  for (int my = 0; my < SCOPE_DEC_SLICES; my++) for (int mx = 0; mx < 15; mx++) if (MCU_OK[my][mx]) MCU_N++;
  memset(BLACK_MCU, 0, 256);
  memset(BLACK_MCU + 256, 128, 128);
  for (int b = 0; b < SCOPE_ROT_BANDS; b++) {
    int x0 = W, x1 = -1;
    for (int y = b * SCOPE_ROT_BAND_H; y < (b + 1) * SCOPE_ROT_BAND_H; y++) { if (VXS[y] < x0) x0 = VXS[y]; if (VXE[y] > x1) x1 = VXE[y]; }
    BAND_X[b] = (uint16_t)x0; BAND_W[b] = (uint16_t)(x1 - x0 + 1);
  }
}

// EN: Quantization of the slice: IJG scale, then multiplied by the AAN factors. Only when q changes.
// IT: Quantizzazione della fascia: scala IJG, poi moltiplicata per i fattori AAN. Solo quando cambia q.
static void setQuant(int q) {
  if (q == g_q) return;
  g_q = q;
  const int scale = q < 50 ? 5000 / q : 200 - 2 * q;
  for (int i = 0; i < 64; i++) {
    int l = (Q_LUM[i] * scale + 50) / 100, c = (Q_CHR[i] * scale + 50) / 100;
    l = l < 1 ? 1 : l > 255 ? 255 : l;
    c = c < 1 ? 1 : c > 255 ? 255 : c;
    QSL[i] = (int32_t)(((int64_t)l * AAN[i] + (1 << 11)) >> 12);
    QSC[i] = (int32_t)(((int64_t)c * AAN[i] + (1 << 11)) >> 12);
  }
}

// ===================================================================================
// BIT READER / LETTORE DI BIT
// ===================================================================================
// EN: 32 bits refilled 16 at a time. The end is checked only at the refill: both bytes must be inside data plus
//     margin, otherwise the slice is corrupted and from there only bits 1 are supplied, without reading memory.
// IT: 32 bit ricaricati 16 alla volta. La fine si controlla solo alla ricarica: entrambi i byte devono stare dentro
//     dati piu' margine, altrimenti la fascia e' corrotta e da li' si forniscono solo bit 1, senza leggere memoria.
typedef struct {
  const uint8_t *p;
  const uint8_t *lim;
  uint32_t buf;
  int      cnt;
  int      bad;
} ScopeBits;

SCOPE_INLINE void brFill(ScopeBits *b) {
  if (b->cnt <= 16) {
    uint32_t w;
    if (b->p + 2 <= b->lim) { w = ((uint32_t)b->p[0] << 8) | b->p[1]; b->p += 2; }
    else { w = 0xFFFF; b->bad = 1; }
    b->buf |= w << (16 - b->cnt);
    b->cnt += 16;
  }
}
SCOPE_INLINE uint32_t brPeek(const ScopeBits *b, int n) { return b->buf >> (32 - n); }
SCOPE_INLINE void brSkip(ScopeBits *b, int n) { b->buf <<= n; b->cnt -= n; }
SCOPE_INLINE int brGet(ScopeBits *b, int n) { const uint32_t v = b->buf >> (32 - n); b->buf <<= n; b->cnt -= n; return (int)v; }
SCOPE_INLINE int extendBits(int v, int s) { return v < (1 << (s - 1)) ? v - (1 << s) + 1 : v; }

SCOPE_INLINE int huffDecode(ScopeBits *b, const ScopeHuff *h) {
  brFill(b);
  const uint32_t e = h->lut[brPeek(b, LUTB)];
  if (e >> 8) { brSkip(b, e >> 8); return e & 0xFF; }
  int l = LUTB + 1;
  int32_t code = (int32_t)brPeek(b, l);
  while (code > h->maxcode[l]) { l++; code = (int32_t)brPeek(b, l); }
  if (l > 16) { b->bad = 1; return 0; }   // no code of 16 bits or less: corrupted / nessun codice di 16 bit o meno: corrotta
  brSkip(b, l);
  return h->vals[h->valoff[l] + code];
}

// ===================================================================================
// INTEGER AAN IDCT / IDCT AAN INTERA (IJG jidctfst)
// ===================================================================================
#define CONST_BITS 8
#define PASS1_BITS 2
#define FIX_1_082392200 277
#define FIX_1_414213562 362
#define FIX_1_847759065 473
#define FIX_2_613125930 669
#define MUL(v, c) (((v) * (c)) >> CONST_BITS)
#define SH   (PASS1_BITS + 3)
#define RND  (1 << (SH - 1))
#define BIAS (((128 + 384) << SH) + RND)                         // rounding, +128 and clamp offset, once per row / arrotondamento, +128 e offset del clamp, una volta per riga
#define OUTB(x) CLAMP_TAB[((x) >> SH) & 1023]
#define OUTDC(x) CLAMP_TAB[((((x) + RND) >> SH) + 128 + 384) & 1023]

// EN: maxrow/maxcol: last row/column with non-zero coefficients (reduced variants up to 4x4).
// IT: maxrow/maxcol: ultima riga/colonna con coefficienti non nulli (varianti ridotte fino a 4x4).
static IRAM_ATTR void idctStore(const int32_t *in, uint8_t *out, int stride, int maxrow, int maxcol) {
  int32_t ws[64];
  const int ncol = maxcol <= 3 ? 4 : 8;
  for (int c = 0; c < ncol; c++) {
    const int32_t *ip = in + c;
    int32_t *wp = ws + c;
    if (c > maxcol) { for (int r = 0; r < 8; r++) wp[r * 8] = 0; continue; }
    if (maxrow == 0 || (ip[8] | ip[16] | ip[24] | ip[32] | ip[40] | ip[48] | ip[56]) == 0) {
      const int32_t d = ip[0];
      for (int r = 0; r < 8; r++) wp[r * 8] = d;
      continue;
    }
    if (maxrow <= 3) {   // reduced column: rows 0-3 only / colonna ridotta: solo righe 0-3
      const int32_t w0 = ip[0], w1 = ip[8], w2 = ip[16], w3 = ip[24];
      const int32_t t12 = MUL(w2, FIX_1_414213562) - w2;
      const int32_t t0 = w0 + w2, t3 = w0 - w2, t1 = w0 + t12, t2 = w0 - t12;
      const int32_t t7 = w1 + w3, t11 = MUL(w1 - w3, FIX_1_414213562);
      const int32_t z5 = MUL(w1 - w3, FIX_1_847759065);
      const int32_t t10 = MUL(w1, FIX_1_082392200) - z5, t12b = MUL(-w3, -FIX_2_613125930) + z5;
      const int32_t t6 = t12b - t7, t5 = t11 - t6, t4 = t10 + t5;
      wp[0] = t0 + t7; wp[56] = t0 - t7; wp[8] = t1 + t6; wp[48] = t1 - t6;
      wp[16] = t2 + t5; wp[40] = t2 - t5; wp[32] = t3 + t4; wp[24] = t3 - t4;
      continue;
    }
    int32_t t0 = ip[0], t1 = ip[16], t2 = ip[32], t3 = ip[48];
    int32_t t10 = t0 + t2, t11 = t0 - t2, t13 = t1 + t3, t12 = MUL(t1 - t3, FIX_1_414213562) - t13;
    t0 = t10 + t13; t3 = t10 - t13; t1 = t11 + t12; t2 = t11 - t12;
    int32_t t4 = ip[8], t5 = ip[24], t6 = ip[40], t7 = ip[56];
    const int32_t z13 = t6 + t5, z10 = t6 - t5, z11 = t4 + t7, z12 = t4 - t7;
    t7 = z11 + z13; t11 = MUL(z11 - z13, FIX_1_414213562);
    const int32_t z5 = MUL(z10 + z12, FIX_1_847759065);
    t10 = MUL(z12, FIX_1_082392200) - z5;
    t12 = MUL(z10, -FIX_2_613125930) + z5;
    t6 = t12 - t7; t5 = t11 - t6; t4 = t10 + t5;
    wp[0] = t0 + t7; wp[56] = t0 - t7; wp[8] = t1 + t6; wp[48] = t1 - t6;
    wp[16] = t2 + t5; wp[40] = t2 - t5; wp[32] = t3 + t4; wp[24] = t3 - t4;
  }
  for (int r = 0; r < 8; r++) {
    int32_t *wp = ws + r * 8;
    uint8_t *op = out + r * stride;
    wp[0] += BIAS;
    if (maxcol == 0) {
      const uint32_t v = OUTB(wp[0]) * 0x01010101u;
      ((uint32_t *)op)[0] = v; ((uint32_t *)op)[1] = v;
      continue;
    }
    if (maxcol <= 3) {
      const int32_t w0 = wp[0], w1 = wp[1], w2 = wp[2], w3 = wp[3];
      const int32_t t12 = MUL(w2, FIX_1_414213562) - w2;
      const int32_t t0 = w0 + w2, t3 = w0 - w2, t1 = w0 + t12, t2 = w0 - t12;
      const int32_t t7 = w1 + w3, t11 = MUL(w1 - w3, FIX_1_414213562);
      const int32_t z5 = MUL(w1 - w3, FIX_1_847759065);
      const int32_t t10 = MUL(w1, FIX_1_082392200) - z5, t12b = MUL(-w3, -FIX_2_613125930) + z5;
      const int32_t t6 = t12b - t7, t5 = t11 - t6, t4 = t10 + t5;
      op[0] = OUTB(t0 + t7); op[7] = OUTB(t0 - t7); op[1] = OUTB(t1 + t6); op[6] = OUTB(t1 - t6);
      op[2] = OUTB(t2 + t5); op[5] = OUTB(t2 - t5); op[4] = OUTB(t3 + t4); op[3] = OUTB(t3 - t4);
      continue;
    }
    const int32_t t10 = wp[0] + wp[4], t11 = wp[0] - wp[4], t13 = wp[2] + wp[6];
    const int32_t t12 = MUL(wp[2] - wp[6], FIX_1_414213562) - t13;
    const int32_t t0 = t10 + t13, t3 = t10 - t13, t1 = t11 + t12, t2 = t11 - t12;
    const int32_t z13 = wp[5] + wp[3], z10 = wp[5] - wp[3], z11 = wp[1] + wp[7], z12 = wp[1] - wp[7];
    const int32_t t7 = z11 + z13, t11b = MUL(z11 - z13, FIX_1_414213562);
    const int32_t z5 = MUL(z10 + z12, FIX_1_847759065);
    const int32_t t10b = MUL(z12, FIX_1_082392200) - z5, t12b = MUL(z10, -FIX_2_613125930) + z5;
    const int32_t t6 = t12b - t7, t5 = t11b - t6, t4 = t10b + t5;
    op[0] = OUTB(t0 + t7); op[7] = OUTB(t0 - t7); op[1] = OUTB(t1 + t6); op[6] = OUTB(t1 - t6);
    op[2] = OUTB(t2 + t5); op[5] = OUTB(t2 - t5); op[4] = OUTB(t3 + t4); op[3] = OUTB(t3 - t4);
  }
}

// ===================================================================================
// BLOCK 8x8 / BLOCCO 8x8
// ===================================================================================
// EN: Huffman, dequantization and IDCT of one block. Only-DC blocks without IDCT. coef is all zero at the exit.
// IT: Huffman, dequantizzazione e IDCT di un blocco. Blocchi con solo DC senza IDCT. coef e' tutto a zero all'uscita.
SCOPE_INLINE void decodeBlock(ScopeBits *b, const ScopeHuff *hdc, const int32_t *fac, const ScopeHuff *hac,
                                                              const int32_t *qs, int32_t *pred, int32_t *coef, uint8_t *out, int stride) {
  int s = huffDecode(b, hdc);
  int diff = 0;
  if (s) { brFill(b); diff = extendBits(brGet(b, s), s); }
  *pred += diff;
  uint8_t pos[64];
  int np = 0, maxrow = 0, maxcol = 0;
  coef[0] = *pred * qs[0];
  for (int k = 1; k < 64; k++) {
    brFill(b);
    const int32_t fe = fac[brPeek(b, LUTB)];
    int z;
    if (fe) {
      brSkip(b, fe & 15);
      k += (fe >> 4) & 15;
      if (k > 63) { b->bad = 1; break; }
      z = ZZ[k];
      coef[z] = (fe >> 8) * qs[z];
    } else {
      const int rs = huffDecode(b, hac);
      const int r = rs >> 4;
      s = rs & 15;
      if (s == 0) { if (r == 15) { k += 15; continue; } break; }   // ZRL / EOB
      k += r;
      if (k > 63) { b->bad = 1; break; }
      brFill(b);
      const int v = extendBits(brGet(b, s), s);
      z = ZZ[k];
      coef[z] = v * qs[z];
    }
    pos[np++] = (uint8_t)z;
    if ((z & 7) > maxcol) maxcol = z & 7;
    if ((z >> 3) > maxrow) maxrow = z >> 3;
  }
  if (np == 0) {
    const uint32_t w = OUTDC(coef[0]) * 0x01010101u;
    for (int r = 0; r < 8; r++) { ((uint32_t *)(out + r * stride))[0] = w; ((uint32_t *)(out + r * stride))[1] = w; }
    coef[0] = 0;
    return;
  }
  idctStore(coef, out, stride, maxrow, maxcol);
  coef[0] = 0;
  for (int i = 0; i < np; i++) coef[pos[i]] = 0;
}

// ===================================================================================
// COLOUR / COLORE
// ===================================================================================
// EN: YCbCr 4:2:0 -> RGB565 in display order, 2x2 squares and 32-bit writes. Partial MCU: only the written span.
// IT: YCbCr 4:2:0 -> RGB565 in ordine del display, quadrati 2x2 e scritture a 32 bit. MCU parziale: solo il tratto scritto.
static IRAM_ATTR void colorMcu(uint8_t sy, uint8_t mx, const uint8_t *Y, const uint8_t *Cb, const uint8_t *Cr, uint16_t *work) {
  const int x0 = mx * 16;
  if (MCU_FULL[sy][mx]) {
    uint16_t *d0 = work + x0;
    for (int cy = 0; cy < 8; cy++, d0 += 2 * W, Y += 32, Cb += 8, Cr += 8) {
      uint32_t *o0 = (uint32_t *)d0, *o1 = (uint32_t *)(d0 + W);
      for (int c = 0; c < 8; c++) {
        const int cb = Cb[c], cr = Cr[c];
        const uint16_t *a = R5 + CRR16[cr], *bb = B5 + CBB16[cb], *g = G6 + ((CBG[cb] + CRG[cr]) >> 16);
        const int p0 = Y[2 * c], p1 = Y[2 * c + 1], p2 = Y[16 + 2 * c], p3 = Y[17 + 2 * c];
        o0[c] = (uint32_t)(a[p0] | g[p0] | bb[p0]) | ((uint32_t)(a[p1] | g[p1] | bb[p1]) << 16);
        o1[c] = (uint32_t)(a[p2] | g[p2] | bb[p2]) | ((uint32_t)(a[p3] | g[p3] | bb[p3]) << 16);
      }
    }
    return;
  }
  for (int cr2 = 0; cr2 < 8; cr2++) {
    const int ys = sy * 16 + cr2 * 2;
    if (WXS[ys] > WXE[ys] && WXS[ys + 1] > WXE[ys + 1]) continue;
    const uint16_t *pr[8], *pg[8], *pb[8];
    for (int c = 0; c < 8; c++) {
      const int cb = Cb[cr2 * 8 + c], cr = Cr[cr2 * 8 + c];
      pr[c] = R5 + CRR16[cr]; pb[c] = B5 + CBB16[cb]; pg[c] = G6 + ((CBG[cb] + CRG[cr]) >> 16);
    }
    for (int rr = 0; rr < 2; rr++) {
      const int r = cr2 * 2 + rr, yg = sy * 16 + r;
      int xs = WXS[yg] - x0, xe = WXE[yg] - x0;
      if (xs < 0) xs = 0;
      if (xe > 15) xe = 15;
      if (xs > xe) continue;
      const uint8_t *yp = Y + r * 16;
      uint16_t *dst = work + r * W + x0;
      int x = xs;
      if (x & 1) { const int c = x >> 1, yv = yp[x]; dst[x] = pr[c][yv] | pg[c][yv] | pb[c][yv]; x++; }
      for (; x < xe; x += 2) {
        const int c = x >> 1;
        const uint16_t *a = pr[c], *g = pg[c], *bb = pb[c];
        const int y0 = yp[x], y1 = yp[x + 1];
        *(uint32_t *)(dst + x) = (uint32_t)(a[y0] | g[y0] | bb[y0]) | ((uint32_t)(a[y1] | g[y1] | bb[y1]) << 16);
      }
      if (x == xe) { const int c = x >> 1, yv = yp[x]; dst[x] = pr[c][yv] | pg[c][yv] | pb[c][yv]; }
    }
  }
}

// ===================================================================================
// SLICE / FASCIA
// ===================================================================================
IRAM_ATTR bool scopeDecSlice(const uint8_t *data, uint16_t len, uint8_t sy, uint8_t q, uint16_t *work, ScopeDecIdleFn idle) {
  if (sy >= SCOPE_DEC_SLICES || q < 1 || q > 100) return false;
  setQuant(q);
  static int32_t coef[64];   // always zero between blocks / sempre a zero tra un blocco e l'altro
  alignas(4) uint8_t Y[16 * 16];
  alignas(4) uint8_t Cb[64];
  alignas(4) uint8_t Cr[64];
  ScopeBits b = { data, data + len + SCOPE_DEC_MARGIN, 0, 0, 0 };
  int32_t py = 0, pb = 0, pr = 0;   // DC predictors reset at every slice / predittori DC azzerati a ogni fascia

  for (uint8_t mx = 0; mx < 15; mx++) {
    if (!MCU_OK[sy][mx]) continue;   // MCU outside the circle: not in the bitstream / MCU fuori dal cerchio: non c'e' nel bitstream
    decodeBlock(&b, &HDL, FACL, &HAL, QSL, &py, coef, Y,       16);
    decodeBlock(&b, &HDL, FACL, &HAL, QSL, &py, coef, Y + 8,   16);
    decodeBlock(&b, &HDL, FACL, &HAL, QSL, &py, coef, Y + 128, 16);
    decodeBlock(&b, &HDL, FACL, &HAL, QSL, &py, coef, Y + 136, 16);
    decodeBlock(&b, &HDC, FACC, &HAC, QSC, &pb, coef, Cb,      8);
    decodeBlock(&b, &HDC, FACC, &HAC, QSC, &pr, coef, Cr,      8);
    if (b.bad) return false;   // corrupted: the rest is not needed / corrotta: il resto non serve
    colorMcu(sy, mx, Y, Cb, Cr, work);
    if (idle) idle();
  }
  // EN: Logical limit: the bits consumed must not exceed those of the received data.
  // IT: Limite logico: i bit consumati non devono superare quelli dei dati ricevuti.
  const int32_t consumed = (int32_t)(b.p - data) * 8 - b.cnt;
  return !b.bad && consumed <= (int32_t)len * 8;
}

void IRAM_ATTR scopeDecCopySlice(const uint16_t *work, uint8_t sy, uint16_t *block) {
  for (int r = 0; r < SCOPE_DEC_SLICE_H; r++) {
    const int y = sy * SCOPE_DEC_SLICE_H + r, xs = WXS[y], xe = WXE[y];
    if (xs <= xe) memcpy(block + r * W + xs, work + r * W + xs, (size_t)(xe - xs + 1) * 2);
  }
}

void scopeDecBlackSlice(uint8_t sy, uint16_t *block) {
  for (int r = 0; r < SCOPE_DEC_SLICE_H; r++) {
    const int y = sy * SCOPE_DEC_SLICE_H + r, xs = WXS[y], xe = WXE[y];
    if (xs <= xe) memset(block + r * W + xs, 0, (size_t)(xe - xs + 1) * 2);
  }
}

// ===================================================================================
// ROTATION / ROTAZIONE
// ===================================================================================
uint16_t scopeRotWindowX(uint8_t sy) { return OUT_X[sy]; }
uint16_t scopeRotWindowW(uint8_t sy) { return OUT_W[sy]; }

void scopeRotAngle(int angle_deg, int32_t *c, int32_t *s) {
  const float a = (float)angle_deg * 3.14159265f / 180.0f;
  *c = (int32_t)lroundf(cosf(a) * 65536.0f);
  *s = (int32_t)lroundf(sinf(a) * 65536.0f);
}

// EN: Output pixel (x, y) shows the frame pixel u = 119.5 + c*dx + s*dy, v = 119.5 - s*dx + c*dy (dx = x - 119.5,
//     dy = y - 119.5), rounded to the nearest and limited to 0-239: on the edge the rounding can give -1 or 240.
//     Three faster paths give the same pixels: at 0 degrees each row is a plain copy; a row whose source stays
//     inside 0-239 needs no limit (u and v change linearly along the row, so its two ends decide); up to about
//     8 degrees such a row is copied in runs of pixels with the same source row and consecutive source columns.
// IT: Il pixel di uscita (x, y) mostra il pixel del frame u = 119,5 + c*dx + s*dy, v = 119,5 - s*dx + c*dy (dx = x - 119,5,
//     dy = y - 119,5), arrotondato al piu' vicino e limitato a 0-239: sul bordo l'arrotondamento puo' dare -1 o 240.
//     Tre percorsi piu' veloci danno gli stessi pixel: a 0 gradi ogni riga e' una copia semplice; una riga la cui
//     sorgente resta dentro 0-239 non ha bisogno del limite (u e v cambiano in modo lineare lungo la riga, quindi
//     decidono i suoi due estremi); fino a circa 8 gradi una riga cosi' e' copiata a tratti di pixel con la stessa
//     riga sorgente e colonne sorgente consecutive.
#define SCOPE_ROT_RUN_SIN 9200   // sin(angle) * 65536 up to which the runs are used (about 8 degrees) / fino a cui si usano i tratti (circa 8 gradi)
IRAM_ATTR void scopeRotSlice(const uint16_t *const *rows, uint8_t sy, int32_t c, int32_t s, uint16_t *out) {
  const int ox = OUT_X[sy], ow = OUT_W[sy];
  const bool same = (s == 0 && c == 65536);   // 0 degrees / 0 gradi
  const int32_t d = 65536 - c, as = s < 0 ? -s : s;
  // EN: The run lengths need fu / d and fv / |s|, done as (f * m) >> 32 with m = 2^32 / divisor + 1: exact for
  //     f < 65536 and divisor >= 2 (checked on every pair). With whole degrees d >= 10 and |s| >= 1144.
  // IT: Le lunghezze dei tratti richiedono fu / d e fv / |s|, fatte come (f * m) >> 32 con m = 2^32 / divisore + 1: esatto per
  //     f < 65536 e divisore >= 2 (verificato su ogni coppia). Con gradi interi d >= 10 e |s| >= 1144.
  const bool runs = !same && c > 32768 && as >= 2 && as <= SCOPE_ROT_RUN_SIN && d != 1;
  const uint32_t md = (runs && d >= 2) ? (uint32_t)((((uint64_t)1 << 32) / (uint32_t)d) + 1) : 0;   // 0: d = 0, u never changes row / u non cambia mai
  const uint32_t ms = runs ? (uint32_t)((((uint64_t)1 << 32) / (uint32_t)as) + 1) : 0;
  for (int r = 0; r < SCOPE_DEC_SLICE_H; r++) {
    const int y = sy * SCOPE_DEC_SLICE_H + r, xs = VXS[y], xe = VXE[y];
    uint16_t *o = out + r * ow - ox;   // o[x] = pixel x of the display / o[x] = pixel x del display
    for (int x = ox; x < xs; x++) o[x] = 0;               // black outside the circle / nero fuori dal cerchio
    for (int x = xe + 1; x < ox + ow; x++) o[x] = 0;
    if (same) { memcpy(o + xs, rows[y] + xs, (size_t)(xe - xs + 1) * 2); continue; }
    const int32_t dx = ((xs << 1) - 239) * 32768, dy = ((y << 1) - 239) * 32768;   // 16.16
    int32_t u = (120 << 16) + (int32_t)(((int64_t)c * dx + (int64_t)s * dy) >> 16);   // 119.5 + 0.5 for the rounding / 119,5 + 0,5 per l'arrotondamento
    int32_t v = (120 << 16) + (int32_t)(((int64_t)c * dy - (int64_t)s * dx) >> 16);
    const int32_t n = xe - xs, u1 = u + n * c, v1 = v - n * s;   // source at the end of the row / sorgente a fine riga
    const int32_t umin = u < u1 ? u : u1, umax = u < u1 ? u1 : u;
    const int32_t vmin = v < v1 ? v : v1, vmax = v < v1 ? v1 : v;
    if (umin >= 0 && vmin >= 0 && umax < (W << 16) && vmax < (H << 16)) {
      if (runs) {
        int x = xs;
        while (x <= xe) {
          // EN: the run ends when u falls one column behind (after fu / d + 1 pixels) or v changes row (fv / |s| + 1).
          // IT: il tratto finisce quando u resta indietro di una colonna (dopo fu / d + 1 pixel) o v cambia riga (fv / |s| + 1).
          const uint32_t fu = (uint32_t)(u & 0xFFFF);
          const uint32_t fv = (uint32_t)(s > 0 ? (v & 0xFFFF) : (0xFFFF - (v & 0xFFFF)));
          int32_t k = xe - x + 1;
          if (md) { const int32_t ku = (int32_t)(((uint64_t)fu * md) >> 32) + 1; if (ku < k) k = ku; }
          const int32_t kv = (int32_t)(((uint64_t)fv * ms) >> 32) + 1;
          if (kv < k) k = kv;
          const uint16_t *src = rows[v >> 16] + (u >> 16);
          uint16_t *dst = o + x;
          for (int32_t j = 0; j < k; j++) dst[j] = src[j];
          x += k;
          u += k * c;
          v -= k * s;
        }
      } else {
        for (int x = xs; x <= xe; x++) {
          o[x] = rows[v >> 16][u >> 16];
          u += c;
          v -= s;
        }
      }
      continue;
    }
    for (int x = xs; x <= xe; x++) {
      int ui = u >> 16, vi = v >> 16;
      ui = ui < 0 ? 0 : (ui > W - 1 ? W - 1 : ui);
      vi = vi < 0 ? 0 : (vi > H - 1 ? H - 1 : vi);
      o[x] = rows[vi][ui];
      u += c;
      v -= s;
    }
  }
}

// ===================================================================================
// FRAME IN YCbCr / FRAME IN YCbCr
// ===================================================================================
uint16_t scopeDecMcuCount() { return MCU_N; }

// EN: Same bitstream and same blocks as scopeDecSlice: every coded MCU goes into its place of the frame as it is decoded.
// IT: Stesso bitstream e stessi blocchi di scopeDecSlice: ogni MCU codificato va nel suo posto del frame appena decodificato.
IRAM_ATTR bool scopeDecSliceYcc(const uint8_t *data, uint16_t len, uint8_t sy, uint8_t q, uint8_t *const *mcu, ScopeDecIdleFn idle) {
  if (sy >= SCOPE_DEC_SLICES || q < 1 || q > 100) return false;
  setQuant(q);
  static int32_t coef[64];   // always zero between blocks / sempre a zero tra un blocco e l'altro
  ScopeBits b = { data, data + len + SCOPE_DEC_MARGIN, 0, 0, 0 };
  int32_t py = 0, pb = 0, pr = 0;   // DC predictors reset at every slice / predittori DC azzerati a ogni fascia
  uint8_t *const *row = mcu + sy * 16;

  for (uint8_t mx = 0; mx < 15; mx++) {
    if (!MCU_OK[sy][mx]) continue;   // MCU outside the circle: not in the bitstream / MCU fuori dal cerchio: non c'e' nel bitstream
    uint8_t *m = row[mx];
    decodeBlock(&b, &HDL, FACL, &HAL, QSL, &py, coef, m,       16);
    decodeBlock(&b, &HDL, FACL, &HAL, QSL, &py, coef, m + 8,   16);
    decodeBlock(&b, &HDL, FACL, &HAL, QSL, &py, coef, m + 128, 16);
    decodeBlock(&b, &HDL, FACL, &HAL, QSL, &py, coef, m + 136, 16);
    decodeBlock(&b, &HDC, FACC, &HAC, QSC, &pb, coef, m + 256, 8);
    decodeBlock(&b, &HDC, FACC, &HAC, QSC, &pr, coef, m + 320, 8);
    if (b.bad) return false;   // corrupted: the rest is not needed / corrotta: il resto non serve
    if (idle) idle();
  }
  // EN: Logical limit: the bits consumed must not exceed those of the received data.
  // IT: Limite logico: i bit consumati non devono superare quelli dei dati ricevuti.
  const int32_t consumed = (int32_t)(b.p - data) * 8 - b.cnt;
  return !b.bad && consumed <= (int32_t)len * 8;
}

void scopeDecBlackSliceYcc(uint8_t sy, uint8_t *const *mcu) {
  for (uint8_t mx = 0; mx < 15; mx++) if (MCU_OK[sy][mx]) memcpy(mcu[sy * 16 + mx], BLACK_MCU, SCOPE_DEC_MCU_BYTES);
}

void scopeDecMcuTable(uint8_t *const *piece, const uint8_t *count, uint8_t nPiece, uint8_t **mcu) {
  uint8_t pi = 0, k = 0;
  for (int my = 0; my < SCOPE_DEC_SLICES; my++) {
    for (int mx = 0; mx < 15; mx++) {
      if (!MCU_OK[my][mx]) { mcu[my * 16 + mx] = BLACK_MCU; continue; }
      while (pi < nPiece && k >= count[pi]) { pi++; k = 0; }
      mcu[my * 16 + mx] = (pi < nPiece) ? piece[pi] + k++ * SCOPE_DEC_MCU_BYTES : BLACK_MCU;
    }
    mcu[my * 16 + 15] = BLACK_MCU;   // never read / mai letto
  }
}

uint16_t scopeRotBandX(uint8_t band) { return BAND_X[band]; }
uint16_t scopeRotBandW(uint8_t band) { return BAND_W[band]; }

// EN: Pixel (u, v) of the frame in RGB565, display order: same tables and same sums as colorMcu. The three tables
//     of the colour depend only on the chroma sample (2x2 pixels): they are computed again only when it changes.
// IT: Pixel (u, v) del frame in RGB565, ordine del display: stesse tabelle e stesse somme di colorMcu. Le tre tabelle
//     del colore dipendono solo dal campione di croma (2x2 pixel): si ricalcolano solo quando cambia.
typedef struct {
  uint32_t key;   // chroma sample (v / 2, u / 2) of the tables below / campione di croma (v / 2, u / 2) delle tabelle qui sotto
  const uint16_t *r, *g, *b;
} ScopeYccCache;

SCOPE_INLINE uint16_t yccPixel(const uint8_t *const *mcu, uint32_t u, uint32_t v, ScopeYccCache *k) {
  const uint8_t *m = mcu[((v >> 4) << 4) | (u >> 4)];
  const int yy = m[((v & 15) << 4) | (u & 15)];
  const uint32_t key = ((v >> 1) << 8) | (u >> 1);
  if (key != k->key) {
    k->key = key;
    const uint8_t *cp = m + 256 + ((v & 14) << 2) + ((u & 15) >> 1);
    const int cb = cp[0], cr = cp[64];
    k->r = R5 + CRR16[cr]; k->g = G6 + ((CBG[cb] + CRG[cr]) >> 16); k->b = B5 + CBB16[cb];
  }
  return (uint16_t)(k->r[yy] | k->g[yy] | k->b[yy]);
}

// EN: n consecutive pixels of row v of the frame from column u: the colour of a chroma pair is computed once.
// IT: n pixel consecutivi della riga v del frame dalla colonna u: il colore di una coppia di croma si calcola una volta.
SCOPE_INLINE void yccSpan(const uint8_t *const *mcu, int u, int v, int n, uint16_t *dst) {
  const uint8_t *const *mrow = mcu + ((v >> 4) << 4);
  const int ry = (v & 15) << 4, rc = 256 + (((v & 15) >> 1) << 3);
  const int ue = u + n;
  while (u < ue) {
    const uint8_t *m = mrow[u >> 4];
    const uint8_t *yp = m + ry, *cp = m + rc;
    int lx = u & 15;
    int end = 16 - lx;                       // pixels left in this MCU / pixel rimasti in questo MCU
    if (end > ue - u) end = ue - u;
    const int le = lx + end;
    while (lx < le) {
      const int c2 = lx >> 1, cb = cp[c2], cr = cp[c2 + 64];
      const uint16_t *a = R5 + CRR16[cr], *g = G6 + ((CBG[cb] + CRG[cr]) >> 16), *bb = B5 + CBB16[cb];
      int yv = yp[lx];
      *dst++ = (uint16_t)(a[yv] | g[yv] | bb[yv]);
      lx++;
      if ((lx & 1) && lx < le) { yv = yp[lx]; *dst++ = (uint16_t)(a[yv] | g[yv] | bb[yv]); lx++; }
    }
    u += end;
  }
}

// EN: Two rows y (even) and y + 1, n pixels from column x (x even, n even): they share the chroma samples, so the
//     colour tables are computed once every 4 pixels, as in colorMcu.
// IT: Due righe y (pari) e y + 1, n pixel dalla colonna x (x pari, n pari): condividono i campioni di croma, quindi le
//     tabelle del colore si calcolano una volta ogni 4 pixel, come in colorMcu.
SCOPE_INLINE void yccQuad(const uint8_t *const *mcu, int x, int y, int n, uint16_t *da, uint16_t *db) {
  const uint8_t *const *mrow = mcu + ((y >> 4) << 4);
  const int ry = (y & 15) << 4, rc = 256 + (((y & 15) >> 1) << 3);
  const int xe = x + n;
  while (x < xe) {
    const uint8_t *m = mrow[x >> 4];
    const uint8_t *yp = m + ry + (x & 15), *cp = m + rc + ((x & 15) >> 1);
    int cnt = 16 - (x & 15);
    if (cnt > xe - x) cnt = xe - x;
    x += cnt;
    for (; cnt > 0; cnt -= 2, yp += 2, cp++, da += 2, db += 2) {
      const int cb = cp[0], cr = cp[64];
      const uint16_t *a = R5 + CRR16[cr], *g = G6 + ((CBG[cb] + CRG[cr]) >> 16), *bb = B5 + CBB16[cb];
      int p = yp[0];  da[0] = (uint16_t)(a[p] | g[p] | bb[p]);
      p = yp[1];      da[1] = (uint16_t)(a[p] | g[p] | bb[p]);
      p = yp[16];     db[0] = (uint16_t)(a[p] | g[p] | bb[p]);
      p = yp[17];     db[1] = (uint16_t)(a[p] | g[p] | bb[p]);
    }
  }
}

// EN: Same geometry and same paths as scopeRotSlice, on bands of 8 rows; only the reading of the pixel changes.
// IT: Stessa geometria e stessi percorsi di scopeRotSlice, su fasce da 8 righe; cambia solo la lettura del pixel.
IRAM_ATTR void scopeRotBandYcc(const uint8_t *const *mcu, uint8_t band, int32_t c, int32_t s, uint16_t *out) {
  const int ox = BAND_X[band], ow = BAND_W[band];
  const bool same = (s == 0 && c == 65536);   // 0 degrees / 0 gradi
  const int32_t d = 65536 - c, as = s < 0 ? -s : s;
  const bool runs = !same && c > 32768 && as >= 2 && as <= SCOPE_ROT_RUN_SIN && d != 1;
  const uint32_t md = (runs && d >= 2) ? (uint32_t)((((uint64_t)1 << 32) / (uint32_t)d) + 1) : 0;   // 0: d = 0, u never changes row / u non cambia mai
  const uint32_t ms = runs ? (uint32_t)((((uint64_t)1 << 32) / (uint32_t)as) + 1) : 0;
  if (same) {
    // EN: 0 degrees: the rows two by two; the common part of the two spans 4 pixels at a time, the ends one row at a time.
    // IT: 0 gradi: le righe a due a due; la parte comune dei due tratti 4 pixel alla volta, le estremita' una riga alla volta.
    for (int r = 0; r < SCOPE_ROT_BAND_H; r += 2) {
      const int y = band * SCOPE_ROT_BAND_H + r;
      uint16_t *oa = out + r * ow - ox, *ob = oa + ow;
      const int xsa = VXS[y], xea = VXE[y], xsb = VXS[y + 1], xeb = VXE[y + 1];
      for (int x = ox; x < xsa; x++) oa[x] = 0;               // black outside the circle / nero fuori dal cerchio
      for (int x = xea + 1; x < ox + ow; x++) oa[x] = 0;
      for (int x = ox; x < xsb; x++) ob[x] = 0;
      for (int x = xeb + 1; x < ox + ow; x++) ob[x] = 0;
      const int cs = ((xsa > xsb ? xsa : xsb) + 1) & ~1;        // even / pari
      const int ce = ((((xea < xeb ? xea : xeb) + 1) & ~1)) - 1; // odd / dispari
      if (cs > ce) {
        yccSpan(mcu, xsa, y, xea - xsa + 1, oa + xsa);
        yccSpan(mcu, xsb, y + 1, xeb - xsb + 1, ob + xsb);
        continue;
      }
      if (xsa < cs) yccSpan(mcu, xsa, y, cs - xsa, oa + xsa);
      if (xsb < cs) yccSpan(mcu, xsb, y + 1, cs - xsb, ob + xsb);
      yccQuad(mcu, cs, y, ce - cs + 1, oa + cs, ob + cs);
      if (xea > ce) yccSpan(mcu, ce + 1, y, xea - ce, oa + ce + 1);
      if (xeb > ce) yccSpan(mcu, ce + 1, y + 1, xeb - ce, ob + ce + 1);
    }
    return;
  }
  for (int r = 0; r < SCOPE_ROT_BAND_H; r++) {
    const int y = band * SCOPE_ROT_BAND_H + r, xs = VXS[y], xe = VXE[y];
    uint16_t *o = out + r * ow - ox;   // o[x] = pixel x of the display / o[x] = pixel x del display
    for (int x = ox; x < xs; x++) o[x] = 0;               // black outside the circle / nero fuori dal cerchio
    for (int x = xe + 1; x < ox + ow; x++) o[x] = 0;
    const int32_t dx = ((xs << 1) - 239) * 32768, dy = ((y << 1) - 239) * 32768;   // 16.16
    int32_t u = (120 << 16) + (int32_t)(((int64_t)c * dx + (int64_t)s * dy) >> 16);   // 119.5 + 0.5 for the rounding / 119,5 + 0,5 per l'arrotondamento
    int32_t v = (120 << 16) + (int32_t)(((int64_t)c * dy - (int64_t)s * dx) >> 16);
    const int32_t n = xe - xs, u1 = u + n * c, v1 = v - n * s;   // source at the end of the row / sorgente a fine riga
    const int32_t umin = u < u1 ? u : u1, umax = u < u1 ? u1 : u;
    const int32_t vmin = v < v1 ? v : v1, vmax = v < v1 ? v1 : v;
    if (umin >= 0 && vmin >= 0 && umax < (W << 16) && vmax < (H << 16)) {
      if (runs) {
        int x = xs;
        while (x <= xe) {
          const uint32_t fu = (uint32_t)(u & 0xFFFF);
          const uint32_t fv = (uint32_t)(s > 0 ? (v & 0xFFFF) : (0xFFFF - (v & 0xFFFF)));
          int32_t k = xe - x + 1;
          if (md) { const int32_t ku = (int32_t)(((uint64_t)fu * md) >> 32) + 1; if (ku < k) k = ku; }
          const int32_t kv = (int32_t)(((uint64_t)fv * ms) >> 32) + 1;
          if (kv < k) k = kv;
          yccSpan(mcu, u >> 16, v >> 16, k, o + x);
          x += k;
          u += k * c;
          v -= k * s;
        }
      } else {
        ScopeYccCache k = { 0xFFFFFFFFu, R5, G6, B5 };
        for (int x = xs; x <= xe; x++) {
          o[x] = yccPixel(mcu, (uint32_t)u >> 16, (uint32_t)v >> 16, &k);
          u += c;
          v -= s;
        }
      }
      continue;
    }
    ScopeYccCache k = { 0xFFFFFFFFu, R5, G6, B5 };
    for (int x = xs; x <= xe; x++) {
      int ui = u >> 16, vi = v >> 16;
      ui = ui < 0 ? 0 : (ui > W - 1 ? W - 1 : ui);
      vi = vi < 0 ? 0 : (vi > H - 1 ? H - 1 : vi);
      o[x] = yccPixel(mcu, (uint32_t)ui, (uint32_t)vi, &k);
      u += c;
      v -= s;
    }
  }
}
