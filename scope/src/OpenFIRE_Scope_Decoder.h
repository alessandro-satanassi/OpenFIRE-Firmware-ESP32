/*!
 * @file OpenFIRE_Scope_Decoder.h
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

#pragma once

#include <stdint.h>

// ===================================================================================
// VIDEO FORMAT OF THE SCOPE / FORMATO VIDEO DEL MIRINO
// ===================================================================================
// EN: Frame 240x240 in 15 slices of 16 rows, each slice decoded by itself. Bitstream derived from baseline JPEG:
//     4:2:0, MCU 16x16 (Y Y Y Y Cb Cr), standard tables K.1-K.6 of ITU-T T.81, DC predictors reset at every slice,
//     MCU outside the circle omitted, no header, marker or byte stuffing, end of slice padded with bits 1.
//     Output RGB565 with the bytes already in the order of the display (big-endian).
// IT: Frame 240x240 in 15 fasce da 16 righe, ogni fascia si decodifica da sola. Bitstream derivato dal JPEG baseline:
//     4:2:0, MCU 16x16 (Y Y Y Y Cb Cr), tabelle standard K.1-K.6 di ITU-T T.81, predittori DC azzerati a ogni fascia,
//     MCU fuori dal cerchio omessi, niente intestazione, marker o byte stuffing, fine fascia riempita con bit 1.
//     Uscita RGB565 con i byte gia' nell'ordine del display (big-endian).
#define SCOPE_DEC_W          240
#define SCOPE_DEC_H          240
#define SCOPE_DEC_SLICES     15
#define SCOPE_DEC_SLICE_H    16
#define SCOPE_DEC_MARGIN     4     // bytes 0xFF that the receiver puts after the data / byte 0xFF che chi riceve mette dopo i dati

// EN: Called after every MCU: the video module uses it to start the next output slice on the DMA.
// IT: Chiamata dopo ogni MCU: il modulo video la usa per far partire sul DMA la fascia di uscita successiva.
typedef void (*ScopeDecIdleFn)(void);

// EN: Tables (Huffman, fast AC, colour, clamp, circle). Once at power on.
// IT: Tabelle (Huffman, AC veloce, colore, clamp, cerchio). Una volta all'accensione.
void scopeDecInit();

// EN: Decodes the slice sy (0-14) with quality q (1-100) into work (240x16 RGB565, display order).
//     data: len bytes followed by SCOPE_DEC_MARGIN bytes 0xFF. Returns true if the slice is valid: every MCU complete
//     and no bit read beyond the data. A false slice must not be copied into the frame.
// IT: Decodifica la fascia sy (0-14) con qualita' q (1-100) in work (240x16 RGB565, ordine del display).
//     data: len byte seguiti da SCOPE_DEC_MARGIN byte 0xFF. Torna true se la fascia e' valida: tutti gli MCU completi
//     e nessun bit letto oltre i dati. Una fascia false non va copiata nel frame.
bool scopeDecSlice(const uint8_t *data, uint16_t len, uint8_t sy, uint8_t q, uint16_t *work, ScopeDecIdleFn idle);

// EN: The frame is kept as 15 separate blocks of 240x16 pixels, one per slice (no need of 115 KB in one piece).
//     Copies into the block of slice sy only the pixels written by the decoder (circle and coded MCU).
// IT: Il frame si tiene in 15 blocchi separati da 240x16 pixel, uno per fascia (non servono 115 KB in un pezzo solo).
//     Copia nel blocco della fascia sy solo i pixel scritti dal decoder (cerchio e MCU codificati).
void scopeDecCopySlice(const uint16_t *work, uint8_t sy, uint16_t *block);

// EN: Black on the written pixels of the block of slice sy (slice missing for too long).
// IT: Nero sui pixel scritti del blocco della fascia sy (fascia mancante da troppo tempo).
void scopeDecBlackSlice(uint8_t sy, uint16_t *block);

// ===================================================================================
// ROLL ROTATION / ROTAZIONE ROLL
// ===================================================================================
// EN: Window of the output slice sy on the display: only as wide as the visible circle in those 16 rows.
// IT: Finestra della fascia di uscita sy sul display: larga solo quanto il cerchio visibile in quelle 16 righe.
uint16_t scopeRotWindowX(uint8_t sy);
uint16_t scopeRotWindowW(uint8_t sy);

// EN: cos and sin in 16.16 of the rotation of the image on the display, angle in degrees, positive = clockwise.
// IT: cos e sin in 16.16 della rotazione dell'immagine sul display, angolo in gradi, positivo = orario.
void scopeRotAngle(int angle_deg, int32_t *c, int32_t *s);

// EN: Nearest-neighbour rotation around (119.5, 119.5) of the frame (rows[y] = row y of the frame) into the
//     output slice sy (scopeRotWindowW(sy) x 16 pixels, black outside the visible circle).
// IT: Rotazione nearest-neighbour attorno a (119,5; 119,5) del frame (rows[y] = riga y del frame) nella fascia
//     di uscita sy (scopeRotWindowW(sy) x 16 pixel, nero fuori dal cerchio visibile).
void scopeRotSlice(const uint16_t *const *rows, uint8_t sy, int32_t c, int32_t s, uint16_t *out);
