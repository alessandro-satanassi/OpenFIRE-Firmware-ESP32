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

// ===================================================================================
// FRAME IN YCbCr / FRAME IN YCbCr
// ===================================================================================
// EN: Frame kept as it arrives: only the coded MCU, in YCbCr 4:2:0 (201 MCU x 384 bytes = 75 KB instead of the
//     115 KB of the RGB565 frame), in pieces of any size. Every MCU = Y 16x16, Cb 8x8, Cr 8x8. The colour is made
//     by the rotation, with the same tables and the same pixels as the RGB565 path.
// IT: Frame tenuto come arriva: solo gli MCU codificati, in YCbCr 4:2:0 (201 MCU x 384 byte = 75 KB invece dei
//     115 KB del frame RGB565), a pezzi di qualunque misura. Ogni MCU = Y 16x16, Cb 8x8, Cr 8x8. Il colore lo fa
//     la rotazione, con le stesse tabelle e gli stessi pixel del percorso RGB565.
#define SCOPE_DEC_MCU_BYTES  384
#define SCOPE_DEC_MCU_TABLE  (SCOPE_DEC_SLICES * 16)

// EN: Coded MCU of the frame (SCOPE_DEC_MCU_CODED with the 240x240 circle): the frame needs scopeDecMcuCount() x
//     SCOPE_DEC_MCU_BYTES bytes.
// IT: MCU codificati del frame (SCOPE_DEC_MCU_CODED col cerchio 240x240): il frame vuole scopeDecMcuCount() x
//     SCOPE_DEC_MCU_BYTES byte.
#define SCOPE_DEC_MCU_CODED  201
uint16_t scopeDecMcuCount();

// EN: Table of the 15 x 15 MCU of the frame (mcu[my * 16 + mx], rows of 16 for a shift instead of a multiplication).
//     The frame can be in pieces (any free memory): piece i holds count[i] coded MCU, in order of slice and column.
//     A not coded MCU points to a black MCU, never written. Once, after the pieces are allocated.
// IT: Tabella dei 15 x 15 MCU del frame (mcu[my * 16 + mx], righe da 16 per uno scorrimento invece di una moltiplicazione).
//     Il frame puo' essere a pezzi (qualunque memoria libera): il pezzo i tiene count[i] MCU codificati, in ordine di
//     fascia e colonna. Un MCU non codificato punta a un MCU nero, mai scritto. Una volta, dopo l'allocazione dei pezzi.
void scopeDecMcuTable(uint8_t *const *piece, const uint8_t *count, uint8_t nPiece, uint8_t **mcu);

// EN: As scopeDecSlice, but straight into the MCU of the frame in YCbCr, without colour. If not valid the slice is
//     partly written: the caller makes it black (scopeDecBlackSliceYcc).
// IT: Come scopeDecSlice, ma direttamente negli MCU del frame in YCbCr, senza colore. Se non valida la fascia e'
//     scritta in parte: chi chiama la fa nera (scopeDecBlackSliceYcc).
bool scopeDecSliceYcc(const uint8_t *data, uint16_t len, uint8_t sy, uint8_t q, uint8_t *const *mcu, ScopeDecIdleFn idle);

// EN: Black slice sy of the frame (Y 0, Cb and Cr 128: RGB565 0 as in the RGB565 path).
// IT: Fascia sy del frame nera (Y 0, Cb e Cr 128: RGB565 0 come nel percorso RGB565).
void scopeDecBlackSliceYcc(uint8_t sy, uint8_t *const *mcu);

// EN: Output bands of 8 rows (half the memory of the 16-row slices for the DMA buffers), window as wide as the
//     visible circle in those rows.
// IT: Fasce di uscita da 8 righe (meta' memoria delle fasce da 16 righe per i buffer del DMA), finestra larga quanto
//     il cerchio visibile in quelle righe.
#define SCOPE_ROT_BAND_H     8
#define SCOPE_ROT_BANDS      (SCOPE_DEC_H / SCOPE_ROT_BAND_H)   // 30
uint16_t scopeRotBandX(uint8_t band);
uint16_t scopeRotBandW(uint8_t band);

// EN: As scopeRotSlice (same pixels, same three faster paths), from the YCbCr frame into the output band
//     (scopeRotBandW(band) x 8 pixels, RGB565 in display order, black outside the visible circle).
// IT: Come scopeRotSlice (stessi pixel, stessi tre percorsi piu' veloci), dal frame YCbCr nella fascia di uscita
//     (scopeRotBandW(band) x 8 pixel, RGB565 in ordine del display, nero fuori dal cerchio visibile).
void scopeRotBandYcc(const uint8_t *const *mcu, uint8_t band, int32_t c, int32_t s, uint16_t *out);
