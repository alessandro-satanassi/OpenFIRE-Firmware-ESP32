#ifdef USE_MULTI_ONE_EURO_FILTER
/*!
 * @file OpenFIRE_Multi_One_Euro_Filter.h
 * @brief Library for One Euro Filter for 4 LED
 * @n CPP Library for One Euro Filter for 4 LED
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V2.0
 * @date 2026
 */

#ifndef OpenFIRE_Multi_One_Euro_Filter_h
#define OpenFIRE_Multi_One_Euro_Filter_h

#include <Arduino.h>
#include "OpenFIREConst.h"
#include <OpenFIRECameraProfile.h>

class OpenFIRE_One_Euro_Multi {
private:
    // Struttura fusa: mantiene sia lo storico posizionale (prev/hat) che vettoriale (vel).
    // Usare una singola struct garantisce che la cache L1 della CPU legga tutto lo stato
    // di un singolo asse in una sola operazione (Locality of Reference).
    struct FilterState {
        float x_prev, x_hat;
        float y_prev, y_hat;
        float vel_x_hat, vel_y_hat;
    };

    FilterState states[4];
    unsigned long lastMicros;
    bool initialized = false;

    // ==========================================
    // --- PARAMETRI DI TUNING E-SPORTS (BILANCIAMENTO DEFINITIVO) ---

    // min_cutoff: La "lentezza" del mirino quando ti muovi pochissimo o sei fermo.
    // Impostato a 0.1f: il punto di equilibrio perfetto. 0.2f era leggermente
    // scivoloso, 0.05f era troppo rigido. 0.1f garantisce mira da cecchino solida.
    static constexpr float min_cutoff = 0.1f;

    // --- GESTIONE ASIMMETRICA DELLA VELOCITÀ ---
    // d_cutoff_base: Reattività per i movimenti di precisione.
    // Il tuning temporale resta invariato e viene gestito tramite il dt reale.
    static constexpr float d_cutoff_base = 1.0f;

    // d_cutoff_snap: Reattività per scatti violenti e frenate brusche.
    static constexpr float d_cutoff_snap = 25.0f;

    float sensor_noise_unit_x = 1.0f;
    float sensor_noise_unit_y = 1.0f;

    // Coefficienti di tuning normalizzati sul pixel fisico del sensore.
    // Sono ricavati dal comportamento DFRobot originale, ma non contengono
    // risoluzioni o dimensioni hardcoded di una specifica CAM.
    // Con DFRobot: sensor_noise_unit X/Y = 32, quindi si ottengono
    // esattamente snap=1000, edge=2000 e beta=0.0077 come prima.
    static constexpr float snap_sensor_gain = 31.25f;
    static constexpr float snap_edge_sensor_gain = 62.5f;
    static constexpr float beta_sensor_gain = 0.2464f;

    // Soglie di snap adattate alla granularita fisica/rumore della CAM.
    float snap_base_x = 0.0f;
    float snap_base_y = 0.0f;

    float snap_edge_multiplier_x = 0.0f;
    float snap_edge_multiplier_y = 0.0f;

    // max_cutoff: Il limite di banda passante superiore.
    static constexpr float max_cutoff = 30.0f;

    // Beta adattivo normalizzato sulla granularita fisica del sensore.
    // Una CAM con pixel fisici relativamente piu grossi (o CamNoiseFactor > 1)
    // parte piu conservativa; una CAM piu pulita puo essere resa piu reattiva
    // diminuendo CamNoiseFactor.
    float beta_base_x = 0.0f;
    float beta_base_y = 0.0f;

    // MICRO-SNAP: resta volutamente invariato. Non e una soglia di rumore ottico:
    // elimina soltanto la coda frazionaria quando l'input intero e gia fermo.
    static constexpr float micro_snap_threshold = 0.5f;

    // ==========================================


    float inv_center_x = 0.0f;
    float inv_center_y = 0.0f;
    int mouseMaxX = 0;
    int mouseMaxY = 0;
    float fallbackDt = 1.0f / 209.0f;

    static constexpr float OEF_TWO_PI = 6.28318530718f;

    // Sostituisce l'implementazione classica dell'Exponential Moving Average.
    // Espande la formula matematica per evitare chiamate di funzione ricorsive.
    static inline float fast_alpha(float cutoff, float dt_two_pi) {
        float te = cutoff * dt_two_pi;
        return te / (te + 1.0f);
    }

public:
    OpenFIRE_One_Euro_Multi();
    void configure(const CameraProfile& profile);
    
    // Il passaggio tramite puntatori permette di leggere i raw (int) e restituire le coordinate
    // sub-pixel perfette (float) senza copie di array.
    void process(int* x_in, int* y_in, float* x_out, float* y_out);
};

#endif // OpenFIRE_Multi_One_Euro_Filter_h

#endif // USE_MULTI_ONE_EURO_FILTER