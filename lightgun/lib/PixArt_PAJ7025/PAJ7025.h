#ifdef PAJ7025_CAM

/**
 * @file PAJ7025.h
 * @brief Libreria C++ completa e unificata per il sensore IR PixArt PAJ7025R2 / PAJ7025R3 (Multiple Objects Tracking)
 * 
 * Questa libreria permette di interfacciarsi con il sensore di tracciamento IR PixArt tramite bus SPI.
 * È compatibile sia con ESP32 che con RP2040 (earlephilhower core).
 * 
 * ====================================================================================================
 * NOZIONI FONDAMENTALI SUL SENSORE PAJ7025
 * ====================================================================================================
 * 
 * 1) NUMERO DI PUNTI TRACCIATI (I 16 SLOT)
 * Il DSP della telecamera elabora l'immagine e ordina i punti IR in 16 "slot" di memoria fissi.
 * Non esiste un comando per forzare il sensore a cercare "solo 4 punti". Calcolerà sempre 16 oggetti.
 * Se ci sono meno di 16 oggetti nel campo visivo, gli slot inutilizzati avranno il parametro Area = 0.
 * IMPORTANTE: I punti non mantengono un "Tracking ID" persistente hardware. Se un oggetto scompare,
 * i punti successivi non traslano necessariamente in modo prevedibile. L'unico modo matematicamente
 * sicuro per trovare i punti attivi è scorrere TUTTI i 16 slot e scartare quelli con Area == 0.
 * 
 * 2) OTTIMIZZAZIONE DELLA VELOCITÀ DI LETTURA E FORMATI (FORMAT 1-4)
 * Dato che bisogna leggere 16 slot, il tempo perso sul bus SPI dipende dal Formato Dati scelto:
 * - FORMAT 1 (Completo): 16 bytes per oggetto. Totale 256 bytes per lettura. (Raccomandato se serve tutto).
 * - FORMAT 2 (Compatto): 6 bytes per oggetto. Totale 96 bytes per lettura. Contiene SOLO Area, X e Y.
 *   Usa il Formato 2 se hai bisogno di latenza zero su microcontrollori lenti e non ti servono Raggio o Bound.
 * - FORMAT 3 (Standard): 9 bytes per oggetto (144 bytes).
 * - FORMAT 4 (Avanzato): 13 bytes per oggetto (208 bytes).
 * Esiste un trucco per esperti: se ti interessano solo i primi 4 punti elaborati dal DSP (i più luminosi/grandi),
 * puoi interrompere la lettura SPI ("burst read") dopo aver letto solo i primi (4 * size) byte, risparmiando tempo.
 * 
 * 3) RISOLUZIONE, FRAMERATE ED ESPOSIZIONE (OVERCLOCKING)
 * - Risoluzione: Il sensore fisico (CMOS) è molto piccolo (~98x98 px), ma il DSP interpola a livello sub-pixel.
 *   Puoi impostare una risoluzione di uscita (es. 2940x2940) fino ad un massimo di 4095x4095.
 * - Frame Rate (FPS): Programmabile da 10 FPS a 200 FPS certificati (datasheet). Puoi forzare matematicamente
 *   valori come 209 FPS, mandando il sensore in lieve "overclock" (spesso retto benissimo dal chip).
 * - Esposizione: Se alzi l'FPS, il tempo fisico per scattare la "foto" si riduce! La formula vitale è:
 *   Esposizione Massima Consentita in uSec = (1 / FPS * 1.000.000) - 2700 uSec.
 *   Se superi questo limite (es. chiedendo esposizione di 3000 uSec a 200 FPS), il chip sballerà.
 * ====================================================================================================
 */

#ifndef PAJ7025_H
#define PAJ7025_H

#include <Arduino.h>
#include <SPI.h>

// Formati di lettura dei banchi di memoria interni. 
// Ogni formato restituisce dati diversi per ognuno dei 16 oggetti tracciati.
//
// FORMATO 1: 16 bytes per oggetto -> [Area 0-1] [X 2-3] [Y 4-5] [AvgBright 6] [MaxBright 7] [Radius/Range 8] [Bound L-R-U-D 9-12] [Aspect 13] [Vx 14] [Vy 15]
#define PAJ7025_FORMAT_1_256_BYTE 1

// FORMATO 2: 6 bytes per oggetto -> [Area 0-1] [X 2-3] [Y 4-5] (Solo dati base, ultra-veloce)
#define PAJ7025_FORMAT_2_96_BYTE  2

// FORMATO 3: 9 bytes per oggetto -> [Area 0-1] [X 2-3] [Y 4-5] [AvgBright 6] [MaxBright 7] [Radius/Range 8]
#define PAJ7025_FORMAT_3_144_BYTE 3

// FORMATO 4: 13 bytes per oggetto -> [Area 0-1] [X 2-3] [Y 4-5] [Bound L-R-U-D 6-9] [Aspect 10] [Vx 11] [Vy 12]
#define PAJ7025_FORMAT_4_208_BYTE 4

/**
 * @brief Struttura dati per rappresentare un singolo oggetto IR tracciato (uno dei 16 slot).
 * 
 * L'ordine dei byte nel chip è: [Area(2)], [CX(2)], [CY(2)], [AvgBright(1)], [MaxBright(1)], [Radius(1)], [Bounds(4)]...
 */
struct PAJ7025_Object {
    bool is_valid;              ///< true se l'oggetto esiste effettivamente in questo slot (Area > 0).
    uint16_t area;              ///< Dimensione del blob IR. (Formati 1, 2, 3, 4)
    uint16_t cx;                ///< Coordinata X interpolata del centroide (0 - Risoluzione Max impostata). (Formati 1, 2, 3, 4)
    uint16_t cy;                ///< Coordinata Y interpolata del centroide (0 - Risoluzione Max impostata). (Formati 1, 2, 3, 4)
    
    uint8_t average_brightness; ///< Luminosità media dell'intero blob (0-255). (Formati 1, 3)
    uint8_t max_brightness;     ///< Luminosità massima del pixel più brillante al centro (0-255). (Formati 1, 3)
    uint8_t range;              ///< Range stimato dell'oggetto. (Formati 1, 3)
    uint8_t radius;             ///< Raggio del blob IR. (Formati 1, 3)
    
    uint8_t boundary_left;      ///< Limite sinistro del bounding box. (Formati 1, 4)
    uint8_t boundary_right;     ///< Limite destro del bounding box. (Formati 1, 4)
    uint8_t boundary_up;        ///< Limite superiore del bounding box. (Formati 1, 4)
    uint8_t boundary_down;      ///< Limite inferiore del bounding box. (Formati 1, 4)
    uint8_t aspect_ratio;       ///< Rapporto d'aspetto larghezza/altezza dell'oggetto. (Formati 1, 4)
    uint8_t vx;                 ///< Velocità sull'asse X stimata dal chip. (Formati 1, 4)
    uint8_t vy;                 ///< Velocità sull'asse Y stimata dal chip. (Formati 1, 4)
};

class PAJ7025 {
  public:
    PAJ7025();

    /**
     * @brief Inizializza il sensore. Va chiamata nel setup().
     * @param spiPort Puntatore all'oggetto SPIClass (&SPI, &SPI1, ecc.). Lascia a te il pieno controllo sui pin hardware.
     * @param csPin Pin usato per il Chip Select.
     * @return true se la telecamera risponde col Model ID corretto (0x7025), false se scollegata/guasta.
     */
    bool begin(SPIClass* spiPort, uint8_t csPin);

    /**
     * @brief Imposta gli FPS (Frame per secondo).
     * @param fps Valori da 10 a 200 fps. Valori > 200 (es. 209) spingono il sensore in Overclocking.
     */
    void setFrameRate(int fps);

    /**
     * @brief Imposta il tempo di esposizione alla luce infrarossa.
     * @param exposure_uSec Tempo in microsecondi. Valore di default interno: ~1638 uSec.
     * @warning Se imposti l'FPS a 200, l'esposizione DEVE essere inferiore a ~2300 uSec.
     *          Formula vitale: (1 / FPS * 1.000.000) - 2700 uSec.
     */
    void setExposure(int exposure_uSec);

    /**
     * @brief Imposta il guadagno analogico (Gain) del sensore CMOS.
     * @param global Guadagno globale. Default = 15. Aumentalo per rilevare led deboli, abbassalo per ridurre il rumore.
     * @param ggh Guadagno extra addizionale (default 0).
     */
    void setGain(uint8_t global = 15, uint8_t ggh = 0);

    /**
     * @brief Definisce il tetto massimo della risoluzione in output. Il sensore scalerà i centroidi.
     * @param x_res Risoluzione asse X (Max 4095). Default 2940.
     * @param y_res Risoluzione asse Y (Max 4095). Default 2940.
     */
    void setResolution(int x_res = 2940, int y_res = 2940);

    /**
     * @brief Calibra i filtri del Digital Signal Processor (DSP) per riconoscere o scartare blob di luce.
     * @param area_min Dimensione minima del blob (pixel). Oggetti più piccoli vengono ignorati. Default = 3.
     * @param brightness_th Soglia di luminosità minima (0-255). Pixel più scuri di così non formano blob. Default = 151 (0x97).
     * @param area_max Dimensione massima del blob. Se supera questo valore (es. grossi riflessi del sole), lo scarta. Default = 9605 (0x2585).
     * @param noise_th Soglia di soppressione rumore di fondo. Default = 10 (0x0A).
     */
    void setDSP(uint8_t area_min = 3, uint8_t brightness_th = 0x97, uint16_t area_max = 9605, uint8_t noise_th = 10);

    /**
     * @brief Se true, il sensore calcola l'immagine sottraendo costantemente i fotogrammi senza illuminazione.
     * Utile per eliminare completamente le fonti di luce IR ambientali statiche (sole, lampadine) se usato con LED stroboscopico.
     */
    void enableFrameSubtraction(bool enable);

    /**
     * @brief Permette alla telecamera di pilotare un illuminatore LED (spesso pin G13 del modulo).
     * @param led_on Attiva fisicamente il pin GPIO dedicato ai LED.
     * @param led_frame_subtraction Se attivato, il LED lampeggerà a frame alterni sincronizzato con la Frame Subtraction interna.
     */
    void setLedGPIO(bool led_on, bool led_frame_subtraction);

    /**
     * @brief Abilita l'emissione del segnale di sync dell'esposizione (spesso pin G7). Utile per trigger hardware.
     */
    void setExposureGPIO(bool exposure_on);

    /**
     * @brief Funzione di correzione e calibrazione difetti fissi ottici (Fix Object Defect).
     */
    void setFOD(bool enable);

    /**
     * @brief Abilita pattern virtuali iniettati dal sensore, scavalcando la lente. Usato per debug SPI.
     * @param mode 0=OFF, 5=16 fissi, 6=4 fissi, 7=2 rotanti, 11=4 punti limite di un pixel.
     */
    void setDebugMode(uint8_t mode);

    /**
     * @brief Mette in standby totale il sensore (Power Down mode).
     */
    void powerDown();

    /**
     * @brief Azzera il sensore con un soft-reset dei registri.
     */
    void reset();
    
    /**
     * @brief Legge la firma hardware (Model ID).
     * @return Dovrebbe ritornare sempre 28709 (che in esadecimale è 0x7025).
     */
    uint16_t getModel();

    /**
     * @brief Legge direttamente dai registri quanto tempo passa per il calcolo di un frame.
     * @return Tempo in microsecondi (uSec).
     */
    uint32_t getFramePeriodMicroseconds();

    /**
     * @brief Lettura grezza dei byte dal banco di memoria specificato dal formato.
     * @param buffer Array raw. Deve essere capiente abbastanza (max 256 bytes).
     * @param format Uno dei PAJ7025_FORMAT_ (1, 2, 3 o 4).
     * @return Il numero effettivo di bytes letti.
     */
    int readDataRaw(uint8_t* buffer, uint8_t format = PAJ7025_FORMAT_1_256_BYTE);

    /**
     * @brief Funzione principe: legge, decodifica e compila gli struct PAJ7025_Object.
     * Attenzione: la PAJ7025 alloca sempre slot per 16 oggetti. Tu passa un array di 16 oggetti.
     * @param objects Array di minimo 16 elementi `PAJ7025_Object`.
     * @param format Quale banco leggere. Se la latenza SPI è critica, usa PAJ7025_FORMAT_2_96_BYTE per dimezzare il tempo di lettura.
     * @return Il numero di oggetti effettivamente avvistati con is_valid = true (da 0 a 16).
     */
    int readData(PAJ7025_Object* objects, uint8_t format = PAJ7025_FORMAT_1_256_BYTE);

  private:
    SPIClass* _spi;
    uint8_t _csPin;

    void select();
    void deselect();
    
    void writeRegister(uint8_t reg, uint8_t data);
    uint8_t readRegister(uint8_t reg);
    void burstRead(uint8_t reg_base, uint8_t* buffer, uint16_t num_bytes);
    
    void changeBank(uint8_t bank);
    void applyCommand(uint8_t bank);
    void loadInitialSettings();
    void parseObject(const uint8_t* rawData, PAJ7025_Object& obj, uint8_t format);
};

#endif // PAJ7025_H

#endif //PAJ7025_CAM