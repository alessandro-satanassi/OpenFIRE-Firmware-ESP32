#ifdef PAJ7025_CAM


#include "PAJ7025.h"

PAJ7025::PAJ7025() : _spi(nullptr), _csPin(0) {
}

bool PAJ7025::begin(SPIClass* spiPort, uint8_t csPin) {
    if (spiPort == nullptr) return false;
    
    _spi = spiPort;
    _csPin = csPin;
    
    if (_csPin != 255) {
        pinMode(_csPin, OUTPUT);
    }
    deselect();
    
    // Configura i settaggi iniziali di PixArt
    select();
    loadInitialSettings();
    deselect();
    
    // Controlla che il sensore risponda correttamente tramite lettura del product ID
    uint16_t model = getModel();
    if (model != 28709) { // 28709 in esadecimale è 0x7025 (ovvero PAJ 7025)
        return false;
    }
    
    return true;
}

void PAJ7025::select() {
    // Inizializza SPI transaction con la velocità raccomandata (max 14 MHz per PAJ7025)
    // Modalità SPI: LSB First, Mode 3 (CPOL=1, CPHA=1)
    _spi->beginTransaction(SPISettings(8000000, LSBFIRST, SPI_MODE3));
    if (_csPin != 255) {
        digitalWrite(_csPin, LOW);
    }
}

void PAJ7025::deselect() {
    if (_csPin != 255) {
        digitalWrite(_csPin, HIGH);
    }
    _spi->endTransaction();
}

void PAJ7025::writeRegister(uint8_t reg, uint8_t data) {
    // Write format: Bit 7 = 0 (Write), Bit 6-0 = 0 (Single byte) => 0x00
    _spi->transfer(0x00); 
    _spi->transfer(reg);
    _spi->transfer(data);
}

uint8_t PAJ7025::readRegister(uint8_t reg) {
    // Read format: Bit 7 = 1 (Read), Bit 6-0 = 0 (Single byte) => 0x80
    _spi->transfer(0x80); 
    _spi->transfer(reg);
    return _spi->transfer(0x00); // Send dummy byte to clock data in
}

void PAJ7025::burstRead(uint8_t reg_base, uint8_t* buffer, uint16_t num_bytes) {
    // Burst Read format: Bit 7 = 1, Bit 6-0 = 1 => 0x81 (Sequential read)
    _spi->transfer(0x81); 
    _spi->transfer(reg_base); // Usually starts at register 0x00 in the target bank
    for (uint16_t i = 0; i < num_bytes; i++) {
        buffer[i] = _spi->transfer(0x00);
    }
}

void PAJ7025::changeBank(uint8_t bank) {
    // Il Bank register è 0xEF
    writeRegister(0xEF, bank);
}

void PAJ7025::applyCommand(uint8_t bank) {
    // Il cambio di alcuni registri richiede di lanciare un "Apply Command" 
    // scrivendo 0x01 nel registro 0x01 del bank target.
    changeBank(bank);
    writeRegister(0x01, 0x01);
}

void PAJ7025::loadInitialSettings() {
    // Questa sequenza "magica" abilita l'elaborazione interna (DSP) 
    // della telecamera PixArt e configura i pin di output.
    // Estratta direttamente dal datasheet del PAJ7025R2 (Pagina 16-17).
    
    changeBank(0x00);
    writeRegister(0xDC, 0x00); // disable internal system control
    writeRegister(0xFB, 0x04); // LEDDAC disable [2]
    writeRegister(0x2F, 0x05); // Sensor ON (wake up from power down)
    writeRegister(0x30, 0x00);
    writeRegister(0x30, 0x01); // Power control updated
    writeRegister(0x1F, 0x00); // freerun_irtx_disable
    
    applyCommand(0x01);
    
    changeBank(0x01);
    writeRegister(0x2D, 0x00); // Vertical Flip OFF (0x00). Scrivere 0x01 per ribaltare Y.
    
    changeBank(0x0C);
    // Azzera impostazioni di mode
    writeRegister(0x64, 0x00); 
    writeRegister(0x65, 0x00);
    writeRegister(0x66, 0x00);
    writeRegister(0x67, 0x00);
    writeRegister(0x68, 0x00);
    writeRegister(0x69, 0x00);
    
    // Default GPIO e modalità opzionali
    setFOD(false);
    setExposureGPIO(true);
    
    writeRegister(0x6C, 0x00);
    writeRegister(0x72, 0x00);
    writeRegister(0x12, 0x00); // disable keyscan
    writeRegister(0x13, 0x00); // disable keyscan
    
    applyCommand(0x00);
}

void PAJ7025::powerDown() {
    select();
    changeBank(0x00);
    writeRegister(0x2F, 0x00); // Sensor power down
    writeRegister(0x30, 0x00);
    writeRegister(0x31, 0x01);
    applyCommand(0x00);
    deselect();
}

void PAJ7025::reset() {
    select();
    changeBank(0x00);
    writeRegister(0x64, 0x00); // Reset soft
    deselect();
}

void PAJ7025::setFrameRate(int fps) {
    if(fps <= 0) return;
    
    // Calcola il periodo in nanosecondi e scala per le unità del registro
    // Formula dal codice driver originale: frame_time = 1/fps. frame_long = frame_time * 10.000.000
    float frame_time = 1.0f / fps;
    uint32_t frame_long = (uint32_t)(frame_time * 10000000);

    select();
    changeBank(0x0C);
    writeRegister(0x07, frame_long & 0xFF);         // frame period Low Byte
    writeRegister(0x08, (frame_long >> 8) & 0xFF);  // frame period High Byte
    writeRegister(0x09, (frame_long >> 16) & 0xFF); // frame period Highest Byte
    deselect();
}

void PAJ7025::setExposure(int exposure_uSec) {
    // Assicurati di non chiamare questa funzione con un valore in uSec superiore al frame_rate.
    // L'esposizione massima ammessa è ~ (1 / fps * 1000000) - 2700 uSec.
    
    // Conversione: 1 unita di registro = 0.2 uSec (200 nSec)
    int exposure_time = (int)(exposure_uSec / 0.2f);
    
    select();
    changeBank(0x0C);
    writeRegister(0x0F, exposure_time & 0xFF);          // exposure LB
    writeRegister(0x10, (exposure_time >> 8) & 0xFF);   // exposure HB
    applyCommand(0x01);
    deselect();
}

void PAJ7025::setGain(uint8_t global, uint8_t ggh) {
    select();
    changeBank(0x0C);
    writeRegister(0x0B, global); // B_global (0-15 di default, può spingersi oltre)
    writeRegister(0x0C, ggh);    // B_ggh
    applyCommand(0x01);
    applyCommand(0x00);
    deselect();
}

void PAJ7025::setResolution(int x_res, int y_res) {
    select();
    changeBank(0x0C);
    writeRegister(0x60, x_res & 0xFF);          // x res LB
    writeRegister(0x61, (x_res >> 8) & 0xFF);   // x res HB
    writeRegister(0x62, y_res & 0xFF);          // y res LB
    writeRegister(0x63, (y_res >> 8) & 0xFF);   // y res HB
    applyCommand(0x01);
    deselect();
}

void PAJ7025::setDSP(uint8_t area_min, uint8_t brightness_th, uint16_t area_max, uint8_t noise_th) {
    select();
    changeBank(0x0C);
    writeRegister(0x46, area_min);        // Soglia dimensione area minima per accettare l'oggetto
    writeRegister(0x47, brightness_th);   // Soglia di luminosità
    
    changeBank(0x00);
    writeRegister(0x0B, area_max & 0xFF);          // max area threshold Low Byte
    writeRegister(0x0C, (area_max >> 8) & 0xFF);   // max area threshold High Byte
    writeRegister(0x0F, noise_th);                 // threshold del rumore
    writeRegister(0x2B, 0x01);                     // configurazione interna DSP PixArt
    applyCommand(0x00);
    deselect();
}

void PAJ7025::enableFrameSubtraction(bool enable) {
    select();
    changeBank(0x00);
    writeRegister(0x28, enable ? 0x01 : 0x00); // 0x28 -> Frame subtraction enable/disable
    applyCommand(0x00);
    deselect();
}

void PAJ7025::setLedGPIO(bool led_on, bool led_frame_subtraction) {
    select();
    changeBank(0x0C);
    // 0x71 gestisce G13 / LED_SIDE. Scrive 0x08 per attivarlo
    writeRegister(0x71, led_on ? 0x08 : 0x00);
    
    changeBank(0x00);
    // Se la sottrazione del frame è on, salta frame per leggere il background senza LED
    writeRegister(0x4F, led_frame_subtraction ? 0x5C : 0x2C); 
    applyCommand(0x00);
    deselect();
}

void PAJ7025::setExposureGPIO(bool exposure_on) {
    select();
    changeBank(0x0C);
    writeRegister(0x6B, exposure_on ? 0x08 : 0x00); // 0x6B gestisce Exposure GPIO G7
    applyCommand(0x00);
    deselect();
}

void PAJ7025::setFOD(bool enable) {
    select();
    changeBank(0x0C);
    writeRegister(0x6A, enable ? 0x07 : 0x00); // G6 GPIO settings
    applyCommand(0x00);
    deselect();
}

void PAJ7025::setDebugMode(uint8_t mode) {
    select();
    changeBank(0x01);
    writeRegister(0x2B, mode); 
    // Modi noti:
    // 0 = OFF (Immagine reale dal sensore)
    // 5 = Pattern di 16 oggetti fissi
    // 6 = Pattern di 4 oggetti fissi
    // 7 = Pattern di 2 oggetti mobili in cerchio
    // 11 = Pattern di 4 oggetti limiti a singolo pixel
    deselect();
}

uint16_t PAJ7025::getModel() {
    select();
    changeBank(0x00);
    uint8_t model_lb = readRegister(0x02);
    uint8_t model_hb = readRegister(0x03);
    deselect();
    return (model_hb << 8) | model_lb;
}

uint32_t PAJ7025::getFramePeriodMicroseconds() {
    select();
    changeBank(0x0C);
    uint32_t cmd_frame_period = readRegister(0x07);
    cmd_frame_period |= readRegister(0x08) << 8;
    cmd_frame_period |= readRegister(0x09) << 16;
    deselect();

    // Lettura è in unità di 100ns -> convertiamo in us, arrotondando
    uint32_t microseconds = (cmd_frame_period / 10) + ((cmd_frame_period % 10) >= 5 ? 1 : 0);
    return microseconds;
}

int PAJ7025::readDataRaw(uint8_t* buffer, uint8_t format) {
    int num_bytes = 256;
    uint8_t format_code = 5;
    
    // Mappatura tra formato utente (1..4) e codice banco di lettura (5, 9, 10, 11)
    switch (format) {
        case PAJ7025_FORMAT_1_256_BYTE: num_bytes = 256; format_code = 5;  break;
        case PAJ7025_FORMAT_2_96_BYTE:  num_bytes = 96;  format_code = 9;  break;
        case PAJ7025_FORMAT_3_144_BYTE: num_bytes = 144; format_code = 10; break;
        case PAJ7025_FORMAT_4_208_BYTE: num_bytes = 208; format_code = 11; break;
        default: break;
    }

    select();
    writeRegister(0xEF, format_code); // Seleziona il bank corretto per il burst
    burstRead(0x00, buffer, num_bytes); // I dati dei blob partono sempre dal registro 0x00
    deselect();
    
    return num_bytes;
}

int PAJ7025::readData(PAJ7025_Object* objects, uint8_t format) {
    uint8_t buffer[256];
    int bytes_read = readDataRaw(buffer, format);
    
    int obj_size = bytes_read / 16;
    int valid_objects = 0;
    
    for (int i = 0; i < 16; i++) {
        // Estrapola il singolo oggetto in base al formato
        parseObject(&buffer[i * obj_size], objects[i], format);
        if (objects[i].is_valid) {
            valid_objects++;
        }
    }
    
    return valid_objects; // Ritorna quanti blob risultano attivi in questo momento
}

void PAJ7025::parseObject(const uint8_t* data, PAJ7025_Object& obj, uint8_t format) {
    // Azzera i dati per prevenire valori "sporchi" di passate letture
    memset(&obj, 0, sizeof(PAJ7025_Object));
    
    // Tutti i formati (1, 2, 3, 4) condividono i primi 6 byte per Area e Centroide
    // Byte 0-1: Area
    obj.area = data[0] | ((data[1] & 0x3f) << 8);
    // Byte 2-3: X Centroide (Coordinate Native 0-4095)
    obj.cx = data[2] | ((data[3] & 0x0f) << 8);
    // Byte 4-5: Y Centroide (Coordinate Native 0-4095)
    obj.cy = data[4] | ((data[5] & 0x0f) << 8);
    
    // Determina la validità dell'oggetto: di base se l'area è > 0
    obj.is_valid = (obj.area > 0);

    // I formati 1 e 3 riportano luminosità e raggio
    if (format == PAJ7025_FORMAT_1_256_BYTE || format == PAJ7025_FORMAT_3_144_BYTE) {
        obj.average_brightness = data[6];
        obj.max_brightness     = data[7];
        obj.range              = data[8] >> 4;
        obj.radius             = data[8] & 0x0f;
        
        if (!obj.is_valid && obj.max_brightness > 0) {
            obj.is_valid = true; // Safe fallback: se c'è luminosità l'oggetto c'è
        }
    }

    // I formati 1 e 4 riportano Bound, Aspect Ratio e Velocità
    if (format == PAJ7025_FORMAT_1_256_BYTE || format == PAJ7025_FORMAT_4_208_BYTE) {
        // Nel formato 4 i bound partono dal byte 6 (dato che mancano luminosità/raggio)
        int offset = (format == PAJ7025_FORMAT_4_208_BYTE) ? 3 : 0;
        
        obj.boundary_left  = data[9 - offset] & 0x7f;
        obj.boundary_right = data[10 - offset] & 0x7f;
        obj.boundary_up    = data[11 - offset] & 0x7f;
        obj.boundary_down  = data[12 - offset] & 0x7f;
        obj.aspect_ratio   = data[13 - offset];
        obj.vx             = data[14 - offset];
        obj.vy             = data[15 - offset];
    }
}

#endif // PAJ7025_CAM