#ifdef PAJ7025_CAM

/*!
 * @file PixArt_PAJ7025.cpp
 *
 * @copyright Alessandro Satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author Alessandro Satanassi
 * @version V1.0
 * @date 2026
 */

#include "PixArt_PAJ7025.h"

#ifdef ARDUINO_ARCH_ESP32
    #define delay(ms) vTaskDelay(pdMS_TO_TICKS(ms))
#endif


PAJ7025::PAJ7025() : _spi(nullptr), _csPin(0xFF), bank_global(0xFF) {
}

bool PAJ7025::begin(SPIClass* spiPort, uint8_t csPin, uint32_t clock) {
    if (spiPort == nullptr) return false;

    _spi = spiPort;
    _csPin = csPin;

    paj7025_spi_clock = clock;

    _spiSettings = SPISettings(paj7025_spi_clock, LSBFIRST, SPI_MODE3);

    pinMode(_csPin, OUTPUT);
    digitalWrite(_csPin, HIGH);


    uint16_t model;
    uint8_t i = 0;
    do {
        delay(10);

        model = getModel();
        i++;
    } while ((model != 28709) && (i < 10));

    if (model != 28709) {
        return false;
    }


    loadInitialSettings();

    changeBank(0x00);
    writeRegister(0x19, 4);

    return true;
}

void PAJ7025::select() {


    _spi->beginTransaction(_spiSettings);


    digitalWrite(_csPin, LOW);


}

void PAJ7025::deselect() {
    digitalWrite(_csPin, HIGH);

    _spi->endTransaction();
}

void PAJ7025::writeRegister(uint8_t reg, uint8_t data) {


    uint8_t buffer[3] = { 0x00, reg, data };

    select();
    _spi->transfer(buffer, 3);
    deselect();
}


const uint8_t* PAJ7025::burstRead(uint8_t reg_base, uint16_t num_bytes) {


    uint32_t total_length = static_cast<uint32_t>(num_bytes) + 2U;
    memset(_transferBuffer, 0x00, total_length);
    _transferBuffer[0] = 0x81;
    _transferBuffer[1] = reg_base;


    select();
    _spi->transfer(_transferBuffer, total_length);
    deselect();


    return &_transferBuffer[2];
}

void PAJ7025::universalReadRegister(uint8_t reg, uint8_t* buffer, uint16_t num_bytes) {

    if (buffer == nullptr) return;

    uint16_t total_length = num_bytes + 2;
    memset(_transferBuffer, 0x00, total_length);
    _transferBuffer[0] = (num_bytes > 1) ? 0x81 : 0x80;
    _transferBuffer[1] = reg;


    select();
    _spi->transfer(_transferBuffer, total_length);
    deselect();

    memcpy(buffer, &_transferBuffer[2], num_bytes);
}


void PAJ7025::changeBank(uint8_t bank) {

    if (bank == bank_global) return;
    writeRegister(0xEF, bank);
    bank_global = bank;
}

void PAJ7025::applyCommand(uint8_t bank) {


    changeBank(bank);
    writeRegister(0x01, 0x01);
}

void PAJ7025::loadInitialSettings() {


    writeRegister(0xEF, 0x00);
    writeRegister(0xDC, 0x00);
    writeRegister(0xFB, 0x04);
    writeRegister(0xEF, 0x00);
    writeRegister(0x2F, 0x05);
    writeRegister(0x30, 0x00);
    writeRegister(0x30, 0x01);
    writeRegister(0x1F, 0x00);


    writeRegister(0xEF, 0x01);
    writeRegister(0x2D, 0x00);


    writeRegister(0xEF, 0x0C);

    writeRegister(0x64, 0x00);
    writeRegister(0x65, 0x00);
    writeRegister(0x66, 0x00);
    writeRegister(0x67, 0x00);
    writeRegister(0x68, 0x00);
    writeRegister(0x69, 0x00);


    writeRegister(0x6A, 0x00);
    writeRegister(0x6B, 0x00);

    writeRegister(0x6C, 0x00);
    writeRegister(0x71, 0x00);
    writeRegister(0x72, 0x00);
    writeRegister(0x12, 0x00);
    writeRegister(0x13, 0x00);

    writeRegister(0xEF, 0x00);
    writeRegister(0x01, 0x01);
}


void PAJ7025::setFrameRate(uint16_t fps) {


    uint32_t frame_long = uint32_t{10000000} / fps;

    changeBank(0x0C);


    universalWriteRegister(0x07, frame_long & 0xFF, (frame_long >> 8) & 0xFF, (frame_long >> 16) & 0xFF);
}

void PAJ7025::setExposure(uint16_t exposure_uSec) {


    uint32_t exposure_time = static_cast<uint32_t>(exposure_uSec) * 5U;

    changeBank(0x0C);


    universalWriteRegister(0x0F, exposure_time & 0xFF, (exposure_time >> 8) & 0xFF);
    applyCommand(0x01);
}

void PAJ7025::setGain(uint8_t global, uint8_t ggh) {
    changeBank(0x0C);


    universalWriteRegister(0x0B, global, ggh);
    applyCommand(0x01);
}

void PAJ7025::setResolution(uint16_t x_res, uint16_t y_res) {
    changeBank(0x0C);


    universalWriteRegister(0x60, x_res & 0xFF, (x_res >> 8) & 0xFF, y_res & 0xFF, (y_res >> 8) & 0xFF);
}

void PAJ7025::setDSP(uint8_t area_min, uint8_t brightness_th, uint16_t area_max, uint8_t noise_th) {
    changeBank(0x0C);
    writeRegister(0x46, area_min);
    writeRegister(0x47, brightness_th);

    changeBank(0x00);


    universalWriteRegister(0x0B, area_max & 0xFF, (area_max >> 8) & 0xFF);
    writeRegister(0x0F, noise_th);

}


uint16_t PAJ7025::getModel() {
    changeBank(0x00);


    uint16_t model = 0;
    universalReadRegister(0x02, &model);
    return model;
}

static inline uint8_t PAJ7025_getFormatInfo(uint8_t format, uint8_t& register_bank) {
    if (format == PAJ7025_FORMAT_BASIC) {
        register_bank = 0x09;
        return 6;
    }
    register_bank = 0x05;
    return 16;
}

int32_t PAJ7025::readData(PAJ7025_Object* objects, uint8_t format) {
    uint8_t register_bank;
    const uint8_t obj_size = PAJ7025_getFormatInfo(format, register_bank);
    const uint16_t bytes_to_read = 4U * obj_size;

    changeBank(register_bank);
    const uint8_t* buffer = burstRead(0x00, bytes_to_read);

    int32_t valid_objects = 0;
    for (int32_t i = 0; i < 4; i++) {

        parseObject(&buffer[i * obj_size], objects[i], format);
        if (objects[i].is_valid) {
            valid_objects++;
        }
    }

    return valid_objects;
}

void PAJ7025::parseObject(const uint8_t* data, PAJ7025_Object& obj, uint8_t format) {

    memset(&obj, 0, sizeof(PAJ7025_Object));

    obj.area = data[0] | ((data[1] & 0x3f) << 8);
    obj.cx   = data[2] | ((data[3] & 0x0f) << 8);
    obj.cy   = data[4] | ((data[5] & 0x0f) << 8);

    obj.is_valid = (obj.area > 0);

    if (format == PAJ7025_FORMAT_EXTENDED) {
        obj.average_brightness = data[6];
        obj.max_brightness     = data[7];
        obj.range              = data[8] >> 4;
        obj.radius             = data[8] & 0x0f;

        obj.boundary_left      = data[9] & 0x7f;
        obj.boundary_right     = data[10] & 0x7f;
        obj.boundary_up        = data[11] & 0x7f;
        obj.boundary_down      = data[12] & 0x7f;
        obj.aspect_ratio       = data[13];
        obj.vx                 = data[14];
        obj.vy                 = data[15];
    }
}

#endif