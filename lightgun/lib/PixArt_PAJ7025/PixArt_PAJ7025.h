#ifdef PAJ7025_CAM

/*!
 * @file PixArt_PAJ7025.h
 *
 * @copyright Alessandro Satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author Alessandro Satanassi
 * @version V1.0
 * @date 2026
 */


#ifndef PixArt_PAJ7025_H
#define PixArt_PAJ7025_H

#include <Arduino.h>
#include <SPI.h>


#define PAJ7025_FORMAT_EXTENDED 1


#define PAJ7025_FORMAT_BASIC     2


struct PAJ7025_Object {
    bool is_valid;
    uint16_t area;
    uint16_t cx;
    uint16_t cy;

    uint8_t average_brightness;
    uint8_t max_brightness;
    uint8_t range;
    uint8_t radius;

    uint8_t boundary_left;
    uint8_t boundary_right;
    uint8_t boundary_up;
    uint8_t boundary_down;
    uint8_t aspect_ratio;
    uint8_t vx;
    uint8_t vy;
};

class PAJ7025 {
  public:
    PAJ7025();


    bool begin(SPIClass* spiPort, uint8_t csPin, uint32_t clock);


    void setFrameRate(uint16_t fps);


    void setExposure(uint16_t exposure_uSec);


    void setGain(uint8_t global, uint8_t ggh);


    void setResolution(uint16_t x_res, uint16_t y_res);


    void setDSP(uint8_t area_min, uint8_t brightness_th, uint16_t area_max, uint8_t noise_th);


    uint16_t getModel();


    int32_t readData(PAJ7025_Object* objects, uint8_t format);

  private:
    SPIClass* _spi;
    uint8_t _csPin;
    SPISettings _spiSettings;
    uint8_t bank_global = 0xFF;
    uint32_t paj7025_spi_clock = 2000000 ;
    uint8_t _transferBuffer[258];

    void select();
    void deselect();

    void writeRegister(uint8_t reg, uint8_t data);

    template <typename... Args>
    void universalWriteRegister(uint8_t reg, Args... args) {

      constexpr uint8_t length = sizeof...(args);


      if (length == 0) return;


      constexpr uint8_t codice = (length > 1) ? 0x01 : 0x00;


      uint8_t buffer[length + 2] = { codice, reg, static_cast<uint8_t>(args)... };


      select();
      _spi->transfer(buffer, length + 2);
      deselect();
    }

    const uint8_t* burstRead(uint8_t reg_base, uint16_t num_bytes);


    void universalReadRegister(uint8_t reg, uint8_t* buffer, uint16_t num_bytes);


    inline void universalReadRegister(uint8_t reg, uint16_t* dest, uint16_t num_byte = 2) { universalReadRegister(reg, (uint8_t*)dest, num_byte > 2 ? 2 : num_byte); }


    void changeBank(uint8_t bank);
    void applyCommand(uint8_t bank);
    void loadInitialSettings();
    void parseObject(const uint8_t* rawData, PAJ7025_Object& obj, uint8_t format);
};

#endif

#endif