#pragma once

class LGFX : public lgfx::LGFX_Device
    {
      lgfx::Panel_ST7735S _panel_instance;
      lgfx::Bus_SPI       _bus_instance;
      lgfx::Light_PWM     _light_instance;

      public:
      LGFX(void)
      {
        {
          auto cfg = _bus_instance.config();
          cfg.spi_host   = SPI_HOST; // utilizzabili SPI2_HOST -> HSPI   SPI3_HOST -> VSPI
          cfg.spi_mode   = SPI_MODE;
          cfg.freq_write = SPI_FREQUENCY; //80000000;  // 50000000
          cfg.spi_3wire  = true; // non ha MISO (solo MOSI e SCK)
          cfg.use_lock   = false; // true se più dispositivi utilizzano il bus SPI - false se è l'unico sul bus SPI 
          cfg.dma_channel = SPI_DMA_CH_AUTO;
          cfg.pin_sclk   = TFT_SCLK;
          cfg.pin_miso   = TFT_MISO;
          cfg.pin_mosi   = TFT_MOSI;
          cfg.pin_dc     = TFT_DC;
          _bus_instance.config(cfg);
          _panel_instance.setBus(&_bus_instance);
        }
        {
          auto cfg = _panel_instance.config();
          cfg.pin_cs       = TFT_CS;
          cfg.pin_rst      = TFT_RST;
          cfg.panel_width  = TFT_WIDTH;
          cfg.panel_height = TFT_HEIGHT;
          cfg.offset_x     = TFT_OFFSET_X; //26
          cfg.offset_y     = TFT_OFFSET_Y; //1
          cfg.offset_rotation = 2;   // la rotazione del displai - 0 -> nessunoa, 1 90°, 2 -> 180 °, 3 -> 270°
          cfg.invert       = TFT_INVERT; //true;
          cfg.rgb_order    = TFT_RGB_ORDER; //false;
          _panel_instance.config(cfg);
        }
        {
          auto cfg = _light_instance.config();
          cfg.pin_bl = TFT_PIN_BL; // pin per retoilluminazione
          cfg.pwm_channel = TFT_PIN_BL_PWM_CHANNEL; // canale pwm se display supporta vari livelli di retroilluminazione
          _light_instance.config(cfg);
          _panel_instance.setLight(&_light_instance);
        }
      setPanel(&_panel_instance);
      }
    };
