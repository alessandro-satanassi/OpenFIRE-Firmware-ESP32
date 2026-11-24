/*!
 * @file LGFX_128_GC9A01_240x240.hpp
 * @brief Library for LGFX_128_GC9A01_240x240
 * @n CPP Library for LGFX_128_GC9A01_240x240
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#pragma once

// ===================================================================================
// DEFINIZIONE DEVICE: APPROCCIO MACRO-DRIVEN (SPI)
// ===================================================================================
// Per garantire prestazioni estreme (es. 80MHz) necessarie per animazioni a colori,
// questa classe usa un'architettura statica. I parametri non vengono passati a 
// runtime, ma risolti dal preprocessore C++ tramite le macro del platformio.ini.
// Questo azzera l'overhead di memoria e rende il binario compilato strettamente
// accoppiato alla board target (es. SCOPE su Waveshare ESP32-C5-Zero).
class LGFX : public lgfx::LGFX_Device
    {
      lgfx::Panel_GC9A01  _panel_instance;   // display rotondo 1.28" 240x240
      lgfx::Bus_SPI       _bus_instance;
      lgfx::Light_PWM     _light_instance;

      public:
      LGFX(void)
      {
        // ===================================================================================
        // HARDWARE LAYER: BUS SPI & DMA
        // ===================================================================================
        {
          auto cfg = _bus_instance.config();
          cfg.spi_host   = SPI_HOST; // sull'ESP32-C5 l'unico SPI libero e' SPI2_HOST
          cfg.spi_mode   = SPI_MODE;
          cfg.freq_write = SPI_FREQUENCY; //80000000;  // 50000000
          
          // Ottimizzazione segnale: Disattivando il MISO (poiché il display non invia 
          // dati rilevanti indietro) si risparmia un pin GPIO e si semplifica il routing.
          cfg.spi_3wire  = true; // non ha MISO (solo MOSI e SCK)
          
          cfg.use_lock   = false; // true se più dispositivi utilizzano il bus SPI - false se è l'unico sul bus SPI 
          
          // Delegando la scelta del canale DMA (Direct Memory Access) in AUTO, si evita
          // il rischio di "core panic" sull'ESP32, permettendo al framework di allocare
          // la risorsa hardware libera, garantendo refresh rate elevatissimi.
          cfg.dma_channel = SPI_DMA_CH_AUTO;
          
          cfg.pin_sclk   = TFT_SCLK;
          cfg.pin_miso   = TFT_MISO;
          cfg.pin_mosi   = TFT_MOSI;
          cfg.pin_dc     = TFT_DC;
          _bus_instance.config(cfg);
          _panel_instance.setBus(&_bus_instance);
        }
        
        // ===================================================================================
        // PRESENTATION LAYER: PANNELLO TFT E VARIANTI CLONI
        // ===================================================================================
        {
          auto cfg = _panel_instance.config();
          cfg.pin_cs       = TFT_CS;
          cfg.pin_rst      = TFT_RST;
          cfg.panel_width  = TFT_WIDTH;
          cfg.panel_height = TFT_HEIGHT;
          
          // GC9A01: la memoria del controller e' 240x240 come il pannello, quindi
          // nessun offset. I moduli GC9A01 vogliono i colori invertiti (TFT_INVERT=1).
          cfg.offset_x     = TFT_OFFSET_X; //0
          cfg.offset_y     = TFT_OFFSET_Y; //0
          cfg.offset_rotation = 0;   // la rotazione del display - 0 -> nessuna, 1 90°, 2 -> 180 °, 3 -> 270° (da regolare a montaggio fatto)
          cfg.invert       = TFT_INVERT; //true;
          cfg.rgb_order    = TFT_RGB_ORDER; //false;
          cfg.readable     = false;  // senza MISO non si legge dal display
          cfg.bus_shared   = false;  // il display e' l'unico dispositivo sul bus SPI
          _panel_instance.config(cfg);
        }
        
        // ===================================================================================
        // PERIPHERAL LAYER: GESTIONE RETROILLUMINAZIONE PWM
        // ===================================================================================
        {
          auto cfg = _light_instance.config();
          cfg.pin_bl = TFT_PIN_BL; // pin per retoilluminazione
          
          // L'assegnazione flessibile del canale PWM previene conflitti hardware 
          // (es. LEDC timer clash) con altri componenti del sistema (es. solenoidi/audio).
          cfg.pwm_channel = TFT_PIN_BL_PWM_CHANNEL; // canale pwm se display supporta vari livelli di retroilluminazione
          _light_instance.config(cfg);
          _panel_instance.setLight(&_light_instance);
        }
      setPanel(&_panel_instance);
      }
    };