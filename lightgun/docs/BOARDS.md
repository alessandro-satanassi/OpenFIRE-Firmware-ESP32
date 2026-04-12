# Supported Boards Layouts:
 - [Waveshare esp32-s3-pico](#waveshare-esp32-s3-pico)
 - [ESP32_S3_WROOM1_DevKitC_1_N16R8](#esp32_s3_wroom1_devkitc_1_n16r8)
 - [Raspberry Pi Pico (Non/W)](#raspberry-pi-pico-nonw)
 - [Adafruit ItsyBitsy RP2040](#adafruit-itsybitsy-rp2040)
 - [Adafruit Keeboar KB2040](#adafruit-keeboar-kb2040)
 - [Arduino Nano RP2040 Connect](#arduino-nano-rp2040-connect)
 - [Waveshare RP2040 Zero](#waveshare-rp2040-zero)
 - [VCC-GND YD RP2040 (Not Implemented)](#vcc-gnd-yd-rp2040)

```
                                        Symbol Legend:
                       (x) = GND/No Connect | (-) = GPIO | (p) = Power
```

> [!NOTE]
> Any button/function GPIO pin assignments written here are the ***defaults***, and can be changed for any other function from within the OpenFIRE App.
> 
> All of these layouts can also be previewed from within the OpenFIRE App (v3.0+) under *Help->View Compatible Boards*.

### 

## waveshare esp32-s3-pico
```
                                        (_____)
                      A Button   11  |-) *USB* (-| VBUS  5v (USB voltage)
                      B Button   12  |-)       (-| VSYS  (Input from Battery/Output to NeoPixels)
                                 GND |x)       (x| GND
                      C Button   13  |-)       (x| 3V3   En
                         Start   14  |-)       (p| 3V3   Out (to Display/Cam/Analog Inputs)
                        Select   15  |-)       (x|  10   
                   Home Button   16  |-)       (-|  09   Temp Sensor
                                 GND |x)       (x| GND 
                      D-Pad Up   17  |-)       (-|  08   *Unmapped* Analop X
                    D-Pad Down   18  |-)       (-|  07   *Unmapped* Analog Y
                    D-Pad Left   33  |-)       (x| RUN
                   D-Pad Right   34  |-)       (-|  06   *Unmapped*
                                 GND |x)       (x| GND
Peripherals SDA (OLED) RGB Red   35  |-)       (-|  05   Camera SCL
Peripherals SCL OLED RGB Green   36  |-)       (-|  04   Camera SDA
                      RGB Blue   37  |-)       (-|  02   *Unmapped* Peripherals SCL (OLED)
                   Pump Action   38  |-)       (-|  01   *Unmapped* Peripherals SDA (OLED)
                                 GND |x)       (x| GND
                         Pedal   39  |-)       (-|  41   Rumble Signal
                       Trigger   40  |-) _|_|_ (-|  42   Solenoid Signal

* GPIO 21 at RGB LED 

```

## ESP32_S3_WROOM1_DevKitC_1_N16R8
```

                                         (__)  (__)  
                                 GND |-)  USB  USB  (-| GND
                                 GND |-)  SER  OTG  (-| 5V IN-OUT
    COLLEGATO A USB OTG ???  GPIO 19 |x)            (x| GPIO 14
    COLLEGATO A USB OTG ???  GPIO 20 |-)            (x| GPIO 13  MISO -> RX  (non usato per momento) SPI
     2DOWN     GUNC          GPIO 21 |-)            (p| GPIO 12  SCK  -> SCK (non usato per momento) SPI
     2LEFT      PEDAL        GPIO 47 |-)            (x| GPIO 11  MOSI -> TX  (non usato per momento) SPI
       COLLEGATO AL LED RGB  GPIO 48 |-)            (-| GPIO 10  SS   -> CSn (non usato per momento) SPI
     2RIGHT     PUMP         GPIO 45 |x)            (x| GPIO  9  SCL (sono di default li riusiamo in OpenFIRE) CAM SCL
  COLLEGATO A PULSANTE BOOT  GPIO  0 |-)            (-| GPIO 46
    2MIDDLE    HOME          GPIO 35 |-)            (-| GPIO  3 USATO PER USB SERIALE JTAG ????? (comunque non usato)
       2SET    GUNA          GPIO 36 |-)            (x| GPIO  8  SDA (sono di default li riusiamo in OpenFire) CAM SDA
       2RST    GUNB          GPIO 37 |-)            (-| GPIO 18  SDA1 OLED SDA (OpenFIRE)
       1SET   SELECT         GPIO 38 |x)            (x| GPIO 17  SOLENOID signal (OpenFIRE)
       1RST   START          GPIO 39 |-)            (-| GPIO 16  RUMBLE signal (OpenFIRE)
            1DPAD UP         GPIO 40 |-)            (-| GPIO 15  SCL1 OLED SCL (OpenFIRE)
            1DPAD DOWN       GPIO 41 |-)            (-| GPIO  7 LIVELLO BATTERIA (non usato per il momento)
            1DPAD LEFT       GPIO 42 |-)            (-| GPIO  6 TEMP SENSOR (OpenFIRE)
            1DPAD RIGHT      GPIO  2 |x)            (x| GPIO  5 ANALOG JOYSTICK Y
    1MIDDLE     TRIGGER      GPIO  1 |-)            (-| GPIO  4 ANALOG JOYSTICK X
COLLEG 44 A USB SERIAL JTAG U0RXD RX |x)            (x| RST CHIP PU - COLLEGATO A PULSANTE RST
COLLEG 43 A USB SERIAL JTAG U0TXD TX |-)            (-| 3V3 VDD33
                                 GND |-)            (-| 3V3 VDD33
                                     |____ANTENNA_____|  

```
## Raspberry Pi Pico (Non/W)
```

                                        (_____)
                     A Button     0  |-) *USB* (-| VBUS 5v (USB voltage)
                     B Button     1  |-)       (-| VSYS (Input from Battery/Output to NeoPixels)
                                 GND |x)       (x| GND
                     C Button     2  |-)       (x| 3V3 En
                     Start        3  |-)       (p| 3V3 Out (to Display/Cam/Analog Inputs)
                     Select       4  |-)       (x| ADCVREF
                     Home Button  5  |-)       (-|  A2 Temp Sensor 28
                                 GND |x)       (x| AGND (for ADC VREF)
                     D-Pad Up     6  |-)       (-|  A1 *Unmapped* Analop X 27
                     D-Pad Down   7  |-)       (-|  A0 *Unmapped* Analog Y 26
                     D-Pad Left   8  |-)       (x| RUN
                     D-Pad Right  9  |-)       (-|  22 *Unmapped*
                                 GND |x)       (x| GND
Peripherals SDA OLED RGB Red     10  |-)       (-|  21 Camera SCL
Peripherals SCL OLED RGB Green   11  |-)       (-|  20 Camera SDA
                     RGB Blue    12  |-)       (-|  19 *Unmapped* Peripherals SCL (impostato da App)
                     Pump Action 13  |-)       (-|  18 *Unmapped* Peripherals SDA (impostato da App)
                                 GND |x)       (x| GND
                     Pedal       14  |-)       (-|  17 Rumble Signal
                     Trigger     15  |-) _|_|_ (-|  16 Solenoid Signal

```
## Adafruit ItsyBitsy RP2040
```

                                        (_____)
                               RST |x)   *USB*   (p| BAT (Battery Input)
                               3V3 |p)           (x| GND
   (Display/Cam/Rumble/Analog) 3V3 |p)           (p| USB (5V USB voltage, to NeoPixels)
            (VSYS-like Output) VHi |p)           (-|  11 C Button
               B Button        A0  |-)           (-|  10 D-Pad Right
               A Button        A1  |-)           (-|  9  D-Pad Up
               Start           A2  |-)           (-|  8  D-Pad Left
               Select          A3  |-)           (-|  7  D-Pad Down
               Rumble Signal   24  |-)           (-|  6  Trigger
               Solenoid Signal 25  |-)           (x| !5  (Output Only)
               *Unmapped*      18  |-)           (-|  3  Camera SCL
               *Unmapped*      19  |-)           (-|  2  Camera SDA
               *Unmapped*      20  |-)           (-|  0  *Unmapped*
               *Unmapped*      12  |-) x|x|x|-|- (-|  1  *Unmapped*
                                             4 5
                                     4 - *Unmapped*
                                     5 - Pedal

```
## Adafruit Keeboar KB2040
```
                                        (_____)
                   (USB Data+)  D+ |x)   *USB*   (x| D-  (USB Data-)
               *Unmapped*       0  |-)           (p| RAW (5V USB voltage, to NeoPixels)
               *Unmapped*       1  |-)           (x| GND
                               GND |x)           (x| RST
                               GND |x)           (p| 3V3 (Display/Cam/Rumble/Analog)
               Camera SDA       2  |-)           (-|  A3 A Button
               Camera SCL       3  |-)           (-|  A2 Trigger
               B Button         4  |-)           (-|  A1 Home Button
               Rumble Signal    5  |-)           (-|  A0 Temp Sensor
               Button C         6  |-)           (-|  18 D-Pad Up
               Solenoid Signal  7  |-)           (-|  20 D-Pad Down
               Select           8  |-)           (-|  19 D-Pad Left
               Start            9  |-)___________(-|  10 D-Pad Right
```
## Arduino Nano RP2040 Connect
##### *Note: A4/A5/A6/A7 are handled by the NiNa WiFi chip, but are not yet exposed to OpenFIRE.
```
                                         (_____)
                                    |     *USB*     |
                *Unmapped*       6  |-)           (-|  4  A Button
   (Display/Cam/Rumble/Analog)  3V3 |p)           (-|  7  B Button
                               AREF |x)           (-|  5  C Button
                *Unmapped*      A0  |-)           (-|  21 *Unmapped*
                *Unmapped*      A1  |-)           (-|  20 *Unmapped*
                *Unmapped*      A2  |-)           (-|  19 *Unmapped*
                *Unmapped*      A3  |-)           (-|  18 *Unmapped*
                Camera SDA      12  |x)           (-|  17 *Unmapped*
                Camera SCL      13  |x)           (-|  16 *Unmapped*
                *N/C*           A6  |x)           (-|  15 *Unmapped*
                *N/C*           A7  |x)           (-|  25 *Unmapped*
    (USB voltage, to NeoPixels)  5V |p)           (x| GND
                *Unmapped*      REC |x)           (x| RST
                *Unmapped*      GND |x)           (-|  1  Pedal
                *Unmapped*      VIN |p)           (-|  0  Trigger
                                    |_______________|
```

## Waveshare RP2040 Zero
##### *Note: The underside GPIO pads *17-25* are not yet exposed to the OpenFIRE App.
```

                                         (_____)
   (USB voltage, to NeoPixels)  5V |p)    *USB*    (-|  0  Trigger
                               GND |x)             (-|  1  A Button
   (Display/Cam/Rumble/Analog) 3V3 |p)             (-|  2  B Button
               Temp Sensor     A3  |-)             (-|  3  C Button
               *Unmapped*      A2  |-)             (-|  4  Start
               Camera SCL      A1  |-)             (-|  5  Select
               Camera SDA      A0  |-)             (-|  6  *Unmapped*
               *Unmapped*      15  |-)             (-|  7  *Unmapped*
               *Unmapped*      14  |-)-| -| - |- |-(-|  8  *Unmapped*
                                     13 12 11  10 9
                                      9 - *Unmapped*
                                     10 - *Unmapped*
                                     11 - *Unmapped*
                                     12 - *Unmapped*
                                     13 - *Unmapped*

```
## VCC-GND YD RP2040
```

                          *no default layout yet, to be added*

```
