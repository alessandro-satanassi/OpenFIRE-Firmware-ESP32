# Default board layouts

These are the pins each supported board uses after a clean installation. Every pin can be changed in the **Board Layout** tab of the WebApp, which also shows the board picture with each pin highlighted; the WebApp is always the reference for your firmware version. The illustrated pinouts are in the [lightgun guide](../README.md#english-version).

> **Italiano:** questa pagina elenca i pin predefiniti di ogni scheda supportata dopo un'installazione pulita. Ogni pin si può cambiare nella sezione **Layout Scheda** della WebApp, che mostra anche l'immagine della scheda con ogni pin evidenziato. Le tabelle usano i nomi delle funzioni in inglese, come nella WebApp in inglese; i pinout illustrati sono nella [guida lightgun](../README.md#versione-italiana).

## ESP32-S3 boards

The PixArt PAJ7025R2/R3 cameras use four SPI pins (**Camera SPI MISO, MOSI, SCK, CS**): the DevKitC-1 has them by default on GPIO 13, 11, 12 and 10; on the PICO and the ZERO assign them to free pins. The DFRobot/Wii camera uses **Camera SDA** and **Camera SCL**. An OLED display uses **Peripherals SDA** and **Peripherals SCL** (default on the DevKitC-1 and the PICO; choose free pins on the ZERO). On the DevKitC-1, GPIO 0, 19, 20, 43, 44 and 48 are reserved (BOOT button, USB, serial port, built-in RGB LED); pins that are not listed on any board are not available.

### ESP32-S3-DevKitC-1 (N16R8 / N8R2)

| GPIO | Default function |
| ---: | --- |
| 1 | Trigger |
| 2 | D-Pad Right |
| 4 | Analog Stick X |
| 5 | Analog Stick Y |
| 6 | Temperature Sensor |
| 8 | Camera SDA |
| 9 | Camera SCL |
| 10 | Camera SPI CS |
| 11 | Camera SPI MOSI |
| 12 | Camera SPI SCK |
| 13 | Camera SPI MISO |
| 14 | External NeoPixel |
| 15 | Peripherals SCL |
| 16 | Rumble Signal |
| 17 | Solenoid Signal |
| 18 | Peripherals SDA |
| 21 | Button C |
| 35 | Home Button |
| 36 | Button A |
| 37 | Button B |
| 38 | Select |
| 39 | Start |
| 40 | D-Pad Up |
| 41 | D-Pad Down |
| 42 | D-Pad Left |
| 45 | Pump Action |

Free GPIO, available in Board Layout: 3, 7, 46, 47.

### Waveshare ESP32-S3-PICO

| GPIO | Default function |
| ---: | --- |
| 4 | Camera SDA |
| 5 | Camera SCL |
| 7 | Analog Stick Y |
| 8 | Analog Stick X |
| 9 | Temperature Sensor |
| 11 | Button A |
| 12 | Button B |
| 13 | Button C |
| 14 | Start |
| 15 | Select |
| 16 | Home Button |
| 17 | D-Pad Up |
| 18 | D-Pad Down |
| 33 | D-Pad Left |
| 34 | D-Pad Right |
| 35 | Peripherals SDA |
| 36 | Peripherals SCL |
| 37 | External NeoPixel |
| 38 | Pump Action |
| 40 | Trigger |
| 41 | Rumble Signal |
| 42 | Solenoid Signal |

Free GPIO, available in Board Layout: 1, 2, 6, 10, 39.

### Waveshare ESP32-S3-ZERO (N8R8 / N4R2)

| GPIO | Default function |
| ---: | --- |
| 4 | Camera SDA |
| 5 | Camera SCL |
| 6 | Trigger |
| 7 | Button A |
| 8 | Button B |
| 9 | Button C |
| 10 | Start |
| 11 | Select |
| 12 | Rumble Signal |
| 13 | Solenoid Signal |

Free GPIO, available in Board Layout: 1, 2, 3, 14, 15, 16, 17, 18, 38, 39, 40, 41, 42, 45.


---

## RP2040 boards (wired only)

These layouts come from the original OpenFIRE project. On RP2040 boards the firmware works by USB cable only, without the wireless features.

```
                                        Symbol Legend:
                       (x) = GND/No Connect | (-) = GPIO | (p) = Power
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
                     D-Pad Up     6  |-)       (-|  A1 *Unmapped* Analog X 27
                     D-Pad Down   7  |-)       (-|  A0 *Unmapped* Analog Y 26
                     D-Pad Left   8  |-)       (x| RUN
                     D-Pad Right  9  |-)       (-|  22 *Unmapped*
                                 GND |x)       (x| GND
Peripherals SDA OLED RGB Red     10  |-)       (-|  21 Camera SCL
Peripherals SCL OLED RGB Green   11  |-)       (-|  20 Camera SDA
                     RGB Blue    12  |-)       (-|  19 *Unmapped* Peripherals SCL (assign in Board Layout)
                     Pump Action 13  |-)       (-|  18 *Unmapped* Peripherals SDA (assign in Board Layout)
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
