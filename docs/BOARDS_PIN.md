# Supported Boards Layouts:
 - [Waveshare esp32-s3-pico](#waveshare-esp32-s3-pico)
 - [ESP32_S3_WROOM1_DevKitC_1_N16R8](#esp32_s3_wroom1_devkitc_1_n16r8)
 - [Waveshare esp32-s3-zero](#waveshare-esp32-s3-zero)


```
                                        Symbol Legend:
                       (x) = GND/No Connect | (-) = GPIO | (p) = Power
```

### 

## waveshare-esp32-s3-pico
```

;                                     ______(*USB*)______     
;                           (GPIO 11) |-)  1  PIN  40 (p| (VBUS)
;                           (GPIO 12) |-)  2       39 (p| (VSYS)
;                               (GND) |x)  3       38 (x| (GND)
;                           (GPIO 13) |-)  4       37 (x| (EN)
;                           (GPIO 14) |-)  5       36 (p| (3V3 OUT)                    
;                           (GPIO 15) |-)  6       35 (-| (GPIO 10)   
;                           (GPIO 16) |-)  7       34 (-| (GPIO 09)
;                               (GND) |x)  8       33 (x| (GND) 
;                           (GPIO 17) |-)  9       32 (-| (GPIO 08)
;                           (GPIO 18) |-) 10       31 (-| (GPIO 07)
;                           (GPIO 33) |-) 11       30 (x| (RUN)
;                           (GPIO 34) |-) 12       29 (-| (GPIO 06)
;                               (GND) |x) 13       28 (x| (GND)
;                           (GPIO 35) |-) 14       27 (-| (GPIO 05)
;                           (GPIO 36) |-) 15       26 (-| (GPIO 04)
;                           (GPIO 37) |-) 16       25 (-| (GPIO 02)
;                           (GPIO 38) |-) 17       24 (-| (GPIO 01)
;                               (GND) |x) 18       23 (x| (GND)
;                           (GPIO 39) |-) 19       22 (-| (GPIO 41)
;                           (GPIO 40) |-) 20  PIN  21 (-| (GPIO 42)
;                                     |_________________|

```

## esp32_s3_wroom1_devkitc_1_n16r8
```

;                                          SER     OTG
;                                      ___(USB)___(USB)___  
;                                (GND) |x)  1  PIN  44 (x| (GND)
;                                (GND) |x)  2       43 (p| (5V OUT-IN)
;                            (GPIO 19) |-)  3       42 (-| (GPIO 14)
;                            (GPIO 20) |-)  4       41 (-| (GPIO 13) 
;                            (GPIO 21) |-)  5       40 (-| (GPIO 12) 
;                            (GPIO 47) |-)  6       39 (-| (GPIO 11) 
;                            (GPIO 48) |-)  7       38 (-| (GPIO 10) 
;                            (GPIO 45) |-)  8       37 (-| (GPIO 09) 
;                            (GPIO 00) |-)  9       36 (-| (GPIO 46)
;                            (GPIO 35) |-) 10       35 (-| (GPIO 03)
;                            (GPIO 36) |-) 11       34 (-| (GPIO 08) 
;                            (GPIO 37) |-) 12       33 (-| (GPIO 18) 
;                            (GPIO 38) |-) 13       32 (-| (GPIO 17) 
;                            (GPIO 39) |-) 14       31 (-| (GPIO 16) 
;                            (GPIO 40) |-) 15       30 (-| (GPIO 15) 
;                            (GPIO 41) |-) 16       29 (-| (GPIO 07)
;                            (GPIO 42) |-) 17       28 (-| (GPIO 06)
;                            (GPIO 02) |-) 18       27 (-| (GPIO 05)
;                            (GPIO 01) |-) 19       26 (-| (GPIO 04)
;                                 (RX) |-) 20       25 (x| (RST)
;                                 (TX) |-) 21       24 (p| (3V3 OUT)
;                                (GND) |x) 22  PIN  23 (p| (3V3 OUT)
;                                      |_________________|  
;                                          | ANTENNA |
;                                          |_________|

```

## waveshare-esp32-s3-zero
```

;                                     ______(*USB*)______     
;                                (5V) |p)  1  PIN  18 (-| (TX)
;                               (GND) |x)  2       17 (-| (RX)
;                           (3V3 OUT) |p)  3       16 (-| (GPIO 13)
;                           (GPIO 01) |-)  4       15 (-| (GPIO 12)
;                           (GPIO 02) |-)  5       14 (-| (GPIO 11)    
;                           (GPIO 03) |-)  6       13 (-| (GPIO 10)   
;                           (GPIO 04) |-)  7       12 (-| (GPIO 09)                  
;                           (GPIO 05) |-)  8       11 (-| (GPIO 08)                 
;                           (GPIO 06) |-)  9  PIN  10 (-| (GPIO 07)         
;                                     |_________________|


```
