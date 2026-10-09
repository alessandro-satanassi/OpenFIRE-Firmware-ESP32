/*!
 * @file OpenFIRECamera.cpp
 * @brief Runtime camera manager and common camera API for OpenFIRE.
 *
 * @copyright Alessandro Satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author Alessandro Satanassi
 * @version V1.0
 * @date 2026
 */

#include <Arduino.h>
#include <Wire.h>
#include <SPI.h>
#include <DFRobotIRPositionEx.h>
#include <PixArt_PAJ7025.h>
#include "OpenFIRECamera.h"

#include "OpenFIREprefs.h"

#ifdef ARDUINO_ARCH_RP2040
#include <pico/stdlib.h>
#include <hardware/clocks.h>
#include <hardware/pwm.h>
#endif

namespace {
TwoWire* activeI2C = nullptr;
SPIClass* activeSPI = nullptr;
int8_t activeWiiClockPin = -1;
DFRobotIRPositionEx* dfrCamera = nullptr;
PAJ7025* pajCamera = nullptr;
int pajX[4] = {0};
int pajY[4] = {0};

#ifdef ARDUINO_ARCH_ESP32
SPIClass pajSPI(FSPI);
#endif
}

const CameraProfile* OpenFIRECamera::activeProfile = nullptr;
const OpenFIRECamera::CameraOps* OpenFIRECamera::activeOps = nullptr;
OpenFIRECamera::ReadFn OpenFIRECamera::activeRead = nullptr;
const int* OpenFIRECamera::activeX = nullptr;
const int* OpenFIRECamera::activeY = nullptr;
unsigned int OpenFIRECamera::activeSeen = 0;
OpenFIRECamera::ObjectData OpenFIRECamera::objectData[4] = {};
uint16_t OpenFIRECamera::activeExtendedCapabilities = 0;
OpenFIRECamera::DataFormat_e OpenFIRECamera::activeFormat = OpenFIRECamera::DataFormat_Basic;
bool OpenFIRECamera::ready = false;

void OpenFIRECamera::ClearObjectData() {
    for (int i = 0; i < 4; i++) {
        objectData[i] = {};
    }
}

bool OpenFIRECamera::Select() {
    // The DFRobot Extended format is read with the camera's Full format (ReadDFRobotFull).
    static constexpr uint16_t DFRobotExtendedCapabilities =
        ExtendedData_Size |
        ExtendedData_MaxBrightness |
        ExtendedData_Boundaries;

    static constexpr uint16_t PAJ7025ExtendedCapabilities =
        ExtendedData_Size |
        ExtendedData_Area |
        ExtendedData_AverageBrightness |
        ExtendedData_MaxBrightness |
        ExtendedData_Range |
        ExtendedData_Radius |
        ExtendedData_Boundaries |
        ExtendedData_AspectRatio |
        ExtendedData_Velocity;

    static const CameraOps DFRobotOps = {
        &OpenFIRECamera::BeginDFRobot,
        &OpenFIRECamera::ReadDFRobotBasic,
        &OpenFIRECamera::ReadDFRobotFull,
        &OpenFIRECamera::DataFormatDFRobot,
        &OpenFIRECamera::SensitivityDFRobot,
        &OpenFIRECamera::EndDFRobot,
        DFRobotExtendedCapabilities
    };

    // R2 and R3 share the same backend; only the sensitivity presets differ.
    static const CameraOps PAJ7025R2Ops = {
        &OpenFIRECamera::BeginPAJ7025,
        &OpenFIRECamera::ReadPAJ7025Basic,
        &OpenFIRECamera::ReadPAJ7025Extended,
        &OpenFIRECamera::DataFormatPAJ7025,
        &OpenFIRECamera::SensitivityPAJ7025R2,
        &OpenFIRECamera::EndPAJ7025,
        PAJ7025ExtendedCapabilities
    };

    static const CameraOps PAJ7025R3Ops = {
        &OpenFIRECamera::BeginPAJ7025,
        &OpenFIRECamera::ReadPAJ7025Basic,
        &OpenFIRECamera::ReadPAJ7025Extended,
        &OpenFIRECamera::DataFormatPAJ7025,
        &OpenFIRECamera::SensitivityPAJ7025R3,
        &OpenFIRECamera::EndPAJ7025,
        PAJ7025ExtendedCapabilities
    };

    if (ready) {
        return false;
    }

    activeSeen = 0;
    activeX = nullptr;
    activeY = nullptr;
    activeFormat = DataFormat_Basic;
    ClearObjectData();

    const CameraModel model = (CameraModel)OF_Prefs::settings[OF_Const::cameraModel];

    switch (model) {
        // An unknown saved model (damaged settings, or saved by a newer firmware) uses the
        // DFRobot, the model of a clean installation, so that the profile is never empty.
        default:
        case OF_Const::DFRobot_SEN0158:
            activeProfile = &OpenFIRE_CameraProfiles::DFRobot_SEN0158;
            activeOps = &DFRobotOps;
            break;
        case OF_Const::PixArt_PAJ7025R2:
            activeProfile = &OpenFIRE_CameraProfiles::PixArt_PAJ7025R2;
            activeOps = &PAJ7025R2Ops;
            break;
        case OF_Const::PixArt_PAJ7025R3:
            activeProfile = &OpenFIRE_CameraProfiles::PixArt_PAJ7025R3;
            activeOps = &PAJ7025R3Ops;
            break;
    }

    activeRead = activeOps->readBasic;
    activeExtendedCapabilities = activeOps->extendedCapabilities;
    return true;
}

const CameraProfile& OpenFIRECamera::Profile() {
    return *activeProfile;
}

CameraModel OpenFIRECamera::Model() {
    return activeProfile->model;
}

bool OpenFIRECamera::Begin()
{
    if (ready)
        return false;

    if (!Select())
        return false;

    const uint8_t sensitivity = ClampSensitivity(
        (uint8_t)OF_Prefs::profiles[OF_Prefs::currentProfile].irSens
    );

    ready = activeOps->begin(sensitivity);
    return ready;
}

void OpenFIRECamera::End() {
    if (activeOps != nullptr) {
        activeOps->end();
    }

    ready = false;
    activeX = nullptr;
    activeY = nullptr;
    activeSeen = 0;
    ClearObjectData();
}

bool OpenFIRECamera::SetDataFormat(DataFormat_e format) {
    if (activeOps == nullptr) {
        return false;
    }

    if (format != DataFormat_Basic && format != DataFormat_Extended) {
        return false;
    }

    activeFormat = format;
    activeRead = (format == DataFormat_Extended) ? activeOps->readExtended : activeOps->readBasic;
    ClearObjectData();

    if (ready) {
        activeOps->dataFormat(format);
    }

    return true;
}

void OpenFIRECamera::SetSensitivity(Sensitivity_e sensitivity) {
    if (activeOps != nullptr && ready) {
        activeOps->sensitivity(ClampSensitivity((uint8_t)sensitivity));
    }
}

bool OpenFIRECamera::BeginDFRobot(uint8_t sensitivity) {
    const int8_t pin_sda = OF_Prefs::pins[OF_Const::camSDA];
    const int8_t pin_scl = OF_Prefs::pins[OF_Const::camSCL];
    const int8_t pin_wiiClock = OF_Prefs::pins[OF_Const::wiiClockGen];

    if (pin_sda < 0 || pin_scl < 0) {
        return false;
    }
    
    #ifdef ARDUINO_ARCH_ESP32
    Wire.setPins(pin_sda, pin_scl);
    dfrCamera = new DFRobotIRPositionEx(Wire);

    #ifdef CLOCK_CAM_WII
    constexpr uint32_t WII_CLOCK_FREQUENCY_HZ = 20000000U;
    constexpr uint32_t WII_CLOCK_DUTY_CYCLE = 1U;

    if (pin_wiiClock > -1) {
        if (!ledcSetClockSource(LEDC_AUTO_CLK)) {
            log_e("ERRORE: ledcSetClockSource fallita!");
        }
        if (!ledcAttach(pin_wiiClock, WII_CLOCK_FREQUENCY_HZ, 1)) {
            log_e("ERRORE: ledcAttach fallita!");
        } else if (!ledcWrite(pin_wiiClock, WII_CLOCK_DUTY_CYCLE)) {
            log_e("ERRORE: ledcWrite fallita!");
            ledcDetach(pin_wiiClock);
        } else {
            activeWiiClockPin = pin_wiiClock;
        }
    } else {
        log_e("Clock non attivato (GPIO %d non valido: <= -1).\n", pin_wiiClock);
    }
    #endif
#else  // rp2040
    if (bitRead(pin_scl, 1) && bitRead(pin_sda, 1)) {
        if (bitRead(pin_scl, 0) && !bitRead(pin_sda, 0)) {
            Wire1.setSDA(pin_sda);
            Wire1.setSCL(pin_scl);
            activeI2C = &Wire1;
            dfrCamera = new DFRobotIRPositionEx(Wire1);
        }
    } else if (!bitRead(pin_scl, 1) && !bitRead(pin_sda, 1)) {
        if (bitRead(pin_scl, 0) && !bitRead(pin_sda, 0)) {
            Wire.setSDA(pin_sda);
            Wire.setSCL(pin_scl);
            activeI2C = &Wire;
            dfrCamera = new DFRobotIRPositionEx(Wire);
        }
    }

    #ifdef CLOCK_CAM_WII
    if (pin_wiiClock > -1) {
        set_sys_clock_khz(125000, true);
        gpio_set_function(pin_wiiClock, GPIO_FUNC_PWM);
        const uint slice_num = pwm_gpio_to_slice_num(pin_wiiClock);
        pwm_set_clkdiv(slice_num, 1.0f);
        pwm_set_wrap(slice_num, 4);
        pwm_set_chan_level(slice_num, pwm_gpio_to_channel(pin_wiiClock), 2);
        pwm_set_enabled(slice_num, true);
        activeWiiClockPin = pin_wiiClock;
    } 
   #endif
#endif
    
    if (dfrCamera == nullptr) {
        EndDFRobot();
        return false;
    }

    const DFRobotIRPositionEx::DataFormat_e format =
        (activeFormat == DataFormat_Extended)
            ? DFRobotIRPositionEx::DataFormat_Full
            : DFRobotIRPositionEx::DataFormat_Basic;

    if (!dfrCamera->begin(activeProfile->busClock,
                          format,
                          (DFRobotIRPositionEx::Sensitivity_e)sensitivity)) {
        EndDFRobot();
        return false;
    }

    activeX = dfrCamera->xPositions();
    activeY = dfrCamera->yPositions();
    activeSeen = dfrCamera->seen();
    return true;
}

int OpenFIRECamera::ReadDFRobotBasic() {
    const int error = dfrCamera->basicAtomic(DFRobotIRPositionEx::Retry_2);
    if (error >= DFRobotIRPositionEx::Error_Success) {
        activeSeen = dfrCamera->seen();
    }
    return error;
}

// Reads the camera's Extended format. Not bound at the moment: the OpenFIRE Extended
// format is read with the camera's Full format (ReadDFRobotFull), which also gives the
// intensity. Kept to go back to it by binding it in DFRobotOps (and DataFormat_Extended
// in BeginDFRobot and DataFormatDFRobot).
int OpenFIRECamera::ReadDFRobotExtended() {
    // The DFRobot Extended format measures only x, y and size (0..15).
    // The remaining ObjectData fields are filled here, so that consumers
    // (e.g. the IR test view) handle every camera in the same way.
    //
    // Area: size is proportional to the blob diameter, so the equivalent area
    // grows with size^2. A blob of DFR_REF_SIZE gets DFR_REF_AREA, which the
    // web app IR test draws with the radius of the former fixed circle
    // (25 px in 1920x1080 test space, see
    // IRTEST_BLOB_AREA_MAX and IRTEST_BLOB_RADIUS_MAX in fullscreen.js):
    // 300 * (25 / 60)^2 = 52. Other sizes scale in proportion.
    // IRTEST_BLOB_RADIUS_SCALE in fullscreen.js then enlarges every circle of
    // the IR test view alike, so this reference does not depend on it.
    // Brightness is not measured: fixed neutral values are used.
    // Range, radius, boundaries, aspect ratio and velocity are not measured
    // and not used by any consumer: they stay at zero.
    static constexpr uint32_t DFR_REF_SIZE = 10U;
    static constexpr uint32_t DFR_REF_AREA = 52U;
    static constexpr uint8_t DFR_AVERAGE_BRIGHTNESS = 190U;
    static constexpr uint8_t DFR_MAX_BRIGHTNESS = 255U;

    const int error = dfrCamera->extendedAtomic(DFRobotIRPositionEx::Retry_2);

    if (error >= DFRobotIRPositionEx::Error_Success) {
        activeSeen = dfrCamera->seen();

        for (int i = 0; i < 4; i++) {
            if ((activeSeen & (1U << i)) != 0U) {
                const int size = dfrCamera->size(i);
                const uint32_t size2 = (uint32_t)(size & 0x0F) * (uint32_t)(size & 0x0F);
                // Rounded; a seen blob never reports area 0 (0 = no object on PAJ7025).
                uint32_t area = (DFR_REF_AREA * size2 + (DFR_REF_SIZE * DFR_REF_SIZE) / 2U) /
                                (DFR_REF_SIZE * DFR_REF_SIZE);
                if (area == 0U) area = 1U;

                objectData[i].valid = true;
                objectData[i].x = dfrCamera->x(i);
                objectData[i].y = dfrCamera->y(i);
                objectData[i].size = size;
                objectData[i].area = (uint16_t)area;
                objectData[i].averageBrightness = DFR_AVERAGE_BRIGHTNESS;
                objectData[i].maxBrightness = DFR_MAX_BRIGHTNESS;
            } else {
                objectData[i].valid = false;
            }
        }
    }

    return error;
}

int OpenFIRECamera::ReadDFRobotFull() {
    // The DFRobot Full format measures x, y, size (0..15), the bounding box of the blob
    // in the sensor's 128x96 array and an 8 bit intensity. The fields are converted here
    // to the common ObjectData scale, so that consumers handle every camera alike.
    //
    // Area: pixels of the bounding box, the same unit as the PAJ7025 area (pixels of
    // the sensor); a box is about 27% larger than a round blob. To check on the hardware.
    // Boundaries: the bounding box as the camera gives it.
    //
    // Brightness, PHASE 1 (measurements): the DFRobot intensity is on a much lower scale
    // than the PAJ7025 brightness (LEDs about 4..30, 4 at the detection limit, measured by
    // the LIGHTGUN-STUDIO project). averageBrightness carries the raw intensity, to be read
    // with the WebApp IR test (?irdebug in its address); maxBrightness a provisional linear
    // conversion: DFR_INTENSITY_MIN -> PAJ_BRIGHTNESS_MIN (the PAJ7025 detection limit),
    // DFR_INTENSITY_FULL -> 255. PHASE 2: both get the conversion measured on the hardware.
    static constexpr int32_t DFR_INTENSITY_MIN = 4;
    static constexpr int32_t DFR_INTENSITY_FULL = 30;
    static constexpr int32_t PAJ_BRIGHTNESS_MIN = 130;
    static constexpr int32_t BRIGHTNESS_MAX = 255;

    const int error = dfrCamera->fullAtomic(DFRobotIRPositionEx::Retry_2);

    if (error >= DFRobotIRPositionEx::Error_Success) {
        activeSeen = dfrCamera->seen();

        for (int i = 0; i < 4; i++) {
            if ((activeSeen & (1U << i)) != 0U) {
                const DFRobotIRPositionEx::Box_t& box = dfrCamera->box(i);
                const uint32_t width = (box.xMax > box.xMin) ? (uint32_t)(box.xMax - box.xMin) : 0U;
                const uint32_t height = (box.yMax > box.yMin) ? (uint32_t)(box.yMax - box.yMin) : 0U;

                const int32_t intensity = dfrCamera->intensity(i);
                int32_t brightness = PAJ_BRIGHTNESS_MIN +
                    (intensity - DFR_INTENSITY_MIN) * (BRIGHTNESS_MAX - PAJ_BRIGHTNESS_MIN) /
                    (DFR_INTENSITY_FULL - DFR_INTENSITY_MIN);
                if (brightness < PAJ_BRIGHTNESS_MIN) brightness = PAJ_BRIGHTNESS_MIN;
                if (brightness > BRIGHTNESS_MAX) brightness = BRIGHTNESS_MAX;

                objectData[i].valid = true;
                objectData[i].x = dfrCamera->x(i);
                objectData[i].y = dfrCamera->y(i);
                objectData[i].size = dfrCamera->size(i);
                objectData[i].area = (uint16_t)((width + 1U) * (height + 1U)); // 7 bit box: fits 16 bits (the sensor gives at most 128 * 96)
                objectData[i].averageBrightness = (uint8_t)intensity;            // PHASE 1: raw value
                objectData[i].maxBrightness = (uint8_t)brightness;
                objectData[i].boundaryLeft = box.xMin;
                objectData[i].boundaryRight = box.xMax;
                objectData[i].boundaryUp = box.yMin;
                objectData[i].boundaryDown = box.yMax;
            } else {
                objectData[i].valid = false;
            }
        }
    }

    return error;
}

void OpenFIRECamera::DataFormatDFRobot(DataFormat_e format) {
    dfrCamera->dataFormat(
        (format == DataFormat_Extended)
            ? DFRobotIRPositionEx::DataFormat_Full
            : DFRobotIRPositionEx::DataFormat_Basic
    );
}

void OpenFIRECamera::SensitivityDFRobot(uint8_t sensitivity) {
    dfrCamera->sensitivityLevel((DFRobotIRPositionEx::Sensitivity_e)sensitivity);
}

void OpenFIRECamera::EndDFRobot() {
    if (dfrCamera != nullptr) {
        delete dfrCamera;
        dfrCamera = nullptr;
    }

#ifdef ARDUINO_ARCH_ESP32
    Wire.end();
    #ifdef CLOCK_CAM_WII
    if (activeWiiClockPin > -1) {
        ledcDetach(activeWiiClockPin);
        activeWiiClockPin = -1;
    }
    #endif
#elif defined(ARDUINO_ARCH_RP2040)
    if (activeI2C) {
        activeI2C->end();
        activeI2C = nullptr;
    }
    #ifdef CLOCK_CAM_WII
    if (activeWiiClockPin > -1) {
        // set_sys_clock_khz(133000, true);
        pwm_set_enabled(pwm_gpio_to_slice_num(activeWiiClockPin), false);
        gpio_set_function(activeWiiClockPin, GPIO_FUNC_SIO);
        gpio_put(activeWiiClockPin, 0);
        gpio_set_dir(activeWiiClockPin, GPIO_IN);
        activeWiiClockPin = -1;
        
        set_sys_clock_khz((uint32_t)(F_CPU / 1000UL), true);  // ripristina il clock originale funzio sia su rp2040 cherp2350/pico 2
    }
    #endif
#endif
}

// ============================================================================
// [3] SEZIONE: PixArt PAJ7025 (Driver Nativo)
// ============================================================================
bool OpenFIRECamera::BeginPAJ7025(uint8_t sensitivity) {
    const int8_t pin_spiSck = OF_Prefs::pins[OF_Const::cam_SPI_SCK];
    const int8_t pin_spiMiso = OF_Prefs::pins[OF_Const::cam_SPI_MISO];
    const int8_t pin_spiMosi = OF_Prefs::pins[OF_Const::cam_SPI_MOSI];
    const int8_t pin_spiCs = OF_Prefs::pins[OF_Const::cam_SPI_CS];
    
    if (pin_spiSck < 0 || pin_spiMiso < 0 || pin_spiMosi < 0 || pin_spiCs < 0) {
        return false;
    }

#ifdef ARDUINO_ARCH_ESP32
    if (!pajSPI.begin(pin_spiSck, pin_spiMiso, pin_spiMosi, pin_spiCs)) {
        pajSPI.end();
        return false;
    }
    if (pajCamera == nullptr) pajCamera = new PAJ7025();

    if (pajCamera == nullptr || !pajCamera->begin(&pajSPI, pin_spiCs, activeProfile->busClock)) {
        EndPAJ7025();
        return false;
    }
#else // rp2040
    // RP2040/RP2350:
    // bits 0-1 identify the SPI function:
    // 0 = MISO/RX, 1 = CS, 2 = SCK, 3 = MOSI/TX
    // bit 3 identifies the controller: 0 = SPI0, 1 = SPI1.
    const uint8_t spiBus = bitRead(pin_spiSck, 3);

    if ((uint8_t)pin_spiSck  >= (uint8_t)NUM_BANK0_GPIOS ||
            (uint8_t)pin_spiMiso >= (uint8_t)NUM_BANK0_GPIOS ||
            (uint8_t)pin_spiMosi >= (uint8_t)NUM_BANK0_GPIOS ||
            (uint8_t)pin_spiCs   >= (uint8_t)NUM_BANK0_GPIOS ||
            (pin_spiMiso & 0x03) != 0x00 ||
            (pin_spiCs   & 0x03) != 0x01 ||
            (pin_spiSck  & 0x03) != 0x02 ||
            (pin_spiMosi & 0x03) != 0x03 ||
            bitRead(pin_spiMiso, 3) != spiBus ||
            bitRead(pin_spiMosi, 3) != spiBus ||
            bitRead(pin_spiCs,   3) != spiBus) {
        return false;
    }

    SPIClassRP2040* spiPort = spiBus ? &SPI1 : &SPI;

    if (!spiPort->setSCK(pin_spiSck) ||
            !spiPort->setMISO(pin_spiMiso) ||
            !spiPort->setMOSI(pin_spiMosi) ||
            !spiPort->setCS(pin_spiCs)) {
        return false;
    }

    activeSPI = spiPort;
    activeSPI->begin();

    if (pajCamera == nullptr) pajCamera = new PAJ7025();

    if (pajCamera == nullptr || !pajCamera->begin(activeSPI, pin_spiCs, activeProfile->busClock)) {
        EndPAJ7025();
        return false;
    }

#endif



    pajCamera->setFrameRate(activeProfile->fps);
    pajCamera->setExposure(300);
    // Use the selected model's preset also at startup, before ready is set.
    activeOps->sensitivity(sensitivity);
    pajCamera->setResolution((uint16_t)activeProfile->camMaxX, (uint16_t)activeProfile->camMaxY);

    activeX = pajX;
    activeY = pajY;
    activeSeen = 0;
    return true;
}

int OpenFIRECamera::ReadPAJ7025Basic() {
    PAJ7025_Object rawObjects[4];
    pajCamera->readData(rawObjects, PAJ7025_FORMAT_BASIC);

    activeSeen = 0;
    for (int i = 0; i < 4; i++) {
        if (rawObjects[i].is_valid) {
            pajX[i] = rawObjects[i].cx;
            pajY[i] = rawObjects[i].cy;
            activeSeen |= (1U << i);
        }
    }

    return Error_Success;
}

int OpenFIRECamera::ReadPAJ7025Extended() {
    PAJ7025_Object rawObjects[4];
    pajCamera->readData(rawObjects, PAJ7025_FORMAT_EXTENDED);

    activeSeen = 0;
    for (int i = 0; i < 4; i++) {
        if (rawObjects[i].is_valid) {
            // Array di X e Y di base per i ritorni costanti (activeX / activeY)
            pajX[i] = rawObjects[i].cx;
            pajY[i] = rawObjects[i].cy;

            // Direct struct copy - evita del tutto i wrapper intermedi e non fa double copy!
            objectData[i].valid = true;
            objectData[i].x = rawObjects[i].cx;
            objectData[i].y = rawObjects[i].cy;
            objectData[i].size = (rawObjects[i].area > 15U) ? 15 : (int)rawObjects[i].area;
            objectData[i].area = rawObjects[i].area;
            objectData[i].averageBrightness = rawObjects[i].average_brightness;
            objectData[i].maxBrightness = rawObjects[i].max_brightness;
            objectData[i].range = rawObjects[i].range;
            objectData[i].radius = rawObjects[i].radius;
            objectData[i].boundaryLeft = rawObjects[i].boundary_left;
            objectData[i].boundaryRight = rawObjects[i].boundary_right;
            objectData[i].boundaryUp = rawObjects[i].boundary_up;
            objectData[i].boundaryDown = rawObjects[i].boundary_down;
            objectData[i].aspectRatio = rawObjects[i].aspect_ratio;
            objectData[i].vx = rawObjects[i].vx;
            objectData[i].vy = rawObjects[i].vy;
            
            activeSeen |= (1U << i);
        } else {
            // Ottimizzazione Hot-Path 1: invalidiamo il punto senza usare uno sprecone memset / {}
            // Le coordinate X e Y sporche verranno semplicemente ignorate dal sistema di mira grazie a valid=false e a seenFlags.
            objectData[i].valid = false;
        }
    }

    return Error_Success;
}

void OpenFIRECamera::DataFormatPAJ7025(DataFormat_e format) {
    (void)format;
    // PAJ7025 format is selected directly by the bound read function.
}

void OpenFIRECamera::SensitivityPAJ7025R2(uint8_t sensitivity) {
    if (pajCamera == nullptr) return;
    
    if (sensitivity == 0U) {
        pajCamera->setGain(0x10, 0x00);
        pajCamera->setDSP(2, 130, 150, 40);
    }
    else if (sensitivity == 1U) {
        pajCamera->setGain(0x10, 0x02);
        pajCamera->setDSP(2, 130, 200, 50);
    }
    else {
        pajCamera->setGain(0x10, 0x03);
        pajCamera->setDSP(1, 150, 300, 60);
    }
}

void OpenFIRECamera::SensitivityPAJ7025R3(uint8_t sensitivity) {
    if (pajCamera == nullptr) return;

    // Initial R3 presets to validate on hardware, keeping 300 us exposure.
    // setDSP arguments: minimum area, brightness threshold, maximum area, noise threshold.
    if (sensitivity == 0U) {
        pajCamera->setGain(0x10, 0x02); // 4x
        pajCamera->setDSP(1, 130, 150, 40);
    }
    else if (sensitivity == 1U) {
        pajCamera->setGain(0x08, 0x03); // 6x
        pajCamera->setDSP(1, 130, 200, 50);
    }
    else {
        pajCamera->setGain(0x10, 0x03); // 8x
        pajCamera->setDSP(1, 130, 300, 60);
    }
}

void OpenFIRECamera::EndPAJ7025() {
    if (pajCamera != nullptr) {
        delete pajCamera;
        pajCamera = nullptr;
    }
#ifdef ARDUINO_ARCH_ESP32
    pajSPI.end();
#else
    if (activeSPI) {
        activeSPI->end();
        activeSPI = nullptr;
    }
#endif
}
