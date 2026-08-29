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
#include <OpenFIRE_PAJ7025.h>
#include "OpenFIRECamera.h"

#include "OpenFIREprefs.h"

#ifdef ARDUINO_ARCH_RP2040
#include <pico/stdlib.h>
#include <hardware/clocks.h>
#include <hardware/pwm.h>
#endif

namespace {
DFRobotIRPositionEx* dfrCamera = nullptr;
OpenFIRE_PAJ7025* pajCamera = nullptr;

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

bool OpenFIRECamera::Select(CameraModel model) {
    static constexpr uint16_t DFRobotExtendedCapabilities =
        ExtendedData_Size;

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
        &OpenFIRECamera::ReadDFRobotExtended,
        &OpenFIRECamera::DataFormatDFRobot,
        &OpenFIRECamera::SensitivityDFRobot,
        &OpenFIRECamera::EndDFRobot,
        DFRobotExtendedCapabilities
    };

    static const CameraOps PAJ7025Ops = {
        &OpenFIRECamera::BeginPAJ7025,
        &OpenFIRECamera::ReadPAJ7025Basic,
        &OpenFIRECamera::ReadPAJ7025Extended,
        &OpenFIRECamera::DataFormatPAJ7025,
        &OpenFIRECamera::SensitivityPAJ7025,
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

    switch (model) {
        case OF_Const::DFRobot_SEN0158:
            activeProfile = &OpenFIRE_CameraProfiles::DFRobot_SEN0158;
            activeOps = &DFRobotOps;
            break;
        case OF_Const::PixArt_PAJ7025R2:
            activeProfile = &OpenFIRE_CameraProfiles::PixArt_PAJ7025R2;
            activeOps = &PAJ7025Ops;
            break;
        case OF_Const::PixArt_PAJ7025R3:
            activeProfile = &OpenFIRE_CameraProfiles::PixArt_PAJ7025R3;
            activeOps = &PAJ7025Ops;
            break;
        default:
            activeProfile = nullptr;
            activeOps = nullptr;
            activeRead = nullptr;
            activeExtendedCapabilities = 0;
            return false;
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

    const CameraModel model =
        (CameraModel)OF_Prefs::settings[OF_Const::cameraModel];

    if (!Select(model))
        return false;

    const uint8_t sensitivity = ClampSensitivity(
        (uint8_t)OF_Prefs::profiles[OF_Prefs::currentProfile].irSens
    );

    ready = activeOps->begin(sensitivity);
    return ready;
}
/*
bool OpenFIRECamera::Begin(Sensitivity_e sensitivity,
                           DataFormat_e format) {
    if (activeOps == nullptr || activeProfile == nullptr) {
        return false;
    }

    if (format != DataFormat_Basic && format != DataFormat_Extended) {
        return false;
    }

    activeFormat = format;
    activeRead = (format == DataFormat_Extended) ? activeOps->readExtended : activeOps->readBasic;
    ClearObjectData();

    ready = activeOps->begin(ClampSensitivity((uint8_t)sensitivity));
    return ready;
}
*/

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
            dfrCamera = new DFRobotIRPositionEx(Wire1);
        }
    } else if (!bitRead(pin_scl, 1) && !bitRead(pin_sda, 1)) {
        if (bitRead(pin_scl, 0) && !bitRead(pin_sda, 0)) {
            Wire.setSDA(pin_sda);
            Wire.setSCL(pin_scl);
            dfrCamera = new DFRobotIRPositionEx(Wire);
        }
    }

    if (pin_wiiClock > -1) {
        set_sys_clock_khz(125000, true);
        gpio_set_function(pin_wiiClock, GPIO_FUNC_PWM);
        const uint slice_num = pwm_gpio_to_slice_num(pin_wiiClock);
        pwm_set_clkdiv(slice_num, 1.0f);
        pwm_set_wrap(slice_num, 4);
        pwm_set_chan_level(slice_num, pwm_gpio_to_channel(pin_wiiClock), 2);
        pwm_set_enabled(slice_num, true);
    } else {
        set_sys_clock_khz(133000, true);
        //pwm_set_enabled(pwm_gpio_to_slice_num(pin_wiiClock), false);
        //gpio_set_function(pin_wiiClock, GPIO_FUNC_SIO);
        //gpio_put(pin_wiiClock, 0);
        //gpio_set_dir(pin_wiiClock, GPIO_IN);
    }
#endif

    if (dfrCamera == nullptr) {
        return false;
    }

    const DFRobotIRPositionEx::DataFormat_e format =
        (activeFormat == DataFormat_Extended)
            ? DFRobotIRPositionEx::DataFormat_Extended
            : DFRobotIRPositionEx::DataFormat_Basic;

    if (!dfrCamera->begin(activeProfile->busClock,
                          format,
                          (DFRobotIRPositionEx::Sensitivity_e)sensitivity)) {
        delete dfrCamera;
        dfrCamera = nullptr;
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

int OpenFIRECamera::ReadDFRobotExtended() {
    const int error = dfrCamera->extendedAtomic(DFRobotIRPositionEx::Retry_2);

    if (error >= DFRobotIRPositionEx::Error_Success) {
        activeSeen = dfrCamera->seen();

        for (int i = 0; i < 4; i++) {
            objectData[i] = {};

            if ((activeSeen & (1U << i)) != 0U) {
                objectData[i].valid = true;
                objectData[i].x = dfrCamera->x(i);
                objectData[i].y = dfrCamera->y(i);
                objectData[i].size = dfrCamera->size(i);
            }
        }
    }

    return error;
}

void OpenFIRECamera::DataFormatDFRobot(DataFormat_e format) {
    dfrCamera->dataFormat(
        (format == DataFormat_Extended)
            ? DFRobotIRPositionEx::DataFormat_Extended
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
}

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
    pajCamera = new OpenFIRE_PAJ7025(&pajSPI, pin_spiCs, *activeProfile);
#else
    SPI.setSCK(pin_spiSck);
    SPI.setMISO(pin_spiMiso);
    SPI.setMOSI(pin_spiMosi);
    SPI.setCS(pin_spiCs);
    SPI.begin();
    pajCamera = new OpenFIRE_PAJ7025(&SPI, pin_spiCs, *activeProfile);
#endif

    if (!pajCamera->begin(activeProfile->busClock, sensitivity)) {
        delete pajCamera;
        pajCamera = nullptr;
        return false;
    }

    activeX = pajCamera->xPositions();
    activeY = pajCamera->yPositions();
    activeSeen = pajCamera->seen();
    return true;
}

int OpenFIRECamera::ReadPAJ7025Basic() {
    const int error = pajCamera->readBasic();
    activeSeen = pajCamera->seen();
    return error;
}

int OpenFIRECamera::ReadPAJ7025Extended() {
    const int error = pajCamera->readExtended();
    activeSeen = pajCamera->seen();

    if (error >= Error_Success) {
        const PAJ7025_Object* source = pajCamera->objects();

        for (int i = 0; i < 4; i++) {
            objectData[i] = {};

            if ((activeSeen & (1U << i)) != 0U) {
                objectData[i].valid = source[i].is_valid;
                objectData[i].x = pajCamera->x(i);
                objectData[i].y = pajCamera->y(i);
                objectData[i].size = pajCamera->size(i);
                objectData[i].area = source[i].area;
                objectData[i].averageBrightness = source[i].average_brightness;
                objectData[i].maxBrightness = source[i].max_brightness;
                objectData[i].range = source[i].range;
                objectData[i].radius = source[i].radius;
                objectData[i].boundaryLeft = source[i].boundary_left;
                objectData[i].boundaryRight = source[i].boundary_right;
                objectData[i].boundaryUp = source[i].boundary_up;
                objectData[i].boundaryDown = source[i].boundary_down;
                objectData[i].aspectRatio = source[i].aspect_ratio;
                objectData[i].vx = source[i].vx;
                objectData[i].vy = source[i].vy;
            }
        }
    }

    return error;
}

void OpenFIRECamera::DataFormatPAJ7025(DataFormat_e format) {
    (void)format;
    // PAJ7025 format is selected directly by the bound read function.
}

void OpenFIRECamera::SensitivityPAJ7025(uint8_t sensitivity) {
    pajCamera->sensitivityLevel(sensitivity);
}

void OpenFIRECamera::EndPAJ7025() {
    if (pajCamera != nullptr) {
        delete pajCamera;
        pajCamera = nullptr;
    }
#ifdef ARDUINO_ARCH_ESP32
    pajSPI.end();
#endif
}
