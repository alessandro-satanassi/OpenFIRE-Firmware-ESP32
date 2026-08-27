/*!
 * @file OpenFIRECamera.h
 * @brief Runtime camera manager and common camera API for OpenFIRE.
 *
 * @copyright Alessandro Satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author Alessandro Satanassi
 * @version V1.0
 * @date 2026
 */

#ifndef OPENFIRE_CAMERA_H
#define OPENFIRE_CAMERA_H

#include <stdint.h>
#include "OpenFIRECameraProfile.h"

struct OpenFIRECameraPins {
    int8_t sda;
    int8_t scl;
    int8_t wiiClock;
    int8_t spiSck;
    int8_t spiMiso;
    int8_t spiMosi;
    int8_t spiCs;
};

class OpenFIRECamera {
public:
    // Common data format exposed to the OpenFIRE firmware.
    enum DataFormat_e : uint8_t {
        DataFormat_Basic = 0,
        DataFormat_Extended = 1
    };

    // Common sensitivity levels exposed to the OpenFIRE firmware.
    enum Sensitivity_e : uint8_t {
        Sensitivity_Min = 0,
        Sensitivity_Default = 0,
        Sensitivity_High = 1,
        Sensitivity_Max = 2
    };

    // Common camera read results exposed to the OpenFIRE firmware.
    enum Errors_e : int8_t {
        Error_SuccessMismatch = 1,
        Error_Success = 0,
        Error_Communication = -1,
        Error_IICerror = Error_Communication, // Legacy-compatible alias.
        Error_DataMismatch = -2,
    };

    // Fields available through ObjectData when Extended format is active.
    // Unsupported camera-specific fields are returned as zero.
    enum ExtendedData_e : uint16_t {
        ExtendedData_Size              = 1U << 0,
        ExtendedData_Area              = 1U << 1,
        ExtendedData_AverageBrightness = 1U << 2,
        ExtendedData_MaxBrightness     = 1U << 3,
        ExtendedData_Range             = 1U << 4,
        ExtendedData_Radius            = 1U << 5,
        ExtendedData_Boundaries        = 1U << 6,
        ExtendedData_AspectRatio       = 1U << 7,
        ExtendedData_Velocity          = 1U << 8
    };

    // Common extended object data. X/Y/valid/size are normalized to the
    // selected backend's OpenFIRE output slots; the remaining fields retain
    // the native values supplied by cameras that support them.
    struct ObjectData {
        int x;
        int y;
        int size;
        uint16_t area;
        bool valid;
        uint8_t averageBrightness;
        uint8_t maxBrightness;
        uint8_t range;
        uint8_t radius;
        uint8_t boundaryLeft;
        uint8_t boundaryRight;
        uint8_t boundaryUp;
        uint8_t boundaryDown;
        uint8_t aspectRatio;
        uint8_t vx;
        uint8_t vy;
    };

    static bool Select(CameraModel model);
    static const CameraProfile& Profile();
    static CameraModel Model();

    static bool Begin(const OpenFIRECameraPins& pins,
                      Sensitivity_e sensitivity,
                      DataFormat_e format = DataFormat_Basic);
    static void End();

    // Hot path: activeRead is bound to the selected backend and data format.
    static inline int Read() { return activeRead(); }

    static bool SetDataFormat(DataFormat_e format);
    static DataFormat_e DataFormat() { return activeFormat; }
    static void SetSensitivity(Sensitivity_e sensitivity);

    static bool IsReady() { return ready; }
    static const int* XPositions() { return activeX; }
    static const int* YPositions() { return activeY; }
    static unsigned int Seen() { return activeSeen; }

    // Extended data is updated only by successful Extended reads. Switching
    // back to Basic clears this buffer so stale Extended values are never
    // exposed as current data.
    static const ObjectData* Objects() { return objectData; }
    static const ObjectData& Object(uint8_t index) { return objectData[index]; }
    static uint16_t ExtendedCapabilities() { return activeExtendedCapabilities; }
    static bool SupportsExtendedData(ExtendedData_e field) {
        return (activeExtendedCapabilities & (uint16_t)field) != 0U;
    }

private:
    using BeginFn = bool (*)(const OpenFIRECameraPins&, uint8_t);
    using ReadFn = int (*)();
    using DataFormatFn = void (*)(DataFormat_e);
    using SensitivityFn = void (*)(uint8_t);
    using EndFn = void (*)();

    struct CameraOps {
        BeginFn begin;
        ReadFn readBasic;
        ReadFn readExtended;
        DataFormatFn dataFormat;
        SensitivityFn sensitivity;
        EndFn end;
        uint16_t extendedCapabilities;
    };

    static const CameraProfile* activeProfile;
    static const CameraOps* activeOps;
    static ReadFn activeRead;
    static const int* activeX;
    static const int* activeY;
    static unsigned int activeSeen;
    static ObjectData objectData[4];
    static uint16_t activeExtendedCapabilities;
    static DataFormat_e activeFormat;
    static bool ready;

    static void ClearObjectData();
    static inline uint8_t ClampSensitivity(uint8_t sensitivity) {
        return (sensitivity > (uint8_t)Sensitivity_Max) ? (uint8_t)Sensitivity_Max : sensitivity;
    }


    // DFRobot SEN0158 / WiiCam
    static bool BeginDFRobot(const OpenFIRECameraPins& pins, uint8_t sensitivity);
    static int ReadDFRobotBasic();
    static int ReadDFRobotExtended();
    static void DataFormatDFRobot(DataFormat_e format);
    static void SensitivityDFRobot(uint8_t sensitivity);
    static void EndDFRobot();

    // PixArt PAJ7025R2 / PixArt PAJ7025R3
    static bool BeginPAJ7025(const OpenFIRECameraPins& pins, uint8_t sensitivity);
    static int ReadPAJ7025Basic();
    static int ReadPAJ7025Extended();
    static void DataFormatPAJ7025(DataFormat_e format);
    static void SensitivityPAJ7025(uint8_t sensitivity);
    static void EndPAJ7025();
};

#endif // OPENFIRE_CAMERA_H
