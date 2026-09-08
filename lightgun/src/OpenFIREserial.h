 /*!
 * @file OpenFIREserial.h
 * @brief Serial RX buffer reading routines.
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V2.0
 * @date 2026
 *
 * I thank you for producing the first original code:
 * 
 * @copyright That One Seong, 2025
 * @copyright GNU Lesser General Public License
 */ 

#ifndef _OPENFIRESERIAL_H_
#define _OPENFIRESERIAL_H_

#include <Arduino.h>
#include <unordered_map>
#include <string_view>
#include "OpenFIREDefines.h"

class OF_Serial
{
public:
    /// @brief    Method for processing the Serial buffer when docked to the Desktop App
    /// @details  Only method that allows for reading/writing to system settings.
    static void SerialProcessingDocked();

    static bool AppSerialSessionIsActive() { return appSerialSessionActive; }

    /// @brief    Sends an acknowledged response to the Desktop App.
    static bool AppSerialSendResponse(uint8_t command,
                                      const void *payload = nullptr,
                                      uint8_t length = 0,
                                      bool final = false);

    /// @brief    Sends a real-time event without acknowledgement or retry.
    static bool AppSerialSendEvent(uint8_t command, const void *payload = nullptr, uint8_t length = 0);

    /// @brief    Returns and clears a calibration-cancel request received from the App.
    static bool AppSerialTakeCalibrationCancel();

    /// @brief Send an existing error code (zero = fatal operation error).
    static void AppSerialSendError(uint8_t = 0);

    // Main routine that prints information to connected serial monitor when the gun enters Pause Mode.
    static void PrintResults();

    // utility function to wait for n bytes with timeout //[ESP32_PORT] inserita da me
    static bool Serial_available(uint8_t min = 1);

    #ifdef DEBUG_SERIAL
    static void PrintDebugSerial();
    #endif // DEBUG_SERIAL

    #ifdef MAMEHOOKER
    /// @brief    Main method processing the Serial buffer.
    static void SerialProcessing();

    /// @brief    Handling gun events that may have been processed in SerialProcessing
    static void SerialHandling();

    // For serial mode:
    enum SerialQueueBits {
        SerialQueue_Solenoid = 0,
        SerialQueue_SolPulse,
        SerialQueue_Rumble,
        SerialQueue_RumbPulse,
        SerialQueue_Red,
        SerialQueue_Green,
        SerialQueue_Blue,
        SerialQueue_LEDPulse,
        SerialQueueBitsCount
    };

    static inline bool serialMode = false;                         // Set if we're prioritizing force feedback over serial commands or not.
    static inline bool serialQueue[SerialQueueBitsCount] = {false};// Array of events we've queued from the serial receipt.
    static inline bool serialARcorrection = false;                 // 4:3 AR correction mode flag
    static inline bool serialMappingsOffscreenShot = false;        // Marker if Offscreen Shot Mode's been enabled, for FW_Common::UpdateBindings
    static inline int  serialMappingsPedalMode = 0;                // Marker if Pedal has been remapped, for FW_Common::UpdateBindings
    // from least to most significant bit: solenoid digital, solenoid pulse, rumble digital, rumble pulse, R/G/B direct, RGB (any) pulse.

    // These do get addressed by the main code
    #ifdef USES_DISPLAY
    static inline bool serialDisplayChange = false;                // Signal of pending display update, sent by Core 2 to be used by Core 1 in dual core configs
    static inline uint serialLifeCount = 0;		                   // Changed from uint8_t for games with life values > 255
    static inline uint serialAmmoCount = 0;
    #endif // USES_DISPLAY

    #endif // MAMEHOOKER

private:
    // Desktop App serial framing. The two start bytes are used only to find
    // frame boundaries and are deliberately excluded from the CRC.
    static constexpr uint8_t  APP_SERIAL_START_1       = 0xA5;
    static constexpr uint8_t  APP_SERIAL_START_2       = 0x5A;
    static constexpr uint8_t  APP_SERIAL_MAX_PAYLOAD   = 200;
    static constexpr uint8_t  APP_SERIAL_OVERHEAD      = 7;
    static constexpr uint16_t APP_SERIAL_MAX_FRAME     = APP_SERIAL_MAX_PAYLOAD + APP_SERIAL_OVERHEAD;
    static constexpr uint16_t APP_SERIAL_FRAME_TIMEOUT = 250;
    static constexpr uint16_t APP_SERIAL_ACK_TIMEOUT   = 500;
    // Additional attempts after the initial transmission.
    static constexpr uint8_t  APP_SERIAL_MAX_RETRIES   = 3;

    // sError payload: [retry required] [originating request command].
    static constexpr uint8_t APP_SERIAL_ERR_COMMIT_RETRY = 0x82;

    enum AppSerialType_e : uint8_t {
        APP_SERIAL_TYPE_REQUEST  = 0x00,
        APP_SERIAL_TYPE_RESPONSE = 0x01,
        APP_SERIAL_TYPE_EVENT    = 0x02,
        APP_SERIAL_TYPE_ACK      = 0x03,
        APP_SERIAL_TYPE_MASK     = 0x03,
        APP_SERIAL_FLAG_FINAL    = 0x80
    };

    typedef struct AppSerialFrame_t {
        uint8_t typeFlags;
        uint8_t command;
        uint8_t sequence;
        uint8_t length;
        uint8_t payload[APP_SERIAL_MAX_PAYLOAD];
        uint8_t crc;
    } AppSerialFrame_s;

    static uint8_t AppSerialCRC8(const uint8_t *data, uint16_t length);
    static bool AppSerialWriteFrame(uint8_t typeFlags, uint8_t command,
                                    uint8_t sequence, const void *payload,
                                    uint8_t length);
    static bool AppSerialSendReliable(uint8_t typeFlags, uint8_t command,
                                      const void *payload, uint8_t length,
                                      uint8_t sequence = 0);
    static bool AppSerialReadFrame(AppSerialFrame_s &frame);
    static void AppSerialHandleFrame(const AppSerialFrame_s &frame);
    static bool AppSerialRequestMatchesLast(const AppSerialFrame_s &frame);
    static void AppSerialProcessDeferredRequest();
    static void AppSerialDispatchRequest(const AppSerialFrame_s &frame);
    // Returns true when the commit state machine consumed the request.
    static bool AppSerialDispatchCommit(const AppSerialFrame_s &frame,
                                        bool fullDisconnect);
    static void AppSerialRestoreRunState();
    static void AppSerialSendAck(uint8_t command, uint8_t sequence);
    static void AppSerialSendCommitError(const AppSerialFrame_s &frame);
    static void AppSerialSessionBegin();
    static void AppSerialSessionEnd();
    static uint8_t AppSerialNextSequence();

    static bool AppSerialSendRecords(
        uint8_t command, const void *data,
        const std::unordered_map<std::string_view, int> &fields,
        size_t dataSize, int profile = -1, bool sendFinal = true);
    static bool AppSerialReceiveRecord(
        const uint8_t *payload, uint8_t length, void *data,
        const std::unordered_map<std::string_view, int> &fields,
        size_t dataSize, bool profileData = false);

    // Frame buffers and receive-parser state.
    static inline uint8_t appSerialRxBuffer[APP_SERIAL_MAX_FRAME];
    static inline uint8_t appSerialTxBuffer[APP_SERIAL_MAX_FRAME];
    static inline uint16_t appSerialRxLength = 0;
    static inline unsigned long appSerialRxTimestamp = 0;

    // Session and outgoing-response state.
    static inline uint8_t appSerialTxSequence = 0;
    static inline uint8_t appSerialRawDockState = 0;
    static inline bool appSerialSessionActive = false;

    // A request received while a reliable response is waiting for its ACK
    // is retained and processed after the current transaction has completed.
    static inline bool appSerialWaitingForAck = false;
    static inline bool appSerialDispatching = false;
    static inline bool appSerialProcessingDeferred = false;
    static inline bool appSerialDeferredValid = false;
    static inline AppSerialFrame_s appSerialDeferredFrame = {};

    // Identity of the last request, used to acknowledge retries without
    // executing the same command twice.
    static inline bool appSerialLastRxValid = false;
    static inline uint8_t appSerialLastRxCommand = 0;
    static inline uint8_t appSerialLastRxSequence = 0;
    static inline uint8_t appSerialLastRxLength = 0;
    static inline uint8_t appSerialLastRxCRC = 0;

    // Configuration transaction and calibration state.
    static inline bool appSerialCommitActive = false;
    static inline bool appSerialCommitFailed = false;
    static inline bool appSerialCalibrationCancel = false;

    #ifdef MAMEHOOKER

    #ifdef LED_ENABLE
    static inline unsigned long serialLEDPulsesLastUpdate = 0;     // The timestamp of the last serial-invoked LED pulse update we iterated.
    static inline unsigned int serialLEDPulsesLength = 2;          // How long each stage of a serial-invoked pulse rumble is, in ms.
    static inline bool serialLEDChange = false;                    // Set on if we set an LED command this cycle.
    static inline bool serialLEDPulseRising = true;                // In LED pulse events, is it rising now? True to indicate rising, false to indicate falling; default to on for very first pulse.
    static inline uint serialLEDPulses = 0;                        // How many LED pulses are we being told to do?
    static inline uint serialLEDPulsesLast = 0;                    // What LED pulse we've processed last.
    static inline uint8_t serialLEDR = 0;                          // For the LED, how strong should it be?
    static inline uint8_t serialLEDG = 0;                          // Each channel is defined as three brightness values
    static inline uint8_t serialLEDB = 0;                          // So yeah.
    static inline uint8_t serialLEDPulseColorMap = 0b00000000;     // The map of what LEDs should be pulsing (we use the rightmost three of this bitmask for R, G, or B).
    #endif // LED_ENABLE

    #ifdef USES_RUMBLE
    static inline unsigned long serialRumbPulsesLastUpdate = 0;    // The timestamp of the last serial-invoked pulse rumble we updated.
    static constexpr uint serialRumbPulsesLength = 60;             // How long each stage of a serial-invoked pulse rumble is, in ms.
    static inline uint serialRumbPulseStage = 0;                   // 0 = start/rising, 1 = peak, 2 = falling, 3 = final check/reset to start
    static inline uint serialRumbPulses = 0;                       // If rumble is commanded to do pulse responses, how many?
    static inline uint serialRumbPulsesLast = 0;                   // Counter of how many pulse rumbles we did so far.
    static inline uint serialRumbCustomHoldLength = 0;             // Determines custom solenoid ON state length for sol "pulse" commands - 0 = use system settings
    static inline uint serialRumbCustomPauseLength = 0;            // Determines custom solenoid OFF state length for sol "pulse" commands - 0 = use system settings
    #endif // USES_RUMBLE

    #ifdef USES_SOLENOID
    static inline unsigned long serialSolPulsesLastUpdate = 0;     // The timestamp of the last serial-invoked pulse solenoid event we updated.
    static inline uint serialSolPulses = 0;                        // How many solenoid pulses are we being told to do?
    static inline uint serialSolPulsesLast = 0;                    // What solenoid pulse we've processed last.
    static inline unsigned long serialSolTimestamp = 0;            // Timestamp of the last solenoid static on command (for safety)
    static inline uint serialSolCustomHoldLength = 0;              // Determines custom solenoid ON state length for sol "pulse" commands - 0 = use system settings
    static inline uint serialSolCustomPauseLength = 0;             // Determines custom solenoid OFF state length for sol "pulse" commands - 0 = use system settings
    #ifdef USES_TEMP
    // When tempStatus is above Temp_Safe, new static solenoid ON commands toggles this.
    // False = disable solenoid ON command, True = allow solenoid ON command
    static inline bool serialSolTempBuffer = false;
    #endif // USES_TEMP
    #define SERIAL_SOLENOID_MAXSHUTOFF 2000
    #endif // USES_SOLENOID

    #endif // MAMEHOOKER

    //// Printing
    // used for periodic serial prints
    static inline unsigned long lastPrintMillis = 0;

    // used for debug prints
    #ifdef DEBUG_SERIAL
    static inline unsigned long serialDbMs = 0;
    static inline unsigned long frameCount = 0;
    static inline unsigned long irPosCount = 0;
    #endif
};

#endif // _OPENFIRESERIAL_H_
