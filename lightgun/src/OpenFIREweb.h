#pragma once
#include <stdint.h>
#include <stddef.h>

/// @brief Set when the user requested the web configuration mode at boot.
///        While true the lightgun also serves the App over WiFi (WebSocket):
///        an App can dock on the serial port or on the WebSocket, one at a time.
extern bool OF_WebConfigModeActive;

/// @brief Starts the access point, the HTTP server and the WebSocket (ESP32 only).
void WebApp_Init();
void WebApp_Loop();

/// @brief Returns and clears the notification that the WebSocket App client
///        disconnected or was replaced by a new page (always false on RP2040).
bool WebApp_TakeClientLost();

/// @brief True while a client-lost notification is waiting (does not clear it).
bool WebApp_ClientLostPending();

namespace OF_WebSerialWrapper {
    using AvailableFn = int (*)();
    using ReadFn = int (*)();
    using WriteByteFn = size_t (*)(uint8_t);
    using WriteBufFn = size_t (*)(const uint8_t*, size_t);
    using FlushFn = void (*)();
    using ReadBytesFn = size_t (*)(char*, size_t);

    struct SerialOps {
        AvailableFn available;
        ReadFn read;
        WriteByteFn writeByte;
        WriteBufFn writeBuf;
        FlushFn flush;
        ReadBytesFn readBytes;
    };
}

/// @brief Links that can carry the App protocol.
// (No member named "Serial": on ESP32 Serial is a macro.)
enum class WebAppLink : uint8_t {
    SerialPort, // USB serial port or wireless dongle stream
    WebSocket   // web configuration mode (ESP32 only)
};

/// @brief Byte stream used by the App protocol. It is the link on which the
///        current App docked (the serial link until an App docks on the WebSocket).
class WebAppSerial {
public:
    static OF_WebSerialWrapper::SerialOps ops;
    static WebAppLink link;

    /// @brief Selects the link used by the App protocol (called when an App docks).
    static void Use(WebAppLink newLink);
    /// @brief Direct access to one link, whatever link the App session uses.
    static const OF_WebSerialWrapper::SerialOps& Ops(WebAppLink which);
    static inline bool IsWebSocket() { return link == WebAppLink::WebSocket; }

    static inline int available() { return ops.available(); }
    static inline int read() { return ops.read(); }
    static inline size_t write(uint8_t c) { return ops.writeByte(c); }
    static inline size_t write(const uint8_t* buf, size_t size) { return ops.writeBuf(buf, size); }
    static inline void flush() { ops.flush(); }
    static inline size_t readBytes(char* buf, size_t size) { return ops.readBytes(buf, size); }
};
