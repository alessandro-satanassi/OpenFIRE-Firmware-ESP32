#pragma once
#include <stdint.h>
#include <stddef.h>

extern bool OF_WebConfigModeActive;

void WebApp_Init();
void WebApp_Loop();
void WebApp_SendToApp(const uint8_t* buffer, size_t length);

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

class WebAppSerial {
public:
    static OF_WebSerialWrapper::SerialOps ops;

    static inline int available() { return ops.available(); }
    static inline int read() { return ops.read(); }
    static inline size_t write(uint8_t c) { return ops.writeByte(c); }
    static inline size_t write(const uint8_t* buf, size_t size) { return ops.writeBuf(buf, size); }
    static inline void flush() { ops.flush(); }
    static inline size_t readBytes(char* buf, size_t size) { return ops.readBytes(buf, size); }
};

