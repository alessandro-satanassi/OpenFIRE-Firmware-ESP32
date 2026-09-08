#pragma once
#include <stdint.h>
#include <stddef.h>

extern bool OF_WebConfigModeActive;

void WebApp_Init();
void WebApp_Loop();
void WebApp_SendToApp(const uint8_t* buffer, size_t length);
