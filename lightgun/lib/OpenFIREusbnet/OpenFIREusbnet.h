#pragma once

// Boot-only USB selection: normal HID + CDC, or web configuration HID + NCM.
// Wi-Fi, the dongle's serial stream and the App protocol remain independent.
#if defined(ARDUINO_ARCH_ESP32) && defined(OPENFIRE_USB_NCM)
#include <Arduino.h>

// Call once after CheckBootRequests(), after setting the USB identity.
// Adds all interfaces, then starts USB. A failed NCM preparation falls back
// to CDC; OpenFIREUsbNetActive() reports which configuration was selected.
bool OpenFIREUsbBegin(bool webConfig, uint8_t pollRate);
Stream& OpenFIREUsbSerial(); // Real CDC, or a nonblocking empty stream in web mode.
bool OpenFIREUsbNetActive();
void OpenFIREUsbNetMDNSReady(); // Call after the shared MDNS.begin() succeeds.
void OpenFIREUsbDetach();   // Main task: dongle fallback or before ROM bootloader.
#else
static inline bool OpenFIREUsbNetActive() { return false; }
#endif
