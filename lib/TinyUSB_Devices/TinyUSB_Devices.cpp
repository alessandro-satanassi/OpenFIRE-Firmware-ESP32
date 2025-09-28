/*
 * This module (hastily) combines Mike Lynch's modified "AbsMouse5" library, and
 * Chris Young's TinyUSB Mouse and Keyboard library (the Keyboard half, anyways),
 * which in itself uses pieces of Arduino's basic Keyboard library.
 */
#ifdef USE_TINYUSB
  #include <Adafruit_TinyUSB.h>
#elifdef CFG_TUSB_MCU
  #error Incompatible USB stack. Use Adafruit TinyUSB.
#else
  #include <HID.h>
#endif

#include "TinyUSB_Devices.h"

#if defined(PIO_FRAMEWORK_ARDUINO_ENABLE_BLUETOOTH) && defined(ARDUINO_RASPBERRY_PI_PICO_W) && defined(ENABLE_CLASSIC)  
  #include <HID_Bluetooth.h>
  #include <PicoBluetoothHID.h>
#endif // ARDUINO_RASPBERRY_PI_PICO_W

/*****************************
 *   GLOBAL SECTION
 *****************************/
Adafruit_USBD_HID usbHid;

#ifdef COMMENTO
/* // 696969 spostato in TinyUSB_Devices
enum HID_RID_e{
  HID_RID_KEYBOARD = 1,
  HID_RID_MOUSE,   //2  
  HID_RID_GAMEPAD // 3  
};

uint8_t const desc_hid_report[] = {
  TUD_HID_REPORT_DESC_KEYBOARD(HID_REPORT_ID(HID_RID_e::HID_RID_KEYBOARD)),
  TUD_HID_REPORT_DESC_ABSMOUSE(HID_REPORT_ID(HID_RID_e::HID_RID_MOUSE)),
  TUD_HID_REPORT_DESC_GAMEPAD16(HID_REPORT_ID(HID_RID_e::HID_RID_GAMEPAD))
};
*/
#endif // COMMENTO

#if defined(PIO_FRAMEWORK_ARDUINO_ENABLE_BLUETOOTH) && defined(ARDUINO_RASPBERRY_PI_PICO_W) && defined(ENABLE_CLASSIC)
enum HID_BT_e {
    HID_BT_KEYBOARD = 1,
    HID_BT_MOUSE,  // 2
    HID_BT_GAMEPAD // 3
};

uint8_t desc_bt_report[] = {
    TUD_HID_REPORT_DESC_KEYBOARD(HID_REPORT_ID(HID_BT_e::HID_BT_KEYBOARD)),
    TUD_HID_REPORT_DESC_ABSMOUSE(HID_REPORT_ID(HID_BT_e::HID_BT_MOUSE)),
    TUD_HID_REPORT_DESC_GAMEPAD16(HID_REPORT_ID(HID_BT_e::HID_BT_GAMEPAD))
};
#endif // ARDUINO_RASPBERRY_PI_PICO_W

void TinyUSBDevices_::begin(int polRate) {
    usbHid.setPollInterval(polRate);
    usbHid.setReportDescriptor(desc_hid_report, sizeof(desc_hid_report));
    usbHid.begin();
    // 696969 aggiunto da me If already enumerated, additional class driverr begin() e.g msc, hid, midi won't take effect until re-enumeration
    if (TinyUSBDevice.mounted()) {
      TinyUSBDevice.detach();
      delay(10);
      TinyUSBDevice.attach();
    }
    onBattery = false;
}

#if defined(PIO_FRAMEWORK_ARDUINO_ENABLE_BLUETOOTH) && defined(ARDUINO_RASPBERRY_PI_PICO_W) && defined(ENABLE_CLASSIC)
void TinyUSBDevices_::beginBT(const char *localName, const char *hidName) {
    // third arg is the type of device that this is exposed as, i.e. the icon displayed on the PC.
    // for BLE: 0x03C2 is mouse, 0x03C1 is keyboard, 0x03C4 is gamepad, 0x03C0 is "generic" bluetooth icon
    // for BT: 0x2580 is mouse, 0x2540 is keyboard, 0x2508 is gamepad, 0x25C0 is "combo".
    // also bluetooth classic for some reason has a "subclass"?
    PicoBluetoothHID.startHID(localName, hidName, 0x2580, 33, desc_bt_report, sizeof(desc_bt_report));
    onBattery = true;
}
#endif // ARDUINO_RASPBERRY_PI_PICO_W

TinyUSBDevices_ TinyUSBDevices;
  
/*****************************
 *   MOUSE SECTION
 *****************************/

AbsMouse5_::AbsMouse5_() {}

void AbsMouse5_::report(void)
{
  if(!TinyUSBDevices.onBattery) {
    while(!usbHid.ready()) yield();  
    usbHid.sendReport(HID_RID_e::HID_RID_MOUSE, &absmouse5Report, sizeof(absmouse5Report));
  }
  #ifdef OPENFIRE_WIRELESS_ENABLE
    else {
    switch (TinyUSBDevices.wireless_mode)
    {
      #if defined(ENABLE_BLUETOOTH_OPENFIRE)
      case ENABLE_BLUETOOTH_TO_PC:
        BT_OpenFIRE.send(HID_BT_MOUSE, &absmouse5Report, sizeof(absmouse5Report));
      case ENABLE_BLUETOOTH_TO_DONGLE:
        /* code */
        break;
      #endif // ENABLE_BLUETOOTH_OPENFIRE
      #ifdef OPENFIRE_WIRELESS_DEVICE_ESPNOW
      case ENABLE_ESP_NOW_TO_DONGLE:  
        SerialWireless.SendPacket((const uint8_t *)&absmouse5Report, sizeof(absmouse5Report), PACKET_TX::MOUSE_TX);
        break;
      #endif //OPENFIRE_WIRELESS_DEVICE_ESPNOW
      case ENABLE_WIFI_TO_DONGLE:
        /* code */
        break;
      case ENABLE_NRF_24L01PLUS_TO_DONGLE:
        /* code */
        break;
      default:
        break;
    }  
  }
  #endif //OPENFIRE_WIRELESS_ENABLE
  TinyUSBDevices.newReport[TinyUSBDevices_::reportMouse] = false;
}

void AbsMouse5_::move_wheel_pan(int8_t wheel, int8_t pan) {
  if(pan != absmouse5Report.pan || wheel != absmouse5Report.wheel) {
		absmouse5Report.wheel = wheel;
    absmouse5Report.pan = pan;
		report();

  }
}

void AbsMouse5_::move(int16_t x, int16_t y)
{
	if(x != absmouse5Report.x || y != absmouse5Report.y) {
		absmouse5Report.x = x;
		absmouse5Report.y = y;
    TinyUSBDevices.newReport[TinyUSBDevices_::reportMouse] = true;
	  //report();
		
	}
}

void AbsMouse5_::press(uint8_t button)
{
	if(!(absmouse5Report.buttons & button)) {
    absmouse5Report.buttons |= button;
    TinyUSBDevices.newReport[TinyUSBDevices_::reportMouse] = true;
	  //report();
  }
}

void AbsMouse5_::release(uint8_t button)
{
	if(absmouse5Report.buttons & button) {
    absmouse5Report.buttons &= ~button;
	  TinyUSBDevices.newReport[TinyUSBDevices_::reportMouse] = true;
    //report();
  }
}

void AbsMouse5_::click(uint8_t b) {
  absmouse5Report.buttons = b;
  report();
  absmouse5Report.buttons = 0;
  report();
}

void AbsMouse5_::buttons(uint8_t b) {
  if (b != absmouse5Report.buttons) {
    absmouse5Report.buttons = b;
    report();
  }
}

bool AbsMouse5_::isPressed(uint8_t b) {
  if ((b & absmouse5Report.buttons) > 0) return true;
  else return false;
}

AbsMouse5_ AbsMouse5;


 /*****************************
 *   KEYBOARD SECTION
 *****************************/ 

Keyboard_::Keyboard_(void) {}
  

void Keyboard_::report(void) 
{
  if(!TinyUSBDevices.onBattery) {
    if (TinyUSBDevice.suspended())  { TinyUSBDevice.remoteWakeup(); }
    while(!usbHid.ready()) yield();   
    usbHid.sendReport(HID_RID_e::HID_RID_KEYBOARD,&_keyReport, sizeof(_keyReport));
  }
  #ifdef OPENFIRE_WIRELESS_ENABLE
  else {
    switch (TinyUSBDevices.wireless_mode)
    {
      #if defined(ENABLE_BLUETOOTH_OPENFIRE)
      case ENABLE_BLUETOOTH_TO_PC:
        BT_OpenFIRE.send(HID_BT_KEYBOARD, &_keyReport, sizeof(_keyReport));
        break;
      case ENABLE_BLUETOOTH_TO_DONGLE:
        /* code */
        break;
      #endif // ENABLE_BLUETOOTH_OPENFIRE
      #ifdef OPENFIRE_WIRELESS_DEVICE_ESPNOW
      case ENABLE_ESP_NOW_TO_DONGLE:
        SerialWireless.SendPacket((const uint8_t *)&_keyReport, sizeof(_keyReport), PACKET_TX::KEYBOARD_TX);
        break;
        #endif //OPENFIRE_WIRELESS_DEVICE_ESPNOW
      case ENABLE_WIFI_TO_DONGLE:
        /* code */
        break;
      case ENABLE_NRF_24L01PLUS_TO_DONGLE:
        /* code */
        break;
      default:
        break;
    }  
  }
  #endif //OPENFIRE_WIRELESS_ENABLE
  TinyUSBDevices.newReport[TinyUSBDevices_::reportKeyboard] = false;
}
 
size_t Keyboard_::press(uint8_t k)
{  
	uint8_t i;
	if (k >= 136) {			// it's a non-printing key (not a modifier)
		k = k - 136;
	} else if (k >= 128) {	// it's a modifier key
		_keyReport.modifier |= (1<<(k-128));
		k = 0;
	} else {				// it's a printing key
		//k = pgm_read_byte(_asciimap + k);
		k = pgm_read_byte(KeyboardLayout_it_IT + k);
    //k = conv_table[k][1];
    if (!k) {
			setWriteError();
			return 0;
		}
		if ((k & ALT_GR) == ALT_GR) {
			_keyReport.modifier |= 0x40;   // AltGr = right Alt
			k &= 0x3F; // k &= ~ALT_GR;
		} else if ((k & SHIFT) == SHIFT) {
			_keyReport.modifier |= 0x02;	// the left shift modifier
			k &= 0x7F; // k &= ~SHIFT;
		}
		if (k == ISO_REPLACEMENT) {
			k = ISO_KEY;
		}
	}

	
  if (k >= 0xE0 && k < 0xE8) {
    // it's a modifier key
    _keyReport.modifier |= (1 << (k - 0xE0));
  } else if (k && k < 0xA5) {
  // Add k to the key report only if it's not already present
	// and if there is an empty slot.
	  if (_keyReport.keycode[0] != k && _keyReport.keycode[1] != k &&
		  _keyReport.keycode[2] != k && _keyReport.keycode[3] != k &&
		  _keyReport.keycode[4] != k && _keyReport.keycode[5] != k) {

		  for (i=0; i<6; i++) {
			  if (_keyReport.keycode[i] == 0x00) {
				  _keyReport.keycode[i] = k;
				  break;
			  }
		  }
		  if (i == 6) {
			  setWriteError();
			  return 0;
		  }
	  }
  }  else if (_keyReport.modifier == 0) {
    //not a modifier and not a key
    return 0;
  }
  TinyUSBDevices.newReport[TinyUSBDevices_::reportKeyboard] = true;
	//report();
	return 1;
}
  
size_t Keyboard_::release(uint8_t k)
{
  uint8_t i;
	if (k >= 136) {			// it's a non-printing key (not a modifier)
		k = k - 136;
	} else if (k >= 128) {	// it's a modifier key
		_keyReport.modifier &= ~(1<<(k-128));
		k = 0;
	} else {				// it's a printing key
		//k = pgm_read_byte(_asciimap + k);
    k = pgm_read_byte(KeyboardLayout_it_IT + k);
    //k = conv_table[k][1];
		if (!k) {
			return 0;
		}
		if ((k & ALT_GR) == ALT_GR) {
			_keyReport.modifier &= ~(0x40);   // AltGr = right Alt
			k &= 0x3F; // k &= ~ALT_GR;
		} else if ((k & SHIFT) == SHIFT) {
			_keyReport.modifier &= ~(0x02);	// the left shift modifier
			k &= 0x7F; // k &= ~SHIFT;
		}
		if (k == ISO_REPLACEMENT) {
			k = ISO_KEY;
		}
	}

	if (k >= 0xE0 && k < 0xE8) {
    // it's a modifier key
    _keyReport.modifier &= ~(1 << (k - 0xE0));
  } else if (k && k < 0xA5) {
    // Test the key report to see if k is present.  Clear it if it exists.
    // Check all positions in case the key is present more than once (which it shouldn't be)  
    for (i=0; i<6; i++) {
		  if (0 != k && _keyReport.keycode[i] == k) {
			  _keyReport.keycode[i] = 0x00;
		  }
	  }
  }
  TinyUSBDevices.newReport[TinyUSBDevices_::reportKeyboard] = true;
  //report();
	return 1;
}
  
void Keyboard_::releaseAll(void)
{
  tu_memclr(&_keyReport.keycode, 6);
  _keyReport.modifier = 0;
  TinyUSBDevices.newReport[TinyUSBDevices_::reportKeyboard] = true;
  //report();
}
  
size_t Keyboard_::write(uint8_t c)
{
  uint8_t p = press(c);  // Keydown
  release(c);            // Keyup
  return p;              // just return the result of press() since release() almost always returns 1
}
  
size_t Keyboard_::write(const uint8_t *buffer, size_t size) 
{
  size_t n = 0;
  while (size--) {
    if (*buffer != '\r') {
      if (write(*buffer)) {
        n++;
      } else {
        break;
      }
    }
    buffer++;
  }
  return n;
}

  Keyboard_ Keyboard;  //create an instance of the Keyboard object

/*****************************
 *   GAMEPAD SECTION 16
 *****************************/

 Gamepad16_::Gamepad16_(void) {}

void Gamepad16_::moveCam(uint16_t origX, uint16_t origY) {
  if(stickRight) {
      gamepad16Report.x = map(origX, 0, 32767, -32767, 32767);
      gamepad16Report.y = map(origY, 0, 32767, -32767, 32767);
  } else {
      gamepad16Report.rx = map(origX, 0, 32767, -32767, 32767);
      gamepad16Report.ry = map(origY, 0, 32767, -32767, 32767);
  }
  TinyUSBDevices.newReport[TinyUSBDevices_::reportGamepad] = true;
  //report(); 
}

void Gamepad16_::moveStick(uint16_t origX, uint16_t origY) {
  // TODO: inverted output for Cabela's Top Shot Elite sticks, but it might be backwards for others.
  if(origX != _x || origY != _y) {
    _x = origX, _y = origY;
    if(stickRight) {
      gamepad16Report.rx = map(_x, 0, 4095, 32767, -32767);
      gamepad16Report.ry = map(_y, 0, 4095, 32767, -32767);
    } else {
      gamepad16Report.x = map(_x, 0, 4095, 32767, -32767);
      gamepad16Report.y = map(_y, 0, 4095, 32767, -32767);
    }
    TinyUSBDevices.newReport[TinyUSBDevices_::reportGamepad] = true;
    //report();
  }
}

void Gamepad16_::press(uint8_t buttonNum) {
  if(!(gamepad16Report.buttons & (1 << buttonNum))) {
    bitSet(gamepad16Report.buttons, buttonNum);
    TinyUSBDevices.newReport[TinyUSBDevices_::reportGamepad] = true;
    //report();
  }
}

void Gamepad16_::release(uint8_t buttonNum) {
  if(gamepad16Report.buttons & (1 << buttonNum)) {
    bitClear(gamepad16Report.buttons, buttonNum);
    TinyUSBDevices.newReport[TinyUSBDevices_::reportGamepad] = true;
    //report();
  }
}

void Gamepad16_::padUpdate(uint8_t padMask) {
  if(gamepad16Report.hat != padMask) {
    gamepad16Report.hat = padMask;
    TinyUSBDevices.newReport[TinyUSBDevices_::reportGamepad] = true;
    //report();
  }
}

void Gamepad16_::report()
{  
  if(!TinyUSBDevices.onBattery) {
    if (TinyUSBDevice.suspended())  { TinyUSBDevice.remoteWakeup(); }
    while(!usbHid.ready()) yield();   
    usbHid.sendReport(HID_RID_e::HID_RID_GAMEPAD,&gamepad16Report,sizeof(gamepad16Report));
  }
  #ifdef OPENFIRE_WIRELESS_ENABLE
  else {
    switch (TinyUSBDevices.wireless_mode)
    {
      #if defined(ENABLE_BLUETOOTH_OPENFIRE)
      case ENABLE_BLUETOOTH_TO_PC:
        BT_OpenFIRE.send(HID_BT_GAMEPAD, &gamepad16Report, sizeof(gamepad16Report));
        break;
      case ENABLE_BLUETOOTH_TO_DONGLE:
        /* code */
        break;
      #endif // ENABLE_BLUETOOTH_OPENFIRE
      #ifdef OPENFIRE_WIRELESS_DEVICE_ESPNOW
      case ENABLE_ESP_NOW_TO_DONGLE:
        SerialWireless.SendPacket((const uint8_t *)&gamepad16Report, sizeof(gamepad16Report), PACKET_TX::GAMEPADE_TX);
        break;
        #endif //OPENFIRE_WIRELESS_DEVICE_ESPNOW
      case ENABLE_WIFI_TO_DONGLE:
        /* code */
        break;
      case ENABLE_NRF_24L01PLUS_TO_DONGLE:
        /* code */
        break;
      default:
        break;
    }  
  }
  #endif // OPENFIRE_WIRELESS_ENABLE
  TinyUSBDevices.newReport[TinyUSBDevices_::reportGamepad] = false;
}

void Gamepad16_::releaseAll() {
  tu_memclr(&gamepad16Report, sizeof(gamepad16Report));
  _x = 2048, _y = 2048;
  TinyUSBDevices.newReport[TinyUSBDevices_::reportGamepad] = true;
  //report();
}


  Gamepad16_ Gamepad16;
