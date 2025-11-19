/*
 * This module simulates the standard Arduino "Mouse.h" and
 * "Keyboard.h" API for use with the TinyUSB HID API. Instead of doing
 *  #include <HID.h>
 *  #include <Mouse.h>
 *  #include <Keyboard.h>
 *  
 *  Simply do
 *  
 *  #include <TinyUSB_Mouse_Keyboard.h>
 *  
 *  and this module will automatically select whether or not to use the
 *  standard Arduino mouse and keyboard API or the TinyUSB API. We had to
 *  combine them into a single library because of the way TinyUSB handles
 *  descriptors.
 *  
 *  For details on Arduino Mouse.h see
 *   https://www.arduino.cc/reference/en/language/functions/usb/mouse/
 *  For details on Arduino Keyboard.h see
 *   https://www.arduino.cc/reference/en/language/functions/usb/keyboard/
 *
 *  NOTE: This code is derived from the standard Arduino Mouse.h, Mouse.cpp,
 *    Keyboard.h, and Keyboard.cpp code. The copyright on that original code
 *    is as follows.
 *   
 *  Copyright (c) 2015, Arduino LLC
 *  Original code (pre-library): Copyright (c) 2011, Peter Barrett
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#pragma once

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)

#include <Arduino.h>


#include <Adafruit_TinyUSB.h> // serve per alcune intestazioni
#include "TinyUSB_Devices.h"
#include "../OpenFIRE_Packet/OpenFIRE_Packet.h"

#if defined(ARDUINO_ARCH_ESP32)
  #include <WiFi.h> // serve per leggere max addres
  #include <esp_now.h>
  //#include <WiFi.h> // serve per leggere max addres
  #include <esp_mac.h>  // For the MAC2STR and MACSTR macros
  #include "MacAddress.h"
  #include "esp_wifi.h"
  #include "freertos/semphr.h"
  #include "esp_timer.h"
#elif defined(ARDUINO_ARCH_RP2040)
  // vediamo
#endif

/*****************************
 *   ESP_NOW
 *****************************/

 
extern uint8_t espnow_wifi_channel;
extern uint8_t espnow_wifi_power;
extern uint8_t lastDongleAddress[6];
extern uint8_t lastDongleChannel;
extern uint8_t peerAddress[6];
extern bool lastDongleSave;
 
 

enum WIRELESS_MODE {
  NONE_WIRELESS =  0,                
  ENABLE_BLUETOOTH_TO_PC,         //1
  ENABLE_BLUETOOTH_TO_DONGLE,     //2
  ENABLE_ESP_NOW_TO_DONGLE,       //3
  ENABLE_WIFI_TO_DONGLE,          //4
  ENABLE_NRF_24L01PLUS_TO_DONGLE  //5
};

enum PACKET_TX {
  SERIAL_TX = 1, // DATI SERIALI
  KEYBOARD_TX,
  MOUSE_TX,
  GAMEPADE_TX,
  CONNECTION,  // CONNESSIONE E ASSOCIAZIONE DONGLE CON GUN
  CHECK_CONNECTION_LAST_DONGLE, //VERIFICA LA CONNESSIONE WIRELESS TRA GUN E DONGLE E VICEVERSA
  // ========= VALUTARE
  RESET_DATA_USB,
  REBOOT,      // REBOOT DEL DONGLE
  TX_DATA_USB, // TRASMISSIONE DATI DELL'USB DA GUN A DONGLE
  CHANGE_DATA_USB, // CAMBIO DEI DATI DELL'USB, QUINDI ?????
  CONN_REQ_TO_GUN,      // RICHIESTA INVIATA DA DONGLE VERSO GUN (INVIA MAC ADDRESS DEL DONGLE ed il MAC ADDRES A CUI E' STATA INVIATA LA RICHIESTA (PUO' ESSERE BROADCAST))
  CONN_REQ_TO_DONGLE,   // RICHIESTA INVIATA DA GUN A DONGLE, PER ACCETTARE RICHIESTA DI CONNESSIONE (INVIA MAC ADDRESS DEL GUN E MAC ADDRES A CUI E' STATA INVIATA RICHIESTA)
  CONN_ACCEPT,          // IL DONGLE HA ACCETTATO CONNESSIONE (INVIA ANCHE MAC ADDRES DEL DONGLE E DELLA GUN A CUI E' STATA ACCETTATA LA RICHIESTA)
};

enum CONNECTION_STATE {
  NONE_CONNECTION = 0,
  TX_DONGLE_SEARCH_GUN_BROADCAST,  //1
  TX_GUN_TO_DONGLE_PRESENCE,  //2
  TX_DONGLE_TO_GUN_ACCEPT,  //3
  TX_GUN_TO_DONGLE_CONFERM,  //4
  DEVICES_CONNECTED,  //5
  // =============================
  // === per controllo connessione
  NONE_CONNECTION_DONGLE,
  TX_CHECK_CONNECTION_LAST_DONGLE, //6  GUN -> DONGLE
  TX_CONFERM_CONNECTION_LAST_DONGLE, //7 DONGLE -> GUN
  DEVICES_CONNECTED_WITH_LAST_DONGLE, // 8
};

typedef struct __attribute__ ((packed)) {
  char deviceManufacturer[21];
  char deviceName[16];
  uint16_t deviceVID;
  uint16_t devicePID;
  uint8_t devicePlayer;
  uint8_t channel; // da 0 a 13
  //char deviceSerial[11];
} USB_Data_GUN_Wireless;


extern USB_Data_GUN_Wireless usb_data_wireless;


/*****************************
 *   SERIAL WIRELESS SECTION
 *****************************/
#ifndef _SERIAL_STREAM_H_
#define _SERIAL_STREAM_H_

class SerialWireless_ : public Stream       
  {
    
  public:
  
  uint8_t mac_esp_inteface[6];
  uint8_t mac_esp_another_card[6];
  uint8_t stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION; //0;

  // ======= per FIFO SERIAL ===============
  // ===== per write === buffer lineare ====
  #define FIFO_SIZE_WRITE_SERIAL 128
  #define TIME_OUT_SERIAL_WRITE 3  // in seguito rimuovere o rideterminare in microsecondi
  #define TIMER_HANDLE_SERIAL_DURATION_MICROS 3000 // 3000 microsecondi un 3 millisecondi
  uint8_t bufferSerialWrite[FIFO_SIZE_WRITE_SERIAL];
  unsigned long startTimeSerialWrite = 0; // = millis(); // poi convertire in microsecondi // in seguito rimovere
  volatile uint16_t lenBufferSerialWrite = 0;
  esp_timer_handle_t timer_handle_serial;
  // ====== per read ====== buffer circolare =====
  #define FIFO_SIZE_READ_SERIAL 200
  uint8_t bufferSerialRead[FIFO_SIZE_READ_SERIAL];
  volatile uint16_t lenBufferSerialRead = 0;
  volatile uint16_t _writerSerialRead = 0;
  volatile uint16_t _readerSerialRead = 0;
  bool _overflow_bufferSerialRead = false; 
  void write_on_rx_serialBuffer(const uint8_t *data, int len);
  // ============================================

  Packet  packet;
  
  // per buffer lettura
  volatile uint16_t _readLen = 0;
  volatile uint16_t _writer = 0;
  volatile uint16_t _reader = 0;
  #define FIFO_SIZE_READ 1024 // buffer lettura
  uint8_t _queue[FIFO_SIZE_READ];
  bool _overflow_read = false; 
  // fine per buffer lettura

  // per buffer scrittura 
  volatile uint16_t writeIndex = 0;
  volatile uint16_t readIndex = 0;
  volatile uint16_t _writeLen = 0;
  #define BUFFER_SIZE 1024 //buffer scrittura
  uint8_t buffer[BUFFER_SIZE];
  bool _overflow_write = false; 
  // fine per buffer scrittura
 
  SerialWireless_() : Stream() {}
  ~SerialWireless_() {}
  
  void init_wireless();
  void begin();
  bool end();
  
  // overrade da ::Stream
  int peek() override;
  int read() override;
  int available() override;
  // fine override da ::stream
  
  // overraid da ::Print
  int availableForWrite() override;
  void flush() override;
  size_t write(uint8_t c) override;
  size_t write(const uint8_t *data, size_t len) override;
  using Print::write;
  // ==== fine override da :: Print

  operator bool();

  // inserire da me per gestione buffer uscita
  size_t writeBin(uint8_t c);
  size_t writeBin(const uint8_t *data, size_t len);
  void flushBin();
  int availableForWriteBin();
  bool flush_sem();
  // inserire da me per gestione buffer ingresso
  int peekBin();
  int readBin();
  int availableBin();

  //inserito per gestire Packet
  volatile uint16_t numAvailablePacket = 0;
  int availablePacket();
        
  // ======== generiche ============
  void SendData();  // utilizziamo anche flush
  void SendData_sem();
  void SendPacket(const uint8_t *data, const uint8_t &len,const uint8_t &packetID); // non penso lo utilizzeremo

  bool checkForRxPacket(); // andrà nel ciclo principale .. messo nella callback

  bool connection_dongle();
  bool connection_gun();
  bool connection_gun_at_last_dongle();

  // ===============================
  // ===== per i timer ================

//void setupTimer(uint64_t duration_us);
void setupTimer();
void stopTimer_serial();
void resetTimer_serial(uint64_t duration_us);

private:

};
extern SerialWireless_ SerialWireless;

#endif // _SERIAL_STREAM_H_

#endif //OPENFIRE_WIRELESS_ENABLE
