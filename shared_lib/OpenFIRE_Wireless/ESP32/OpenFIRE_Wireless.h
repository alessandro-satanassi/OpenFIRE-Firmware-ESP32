/*!
 * @file OpenFIRE_Wireless.h
 * @brief Library for OpenFIRE wireless with Esp-Now
 * @n CPP Library for OpenFIRE wireless with Esp-Now
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#pragma once

// ===================================================================================
// ASTRAZIONE MULTI-PIATTAFORMA (WIP) E DIPENDENZE
// ===================================================================================
// L'obiettivo è mantenere l'hardware sottostante il più possibile astratto,
// permettendo in futuro porting su altre architetture (es. RP2040) se le librerie
// radio lo consentiranno. Attualmente la rete è fortemente legata all'ecosistema 
// Espressif (ESP-NOW) per questioni di latenza ultra-bassa.

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)

#include <Arduino.h>


#include <Adafruit_TinyUSB.h> // serve per alcune intestazioni
#include "TinyUSB_Devices.h"
#include "../OpenFIRE_Packet/OpenFIRE_Packet.h"

#if defined(ARDUINO_ARCH_ESP32)
  #include <WiFi.h> // serve per leggere max addres
  #include <esp_now.h>
  // #include <WiFi.h> // serve per leggere max addres
  #include <esp_mac.h>  // For the MAC2STR and MACSTR macros
  #include "MacAddress.h"
  #include "esp_wifi.h"
  #include "freertos/semphr.h"
  #include "esp_timer.h"
#elif defined(ARDUINO_ARCH_RP2040)
  // vediamo
#endif

// numero di secondi massimi a disposizione per la connessione del Pedal
#define MAX_SECONDI_CONNESSIONE_PEDAL 10

/*****************************
 *   ESP_NOW E VARIABILI DI STATO GLOBALI
 *****************************/
// ===================================================================================
// GESTIONE DELLO STATO DISTRIBUITO (EXTERN)
// ===================================================================================
// L'uso estensivo di extern permette di condividere lo stato della rete (MAC address, 
// canali, report USB) tra diversi moduli senza overhead di passaggio parametri. 
// Le struct report (es. absmouse5Report_last_wifi) fungono da "Shadow State" per 
// inviare i dati all'host USB solo quando c'è un reale cambiamento hardware.

extern uint8_t espnow_wifi_channel;
extern uint8_t espnow_wifi_power;
extern hid_abs_mouse_report_t absmouse5Report_last_wifi;
extern hid_keyboard_report_t  keyReport_last_wifi;
extern hid_gamepad16_report_t gamepad16Report_last_wifi;
extern uint8_t lastDongleAddress[6];
extern uint8_t lastPedalAddress[6];
extern uint8_t lastDongleChannel;
extern uint8_t lastPedalChannel;
extern uint8_t peerAddress[6];
extern bool lastDongleSave;
extern bool lastPedalSave; 
extern const uint8_t BROADCAST_ADDR[6];

extern uint8_t peerAddress_pedal[6];
 

// ===================================================================================
// DEFINIZIONI DEL PROTOCOLLO (ENUMERATORI) E MACCHINA A STATI DI RETE
// ===================================================================================
// Definiscono il vocabolario del nostro livello di trasporto proprietario.
// PACKET_TX: Identifica il TIPO di payload trasportato (dati USB, force feedback seriale, o comandi di rete).
// CONNECTION_STATE: Una macchina a stati finiti (FSM) rigorosa che guida la fase di "Handshake" 
// (ricerca, associazione, conferma) garantendo che i dispositivi (Gun, Dongle, Pedal) non si incrocino.


enum WIRELESS_MODE {
  NONE_WIRELESS = 0,   // WIFI SPENTO
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
  PEDAL_TX,
  CONNECTION,  // CONNESSIONE E ASSOCIAZIONE DONGLE CON GUN
  CHECK_CONNECTION_LAST_DONGLE, //VERIFICA LA CONNESSIONE WIRELESS TRA GUN E DONGLE E VICEVERSA
  CONNECTION_PEDAL,  // CONNESSIONE E ASSOCIAZIONE PEDAL CON GUN
  CHECK_CONNECTION_LAST_PEDAL, //VERIFICA LA CONNESSIONE WIRELESS TRA GUN E PEDAL E VICEVERSA
  #ifdef OPENFIRE_USE_ESPNOW_UNIFIED_PACKET
  MOUSE_KEY_PAD_TX,
  #endif // OPENFIRE_USE_ESPNOW_UNIFIED_PACKET
  // ========= VALUTARE
  RESET_DATA_USB,
  REBOOT,      // REBOOT DEL DONGLE
  TX_DATA_USB, // TRASMISSIONE DATI DELL'USB DA GUN A DONGLE
  CHANGE_DATA_USB, // CAMBIO DEI DATI DELL'USB
  CONN_REQ_TO_GUN,      // RICHIESTA INVIATA DA DONGLE VERSO GUN (INVIA MAC ADDRESS DEL DONGLE ed il MAC ADDRES A CUI E' STATA INVIATA LA RICHIESTA (PUO' ESSERE BROADCAST))
  CONN_REQ_TO_DONGLE,   // RICHIESTA INVIATA DA GUN A DONGLE, PER ACCETTARE RICHIESTA DI CONNESSIONE (INVIA MAC ADDRESS DEL GUN E MAC ADDRES A CUI E' STATA INVIATA RICHIESTA)
  CONN_ACCEPT,          // IL DONGLE HA ACCETTATO CONNESSIONE (INVIA ANCHE MAC ADDRES DEL DONGLE E DELLA GUN A CUI E' STATA ACCETTATA LA RICHIESTA)
};

enum CONNECTION_STATE {
  NONE_CONNECTION = 0,
  DEVICES_CONNECTED,  //5
  // =============================
  // === per controllo connessione
  NONE_CONNECTION_DONGLE,
  TX_CHECK_CONNECTION_LAST_DONGLE, //6  GUN -> DONGLE
  TX_CONFERM_CONNECTION_LAST_DONGLE, //7 DONGLE -> GUN
  DEVICES_CONNECTED_WITH_LAST_DONGLE, // 8
  // ========================================
  NONE_CONNECTION_PEDAL,
  TX_CHECK_CONNECTION_LAST_PEDAL,
  TX_CONFERM_CONNECTION_LAST_PEDAL, //7 DONGLE -> GUN
  DEVICES_CONNECTED_WITH_LAST_PEDAL, // 8
  //======CONNESSIONE GUN - DONGLE =========================
  TX_GUN_SEARCH_DONGLE_BROADCAST,  //9
  TX_DONGLE_TO_GUN_PRESENCE,  //10
  TX_GUN_TO_DONGLE_ACCEPT,  //11
  TX_DONGLE_TO_GUN_CONFERM,  //12
  //====== CONNESSIONE GUN - PEDAL =========================
  TX_GUN_SEARCH_PEDAL_BROADCAST,  //9
  TX_PEDAL_TO_GUN_PRESENCE,  //10
  TX_GUN_TO_PEDAL_ACCEPT,  //11
  TX_PEDAL_TO_GUN_CONFERM,  //12
  // ======= CONNESSIONE ULTIMO PEDAL


};

// Struttura vitale per l'USB Spoofing. 
// Passata via rete all'avvio, permette al Dongle di clonare dinamicamente 
// l'identità (VID/PID) della specifica pistola connessa in quel momento.
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
// ===================================================================================
// CLASSE: SerialWireless_ (INTERFACE PATTERN)
// ===================================================================================
// Architettura Bridge: Ereditando dalla classe standard `Stream`, facciamo 
// credere all'ecosistema Arduino (e al modulo Mamehooker) di stare parlando con una 
// normale porta seriale hardware. Sotto il cofano, i dati vengono manipolati 
// in modo asincrono, impacchettati via COBS e sparati sull'etere tramite ESP-NOW.

#ifndef _SERIAL_STREAM_H_
#define _SERIAL_STREAM_H_

class SerialWireless_ : public Stream       
  {
    
  public:
  
  uint8_t mac_esp_inteface[6];
  uint8_t mac_esp_another_card[6];
  volatile uint8_t stato_connessione_wireless = CONNECTION_STATE::NONE_CONNECTION; //0;

  // ===================================================================================
  // STRUTTURE DATI LOCK-FREE: GESTIONE BUFFER (TX/RX)
  // ===================================================================================
  // Per garantire operazioni in tempo reale durante gli interrupt (ISR) del WiFi, 
  // si utilizzano Ring Buffers (Buffer Circolari) basati su potenze di 2 (es. 128, 256).
  // L'uso di maschere bit a bit (es. & MASK) sostituisce i lenti operatori modulo (%), 
  // e la logica "writer/reader" evita l'uso di semafori bloccanti nelle fasi critiche.
  
  // ======= per FIFO SERIAL ===============
  // ===== per write === buffer lineare ====
  #define FIFO_SIZE_WRITE_SERIAL 128 // deve essere una potenza di 2
  #define TIME_OUT_SERIAL_WRITE 3  // in seguito rimuovere o rideterminare in microsecondi
  #define TIMER_HANDLE_SERIAL_DURATION_MICROS 3000 // 3000 microsecondi un 3 millisecondi
  uint8_t bufferSerialWrite[FIFO_SIZE_WRITE_SERIAL];
  
  esp_timer_handle_t timer_handle_serial;
  // Nuove variabili TX Serial (Circolare)
  volatile uint16_t _writerSerialWrite = 0;
  volatile uint16_t _readerSerialWrite = 0;
  const uint16_t MASK_WRITE_SERIAL = FIFO_SIZE_WRITE_SERIAL - 1;

  // ====== per read ====== buffer circolare =====
  #define FIFO_SIZE_READ_SERIAL 256 // deve essere una potenza di due
  uint8_t bufferSerialRead[FIFO_SIZE_READ_SERIAL];
  volatile uint16_t _writerSerialRead = 0;
  volatile uint16_t _readerSerialRead = 0;
  bool _overflow_bufferSerialRead = false; 
  void write_on_rx_serialBuffer(const uint8_t *data, int len);

  // Maschera per RX Serial
  const uint16_t MASK_READ_SERIAL = FIFO_SIZE_READ_SERIAL - 1;
  // ============================================
  volatile bool _serial_needs_recovery = false;

  Packet  packet;
  
  // per buffer lettura
  //volatile uint16_t _readLen = 0;
  volatile uint16_t _writer = 0;
  volatile uint16_t _reader = 0;
  #define FIFO_SIZE_READ 1024 // 1024, 2048, 4096 deve essere una potenza di 2 buffer lettura e si vuole ottimizzazione
  uint8_t _queue[FIFO_SIZE_READ];
  bool _overflow_read = false; 
  // fine per buffer lettura

  // per buffer scrittura 
  volatile uint16_t writeIndex = 0;
  volatile uint16_t readIndex = 0;
  //volatile uint16_t _writeLen = 0;
  #define BUFFER_SIZE 1024 //buffer scrittura  -- 1024, 2048, 4096 deve essere una potenza di 2
  uint8_t buffer[BUFFER_SIZE];
  bool _overflow_write = false; 
  // fine per buffer scrittura
 
  // ===================================================================================
  // CONTRATTO API: FUNZIONI DELLA CLASSE STREAM E FUNZIONI CUSTOM
  // ===================================================================================
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
  int availableForWriteBin();
  bool flush_sem();
  // inserire da me per gestione buffer ingresso
  int peekBin();
  int readBin(); // mai usata
  int availableBin();
  int availableBufferSerialWrite();

  
        
  // ======== generiche ============
  void SendData_sem();
  void SendPacket(const uint8_t *data, uint8_t len,uint8_t packetID); // non penso lo utilizzeremo

  bool checkForRxPacket(); // andrà nel ciclo principale .. messo nella callback

  bool connection_dongle();
  bool connection_gun();
  bool connection_pedal();
  bool connection_gun_at_last_dongle();
  bool connection_gun_at_pedal();
  bool connection_gun_at_last_pedal();
  void tx_gun_at_dongle_pedal_ready();

  // ===============================
  // ===== per i timer ================

void setupTimerSerial();


volatile bool is_pedal_wireless_comunication; // serve per sincronizzazione col dongle, ma andrebbe rivisto e poi probabilmente tolto

private:

};
extern SerialWireless_ SerialWireless;

#endif // _SERIAL_STREAM_H_

#endif //OPENFIRE_WIRELESS_ENABLE