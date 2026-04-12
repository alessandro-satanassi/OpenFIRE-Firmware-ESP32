/*
11111111 00000000 11111111 00000000 00000000 00000000 ... 00000000 10000001
|      | |      | |      | |      | |      | |      | | | |      | |______|__Stop byte
|      | |      | |      | |      | |      | |      | | | |______|___________8-bit CRC
|      | |      | |      | |      | |      | |      | |_|____________________Rest of payload
|      | |      | |      | |      | |      | |______|________________________2nd payload byte
|      | |      | |      | |      | |______|_________________________________1st payload byte
|      | |      | |      | |______|__________________________________________# len  of payload bytes (lunghezza del contenuto del paccchetto/payload sempre minore di START_BYTE)
|      | |      | |______|___________________________________________________Packet ID (0 by default) (non deve però essere mai come Start byte)
|      | |______|____________________________________________________________COBS Overhead byte (sempre minore di START_BYTE)
|______|_____________________________________________________________________Start byte (constant)
*/

#ifdef OPENFIRE_WIRELESS_ENABLE

#pragma once
#include "Arduino.h"
//#include "PacketCRC.h"


typedef void (*functionPtr)();


const int8_t CONTINUE           = 3;
const int8_t NEW_DATA           = 2;
const int8_t NO_DATA            = 1;
const int8_t CRC_ERROR          = 0;
const int8_t PAYLOAD_ERROR      = -1;
const int8_t STOP_BYTE_ERROR    = -2;
const int8_t STALE_PACKET_ERROR = -3;

const uint8_t START_BYTE = 0xFF; // 255 -> 11111111 ==== // vecchio valore 0x7E; // 126 ->  01111110
const uint8_t STOP_BYTE  = 0x81; // 129 ->  10000001

const uint8_t PREAMBLE_SIZE   = 4;    // START_BYTE + COBS + PACKET_ID + LEN
const uint8_t POSTAMBLE_SIZE  = 2;    // CRC + STOP_BYTE
const uint8_t MAX_PACKET_SIZE = 0xF4; // 244 = 250-6  =====   //0xFE; // 254 Maximum allowed payload bytes per packet
const uint8_t MAX_COBS        = 0xFE; // 254 aggiunto da me

const uint8_t DEFAULT_TIMEOUT = 50;


struct configST
{
	Stream*            port         = &Serial;
	Stream*            debugPort    = &Serial;
	bool               debug        = true;
	const functionPtr* callbacks    = NULL;
	uint8_t            callbacksLen = 0;
	uint32_t           timeout      = __UINT32_MAX__;
};


class Packet
{
  public: // <<---------------------------------------//public
	uint8_t txBuff[MAX_PACKET_SIZE+PREAMBLE_SIZE+POSTAMBLE_SIZE];
	uint8_t rxBuff[MAX_PACKET_SIZE+PREAMBLE_SIZE+POSTAMBLE_SIZE];
	uint8_t preamble[PREAMBLE_SIZE]   = {START_BYTE, 0, 0, 0};
	uint8_t postamble[POSTAMBLE_SIZE] = {0, STOP_BYTE};

	uint8_t bytesRead = 0;
	//int8_t  status    = 0;


	void    begin(const configST& configs);
	void    begin(Stream& _port, const bool& _debug = true, Stream& _debugPort = Serial, const uint32_t& _timeout = DEFAULT_TIMEOUT);
	uint8_t constructPacket(const uint8_t& messageLen, const uint8_t& packetID = 0);
	uint8_t parse(const uint8_t& recChar, const bool& valid = true);
	uint8_t currentPacketID();
	void    reset();

	template <typename T>
	uint8_t txObj(const T* val, const uint8_t& index = 0, const uint8_t& len = sizeof(T))
	{
		memcpy(&txBuff[PREAMBLE_SIZE + index], val, len);
		return len + index;
	}

	template <typename T>
	uint16_t rxObj(const T* val, uint8_t& index = 0, const uint8_t& len = sizeof(T))
	{
		memcpy(val, &txBuff[PREAMBLE_SIZE + index], len);
		return len + index;
	}


  private: // <<---------------------------------------//private
	enum fsm
	{
		find_start_byte,
		find_id_byte,
		find_overhead_byte,
		find_payload_len,
		find_payload,
		find_crc,
		find_end_byte
	};
	fsm state = find_start_byte;

	const functionPtr* callbacks    = NULL;
	uint8_t            callbacksLen = 0;

	Stream* port;  // aggiunto da mw
	Stream* debugPort;
	bool    debug = false;

	uint8_t bytesToRec      = 0;
	uint8_t payIndex        = 0;
	uint8_t idByte          = 0;
	uint8_t overheadByte    = 0;
	uint8_t recOverheadByte = 0;

	uint32_t packetStart    = 0;
	uint32_t timeout = DEFAULT_TIMEOUT;

	uint8_t csTable_crc[256]; 
	void generateTable_crc();
	uint8_t calculate_crc(const uint8_t& val);
	uint8_t calculate_crc(uint8_t arr[], uint8_t len);

	void    calcOverhead(uint8_t arr[], const uint8_t& len);
	int16_t findLast(uint8_t arr[], const uint8_t& len);
	void    stuffPacket(uint8_t arr[], const uint8_t& len);
	void    unpackPacket(uint8_t arr[], const uint8_t &len);
};

#endif //OPENFIRE_WIRELESS_ENABLE
