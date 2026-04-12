#ifdef OPENFIRE_WIRELESS_ENABLE

#include "OpenFIRE_Packet.h"

//PacketCRC crc;  0x9B = 10011011 = 155
//#define POLY_CRC 0x9B
const uint8_t poly = 0x9B;
void Packet::generateTable_crc()
{
	for (uint16_t i = 0; i < 256; ++i) {
		int curr = i;
		for (int j = 0; j < 8; ++j) {
			if ((curr & 0x80) != 0) curr = (curr << 1) ^ (int)poly;
			else curr <<= 1;
		}
		csTable_crc[i] = (byte)curr;
	}
}
	
uint8_t Packet::calculate_crc(const uint8_t& val)
{
	if (val < 256) return csTable_crc[val];
	return 0;
}
	
uint8_t Packet::calculate_crc(uint8_t arr[], uint8_t len)
{
	uint8_t crc = 0;
	for (uint16_t i = 0; i < len; i++) crc = csTable_crc[crc ^ arr[i]];
	return crc;
}

void Packet::begin(const configST& configs)
{
	port         = configs.port; 
	debugPort    = configs.debugPort; 
	debug        = configs.debug;
	callbacks    = configs.callbacks;
	callbacksLen = configs.callbacksLen;
	timeout 	 = configs.timeout;
	
	txBuff[0] = START_BYTE;
	generateTable_crc();
}

void Packet::begin(Stream& _port, const bool& _debug, Stream& _debugPort, const uint32_t& _timeout)
{
	port      = &_port;
	debugPort = &_debugPort;
	debug     = _debug;
	timeout   = _timeout;

	txBuff[0] = START_BYTE;
	generateTable_crc();
}

uint8_t Packet::constructPacket(const uint8_t& messageLen, const uint8_t& packetID)
{
	txBuff[2] = packetID;
	txBuff[3] = messageLen;

	uint8_t crcVal = calculate_crc(&txBuff[2], messageLen+2);
	txBuff[messageLen + PREAMBLE_SIZE] = crcVal;
	
	calcOverhead(&txBuff[PREAMBLE_SIZE], (uint8_t)messageLen+1);
	stuffPacket(&txBuff[PREAMBLE_SIZE], (uint8_t)messageLen+1);

	txBuff[1] = overheadByte;
	txBuff[messageLen + PREAMBLE_SIZE + 1] = STOP_BYTE;
	return (uint8_t)messageLen;
}

uint8_t Packet::parse(const uint8_t& recChar, const bool& valid)
{
	bool packet_fresh = (packetStart == 0) || ((millis() - packetStart) < timeout);
	if (!packet_fresh) {
		//if (debug) debugPort->println("ERROR: STALE PACKET");
		reset();
		//bytesRead   = 0;
		//state       = find_start_byte;
		//status      = STALE_PACKET_ERROR;
		//packetStart = 0;
		return 0; //bytesRead;
	}

	if (valid) {
		switch (state) {
			case find_start_byte: /////////////////////////////////////////
			{
				if (recChar == START_BYTE) {
					state       = find_overhead_byte; //find_id_byte;
					packetStart = millis();	//start the timer
					rxBuff[0] = recChar;
				}
				break;
			}

			case find_overhead_byte: //////////////////////////////////////
			{
				recOverheadByte = recChar;
				state           = find_id_byte; //find_payload_len;
				rxBuff[1] = recChar;
				break;
			}

			case find_id_byte: ////////////////////////////////////////////
			{
				idByte = recChar;
				state  = find_payload_len; // find_overhead_byte;
				rxBuff[2] = recChar;
				break;
			}

			case find_payload_len: ////////////////////////////////////////
			{
				if ((recChar > 0) && (recChar <= MAX_PACKET_SIZE)) {
					bytesToRec = recChar;
					payIndex   = 0;
					state      = find_payload;
					rxBuff[3] = recChar;
				}
				else {
					//bytesRead = 0;
					//state     = find_start_byte;
					//status    = PAYLOAD_ERROR;
					//if (debug) debugPort->println("ERROR: PAYLOAD_ERROR");
					reset();
					return 0; //bytesRead;
				}
				break;
			}

			case find_payload: ////////////////////////////////////////////
			{
				if (payIndex < bytesToRec) {
					rxBuff[payIndex + PREAMBLE_SIZE] = recChar;
					payIndex++;
					if (payIndex == bytesToRec) state    = find_crc;
				}
				break;
			}

			case find_crc: ///////////////////////////////////////////
			{		
				rxBuff[bytesToRec + PREAMBLE_SIZE] = recChar;						
				unpackPacket(&rxBuff[PREAMBLE_SIZE], bytesToRec + 1);
				uint8_t calcCrc = calculate_crc(&rxBuff[2], bytesToRec+2);
				if (calcCrc == rxBuff[bytesToRec + PREAMBLE_SIZE]) state = find_end_byte;
				else {
					//bytesRead = 0;
					//state     = find_start_byte;
					//status    = CRC_ERROR;
					//if (debug) debugPort->println("ERROR: CRC_ERROR");
					reset();
					return 0; //bytesRead;
				}
				break;
			}

			case find_end_byte: ///////////////////////////////////////////
			{
				state = find_start_byte;
				if (recChar == STOP_BYTE) {
					bytesRead = bytesToRec;
					//status    = NEW_DATA;
					if (callbacks) callbacks[0]();
					packetStart = 0;	// reset the timer
					return bytesToRec;
				}
				//bytesRead = 0;
				//status    = STOP_BYTE_ERROR;
				//if (debug) debugPort->println("ERROR: STOP_BYTE_ERROR");
				reset();
				return 0; //bytesRead;
				break;
			}

			default:
			{
				/*
				if (debug) {
					debugPort->print("ERROR: Undefined state ");
					debugPort->println(state);
				}
				*/
				reset();
				//bytesRead = 0;
				//state     = find_start_byte;
				break;
			}
		}
	}
	else {
		bytesRead = 0;
		//status    = NO_DATA;
		return 0; //bytesRead;
	}
	bytesRead = 0;
	//status    = CONTINUE;
	return 0; //bytesRead;
}

uint8_t Packet::currentPacketID()
{
	return idByte;
}

void Packet::calcOverhead(uint8_t arr[], const uint8_t& len)
{
	overheadByte = MAX_COBS;  // 254 ===== //0xFF;

	for (uint8_t i = 0; i < len; i++) {
		if (arr[i] == START_BYTE) {
			overheadByte = i;
			break;
		}
	}
}

int16_t Packet::findLast(uint8_t arr[], const uint8_t& len)
{
	for (uint8_t i = (len - 1); i != 0xFF; i--) if (arr[i] == START_BYTE) return i;
	return -1;
}

void Packet::stuffPacket(uint8_t arr[], const uint8_t& len)
{
	int16_t refByte = findLast(arr, len);

	if (refByte != -1) {
		for (uint8_t i = (len - 1); i != 0xFF; i--) {
			if (arr[i] == START_BYTE) {
				arr[i]  = refByte - i;
				refByte = i;
			}
		}
	}
}

void Packet::unpackPacket(uint8_t arr[], const uint8_t &len)
{
	uint8_t testIndex = recOverheadByte;
	uint8_t delta     = 0;

	if (testIndex != MAX_COBS) {
		while (arr[testIndex]) {
			delta          = arr[testIndex];
			arr[testIndex] = START_BYTE;
			testIndex += delta;
		}
		arr[testIndex] = START_BYTE;
	}
}

void Packet::reset()
{
	bytesRead   = 0;
	packetStart = 0;
	state       = find_start_byte;
}

#endif //OPENFIRE_WIRELESS_ENABLE