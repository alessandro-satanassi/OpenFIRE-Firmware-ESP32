/*!
 * @file OpenFIRE_Packet.cpp
 * @brief Library for OpenFIRE Packet wireless with Esp-Now
 * @n CPP Library for OpenFIRE Packet wireless with Esp-Now
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 * 
 * This code is partially derived from the SerialTransfer library by PowerBroker2
 * (https://github.com/PowerBroker2/SerialTransfer), and I would like to express my gratitude to him.
 * 
 */
#ifdef OPENFIRE_WIRELESS_ENABLE

#include "OpenFIRE_Packet.h"

// ===================================================================================
// INIZIALIZZAZIONE E GENERAZIONE CRC 
// ===================================================================================

const uint8_t poly = 0x9B;

// La tabella CRC viene generata dinamicamente in RAM all'avvio (boot).
// Questo trade-off sacrifica una frazione di millisecondo al boot per risparmiare 
// preziosi cicli di CPU durante il calcolo CRC bit-a-bit in esecuzione (runtime).
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
	
uint8_t Packet::calculate_crc(uint8_t val)
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
	callbacks    = configs.callbacks;
	callbacksLen = configs.callbacksLen;
	timeout 	 = configs.timeout;
	
	txBuff[0] = START_BYTE;
	generateTable_crc();
}

void Packet::begin(Stream& _port, const bool& _debug, Stream& _debugPort, const uint32_t& _timeout)
{
	timeout   = _timeout;

	txBuff[0] = START_BYTE;
	generateTable_crc();
}

// ===================================================================================
// COSTRUZIONE DEL PACCHETTO (TX)
// ===================================================================================

uint8_t Packet::constructPacket(uint8_t messageLen, uint8_t packetID)
{
	txBuff[2] = packetID;
	txBuff[3] = messageLen;

	// Il CRC copre l'intestazione e il payload. Deve essere calcolato PRIMA di 
	// applicare il COBS, in modo che il ricevitore possa validare il pacchetto
	// solo dopo averlo riportato al suo stato originale.
	uint8_t crcVal = calculate_crc(&txBuff[2], messageLen+2);
	txBuff[messageLen + PREAMBLE_SIZE] = crcVal;
	
	// Il COBS "nasconde" gli 0xFF trasformandoli in puntatori al byte successivo.
	calcOverhead(&txBuff[PREAMBLE_SIZE], (uint8_t)messageLen+1);
	stuffPacket(&txBuff[PREAMBLE_SIZE], (uint8_t)messageLen+1);

	txBuff[1] = overheadByte;
	txBuff[messageLen + PREAMBLE_SIZE + 1] = STOP_BYTE;
	return (uint8_t)messageLen;
}

// ===================================================================================
// MOTORE FSM DI RICEZIONE E PARSING (RX)
// ===================================================================================
// Macchina a stati finiti reattiva. Legge un carattere alla volta uscendo subito 
// dalla funzione. Nessun ciclo bloccante per rispettare i tempi di ESP-NOW.

uint8_t Packet::parse(uint8_t recChar, bool valid)
{
	// Valutazione a cortocircuito: se packetStart == 0 (siamo in stato di attesa), 
	// la condizione è vera e millis() NON viene valutata. Risparmia chiamate CPU.
	bool packet_fresh = (packetStart == 0) || ((millis() - packetStart) < timeout);
	if (!packet_fresh) {
		reset();
		// RIMOSSO: return 0; 
		// Il byte arrivato DEVE essere passato al blocco valid qui sotto, 
		// perché quasi certamente è lo START_BYTE di un nuovo pacchetto in ritardo!
	}

	if (valid) {
		// FIX FONDAMENTALE: Hard-Sync. 
		// Nel protocollo COBS, lo START_BYTE non compare MAI nei dati. 
		// Se ne ricevi uno, significa matematicamente che sta iniziando un nuovo pacchetto.
		if (recChar == START_BYTE && state != find_start_byte) {
			reset(); // Azzera tutto istantaneamente
			// Lasciamo che lo switch qui sotto elabori questo START_BYTE in modo pulito!
		}

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
				
				// Rimuove l'overhead COBS per riportare i dati alla forma originale 
				// necessaria per la validazione crittografica.
				// Un overhead o un collegamento COBS che punta fuori dal pacchetto (pacchetto
				// corrotto o estraneo) lo scarta subito, senza scrivere oltre i suoi limiti.
				if (!unpackPacket(&rxBuff[PREAMBLE_SIZE], bytesToRec + 1)) {
					reset();
					return 0;
				}
				
				uint8_t calcCrc = calculate_crc(&rxBuff[2], bytesToRec+2);
				if (calcCrc == rxBuff[bytesToRec + PREAMBLE_SIZE]) state = find_end_byte;
				else {
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
					
					// Esecuzione callback se presente. Implementazione a puntatore a funzione
					// per minimizzare la latenza di notifica al layer superiore.
					if (callbacks) callbacks[0]();
					packetStart = 0;	// reset the timer
					return bytesToRec;
				}
				reset();
				return 0; //bytesRead;
				break;
			}

			default:
			{
				reset();
				break;
			}
		}
	}
	else {
		bytesRead = 0;
		return 0; //bytesRead;
	}
	bytesRead = 0;
	return 0; //bytesRead;
}

uint8_t Packet::currentPacketID()
{
	return idByte;
}

// ===================================================================================
// FUNZIONI SUPPORTO AL COBS (Consistent Overhead Byte Stuffing)
// ===================================================================================
// Sostituisce gli 0xFF nel flusso dati con l'offset (la distanza) per raggiungere 
// lo 0xFF successivo, creando una sorta di "linked list" all'interno del pacchetto.

void Packet::calcOverhead(uint8_t arr[], uint8_t len)
{
	overheadByte = MAX_COBS;  // 254 ===== //0xFF;

	for (uint8_t i = 0; i < len; i++) {
		if (arr[i] == START_BYTE) {
			overheadByte = i;
			break;
		}
	}
}

int16_t Packet::findLast(uint8_t arr[], uint8_t len)
{
	for (uint8_t i = (len - 1); i != 0xFF; i--) if (arr[i] == START_BYTE) return i;
	return -1;
}

void Packet::stuffPacket(uint8_t arr[], uint8_t len)
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

// EN: Bounded: an overhead byte or a link that points past the packet (corrupted or foreign
//     frame) returns false before anything is written outside it, and parse() drops the packet.
//     A valid packet is decoded exactly as before, with one comparison per stuffed byte.
// IT: Con limiti: un overhead o un collegamento che punta oltre il pacchetto (pacchetto corrotto
//     o estraneo) restituisce false prima di scrivere fuori, e parse() scarta il pacchetto.
//     Un pacchetto valido si decodifica come prima, con un confronto in piu' per byte sostituito.
bool Packet::unpackPacket(uint8_t arr[], uint8_t len)
{
	uint16_t testIndex = recOverheadByte; // 16 bit: the sum below cannot wrap back inside / la somma non puo' tornare dentro

	if (testIndex == MAX_COBS) return true; // no START_BYTE in the data / nessuno START_BYTE nei dati

	while (testIndex < len) {
		const uint8_t delta = arr[testIndex];
		arr[testIndex] = START_BYTE;
		if (delta == 0) return true; // last link of the chain / ultimo collegamento della catena
		testIndex += delta;
	}
	return false;
}

void Packet::reset()
{
	bytesRead   = 0;
	packetStart = 0;
	state       = find_start_byte;
}

#endif //OPENFIRE_WIRELESS_ENABLE