// ============================================================================
// Costanti del protocollo OpenFIRE (Da appserial.h e OpenFIREshared.h)
// ============================================================================
const OF_CONST = {
    // Framing
    APP_SERIAL_START_1: 0xA5,
    APP_SERIAL_START_2: 0x5A,
    APP_SERIAL_MAX_PAYLOAD: 200,
    APP_SERIAL_OVERHEAD: 7,

    // Types
    APP_SERIAL_TYPE_REQUEST:  0x00,
    APP_SERIAL_TYPE_RESPONSE: 0x01,
    APP_SERIAL_TYPE_EVENT:    0x02,
    APP_SERIAL_TYPE_ACK:      0x03,
    APP_SERIAL_TYPE_MASK:     0x03,
    APP_SERIAL_FLAG_FINAL:    0x80,

    // Commands
    sDock1: 1,
    sDock2: 2,
    serialTerminator: 0 
};

// ============================================================================
// Utilities: CRC8 e Assemblatore Frame (Porting di appserial.cpp)
// ============================================================================

function AppSerialCRC8(dataArray) {
    let crc = 0;
    for (let i = 0; i < dataArray.length; i++) {
        crc ^= dataArray[i];
        for (let bit = 0; bit < 8; bit++) {
            crc = (crc & 0x80) ? ((crc << 1) ^ 0x9B) : (crc << 1);
        }
    }
    return crc & 0xFF;
}

function BuildFrame(typeFlags, command, sequence, payloadUint8) {
    const length = payloadUint8 ? payloadUint8.length : 0;
    const frameLength = OF_CONST.APP_SERIAL_OVERHEAD + length;
    const buffer = new Uint8Array(frameLength);
    
    buffer[0] = OF_CONST.APP_SERIAL_START_1;
    buffer[1] = OF_CONST.APP_SERIAL_START_2;
    buffer[2] = typeFlags;
    buffer[3] = command;
    buffer[4] = sequence;
    buffer[5] = length;
    
    if (length > 0) {
        buffer.set(payloadUint8, 6);
    }
    
    const crcData = buffer.subarray(2, 6 + length);
    buffer[6 + length] = AppSerialCRC8(crcData);
    
    return buffer;
}

// ============================================================================
// Gestore Connessione Ibrido (WebSocket / WebSerial) e Parser
// ============================================================================

class OpenFIREConnection {
    constructor() {
        this.socket = null;
        this.serialPort = null;
        this.serialReader = null;
        this.serialWriter = null;
        this.isWebSocket = false;
        
        // Stato del parser e protocollo
        this.rxBuffer = new Uint8Array(0);
        this.appSerialResponses = [];
        this.onEventReceived = null; // Callback per gli EVENT
    }

    async connect() {
        const hostname = window.location.hostname;
        // Auto-Discovery per WebSocket
        const isLocalHost = hostname === "openfire.local" || /^192\.168\./.test(hostname) || hostname === "10.0.0.1" || hostname === "localhost";
        
        if (isLocalHost) {
            try {
                await this.connectWebSocket(`ws://${hostname}/ws`);
                return true;
            } catch (e) {
                console.warn("WebSocket fallito, provo fallback USB...");
            }
        }
        
        return await this.connectWebSerial();
    }

    connectWebSocket(url) {
        return new Promise((resolve, reject) => {
            this.socket = new WebSocket(url);
            this.socket.binaryType = "arraybuffer"; 

            this.socket.onopen = () => {
                this.isWebSocket = true;
                console.log("Connesso via Wi-Fi (WebSocket)!");
                resolve();
            };

            this.socket.onerror = (err) => reject(err);

            this.socket.onmessage = (event) => {
                this.processIncoming(new Uint8Array(event.data));
            };
        });
    }

    async connectWebSerial() {
        if (!navigator.serial) {
            alert("Web Serial API non supportata da questo browser. Usa Chrome/Edge o la modalità Wi-Fi.");
            return false;
        }

        try {
            this.serialPort = await navigator.serial.requestPort();
            await this.serialPort.open({ baudRate: 9600 });
            this.isWebSocket = false;
            console.log("Connesso via Cavo USB (Web Serial)!");

            this.serialWriter = this.serialPort.writable.getWriter();
            this.serialReader = this.serialPort.readable.getReader();
            
            this.readSerialLoop();
            return true;
        } catch (e) {
            console.error("Connessione Seriale interrotta o rifiutata:", e);
            return false;
        }
    }

    async readSerialLoop() {
        try {
            while (true) {
                const { value, done } = await this.serialReader.read();
                if (done) break;
                if (value) {
                    this.processIncoming(new Uint8Array(value));
                }
            }
        } catch (error) {
            console.error("Errore di lettura Seriale:", error);
        } finally {
            this.serialReader.releaseLock();
        }
    }

    async write(uint8Array) {
        if (this.isWebSocket && this.socket) {
            this.socket.send(uint8Array);
        } else if (this.serialWriter) {
            await this.serialWriter.write(uint8Array);
        } else {
            console.error("Nessuna connessione attiva per scrivere!");
        }
    }

    // ========================================================================
    // Protocollo Binario
    // ========================================================================

    processIncoming(newBytes) {
        const combined = new Uint8Array(this.rxBuffer.length + newBytes.length);
        combined.set(this.rxBuffer);
        combined.set(newBytes, this.rxBuffer.length);
        this.rxBuffer = combined;

        while (this.rxBuffer.length >= 2) {
            let start = 0;
            while (start + 1 < this.rxBuffer.length && 
                  (this.rxBuffer[start] !== OF_CONST.APP_SERIAL_START_1 || 
                   this.rxBuffer[start + 1] !== OF_CONST.APP_SERIAL_START_2)) {
                start++;
            }

            if (start + 1 >= this.rxBuffer.length) {
                if (this.rxBuffer[this.rxBuffer.length - 1] === OF_CONST.APP_SERIAL_START_1) {
                    this.rxBuffer = new Uint8Array([OF_CONST.APP_SERIAL_START_1]);
                } else {
                    this.rxBuffer = new Uint8Array(0);
                }
                break;
            }

            if (start > 0) {
                this.rxBuffer = this.rxBuffer.subarray(start);
            }

            if (this.rxBuffer.length < OF_CONST.APP_SERIAL_OVERHEAD) {
                break; 
            }

            const length = this.rxBuffer[5];
            const totalFrameSize = OF_CONST.APP_SERIAL_OVERHEAD + length;

            if (this.rxBuffer.length < totalFrameSize) {
                break;
            }

            const crcData = this.rxBuffer.subarray(2, 6 + length);
            const expectedCRC = AppSerialCRC8(crcData);
            const actualCRC = this.rxBuffer[6 + length];

            if (expectedCRC === actualCRC) {
                const frame = {
                    typeFlags: this.rxBuffer[2],
                    command: this.rxBuffer[3],
                    sequence: this.rxBuffer[4],
                    payload: this.rxBuffer.subarray(6, 6 + length)
                };

                this.handleFrame(frame);
                this.rxBuffer = this.rxBuffer.subarray(totalFrameSize);
            } else {
                console.warn("CRC Errato! Scarto il pacchetto corroto.");
                this.rxBuffer = this.rxBuffer.subarray(2);
            }
        }
    }

    handleFrame(frame) {
        const type = frame.typeFlags & OF_CONST.APP_SERIAL_TYPE_MASK;
        if (type === OF_CONST.APP_SERIAL_TYPE_EVENT) {
            if (this.onEventReceived) this.onEventReceived(frame);
        } else {
            this.appSerialResponses.push(frame);
        }
    }

    // ========================================================================
    // Operazioni ad alto livello (Handshake)
    // ========================================================================

    async beginDock() {
        this.appSerialResponses = []; // Svuota lo stato
        
        // AppSerial::BeginDock() in C++ invia nudi i byte sDock1 e sDock2
        const dockHandshake = new Uint8Array([OF_CONST.sDock1, OF_CONST.sDock2]);
        await this.write(dockHandshake);
        console.log("Handshake sDock inviato, in attesa di boardInfo...");

        // Ora aspettiamo il pacchetto RESPONSE con command == sDock2
        // Simuliamo un timeout asincrono (polling per semplicità)
        for(let i = 0; i < 30; i++) { // Timeout ~3 secondi
            await new Promise(r => setTimeout(r, 100));
            const responseIndex = this.appSerialResponses.findIndex(f => f.command === OF_CONST.sDock2);
            if (responseIndex !== -1) {
                const deviceInfoFrame = this.appSerialResponses.splice(responseIndex, 1)[0];
                console.log("Docked con successo! BoardInfo ricevuto:", new TextDecoder().decode(deviceInfoFrame.payload));
                return deviceInfoFrame.payload;
            }
        }
        throw new Error("Timeout durante l'handshake Docked!");
    }

    async sendCommand(command, payload = null) {
        if (this.appSerialTxSequence === undefined) this.appSerialTxSequence = 0;
        this.appSerialTxSequence = (this.appSerialTxSequence + 1) % 256;
        if (this.appSerialTxSequence === 0) this.appSerialTxSequence = 1;
        
        const frame = BuildFrame(OF_CONST.APP_SERIAL_TYPE_REQUEST, command, this.appSerialTxSequence, payload);
        await this.write(frame);
    }

    async receiveSettingsRecords(command) {
        const records = [];
        const timeoutMs = 3000;
        
        while (true) {
            let foundFrame = null;
            for (let i = 0; i < timeoutMs / 10; i++) {
                const idx = this.appSerialResponses.findIndex(f => f.command === command);
                if (idx !== -1) {
                    foundFrame = this.appSerialResponses.splice(idx, 1)[0];
                    break;
                }
                await new Promise(r => setTimeout(r, 10));
            }

            if (!foundFrame) throw new Error(`Timeout attesa dati per comando 0x${command.toString(16)}`);

            if (foundFrame.typeFlags & OF_CONST.APP_SERIAL_FLAG_FINAL) {
                break;
            }

            const payload = foundFrame.payload;
            if (payload.length === 0) continue; // Salta gli ACK vuoti iniziali

            const nullPos = payload.indexOf(0);
            if (nullPos < 0) continue;

            const fieldName = new TextDecoder().decode(payload.subarray(0, nullPos));
            let pos = nullPos + 1;
            const valueSize = payload[pos++];
            
            let profNum = null;
            if (command === OpenFIREshared.serialCmdTypes_e.sGetProfile && fieldName !== "CurrentProf") {
                profNum = payload[pos++];
            }
            
            const valueBytes = payload.subarray(pos, pos + valueSize);
            const dataView = new DataView(valueBytes.buffer, valueBytes.byteOffset, valueBytes.byteLength);

            let value = 0;
            if (valueSize === 1) {
                value = dataView.getInt8(0);
            } else if (valueSize === 2) {
                value = dataView.getInt16(0, true);
            } else if (valueSize === 4) {
                if (["TLled", "TRled", "AdjX", "AdjY"].includes(fieldName)) {
                    value = dataView.getFloat32(0, true);
                } else if (["TopOffset", "BottomOffset", "LeftOffset", "RightOffset"].includes(fieldName)) {
                    value = dataView.getInt32(0, true);
                } else {
                    value = dataView.getUint32(0, true);
                }
            } else {
                value = new TextDecoder().decode(valueBytes).replace(/\0/g, '');
            }

            records.push({ fieldName, profNum, value });
        }
        return records;
    }

    async syncSettings() {
        const config = { toggles: {}, pins: {}, settings: {}, buttons: {}, profiles: [] };
        const cmds = OpenFIREshared.serialCmdTypes_e;

        console.log("Richiesta Toggles...");
        await this.sendCommand(cmds.sGetToggles);
        const toggles = await this.receiveSettingsRecords(cmds.sGetToggles);
        toggles.forEach(r => config.toggles[r.fieldName] = (r.value !== 0));

        if (config.toggles["CustomPins"]) {
            console.log("Richiesta Pins...");
            await this.sendCommand(cmds.sGetPins);
            const pins = await this.receiveSettingsRecords(cmds.sGetPins);
            pins.forEach(r => config.pins[r.fieldName] = r.value);
        }

        console.log("Richiesta Settings...");
        await this.sendCommand(cmds.sGetSettings);
        const settings = await this.receiveSettingsRecords(cmds.sGetSettings);
        settings.forEach(r => config.settings[r.fieldName] = r.value);

        console.log("Richiesta Buttons...");
        await this.sendCommand(cmds.sGetBtns);
        const btns = await this.receiveSettingsRecords(cmds.sGetBtns);
        btns.forEach(r => config.buttons[r.fieldName] = r.value);

        console.log("Richiesta Profili...");
        await this.sendCommand(cmds.sGetProfile);
        const profiles = await this.receiveSettingsRecords(cmds.sGetProfile);
        profiles.forEach(r => {
            if (r.fieldName === "CurrentProf") {
                config.currentProfile = r.value;
            } else {
                if (!config.profiles[r.profNum]) config.profiles[r.profNum] = {};
                config.profiles[r.profNum][r.fieldName] = r.value;
            }
        });

        console.log("Sincronizzazione completata con successo!", config);
        return config;
    }
}

// ============================================================================
// Logica UI Base
// ============================================================================
document.addEventListener("DOMContentLoaded", () => {
    const statusText = document.getElementById("status");
    const btnTest = document.getElementById("btn-test");
    
    const connection = new OpenFIREConnection();
    statusText.innerText = "Pronto. Clicca Connetti per avviare l'Handshake e la Sincronizzazione.";

    btnTest.addEventListener("click", async () => {
        statusText.innerText = "Connessione in corso...";
        const success = await connection.connect();
        
        if (success) {
            statusText.innerText = "Connesso! Avvio Handshake...";
            try {
                const boardInfo = await connection.beginDock();
                statusText.innerText = "Docked! Sincronizzazione dei Settings in corso...";
                
                const gunConfig = await connection.syncSettings();
                statusText.innerText = `Sincronizzazione completata! Profilo attivo: ${gunConfig.currentProfile}`;
                
                // Alert o UI update con i dati completi
                console.log("Gun Config Completo: ", gunConfig);
            } catch (err) {
                statusText.innerText = "Errore: " + err.message;
            }
        } else {
            statusText.innerText = "Errore di connessione.";
        }
    });
});
