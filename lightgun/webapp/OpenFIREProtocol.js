// --- Costanti del protocollo OpenFIRE (Copiate dal firmware/App Qt) ---
const OF_CONST = {
    sDock1: 0xC6,
    sDock2: 0x6C,
    APP_SERIAL_START_1: 0xA5,
    APP_SERIAL_START_2: 0x5A,
    APP_SERIAL_TYPE_REQUEST: 0x10,
    APP_SERIAL_TYPE_RESPONSE: 0x20,
    APP_SERIAL_TYPE_ACK: 0x40,
    APP_SERIAL_TYPE_EVENT: 0x80,
    APP_SERIAL_FLAG_FINAL: 0x01,
    APP_SERIAL_OVERHEAD: 7
};

// --- Funzione CRC8 (Porting speculare dal C++ di appserial.cpp) ---
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

function BuildFrame(typeFlags, command, sequence, payload = null) {
    const payloadLen = payload ? payload.length : 0;
    const frame = new Uint8Array(OF_CONST.APP_SERIAL_OVERHEAD + payloadLen);
    frame[0] = OF_CONST.APP_SERIAL_START_1;
    frame[1] = OF_CONST.APP_SERIAL_START_2;
    frame[2] = typeFlags;
    frame[3] = command;
    frame[4] = sequence;
    frame[5] = payloadLen;
    if (payloadLen > 0) frame.set(payload, 6);
    frame[6 + payloadLen] = AppSerialCRC8(frame.subarray(2, 6 + payloadLen));
    return frame;
}

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
        this.isOpen = false;
        this.gunConfig = null;
    }

    async connect() {
        const hostname = window.location.hostname;
        const isLocalHost = hostname === "openfire.local" || /^192\.168\./.test(hostname) || hostname === "10.0.0.1" || hostname === "localhost";
        
        if (isLocalHost) {
            try {
                await this.connectWebSocket(ws:// + hostname + /ws);
                this.isOpen = true;
                return true;
            } catch (e) {
                console.warn("WebSocket fallito, provo fallback USB...");
            }
        }
        
        const success = await this.connectWebSerial();
        if (success) this.isOpen = true;
        return success;
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
            alert("Web Serial API non supportata da questo browser. Usa Chrome/Edge o la modalit Wi-Fi.");
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
                console.warn("CRC Errato! Scarto il pacchetto.");
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
    // Operazioni ad alto livello
    // ========================================================================

    async sendCommand(command, payload = null) {
        if (this.appSerialTxSequence === undefined) this.appSerialTxSequence = 0;
        this.appSerialTxSequence = (this.appSerialTxSequence + 1) % 256;
        if (this.appSerialTxSequence === 0) this.appSerialTxSequence = 1;
        
        const frame = BuildFrame(OF_CONST.APP_SERIAL_TYPE_REQUEST, command, this.appSerialTxSequence, payload);
        await this.write(frame);
        return true;
    }

    async waitForResponse(command, timeoutMs = 3000) {
        for (let i = 0; i < timeoutMs / 10; i++) {
            const idx = this.appSerialResponses.findIndex(f => f.command === command);
            if (idx !== -1) {
                return this.appSerialResponses.splice(idx, 1)[0];
            }
            await new Promise(r => setTimeout(r, 10));
        }
        return null;
    }

    async beginDock() {
        this.appSerialResponses = [];
        
        const dockHandshake = new Uint8Array([OF_CONST.sDock1, OF_CONST.sDock2]);
        await this.write(dockHandshake);
        
        const deviceInfoFrame = await this.waitForResponse(OpenFIREshared.serialCmdTypes_e.sDock2, 3000);
        if (!deviceInfoFrame) throw new Error("Timeout durante l'handshake Docked!");

        const boardInfoRaw = deviceInfoFrame.payload;
        const firstSeparator = boardInfoRaw.indexOf(OpenFIREshared.serialCmdTypes_e.serialTerminator);
        const secondSeparator = boardInfoRaw.indexOf(OpenFIREshared.serialCmdTypes_e.serialTerminator, firstSeparator + 1);

        const versionBytes = boardInfoRaw.subarray(0, firstSeparator);
        const typeBytes = boardInfoRaw.subarray(firstSeparator + 1, secondSeparator);
        
        const usbOffset = secondSeparator + 1;
        if (boardInfoRaw.length >= usbOffset + 18) {
            if (!this.gunConfig) this.gunConfig = { profiles: [] };
            this.gunConfig.tinyUSBtable = new Uint8Array(boardInfoRaw.subarray(usbOffset, usbOffset + 18));
        }

        const tailOffset = usbOffset + 18;
        if (boardInfoRaw.length > tailOffset && boardInfoRaw[tailOffset] === OpenFIREshared.serialCmdTypes_e.serialTerminator) {
            if (!this.gunConfig) this.gunConfig = { profiles: [] };
            this.gunConfig.currentProfile = boardInfoRaw[tailOffset + 1];
        }

        return {
            firmwareVersion: new TextDecoder().decode(versionBytes),
            boardName: new TextDecoder().decode(typeBytes).trim()
        };
    }

    async receiveSettingsRecords(command) {
        const records = [];
        const timeoutMs = 3000;
        
        while (true) {
            const foundFrame = await this.waitForResponse(command, timeoutMs);
            if (!foundFrame) throw new Error(Timeout attesa dati per comando 0x);

            if (foundFrame.typeFlags & OF_CONST.APP_SERIAL_FLAG_FINAL) {
                break;
            }

            const payload = foundFrame.payload;
            if (payload.length === 0) continue; 

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
        if (!this.gunConfig) this.gunConfig = { profiles: [] };
        const config = this.gunConfig;
        if (!config.toggles) config.toggles = {};
        if (!config.pins) config.pins = {};
        if (!config.settings) config.settings = {};
        if (!config.buttons) config.buttons = {};
        if (!config.profiles) config.profiles = [];

        const cmds = OpenFIREshared.serialCmdTypes_e;

        await this.sendCommand(cmds.sGetToggles);
        const toggles = await this.receiveSettingsRecords(cmds.sGetToggles);
        toggles.forEach(r => config.toggles[r.fieldName] = (r.value !== 0));

        if (config.toggles["CustomPins"]) {
            await this.sendCommand(cmds.sGetPins);
            const pins = await this.receiveSettingsRecords(cmds.sGetPins);
            pins.forEach(r => config.pins[r.fieldName] = r.value);
        }

        await this.sendCommand(cmds.sGetSettings);
        const settings = await this.receiveSettingsRecords(cmds.sGetSettings);
        settings.forEach(r => config.settings[r.fieldName] = r.value);

        await this.sendCommand(cmds.sGetBtns);
        const btns = await this.receiveSettingsRecords(cmds.sGetBtns);
        btns.forEach(r => config.buttons[r.fieldName] = r.value);

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

        return config;
    }

    async commitSettings() {
        if (!this.isOpen) return false;

        const cmds = OpenFIREshared.serialCmdTypes_e;

        if (!await this.sendCommand(cmds.sCommitStart)) return false;
        const startAck = await this.waitForResponse(cmds.sCommitStart, 2000);
        if (!startAck || startAck.payload.length > 0 || (startAck.typeFlags & OF_CONST.APP_SERIAL_FLAG_FINAL)) {
            console.error("Failed to start commit.");
            return false;
        }

        const buildPayload = (name, valueSize, profNum, valueDataView) => {
            const nameBytes = new TextEncoder().encode(name);
            const len = nameBytes.length + 1 + 1 + (profNum !== null ? 1 : 0) + valueSize;
            const payload = new Uint8Array(len);
            let offset = 0;
            payload.set(nameBytes, offset); offset += nameBytes.length;
            payload[offset++] = 0; // null terminator
            payload[offset++] = valueSize;
            if (profNum !== null) payload[offset++] = profNum;
            
            payload.set(new Uint8Array(valueDataView.buffer, valueDataView.byteOffset, valueSize), offset);
            return payload;
        };

        const sendRec = async (cmd, name, val, size, profNum = null) => {
            const dv = new DataView(new ArrayBuffer(size));
            if (size === 1) dv.setInt8(0, val);
            else if (size === 4 && typeof val === 'number') dv.setUint32(0, val, true);
            else if (size === 16 && typeof val === 'string') {
                const strBytes = new TextEncoder().encode(val.substring(0, 15));
                for(let i=0; i<strBytes.length; i++) dv.setUint8(i, strBytes[i]);
                dv.setUint8(strBytes.length, 0);
            }
            const payload = buildPayload(name, size, profNum, dv);
            return await this.sendCommand(cmd, payload);
        };

        for (const [key, val] of Object.entries(this.gunConfig.toggles)) {
            await sendRec(cmds.sCommitToggles, key, val ? 1 : 0, 1);
        }

        if (this.gunConfig.toggles["CustomPins"]) {
            for (const [key, val] of Object.entries(this.gunConfig.pins)) {
                if (OpenFIREshared.boardInputs_Strings[key] === undefined) continue;
                await sendRec(cmds.sCommitPins, key, val, 1);
            }
        }

        for (const [key, val] of Object.entries(this.gunConfig.settings)) {
            await sendRec(cmds.sCommitSettings, key, val, 4);
        }

        for (const [key, arr] of Object.entries(this.gunConfig.buttons)) {
            if (OpenFIREshared.boardInputs_Strings[key] === undefined) continue;
            const dv = new DataView(new ArrayBuffer(6));
            for(let i=0; i<6; i++) dv.setUint8(i, arr[i] || 0);
            const payload = buildPayload(key, 6, null, dv);
            await this.sendCommand(cmds.sCommitBtns, payload);
        }

        let currentProfSent = false;
        for (let i = 0; i < 4; i++) {
            if (!this.gunConfig.profiles[i]) continue;
            
            if (!currentProfSent) {
                await sendRec(cmds.sCommitProfile, "CurrentProf", this.gunConfig.currentProfile, 1, null);
                currentProfSent = true;
            }

            for (const [key, val] of Object.entries(this.gunConfig.profiles[i])) {
                const ignoreKeys = ["TopOffset", "BtmOffset", "LftOffset", "RhtOffset", "TLLed", "TRLed", "AdjX", "AdjY"];
                if (ignoreKeys.includes(key)) continue;

                const size = (key === "Name") ? 16 : 4;
                await sendRec(cmds.sCommitProfile, key, val, size, i);
            }
        }

        if (this.gunConfig.tinyUSBtable) {
            await this.sendCommand(cmds.sCommitID, this.gunConfig.tinyUSBtable);
            const r2 = await this.waitForResponse(cmds.sCommitID, 2000);
            if (!r2) return false;
        }

        if (!await this.sendCommand(cmds.sSave)) return false;
        const saveAck = await this.waitForResponse(cmds.sSave, 5000);
        if (!saveAck || saveAck.payload.length === 0 || saveAck.payload[0] !== 1) {
            console.error("Save failed by device.");
            return false;
        }

        console.log("Save successful!");
        return true;
    }
}
