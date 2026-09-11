/*  OpenFIRE Web App - simulated lightgun (development without hardware).

    Serves a folder over HTTP and a WebSocket at /ws connected to the mock firmware,
    like the web server of the ESP32 in web configuration mode.

      node tests/sim-server.js --root ../dist/device --board waveshare-esp32-s3-zero
      node tests/sim-server.js --root .            (unbundled webapp: open /?ws)

    Options: --port 8080  --board <board>  --root <folder>  --latency <ms>

    Build the device folder with:
      python scripts/webapp_build.py device --board waveshare-esp32-s3-zero --out dist/device
*/
'use strict';

const crypto = require('crypto');
const fs = require('fs');
const http = require('http');
const path = require('path');
const vm = require('vm');

const { MockFirmware, MemoryLink } = require('./mock-firmware.js');

const WEBAPP = path.join(__dirname, '..');

function option(name, fallback) {
    const index = process.argv.indexOf(`--${name}`);
    return index > 0 && index + 1 < process.argv.length ? process.argv[index + 1] : fallback;
}

const MIME = {
    '.html': 'text/html; charset=utf-8', '.js': 'application/javascript; charset=utf-8',
    '.css': 'text/css; charset=utf-8', '.svg': 'image/svg+xml', '.json': 'application/json',
};

// ----- Minimal WebSocket server (RFC 6455, binary messages only) -------------

function acceptWebSocket(request, socket, onOpen) {
    const key = request.headers['sec-websocket-key'];
    const accept = crypto.createHash('sha1').update(key + '258EAFA5-E914-47DA-95CA-C5AB0DC85B11').digest('base64');
    socket.write('HTTP/1.1 101 Switching Protocols\r\nUpgrade: websocket\r\nConnection: Upgrade\r\n' +
                 `Sec-WebSocket-Accept: ${accept}\r\n\r\n`);
    socket.setNoDelay(true);

    let buffer = Buffer.alloc(0);
    let closed = false;
    const ws = {
        onMessage: null,
        onClose: null,
        send(bytes) {
            if (closed) return;
            const data = Buffer.from(bytes);
            const header = data.length < 126 ? Buffer.from([0x82, data.length])
                : Buffer.from([0x82, 126, data.length >> 8, data.length & 0xFF]);
            socket.write(Buffer.concat([header, data]));
        },
        close() {
            if (closed) return;
            closed = true;
            try { socket.end(Buffer.from([0x88, 0])); } catch (e) { /* already gone */ }
            socket.destroy();
            if (ws.onClose) ws.onClose();
        },
    };

    socket.on('data', (chunk) => {
        buffer = Buffer.concat([buffer, chunk]);
        while (buffer.length >= 2) {
            const opcode = buffer[0] & 0x0F;
            let length = buffer[1] & 0x7F;
            let offset = 2;
            if (length === 126) { if (buffer.length < 4) return; length = buffer.readUInt16BE(2); offset = 4; }
            else if (length === 127) { if (buffer.length < 10) return; length = Number(buffer.readBigUInt64BE(2)); offset = 10; }
            const masked = (buffer[1] & 0x80) !== 0;
            const total = offset + (masked ? 4 : 0) + length;
            if (buffer.length < total) return;
            const mask = masked ? buffer.subarray(offset, offset + 4) : null;
            const payload = Buffer.from(buffer.subarray(total - length, total));
            if (mask) for (let i = 0; i < payload.length; ++i) payload[i] ^= mask[i & 3];
            buffer = buffer.subarray(total);

            if (opcode === 0x8) { ws.close(); return; }
            if (opcode === 0x9) { socket.write(Buffer.concat([Buffer.from([0x8A, payload.length]), payload])); continue; }
            if ((opcode === 0x2 || opcode === 0x1 || opcode === 0x0) && ws.onMessage) ws.onMessage(new Uint8Array(payload));
        }
    });
    socket.on('close', () => { if (!closed) { closed = true; if (ws.onClose) ws.onClose(); } });
    socket.on('error', () => {});
    onOpen(ws);
}

// ----- Simulated lightgun ----------------------------------------------------

function startServer({ port = 8080, root = WEBAPP, board = 'waveshare-esp32-s3-zero', latency = 2 } = {}) {
    if (!globalThis.OpenFIREshared) {
        vm.runInThisContext(fs.readFileSync(path.join(WEBAPP, 'boards', 'OpenFIREshared.js'), 'utf8') +
                            '\nglobalThis.OpenFIREshared = OpenFIREshared;');
    }
    // The page talks to the WebSocket link; serialLink is the (unused) USB port.
    const serialLink = new MemoryLink({ latency });
    const link = new MemoryLink({ latency });
    const firmware = new MockFirmware(globalThis.OpenFIREshared, serialLink, { boardType: board, webLink: link });
    firmware.loadPresets();
    firmware.start();

    let client = null;
    link.toApp = (bytes) => { if (client) client.send(bytes); };

    const server = http.createServer((request, response) => {
        const url = new URL(request.url, 'http://localhost');
        let file = path.normalize(path.join(root, decodeURIComponent(url.pathname)));
        if (!file.startsWith(path.normalize(root))) { response.writeHead(403); response.end(); return; }
        if (fs.existsSync(file) && fs.statSync(file).isDirectory()) file = path.join(file, 'index.html');
        if (!fs.existsSync(file)) { response.writeHead(404); response.end('not found'); return; }
        response.writeHead(200, { 'Content-Type': MIME[path.extname(file)] || 'application/octet-stream', 'Cache-Control': 'no-cache' });
        fs.createReadStream(file).pipe(response);
    });

    server.on('upgrade', (request, socket) => {
        if (new URL(request.url, 'http://localhost').pathname !== '/ws') { socket.destroy(); return; }
        acceptWebSocket(request, socket, (ws) => {
            // Like OpenFIREweb.cpp: a new page replaces the previous client, which is closed.
            if (client) { const old = client; client = null; old.close(); firmware.clientLost = true; }
            client = ws;
            ws.onMessage = (bytes) => link.appWrite(bytes);
            ws.onClose = () => {
                if (client === ws) { client = null; firmware.clientLost = true; }
            };
        });
    });

    return new Promise((resolve) => server.listen(port, () => resolve({
        server, firmware, link, serialLink,
        dropClient() { if (client) client.close(); },
        async close() {
            if (client) client.close();
            await firmware.stop();
            await new Promise((done) => server.close(done));
        },
    })));
}

module.exports = { startServer };

if (require.main === module) {
    const port = Number(option('port', 8080));
    const root = path.resolve(option('root', WEBAPP));
    const board = option('board', 'waveshare-esp32-s3-zero');
    startServer({ port, root, board, latency: Number(option('latency', 2)) }).then(() => {
        console.log(`Simulated ${board} on http://localhost:${port}/ (files from ${root}, WebSocket /ws)`);
    });
}
