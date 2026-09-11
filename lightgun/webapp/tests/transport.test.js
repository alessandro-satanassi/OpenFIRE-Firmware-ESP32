/*  OpenFIRE Web App: transport tests with fake browser APIs (Web Serial, WebSocket).
    Run from lightgun/webapp:   node --test tests/*.test.js
*/
'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('fs');
const path = require('path');
const vm = require('vm');

vm.runInThisContext(fs.readFileSync(path.join(__dirname, '..', 'js/core/transport.js'), 'utf8'));
const OF = globalThis.OF;

if (!process.env.OF_TEST_VERBOSE)
    console.warn = () => {};

function setNavigator(value) {
    Object.defineProperty(globalThis, 'navigator', { value, configurable: true, writable: true });
}

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

// ----- Fake Web Serial -------------------------------------------------------

class FakeSerialPort extends EventTarget {
    constructor(info = { usbVendorId: 0xF143, usbProductId: 0x0002 }) {
        super();
        this.info = info;
        this.opened = [];
        this.signals = [];
        this.written = [];
        this.readable = null;
        this.writable = null;
        this.busy = false;
    }
    getInfo() { return this.info; }
    async open(options) {
        if (this.busy) {
            const error = new Error('Failed to open serial port.');
            error.name = 'NetworkError';
            throw error;
        }
        this.opened.push(options);
        this.readable = new ReadableStream({ start: (controller) => { this.controller = controller; } });
        this.writable = new WritableStream({ write: (chunk) => { this.written.push(Array.from(chunk)); } });
    }
    async setSignals(signals) { this.signals.push(signals); }
    async close() { this.closed = (this.closed || 0) + 1; this.readable = null; this.writable = null; }
    push(bytes) { this.controller.enqueue(Uint8Array.from(bytes)); }
    fail(name) { const e = new Error(name); e.name = name; this.controller.error(e); }
}

test('Web Serial: open at 9600 with DTR, ordered writes, incoming data, close', async () => {
    const port = new FakeSerialPort();
    setNavigator({ serial: new EventTarget() });
    navigator.serial.getPorts = async () => [port, new FakeSerialPort({ usbVendorId: 0x2E8A })];

    assert.equal((await OF.WebSerialTransport.getKnownPorts()).length, 1);
    assert.equal(OF.WebSerialTransport.describePort(port).label, 'OpenFIRE (PID 0x0002)');

    const transport = new OF.WebSerialTransport(port);
    const received = [];
    const closes = [];
    transport.onData = (bytes) => received.push(...bytes);
    transport.onClose = (reason) => closes.push(reason);

    await transport.open();
    assert.deepEqual(port.opened[0], { baudRate: 9600, bufferSize: 4096 });
    assert.deepEqual(port.signals[0], { dataTerminalReady: true });

    await Promise.all([transport.write(Uint8Array.of(1, 2)), transport.write(Uint8Array.of(3)), transport.write(Uint8Array.of(4, 5))]);
    assert.deepEqual(port.written, [[1, 2], [3], [4, 5]]);

    port.push([0xA5, 0x5A]);
    port.push([0x01]);
    await sleep(10);
    assert.deepEqual(received, [0xA5, 0x5A, 0x01]);

    await transport.close();
    assert.equal(transport.isOpen, false);
    assert.equal(port.closed, 1);
    assert.deepEqual(closes, ['closed']);
    assert.equal(await transport.write(Uint8Array.of(9)), false);
});

test('Web Serial: busy port and device loss', async () => {
    setNavigator({ serial: new EventTarget() });
    const busyPort = new FakeSerialPort();
    busyPort.busy = true;
    await assert.rejects(new OF.WebSerialTransport(busyPort).open(), (error) => error.code === 'port_busy');

    const port = new FakeSerialPort();
    const transport = new OF.WebSerialTransport(port);
    const closes = [];
    transport.onClose = (reason) => closes.push(reason);
    await transport.open();
    port.fail('NetworkError');
    await sleep(20);
    assert.deepEqual(closes, ['device_lost']);
    assert.equal(transport.isOpen, false);
});

test('Web Serial: recoverable framing errors keep the port open', async () => {
    setNavigator({ serial: new EventTarget() });
    const port = new FakeSerialPort();
    const transport = new OF.WebSerialTransport(port);
    const closes = [];
    transport.onClose = (reason) => closes.push(reason);
    await transport.open();
    // A recoverable error replaces the stream, as Chromium does.
    const oldController = port.controller;
    port.readable = new ReadableStream({ start: (controller) => { port.controller = controller; } });
    const e = new Error('framing'); e.name = 'FramingError';
    oldController.error(e);
    await sleep(20);
    assert.equal(transport.isOpen, true);
    assert.deepEqual(closes, []);
    await transport.close();
});

test('Web Serial: 1200-baud bootloader touch reopens the port and drops DTR', async () => {
    setNavigator({ serial: new EventTarget() });
    const port = new FakeSerialPort();
    const transport = new OF.WebSerialTransport(port);
    await transport.open();
    await transport.touch1200();
    assert.deepEqual(port.opened[1], { baudRate: 1200 });
    assert.deepEqual(port.signals[port.signals.length - 1], { dataTerminalReady: false });
    assert.equal(port.closed, 2);
});

// ----- Fake WebSocket ----------------------------------------------------------

class FakeWebSocket {
    static get OPEN() { return 1; }
    constructor(url) {
        this.url = url;
        this.readyState = 0;
        this.sent = [];
        FakeWebSocket.last = this;
        if (!url.includes('refuse'))
            setTimeout(() => { this.readyState = 1; this.onopen && this.onopen(); }, 5);
        else
            setTimeout(() => { this.onerror && this.onerror(new Error('refused')); this.onclose && this.onclose(); }, 5);
    }
    send(data) { this.sent.push(Array.from(data)); }
    close() { this.readyState = 3; this.onclose && this.onclose(); }
}

test('WebSocket: default URL, binary messages, send and remote close', async () => {
    globalThis.WebSocket = FakeWebSocket;
    assert.equal(OF.WebSocketTransport.defaultUrl({ protocol: 'http:', host: '192.168.4.1' }), 'ws://192.168.4.1/ws');
    assert.equal(OF.WebSocketTransport.defaultUrl({ protocol: 'https:', host: 'openfire.local' }), 'wss://openfire.local/ws');

    const transport = new OF.WebSocketTransport('ws://192.168.4.1/ws');
    const received = [];
    const closes = [];
    transport.onData = (bytes) => received.push(...bytes);
    transport.onClose = (reason) => closes.push(reason);
    await transport.open();
    const socket = FakeWebSocket.last;
    assert.equal(socket.binaryType, 'arraybuffer');

    assert.equal(await transport.write(Uint8Array.of(7, 8)), true);
    assert.deepEqual(socket.sent, [[7, 8]]);

    socket.onmessage({ data: Uint8Array.of(0xA5, 0x5A).buffer });
    assert.deepEqual(received, [0xA5, 0x5A]);

    socket.readyState = 3;
    socket.onclose();
    assert.deepEqual(closes, ['connection_lost']);
    assert.equal(transport.isOpen, false);
});

test('WebSocket: refused connection rejects open', async () => {
    globalThis.WebSocket = FakeWebSocket;
    const transport = new OF.WebSocketTransport('ws://refuse/ws');
    await assert.rejects(transport.open(), (error) => error.code === 'open_failed');
});
