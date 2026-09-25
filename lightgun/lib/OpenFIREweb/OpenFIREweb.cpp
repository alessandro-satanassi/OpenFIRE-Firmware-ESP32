#include <Arduino.h>
#include "OpenFIREweb.h"
#include "OpenFIREusbnet.h"

#if defined(ARDUINO_ARCH_ESP32)
// Included before the Serial redefinition below, as in the original file.
#include "../../include/web_assets.h"
//#include <WiFi.h> //
//#include <esp_wifi.h> //
#include <esp_http_server.h>
#include <DNSServer.h>
#include <ESPmDNS.h>

#include <NetBIOS.h>

#include <unistd.h>
#include <lwip/sockets.h>
#include "../../src/OpenFIREserial.h"   // http://<gun>/status

#include "../../shared_lib/OpenFIRE_Wireless/ESP32/OpenFIRE_Wireless.h"

#endif // ARDUINO_ARCH_ESP32

bool OF_WebConfigModeActive = false;

// ============ [ESP32_PORT] ============
// redefinition of Serial to handle wireless serial connections / redifinizione di Serial per gestire le connessione wireless seriali
#ifdef OPENFIRE_WIRELESS_ENABLE
    extern Stream* Serial_OpenFIRE_Stream;
    #ifdef Serial
        #define AUX_SERIAL Serial
        #undef Serial
    #endif
    #define Serial (*Serial_OpenFIRE_Stream)
#endif // OPENFIRE_WIRELESS_ENABLE
// ============ [ESP32_PORT] ============
// END redefinition of Serial to handle wireless serial connections / fine redifinizione di Serial per gestire le connessione wireless seriali ========


// =================================================================================================
// --- SERIAL STREAM ENGINE (all architectures) / MOTORE SERIALE (tutte le architetture) ---
// Default App protocol link: the USB serial port, or the wireless stream through the dongle.
static int hw_available() { return Serial.available(); }
static int hw_read() { return Serial.read(); }
static size_t hw_writeByte(uint8_t c) { return Serial.write(c); }
static size_t hw_writeBuf(const uint8_t* buf, size_t size) { return Serial.write(buf, size); }
static void hw_flush() { Serial.flush(); }
static size_t hw_readBytes(char* buf, size_t size) { return Serial.readBytes(buf, size); }

static constexpr OF_WebSerialWrapper::SerialOps hardwareOps = {
    hw_available, hw_read, hw_writeByte, hw_writeBuf, hw_flush, hw_readBytes
};

OF_WebSerialWrapper::SerialOps WebAppSerial::ops = hardwareOps;
WebAppLink WebAppSerial::link = WebAppLink::SerialPort;
// =================================================================================================


#if defined(ARDUINO_ARCH_ESP32)


#ifndef WEBAPP_AP_SSID
    #define WEBAPP_AP_SSID "OpenFIRE_Config"
#endif



// Channel used when the radio is not already in use by the ESP-NOW link (cable
// only). Without an explicit build flag it is the channel of that link, so a
// dongle or a wireless pedal connecting later finds the radio already tuned to
// it. OPENFIRE_ESPNOW_WIFI_CHANNEL is defined inside OpenFIRE_Wireless.cpp and
// is not visible here, but the library exports the channel in use as a variable
// (and that one also holds a channel negotiated at run time).
// /
// Canale usato quando la radio non e' gia' impegnata dal collegamento ESP-NOW
// (solo cavo). Senza un flag di compilazione esplicito e' il canale di quel
// collegamento, cosi' un dongle o un pedale wireless che si connette dopo trova
// la radio gia' sintonizzata. OPENFIRE_ESPNOW_WIFI_CHANNEL e' definita dentro
// OpenFIRE_Wireless.cpp e qui non si vede, ma la libreria espone il canale in
// uso come variabile (che contiene anche un canale negoziato a run time).
/*
#ifndef WEBAPP_AP_DEFAULT_CHANNEL
    #define WEBAPP_AP_CHANNEL_FROM_LINK
    #define WEBAPP_AP_DEFAULT_CHANNEL 1
#endif
*/

/*
#if defined(WEBAPP_AP_CHANNEL_FROM_LINK) && defined(OPENFIRE_WIRELESS_ENABLE)
extern uint8_t espnow_wifi_channel;   // OpenFIRE_Wireless.h
#endif
*/

static httpd_handle_t web_server = NULL;
static DNSServer dnsServer;

// =================================================================================================
// --- WEBSOCKET ENGINE / MOTORE WEBSOCKET ---
// One App page at a time. The HTTP server task produces incoming bytes, the
// firmware loop consumes them (single-producer / single-consumer ring buffer).

#define WS_RX_BUFFER_SIZE 1024   // a few App frames (max 207 bytes each)
#define WS_MAX_MESSAGE    512    // larger WebSocket messages are refused
static uint8_t ws_rx_buffer[WS_RX_BUFFER_SIZE];
// Counters of bytes, not positions: they never wrap back to the start of the buffer,
// so "is this mark still ahead of what has been read?" is a plain comparison and stays
// right even when the 32-bit counters roll over (the difference is always < the buffer).
// The slot in the buffer is the counter modulo its size.
//
// Everything the two tasks share goes through the two functions below, and nothing else:
// that is the whole rule. `volatile` would not be enough. What matters here is that the
// byte reaches the buffer BEFORE the counter that announces it, and that the byte is read
// BEFORE its slot is released - and the buffer is an ordinary array, which the compiler is
// free to move around a volatile access. GCC says so in as many words: a volatile object
// cannot be used as a memory barrier to order writes to non-volatile memory. It is not a
// theoretical freedom: with volatile counters, GCC 13 on x86-64 at -O2 emits the store to
// the counter BEFORE the store of the byte, and the load of the byte AFTER the slot has
// been released - and nothing in the language stops the compiler for this target from
// doing the same. A publish/read pair with release/acquire says what is meant, and the
// compiler keeps it.
//
// Only whole aligned words are published, and only with plain loads and stores: never a
// read-modify-write, which on this target could turn into a library call with a lock
// inside. To check that none is left:
//   xtensa-esp32s3-elf-nm <build>/OpenFIREweb.cpp.o | grep __atomic
// Nothing listed means the compiler inlined them all. (Do not assert it with
// __atomic_always_lock_free(4, 0): with a null pointer it answers for the target's typical
// alignment, and it answers about lock-free atomic operations in general - not only the
// plain loads and stores used here. On this toolchain it answers no, so the assertion just
// breaks the build without saying anything about this code: read the object file instead.)

// Widths are spelled out on purpose: int32_t and uint32_t, never a bare int.
// On this target int is already 32 bits wide, so int32_t is its exact equivalent and
// nothing here was narrowed; the explicit spelling only says out loud how wide the
// value is, which matters for the shared counters above. Where a type is imposed from
// outside - SerialOps, an ESP-IDF callback, a printf format - it is left as it is and
// the conversion is written on the spot, so it is visible instead of implied.

// Read a value the other task may be publishing right now.
static inline uint32_t sharedGet(const uint32_t *value) {
    return __atomic_load_n(value, __ATOMIC_ACQUIRE);
}
static inline int32_t sharedGet(const int32_t *value) {
    return __atomic_load_n(value, __ATOMIC_ACQUIRE);
}
// Publish a value for the other task, after whatever it announces is already in place.
// Two plain overloads and not one template: the type the pointer points at picks the
// function, and the value converts to it. A template would have to deduce the same type
// from both, and uint32_t is not spelled the same everywhere - unsigned long here,
// unsigned int elsewhere - so a literal written one way would stop matching.
static inline void sharedPublish(uint32_t *value, uint32_t now) {
    __atomic_store_n(value, now, __ATOMIC_RELEASE);
}
static inline void sharedPublish(int32_t *value, int32_t now) {
    __atomic_store_n(value, now, __ATOMIC_RELEASE);
}

static uint32_t ws_rx_written = 0;   // written by the HTTP server task only
static uint32_t ws_rx_read = 0;      // written by the firmware loop only

// Socket of the current App page (-1 = none). Changed in the server task, read by the
// sender in the firmware loop as well, so it follows the same rule as everything else.
// The server takes it as a plain int, so it is converted back where it is handed over.
static int32_t ws_client_fd = -1;

// Set by the server task, consumed by the firmware loop (0 = no, 1 = yes).
static uint32_t ws_client_lost = 0;     // the page that owned the session is gone
// A page that simply closed (no new page taking over) is reported only once what it
// sent has been read: its last request - the reboot to the bootloader, for instance -
// is still in the buffer, and ending the session first would throw it away.
static uint32_t ws_client_closed = 0;   // the page closed, nobody replaced it
static uint32_t ws_client_closed_at = 0;
#define WS_CLOSE_DRAIN_MS 250                    // ...at the latest (an unfinished frame)
// A new page asks for the bytes left by the previous one to be dropped. The request
// carries its own number instead of a flag to be cleared: reading the number means the
// mark published with it is in place too, and a request arriving while this one is being
// handled keeps a number of its own, so it cannot be swallowed. A flag would have to be
// cleared, and the clear can land after a newer request without ever having seen it -
// leaving the flag down and that request lost.
static uint32_t ws_flush_seq = 0;       // published by the HTTP server task
static uint32_t ws_flush_to = 0;        // ...up to here: what the new page sent is kept
static uint32_t ws_flush_done = 0;      // what the firmware loop has already applied (its own)

// Counters of http://<gun>/status (diagnosis of the link with the App page).
// Plain counters: only read for the diagnosis, an exact count is not needed
// (C++20 deprecates ++ on a volatile).
// The first three are written and read in the server task alone. The last three are
// written in the firmware loop and read by /status in the server task, so they are
// published like everything else: an exact count is not needed, but the accesses still
// have to be proper ones.
static uint32_t ws_stat_handshakes = 0;  // pages that opened the WebSocket
static uint32_t ws_stat_rx = 0;          // bytes received from the page
static uint32_t ws_stat_dropped = 0;     // bytes discarded (buffer full, or not the current page)
static uint32_t ws_stat_tx = 0;          // bytes sent to the page
static uint32_t ws_stat_tx_failed = 0;   // sends the server refused
static uint32_t ws_stat_flushed = 0;     // bytes dropped when a new page arrived

// The request is published by the server task in two steps (the mark, then its number) and
// consumed here in several more: the two tasks can interleave, so this must not rely on the
// number and the mark being written together. It does not. The number is read first, and it
// is published after the mark, so the mark that belongs to it is already in place. The mark
// is then applied only while it is still ahead of what has been read, so an older request
// can only ask to drop bytes that are already gone, which is nothing at all. A request
// arriving while this runs carries a number of its own, so the next call sees it rather
// than losing it - which is what a flag to be cleared could not guarantee.
// Moving the read counter backwards would hand the same byte to the protocol twice.
static void ws_handle_flush() {
    const uint32_t seq = sharedGet(&ws_flush_seq);
    if (seq == ws_flush_done) return;
    // The number is published after the mark, so reading it means the mark is in place.
    const uint32_t mark = sharedGet(&ws_flush_to);
    ws_flush_done = seq;
    const uint32_t alreadyRead = sharedGet(&ws_rx_read);
    if ((int32_t)(mark - alreadyRead) > 0) {
        sharedPublish(&ws_stat_flushed, sharedGet(&ws_stat_flushed) + (mark - alreadyRead));
        sharedPublish(&ws_rx_read, mark);
    }
}

// available() and read() answer with a plain int because that is what SerialOps declares,
// following the Arduino Stream convention (-1 = nothing there). Not our choice to make.
static int ws_available() {
    ws_handle_flush();
    const uint32_t alreadyRead = sharedGet(&ws_rx_read);
    return (int)(sharedGet(&ws_rx_written) - alreadyRead);
}

static int ws_read() {
    ws_handle_flush();
    // Snapshot the counter: this task is its only writer. The slot is released only
    // after the byte has been taken out of it.
    const uint32_t alreadyRead = sharedGet(&ws_rx_read);
    if (sharedGet(&ws_rx_written) == alreadyRead) return -1;
    uint8_t c = ws_rx_buffer[alreadyRead % WS_RX_BUFFER_SIZE];
    sharedPublish(&ws_rx_read, alreadyRead + 1);
    return c;
}

static size_t ws_writeBuf(const uint8_t* buf, size_t size) {
    const int32_t fd = sharedGet(&ws_client_fd);
    if (!web_server || fd < 0) return 0;

    httpd_ws_frame_t ws_pkt;
    memset(&ws_pkt, 0, sizeof(httpd_ws_frame_t));
    ws_pkt.type = HTTPD_WS_TYPE_BINARY; // The whole protocol is binary / tutto il protocollo è binario
    ws_pkt.payload = (uint8_t*)buf;
    ws_pkt.len = size;

    if (httpd_ws_send_data(web_server, (int)fd, &ws_pkt) != ESP_OK) {
        sharedPublish(&ws_stat_tx_failed, sharedGet(&ws_stat_tx_failed) + 1u);
        return 0;
    }
    sharedPublish(&ws_stat_tx, sharedGet(&ws_stat_tx) + (uint32_t)size);
    return size;
}

static size_t ws_writeByte(uint8_t c) {
    return ws_writeBuf(&c, 1);
}

static void ws_flush() {
    // WebSocket frames are sent whole: there is no TX buffer to flush.
}

static size_t ws_readBytes(char* buf, size_t size) {
    size_t count = 0;
    const uint32_t startMillis = (uint32_t)millis();
    // 1000 ms timeout, emulating the standard Arduino Serial behaviour.
    while (count < size && ((uint32_t)millis() - startMillis < 1000u)) {
        if (ws_available() > 0) {
            // A byte counted a moment ago can still fail to arrive: ws_read() handles a
            // pending flush first, and a page that has just taken over can carry the read
            // counter past it. Storing that -1 would invent a 0xFF that nobody sent.
            const int32_t incoming = ws_read();
            if (incoming >= 0)
                buf[count++] = (char)incoming;
        } else {
            vTaskDelay(pdMS_TO_TICKS(1));
        }
    }
    return count;
}

static constexpr OF_WebSerialWrapper::SerialOps websocketOps = {
    ws_available, ws_read, ws_writeByte, ws_writeBuf, ws_flush, ws_readBytes
};

bool WebApp_TakeClientLost() {
    // A new page took over: immediate, its predecessor has nothing left to say.
    if (sharedGet(&ws_client_lost)) {
        sharedPublish(&ws_client_lost, 0u);
        sharedPublish(&ws_client_closed, 0u);
        return true;
    }
    if (!sharedGet(&ws_client_closed)) return false;
    // The page just closed: let the firmware read what it sent before leaving. Reading the
    // flag set means the time of that close is in place too.
    if (ws_available() > 0 &&
        (uint32_t)(millis() - sharedGet(&ws_client_closed_at)) < WS_CLOSE_DRAIN_MS)
        return false;
    sharedPublish(&ws_client_closed, 0u);
    return true;
}

bool WebApp_ClientLostPending() {
    return sharedGet(&ws_client_lost) || sharedGet(&ws_client_closed);
}

void WebApp_RadioState(uint8_t *channel, uint8_t *powerSave) { // DA TOGLIERE
    if (channel) *channel = 0;
    if (powerSave) *powerSave = 0;
    /*
    if (channel) {
        uint8_t primary = 0;
        wifi_second_chan_t second = WIFI_SECOND_CHAN_NONE;
        *channel = esp_wifi_get_channel(&primary, &second) == ESP_OK ? primary : 0;
    }
    if (powerSave) {
        wifi_ps_type_t ps = WIFI_PS_NONE;
        *powerSave = esp_wifi_get_ps(&ps) == ESP_OK ? (uint8_t)ps : 0;
    }
    */
}

// =================================================================================================
// --- HTTP SERVER / SERVER HTTP ---

// Captive portal: every unknown URL (OS connectivity checks included) goes to the App page.
// The socket is closed right after the redirect: phones repeat these checks every few
// seconds and their idle connections would otherwise use up the server's sockets
// (and the App WebSocket would be the one closed to make room).
static esp_err_t captive_portal_handler(httpd_req_t *req, httpd_err_code_t error) {
    const char *location = "http://192.168.4.1/";
    #ifdef OPENFIRE_USB_NCM
    // A request received over USB must not be redirected to the Wi-Fi address.
    struct sockaddr_storage local = {};
    socklen_t localSize = sizeof(local);
    if(OpenFIREUsbNetActive() &&
       getsockname(httpd_req_to_sockfd(req), (struct sockaddr*)&local, &localSize) == 0) {
        uint32_t address = 0;
        if(local.ss_family == AF_INET)
            address = ((struct sockaddr_in*)&local)->sin_addr.s_addr;
        #if LWIP_IPV6
        else if(local.ss_family == AF_INET6) {
            const struct in6_addr *v6 = &((struct sockaddr_in6*)&local)->sin6_addr;
            if(IN6_IS_ADDR_V4MAPPED(v6)) memcpy(&address, &v6->s6_addr[12], sizeof(address));
        }
        #endif
        if(address == inet_addr("192.168.7.1")) location = "http://192.168.7.1/";
    }
    #endif
    httpd_resp_set_status(req, "302 Found");
    httpd_resp_set_hdr(req, "Location", location);
    httpd_resp_set_hdr(req, "Connection", "close");
    httpd_resp_send(req, NULL, 0);
    if (web_server)
        httpd_sess_trigger_close(web_server, httpd_req_to_sockfd(req));
    return ESP_OK;
}

static esp_err_t send_gzip(httpd_req_t *req, const char *type, const uint8_t *data, size_t len) {
    httpd_resp_set_type(req, type);
    httpd_resp_set_hdr(req, "Content-Encoding", "gzip");
    // The assets change with every firmware build.
    httpd_resp_set_hdr(req, "Cache-Control", "no-cache");
    return httpd_resp_send(req, (const char*)data, len);
}

static esp_err_t index_get_handler(httpd_req_t *req) {
    return send_gzip(req, "text/html", web_index_html_gz, web_index_html_gz_len);
}

static esp_err_t style_get_handler(httpd_req_t *req) {
    return send_gzip(req, "text/css", web_style_css_gz, web_style_css_gz_len);
}

static esp_err_t app_js_get_handler(httpd_req_t *req) {
    return send_gzip(req, "application/javascript", web_app_js_gz, web_app_js_gz_len);
}

// http://<gun>/status: state of the link with the App page (diagnosis, no secrets).
//
// "docked" and "link" are the only two values here that are not published: they are read
// straight from the firmware loop's own variables (appSerialSessionActive, and the link
// WebAppSerial::Use() selects), which it writes without synchronising. That is left as it
// is, on purpose. Both are a single byte, so a read can never see anything but one of
// their valid values - but they are read one after the other, so this page can show a
// combination that never existed in one instant, such as docked=1 with link="serial"
// while a page is docking. Nothing else follows from it: the two values are printed here
// and nowhere else, and everything that the protocol actually acts on is published
// properly. Making them agree would mean publishing a snapshot from the firmware loop,
// which is a new path between tasks - more machinery, and more to get wrong, than two
// diagnostic fields being a step apart during a handover are worth.
static esp_err_t status_get_handler(httpd_req_t *req) {
    char body[400];
    uint8_t radioChannel = 0, radioPowerSave = 0;
    WebApp_RadioState(&radioChannel, &radioPowerSave); // DA TOGLIERE
    const int32_t length = snprintf(body, sizeof(body),
        "{\"webConfig\":%d,\"docked\":%d,\"link\":\"%s\",\"wsClient\":%d,"
        "\"handshakes\":%u,\"rx\":%u,\"tx\":%u,\"txFailed\":%u,\"dropped\":%u,\"flushed\":%u,"
        "\"pending\":%d,"
        // Radio: channel in use and power saving (0 = none, as ESP-NOW needs).
        "\"channel\":%u,\"powerSave\":%u,"
        "\"uptimeMs\":%lu}",
        OF_WebConfigModeActive ? 1 : 0,
        OF_Serial::AppSerialSessionIsActive() ? 1 : 0,
        WebAppSerial::IsWebSocket() ? "ws" : "serial",
        (int)sharedGet(&ws_client_fd),
        (unsigned)ws_stat_handshakes, (unsigned)ws_stat_rx, (unsigned)sharedGet(&ws_stat_tx),
        (unsigned)sharedGet(&ws_stat_tx_failed), (unsigned)ws_stat_dropped,
        (unsigned)sharedGet(&ws_stat_flushed),
        (int)(sharedGet(&ws_rx_written) - sharedGet(&ws_rx_read)), // without consuming a pending flush
        (unsigned)radioChannel, (unsigned)radioPowerSave,
        (unsigned long)millis());
    // snprintf returns the length the text WOULD have had: should the diagnosis ever
    // outgrow the buffer, sending that number would read past it and put whatever the
    // stack holds on the wire. Send what was actually written. Nothing changes as long
    // as it fits, which today it does with room to spare.
    const size_t sent = (length > 0)
        ? ((size_t)length < sizeof(body) ? (size_t)length : sizeof(body) - 1)
        : 0;
    httpd_resp_set_type(req, "application/json");
    httpd_resp_set_hdr(req, "Cache-Control", "no-store");
    return httpd_resp_send(req, body, sent);
}

// A new App page owns the link: the previous one (if any) is closed and its
// session abandoned, and the bytes received so far are dropped (not the ones of
// this page, which arrive after the mark).
static void ws_adopt_client(int32_t fd) {
    const int32_t previous = sharedGet(&ws_client_fd);
    ws_stat_handshakes++;
    sharedPublish(&ws_client_fd, fd);
    if (previous >= 0 && previous != fd)
        httpd_sess_trigger_close(web_server, (int)previous);
    // The mark first, its number after: the loop reads the number and knows the mark is there.
    sharedPublish(&ws_flush_to, sharedGet(&ws_rx_written));
    sharedPublish(&ws_flush_seq, sharedGet(&ws_flush_seq) + 1u);
    sharedPublish(&ws_client_lost, 1u);
}

static esp_err_t ws_handler(httpd_req_t *req) {
    const int32_t fd = httpd_req_to_sockfd(req);

    // Handshake of a new page. Some esp_http_server versions do not report it
    // here: the first message of an unknown socket (below) does the same.
    if (req->method == HTTP_GET) {
        ws_adopt_client(fd);
        return ESP_OK;
    }

    httpd_ws_frame_t ws_pkt;
    memset(&ws_pkt, 0, sizeof(httpd_ws_frame_t));
    ws_pkt.type = HTTPD_WS_TYPE_BINARY;

    // 1. Ask the server how long the incoming packet is
    esp_err_t ret = httpd_ws_recv_frame(req, &ws_pkt, 0);
    if (ret != ESP_OK || ws_pkt.len == 0) return ret;

    // Not an App protocol message: refuse it (the server closes the socket)
    if (ws_pkt.len > WS_MAX_MESSAGE) return ESP_ERR_INVALID_SIZE;

    uint8_t *buf = (uint8_t*)malloc(ws_pkt.len);
    if (!buf) return ESP_ERR_NO_MEM;
    ws_pkt.payload = buf;

    // 2. Receive the data, 3. queue it: a message from another socket is a new page
    ret = httpd_ws_recv_frame(req, &ws_pkt, ws_pkt.len);
    if (ret == ESP_OK && fd != sharedGet(&ws_client_fd))
        ws_adopt_client(fd);
    if (ret == ESP_OK && fd == sharedGet(&ws_client_fd)) {
        for (size_t i = 0; i < ws_pkt.len; i++) {
            // Snapshot the counter: this task is its only writer. The byte is announced
            // only after it is in the buffer.
            const uint32_t alreadyWritten = sharedGet(&ws_rx_written);
            if (alreadyWritten - sharedGet(&ws_rx_read) < WS_RX_BUFFER_SIZE) { // if there is room / se c'è spazio
                ws_rx_buffer[alreadyWritten % WS_RX_BUFFER_SIZE] = buf[i];
                sharedPublish(&ws_rx_written, alreadyWritten + 1);
                ws_stat_rx++;
            } else ws_stat_dropped++;
        }
    } else if (ret == ESP_OK) {
        ws_stat_dropped += ws_pkt.len; // frame of a page that is no longer the current one
    }
    free(buf);
    return ret;
}

// Called by the server for every closed socket; it must close the socket itself.
static void web_close_fn(httpd_handle_t hd, int sockfd) {
    if ((int32_t)sockfd == sharedGet(&ws_client_fd)) {
        sharedPublish(&ws_client_fd, -1);
        // The time first, the flag second. The firmware loop reads the time as soon as
        // it sees the flag, to decide how long to wait for the last bytes of the page;
        // catching the flag before the time was written, it would use the time of some
        // earlier close, find the wait long over, and drop what the page just sent -
        // its request to reboot into flashing mode, for instance.
        sharedPublish(&ws_client_closed_at, (uint32_t)millis());
        sharedPublish(&ws_client_closed, 1u);
    }
    close(sockfd);
}

static const httpd_uri_t uri_index = { .uri = "/", .method = HTTP_GET, .handler = index_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_style = { .uri = "/style.css", .method = HTTP_GET, .handler = style_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_app   = { .uri = "/app.js", .method = HTTP_GET, .handler = app_js_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_status = { .uri = "/status", .method = HTTP_GET, .handler = status_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_ws    = { .uri = "/ws", .method = HTTP_GET, .handler = ws_handler, .user_ctx = NULL, .is_websocket = true };

// Background task for the captive portal DNS requests
static void dns_server_task(void *pvParameters) {
    while (true) {
        if (OF_WebConfigModeActive) {
            dnsServer.processNextRequest();
        }
        vTaskDelay(pdMS_TO_TICKS(10)); // 10 ms pause so the CPU is not blocked
    }
}

void WebApp_Init() {
    if (!OF_WebConfigModeActive) return;
        //else return;

    // The App protocol keeps using the serial link until an App docks on the
    // WebSocket (OF_Serial::SerialProcessingWebDock selects the link).

    // 1. Access point. The App page over WiFi and the ESP-NOW link (dongle or
    //    wireless pedal) work together, but a single radio has a single channel:
    //    when that link is already running, the access point joins it on ITS
    //    channel and its radio settings are put back afterwards (starting the
    //    WiFi of the Arduino layer resets some of them). With the cable alone the
    //    radio is free: the channel of the ESP-NOW link is used anyway (see the
    //    define above), so a later connection finds it already tuned.
    //    The channel is asked to the radio itself: the ESP-NOW link starts it
    //    with the IDF API, without going through the Arduino WiFi class, whose
    //    mode would still read OFF.
    
    /*
    uint8_t apChannel = WEBAPP_AP_DEFAULT_CHANNEL;
    #if defined(WEBAPP_AP_CHANNEL_FROM_LINK) && defined(OPENFIRE_WIRELESS_ENABLE)
        if (espnow_wifi_channel >= 1 && espnow_wifi_channel <= 13)
            apChannel = espnow_wifi_channel;
    #endif
    */
    /*
    bool linkRunning = false;
    {
        uint8_t primary = 0;
        wifi_second_chan_t second = WIFI_SECOND_CHAN_NONE;
        if (esp_wifi_get_channel(&primary, &second) == ESP_OK && primary >= 1 && primary <= 13) {
            apChannel = primary;
            linkRunning = true;
        }
    }
    */
    
    // With the ESP-NOW link already configured and running, nothing of the radio
    // is touched: only the access point interface is added, on the channel that
    // link is already on. Protocol, bandwidth and power stay as the wireless
    // library set them (the Arduino layer only rewrites them for Long Range,
    // which is not used here), and the access point configuration is written once.
    //
    // The one exception is power saving: every time the station interface starts,
    // the Arduino layer applies ITS OWN setting (esp_wifi_set_ps(WiFi.getSleep()),
    // WIFI_PS_MIN_MODEM by default) in the STA_START event handler, undoing the
    // WIFI_PS_NONE the wireless library needs - and a modem that sleeps loses the
    // ESP-NOW packets it should receive. Saying it here beforehand is not a change
    // of configuration: it is the same value the library sets, so that the Arduino
    // layer stops putting its own back.
    // /
    // Con il collegamento ESP-NOW gia' configurato e attivo non si tocca nulla
    // della radio: si aggiunge solo l'interfaccia access point, sul canale su cui
    // quel collegamento si trova gia'.
    //
    // L'unica eccezione e' il power save: a ogni avvio dell'interfaccia station il
    // livello Arduino applica la PROPRIA impostazione (esp_wifi_set_ps con
    // WiFi.getSleep(), di default WIFI_PS_MIN_MODEM) nel gestore dell'evento
    // STA_START, annullando il WIFI_PS_NONE che serve alla libreria wireless - e un
    // modem che dorme perde i pacchetti ESP-NOW in ricezione. Dirglielo qui prima
    // non cambia la configurazione: e' lo stesso valore impostato dalla libreria.
    
    //WiFi.setSleep(false);   // WIFI_PS_NONE, come la libreria wireless
    
    
    //WiFi.persistent(false);
    //WiFi.mode(WIFI_AP_STA);
    //WiFi.softAPdisconnect();
    //esp_err_t err;
    //err = esp_wifi_set_protocol(WIFI_IF_AP, WIFI_PROTOCOL_11G);
    //err = esp_wifi_set_bandwidth(WIFI_IF_AP, WIFI_BW_HT20);
    //err = esp_wifi_set_ps(WIFI_PS_NONE); // non dovrebbe servire
    //WiFi.softAP(WEBAPP_AP_SSID, NULL, apChannel);
    //WiFi.softAP(WEBAPP_AP_SSID, NULL, 13);
    //esp_err_t err;
    //err = esp_wifi_set_protocol(WIFI_IF_AP, WIFI_PROTOCOL_11G);
    //err = esp_wifi_set_bandwidth(WIFI_IF_AP, WIFI_BW_HT20);
    SerialWireless.startAccessPoint(WEBAPP_AP_SSID);   // rete aperta, canale di ESP-NOW

    // 2. Lightweight native ESP-IDF HTTP/WebSocket server
    httpd_config_t config = HTTPD_DEFAULT_CONFIG();
    config.max_uri_handlers = 8;
    config.close_fn = web_close_fn;
    // Room for the App page, its files and the connectivity checks of a phone;
    // without it the oldest socket is closed to make room, App WebSocket included.
    config.max_open_sockets = 7;
    config.lru_purge_enable = true;

    if (httpd_start(&web_server, &config) == ESP_OK) {
        httpd_register_uri_handler(web_server, &uri_index);
        httpd_register_uri_handler(web_server, &uri_style);
        httpd_register_uri_handler(web_server, &uri_app);
        httpd_register_uri_handler(web_server, &uri_status);
        httpd_register_uri_handler(web_server, &uri_ws);

        // Captive portal for everything else
        httpd_register_err_handler(web_server, HTTPD_404_NOT_FOUND, captive_portal_handler);
    }

    // 3. mDNS: http://openfire.local
    if (MDNS.begin("openfire")) {
        MDNS.addService("http", "tcp", 80);
        #ifdef OPENFIRE_USB_NCM
            OpenFIREUsbNetMDNSReady(); // Same mDNS service on Wi-Fi and USB.
        #endif
    }

    NBNS.begin("openfire");
    
    // 4. DNS server for the captive portal
    dnsServer.start(53, "*", SerialWireless.ipAddressAP());
    //dnsServer.start(53, "*", WiFi.softAPIP());
    

    /*
    // 4. DNS server for the captive portal
    IPAddress dnsIP = SerialWireless.ipAddressAP();
    if (dnsIP == IPAddress(0, 0, 0, 0)) {
        dnsIP = IPAddress(192, 168, 7, 1);
    }
    dnsServer.start(53, "*", dnsIP);
    */

    // 5. DNS task on Core 0 so it does not interfere with the main loop
    xTaskCreatePinnedToCore(dns_server_task, "dns_task", 2048, NULL, 1, NULL, 0);

    /*
    // With no ESP-NOW link there is nothing to preserve: power saving off, so the
    // App page is as responsive as it is on the cable.
    if (!linkRunning)
        esp_wifi_set_ps(WIFI_PS_NONE);
    */
}

void WebApp_Loop() {
    // Nothing to do: the server and the DNS run in their own tasks.
}

static const OF_WebSerialWrapper::SerialOps& webSocketLinkOps() { return websocketOps; }

#else
// RP2040: no web configuration mode; the App protocol always uses Serial.
void WebApp_Init() {}
void WebApp_Loop() {}
bool WebApp_TakeClientLost() { return false; }
bool WebApp_ClientLostPending() { return false; }
void WebApp_RadioState(uint8_t *channel, uint8_t *powerSave) { //DA TOGLIERE
    if (channel) *channel = 0;
    if (powerSave) *powerSave = 0;
}

static int none_available() { return 0; }
static int none_read() { return -1; }
static size_t none_writeByte(uint8_t) { return 0; }
static size_t none_writeBuf(const uint8_t*, size_t) { return 0; }
static void none_flush() {}
static size_t none_readBytes(char*, size_t) { return 0; }
static constexpr OF_WebSerialWrapper::SerialOps noLinkOps = {
    none_available, none_read, none_writeByte, none_writeBuf, none_flush, none_readBytes
};
static const OF_WebSerialWrapper::SerialOps& webSocketLinkOps() { return noLinkOps; }
#endif

const OF_WebSerialWrapper::SerialOps& WebAppSerial::Ops(WebAppLink which)
{
    return which == WebAppLink::WebSocket ? webSocketLinkOps() : hardwareOps;
}

void WebAppSerial::Use(WebAppLink newLink)
{
    ops = Ops(newLink);
    link = newLink;
}


// ============ [ESP32_PORT] ============
// restore Serial after it was redefined for serial connections / ripristino di Serial dopo definizione per connessione seriali ==============
#ifdef OPENFIRE_WIRELESS_ENABLE
    #undef Serial
    #ifdef AUX_SERIAL
        #define Serial AUX_SERIAL
        #undef AuxSerial
    #endif
#endif // OPENFIRE_WIRELESS_ENABLE
// ============ [ESP32_PORT] ============
// restore Serial after it was redefined for serial connections / fine ripristino di Serial dopo definizione per connessione seriali ==============
