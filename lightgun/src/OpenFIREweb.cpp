#include <Arduino.h>
#include "OpenFIREweb.h"

#if defined(ARDUINO_ARCH_ESP32)
// Included before the Serial redefinition below, as in the original file.
#include "web_assets.h"
#include <WiFi.h>
#include <esp_wifi.h>
#include <esp_http_server.h>
#include <DNSServer.h>
#include <ESPmDNS.h>
#include <unistd.h>
#include "OpenFIREserial.h"   // http://<gun>/status
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
#ifndef WEBAPP_AP_PASSWORD
    #define WEBAPP_AP_PASSWORD "12345678"
#endif
// Channel used when the radio is not already in use by the ESP-NOW link.
#ifndef WEBAPP_AP_DEFAULT_CHANNEL
    #define WEBAPP_AP_DEFAULT_CHANNEL 1
#endif

static httpd_handle_t web_server = NULL;
static DNSServer dnsServer;

// =================================================================================================
// --- WEBSOCKET ENGINE / MOTORE WEBSOCKET ---
// One App page at a time. The HTTP server task produces incoming bytes, the
// firmware loop consumes them (single-producer / single-consumer ring buffer).

#define WS_RX_BUFFER_SIZE 1024   // a few App frames (max 207 bytes each)
#define WS_MAX_MESSAGE    512    // larger WebSocket messages are refused
static uint8_t ws_rx_buffer[WS_RX_BUFFER_SIZE];
static volatile uint16_t ws_rx_head = 0;   // written by the HTTP server task only
static volatile uint16_t ws_rx_tail = 0;   // written by the firmware loop only

// Socket of the current App page (-1 = none). Changed only in the server task.
static volatile int ws_client_fd = -1;

// Set by the server task, consumed by the firmware loop.
static volatile bool ws_client_lost = false;     // the page that owned the session is gone
static volatile bool ws_flush_request = false;   // drop bytes left by a previous page
static volatile uint16_t ws_flush_to = 0;        // ...up to here: what the new page sent is kept

// Counters of http://<gun>/status (diagnosis of the link with the App page).
// Plain counters: only read for the diagnosis, an exact count is not needed
// (C++20 deprecates ++ on a volatile).
static uint32_t ws_stat_handshakes = 0;  // pages that opened the WebSocket
static uint32_t ws_stat_rx = 0;          // bytes received from the page
static uint32_t ws_stat_dropped = 0;     // bytes discarded (buffer full, or not the current page)
static uint32_t ws_stat_tx = 0;          // bytes sent to the page
static uint32_t ws_stat_tx_failed = 0;   // sends the server refused
static uint32_t ws_stat_flushed = 0;     // bytes dropped when a new page arrived

static void ws_handle_flush() {
    if (ws_flush_request) {
        ws_flush_request = false;
        const uint16_t mark = ws_flush_to;
        ws_stat_flushed += (WS_RX_BUFFER_SIZE + mark - ws_rx_tail) % WS_RX_BUFFER_SIZE;
        ws_rx_tail = mark;
    }
}

static int ws_available() {
    ws_handle_flush();
    return (WS_RX_BUFFER_SIZE + ws_rx_head - ws_rx_tail) % WS_RX_BUFFER_SIZE;
}

static int ws_read() {
    ws_handle_flush();
    if (ws_rx_head == ws_rx_tail) return -1;
    uint8_t c = ws_rx_buffer[ws_rx_tail];
    ws_rx_tail = (ws_rx_tail + 1) % WS_RX_BUFFER_SIZE;
    return c;
}

static size_t ws_writeBuf(const uint8_t* buf, size_t size) {
    const int fd = ws_client_fd;
    if (!web_server || fd < 0) return 0;

    httpd_ws_frame_t ws_pkt;
    memset(&ws_pkt, 0, sizeof(httpd_ws_frame_t));
    ws_pkt.type = HTTPD_WS_TYPE_BINARY; // The whole protocol is binary / tutto il protocollo è binario
    ws_pkt.payload = (uint8_t*)buf;
    ws_pkt.len = size;

    if (httpd_ws_send_data(web_server, fd, &ws_pkt) != ESP_OK) {
        ws_stat_tx_failed++;
        return 0;
    }
    ws_stat_tx += size;
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
    unsigned long startMillis = millis();
    // 1000 ms timeout, emulating the standard Arduino Serial behaviour.
    while (count < size && (millis() - startMillis < 1000)) {
        if (ws_available() > 0) {
            buf[count++] = (char)ws_read();
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
    if (!ws_client_lost) return false;
    ws_client_lost = false;
    return true;
}

bool WebApp_ClientLostPending() {
    return ws_client_lost;
}

// =================================================================================================
// --- HTTP SERVER / SERVER HTTP ---

// Captive portal: every unknown URL (OS connectivity checks included) goes to the App page.
// The socket is closed right after the redirect: phones repeat these checks every few
// seconds and their idle connections would otherwise use up the server's sockets
// (and the App WebSocket would be the one closed to make room).
static esp_err_t captive_portal_handler(httpd_req_t *req, httpd_err_code_t error) {
    httpd_resp_set_status(req, "302 Found");
    httpd_resp_set_hdr(req, "Location", "http://192.168.4.1/");
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
static esp_err_t status_get_handler(httpd_req_t *req) {
    char body[320];
    const int length = snprintf(body, sizeof(body),
        "{\"webConfig\":%d,\"docked\":%d,\"link\":\"%s\",\"wsClient\":%d,"
        "\"handshakes\":%u,\"rx\":%u,\"tx\":%u,\"txFailed\":%u,\"dropped\":%u,\"flushed\":%u,"
        "\"pending\":%d,\"uptimeMs\":%lu}",
        OF_WebConfigModeActive ? 1 : 0,
        OF_Serial::AppSerialSessionIsActive() ? 1 : 0,
        WebAppSerial::IsWebSocket() ? "ws" : "serial",
        ws_client_fd,
        (unsigned)ws_stat_handshakes, (unsigned)ws_stat_rx, (unsigned)ws_stat_tx,
        (unsigned)ws_stat_tx_failed, (unsigned)ws_stat_dropped, (unsigned)ws_stat_flushed,
        (int)((WS_RX_BUFFER_SIZE + ws_rx_head - ws_rx_tail) % WS_RX_BUFFER_SIZE), // without consuming a pending flush
        (unsigned long)millis());
    httpd_resp_set_type(req, "application/json");
    httpd_resp_set_hdr(req, "Cache-Control", "no-store");
    return httpd_resp_send(req, body, length > 0 ? length : 0);
}

// A new App page owns the link: the previous one (if any) is closed and its
// session abandoned, and the bytes received so far are dropped (not the ones of
// this page, which arrive after the mark).
static void ws_adopt_client(int fd) {
    const int previous = ws_client_fd;
    ws_stat_handshakes++;
    ws_client_fd = fd;
    if (previous >= 0 && previous != fd)
        httpd_sess_trigger_close(web_server, previous);
    ws_flush_to = ws_rx_head;
    ws_flush_request = true;
    ws_client_lost = true;
}

static esp_err_t ws_handler(httpd_req_t *req) {
    const int fd = httpd_req_to_sockfd(req);

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
    if (ret == ESP_OK && fd != ws_client_fd)
        ws_adopt_client(fd);
    if (ret == ESP_OK && fd == ws_client_fd) {
        for (size_t i = 0; i < ws_pkt.len; i++) {
            uint16_t next_head = (ws_rx_head + 1) % WS_RX_BUFFER_SIZE;
            if (next_head != ws_rx_tail) { // if there is room / se c'è spazio
                ws_rx_buffer[ws_rx_head] = buf[i];
                ws_rx_head = next_head;
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
    if (sockfd == ws_client_fd) {
        ws_client_fd = -1;
        ws_client_lost = true;
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

    // The App protocol keeps using the serial link until an App docks on the
    // WebSocket (OF_Serial::SerialProcessingWebDock selects the link).

    // 1. Access point. If the ESP-NOW link to the dongle already uses the radio,
    //    keep its channel: a single radio cannot serve two channels.
    uint8_t apChannel = WEBAPP_AP_DEFAULT_CHANNEL;
    if (WiFi.getMode() != WIFI_OFF) {
        uint8_t primary = 0;
        wifi_second_chan_t second = WIFI_SECOND_CHAN_NONE;
        if (esp_wifi_get_channel(&primary, &second) == ESP_OK && primary >= 1 && primary <= 13)
            apChannel = primary;
    }

    WiFi.mode(WIFI_AP_STA);
    WiFi.softAP(WEBAPP_AP_SSID, WEBAPP_AP_PASSWORD, apChannel);

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
    }

    // 4. DNS server for the captive portal
    dnsServer.start(53, "*", WiFi.softAPIP());

    // 5. DNS task on Core 0 so it does not interfere with the main loop
    xTaskCreatePinnedToCore(dns_server_task, "dns_task", 2048, NULL, 1, NULL, 0);
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
