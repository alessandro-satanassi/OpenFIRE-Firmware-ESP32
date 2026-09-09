#include <Arduino.h>
#include "OpenFIREweb.h"
#include "web_assets.h"
bool OF_WebConfigModeActive = false;

#if defined(ARDUINO_ARCH_ESP32)

#include <WiFi.h>
#include <esp_http_server.h>
#include <DNSServer.h>
#include <ESPmDNS.h>

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


static httpd_handle_t web_server = NULL;
static DNSServer dnsServer;


// =================================================================================================
// --- GRUPPO: MOTORE SERIALE HARDWARE ---
static int hw_available() { return Serial.available(); }
static int hw_read() { return Serial.read(); }
static size_t hw_writeByte(uint8_t c) { return Serial.write(c); }
static size_t hw_writeBuf(const uint8_t* buf, size_t size) { return Serial.write(buf, size); }
static void hw_flush() { Serial.flush(); }
static size_t hw_readBytes(char* buf, size_t size) { return Serial.readBytes(buf, size); }

OF_WebSerialWrapper::SerialOps hardwareOps = {
    hw_available, hw_read, hw_writeByte, hw_writeBuf, hw_flush, hw_readBytes
};

// ===========================================================================================================

// --- GRUPPO: MOTORE WEBSOCKET ---

// --- LOGICA RING BUFFER ---
#define WS_RX_BUFFER_SIZE 256
static uint8_t ws_rx_buffer[WS_RX_BUFFER_SIZE];
static uint16_t ws_rx_head = 0;
static uint16_t ws_rx_tail = 0;

static int ws_available() {
    return (WS_RX_BUFFER_SIZE + ws_rx_head - ws_rx_tail) % WS_RX_BUFFER_SIZE;
}

static int ws_read() {
    if (ws_rx_head == ws_rx_tail) return -1;
    uint8_t c = ws_rx_buffer[ws_rx_tail];
    ws_rx_tail = (ws_rx_tail + 1) % WS_RX_BUFFER_SIZE;
    return c;
}

static size_t ws_writeBuf(const uint8_t* buf, size_t size) {
    if (!web_server) return 0;

    httpd_ws_frame_t ws_pkt;
    memset(&ws_pkt, 0, sizeof(httpd_ws_frame_t));
    ws_pkt.type = HTTPD_WS_TYPE_BINARY; // Tutto il protocollo è binario!
    ws_pkt.payload = (uint8_t*)buf;
    ws_pkt.len = size;

    size_t clients = 8;
    int client_fds[8];
    // Trova tutti i client connessi e spara il pacchetto
    if (httpd_get_client_list(web_server, &clients, client_fds) == ESP_OK) {
        for (size_t i = 0; i < clients; ++i) {
            if (httpd_ws_get_fd_info(web_server, client_fds[i]) == HTTPD_WS_CLIENT_WEBSOCKET) {
                httpd_ws_send_data(web_server, client_fds[i], &ws_pkt);
            }
        }
    }
    return size;
}

static size_t ws_writeByte(uint8_t c) {
    return ws_writeBuf(&c, 1);
}

static void ws_flush() {
    // I websocket inviano a pacchetti interi, non c'è un buffer TX da svuotare
}

static size_t ws_readBytes(char* buf, size_t size) {
    size_t count = 0;
    unsigned long startMillis = millis();
    // Timeout di 1000ms, emulando il comportamento standard della Seriale Arduino
    while (count < size && (millis() - startMillis < 1000)) {
        if (ws_available() > 0) {
            buf[count++] = (char)ws_read();
        } else {
            vTaskDelay(pdMS_TO_TICKS(1)); // Fa respirare il FreeRTOS durante l'attesa
        }
    }
    return count;
}

// --- GRUPPO: MOTORE WEBSOCKET ---
OF_WebSerialWrapper::SerialOps websocketOps = {
    ws_available, ws_read, ws_writeByte, ws_writeBuf, ws_flush, ws_readBytes
};

// ===============================================================================================================

// Partiamo di base caricando il motore Hardware!
OF_WebSerialWrapper::SerialOps WebAppSerial::ops = hardwareOps;
// ===============================================================================================================

//static httpd_handle_t web_server = NULL;
//static DNSServer dnsServer;

// Handler per il Captive Portal: intercetta tutte le richieste 404 e redireziona alla pagina principale
static esp_err_t captive_portal_handler(httpd_req_t *req, httpd_err_code_t error) {
    httpd_resp_set_status(req, "302 Found");
    httpd_resp_set_hdr(req, "Location", "http://192.168.4.1/");
    httpd_resp_send(req, NULL, 0);
    return ESP_OK;
}

static esp_err_t index_get_handler(httpd_req_t *req) {
    httpd_resp_set_type(req, "text/html");
    httpd_resp_set_hdr(req, "Content-Encoding", "gzip");
    return httpd_resp_send(req, (const char*)web_index_html_gz, web_index_html_gz_len);
}

static esp_err_t style_get_handler(httpd_req_t *req) {
    httpd_resp_set_type(req, "text/css");
    httpd_resp_set_hdr(req, "Content-Encoding", "gzip");
    return httpd_resp_send(req, (const char*)web_style_css_gz, web_style_css_gz_len);
}

static esp_err_t app_js_get_handler(httpd_req_t *req) {
    httpd_resp_set_type(req, "application/javascript");
    httpd_resp_set_hdr(req, "Content-Encoding", "gzip");
    return httpd_resp_send(req, (const char*)web_app_js_gz, web_app_js_gz_len);
}

static esp_err_t board_svg_get_handler(httpd_req_t *req) {
    httpd_resp_set_type(req, "image/svg+xml");
    httpd_resp_set_hdr(req, "Content-Encoding", "gzip");
    return httpd_resp_send(req, (const char*)web_board_svg_gz, web_board_svg_gz_len);
}

static esp_err_t ws_handler(httpd_req_t *req) {
    // Handshake iniziale quando il browser si connette
    if (req->method == HTTP_GET) return ESP_OK; 
    
    httpd_ws_frame_t ws_pkt;
    memset(&ws_pkt, 0, sizeof(httpd_ws_frame_t));
    ws_pkt.type = HTTPD_WS_TYPE_BINARY; 
    
    // 1. Chiediamo al server QUANTO è lungo il pacchetto arrivato
    esp_err_t ret = httpd_ws_recv_frame(req, &ws_pkt, 0);
    if (ret != ESP_OK || ws_pkt.len == 0) return ret;

    // Se è più grande del nostro RingBuffer lo rifiutiamo
    if (ws_pkt.len > WS_RX_BUFFER_SIZE - 1) return ESP_ERR_NO_MEM;
    
    // 2. Prepariamo un cesto per raccogliere i dati
    uint8_t *buf = (uint8_t*)malloc(ws_pkt.len);
    if (!buf) return ESP_ERR_NO_MEM;
    ws_pkt.payload = buf;
    
    // 3. Scarichiamo i dati nel cesto
    ret = httpd_ws_recv_frame(req, &ws_pkt, ws_pkt.len);
    if (ret == ESP_OK) {
        // 4. Li infiliamo nel nostro RingBuffer per ingannare la pistola!
        for (size_t i = 0; i < ws_pkt.len; i++) {
            uint16_t next_head = (ws_rx_head + 1) % WS_RX_BUFFER_SIZE;
            if (next_head != ws_rx_tail) { // Se c'è spazio
                ws_rx_buffer[ws_rx_head] = buf[i];
                ws_rx_head = next_head;
            }
        }
    }
    free(buf);
    return ret;
}

static const httpd_uri_t uri_index = { .uri = "/", .method = HTTP_GET, .handler = index_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_style = { .uri = "/style.css", .method = HTTP_GET, .handler = style_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_app   = { .uri = "/app.js", .method = HTTP_GET, .handler = app_js_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_board = { .uri = "/board.svg", .method = HTTP_GET, .handler = board_svg_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_ws    = { .uri = "/ws", .method = HTTP_GET, .handler = ws_handler, .user_ctx = NULL, .is_websocket = true };

// Task background per gestire le richieste DNS del Captive Portal
static void dns_server_task(void *pvParameters) {
    while (true) {
        if (OF_WebConfigModeActive) {
            dnsServer.processNextRequest();
        }
        vTaskDelay(pdMS_TO_TICKS(10)); // Pausa di 10ms per non bloccare la CPU
    }
}

void WebApp_Init() {      
    if (!OF_WebConfigModeActive) return;

    // Inizializza il reindirizzamento al nuovo motore
    WebAppSerial::ops = websocketOps;

    // 1. Accendiamo forzatamente il Wi-Fi in modalità ibrida (AP + Station)
    WiFi.mode(WIFI_AP_STA);
    WiFi.softAP("OpenFIRE_Config", "12345678");

    // 2. Avviamo il leggerissimo server HTTP/WS nativo dell'ESP-IDF
    httpd_config_t config = HTTPD_DEFAULT_CONFIG();
    config.max_uri_handlers = 8;

    if (httpd_start(&web_server, &config) == ESP_OK) {
        httpd_register_uri_handler(web_server, &uri_index);
        httpd_register_uri_handler(web_server, &uri_style);
        httpd_register_uri_handler(web_server, &uri_app);
        httpd_register_uri_handler(web_server, &uri_board);
        httpd_register_uri_handler(web_server, &uri_ws);
        
        // Registriamo il captive portal per intercettare tutto ciò che non esiste
        httpd_register_err_handler(web_server, HTTPD_404_NOT_FOUND, captive_portal_handler);
    }
    
    // 3. Avviamo mDNS per permettere l'accesso tramite http://openfire.local
    if (MDNS.begin("openfire")) {
        MDNS.addService("http", "tcp", 80);
    }
    
    // 4. Avviamo il DNS Server per il Captive Portal
    dnsServer.start(53, "*", WiFi.softAPIP());

    // 5. Creiamo il task background per il server DNS (Core 0, così non interferisce col loop principale)
    xTaskCreatePinnedToCore(dns_server_task, "dns_task", 2048, NULL, 1, NULL, 0);
}

void WebApp_Loop() {
    // Svuotato: non c'è più bisogno di chiamarlo dal loop principale!
}


#else
// RP2040 STUBS - Occupano zero spazio e compilano sempre
void WebApp_Init() {}
void WebApp_Loop() {}
#endif


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
