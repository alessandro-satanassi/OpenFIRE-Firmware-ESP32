#include <Arduino.h>
#include "OpenFIREweb.h"
#include "web_assets.h"


bool OF_WebConfigModeActive = false;

#if defined(ARDUINO_ARCH_ESP32)

#include <WiFi.h>
#include <esp_http_server.h>

static httpd_handle_t web_server = NULL;

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

// WebSocket handler skeleton
static esp_err_t ws_handler(httpd_req_t *req) {
    if (req->method == HTTP_GET) {
        return ESP_OK; // Handshake
    }
    // TODO: Gestione pacchetti binari ricevuti dalla WebApp
    return ESP_OK;
}

static const httpd_uri_t uri_index = { .uri = "/", .method = HTTP_GET, .handler = index_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_style = { .uri = "/style.css", .method = HTTP_GET, .handler = style_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_app   = { .uri = "/app.js", .method = HTTP_GET, .handler = app_js_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_board = { .uri = "/board.svg", .method = HTTP_GET, .handler = board_svg_get_handler, .user_ctx = NULL };
static const httpd_uri_t uri_ws    = { .uri = "/ws", .method = HTTP_GET, .handler = ws_handler, .user_ctx = NULL, .is_websocket = true };

void WebApp_Init() {
    if (!OF_WebConfigModeActive) return;

    // 1. Accendiamo forzatamente il Wi-Fi in modalità ibrida (AP + Station)
    // Non importa se prima era spento o solo in STA, questa riga sovrascrive tutto.
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
    }
}

void WebApp_Loop() {
    if (!OF_WebConfigModeActive) return;
    // Con esp_http_server non c'è bisogno di un loop manuale! 
    // Gira automaticamente su un FreeRTOS task background.
}

void WebApp_SendToApp(const uint8_t* buffer, size_t length) {
    if (!OF_WebConfigModeActive || !web_server) return;

    // TODO: Inviare i pacchetti binari ai client WebSocket connessi
}

#else
// RP2040 STUBS - Occupano zero spazio e compilano sempre
void WebApp_Init() {}
void WebApp_Loop() {}
void WebApp_SendToApp(const uint8_t* buffer, size_t length) {}
#endif
