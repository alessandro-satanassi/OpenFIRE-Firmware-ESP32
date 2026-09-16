#ifdef COMMENTO


#include "OpenFIREusbnet.h"

#if defined(ARDUINO_ARCH_ESP32) && defined(OPENFIRE_USB_NCM)
#include <Arduino.h>
#include <Adafruit_TinyUSB.h>
#include <USB.h>
#include <USBCDC.h>
#include <TinyUSB_Devices.h>
#include <atomic>
#include <esp_event.h>
#include <esp_mac.h>
#include <esp_netif.h>
#include <esp_netif_defaults.h>
#include <mdns.h>
#include <freertos/FreeRTOS.h>
#include <freertos/queue.h>
#include "class/net/net_device.h"
#include "device/dcd.h"
#include "device/usbd_pvt.h"

extern "C" uint8_t tud_network_mac_address[6];
//extern "C" uint8_t tud_network_mac_address[6] = {0, 0, 0, 0, 0, 0};
//////////////uint8_t tud_network_mac_address[6] = {0, 0, 0, 0, 0, 0};


#if !CONFIG_IDF_TARGET_ESP32S3 || ARDUINO_USB_MODE || ARDUINO_USB_ON_BOOT
#error OpenFIRE USB selection requires ESP32-S3 OTG with automatic USB startup disabled.
#endif
#if !CONFIG_TINYUSB_NCM_ENABLED || !CFG_TUD_CDC
#error Enable TinyUSB NCM and CDC drivers in platformio.ini.
#endif
#if !CONFIG_ETH_ENABLED || !CONFIG_MDNS_PREDEF_NETIF_ETH
#error OpenFIRE USB mDNS requires the predefined Ethernet interface in the framework.
#endif

// A separate local subnet; this device is not the PC's Internet gateway.
static constexpr uint8_t USB_NET_TX_SLOTS = 4;
static constexpr uint16_t USB_NET_FRAME_MAX = 1514; // Ethernet header + MTU 1500.
static constexpr uint32_t USB_NET_TX_MAX_AGE_MS = 250;
static constexpr uint32_t USB_NET_MDNS_RETRY_MS = 250;

struct UsbNetPacket {
    uint32_t epoch;
    uint32_t timestamp;
    uint16_t length;
    uint8_t data[USB_NET_FRAME_MAX];
};

static std::atomic<bool> usbNetActive(false);
static std::atomic<bool> usbNetAttached(false);
static std::atomic<bool> usbNetDetached(false);
static std::atomic<bool> usbNetLinkUp(false);
static std::atomic<bool> usbNetPumpPending(false);
static std::atomic<bool> usbNetRxPending(false);
static std::atomic<bool> usbNetMdnsReady(false);
// ESP-IDF's 32-bit atomics also support ISR context with its PSRAM workaround.
// Keep this object in normal static/internal RAM (not in an external-RAM pool).
static std::atomic<uint32_t> usbNetEpoch(0);
static esp_netif_t *usbNetif = nullptr;
static esp_netif_driver_base_t usbNetDriver = {};
static QueueHandle_t usbNetFree = nullptr;
static QueueHandle_t usbNetPending = nullptr;
static UsbNetPacket *usbNetPackets = nullptr;
static TaskHandle_t usbNetTask = nullptr;
static char usbNetMacString[13];
static char usbNetSerialString[17];

bool OpenFIREUsbNetActive() { return usbNetActive.load(); }

void OpenFIREUsbNetMDNSReady() {
    // WebApp_Init() owns the single MDNS.begin(); the USB worker only updates
    // its interface after that initialization has completed successfully.
    if (usbNetActive.load() && usbNetAttached.load()) {
        usbNetMdnsReady.store(true);
        xTaskNotifyGive(usbNetTask);
    }
}

// On ESP32, Adafruit_USBD_CDC is an alias of Arduino's USBCDC. Its begin()
// enables I/O but does not add descriptors to Adafruit's configuration.
// Preserve the original CDC interface order, endpoints and string exactly.
class OpenFIRE_USBD_CDC : public Adafruit_USBD_Interface {
public:
    uint16_t getInterfaceDescriptor(uint8_t, uint8_t *buffer, uint16_t size) override {
        if (!buffer) return TUD_CDC_DESC_LEN;
        if (size < TUD_CDC_DESC_LEN) return 0;
        const uint8_t itf = TinyUSBDevice.allocInterface(2);
        const uint8_t str = TinyUSBDevice.addStringDescriptor("TinyUSB Serial");
        const uint8_t desc[] = {
            TUD_CDC_DESCRIPTOR(itf, str, 0x85, 64, 0x03, 0x84, 64)
        };
        memcpy(buffer, desc, sizeof(desc));
        return sizeof(desc);
    }
};

// No CDC instance or endpoint is created for this stream. Residual diagnostic
// prints/reads in web mode return immediately, without using the UART pins.
class OpenFIRE_EmptyStream : public Stream {
public:
    OpenFIRE_EmptyStream() { setTimeout(0); }
    int available() override { return 0; }
    int read() override { return -1; }
    int peek() override { return -1; }
    void flush() override {}
    size_t readBytes(char*, size_t) override { return 0; }
    size_t readBytes(uint8_t*, size_t) override { return 0; }
    size_t write(uint8_t) override { return 0; }
    size_t write(const uint8_t*, size_t) override { return 0; }
};
static USBCDC *usbSerial = nullptr; // Selected once by setup(), before firmware I/O.

Stream& OpenFIREUsbSerial() {
    static OpenFIRE_EmptyStream empty;
    return usbSerial ? (Stream&)*usbSerial : (Stream&)empty;
}

class OpenFIRE_USBD_NCM : public Adafruit_USBD_Interface {
public:
    uint16_t getInterfaceDescriptor(uint8_t, uint8_t *buffer, uint16_t size) override {
        if (!buffer) return TUD_CDC_NCM_DESC_LEN;
        if (size < TUD_CDC_NCM_DESC_LEN) return 0;
        const uint8_t itf = TinyUSBDevice.allocInterface(2);
        const uint8_t macString = TinyUSBDevice.addStringDescriptor(usbNetMacString);
        const uint8_t notification = TinyUSBDevice.allocEndpoint(TUSB_DIR_IN);
        const uint8_t out = TinyUSBDevice.allocEndpoint(TUSB_DIR_OUT);
        const uint8_t in = TinyUSBDevice.allocEndpoint(TUSB_DIR_IN);
        const uint8_t desc[] = {
            TUD_CDC_NCM_DESCRIPTOR(itf, 0, macString, notification, 64, out, in, 64,
                                   USB_NET_FRAME_MAX, 16, 0)
        };
        memcpy(buffer, desc, sizeof(desc));
        return sizeof(desc);
    }
};
static OpenFIRE_USBD_NCM usbNcm;

// All tud_network_* calls run on the existing TinyUSB task, never on lwIP's
// task. Both producers and the USB consumer exchange owned, bounded buffers.
static void usb_net_pump(void *) {
    if (!usbNetAttached.load() || !tud_mounted()) {
        usbNetRxPending.store(false); // The previous USB connection is gone.
    } else if (usbNetLinkUp.load() && usbNetRxPending.exchange(false)) {
        // A first DHCP packet may beat the worker's mount poll. TinyUSB kept
        // that packet; retry it now that the IP interface can also transmit.
        tud_network_recv_renew();
    }
    UsbNetPacket *packet;
    // Bound each USB-task job even if lwIP refills the queue continuously.
    // Other USB events (especially HID/control transfers) must get their turn.
    for (uint8_t sent = 0; sent < USB_NET_TX_SLOTS &&
         xQueuePeek(usbNetPending, &packet, 0) == pdTRUE; ++sent) {
        const bool stale = !usbNetAttached.load() || !usbNetLinkUp.load() ||
                           packet->epoch != usbNetEpoch.load() ||
                           (uint32_t)(millis() - packet->timestamp) > USB_NET_TX_MAX_AGE_MS;
        if (!stale) {
            if (!tud_ready() || !tud_network_can_xmit(packet->length)) break;
            tud_network_xmit(packet->data, packet->length);
        }
        // TinyUSB copies through tud_network_xmit_cb synchronously.
        xQueueReceive(usbNetPending, &packet, 0);
        xQueueSend(usbNetFree, &packet, 0);
    }
    usbNetPumpPending.store(false);
}

static void usb_net_worker(void *) {
    // No TinyUSB calls before USB.begin(), and no access to resources during
    // a failed boot preparation. The creator releases this gate on success.
    ulTaskNotifyTake(pdTRUE, portMAX_DELAY);
    bool wasUp = false;
    bool mdnsPending = true;
    uint32_t mdnsLastAttempt = millis() - USB_NET_MDNS_RETRY_MS;
    uint32_t lastEpoch = usbNetEpoch.load();
    for (;;) {
        const uint32_t epoch = usbNetEpoch.load();
        const bool up = usbNetAttached.load() && tud_mounted();
        if (wasUp && (!up || epoch != lastEpoch)) {
            usbNetLinkUp.store(false);
            esp_netif_action_disconnected(usbNetif, nullptr, 0, nullptr);
            wasUp = false;
            mdnsPending = true;
        }
        if (up && !wasUp) {
            esp_netif_action_connected(usbNetif, nullptr, 0, nullptr);
            usbNetLinkUp.store(true);
            wasUp = true;
            mdnsPending = true;
        }
        lastEpoch = epoch;
        if ((usbNetRxPending.load() || uxQueueMessagesWaiting(usbNetPending)) &&
            !usbNetPumpPending.exchange(true))
            usbd_defer_func(usb_net_pump, nullptr, false);
        // USB has no Ethernet-driver events: update only its mDNS interface
        // on startup/link changes. ENABLE also reprobes after a USB bus reset.
        // The API queues work to the existing mDNS task; no second server/task.
        if (mdnsPending && usbNetMdnsReady.load() &&
            (uint32_t)(millis() - mdnsLastAttempt) >= USB_NET_MDNS_RETRY_MS) {
            mdnsLastAttempt = millis();
            const mdns_event_actions_t action = wasUp ? MDNS_EVENT_ENABLE_IP4 : MDNS_EVENT_DISABLE_IP4;
            if (mdns_netif_action(usbNetif, action) == ESP_OK)
                mdnsPending = false;
        }
        // Retry a full USB transmitter without spinning or blocking lwIP/HID.
        const TickType_t retry = pdMS_TO_TICKS(up ? 2 : 20);
        ulTaskNotifyTake(pdTRUE, retry ? retry : 1);
    }
}

static esp_err_t usb_net_transmit(void *, void *buffer, size_t length) {
    const uint32_t epoch = usbNetEpoch.load();
    if (!usbNetLinkUp.load()) return ESP_ERR_INVALID_STATE;
    if (!buffer || length < 14 || length > USB_NET_FRAME_MAX) return ESP_ERR_INVALID_SIZE;
    UsbNetPacket *packet;
    if (xQueueReceive(usbNetFree, &packet, 0) != pdTRUE) return ESP_ERR_NO_MEM;
    packet->epoch = epoch;
    packet->timestamp = millis();
    packet->length = (uint16_t)length;
    memcpy(packet->data, buffer, length);
    if (xQueueSend(usbNetPending, &packet, 0) != pdTRUE) {
        xQueueSend(usbNetFree, &packet, 0);
        return ESP_ERR_NO_MEM;
    }
    xTaskNotifyGive(usbNetTask);
    return ESP_OK;
}

static void usb_net_free_rx(void *, void *buffer) { free(buffer); }

static esp_err_t usb_net_post_attach(esp_netif_t *netif, esp_netif_iodriver_handle handle) {
    ((esp_netif_driver_base_t*)handle)->netif = netif;
    return ESP_OK;
}

extern "C" bool tud_network_recv_cb(const uint8_t *src, uint16_t size) {
    if (src && size >= 14 && size <= USB_NET_FRAME_MAX &&
        usbNetAttached.load() && !usbNetLinkUp.load()) {
        // Keep ownership in TinyUSB until the worker raises the link. Do not
        // discard the first DHCP request or block the USB task waiting for it.
        usbNetRxPending.store(true);
        xTaskNotifyGive(usbNetTask);
        return false;
    }
    if (usbNetAttached.load() && usbNetLinkUp.load() &&
        src && size >= 14 && size <= USB_NET_FRAME_MAX) {
        // esp_netif/lwIP owns this copy until driver_free_rx_buffer is called.
        void *copy = malloc(size);
        if (copy) {
            memcpy(copy, src, size);
            // Ethernet input also releases the copy on failure (including a
            // full lwIP queue). Freeing it here again would be a double free.
            esp_netif_receive(usbNetif, copy, size, copy);
        }
    }
    // Also renew after drops: one bad/full receive must not stop the USB link.
    tud_network_recv_renew();
    return true;
}

extern "C" uint16_t tud_network_xmit_cb(uint8_t *dst, void *ref, uint16_t length) {
    memcpy(dst, ref, length);
    return length;
}

extern "C" void tud_network_init_cb(void) {}

// No changes to Arduino's mount/unmount callbacks. A bus reset invalidates
// queued Ethernet packets even when unplug/replug is faster than the worker.
extern "C" void tud_event_hook_cb(uint8_t, uint32_t event, bool) {
    if (usbNetActive.load() && (event == DCD_EVENT_BUS_RESET || event == DCD_EVENT_UNPLUGGED))
        usbNetEpoch.fetch_add(1);
}

// Used only before an interface/task has been published to USB/lwIP.
static void usb_net_cleanup(void) {
    if (usbNetif) { esp_netif_destroy(usbNetif); usbNetif = nullptr; }
    if (usbNetFree) { vQueueDelete(usbNetFree); usbNetFree = nullptr; }
    if (usbNetPending) { vQueueDelete(usbNetPending); usbNetPending = nullptr; }
    free(usbNetPackets);
    usbNetPackets = nullptr;
    usbNetDriver.netif = nullptr;
}

static bool usb_net_prepare() {
    uint8_t hostMac[6], deviceMac[6], serialMac[6];
    if (esp_read_mac(hostMac, ESP_MAC_ETH) != ESP_OK ||
        esp_read_mac(serialMac, ESP_MAC_WIFI_STA) != ESP_OK) return false;
    hostMac[0] = (hostMac[0] | 2) & 0xFE; // Locally administered, unicast.
    memcpy(deviceMac, hostMac, sizeof(deviceMac));
    deviceMac[5] ^= 1; // Host descriptor and MCU netif must NOT share one MAC.

    // --------------------------------
    memcpy(tud_network_mac_address, hostMac, 6);
    // --------------------------------

    snprintf(usbNetMacString, sizeof(usbNetMacString), "%02X%02X%02X%02X%02X%02X",
             hostMac[0], hostMac[1], hostMac[2], hostMac[3], hostMac[4], hostMac[5]);

    
    snprintf(usbNetSerialString, sizeof(usbNetSerialString), "%02X%02X%02X%02X%02X%02X-NCM",
             serialMac[0], serialMac[1], serialMac[2], serialMac[3], serialMac[4], serialMac[5]);
    

    /*
    snprintf(usbNetSerialString, sizeof(usbNetSerialString), "%02X%02X%02X%02X%02X%02X-NCM2",
             serialMac[0], serialMac[1], serialMac[2], serialMac[3], serialMac[4], serialMac[5]);
    */


    if (esp_netif_init() != ESP_OK) return false;
    const esp_err_t eventResult = esp_event_loop_create_default();
    if (eventResult != ESP_OK && eventResult != ESP_ERR_INVALID_STATE) return false;
    usbNetPackets = (UsbNetPacket*)calloc(USB_NET_TX_SLOTS, sizeof(UsbNetPacket));
    usbNetFree = xQueueCreate(USB_NET_TX_SLOTS, sizeof(UsbNetPacket*));
    usbNetPending = xQueueCreate(USB_NET_TX_SLOTS, sizeof(UsbNetPacket*));
    if (!usbNetPackets || !usbNetFree || !usbNetPending) {
        usb_net_cleanup();
        return false;
    }
    for (uint8_t i = 0; i < USB_NET_TX_SLOTS; ++i) {
        UsbNetPacket *packet = &usbNetPackets[i];
        xQueueSend(usbNetFree, &packet, 0);
    }

    esp_netif_ip_info_t ip = {};
    esp_netif_set_ip4_addr(&ip.ip, 192, 168, 7, 1);
    esp_netif_set_ip4_addr(&ip.netmask, 255, 255, 255, 0);
    esp_netif_inherent_config_t base = {};
    base.flags = (esp_netif_flags_t)(ESP_NETIF_DHCP_SERVER | ESP_NETIF_FLAG_AUTOUP);
    base.ip_info = &ip;
    // Reuse mDNS's existing Ethernet slot: this framework has no spare slot
    // for mdns_register_netif(). NCM is the gun's only Ethernet interface.
    // This internal key changes neither the USB identity nor DHCP/routing.
    base.if_key = "ETH_DEF";
    base.if_desc = "OpenFIRE USB";
    base.route_prio = 1; // Wi-Fi retains priority for non-local routes.
    usbNetDriver.post_attach = usb_net_post_attach;
    esp_netif_driver_ifconfig_t driver = {};
    driver.handle = &usbNetDriver;
    driver.transmit = usb_net_transmit;
    driver.driver_free_rx_buffer = usb_net_free_rx;
    esp_netif_config_t config = {};
    config.base = &base;
    config.driver = &driver;
    config.stack = ESP_NETIF_NETSTACK_DEFAULT_ETH;
    usbNetif = esp_netif_new(&config);
    uint8_t noOffer = 0;
    if (!usbNetif || esp_netif_attach(usbNetif, &usbNetDriver) != ESP_OK ||
        esp_netif_set_mac(usbNetif, deviceMac) != ESP_OK ||
        esp_netif_dhcps_option(usbNetif, ESP_NETIF_OP_SET, ESP_NETIF_ROUTER_SOLICITATION_ADDRESS,
                               &noOffer, sizeof(noOffer)) != ESP_OK ||
        esp_netif_dhcps_option(usbNetif, ESP_NETIF_OP_SET, ESP_NETIF_DOMAIN_NAME_SERVER,
                               &noOffer, sizeof(noOffer)) != ESP_OK) {
        usb_net_cleanup();
        return false;
    }
    // Allocate/validate everything before adding a descriptor. On failure the
    // unchanged configuration can still receive CDC and the existing HID.
    esp_netif_action_start(usbNetif, nullptr, 0, nullptr);
    esp_netif_dhcp_status_t dhcpStatus;
    if (esp_netif_dhcps_get_status(usbNetif, &dhcpStatus) != ESP_OK ||
        dhcpStatus != ESP_NETIF_DHCP_STARTED) {
        esp_netif_action_stop(usbNetif, nullptr, 0, nullptr);
        usb_net_cleanup();
        return false;
    }
    if (xTaskCreate(usb_net_worker, "usb_net", 3072, nullptr, 2, &usbNetTask) != pdPASS) {
        esp_netif_action_stop(usbNetif, nullptr, 0, nullptr);
        usb_net_cleanup();
        return false;
    }
    if (!TinyUSBDevice.addInterface(usbNcm)) {
        // The worker is still behind its startup gate; no callback can own
        // these buffers yet. The descriptor rejects lack of space before
        // allocating interface/endpoint numbers.
        vTaskDelete(usbNetTask);
        usbNetTask = nullptr;
        esp_netif_action_stop(usbNetif, nullptr, 0, nullptr);
        usb_net_cleanup();
        return false;
    }
    // Distinguish this instance from the normal COM port in the host cache.
    TinyUSBDevice.setSerialDescriptor(usbNetSerialString);
    usbNetActive.store(true);
    return true;
}

bool OpenFIREUsbBegin(bool webConfig, uint8_t pollRate) {
    // Deliberately not a runtime mode switch. Never edit live descriptors.
    static bool attempted = false;
    static bool started = false;
    if (attempted) return started;
    attempted = true;
    if (TinyUSBDevice.isInitialized()) {
        log_e("OpenFIRE USB was started before boot mode selection");
        return false;
    }

    if (webConfig && !usb_net_prepare())
        log_e("OpenFIRE NCM preparation failed; keeping CDC and Wi-Fi configuration");

    if (!OpenFIREUsbNetActive()) {
        static OpenFIRE_USBD_CDC descriptor;
        static USBCDC serial(0); // Created only on the CDC boot path.
        if (!TinyUSBDevice.addInterface(descriptor)) {
            log_e("OpenFIRE CDC descriptor does not fit");
            return false;
        }
        serial.begin(9600);
        serial.setTimeout(0);
        usbSerial = &serial;
    }

    // The existing single HID uses IN 0x81 with CDC, or IN 0x83 with NCM.
    // Nothing is mounted yet: begin() cannot trigger its detach/attach path.
    TinyUSBDevices.begin(pollRate);
    started = USB.begin();
    if (started && OpenFIREUsbNetActive()) {
        usbNetAttached.store(true);
        xTaskNotifyGive(usbNetTask);
    }
    if (!started) log_e("OpenFIRE USB startup failed");
    return started;
}

/*
bool OpenFIREUsbBegin(bool webConfig, uint8_t pollRate) {
    // Deliberately not a runtime mode switch. Never edit live descriptors.
    static bool attempted = false;
    static bool started = false;
    if (attempted) return started;
    attempted = true;
    if (TinyUSBDevice.isInitialized()) {
        log_e("OpenFIRE USB was started before boot mode selection");
        return false;
    }

    // FONDAMENTALE: Inizializza il device Adafruit per impostare la bDeviceClass a 0xEF (IAD)
    // e resettare correttamente i contatori delle stringhe, PRIMA di aggiungere le interfacce!
    TinyUSBDevice.begin(0);

    if (webConfig && !usb_net_prepare())
        log_e("OpenFIRE NCM preparation failed; keeping CDC and Wi-Fi configuration");

    if (!OpenFIREUsbNetActive()) {
        static OpenFIRE_USBD_CDC descriptor;
        static USBCDC serial(0); // Created only on the CDC boot path.
        if (!TinyUSBDevice.addInterface(descriptor)) {
            log_e("OpenFIRE CDC descriptor does not fit");
            return false;
        }
        serial.begin(9600);
        serial.setTimeout(0);
        usbSerial = &serial;
    }

    // The existing single HID uses IN 0x81 with CDC, or IN 0x83 with NCM.
    // Nothing is mounted yet: begin() cannot trigger its detach/attach path.
    TinyUSBDevices.begin(pollRate);
    started = USB.begin();
    if (started && OpenFIREUsbNetActive()) {
        usbNetAttached.store(true);
        xTaskNotifyGive(usbNetTask);
    }
    if (!started) log_e("OpenFIRE USB startup failed");
    return started;
}
*/

static void usb_net_detach(void *) {
    TinyUSBDevice.detach();
    usbNetDetached.store(true);
}

void OpenFIREUsbDetach() {
    usbNetAttached.store(false);
    if (!TinyUSBDevice.isInitialized()) return;
    if (!OpenFIREUsbNetActive()) {
        TinyUSBDevice.detach();
        return;
    }
    // Main-task only. Finish the current USB job before returning: in
    // particular the ROM-bootloader helper may next disable the USB PHY.
    // Later network jobs only discard packets because attached is false.
    if (!usbNetDetached.load()) {
        usbd_defer_func(usb_net_detach, nullptr, false);
        while (!usbNetDetached.load()) vTaskDelay(1);
    }
    // Descriptors/buffers stay alive for pending callbacks until reboot.
}

/*
#if defined(ARDUINO_ARCH_ESP32) && defined(OPENFIRE_USB_NCM)

extern "C" {
    #include "class/net/net_device.h"

    // Questa è la struttura del driver NCM compilata da Adafruit
    static const usbd_class_driver_t ncm_app_driver = {
        #if CFG_TUSB_DEBUG >= 2
        .name = "NCM",
        #endif
        .init = netd_init,
        .reset = netd_reset,
        .open = netd_open,
        .control_xfer_cb = netd_control_xfer_cb,
        .xfer_cb = netd_xfer_cb,
        .sof = NULL
    };

    // Hook chiamato automaticamente da TinyUSB quando viene fatto USB.begin()
    usbd_class_driver_t const* usbd_app_driver_get_cb(uint8_t* driver_count) {
        // Se la rete virtuale NCM è attiva, inietta il driver
        if (OpenFIREUsbNetActive()) {
            *driver_count = 1;
            return &ncm_app_driver;
        }
        
        *driver_count = 0;
        return NULL;
    }
}

#endif
*/
/*
extern "C" {
    void netd_init(void);
    bool netd_deinit(void);
    void netd_reset(uint8_t rhport);
    uint16_t netd_open(uint8_t rhport, tusb_desc_interface_t const * itf_desc, uint16_t max_len);
    bool netd_control_xfer_cb(uint8_t rhport, uint8_t stage, tusb_control_request_t const * request);
    bool netd_xfer_cb(uint8_t rhport, uint8_t ep_addr, xfer_result_t result, uint32_t xferred_bytes);

    // IL NOSTRO INTERCETTATORE MAGICO
    bool custom_netd_control_xfer_cb(uint8_t rhport, uint8_t stage, tusb_control_request_t const * request) {
        // Intercettiamo la richiesta NCM_GET_NET_ADDRESS (bmRequestType = 0xA1, bRequest = 0x81)
        if (request->bmRequestType == 0xA1 && request->bRequest == 0x81) {
            if (stage == CONTROL_STAGE_SETUP) {
                // Forniamo a Windows il MAC address corretto
                return tud_control_xfer(rhport, request, (void*)tud_network_mac_address, 6);
            }
            return true;
        }
        // Per tutte le altre richieste, passiamo la palla al driver di Adafruit
        return netd_control_xfer_cb(rhport, stage, request);
    }

    static usbd_class_driver_t const ncm_app_driver = {
#if CFG_TUSB_DEBUG >= 2
        .name = "NCM",
#endif
        .init = netd_init,
        .deinit = netd_deinit,
        .reset = netd_reset,
        .open = netd_open,
        .control_xfer_cb = custom_netd_control_xfer_cb, // USIAMO IL NOSTRO INTERCETTATORE!
        .xfer_cb = netd_xfer_cb,
        .sof = NULL
    };

    usbd_class_driver_t const* usbd_app_driver_get_cb(uint8_t* driver_count) {
        *driver_count = 1;
        return &ncm_app_driver;
    }
}
*/

/*
extern "C" {
    void netd_init(void);
    bool netd_deinit(void);
    void netd_reset(uint8_t rhport);
    uint16_t netd_open(uint8_t rhport, tusb_desc_interface_t const * itf_desc, uint16_t max_len);
    bool netd_control_xfer_cb(uint8_t rhport, uint8_t stage, tusb_control_request_t const * request);
    bool netd_xfer_cb(uint8_t rhport, uint8_t ep_addr, xfer_result_t result, uint32_t xferred_bytes);

    // FUNZIONI DUMMY: Salvano il sistema dalla corruzione della memoria (evitano la doppia esecuzione)
    void dummy_netd_init(void) {}
    bool dummy_netd_deinit(void) { return true; }
    void dummy_netd_reset(uint8_t rhport) {}

    // IL NOSTRO INTERCETTATORE MAGICO
    bool custom_netd_control_xfer_cb(uint8_t rhport, uint8_t stage, tusb_control_request_t const * request) {
        // Intercettiamo la richiesta NCM_GET_NET_ADDRESS
        if (request->bmRequestType == 0xA1 && request->bRequest == 0x81) {
            if (stage == CONTROL_STAGE_SETUP) {
                return tud_control_xfer(rhport, request, (void*)tud_network_mac_address, 6);
            }
            return true;
        }
        // Per tutte le altre richieste, passiamo la palla al driver
        return netd_control_xfer_cb(rhport, stage, request);
    }

    static usbd_class_driver_t const ncm_app_driver = {
#if CFG_TUSB_DEBUG >= 2
        .name = "NCM",
#endif
        .init = dummy_netd_init,           // <--- FONDAMENTALE! Evita il bug della doppia init!
        .deinit = dummy_netd_deinit,
        .reset = dummy_netd_reset,
        .open = netd_open,
        .control_xfer_cb = custom_netd_control_xfer_cb,
        .xfer_cb = netd_xfer_cb,
        .sof = NULL
    };

    usbd_class_driver_t const* usbd_app_driver_get_cb(uint8_t* driver_count) {
        // Iniettiamo l'hook SOLO se NCM è attivo. 
        // Questo ripara definitivamente la modalità Normale (CDC) che tornerà a funzionare!
        if (OpenFIREUsbNetActive()) {
            *driver_count = 1;
            return &ncm_app_driver;
        }
        *driver_count = 0;
        return NULL;
    }
}
*/

#ifdef COMMENTO


//volatile bool OF_NcmDataInterfaceReady = false;

extern "C" {
    void netd_init(void);
    bool netd_deinit(void);
    void netd_reset(uint8_t rhport);
    uint16_t netd_open(uint8_t rhport, tusb_desc_interface_t const * itf_desc, uint16_t max_len);
    bool netd_control_xfer_cb(uint8_t rhport, uint8_t stage, tusb_control_request_t const * request);
    bool netd_xfer_cb(uint8_t rhport, uint8_t ep_addr, xfer_result_t result, uint32_t xferred_bytes);

    void dummy_netd_init(void) {}
    bool dummy_netd_deinit(void) { return true; }
    void dummy_netd_reset(uint8_t rhport) {}

    // BUFFER ALLINEATI E FORZATI IN RAM PER IL DMA DELL'ESP32-S3
    alignas(4) uint8_t tud_network_mac_address[6] = {0}; 
    
    /*
    // Parametri NTB (con buffer impostato a 8192 byte) posizionati in RAM
    alignas(4) static uint8_t ntb_parameters_ram[28] = {
        0x1C, 0x00,             // wLength = 28
        0x01, 0x00,             // bmNtbFormatsSupported = 1
        0x00, 0x20, 0x00, 0x00, // dwNtbInMaxSize = 8192
        0x04, 0x00,             // wNdbInDivisor = 4
        0x00, 0x00,             // wNdbInPayloadRemainder = 0
        0x04, 0x00,             // wNdbInAlignment = 4
        0x00, 0x00,             // wReserved = 0
        0x00, 0x20, 0x00, 0x00, // dwNtbOutMaxSize = 8192
        0x04, 0x00,             // wNdbOutDivisor = 4
        0x00, 0x00,             // wNdbOutPayloadRemainder = 0
        0x04, 0x00,             // wNdbOutAlignment = 4
        0x06, 0x00              // wNtbOutMaxDatagrams = 6
    };
    */

    /*
    // Parametri NTB in RAM
    alignas(4) static uint8_t ntb_parameters_ram[28] = {
        0x1C, 0x00,             // wLength
        0x01, 0x00,             // bmNtbFormatsSupported
        0x80, 0x0C, 0x00, 0x00, // dwNtbInMaxSize = 3200 (0x0C80)
        0x04, 0x00,             // wNdpInDivisor
        0x00, 0x00,             // wNdpInPayloadRemainder
        0x04, 0x00,             // wNdpInAlignment
        0x00, 0x00,             // wReserved
        0x80, 0x0C, 0x00, 0x00, // dwNtbOutMaxSize = 3200 (0x0C80)
        0x04, 0x00,             // wNdpOutDivisor
        0x00, 0x00,             // wNdpOutPayloadRemainder
        0x04, 0x00,             // wNdpOutAlignment
        0x06, 0x00              // wNtbOutMaxDatagrams = 6  <--- CAMBIA QUI DA 0x00 A 0x06
        //0x00, 0x00              // wNtbOutMaxDatagrams
    };
    */

    // Aggiungi questa in cima al blocco extern "C"
    volatile bool OF_NcmDataInterfaceReady = false;
    
    /*
    // INTERCETTATORE AVANZATO
    bool custom_netd_control_xfer_cb(uint8_t rhport, uint8_t stage, tusb_control_request_t const * request) {
        if (request->bmRequestType == 0xA1) {
            // Intercetta la richiesta dei parametri NTB
            if (request->bRequest == 0x80) {
                if (stage == CONTROL_STAGE_SETUP) {
                    return tud_control_xfer(rhport, request, ntb_parameters_ram, 28);
                }
                return true;
            }
            // Intercetta la richiesta del MAC Address
            if (request->bRequest == 0x81) {
                if (stage == CONTROL_STAGE_SETUP) {
                    return tud_control_xfer(rhport, request, tud_network_mac_address, 6);
                }
                return true;
            }
        }

        // ---> HAI DIMENTICATO QUESTO PEZZO QUI SOTTO! <---
        // Intercetta il momento esatto in cui Windows apre la scheda di rete
        if (request->bmRequestType == 0x01 && request->bRequest == 0x0B) { // SET_INTERFACE
            if (stage == CONTROL_STAGE_ACK && request->wValue == 1) {
                OF_NcmDataInterfaceReady = true; // WINDOWS È PRONTO, SBLOCCA IL SETUP!
            }
        }
 
        // Per tutto il resto, procedi normalmente
        return netd_control_xfer_cb(rhport, stage, request);
    }
    */

        bool custom_netd_control_xfer_cb(uint8_t rhport, uint8_t stage, tusb_control_request_t const * request) {
        // Intercetta il momento esatto in cui Windows apre la scheda di rete
        if (request->bmRequestType == 0x01 && request->bRequest == 0x0B) { // SET_INTERFACE
            if (stage == CONTROL_STAGE_ACK && request->wValue == 1) {
                OF_NcmDataInterfaceReady = true; // WINDOWS È PRONTO, SBLOCCA IL SETUP!
            }
        }
        
        // Passa SEMPRE il controllo al driver originale di TinyUSB per fare il suo lavoro
        return netd_control_xfer_cb(rhport, stage, request);
    }

    static usbd_class_driver_t const ncm_app_driver = {
#if CFG_TUSB_DEBUG >= 2
        .name = "NCM",
#endif
        .init = dummy_netd_init,
        .deinit = dummy_netd_deinit,
        .reset = dummy_netd_reset,
        .open = netd_open,
        .control_xfer_cb = custom_netd_control_xfer_cb,
        .xfer_cb = netd_xfer_cb,
        .sof = NULL
    };

    usbd_class_driver_t const* usbd_app_driver_get_cb(uint8_t* driver_count) {
        if (OpenFIREUsbNetActive()) {
            *driver_count = 1;
            return &ncm_app_driver;
        }
        *driver_count = 0;
        return NULL;
    }
}

#endif // COMMENTO

extern "C" {
    void netd_init(void);
    bool netd_deinit(void);
    void netd_reset(uint8_t rhport);
    uint16_t netd_open(uint8_t rhport, tusb_desc_interface_t const * itf_desc, uint16_t max_len);
    bool netd_control_xfer_cb(uint8_t rhport, uint8_t stage, tusb_control_request_t const * request);
    bool netd_xfer_cb(uint8_t rhport, uint8_t ep_addr, xfer_result_t result, uint32_t xferred_bytes);

    void dummy_netd_init(void) {}
    bool dummy_netd_deinit(void) { return true; }
    void dummy_netd_reset(uint8_t rhport) {}

    // MAC Address in RAM
    alignas(4) uint8_t tud_network_mac_address[6] = {0}; 
    
    // Parametri NTB in RAM (3200 byte e Datagrams=6)
    alignas(4) static uint8_t ntb_parameters_ram[28] = {
        0x1C, 0x00,             // wLength
        0x01, 0x00,             // bmNtbFormatsSupported
        0x80, 0x0C, 0x00, 0x00, // dwNtbInMaxSize = 3200 (0x0C80)
        0x04, 0x00,             // wNdpInDivisor
        0x00, 0x00,             // wNdpInPayloadRemainder
        0x04, 0x00,             // wNdpInAlignment
        0x00, 0x00,             // wReserved
        0x80, 0x0C, 0x00, 0x00, // dwNtbOutMaxSize = 3200 (0x0C80)
        0x04, 0x00,             // wNdpOutDivisor
        0x00, 0x00,             // wNdpOutPayloadRemainder
        0x04, 0x00,             // wNdpOutAlignment
        0x06, 0x00              // wNtbOutMaxDatagrams = 6  <--- CORRETTO!
    };

    volatile bool OF_NcmDataInterfaceReady = false;
    
    // INTERCETTATORE AVANZATO (Ripristinato)
    bool custom_netd_control_xfer_cb(uint8_t rhport, uint8_t stage, tusb_control_request_t const * request) {
        if (request->bmRequestType == 0xA1) {
            // Intercetta la richiesta dei parametri NTB
            if (request->bRequest == 0x80) {
                if (stage == CONTROL_STAGE_SETUP) {
                    return tud_control_xfer(rhport, request, ntb_parameters_ram, 28);
                }
                return true;
            }
            // Intercetta la richiesta del MAC Address
            if (request->bRequest == 0x81) {
                if (stage == CONTROL_STAGE_SETUP) {
                    return tud_control_xfer(rhport, request, tud_network_mac_address, 6);
                }
                return true;
            }
        }

        // Intercetta il momento esatto in cui Windows apre la scheda di rete
        if (request->bmRequestType == 0x01 && request->bRequest == 0x0B) { // SET_INTERFACE
            if (stage == CONTROL_STAGE_ACK && request->wValue == 1) {
                OF_NcmDataInterfaceReady = true; // WINDOWS È PRONTO, SBLOCCA IL SETUP!
            }
        }
 
        // Per tutto il resto, procedi normalmente
        return netd_control_xfer_cb(rhport, stage, request);
    }

    static usbd_class_driver_t const ncm_app_driver = {
#if CFG_TUSB_DEBUG >= 2
        .name = "NCM",
#endif
        .init = dummy_netd_init,
        .deinit = dummy_netd_deinit,
        .reset = dummy_netd_reset,
        .open = netd_open,
        .control_xfer_cb = custom_netd_control_xfer_cb,
        .xfer_cb = netd_xfer_cb,
        .sof = NULL
    };

    usbd_class_driver_t const* usbd_app_driver_get_cb(uint8_t* driver_count) {
        if (OpenFIREUsbNetActive()) {
            *driver_count = 1;
            return &ncm_app_driver;
        }
        *driver_count = 0;
        return NULL;
    }
}

#endif

#endif // COMMENTO