import os
import gzip
import datetime

Import("env")

# Obtain the current PlatformIO environment name (e.g., WAVESHARE_ESP32_S3_ZERO_N8R8)
env_name = env["PIOENV"]
print(f"\n[WebApp Packer] Building Web Assets for environment: {env_name}")

# Directories
PROJECT_DIR = env.get("PROJECT_DIR")
WEBAPP_DIR = os.path.join(PROJECT_DIR, "webapp")
BOARDS_PICS_DIR = os.path.join(PROJECT_DIR, "..", "shared_boards", "boardPics") # adjust if it's different in the patch
if not os.path.exists(BOARDS_PICS_DIR):
    BOARDS_PICS_DIR = os.path.join(PROJECT_DIR, "src", "boards", "boardPics")

INCLUDE_DIR = os.path.join(PROJECT_DIR, "include")
OUTPUT_FILE = os.path.join(INCLUDE_DIR, "web_assets.h")

# Map environment names to SVG filenames
SVG_MAP = {
    "ESP32_S3_WROOM1_DevKitC_1_N16R8": "esp32-s3-devkitc-1.svg",
    "ESP32_S3_WROOM1_DevKitC_1_N8R2": "esp32-s3-devkitc-1.svg",
    "WAVESHARE_ESP32_S3_PICO": "waveshare-esp32-s3-pico.svg",
    "WAVESHARE_ESP32_S3_ZERO_N8R8": "waveshare-esp32-s3-zero.svg",
    "WAVESHARE_ESP32_S3_ZERO_N4R2": "waveshare-esp32-s3-zero.svg",
    "rpipico": "rpipico.svg",
    "rpipicow": "rpipicow.svg",
    "rpipico2": "rpipico2.svg",
    "rpipico2w": "rpipico2w.svg",
}

svg_filename = SVG_MAP.get(env_name, "generic.svg")
svg_path = os.path.join(BOARDS_PICS_DIR, svg_filename)

# Files to pack (HTML/JS/CSS)
files_to_pack = [
    {"name": "index.html", "var": "web_index_html"},
    {"name": "style.css", "var": "web_style_css"},
    {"name": "app.js", "var": "web_app_js"}
]

if not os.path.exists(INCLUDE_DIR):
    os.makedirs(INCLUDE_DIR)

def generate_c_array(file_path, var_name, is_gzip=True):
    if not os.path.exists(file_path):
        # Create empty placeholder if file doesn't exist yet
        return f"const uint8_t {var_name}_gz[] PROGMEM = {{0x00}};\nconst size_t {var_name}_gz_len = 1;\n"
    
    with open(file_path, "rb") as f:
        data = f.read()
        
    if is_gzip:
        data = gzip.compress(data)
        
    hex_array = ', '.join([f'0x{b:02X}' for b in data])
    
    out = f"// Source: {os.path.basename(file_path)} (GZIPPED: {is_gzip}, Size: {len(data)} bytes)\n"
    out += f"const uint8_t {var_name}_gz[] PROGMEM = {{{hex_array}}};\n"
    out += f"const size_t {var_name}_gz_len = {len(data)};\n\n"
    return out

print(f"[WebApp Packer] Generating {OUTPUT_FILE} ...")

with open(OUTPUT_FILE, "w") as out_f:
    out_f.write(f"// AUTO-GENERATED FILE. DO NOT EDIT.\n")
    out_f.write(f"// Generated on {datetime.datetime.now()}\n")
    out_f.write(f"// Target Environment: {env_name}\n\n")
    out_f.write("#pragma once\n")
    out_f.write("#include <Arduino.h>\n")
    out_f.write("#include <stdint.h>\n")
    out_f.write("#include <stddef.h>\n\n")
    
    # Pack web files
    for f_info in files_to_pack:
        f_path = os.path.join(WEBAPP_DIR, f_info["name"])
        out_f.write(generate_c_array(f_path, f_info["var"], is_gzip=True))
        
    # Pack board specific SVG
    if os.path.exists(svg_path):
        out_f.write(generate_c_array(svg_path, "web_board_svg", is_gzip=True))
    else:
        print(f"[WebApp Packer] WARNING: SVG {svg_path} not found!")
        out_f.write(f"// WARNING: SVG {svg_filename} not found during build\n")
        out_f.write("const uint8_t web_board_svg_gz[] PROGMEM = {0x00};\nconst size_t web_board_svg_gz_len = 1;\n")

print("[WebApp Packer] Done!\n")
