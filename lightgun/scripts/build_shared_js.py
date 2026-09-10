import os
import re

def generate_shared_js(project_dir, webapp_dir):
    header_path = os.path.join(project_dir, "src", "boards", "OpenFIREshared.h")
    out_js_path = os.path.join(webapp_dir, "boards", "OpenFIREshared.js")
    
    if not os.path.exists(os.path.dirname(out_js_path)):
        os.makedirs(os.path.dirname(out_js_path))
        
    with open(header_path, "r", encoding="utf-8") as f:
        content = f.read()
        
    js_content = "// AUTO-GENERATED from OpenFIREshared.h\n"
    js_content += "const OpenFIREshared = {\n"
    
    env_vars = {}
    env_vars["boardArchs[boardRP]"] = "rp2040_235X"
    env_vars["boardArchs[boardESP32_S3]"] = "esp32-s3"
    
    enum_pattern = re.compile(r'enum\s*\{([^}]+)\}\s*([A-Za-z0-9_]+)\s*;')
    for match in enum_pattern.finditer(content):
        enum_body = match.group(1)
        enum_name = match.group(2)
        js_content += f"  {enum_name}: {{\n"
        enum_body = re.sub(r'//.*', '', enum_body)
        enum_body = re.sub(r'/\*.*?\*/', '', enum_body, flags=re.DOTALL)
        
        items = [i.strip() for i in enum_body.split(',') if i.strip()]
        current_val = 0
        for item in items:
            if '=' in item:
                parts = item.split('=')
                k = parts[0].strip()
                v = parts[1].strip()
                try:
                    v_eval = v.replace('0b', '0b')
                    current_val = eval(v_eval, {}, env_vars)
                except Exception:
                    current_val = v
                js_content += f"    {k}: {current_val},\n"
                env_vars[k] = current_val
                if isinstance(current_val, int): current_val += 1
            else:
                js_content += f"    {item}: {current_val},\n"
                env_vars[item] = current_val
                if isinstance(current_val, int): current_val += 1
        js_content += "  },\n"

    # Fix the map pattern
    map_pattern = re.compile(r'const\s+std::(?:unordered_)?map[\s\S]*?\s+([A-Za-z0-9_]+)\s*=\s*\{([\s\S]*?)\};')
    for match in map_pattern.finditer(content):
        map_name = match.group(1)
        map_body = match.group(2)
        if not map_name.startswith('board') and map_name != 'pinCapabilitiesMap' and map_name != 'mcuCapableMaps':
            continue
        js_content += f"  {map_name}: {{\n"
        map_body = re.sub(r'//.*', '', map_body)
        map_body = re.sub(r'/\*.*?\*/', '', map_body, flags=re.DOTALL)
        
        entry_pattern = re.compile(r'\{\s*(?:"([^"]+)"|([A-Za-z0-9_\[\]]+))\s*,\s*([^}]+?)\s*\}')
        for e_match in entry_pattern.finditer(map_body):
            board_name = e_match.group(1) or e_match.group(2)
            if board_name in env_vars:
                board_name = env_vars[board_name]
            val_block = e_match.group(3).strip()
            if val_block.startswith('{'):
                val_block = val_block[1:].strip()
                vals_clean = []
                for v in val_block.split(','):
                    v = v.strip()
                    if not v: continue
                    try:
                        val_int = eval(v.replace('|', '|'), {}, env_vars)
                        vals_clean.append(str(val_int))
                    except Exception as e:
                        vals_clean.append(v)
                js_content += f"    '{board_name}': [{', '.join(vals_clean)}],\n"
            else:
                if val_block.startswith('"') and val_block.endswith('"'):
                    pass # Keep the quotes!
                else:
                    try:
                        val_int = eval(val_block.replace('|', '|'), {}, env_vars)
                        if isinstance(val_int, str):
                            val_block = f"'{val_int}'"
                        else:
                            val_block = str(val_int)
                    except Exception:
                        pass
                js_content += f"    '{board_name}': {val_block},\n"
        js_content += "  },\n"
        
    
    # --- Parse boards.qrc for SVG filenames ---
    qrc_path = os.path.join(project_dir, "src", "boards", "boardPics", "boards.qrc")
    js_content += "  boardImagesMap: {\n"
    images_map_lines = []
    svgs_map_lines = []
    if os.path.exists(qrc_path):
        with open(qrc_path, "r", encoding="utf-8") as f:
            qrc_content = f.read()
        
        svg_base_dir = os.path.dirname(qrc_path)
        qrc_pattern = re.compile(r'<file\s+alias="([^"]+)">([^<]+)</file>')
        for match in qrc_pattern.finditer(qrc_content):
            alias = match.group(1)
            filename = match.group(2)
            js_content += f"    '{alias}': '{filename}',\n"
            
            svg_file_path = os.path.join(svg_base_dir, filename)
            if os.path.exists(svg_file_path):
                with open(svg_file_path, "r", encoding="utf-8") as svg_f:
                    svg_text = svg_f.read()
                    svg_text = svg_text.replace('`', '\\`').replace('$', '\\$')
                    svgs_map_lines.append(f"    '{alias}': `{svg_text}`,\n")
                    
    js_content += "  },\n"
    js_content += "  boardSVGsMap: {\n"
    for line in svgs_map_lines:
        js_content += line
    js_content += "  },\n"

    
    js_content += "};\n"
    
    with open(out_js_path, "w", encoding="utf-8") as f:
        f.write(js_content)

if __name__ == "__main__":
    generate_shared_js("E:/PROGETTI/OpenFIRE-ESP32/OpenFIRE-Firmware-ESP32/lightgun", "F:/OpenFIREFirmware/lightgun/webapp")
