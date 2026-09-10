import os
import shutil
import glob

def build_board_pics(project_dir, webapp_dir):
    src_dir = os.path.join(project_dir, "src", "boards", "boardPics")
    if not os.path.exists(src_dir):
        src_dir = os.path.join(project_dir, "..", "shared_boards", "boardPics")
        
    dest_dir = os.path.join(webapp_dir, "boardPics")
    
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
        
    count = 0
    for svg_file in glob.glob(os.path.join(src_dir, "*.svg")):
        shutil.copy(svg_file, dest_dir)
        count += 1
        
    print(f"[WebApp Packer] Copied {count} SVGs to {dest_dir}")

if __name__ == "__main__":
    build_board_pics("F:/OpenFIREFirmware/lightgun", "F:/OpenFIREFirmware/lightgun/webapp")
