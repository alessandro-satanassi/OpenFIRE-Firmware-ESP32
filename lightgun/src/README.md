<a id="english-version"></a>

[Back to Home](../../README.md#english-version) / [Lightgun Firmware](../README.md#english-version) / **Operational Manual**

<p align="center">
  <a href="#english-version"><img src="../../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# OpenFIRE - The Enclosed Instruction Book!

*... adapted for ESP32 7.0.0 from the original [OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/blob/OpenFIRE-dev/OpenFIREmain/README.md) project repository:*

## Table of Contents:
 - [IR Emitter Setup](#ir-emitter-setup)
 - [Board Configuration](#board-configuration)
   - [Configuration with the WebApp](#configuration-with-the-webapp)
   - [Special boot modes](#special-boot-modes)
   - [Camera, display and startup settings](#camera-display-and-startup-settings)
 - [First-time Setup](#first-time-setup)
 - [Operations Manual](#operations-manual)
   - [Run Modes](#run-modes)
   - [Default Buttons](#default-buttons)
   - [Default Buttons in Pause Mode](#default-buttons-in-pause-mode-hotkey)
   - [Controls for Simple Pause Menu](#controls-for-simple-pause-menu)
   - [How to Calibrate](#how-to-calibrate)
   - [IR Camera Sensitivity](#ir-camera-sensitivity)
   - [Profiles](#profiles)
   - [Software Toggles](#software-toggles)
   - [Saving Settings to Flash](#saving-settings-to-flash)
   - [Test Mode](#test-mode)
 - [Common Problems](#common-problems)
 - [Known Limitations](#known-limitations)
 - [Technical Details & Assorted Errata](#technical-details--assorted-errata)
   - [Serial Handoff (Mame Hooker) Mode](#serial-handoff-mame-hooker-mode)
   - [Multiple Guns and Multiplayer](#multiple-guns-and-multiplayer)

## IR Emitter setup

The IR emitters can be arranged in either two ways:

 - **Square / rectangular layout:** two LEDs at the top and two at the bottom of the display, aligned in two columns. The recommended arrangement uses two pairs centered on the top and bottom edges to form a vertical rectangle (width smaller than height), as illustrated by the alignment assistant. A wider rectangle with the four LEDs at the screen corners is also supported. Avoid a perfect square.
 - **Diamond layout:** one LED at the center of each of the four sides of the display (top, bottom, left and right), not at the corners.

**Moving Square emitters to the screen corners:** keep **Square** selected in the calibration profile and perform a **new calibration**, then save. You do not need to choose a separate vertical/horizontal variant: the firmware derives that distinction from the calibration data. This automatic handling applies within Square; it does not replace the choice between Square and Diamond.

Use **940 nm IR emitters for DFRobot/Wii** and **850 nm for PAJ7025R2/R3**. Match the emitters to your camera; more sensitivity cannot compensate for an unsuitable wavelength or placement.

With a DFRobot/Wii camera and a small PC monitor, you can use 2 Wii sensor bars; one on top of your screen and one below. However, if you're playing on a TV, you should consider building or buying a set of high power black IR LEDs and arranging them like (larger) sensor bars at the top and bottom of the display.

The **OpenFIRE ESP32 WebApp** has an alignment assistant that can be used to help align your emitters to the display (by selecting ***Help->Open IR Emitter Alignment Assistant***, or the **Emitter Alignment** button in the menu bar). 

<table>
  <tr>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_app_001.png">
        <img src="../docs/img/IR_Emitter_app_001.png">
      </a>
    </td>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_app_002.png">
        <img src="../docs/img/IR_Emitter_app_002.png">
      </a>
    </td>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_monitor.jpeg">
        <img src="../docs/img/IR_Emitter_monitor.jpeg">
      </a>
    </td>
  </tr>
</table>

## Board Configuration

The official configuration tool for OpenFIRE ESP32 7.0.0 is the **OpenFIRE ESP32 WebApp**. It configures pins, buttons, camera, calibration profiles and force feedback, and provides input and IR tests. The online WebApp needs a browser with Web Serial, such as Chrome or Edge on a computer; with other browsers (for example Firefox, Safari or phone browsers) use the [offline WebApp](#configuration-with-the-webapp) stored in the gun, which does not need Web Serial and works with current browsers, including on a phone.

If you prefer a program to install, or your browser is too old to open either WebApp, a compatible desktop App, derived from the App of the original OpenFIRE project, can be downloaded from the [Tools page](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-Tools/?lang=en), also reachable from the [project hub](https://alessandro-satanassi.github.io/OpenFIRE-ESP32/?lang=en); it does not depend on the browser, but it does need an operating system supported by the App. Use the **OpenFIRE App CUSTOM for ESP32 7.x firmware** section: the other App builds on that page, including those of the original project, use the previous configuration protocol and are not compatible.

Connecting the App puts the gun into *Docked* configuration mode. Save your changes and wait for confirmation before disconnecting or switching off. If an operation fails, follow the displayed recovery instructions and verify the settings after reconnecting; do not assume an unconfirmed save succeeded.

<a id="configuration-with-the-webapp"></a>

### Configuration with the WebApp

**Online, from a computer:**

1. Start the lightgun normally, without holding a boot combination. Connect its USB OTG port with a data cable, or connect the paired wireless dongle.
2. Open the [WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/) in a browser supporting Web Serial, such as Chrome or Edge on a computer.
3. Connect and authorize the lightgun's serial port (or the paired dongle's port). The launcher identifies the firmware and opens the corresponding published WebApp.
4. Configure, calibrate and test the gun. Save before disconnecting. Close other configuration pages and serial applications that might be using the same device.

**Offline, using the WebApp stored in the lightgun:**

Start with **B held for about 2 seconds**, as explained below. Internet access and a downloaded configuration App are not required.

- **Wi-Fi, including phones:** join **OpenFIRE_Config**. If a welcome page appears, accept using this network even without Internet (on some Android phones: the menu option to use the network as it is). Then open your **normal browser**, not the welcome window, and enter **http://openfire.local/** or **http://192.168.4.1/**. If the name does not resolve, use the IP address. The welcome page is deliberately simple because some phones use a limited captive-portal browser.
- **USB OTG network:** on a computer whose operating system supports USB NCM, connect the lightgun directly and open **http://192.168.7.1/**. **http://openfire** or **http://openfire.local/** may also work, depending on name resolution. Wi-Fi is not required for this path.
- **Wireless gun:** when running on battery, keep the dongle connected and paired; it carries the game inputs while the phone/computer accesses the gun's Wi-Fi network. The dongle's serial port remains available.
- **USB limitation:** in this special mode, the lightgun's USB network replaces its serial port. HID pointing still works, but the gun's USB serial connection, online Web Serial configuration and serial tools such as MAMEHOOKER are unavailable until a normal restart.

Use one configuration page at a time. Once finished, save, disconnect the App and **restart without holding any buttons** to return to normal operation.

<a id="special-boot-modes"></a>

### Special boot modes

Press the lightgun buttons **before powering on or resetting**, and keep them held for about **2 seconds**, until the requested mode starts. These shortcuts require a lightgun already running the 7.0.0 firmware and correctly mapped, working buttons.

| Buttons held at startup | Result |
| --- | --- |
| None | Normal operation. |
| **Trigger + A** | Firmware update mode; the OLED, if fitted, shows **Ready for firmware update**. Connect the lightgun's own USB OTG port and use the [Web Flasher](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/). |
| **B** | Starts the integrated offline WebApp over Wi-Fi and, with a USB cable, USB NCM. On the OLED, the top bar is inverted and shows a gear icon. |

Firmware update takes priority if both combinations are held. **A/B are the lightgun's mapped buttons, not the board's physical BOOT button.** On a blank board, an older firmware or when the shortcuts cannot work, use the board's physical BOOT/RESET procedure instead. On a lightgun running 7.0.0 normally, the Web Flasher can also switch it to flashing mode by itself, and the WebApp can do the same with **Restart Microcontroller in Firmware Update Mode** in the *Gun Tests* tab. After this software restart the serial port changes: start the installation again and select the new port. The firmware cannot be flashed through the wireless dongle or the Wi-Fi configuration page.

For installation, select the image matching the **board and flash/PSRAM variant**, not the camera. A normal update and a clean installation use the same image; a clean installation erases all settings and calibration. When moving from 6.2.1 to 7.0.0, a clean installation is recommended: note your existing settings first, then configure and calibrate again.

### Camera, display and startup settings

- **Camera:** in **Gun Settings → CAMERA Model**, choose **DFRobot SEN0158 / Wii**, **PAJ7025R2** or **PAJ7025R3** to match the installed hardware (enable **Unlock CAMERA Modification** if the controls are locked). Also set the communication pins in **Board Layout**: SDA/SCL for I2C, or RX (MISO), TX (MOSI), SCK and CSn for SPI. These are signal connections, in addition to power and ground. Save to apply the selection; recalibrate after changing the camera, lens or emitter placement. R3's wider lens helps at shorter distances, but is not a guarantee of better accuracy or greater range. After a clean installation the selected camera is **DFRobot SEN0158 / Wii**: with a PAJ7025R2 or R3, select it and set its SPI pins before calibrating, otherwise the WebApp may report **Device Error: Camera not available!**
- **OLED display:** the display is disabled by default. Enable **Enable OLED Display (SSD1306 128x64)** in **Gun Settings → I2C Peripherals** and, in **Board Layout**, assign **Peripherals SDA** and **Peripherals SCL** to the pins the display is wired to; then save and restart. If the display stays blank, try **Use Alternative Device Address**.
- **Startup mode:** choose **Absolute Mouse** (default), **Gamepad (right stick)** or **Gamepad (left stick)** in **Gun Settings → Input / Output**, then save. This selects the output used at the next normal startup; serial game commands can still change it during play.
- **Wireless pedal** (**Gun Settings → Input / Output**): enable this option only if you use that accessory, with both wired pedal inputs unmapped. Otherwise disable it to avoid an unnecessary pedal search at startup. A mapped wired pedal takes precedence. The wireless pedal works only when the gun is connected through the dongle; with the gun connected to the computer by USB cable, use a wired pedal.

## First-time Setup
After flashing a new board, or after a clean installation, the gun has no calibration yet. It waits for the first calibration and does not move the pointer; the board's RGB LED, if present, blinks orange.

1. Connect the gun to the WebApp. If you use custom pins, map at least the *Trigger* and *Button A* in *Board Layout*; select your camera and check the buttons in the *Gun Tests* tab, then save.
2. Start the calibration with any of the *Calibrate Profile* buttons in the *Calibration Profiles* tab, or pull the trigger to run the standalone calibration (see [How to Calibrate](#how-to-calibrate)).
3. Save the calibration.

If the camera is not available (for example because the wrong model is selected), pulling the trigger does not start the calibration: correct the camera settings in the WebApp first.

## Operations Manual
By default, the light gun operates as an absolute positioning mouse (like a stylus!) until the button/combination is pressed to enter pause mode. A different output can be selected with the saved **Startup mode** setting. Alternatively, the gun can be signaled to output using its corresponding HID Gamepad device using a Serial Feedback Distributor program such as MAMEHOOKER - see the [Serial Handoff](#serial-handoff-mame-hooker-mode) section for more info.

**How the computer sees the gun:** the lightgun, or its dongle, is recognised as standard USB devices (an absolute-positioning mouse, a keyboard and a gamepad), so no drivers are needed and games and emulators can use it like any mouse or controller. To make the solenoid, rumble, LEDs and OLED counters follow what happens in the game, you also need a force feedback program such as MAMEHOOKER (see [Serial Handoff](#serial-handoff-mame-hooker-mode)); without it, force feedback reacts only to the trigger.

Any serial terminal (Arduino IDE's Serial Monitor, *PuTTY,* *screen,* etc.) can be used to see information while the gun is paused and during standalone calibration.

Note that the buttons in pause mode (and to enter pause mode) activate when the last button of the combination releases. This is used to detect and differentiate button combinations vs a single button press.

* Note: At its peak, the mouse position updates at 209Hz (both over USB and wireless), or roughly every ~4.8ms, so it is extremely responsive.

### Run modes
The gun has the following modes of operation:
1. Normal - The mouse position updates from each frame from the IR positioning camera (no averaging)
2. Averaging - The position is calculated from a 2 frame moving average (current + previous position)
3. Averaging2 - The position is calculated from a weighted average of the current frame and 2 previous frames
4. Processing - Test mode for use with the WebApp (this mode is prevented from being assigned to a profile)

The averaging modes are subtle but do reduce the motion jitter a bit without adding much if any noticeable lag.

> The ESP32 port of OpenFIRE automatically incorporates advanced anti-jitter algorithms, so it is recommended to use **Normal** mode.

### Default Buttons
- Trigger: Left mouse button
- A: Right mouse button (In low buttons mode, Start if pressed offscreen)
- B: Middle mouse button (In low buttons mode, Select if pressed offscreen)
- C/Reload: Mouse button 4/Side Button 1/Back
- Pump Action (Cabela's or alike): Right mouse button
- Start: 1 key
- Select: 5 key
- Up/Down/Left/Right: Keyboard arrow keys
- Pedal Main: Mouse button 4/Side Button 1/Back
- Alt Pedal: Mouse button 5/Side Button 2/Forward
- C + Start: Esc key

**Low Buttons Mode** (**Gun Settings → Input / Output**) is meant for guns with few buttons: when it is enabled, A and B pressed while pointing off-screen act as Start and Select.

Pause mode can be entered by either pressing C + Select by default, pressing the *Home Button* if used in current pin layout, or *holding the trigger plus the A Button with **no IR points in sight*** if hold-to-pause is enabled - pointing the gun towards the ground is recommended here.

**Which pause mode is active?** By default the gun uses the **Hotkey** pause mode: once paused, each button or combination listed below performs its function directly. To use the **Simple Pause Menu** instead (a list of options scrolled with the buttons, easiest to use with an OLED display), enable **Simple Pause Menu** in **Gun Settings → UI and UX**. Entering pause by holding Trigger + A requires **Hold to Pause Enabled** in the same group, where the hold time is also set. Save after changing these options.

#### Default Buttons in Pause mode (Hotkey)
- A, B, Start, Select: select a profile
- Start + A: Normal gun mode (averaging disabled)
- Start + B: Normal gun with averaging, switch between the 2 averaging modes (use serial monitor to see the setting)
- B + Down: Decrease IR camera sensitivity (use a serial monitor to see the setting)
- B + Up: Increase IR camera sensitivity (use a serial monitor to see the setting)
- C/Reload: Exit pause mode
- Left: Toggle Rumble *(when no rumble switch is detected)*
- Right: Toggle Solenoid *(when no solenoid switch is detected)*
- Trigger: Begin calibration
- Start + Select: save settings to non-volatile flash storage space

#### Controls for Simple Pause Menu
- A or Up: Navigate Cursor Up
- B or Down: Navigate Cursor Down
- Trigger: Select option
- C: Exit pause mode
  - Holding A or B for half the duration of the hold-to-pause time (so ~2s by default) will also exit the simple pause menu.
Available options in simple pause menu are as follows, from first option to last before rolling back:
* Calibrate current profile (always the initial option)
* Switch profiles (submenu)
  * Select from profile 1-4 using the navigation buttons/trigger to select, or press C to back out.
* Save settings to non-volatile memory
* Toggle Rumble *(when rumble is enabled & no switch is detected)*
* Toggle Solenoid *(when solenoid is enabled & no switch is detected)*
* Send escape key signal to the PC

### How to calibrate
##### These instructions apply to the standalone on-board calibration process; the Calibration screens in the OpenFIRE WebApp have a similar procedure with more info to guide the user through the process.

**Before you start:** stand in front of the centre of the screen, hold the gun without rotating it around the barrel, and aim carefully at each target. In the WebApp and desktop App, this posture reminder stays above the bottom instructions throughout calibration, including final aiming verification.

1. Select the profile to calibrate - either through pressing A/B/Start/Select in the Hotkey Pause Mode, or selecting in the Simple Pause Menu - and pull the trigger to begin calibration. Alternatively, calibration can be started from the WebApp (*Calibration Profiles* tab, *Calibrate Profile* buttons). 
2. Aim at the centre of the screen and pull the trigger while keeping a steady aim.
3. The cursor moves to the four edges of the screen: top, bottom, left, right. At each edge it draws a small circle tangent to the edge. Shoot the point where the circle touches the screen edge, **not the centre of that circle**.
4. When the cursor returns to the centre of the screen, aim at the centre of its small circle and shoot to finish calibration.
5. The new calibration profile will be applied and you'll be able to test the tracking. A good sign of a good calibration is maintaining as close to line-of-sight accuracy as possible when aiming at the screen edges and corners.
   - If the calibration is good, pull the trigger to confirm.
   - If you want to start calibration over, press the A or B button in cali verification to restart from the center point in Step 2.
   - Calibration can be canceled outright by pressing C/Reload at any time, or the A/B buttons any time before cali verification (during the first calibration after a clean installation, these buttons do not cancel the procedure; in the WebApp or desktop App, you can still cancel by closing the calibration window).

**Moving cursor in standalone calibration.** By default, the mouse pointer continuously traces a small circle while waiting for each target shot. At the centre of the screen, aim at the centre of the circle; at the four edges, aim at its point of contact with the edge. The animation stops while travelling to the next target and during final aiming verification. It does not change the calibration calculations or the calibration interface in the WebApp or desktop App.

**Calibration guidance on the OLED.** If an SSD1306 display is fitted and enabled, the top status bar remains visible and the area below guides calibration. The display instructions are always in English.

- **Where to aim:** a small TV outline and crosshair indicate the target position on your monitor, not your live aim. Progress runs from **1/6** to **6/6**, in this order: `CENTER`, `TOP`, `BOTTOM`, `LEFT`, `RIGHT`, `CENTER` (final centre). `SHOOT` tells you to shoot the target.
- **Final verification:** after the sixth target, `CHECK AIM` asks you to check your aim on the monitor. `TRIGGER: CONFIRM` means pull the trigger to confirm; `A/B: REPEAT` means press A or B to repeat calibration.
- **IR points:** white points remain visible over the instructions during calibration and verification. They show both detected and reconstructed emitter positions, without distinguishing between them: four points do not necessarily mean that all four emitters are actually visible to the camera.

These are visual aids only: the calibration sequence, buttons and saving procedure described in this guide are unchanged.

**Calibration from the WebApp or the desktop App.** With firmware 7.0.0, the calibration window started from the WebApp (or from the compatible desktop App) also checks the IR emitters before every target shot:
- **Crosshair colour:** red when the camera does not see all four emitters; otherwise it follows the brightness of the weakest emitter on a continuous scale, from red through orange and a light yellow-green to full green, without sudden changes.
- **Crosshair highlight:** a thin dashed outer ring slowly rotates in the same dynamic colour as the crosshair, without pulsing or changing brightness. The original crosshair stays still and keeps its size and position, also at the screen edges. The coloured outer ring is absent during aiming verification and stays still in the WebApp when the system asks for reduced motion.
- **Progress:** six dots above the stage title represent the targets: centre, top, bottom, left, right, final centre. Completed targets are filled grey, the current one is green and the others are outlined grey. This green indicates progress, not IR quality. The titles run from “Calibration: step 1 of 6” through “Calibration: step 6 of 6”, including the initial centre shot; all six dots are filled during “Verify aiming:”.
- **Shot accepted:** a solid, stationary white ring briefly flashes around the target just acquired (about 120 ms). The crosshair stays there until the flash ends, then moves to the next target. The last centre also receives this confirmation before the crosshair follows the aiming cursor in verification. The progress dots do not flash. With reduced motion, the WebApp shows the white ring without fading, for the same duration. Rejected shots and resets do not trigger it. It confirms the target was acquired, not that the calibration was saved: verify your aim, confirm with the trigger, then save the settings as described below.
- **Emitters panel:** a square with rounded corners in the bottom right corner, labelled “IR LEDs” at its centre, showing the four emitters in their layout (Square or Diamond). Size and fill of each circle follow the camera as in the [IR test](#test-mode); each emitter seen has the same colour scale as the crosshair, from its own brightness, and the crosshair follows the weakest one. An emitter that is not seen is a red dashed circle with a red X.
- **Shots refused:** a target shot (centre, the four edges and the final centre) is refused only while the camera does not see all four emitters (red crosshair). Calibration then does not move on and the window explains why: check the emitters, your distance and the [camera sensitivity](#ir-camera-sensitivity), then shoot again. The final confirmation in the verification step is not checked. A weak emitter (orange crosshair) does not stop calibration, but the most reliable result comes with a green crosshair.

The **legend in the bottom left corner**, shown throughout the procedure including aiming verification, explains the IR signal intensity (Weak → Strong), the size of the IR point (Small → Large), and the red dashed circle with an X (IR not detected). The crosshair colour indicates the emitter with the weakest signal; all four emitters must be detected to acquire each target. The legend uses a normal typeface and does not change how calibration works.

The calibration started from the gun itself (pause mode) works as before, without these checks.

<p align="center">
  <img src="../docs/img/webapp_calibration_ir.png" alt="Calibration in the WebApp: red crosshair, panel with the four emitters, one of them not seen, and the message that the shot was refused because the camera does not see all four IR emitters" width="70%">
</p>

Remember to save your calibration and current profile afterwards, either by saving from the WebApp, pressing Start+Select in Hotkey Pause Mode, or selecting the third "Save Settings" option in the Simple Pause Menu.

### IR Camera Sensitivity
The IR camera sensitivity can be adjusted. It is recommended to adjust the sensitivity as high as possible. If the IR sensitivity is too low then the pointer precision can suffer. However, too high of a sensitivity can cause the camera to pick up unwanted reflections that will cause the pointer to jump around. It is impossible to know which setting will work best since it is dependent on the specific setup. It depends on how bright the IR emitters are, the distance, camera lens, and if shiny surfaces may cause reflections.

A sign that the IR sensitivity is too low is if the pointer moves in noticeable coarse steps, as if it has a low resolution to it. If you have the sensitivity level set to max and you notice this, the IR emitters may not be bright enough.

A sign that the IR sensitivity is too high is if the pointer jumps around erratically. If this happens only while aiming at certain areas of the screen, this is a good indication that a reflection is being detected by the camera. If the sensitivity is at max, step it down to high or minimum. Obviously, the best solution is to eliminate the reflective surface. The WebApp's IR test (see [Test Mode](#test-mode)) can help diagnose this problem: it displays the IR points seen by the camera, including how large and bright each one is.

### Profiles
The main OpenFIRE builds are configured with 4 calibration profiles available. Each profile has its own calibration data, run mode, and IR camera sensitivity settings. Each profile can be selected from pause mode by pressing the associated button (A/B/Start/Select), or selecting them via the profiles submenu in simple pause menu. In the WebApp, the **Calibration Profiles** tab shows and changes the settings of every profile, including **Sensitivity** and **Run Mode**: this is the easiest way to check them, since in pause mode they are only reported on a serial monitor. Save after changing them.

### Software Toggles
Hardware features can be toggled at runtime, even without hardware switches defined!

While in pause mode, the toggles are as follows (color indicating what the board's builtin LED lights up with):
- Left D-Pad: **Rumble Toggle** (Salmon) - Enables/disables the rumble functionality. When enabled, the motor will engage for a short period.
- Right D-Pad: **Solenoid Toggle** (Yellow) - Enables/disables the solenoid force feedback. When enabled, the solenoid will engage for a short period.
These can also be done from the respective setting in the Simple Pause Menu.

The current state of these settings is saved when committed to, and pulled from flash storage space at boot.

#### Saving Settings to Flash
The calibration data, profile settings, and extended gun options like custom pins mapping and rumble intensity, can be saved in non-volatile memory by pressing Start + Select while the gun is in Hotkey Pause Mode (or using Save Settings in the Simple Pause Menu), or by saving and receiving confirmation in the WebApp. The currently selected calibration profile is saved as the default for when the light gun is plugged in - gun settings (pins mapping, force feedback, etc.) applies to *all profiles.* In the WebApp and in the desktop App, **Save and Send Settings** pulses while there are unsaved changes, as a reminder to save.

#### Test Mode
Test Mode shows the IR points as seen by the camera. Open it from the WebApp's *Gun Tests* tab with **Open IR Camera Tester...**. It is very useful for aligning the camera when building your light gun, for testing that the camera tracks all 4 points properly, and for spotting possible reflections. The validity of the test points shape (rectangle in Square layout, diamond in Diamond layout) depends on the current profile used and its IR layout setting. Save any change to the IR layout or camera sensitivity before opening the test: the test uses the settings stored in the gun.

With firmware 7.0.0, each emitter circle also shows what the camera measures:
- **Size:** follows the size of the light spot seen by the camera. Moving away from the screen usually makes it only slightly smaller, because the spot size depends mostly on the LED brightness and on the camera optics.
- **Fill:** brightest at the centre (peak brightness of the spot) and fading towards the edge (average brightness). A full circle means a bright, well detected LED; an almost empty circle means a weak LED close to the detection limit.
- **Dashed circle with a red X:** that emitter is not seen by the camera; its position is only estimated.

The gray circle shows where the gun is aiming and the red circle marks the centre of the four emitters. With the DFRobot/Wii camera, size and brightness come from the camera's full data format and are converted to the same scale. With older firmware the classic circles are shown.

<p align="center">
  <img src="../docs/img/webapp_ir_test.png" alt="IR camera test in the WebApp: three emitters seen, with circles of different size and brightness, and one emitter not seen, marked with a red X" width="70%">
</p>

## Common Problems

- **After the installation the gun does nothing and the pointer does not move:** a new or clean-installed gun waits for its first calibration (the board's RGB LED, if present, blinks orange). Pull the trigger to calibrate, or calibrate from the WebApp; see [First-time Setup](#first-time-setup).
- **The WebApp does not find the gun's port:** use a USB **data** cable (some cables only charge) and a computer browser with Web Serial, such as Chrome or Edge. Browsers without Web Serial (for example Firefox, Safari or phone browsers) can use the [offline WebApp](#configuration-with-the-webapp) stored in the gun; for a browser too old for both, see the desktop App in [Board Configuration](#board-configuration). On boards with two USB connectors, use the one connected to the ESP32-S3's own USB (OTG): on the DevKitC-1, seen from above with the connectors pointing forward, it is the one on the right, in front of the built-in NeoPixel LED (the board picture in the WebApp's **Board Layout** tab labels it **USB OTG**). Close other WebApp pages, serial monitors and MAMEHOOKER, which can keep the port busy. If the gun was started holding **B**, its serial port is replaced by the offline WebApp network: restart without holding any buttons.
- **WebApp message "Device Error: Camera not available!":** after a clean installation the camera model is **DFRobot SEN0158 / Wii**. With a PAJ7025R2 or R3, select it in **Gun Settings → CAMERA Model**, set its SPI pins in **Board Layout**, then save and restart. Otherwise check the camera power and that its wires are not swapped or assigned to other functions.
- **No IR points in the test, or some are missing:** check that the emitters are on and suit the camera (940 nm for DFRobot/Wii, 850 nm for PAJ7025R2/R3), that the camera can see the whole emitter layout, and that the camera sensitivity is not too low. Save the sensitivity before opening the test.
- **Calibration does not move on to the next target:** in a calibration started from the WebApp or the desktop App, a target shot is refused while the camera does not see all four emitters (red crosshair). Look at the emitters panel in the bottom right corner: an emitter marked with a red X must come back into the camera view (step back if it leaves the view at the screen edges); check also the emitter, the distance and the [camera sensitivity](#ir-camera-sensitivity). An emitter shown in orange is weak: the shot is accepted, but improving it gives a steadier aim. See [How to calibrate](#how-to-calibrate).
- **The pointer jumps or shakes:** usually a reflection or another light source (window, lamp, shiny surface) seen as an extra IR point. Look for it in the [IR test](#test-mode), remove or cover it, or lower the [camera sensitivity](#ir-camera-sensitivity). Recalibrate after moving the emitters.
- **The gun does not pair with the dongle:** plug in the dongle first and wait about 15 seconds, then switch on the gun while it is not connected to the computer by USB (with USB it works wired). After unplugging or restarting the dongle, restart the gun as well. Use one dongle for each gun, and lightgun and dongle firmware from the same release ([dongle guide](../../dongle/README.md#english-version)).
- **The wireless pedal is not found:** it works only with a gun connected through the dongle. Enable **Wireless pedal**, leave both wired pedal inputs unmapped and save; switch on the pedal before the gun. After restarting the pedal, restart the gun as well.
- **The OLED display stays blank:** the display is disabled by default; enable it as explained in [Camera, display and startup settings](#camera-display-and-startup-settings).
- **The temperature reading is far too high** (for example about 90 °C at room temperature): give the TMP36 its own ground wire to the board, not shared with the camera ([wiring tip](../README.md#optional-modules)).
- **The wireless connection drops or has a short range:** keep a small clear area around the antenna of the gun's board, with no wires over or right next to it and no metal parts covering it ([antenna tip](../README.md#mandatory-components-essentials)); the same applies to a home-made dongle or pedal. Also use solenoid power wires of at least 24 AWG, as thin wires can cause disconnections.

Some behaviours are limitations of the system rather than faults: see [Known Limitations](#known-limitations).

## Known Limitations

- The refined tracking applies to the Square layout; the Diamond layout uses the original OpenFIRE tracking.
- Wireless play requires an ESP32-S3 board; on RP2040 boards the firmware works by USB cable only.
- Up to four guns can be used together; each needs its own dongle, and the wireless pedal works only with a gun connected through the dongle.
- The firmware is installed through each device's own USB port, not through the dongle or the Wi-Fi configuration page.
- In the offline WebApp mode over USB, the gun's USB serial port, and therefore MAMEHOOKER, is unavailable until a normal restart.
- Configuration Apps for firmware older than 7.0.0, including those of the original OpenFIRE project, are not compatible; use the WebApp or the compatible desktop App (see [Board Configuration](#board-configuration)).

## Technical Details & Assorted Errata

### Serial Handoff (Mame Hooker) Mode

Use normal startup for the gun’s USB serial port; it is replaced by USB networking in the special offline WebApp mode. Disconnect the configuration App before using the same port for MAMEHOOKER.
The gun will automatically hand off control to an instance of Mame Hooker that's connected once a start code has been detected! If available, the onboard LED and any *non-static* external NeoPixels will change to a mid-intensity white to signal serial handoff mode (unless any LED events trigger it to change, which will follow those thereafter).

If you aren't already familiar with Mame Hooker, **you'll need compatible inis for each game you play** and **the gun's COM port should be set to match the player number** (COM1 for P1, COM2 for P2, etc.)! COM port assignment can be done in Windows via the Device Manager, or Linux via settings in the Wine registry of the prefix your game/Mame Hooker is started in. [Consult the wiki page on MAMEHOOKER for more information!](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/wiki/MAMEHOOKER_Documentation_EN) For Linux users wanting to use their gun with native emulators' force feedback (currently MAME, Flycast, or their RetroArch ports), consider trying [QMamehook](https://github.com/SeongGino/QMamehook).

### Multiple Guns and Multiplayer
To play with two (or more) OpenFIRE guns on the same computer, each gun must report a different identity. If several devices share the same name and/or Product/Vendor ID (the **USB Implementer's Forum (USB-IF) identifiers**), applications that read individual mouse devices, such as RetroArch and TeknoParrot, cannot tell them apart.

1. Connect one gun at a time to the WebApp and open **Gun Settings → TinyUSB Identifier**.
2. Choose the player number with the **P1**–**P4** buttons: **P1** for the first gun, **P2** for the second, and so on. After a clean installation every gun is **P1**. **Advanced View** lets you enter a custom Product ID and name instead.
3. Save and restart the gun.

The player number also changes the Start and Select keys (P1: 1 and 5, P2: 2 and 6, P3: 3 and 7, P4: 4 and 8; a custom ID keeps the P1 keys), is shown on the OLED displays and by the pedal LEDs, and is the number to use for the MAMEHOOKER COM port (COM1 for P1, COM2 for P2, etc.).

**Wireless play with several guns:** up to four guns can be used wirelessly, each with **its own dongle** and, if you want, its own wireless pedal; each dongle pairs with one gun only. It does not matter which dongle a gun pairs with: the dongle takes on the identity of the gun (name, player number and USB serial number), so the computer always sees P1 as P1. Pairing and the choice of radio channel are automatic, and the dongles may end up on the same channel or on different ones: there is nothing to set.

If you like, set up one gun at a time: plug in the first dongle, wait about 15 seconds, switch on the first gun and its pedal, then do the same for the next one. Each dongle then chooses its channel taking into account any interference, including that produced by the dongles and guns already in use, and each wireless pedal pairs with the intended gun, because a pedal pairs with the first gun that searches for it. This is a recommendation, not a requirement.

---
### Questions or Issues?
For technical support and to join the discussion, please refer to the [Community & Support Section](../../README.md#community-support-english) in the Main Repository.


---

<a id="versione-italiana"></a>

[Torna alla Home](../../README.md#versione-italiana) / [Lightgun Firmware](../README.md#versione-italiana) / **Manuale Operativo**

<p align="center">
  <a href="#english-version"><img src="../../docs/img/gb.png" width="20" alt="English"> English Version</a> &nbsp;•&nbsp; <a href="#versione-italiana"><img src="../../docs/img/it.png" width="20" alt="Italiano"> Versione Italiana</a>
</p>

# OpenFIRE - Il Manuale di utilizzo!

*... adattato per ESP32 7.0.0 dal progetto originale [OpenFIRE](https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/blob/OpenFIRE-dev/OpenFIREmain/README.md):*

## Indice:
 - [Configurazione Emettitori IR](#configurazione-emettitori-ir-italiano)
 - [Configurazione della Scheda](#configurazione-della-scheda-italiano)
   - [Configurazione con la WebApp](#configurazione-con-la-webapp)
   - [Modalità speciali all’avvio](#modalita-speciali-allavvio)
   - [Telecamera, display e impostazioni di avvio](#telecamera-display-e-impostazioni-di-avvio)
 - [Prima Configurazione](#prima-configurazione-italiano)
 - [Manuale Operativo](#manuale-operativo-italiano)
   - [Modalità di Funzionamento](#modalità-di-funzionamento-italiano)
   - [Pulsanti Predefiniti](#pulsanti-predefiniti-italiano)
   - [Pulsanti Predefiniti in Modalità Pausa](#pulsanti-predefiniti-in-modalità-pausa-hotkey-italiano)
   - [Controlli per il Menu di Pausa Semplificato](#controlli-per-il-menu-di-pausa-semplificato-italiano)
   - [Come Calibrare](#come-calibrare-italiano)
   - [Sensibilità della Telecamera IR](#sensibilità-della-telecamera-ir-italiano)
   - [Profili](#profili-italiano)
   - [Interruttori Software (Toggle)](#interruttori-software-toggle-italiano)
   - [Salvataggio delle Impostazioni nella Flash](#salvataggio-delle-impostazioni-nella-flash-italiano)
   - [Modalità di Test](#modalità-di-test-italiano)
 - [Problemi Comuni](#problemi-comuni-italiano)
 - [Limiti noti](#limiti-noti-italiano)
 - [Dettagli Tecnici e Note Varie](#dettagli-tecnici-e-note-varie-italiano)
   - [Modalità Serial Handoff (Mame Hooker)](#modalità-serial-handoff-mame-hooker-italiano)
   - [Più pistole e multigiocatore](#modifica-dell-id-usb-per-pistole-multiple-italiano)


<a id="configurazione-emettitori-ir-italiano"></a>

## Configurazione Emettitori IR

Gli emettitori IR possono essere disposti in due modi:

 - **Layout Square / rettangolare:** due LED in alto e due in basso sullo schermo, allineati su due colonne. La disposizione consigliata usa due coppie centrate sui bordi superiore e inferiore per formare un rettangolo verticale (base minore dell'altezza), come illustrato dall'assistente di allineamento. È supportato anche un rettangolo più largo con i quattro LED agli angoli dello schermo. Evita un quadrato perfetto.
 - **Layout Diamond / a diamante:** un LED al centro di ciascuno dei quattro lati dello schermo (alto, basso, sinistra e destra), non agli angoli.

**Spostando gli emettitori Square agli angoli dello schermo:** mantieni **Square** selezionato nel profilo ed esegui una **nuova calibrazione**, poi salva. Non occorre scegliere una variante verticale/orizzontale separata: il firmware ricava questa distinzione dai dati di calibrazione. Questa gestione automatica vale all'interno di Square; non sostituisce la scelta fra Square e Diamond.

Usa **emettitori IR da 940 nm per DFRobot/Wii** e **da 850 nm per PAJ7025R2/R3**. Abbina gli emettitori alla telecamera: aumentare la sensibilità non compensa una lunghezza d’onda o un posizionamento inadatti.

Con una telecamera DFRobot/Wii e un piccolo monitor per PC, puoi usare 2 barre sensore Wii; una sopra lo schermo e una sotto. Tuttavia, se giochi su una TV, dovresti prendere in considerazione la costruzione o l'acquisto di un set di LED IR neri ad alta potenza e disporli come barre sensore (più grandi) nella parte superiore e inferiore del display.

La WebApp OpenFIRE ESP32 dispone di un assistente di allineamento che può aiutarti ad allineare gli emettitori al display (selezionando *Aiuto → Apri Assistente allineamento emettitori IR*, oppure il pulsante **Allineamento Sensore** nella barra dei menu).

<table>
  <tr>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_app_001.png">
        <img src="../docs/img/IR_Emitter_app_001.png">
      </a>
    </td>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_app_002.png">
        <img src="../docs/img/IR_Emitter_app_002.png">
      </a>
    </td>
    <td valign="middle" width="33%">
      <a href="../docs/img/IR_Emitter_monitor.jpeg">
        <img src="../docs/img/IR_Emitter_monitor.jpeg">
      </a>
    </td>
  </tr>
</table>

<a id="configurazione-della-scheda-italiano"></a>

## Configurazione della Scheda

Lo strumento ufficiale di configurazione per OpenFIRE ESP32 7.0.0 è la **WebApp OpenFIRE ESP32**. Permette di configurare pin, pulsanti, telecamera, profili di calibrazione e force feedback, e di eseguire i test degli ingressi e dei punti IR. La WebApp online richiede un browser con Web Serial, come Chrome o Edge su computer; con altri browser (ad esempio Firefox, Safari o i browser dei telefoni) usa la [WebApp offline](#configurazione-con-la-webapp) contenuta nella pistola, che non richiede Web Serial e funziona con i browser attuali, anche da telefono.

Se preferisci un programma da installare, o il tuo browser è così vecchio da non aprire nessuna delle due WebApp, un'App desktop compatibile, derivata dall'App del progetto OpenFIRE originale, si può scaricare dalla [pagina Tools](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-Tools/?lang=it), raggiungibile anche dal [portale del progetto](https://alessandro-satanassi.github.io/OpenFIRE-ESP32/?lang=it); non dipende dal browser, ma richiede comunque un sistema operativo supportato dall'App. Usa la sezione **OpenFIRE App CUSTOM per firmware ESP32 7.x**: le altre build dell'App presenti in quella pagina, comprese quelle del progetto originale, usano il precedente protocollo di configurazione e non sono compatibili.

Collegando l'App, la pistola entra nello stato *Docked* di configurazione. Salva le modifiche e attendi la conferma prima di disconnettere o spegnere. Se un'operazione fallisce, segui le indicazioni di recupero mostrate e verifica le impostazioni dopo la riconnessione: un salvataggio non confermato non va considerato riuscito.

<a id="configurazione-con-la-webapp"></a>

### Configurazione con la WebApp

**Online, da computer:**

1. Avvia la lightgun normalmente, senza tenere premute combinazioni all'accensione. Collega la sua porta USB OTG con un cavo dati, oppure collega il dongle wireless associato.
2. Apri la [WebApp](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebApp/) con un browser che supporta Web Serial, come Chrome o Edge su computer.
3. Collegati e autorizza la porta seriale della lightgun (oppure quella del dongle associato). Il launcher identifica il firmware e apre la WebApp pubblicata per quella versione.
4. Configura, calibra e prova la pistola. Salva prima di disconnetterti. Chiudi altre pagine di configurazione e programmi seriali che potrebbero utilizzare lo stesso dispositivo.

**Offline, con la WebApp contenuta nella lightgun:**

Avvia tenendo premuto **B per circa 2 secondi**, come spiegato sotto. Non servono Internet né un'App di configurazione da scaricare.

- **Wi-Fi, anche da cellulare:** collegati a **OpenFIRE_Config**. Se compare la pagina di benvenuto, accetta di utilizzare la rete anche senza Internet (su alcuni Android: la voce del menu per utilizzare la rete così com'è). Poi apri il **browser normale**, non la finestra di benvenuto, e digita **http://openfire.local/** oppure **http://192.168.4.1/**. Se il nome non viene risolto, usa l'indirizzo IP. La pagina di benvenuto è volutamente semplice perché alcuni telefoni usano un browser captive portal limitato.
- **Rete USB OTG:** su un computer il cui sistema operativo supporta USB NCM, collega direttamente la lightgun e apri **http://192.168.7.1/**. Possono funzionare anche **http://openfire** o **http://openfire.local/**, secondo la risoluzione dei nomi disponibile. Per questo collegamento non serve il Wi-Fi.
- **Pistola wireless:** se alimentata a batteria, mantieni il dongle collegato e associato: trasporta gli ingressi di gioco mentre il telefono/computer accede alla rete Wi-Fi della pistola. La porta seriale del dongle rimane disponibile.
- **Limite USB:** in questa modalità speciale la rete USB della lightgun sostituisce la sua porta seriale. Il puntamento HID continua a funzionare, ma il collegamento seriale USB della pistola, la configurazione online Web Serial e gli strumenti seriali come MAMEHOOKER non sono disponibili fino al riavvio normale.

Usa una sola pagina di configurazione alla volta. Al termine salva, disconnetti l'App e **riavvia senza tenere premuti pulsanti** per tornare al funzionamento normale.

<a id="modalita-speciali-allavvio"></a>

### Modalità speciali all'avvio

Premi i pulsanti della lightgun **prima di accendere o riavviare** e mantienili premuti per circa **2 secondi**, finché parte la modalità richiesta. Queste scorciatoie richiedono una lightgun su cui sia già presente il firmware 7.0.0 e pulsanti correttamente mappati e funzionanti.

| Pulsanti tenuti premuti all'avvio | Risultato |
| --- | --- |
| Nessuno | Funzionamento normale. |
| **Grilletto + A** | Modalità aggiornamento firmware; sull'OLED, se presente, compare **Ready for firmware update**. Collega la porta USB OTG della lightgun stessa e usa il [Web Flasher](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/). |
| **B** | Avvia la WebApp offline integrata tramite Wi-Fi e, con il cavo USB, tramite USB NCM. Sull'OLED la barra superiore ha i colori invertiti e mostra un ingranaggio. |

L'aggiornamento firmware ha la precedenza se vengono tenute premute entrambe le combinazioni. **A/B sono i pulsanti mappati della lightgun, non il pulsante fisico BOOT della scheda.** Su una scheda vuota, con un firmware precedente o quando le scorciatoie non funzionano, usa invece la procedura BOOT/RESET della scheda. Su una lightgun che esegue normalmente la 7.0.0, il Web Flasher può anche portarla da solo in modalità flashing, e la WebApp può fare lo stesso con **Riavvia il microcontrollore in modalità aggiornamento firmware** nella scheda *Gun Tests*. Dopo questo riavvio software la porta seriale cambia: riavvia l'installazione e seleziona la nuova porta. Non si può installare il firmware attraverso il dongle wireless o la pagina di configurazione Wi-Fi.

Per l'installazione scegli l'immagine corrispondente alla **scheda e alla variante flash/PSRAM**, non alla telecamera. Aggiornamento normale e installazione pulita usano la stessa immagine; l'installazione pulita elimina tutte le impostazioni e calibrazioni. Passando dalla 6.2.1 alla 7.0.0 è consigliata un'installazione pulita: annota prima le impostazioni esistenti, poi configura e calibra nuovamente.

### Telecamera, display e impostazioni di avvio

- **Telecamera:** in **Impostazioni Gun → Modello TELECAMERA**, scegli **DFRobot SEN0158 / Wii**, **PAJ7025R2** o **PAJ7025R3** secondo l'hardware installato (abilita **Sblocca modifica TELECAMERA** se i controlli sono bloccati). Imposta anche i pin di comunicazione nel **Layout Scheda**: SDA/SCL per I2C, oppure RX (MISO), TX (MOSI), SCK e CSn per SPI. Sono collegamenti di segnale, oltre ad alimentazione e massa. Salva per applicare la selezione; ricalibra dopo aver cambiato telecamera, lente o posizione degli emettitori. La lente più ampia della R3 aiuta a distanze ridotte, ma non garantisce maggiore precisione o portata. Dopo un'installazione pulita la telecamera selezionata è **DFRobot SEN0158 / Wii**: con una PAJ7025R2 o R3, selezionala e imposta i suoi pin SPI prima di calibrare, altrimenti la WebApp può segnalare **Errore dispositivo: Fotocamera non disponibile!**
- **Display OLED:** il display è disabilitato in modo predefinito. Attiva **Abilita Display OLED (SSD1306 128x64)** in **Impostazioni Gun → Periferiche I2C** e, nel **Layout Scheda**, assegna **SDA periferiche** e **SCL periferiche** ai pin a cui è collegato il display; poi salva e riavvia. Se il display resta spento, prova **Usa Indirizzo Dispositivo Alternativo**.
- **Modalità all'avvio:** scegli **Mouse assoluto** (predefinito), **Gamepad (stick destro)** o **Gamepad (stick sinistro)** in **Impostazioni Gun → Input / Output**, poi salva. La scelta vale dal successivo avvio normale; i comandi seriali del gioco possono comunque cambiarla durante l'uso.
- **Pedale wireless** (**Impostazioni Gun → Input / Output**): abilita questa opzione solo se usi l'accessorio, lasciando non assegnati entrambi gli ingressi dei pedali cablati. Altrimenti disabilitala per evitare una ricerca inutile del pedale all'avvio. Un pedale cablato mappato ha la precedenza. Il pedale wireless funziona solo quando la pistola è collegata tramite il dongle; con la pistola collegata al computer via cavo USB, usa un pedale cablato.

<a id="prima-configurazione-italiano"></a>

## Prima Configurazione
Dopo l'installazione su una scheda nuova, o dopo un'installazione pulita, la pistola non ha ancora una calibrazione. Resta in attesa della prima calibrazione e non muove il puntatore; il LED RGB della scheda, se presente, lampeggia di arancione.

1. Collega la pistola alla WebApp. Se usi pin personalizzati, assegna almeno il *Grilletto* e il *Pulsante A* nel *Layout Scheda*; scegli la telecamera e prova i pulsanti nella scheda *Gun Tests*, poi salva.
2. Avvia la calibrazione con uno dei pulsanti *Calibra profilo* nella scheda *Profili di calibrazione*, oppure premi il grilletto per la calibrazione autonoma (vedi [Come Calibrare](#come-calibrare-italiano)).
3. Salva la calibrazione.

Se la telecamera non è disponibile (ad esempio perché è selezionato il modello sbagliato), premendo il grilletto la calibrazione non parte: correggi prima le impostazioni della telecamera nella WebApp.

<a id="manuale-operativo-italiano"></a>

## Manuale Operativo
Per impostazione predefinita, la lightgun funziona come un mouse a posizionamento assoluto (come un pennino per tablet!) finché non viene premuto il pulsante/combinazione per entrare in modalità pausa. Si può scegliere un’uscita diversa salvando l’impostazione **Modalità all’avvio**. In alternativa, si può istruire la pistola a inviare output tramite il suo corrispondente dispositivo HID Gamepad utilizzando un programma di distribuzione di feedback seriale come MAMEHOOKER - vedi la sezione [Modalità Serial Handoff](#modalità-serial-handoff-mame-hooker-italiano) per maggiori informazioni.

**Come il computer vede la pistola:** la lightgun, o il suo dongle, viene riconosciuta come normali dispositivi USB (un mouse a posizionamento assoluto, una tastiera e un gamepad): non servono driver e giochi ed emulatori possono usarla come un qualsiasi mouse o controller. Per far seguire a solenoide, rumble, LED e contatori dell'OLED ciò che accade nel gioco serve anche un programma di force feedback come MAMEHOOKER (vedi [Modalità Serial Handoff](#modalità-serial-handoff-mame-hooker-italiano)); senza, il force feedback reagisce solo al grilletto.

Qualsiasi terminale seriale (Monitor Seriale dell'IDE Arduino, *PuTTY*, *screen*, ecc.) può essere utilizzato per visualizzare le informazioni mentre la pistola è in pausa e durante la calibrazione autonoma.

Nota che i pulsanti in modalità pausa (e la combinazione per entrare in modalità pausa) si attivano quando viene rilasciato *l'ultimo pulsante* della combinazione. Questo sistema viene utilizzato per rilevare e differenziare le combinazioni di tasti dalla singola pressione.

* Nota: Al suo picco, la posizione del mouse si aggiorna a 209Hz (sia via cavo USB che via wireless), o all'incirca ogni ~4.8ms, risultando estremamente reattiva.

<a id="modalità-di-funzionamento-italiano"></a>

### Modalità di Funzionamento
La pistola ha le seguenti modalità operative:
1. **Normal** - La posizione del mouse si aggiorna a ogni frame ricevuto dalla telecamera di posizionamento IR (nessuna mediazione).
2. **Averaging** - La posizione è calcolata tramite una media mobile su 2 frame (posizione attuale + precedente).
3. **Averaging2** - La posizione è calcolata tramite una media ponderata del frame attuale e dei 2 frame precedenti.
4. **Processing** - Modalità di test da utilizzare con la WebApp (non è possibile assegnare questa modalità a un profilo).

Le modalità *Averaging* sono sottili ma riducono un po' il jitter (tremolio) del movimento senza aggiungere lag percettibile.

> Il porting di OpenFIRE per ESP32 aggiunge in automatico algoritmi avanzati anti jitter (tremolio), quindi si consiglia di impostare la modalità **Normal**.

<a id="pulsanti-predefiniti-italiano"></a>

### Pulsanti Predefiniti
- **Grilletto (Trigger):** Tasto sinistro del mouse
- **A:** Tasto destro del mouse (In modalità "low buttons", agisce da Start se premuto puntando fuori dallo schermo)
- **B:** Tasto centrale del mouse (In modalità "low buttons", agisce da Select se premuto puntando fuori dallo schermo)
- **C/Reload:** Tasto mouse 4 / Pulsante laterale 1 / Indietro
- **Pump Action (Ricarica a pompa, es. Cabela's):** Tasto destro del mouse
- **Start:** Tasto 1 della tastiera
- **Select:** Tasto 5 della tastiera
- **Su/Giù/Sinistra/Destra:** Frecce direzionali della tastiera
- **Pedale Principale:** Tasto mouse 4 / Pulsante laterale 1 / Indietro
- **Pedale Secondario (Alt):** Tasto mouse 5 / Pulsante laterale 2 / Avanti
- **C + Start:** Tasto Esc della tastiera

La **Modalità pulsanti 'fuori schermo'** (**Impostazioni Gun → Input / Output**; in inglese *Low Buttons Mode*) è pensata per pistole con pochi pulsanti: quando è attiva, A e B premuti puntando fuori dallo schermo funzionano come Start e Select.

Si può entrare in modalità Pausa premendo **C + Select** (impostazione predefinita), premendo il tasto **Home** (se presente nel layout dei pin), oppure **tenendo premuto il grilletto insieme al pulsante A senza alcun punto IR in vista** se l'opzione *hold-to-pause* è abilitata - in quest'ultimo caso si consiglia di puntare la pistola verso il pavimento.

**Quale modalità di pausa è attiva?** In modo predefinito la pistola usa la modalità di pausa **Hotkey**: una volta in pausa, ogni pulsante o combinazione elencata sotto esegue direttamente la sua funzione. Per usare invece il **Menu di Pausa Semplificato** (un elenco di opzioni da scorrere con i pulsanti, più comodo con un display OLED), attiva **Menu pausa semplice** in **Impostazioni Gun → UI e UX**. Per entrare in pausa tenendo premuti Grilletto + A serve **Abilita pausa con pressione prolungata** nello stesso gruppo, dove si imposta anche il tempo di pressione. Salva dopo aver cambiato queste opzioni.

<a id="pulsanti-predefiniti-in-modalità-pausa-hotkey-italiano"></a>

### Pulsanti Predefiniti in Modalità Pausa (Hotkey)
- **A, B, Start, Select:** Seleziona un profilo.
- **Start + A:** Modalità pistola Normal (Averaging disabilitato).
- **Start + B:** Modalità pistola Normal con Averaging, passa da una modalità di media all'altra (usa il monitor seriale per vedere l'impostazione).
- **B + Giù:** Diminuisci la sensibilità della telecamera IR (usa il monitor seriale per vedere l'impostazione).
- **B + Su:** Aumenta la sensibilità della telecamera IR (usa il monitor seriale per vedere l'impostazione).
- **C/Reload:** Esci dalla modalità pausa.
- **Sinistra:** Attiva/Disattiva Rumble *(se non viene rilevato alcuno switch fisico per il rumble)*.
- **Destra:** Attiva/Disattiva Solenoide *(se non viene rilevato alcuno switch fisico per il solenoide)*.
- **Grilletto:** Inizia la calibrazione.
- **Start + Select:** Salva le impostazioni nello spazio di archiviazione flash non volatile.

<a id="controlli-per-il-menu-di-pausa-semplificato-italiano"></a>

#### Controlli per il Menu di Pausa Semplificato
- **A o Su:** Muovi il cursore Su
- **B o Giù:** Muovi il cursore Giù
- **Grilletto:** Seleziona l'opzione
- **C:** Esci dal menu di pausa
  - *Tenendo premuto A o B per metà della durata del tempo di hold-to-pause (circa ~2 secondi di default) si uscirà anche dal menu di pausa semplice.*
  
Le opzioni disponibili nel menu di pausa semplificato sono le seguenti (dalla prima all'ultima, per poi ricominciare):
* Calibra il profilo corrente (sempre la prima opzione iniziale)
* Cambia profilo (sottomenu)
  * Scegli tra i profili 1-4 usando i pulsanti di navigazione/grilletto per selezionare, o premi C per tornare indietro.
* Salva le impostazioni nella memoria non volatile
* Attiva/Disattiva Rumble *(quando abilitato e senza switch fisico)*
* Attiva/Disattiva Solenoide *(quando abilitato e senza switch fisico)*
* Invia segnale del tasto Esc al PC

<a id="come-calibrare-italiano"></a>

### Come Calibrare
*Queste istruzioni si applicano al processo di calibrazione integrato e autonomo della pistola; le schermate di Calibrazione nella WebApp OpenFIRE ESP32 hanno una procedura simile con maggiori informazioni a schermo per guidare l'utente.*

**Prima di iniziare:** posizionati centralmente di fronte allo schermo, tieni la pistola senza ruotarla attorno alla canna e mira con cura ogni bersaglio. Nella WebApp e nell'App desktop questo promemoria resta sopra le istruzioni inferiori per tutta la calibrazione, compresa la verifica finale del puntamento.

1. Seleziona il profilo da calibrare (tramite A/B/Start/Select nella modalità Hotkey Pause, o selezionandolo nel Menu Pausa Semplificato) e premi il grilletto per iniziare la calibrazione. In alternativa, puoi avviarla dalla WebApp (scheda *Profili di calibrazione*, pulsanti *Calibra profilo*).
2. Mira al centro dello schermo e premi il grilletto mantenendo una mira stabile.
3. Il cursore si sposterà sui quattro bordi dello schermo: in alto, in basso, a sinistra e a destra. Su ogni bordo descrive un piccolo cerchio tangente al bordo. Spara nel punto in cui il cerchio tocca il bordo dello schermo, **non al centro di quel cerchio**.
4. Quando il cursore torna al centro dello schermo, mira al centro del piccolo cerchio e spara per terminare la calibrazione.
5. Il nuovo profilo di calibrazione verrà applicato e potrai testare il tracciamento. Un buon indicatore di una calibrazione corretta è il mantenimento di una precisione il più vicino possibile alla linea di vista (line-of-sight) quando si mira ai bordi e agli angoli dello schermo.
   - Se la calibrazione è buona, premi il grilletto per confermare.
   - Se desideri ricominciare la calibrazione, premi il pulsante A o B nella schermata di verifica per ripartire dal punto centrale (Passo 2).
   - La calibrazione può essere annullata del tutto premendo C/Reload in qualsiasi momento, o i pulsanti A/B in qualsiasi momento prima della verifica finale (durante la prima calibrazione dopo un'installazione pulita, questi pulsanti non annullano la procedura; dalla WebApp o dall'App desktop puoi comunque annullarla chiudendo la finestra di calibrazione).

**Cursore in movimento nella calibrazione dalla lightgun.** Di default, il puntatore del mouse descrive continuamente un piccolo cerchio mentre attende ogni tiro. Al centro dello schermo, mira al centro del cerchio; sui quattro bordi, mira al punto di contatto del cerchio con il bordo. L'animazione si interrompe durante lo spostamento verso il bersaglio successivo e nella verifica finale del puntamento. Non cambia i calcoli di calibrazione né l'interfaccia di calibrazione nella WebApp o nell'App desktop.

**Guida alla calibrazione sul display OLED.** Se è installato e abilitato un display SSD1306, la barra di stato superiore resta visibile e lo spazio sottostante guida la calibrazione. Le istruzioni sul display sono sempre in inglese.

- **Dove mirare:** una TV stilizzata e un mirino indicano la posizione del bersaglio sul monitor, non il puntamento attuale della pistola. L'avanzamento va da **1/6** a **6/6**, nell'ordine: `CENTER`, `TOP`, `BOTTOM`, `LEFT`, `RIGHT`, `CENTER` (centro, alto, basso, sinistra, destra, centro finale). `SHOOT` indica di sparare al bersaglio.
- **Verifica finale:** dopo il sesto bersaglio, `CHECK AIM` invita a verificare il puntamento sul monitor. `TRIGGER: CONFIRM` indica di premere il grilletto per confermare; `A/B: REPEAT` indica di premere A o B per ripetere la calibrazione.
- **Punti IR:** i punti bianchi restano visibili sopra le istruzioni durante la calibrazione e la verifica. Mostrano sia le posizioni degli emettitori rilevati sia quelle ricostruite, senza distinguerle: quattro punti non significano necessariamente che la telecamera veda realmente tutti e quattro gli emettitori.

Sono solo aiuti visivi: la sequenza di calibrazione, i pulsanti e la procedura di salvataggio descritti in questa guida restano invariati.

**Calibrazione dalla WebApp o dall'App desktop.** Con il firmware 7.0.0 la finestra di calibrazione avviata dalla WebApp (o dall'App desktop compatibile) controlla anche gli emettitori IR prima di ogni tiro sui bersagli:
- **Colore del mirino:** rosso quando la telecamera non vede tutti e quattro gli emettitori; altrimenti segue la luminosità dell'emettitore più debole su una scala continua, dal rosso all'arancione, poi a un giallo-verde chiaro fino al verde pieno, senza salti.
- **Mirino evidenziato:** un cerchio esterno tratteggiato sottile ruota lentamente con lo stesso colore dinamico del mirino, senza pulsazioni o variazioni di luminosità. Il mirino originale resta fermo e mantiene dimensione e posizione, anche ai bordi dello schermo. Il cerchio esterno colorato non compare durante la verifica del puntamento e resta fermo nella WebApp se il sistema chiede animazioni ridotte.
- **Avanzamento:** sei pallini sopra il titolo del passo rappresentano i bersagli: centro, alto, basso, sinistra, destra, centro finale. Quelli completati sono pieni e grigi, quello corrente è verde e gli altri hanno il contorno grigio. Questo verde indica l’avanzamento, non la qualità IR. I titoli vanno da «Calibrazione: passo 1 di 6» a «Calibrazione: passo 6 di 6», includendo il tiro iniziale al centro; tutti e sei i pallini sono pieni durante «Verifica del puntamento:».
- **Tiro accettato:** un cerchio bianco intero e fermo produce un breve flash intorno al bersaglio appena acquisito (circa 120 ms). Il mirino resta lì fino alla fine del flash, poi si sposta sul bersaglio successivo. Anche il centro finale riceve questa conferma prima che il mirino segua il cursore nella verifica del puntamento. I pallini di avanzamento non lampeggiano. Con animazioni ridotte, la WebApp mostra il cerchio bianco senza sfumarlo, per la stessa durata. Tiri rifiutati e reset non lo attivano. Conferma l’acquisizione del bersaglio, non il salvataggio della calibrazione: verifica la mira, conferma con il grilletto e poi salva le impostazioni come descritto sotto.
- **Riquadro degli emettitori:** un quadrato con angoli arrotondati in basso a destra, con la scritta «LED IR» al centro, che mostra i quattro emettitori nella loro disposizione (Square o Diamond). Dimensione e riempimento di ogni cerchio seguono la telecamera come nel [test IR](#modalità-di-test-italiano); ogni emettitore visto ha la stessa scala di colori del mirino, secondo la propria luminosità, e il mirino segue il più debole. Un emettitore non visto è un cerchio tratteggiato rosso con una X rossa.
- **Tiri rifiutati:** un tiro sui bersagli (centro, i quattro bordi e il centro finale) viene rifiutato solo quando la telecamera non vede tutti e quattro gli emettitori (mirino rosso). In quel caso la calibrazione non va avanti e la finestra spiega il motivo: controlla gli emettitori, la distanza e la [sensibilità della telecamera](#sensibilità-della-telecamera-ir-italiano), poi spara di nuovo. La conferma finale nella fase di verifica non viene controllata. Un emettitore debole (mirino arancione) non blocca la calibrazione, ma il risultato più affidabile si ottiene con il mirino verde.

La **legenda in basso a sinistra**, visibile per tutta la procedura compresa la verifica del puntamento, spiega l'intensità del segnale IR (Debole → Forte), la dimensione del punto IR (Piccolo → Grande) e il cerchio rosso tratteggiato con una X (IR non rilevato). Il colore del mirino indica l'emettitore con il segnale più debole; per acquisire ogni bersaglio devono essere rilevati tutti e quattro gli emettitori. La legenda usa un carattere normale e non cambia il funzionamento della calibrazione.

La calibrazione avviata dalla pistola stessa (modalità pausa) funziona come prima, senza questi controlli.

<p align="center">
  <img src="../docs/img/webapp_calibration_ir.png" alt="Calibrazione nella WebApp: mirino rosso, riquadro con i quattro emettitori di cui uno non visto e il messaggio di tiro rifiutato perché la telecamera non vede tutti e quattro gli emettitori IR" width="70%">
</p>

Ricordati di **salvare la calibrazione** e il profilo corrente subito dopo, salvando dalla WebApp, premendo *Start+Select* nella modalità Hotkey Pause, o scegliendo la terza opzione "Save Settings" nel Menu di Pausa Semplificato.

<a id="sensibilità-della-telecamera-ir-italiano"></a>

### Sensibilità della Telecamera IR
La sensibilità della telecamera IR può essere regolata. Si consiglia di impostarla il più in alto possibile. Se la sensibilità IR è troppo bassa, la precisione del puntatore ne risentirà. Tuttavia, una sensibilità troppo elevata potrebbe far sì che la telecamera rilevi riflessi indesiderati, causando salti improvvisi del puntatore. È impossibile sapere a priori quale impostazione funzionerà meglio, poiché dipende dalle specifiche del tuo setup (luminosità degli emettitori IR, distanza, lente della telecamera ed eventuali superfici lucide che causano riflessi).

Un segno che la sensibilità IR è **troppo bassa** si verifica quando il puntatore si muove in modo "scattoso" e poco fluido, come se avesse una bassa risoluzione. Se noti questo problema nonostante la sensibilità sia impostata al massimo, è probabile che i tuoi emettitori IR non siano abbastanza luminosi.

Un segno che la sensibilità IR è **troppo alta** si verifica quando il puntatore salta in modo irregolare ed erratico. Se ciò accade solo mentre miri a determinate aree dello schermo, è un chiaro indicatore che la telecamera sta rilevando un riflesso. Se la sensibilità è al massimo, riducila su alto o minimo. Ovviamente, la soluzione migliore rimane l'eliminazione della superficie riflettente. Il test IR della WebApp (vedi [Modalità di Test](#modalità-di-test-italiano)) può aiutare a diagnosticare questo problema: mostra i punti IR visti dalla telecamera, indicando anche quanto è grande e luminoso ciascuno.

<a id="profili-italiano"></a>

### Profili
Le build principali di OpenFIRE sono configurate con 4 profili di calibrazione disponibili. Ogni profilo ha i propri dati di calibrazione, la modalità operativa (run mode) e le impostazioni di sensibilità della telecamera IR. Ogni profilo può essere richiamato dalla modalità pausa premendo il pulsante associato (A/B/Start/Select) o selezionandolo dal sottomenu dei profili nel menu di pausa semplificato. Nella WebApp la scheda **Profili di calibrazione** mostra e modifica le impostazioni di ogni profilo, comprese **Sensibilità** e **Modalità**: è il modo più semplice per controllarle, perché in modalità pausa vengono indicate solo su un monitor seriale. Salva dopo averle cambiate.

<a id="interruttori-software-toggle-italiano"></a>

### Interruttori Software (Toggle)
Le funzioni hardware possono essere attivate o disattivate in tempo reale, anche se non hai cablato degli interruttori fisici dedicati!

Mentre sei in modalità pausa, i controlli toggle sono i seguenti (il colore indica come si illumina il LED integrato sulla scheda):
- **D-Pad Sinistra: Rumble Toggle** (Salmone) - Abilita/disabilita la funzione rumble (vibrazione). Quando abilitato, il motore si attiverà per un breve periodo come conferma.
- **D-Pad Destra: Solenoid Toggle** (Giallo) - Abilita/disabilita il force feedback del solenoide. Quando abilitato, il solenoide scatterà per un breve periodo come conferma.

Queste operazioni possono essere eseguite anche dalle rispettive opzioni nel Menu di Pausa Semplificato. Lo stato corrente di queste impostazioni viene salvato nella memoria flash al momento del salvataggio manuale e ricaricato all'avvio.

<a id="salvataggio-delle-impostazioni-nella-flash-italiano"></a>

#### Salvataggio delle Impostazioni nella Flash
I dati di calibrazione, le impostazioni dei profili e le opzioni estese della pistola (come la mappatura personalizzata dei pin e l'intensità del rumble) possono essere salvati nella memoria non volatile premendo **Start + Select** nella modalità Hotkey Pause (oppure usando Save Settings nel Menu di Pausa Semplificato), o salvando e ricevendo conferma nella WebApp. Il profilo di calibrazione attualmente selezionato al momento del salvataggio viene impostato come predefinito per le successive accensioni. Le impostazioni generali della pistola (mappatura pin, force feedback, ecc.) si applicano a *tutti i profili*. Nella WebApp e nell'App desktop il pulsante **Salva e invia impostazioni** pulsa finché ci sono modifiche non salvate, per ricordare di salvare.

<a id="modalità-di-test-italiano"></a>

#### Modalità di Test
La Modalità di Test mostra i punti IR così come li vede la telecamera. Aprila dalla scheda *Gun Tests* della WebApp con **Apri Tester telecamera IR...**. È estremamente utile per allineare la telecamera durante la costruzione della lightgun, per verificare che tracci correttamente tutti e 4 i punti e per individuare eventuali riflessi. La validità della forma dei punti di test (un rettangolo nel layout Square, un diamante nel layout Diamond) dipende dal profilo attualmente in uso e dalla sua impostazione del layout IR. Salva le modifiche al layout IR o alla sensibilità della telecamera prima di aprire il test: il test usa le impostazioni memorizzate nella pistola.

Con il firmware 7.0.0 ogni cerchio degli emettitori mostra anche ciò che la telecamera misura:
- **Dimensione:** segue la grandezza della macchia di luce vista dalla telecamera. Allontanandosi dallo schermo di solito si riduce solo di poco, perché la grandezza della macchia dipende soprattutto dalla luminosità del LED e dall'ottica della telecamera.
- **Riempimento:** più intenso al centro (luminosità di picco della macchia) e sfumato verso il bordo (luminosità media). Un cerchio pieno indica un LED luminoso e ben rilevato; un cerchio quasi vuoto indica un LED debole, vicino al limite di rilevamento.
- **Cerchio tratteggiato con una X rossa:** quell'emettitore non è visto dalla telecamera; la sua posizione è solo stimata.

Il cerchio grigio indica dove sta puntando la pistola e il cerchio rosso il centro dei quattro emettitori. Con la telecamera DFRobot/Wii dimensione e luminosità arrivano dal formato dati completo della telecamera e sono convertite sulla stessa scala. Con firmware precedenti vengono mostrati i cerchi classici.

<p align="center">
  <img src="../docs/img/webapp_ir_test.png" alt="Test della telecamera IR nella WebApp: tre emettitori visti, con cerchi di dimensione e luminosità diverse, e un emettitore non visto, segnato con una X rossa" width="70%">
</p>

<a id="problemi-comuni-italiano"></a>

## Problemi Comuni

- **Dopo l'installazione la pistola non fa nulla e il puntatore non si muove:** una pistola nuova o appena installata da zero attende la prima calibrazione (il LED RGB della scheda, se presente, lampeggia di arancione). Premi il grilletto per calibrare, oppure calibra dalla WebApp; vedi [Prima Configurazione](#prima-configurazione-italiano).
- **La WebApp non trova la porta della pistola:** usa un cavo USB **dati** (alcuni cavi servono solo per la ricarica) e un browser da computer con Web Serial, come Chrome o Edge. Con i browser senza Web Serial (ad esempio Firefox, Safari o i browser dei telefoni) puoi usare la [WebApp offline](#configurazione-con-la-webapp) contenuta nella pistola; per un browser così vecchio da non aprire nessuna delle due, vedi l'App desktop in [Configurazione della Scheda](#configurazione-della-scheda-italiano). Sulle schede con due connettori USB usa quello collegato all'USB nativa dell'ESP32-S3 (OTG): sulla DevKitC-1, guardandola dall'alto con i connettori rivolti in avanti, è quello a destra, davanti al LED NeoPixel integrato (l'immagine della scheda nella sezione **Layout Scheda** della WebApp lo indica come **USB OTG**). Chiudi le altre pagine della WebApp, i monitor seriali e MAMEHOOKER, che possono tenere occupata la porta. Se la pistola è stata avviata tenendo premuto **B**, la sua porta seriale è sostituita dalla rete della WebApp offline: riavvia senza tenere premuti pulsanti.
- **Messaggio della WebApp "Errore dispositivo: Fotocamera non disponibile!":** dopo un'installazione pulita il modello di telecamera è **DFRobot SEN0158 / Wii**. Con una PAJ7025R2 o R3, selezionala in **Impostazioni Gun → Modello TELECAMERA**, imposta i suoi pin SPI nel **Layout Scheda**, poi salva e riavvia. Altrimenti controlla l'alimentazione della telecamera e che i suoi fili non siano invertiti o assegnati ad altre funzioni.
- **Nessun punto IR nel test, o ne manca qualcuno:** controlla che gli emettitori siano accesi e adatti alla telecamera (940 nm per DFRobot/Wii, 850 nm per PAJ7025R2/R3), che la telecamera veda tutta la disposizione degli emettitori e che la sensibilità della telecamera non sia troppo bassa. Salva la sensibilità prima di aprire il test.
- **La calibrazione non passa al bersaglio successivo:** in una calibrazione avviata dalla WebApp o dall'App desktop un tiro sui bersagli viene rifiutato finché la telecamera non vede tutti e quattro gli emettitori (mirino rosso). Guarda il riquadro degli emettitori in basso a destra: un emettitore segnato con una X rossa deve tornare nel campo visivo della telecamera (allontanati se ne esce ai bordi dello schermo); controlla anche l'emettitore, la distanza e la [sensibilità della telecamera](#sensibilità-della-telecamera-ir-italiano). Un emettitore mostrato in arancione è debole: il tiro viene accettato, ma migliorarlo rende la mira più stabile. Vedi [Come Calibrare](#come-calibrare-italiano).
- **Il puntatore salta o trema:** di solito è un riflesso o un'altra sorgente di luce (finestra, lampada, superficie lucida) vista come un punto IR in più. Cercalo nel [test IR](#modalità-di-test-italiano), eliminalo o coprilo, oppure riduci la [sensibilità della telecamera](#sensibilità-della-telecamera-ir-italiano). Ripeti la calibrazione dopo aver spostato gli emettitori.
- **La pistola non si associa al dongle:** inserisci prima il dongle e attendi circa 15 secondi, poi accendi la pistola senza collegarla al computer via USB (con l'USB funziona via cavo). Dopo aver scollegato o riavviato il dongle, riavvia anche la pistola. Usa un dongle per ogni pistola e firmware della stessa release per lightgun e dongle ([guida dongle](../../dongle/README.md#versione-italiana)).
- **Il pedale wireless non viene trovato:** funziona solo con una pistola collegata tramite il dongle. Abilita **Pedale wireless**, lascia non assegnati entrambi gli ingressi dei pedali cablati e salva; accendi il pedale prima della pistola. Dopo aver riavviato il pedale, riavvia anche la pistola.
- **Il display OLED resta spento:** il display è disabilitato in modo predefinito; abilitalo come spiegato in [Telecamera, display e impostazioni di avvio](#telecamera-display-e-impostazioni-di-avvio).
- **La temperatura letta è troppo alta** (ad esempio circa 90 °C a temperatura ambiente): collega il TMP36 alla scheda con un proprio filo di massa, non condiviso con la telecamera ([consiglio di cablaggio](../README.md#moduli-opzionali)).
- **Il collegamento wireless si interrompe o ha poca portata:** lascia libera una piccola area intorno all'antenna della scheda della pistola, senza fili sopra o a ridosso e senza parti metalliche che la coprano ([consiglio sull'antenna](../README.md#componenti-obbligatori-essenziali)); lo stesso vale per un dongle o un pedale autocostruiti. Usa inoltre fili di alimentazione del solenoide di almeno 24 AWG, perché fili sottili possono causare disconnessioni.

Alcuni comportamenti sono limiti del sistema e non guasti: vedi [Limiti noti](#limiti-noti-italiano).

<a id="limiti-noti-italiano"></a>

## Limiti noti

- Il tracciamento affinato vale per il layout Square; il layout Diamond usa il tracciamento originale di OpenFIRE.
- Il gioco senza fili richiede una scheda ESP32-S3; sulle schede RP2040 il firmware funziona solo via cavo USB.
- Si possono usare insieme fino a quattro pistole; ognuna richiede il proprio dongle, e il pedale wireless funziona solo con una pistola collegata tramite il dongle.
- Il firmware si installa dalla porta USB di ciascun dispositivo, non tramite il dongle o la pagina di configurazione Wi-Fi.
- Nella modalità WebApp offline via USB, la porta seriale USB della pistola, e quindi MAMEHOOKER, non è disponibile fino a un riavvio normale.
- Le App di configurazione per firmware precedenti alla 7.0.0, comprese quelle del progetto OpenFIRE originale, non sono compatibili; usa la WebApp o l'App desktop compatibile (vedi [Configurazione della Scheda](#configurazione-della-scheda-italiano)).

<a id="dettagli-tecnici-e-note-varie-italiano"></a>

## Dettagli Tecnici e Note Varie

<a id="modalità-serial-handoff-mame-hooker-italiano"></a>

### Modalità Serial Handoff (Mame Hooker)

Usa l’avvio normale per avere la porta seriale USB della pistola: nella modalità speciale WebApp offline è sostituita dalla rete USB. Disconnetti l’App di configurazione prima di usare la stessa porta con MAMEHOOKER.
La pistola passerà automaticamente il controllo a un'istanza in esecuzione di Mame Hooker (o app similari collegate in seriale) non appena rileverà un codice di avvio compatibile! Se disponibile, il LED integrato sulla scheda e gli eventuali NeoPixel esterni *non statici* diventeranno di colore bianco a media intensità per segnalare l'attivazione della modalità Serial Handoff (salvo modifiche in-game che impongano colori diversi).

Se non hai familiarità con Mame Hooker, **avrai bisogno dei file `.ini` compatibili per ogni gioco che utilizzi** e **la porta COM della pistola dovrebbe essere impostata per corrispondere al numero del giocatore** (COM1 per il P1, COM2 per il P2, ecc.)! L'assegnazione della porta COM può essere fatta su Windows tramite "Gestione dispositivi", o su Linux tramite le impostazioni di registro di Wine nel prefisso in cui avvii il gioco/Mame Hooker. [Consulta la pagina wiki su MAMEHOOKER per maggiori informazioni!](https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/wiki/MAMEHOOKER_Documentation_IT) Per gli utenti Linux che desiderano utilizzare la pistola con il force feedback nativo degli emulatori (attualmente MAME, Flycast o i relativi port su RetroArch), si consiglia di provare [QMamehook](https://github.com/SeongGino/QMamehook).

<a id="modifica-dell-id-usb-per-pistole-multiple-italiano"></a>

### Più pistole e multigiocatore
Per giocare con due (o più) lightgun OpenFIRE sullo stesso computer, ogni pistola deve presentarsi con un'identità diversa. Se più dispositivi condividono lo stesso nome e/o ID Prodotto/Venditore (PID/VID, gli **identificatori USB Implementer's Forum (USB-IF)**), le applicazioni che leggono individualmente i singoli mouse, come RetroArch e TeknoParrot, non riescono a distinguerli.

1. Collega alla WebApp una pistola alla volta e apri **Impostazioni Gun → Identificatore TinyUSB**.
2. Scegli il numero del giocatore con i pulsanti **P1**–**P4**: **P1** per la prima pistola, **P2** per la seconda, e così via. Dopo un'installazione pulita ogni pistola è **P1**. Con **Visualizzazione avanzata** puoi invece inserire un ID prodotto e un nome personalizzati.
3. Salva e riavvia la pistola.

Il numero del giocatore cambia anche i tasti Start e Select (P1: 1 e 5, P2: 2 e 6, P3: 3 e 7, P4: 4 e 8; con un ID personalizzato restano quelli di P1), viene mostrato sui display OLED e dai LED del pedale ed è il numero da usare per la porta COM di MAMEHOOKER (COM1 per P1, COM2 per P2, ecc.).

**Gioco senza fili con più pistole:** si possono usare senza fili fino a quattro pistole, ognuna con **il proprio dongle** e, se vuoi, il proprio pedale wireless; ogni dongle si associa a una sola pistola. Non importa a quale dongle si associ una pistola: il dongle assume l'identità della pistola (nome, numero del giocatore e numero di serie USB), quindi il computer vede sempre P1 come P1. L'associazione e la scelta del canale radio sono automatiche, e i dongle possono trovarsi sullo stesso canale o su canali diversi: non c'è nulla da impostare.

Se vuoi, prepara una pistola alla volta: inserisci il primo dongle, attendi circa 15 secondi, accendi la prima pistola e il suo pedale, poi fai lo stesso con la successiva. Così ogni dongle sceglie il canale valutando anche le eventuali interferenze, comprese quelle prodotte dai dongle e dalle pistole già in uso, e ogni pedale wireless si associa alla pistola giusta, perché un pedale si associa alla prima pistola che lo cerca. È un consiglio, non un obbligo.

---
### Domande o Problemi?
Per supporto tecnico e per unirti alla community, consulta la [Sezione Community e Supporto](../../README.md#community-support-italiano) nella Home del progetto.
