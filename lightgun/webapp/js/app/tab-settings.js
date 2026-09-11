/*  OpenFIRE Web App - "Gun Settings" tab (Qt: settingsTab).

    Texts (titles, labels, descriptions) are those of the Qt form, keyed by the Qt widget name.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    const UI_TEXTS = {
            "forceFeedbackBox": {
                    "title": "Force Feedback"
            },
            "rumbleFFBox": {
                    "title": "Rumble"
            },
            "label_3": {
                    "text": "Rumble Intensity:"
            },
            "rumbleIntensityBox": {
                    "whatsThis": "<html><head/><body><p>This setting determines the intensity of rumble force feedback events. Scales from 0 (disabled) to 255 (maximum).</p><p><span style=\" font-weight:700; font-style:italic;\">Default:</span><span style=\" font-style:italic;\"> 255 (100%)</span></p></body></html>",
                    "accessibleName": "Rumble Intensity Amount",
                    "suffix": " / 255"
            },
            "label_4": {
                    "text": "Rumble Length:"
            },
            "rumbleLengthBox": {
                    "whatsThis": "<html><head/><body><p>When Rumble is enabled, this setting determines the length of how long the rumble motor is activated when pulling the trigger offscreen, in milliseconds.</p><p><span style=\" font-weight:700; font-style:italic;\">Default:</span><span style=\" font-style:italic;\"> 150 (ms)</span></p></body></html>",
                    "accessibleName": "Length of Rumble Events",
                    "suffix": " ms"
            },
            "rumbleFFToggle": {
                    "whatsThis": "<html><head/><body><p>This setting determines if the rumble motor will be used for onscreen trigger feedback; taking the place of the solenoid feedback, rather than using the default off-screen motor feedback. Naturally, this requires <span style=\" font-style:italic;\">Rumble</span> to be enabled.</p><p>Note that <span style=\" font-style:italic;\">Solenoid</span> takes priority over this setting, and will be disabled if this is set.</p></body></html>",
                    "accessibleName": "Rumble-based Force Feedback Toggle",
                    "text": "Use Rumble as primary Force Feedback"
            },
            "rumbleToggle": {
                    "whatsThis": "<html><head/><body><p>Set if the rumble motor is enabled.</p><p>Under normal conditions, the rumble motor will actuate when shooting offscreen <span style=\" font-style:italic;\">(i.e. when the cursor is at the edge of the display)</span>, and can also provide feedback in various sections of the Pause Menu interface.</p></body></html>",
                    "accessibleName": "Rumble Toggle",
                    "text": "Rumble Enabled"
            },
            "autofireToggle": {
                    "whatsThis": "<html><head/><body><p>Enable for rapid solenoid feedback while the trigger is pressed.</p><p>When enabled, the solenoid will always fire off automatically, as if the trigger has been held when using normal solenoid force feedback.</p><p>Note: if <span style=\" font-style:italic;\">Rumble FF</span> is enabled, this option will provide continuous rumble motor feedback instead.</p></body></html>",
                    "accessibleName": "Autofire Toggle",
                    "text": "Autofire Enabled"
            },
            "solenoidFFBox": {
                    "title": "Solenoid"
            },
            "solenoidToggle": {
                    "whatsThis": "<html><head/><body><p>Set if the solenoid is enabled.</p><p>Under normal circumstances, the solenoid actuates when shooting at the display; the &quot;engaged&quot; duration is set in <span style=\" font-style:italic;\">Solenoid Interval</span>. When depressed for the length of <span style=\" font-style:italic;\">Solenoid Hold Length</span>, the solenoid will fire rapidly for as long as the trigger is held, with <span style=\" font-style:italic;\">Solenoid Fast Interval</span> determining state of the &quot;engaged&quot; duration and pause between engagements calculated by <span style=\" font-style:italic;\">Solenoid Fast Interval * Autofire Wait Factor.</span></p><p>Note that the <span style=\" font-style:italic;\">Rumble FF</span> setting conflicts with solenoid functionality, and will be disabled if this is set.</p></body></html>",
                    "accessibleName": "Solenoid Toggle",
                    "text": "Solenoid Enabled"
            },
            "label_14": {
                    "text": "Temperature Warning Threshold:"
            },
            "label_16": {
                    "text": "Temperature Shutoff Threshold:"
            },
            "tempWarningBox": {
                    "whatsThis": "<html><head/><body><p>This determines the point at which solenoid activations will force longer OFF periods between sustained fire/autofire shots, in <span style=\" font-style:italic;\">degrees Celsius.</span> Generally speaking, the higher this is set, the longer the solenoid will maintain full sustained fire speeds.</p><p>Once this point is hit, the solenoid firing rate is reduced until the temperature monitor reads at least five (5)°C <span style=\" font-weight:700;\">below</span> this threshold. The &quot;reduced&quot; firing rate is automatically calculated as <span style=\" font-weight:700; font-style:italic;\">Solenoid ON Fast Length * 4</span><span style=\" font-style:italic;\">.</span></p><p><span style=\" font-weight:700;\">Warning: Setting this value too high may reduce the lifespan of the connected solenoid!</span></p></body></html>",
                    "accessibleName": "Temperature Warning Threshold",
                    "suffix": "°C"
            },
            "tempShutoffBox": {
                    "whatsThis": "<html><head/><body><p>This setting determines the point at which the solenoid will stop activating, in <span style=\" font-style:italic;\">degrees Celsius.</span> Generally speaking, the higher the temperature threshold, the longer the solenoid will activate at reduced sustained fire speeds.</p><p>Once this point is triggered, solenoid activations will not be allowed for either single or sustained shots until the temperature monitor reads at least five (5)°C <span style=\" font-weight:700;\">below</span> this threshold.</p><p><span style=\" font-weight:700;\">Warning: Setting this value too high WILL reduce the lifespan of the connected solenoid!</span></p></body></html>",
                    "accessibleName": "Temperature Shutoff Threshold",
                    "suffix": "°C"
            },
            "label_7": {
                    "text": "Solenoid ON Length:"
            },
            "solenoidOnLengthBox": {
                    "whatsThis": "<html><head/><body><p>When <span style=\" font-style:italic;\">Solenoid</span> is enabled, this setting determines the amount of time the solenoid is engaged when the trigger is held, in milliseconds.</p><p><span style=\" font-weight:700; font-style:italic;\">Default:</span><span style=\" font-style:italic;\"> 45 (ms)</span></p></body></html>",
                    "accessibleName": "Solenoid Engagement Length",
                    "suffix": " ms"
            },
            "solenoidOffLengthBox": {
                    "whatsThis": "<html><head/><body><p>When <span style=\" font-style:italic;\">Solenoid</span> is enabled, this setting determines the amount of time between solenoid activations when the trigger is held in <span style=\" font-style:italic;\">Autofire</span> mode (or after holding the trigger for <span style=\" font-style:italic;\">Solenoid Hold Length),</span> in milliseconds.</p><p><span style=\" font-weight:700; font-style:italic;\">Default:</span><span style=\" font-style:italic;\"> 80 (ms)</span></p></body></html>",
                    "accessibleName": "Solenoid Autofire Wait Length",
                    "suffix": " ms"
            },
            "label_8": {
                    "text": "Solenoid Autofire Wait Time:"
            },
            "label_9": {
                    "text": "Solenoid Hold-to-Autofire Length:"
            },
            "solenoidHoldLengthBox": {
                    "whatsThis": "<html><head/><body><p>When <span style=\" font-style:italic;\">Solenoid</span> is enabled and <span style=\" font-style:italic;\">Autofire</span> is disabled, this setting determines the time it takes to hold the trigger after a <span style=\" font-style:italic;\">single shot</span> solenoid activation before transitioning to a sustained fire feedback, in milliseconds.</p><p><span style=\" font-weight:700; font-style:italic;\">Default:</span><span style=\" font-style:italic;\"> 500 (ms)</span></p></body></html>",
                    "accessibleName": "Solenoid Trigger Hold to Autofire Length",
                    "suffix": " ms"
            },
            "lightingBox": {
                    "title": "Lighting"
            },
            "commonAnodeToggle": {
                    "whatsThis": "<html><head/><body><p>This setting determines if the 4-pin RGB LED module uses a <span style=\" font-style:italic;\">Common Anode</span> or <span style=\" font-style:italic;\">Common Cathode.</span></p><p><span style=\" font-weight:700;\">Enable</span> if your 4-pin LED is a Common Anode type (with common connected to 5V).<br/><span style=\" font-weight:700;\">Disable</span> if your 4-pin LED is a Common Cathode (with common connected to GND).</p></body></html>",
                    "accessibleName": "4-Pin RGB LED Common Anode Toggle",
                    "text": "4-Pin RGB Common Anode"
            },
            "neopixelGroupBox": {
                    "title": "NeoPixels"
            },
            "label": {
                    "text": "NeoPixel Strand Length: "
            },
            "customLEDstaticBtn1": {
                    "whatsThis": "<html><head/><body><p>Click to set the color for this Pixel.</p></body></html>",
                    "accessibleName": "First NeoPixel Static Color",
                    "text": "Static Color 1"
            },
            "customLEDstaticBtn2": {
                    "whatsThis": "<html><head/><body><p>Click to set the color for this Pixel.</p></body></html>",
                    "accessibleName": "Second NeoPixel Static Color",
                    "text": "Static Color 2"
            },
            "customLEDstaticBtn3": {
                    "whatsThis": "<html><head/><body><p>Click to set the color for this Pixel.</p></body></html>",
                    "accessibleName": "Third NeoPixel Static Color",
                    "text": "Static Color 3"
            },
            "neopixelStrandLengthBox": {
                    "whatsThis": "<html><head/><body><p>This setting determines how many LEDs are in the NeoPixel chain. More than one NeoPixel can be daisychained together, starting from the Pixel that's directly connected to <span style=\" font-style:italic;\">NeoPixel Pin.</span></p><p><span style=\" font-weight:700; font-style:italic;\">Default:</span><span style=\" font-style:italic;\"> 1</span></p></body></html>",
                    "accessibleName": "NeoPixel Strand Length",
                    "suffix": " Pixel(s)"
            },
            "customLEDstaticSpinbox": {
                    "whatsThis": "<html><head/><body><p>This setting determines how many LEDs in the NeoPixel chain are statically colored. </p><p>These first (x) Pixels <span style=\" font-weight:700;\">won't be affected</span> by menu elements or Serial LED events. Which colors for which Pixels can be set using the <span style=\" font-style:italic;\">Static Color</span> boxes.</p><p><span style=\" font-weight:700; font-style:italic;\">Default:</span><span style=\" font-style:italic;\"> 0</span></p></body></html>",
                    "accessibleName": "Amount of NeoPixels using Static Colors",
                    "suffix": " Pixel(s)",
                    "prefix": "First "
            },
            "pixelChangeNotice": {
                    "text": "Changes to NeoPixel settings may require a power cycle to update properly!"
            },
            "label_2": {
                    "text": "Static NeoPixels: "
            },
            "invertStaticPixelsBox": {
                    "whatsThis": "<html><head/><body><p>Enabling this setting will have statically illuminated pixels set to the <span style=\" font-style:italic;\">last</span> X Pixel(s), rather than the default <span style=\" font-style:italic;\">first</span> X Pixel(s).</p></body></html>",
                    "accessibleName": "Set Static Pixels to First/Last Pixels Toggle",
                    "text": "Invert Static Colors"
            },
            "inputBox": {
                    "title": "Input"
            },
            "lowButtonsToggle": {
                    "whatsThis": "<html><head/><body><p>This setting determines how <span style=\" font-style:italic;\">Button A</span> &amp; <span style=\" font-style:italic;\">Button B</span> behaves during normal use.</p><p><span style=\" font-weight:700;\">When Enabled, </span><span style=\" font-style:italic;\">Button A</span> &amp; <span style=\" font-style:italic;\">Button B</span> will perform different functions when aiming off-screen, instead actuating the functions of <span style=\" font-style:italic;\">Start</span> &amp; <span style=\" font-style:italic;\">Select</span> respectively.<br/><span style=\" font-weight:700;\">When Disabled,</span> all buttons will perform the same functions, regardless of aiming off-screen or not.</p><p>Enabled is recommended for lightguns <span style=\" font-weight:700;\">with two or less sub buttons.</span></p></body></html>",
                    "accessibleName": "Low Buttons Mode Toggle",
                    "text": "Low Buttons Mode"
            },
            "uiuxBox": {
                    "title": "UI and UX"
            },
            "holdToPauseLengthBox": {
                    "whatsThis": "<html><head/><body><p>When <span style=\" font-style:italic;\">Hold To Pause</span> mode is enabled, this setting determines how long (in milliseconds) should the buttons be held before Pause Mode activates.</p><p><span style=\" font-weight:700; font-style:italic;\">Default:</span><span style=\" font-style:italic;\"> 2500 (ms)</span></p></body></html>",
                    "accessibleName": "Length to activate Hold-to-Pause Menu",
                    "suffix": " ms"
            },
            "label_5": {
                    "text": "Hold-to-Pause Length:"
            },
            "holdToPauseToggle": {
                    "whatsThis": "<html><head/><body><p>This setting determines how to trigger Pause Mode.</p><p><span style=\" font-weight:700;\">When Enabled,</span> Pause Mode can be activated by holding <span style=\" font-style:italic;\">Trigger</span> and <span style=\" font-style:italic;\">Button A </span><span style=\" text-decoration: underline;\">while no IR points are visible.</span><br/><span style=\" font-weight:700;\">When Disabled,</span> Pause Mode can be activated by pressing the Pause Mode hotkey (Button C + Select).</p><p>Enabled is recommended for guns <span style=\" font-weight:700;\">with less than two sub buttons.</span> Note that regardless of setting, Pause Mode can always be accessed by pressing the Home Button (if set).</p></body></html>",
                    "accessibleName": "Pause Mode Hold-to-Activate Toggle",
                    "text": "Hold to Pause Enabled"
            },
            "simplePauseToggle": {
                    "whatsThis": "<html><head/><body><p>This setting determines how Pause Mode functions.</p><p><span style=\" font-weight:700;\">When Enabled,</span> Pause Mode will use a scrolling menu layout - navigate by using <span style=\" font-style:italic;\">Button A/B</span> to select, and press Trigger to activate.<br/><span style=\" font-weight:700;\">When Disabled,</span> Pause Menu options are activated with hotkeys.</p><p>Enabled is recommended for guns <span style=\" font-weight:700;\">with less than two sub buttons,</span> but <span style=\" font-weight:700;\">requires an LED or OLED Display for visual feedback!</span></p></body></html>",
                    "accessibleName": "Simple Pause Menu Toggle",
                    "text": "Simple Pause Menu"
            },
            "i2cGroup": {
                    "title": "I2C Peripherals"
            },
            "i2cOLEDtoggle": {
                    "whatsThis": "<html><head/><body><p>If <span style=\" font-style:italic;\">Peripheral I2C (SDA)</span> and <span style=\" font-style:italic;\">Peripheral I2C (SCL)</span> pins are set, enable use of an SSD1306 OLED display.</p><p>The OLED display can be used to navigate on-gun settings, visualize special modes, and display Ammo &amp; Life counts in compatible software.</p><p><span style=\" font-weight:700;\">NOTE:</span> Because this device does not communicate back to the board, the firmware has no way of knowing whether the device successfully initialized or not. If the display is connected properly, but does <span style=\" font-weight:700;\">not</span> seem to power on on startup/after syncing settings, then try the <span style=\" font-style:italic;\">Use Alternative Device Address</span> setting before attempting to rewire the device.</p></body></html>",
                    "accessibleName": "I2C OLED Display Toggle",
                    "text": "Enable OLED Display (SSD1306 128x64)"
            },
            "oledGroup": {
                    "title": "SSD1306 OLED Display"
            },
            "oledAltAddrsToggle": {
                    "whatsThis": "<html><head/><body><p>If the display doesn't seem to power on when <span style=\" font-style:italic;\">Enable OLED Display</span> is set on startup/after syncing settings to the device (and you have confirmed that it is wired correctly), enabling this will try a different start address.</p><p>The firmware's display component normally defaults to I2C address <span style=\" font-weight:700;\">0x3C</span>, but some 128x64 screens may have device address <span style=\" font-weight:700;\">0x3D</span> instead, which this setting will inform the firmware to try instead.</p></body></html>",
                    "accessibleName": "I2C OLED Device Alternate Address Toggle",
                    "text": "Use Alternative Device Address"
            },
            "cameraModelGroupbox": {
                    "whatsThis": "<html><head/><body><p>Select the camera hardware physically installed in your lightgun.</p><p>The active camera hardware will be swapped automatically when the settings are saved.</p><p>A <span style=\" font-weight:700;\">manual reboot is only recommended</span> if you experience issues.</p><p>Remember to also configure the corresponding communication pins in the Pins Layout tab.</p></body></html>",
                    "accessibleName": "Infrared Camera Model",
                    "title": "CAMERA Model"
            },
            "cameraModelUnlockToggle": {
                    "whatsThis": "<html><head/><body><p>Select the camera hardware physically installed in your lightgun.</p><p>The active camera hardware will be swapped automatically when the settings are saved.</p><p>A <span style=\" font-weight:700;\">manual reboot is only recommended</span> if you experience issues.</p><p>Remember to also configure the corresponding communication pins in the Pins Layout tab.</p></body></html>",
                    "accessibleName": "Infrared Camera Model",
                    "text": "Unlock CAMERA Modification",
                    "toolTip": "Warning: Modifying this is only required once when building the gun."
            },
            "camModel_DFRobot": {
                    "whatsThis": "<html><head/><body><p>This is an <span style=\" font-weight:700;\">I2C</span> camera module.</p><p>When selecting this model, ensure you have correctly mapped the <span style=\" font-style:italic;\">cam_SDA</span> and <span style=\" font-style:italic;\">cam_SCL</span> pins in the Layout tab.</p></body></html>",
                    "accessibleName": "DFRobot SEN0158 / WiiCam",
                    "text": "DFRobot SEN0158 / WiiCam"
            },
            "camModel_PAJ7025R2": {
                    "whatsThis": "<html><head/><body><p>This is an <span style=\" font-weight:700;\">SPI</span> camera module.</p><p>When selecting this model, ensure you have correctly mapped the <span style=\" font-style:italic;\">cam_SPI_MOSI, cam_SPI_MISO, cam_SPI_SCK,</span> and <span style=\" font-style:italic;\">cam_SPI_CS</span> pins in the Layout tab.</p></body></html>",
                    "accessibleName": "PixArt PAJ7025 R2",
                    "text": "PixArt PAJ7025 R2"
            },
            "camModel_PAJ7025R3": {
                    "whatsThis": "<html><head/><body><p>This is an <span style=\" font-weight:700;\">SPI</span> camera module.</p><p>When selecting this model, ensure you have correctly mapped the <span style=\" font-style:italic;\">cam_SPI_MOSI, cam_SPI_MISO, cam_SPI_SCK,</span> and <span style=\" font-style:italic;\">cam_SPI_CS</span> pins in the Layout tab.</p></body></html>",
                    "accessibleName": "PixArt PAJ7025 R3",
                    "text": "PixArt PAJ7025 R3"
            },
            "tinyUSBGroupbox": {
                    "title": "TinyUSB Identifier"
            },
            "label_6": {
                    "text": "For Multiplayer, make sure each gun has its own unique identifier!"
            },
            "label_15": {
                    "text": "<html><head/><body><p>Start/Select key mappings will default to <span style=\" font-style:italic;\">P1</span> binds (<span style=\" font-weight:700;\">Key_1</span> &amp; <span style=\" font-weight:700;\">Key_5</span>) for custom identifiers set here:</p></body></html>"
            },
            "label_11": {
                    "text": "Product ID:"
            },
            "productIdInput": {
                    "whatsThis": "<html><head/><body><p>This determines the custom Product ID that the microcontroller reports when connected to any given device, in base 16 hexadecimal (0-F)</p><p>Custom Product IDs <span style=\" font-weight:700;\">are required</span> to distinguish different lightguns in applications that can individually address multiple mouse devices.</p><p>Note: The USB Vendor ID is hard-coded into the firmware, <span style=\" font-weight:700; font-style:italic;\">0xF143.</span></p></body></html>",
                    "accessibleName": "Custom USB Product ID",
                    "prefix": "0x"
            },
            "label_12": {
                    "text": "Product Name:"
            },
            "productNameInput": {
                    "whatsThis": "<html><head/><body><p>This dictates the custom Product Name that the microcontroller reports when connected to any given device.</p><p>Custom Product Names are helpful (and in some cases, <span style=\" font-weight:700;\">required</span>) to distinguish different lightguns in applications that can individually address multiple mouse devices.</p></body></html>",
                    "accessibleName": "Custom USB Product Name",
                    "placeholderText": "(15 Characters)"
            },
            "label_13": {
                    "text": "Start/Select key mappings will change according to the player identifier set here:"
            },
            "tUSB_p1": {
                    "whatsThis": "<html><head/><body><p>This determines the Product ID and Name that the microcontroller reports when connected to any given device, in decimal (with its hexadecimal equivalent).</p><p>Different Identifiers <span style=\" font-weight:700;\">are required</span> to distinguish different lightguns in applications that can individually address multiple mouse devices.</p><p>When using any <span style=\" font-style:italic;\">Simple USB ID Preset,</span> the <span style=\" font-style:italic;\">Start</span> &amp; <span style=\" font-style:italic;\">Select</span> keyboard mappings (used by default) will be changed to reflect the player slot preset used (e.g. <span style=\" font-style:italic;\">P1</span> will emit <span style=\" font-style:italic;\">1 &amp; 5,</span><span style=\" font-style:italic;\">P2</span> will emit <span style=\" font-style:italic;\">2 &amp; 6,</span> etc.).</p></body></html>",
                    "accessibleName": "Simple USB Identifier Presets",
                    "text": "P1"
            },
            "tinyUSBLayoutToggle": {
                    "whatsThis": "<html><head/><body><p>This determines whether to show/use <span style=\" font-style:italic;\">Simple USB ID Presets</span> or <span style=\" font-style:italic;\">Custom USB ID Settings.</span></p><p>This can be useful for users with either more than four OpenFIRE lightguns, or simply want to personalize each lightgun, but the <span style=\" font-style:italic;\">Simple Presets</span> may be preferable to some.</p><p>All lightguns default to a <span style=\" font-style:italic;\">P1 Simple USB ID Preset.</span></p></body></html>",
                    "accessibleName": "Toggle Advanced USB ID Settings",
                    "text": "Advanced View"
            },
            "settingsDescText": {
                    "whatsThis": "<html><head/><body><p>This tab has a variety of settings to tweak about the lightgun's force feedback, lighting effects, and how it presents itself to the device.</p><p>Hover over an option to view detailed info about it here.</p></body></html>",
                    "text": "<html><head/><body><p>This tab has a variety of settings to tweak about the lightgun's force feedback, lighting effects, and how it presents itself to the device.</p><p>Hover over an option to view detailed info about it here.</p></body></html>"
            },
            "profilesDescText": {
                    "whatsThis": "<html><head/><body><p>This tab shows information and settings to tweak about your current calibration profiles.</p><p>The table above represents the amount of profiles the current board can store, including currently selected profile, names shown for each profile when using a compatible <span style=\" font-style:italic;\">I2C Display,</span> and colors used for lighting devices when switching to them from <span style=\" font-style:italic;\">Pause Mode.</span></p><p>Hover over an option to view detailed info about it here.</p></body></html>",
                    "text": "<html><head/><body><p>This tab shows information and settings to tweak about your current calibration profiles.</p><p>The table above represents the amount of profiles the current board can store, including currently selected profile, names shown for each profile when using a compatible <span style=\" font-style:italic;\">I2C Display,</span> and colors used for lighting devices when switching to them from <span style=\" font-style:italic;\">Pause Mode.</span></p><p>Hover over an option to view detailed info about it here.</p></body></html>"
            }
    };

    function build(app) {
        const { el, t, spin, checkbox, radio } = OF.UI;
        const state = app.state;
        const B = state.B;
        const T = state.T;
        const E = state.E;

        const text = (name, key) => {
            const value = (UI_TEXTS[name] || {})[key];
            return value ? t(value) : '';
        };
        const desc = new OF.UI.DescriptionBox(text('settingsDescText', 'whatsThis'));
        const track = (node, name) => desc.track(node, text(name, 'accessibleName'), text(name, 'whatsThis'));

        const controls = {};

        /** Checkbox bound to a boolean setting. */
        const toggle = (name, index) => {
            const node = checkbox(text(name, 'text'), false, (on) => { state.setToggle(index, on); app.refresh(); },
                { title: text(name, 'toolTip') || undefined });
            controls[name] = node;
            return track(node, name);
        };

        /** Label + spin box bound to a numeric setting (range of the Qt spin box). */
        const ranges = OF.AppState.SETTING_RANGES;
        const settingName = (index) => Object.keys(T).find((key) => T[key] === index && ranges[key]);
        const number = (labelName, name, index) => {
            const [min, max] = ranges[settingName(index)];
            const box = spin({
                min,
                max,
                prefix: text(name, 'prefix'),
                suffix: text(name, 'suffix'),
                onChange: (value) => { state.setSetting(index, value); app.refresh(); },
            });
            controls[name] = box;
            track(box.root, name);
            return [el('label', { class: 'field-label', text: text(labelName, 'text') }), box.root];
        };

        const fieldset = (name, ...children) => {
            const node = el('fieldset', { class: 'group' }, el('legend', { text: text(name, 'title') }), ...children);
            controls[name] = node;
            return node;
        };
        const plain = (name, className, ...children) => {
            const node = el('fieldset', { class: 'plain ' + (className || '') }, ...children);
            controls[name] = node;
            return node;
        };

        // ----- Force Feedback ----------------------------------------------------------
        const solenoidTemp = plain('solenoidTempBox', 'grid-4',
            ...number('label_14', 'tempWarningBox', T.tempWarning),
            ...number('label_16', 'tempShutoffBox', T.tempShutdown));

        const solenoidSettings = plain('solenoidSettingsBox', 'grid-4',
            ...number('label_7', 'solenoidOnLengthBox', T.solenoidOnLength),
            ...number('label_8', 'solenoidOffLengthBox', T.solenoidOffLength),
            el('div', { class: 'grid-row-center' }, ...number('label_9', 'solenoidHoldLengthBox', T.solenoidHoldLength)),
            solenoidTemp);

        const rumbleSettings = plain('rumbleSettingsBox', 'grid-4',
            el('div', { class: 'grid-row-center' }, toggle('rumbleFFToggle', B.rumbleFF)),
            ...number('label_3', 'rumbleIntensityBox', T.rumbleStrength),
            ...number('label_4', 'rumbleLengthBox', T.rumbleInterval));

        const forceFeedback = fieldset('forceFeedbackBox',
            el('div', { class: 'row center' }, toggle('autofireToggle', B.autofire)),
            fieldset('solenoidFFBox', el('div', { class: 'row center' }, toggle('solenoidToggle', B.solenoid)), solenoidSettings),
            fieldset('rumbleFFBox', el('div', { class: 'row center' }, toggle('rumbleToggle', B.rumble)), rumbleSettings));

        // ----- Lighting ----------------------------------------------------------------------
        const colorButtons = [T.customLEDcolor1, T.customLEDcolor2, T.customLEDcolor3].map((index, i) => {
            const name = `customLEDstaticBtn${i + 1}`;
            const swatch = el('span', { class: 'color-swatch' });
            const button = el('button', { class: 'color-button', on: { click: async () => {
                const session = app.session;
                const color = await OF.UI.pickColor(text(name, 'accessibleName'), state.setting(index), { scope: 'session' });
                if (session !== app.session || !state.loaded || app.busy) return;
                if (color !== null) { state.setSetting(index, color); app.refresh(); }
            } } }, swatch, el('span', { text: text(name, 'text') }));
            button.swatch = swatch;
            controls[name] = button;
            return track(button, name);
        });

        // Qt grid: strand length on the first row; static count, colours and invert on the second.
        const pixelNotice = el('p', { class: 'notice', text: text('pixelChangeNotice', 'text') });
        const neopixels = fieldset('neopixelGroupBox',
            el('div', { class: 'row center wrap' }, ...number('label', 'neopixelStrandLengthBox', T.customLEDcount)),
            el('div', { class: 'row center wrap' },
                ...number('label_2', 'customLEDstaticSpinbox', T.customLEDstatic),
                ...colorButtons,
                toggle('invertStaticPixelsBox', B.invertStaticPixels)),
            pixelNotice);

        const lighting = fieldset('lightingBox',
            el('div', { class: 'row center' }, toggle('commonAnodeToggle', B.commonAnode)),
            neopixels);

        // ----- Input / UI and UX ------------------------------------------------------------
        const input = fieldset('inputBox', el('div', { class: 'row center' }, toggle('lowButtonsToggle', B.lowButtonsMode)));
        const uiux = fieldset('uiuxBox',
            el('div', { class: 'row center' }, toggle('simplePauseToggle', B.simplePause)),
            el('div', { class: 'row center wrap' }, toggle('holdToPauseToggle', B.holdToPause),
                ...number('label_5', 'holdToPauseLengthBox', T.holdToPauseLength)));

        // ----- I2C peripherals ------------------------------------------------------------------
        const oled = fieldset('oledGroup', el('div', { class: 'row center' }, toggle('oledAltAddrsToggle', B.i2cOLEDaltAddr)));
        const i2c = fieldset('i2cGroup', el('div', { class: 'row center' }, toggle('i2cOLEDtoggle', B.i2cOLED)), oled);

        // ----- Camera model --------------------------------------------------------------------
        let cameraUnlocked = false;
        const unlock = checkbox(text('cameraModelUnlockToggle', 'text'), false, (on) => { cameraUnlocked = on; update(); },
            { title: text('cameraModelUnlockToggle', 'toolTip') });
        track(unlock, 'cameraModelUnlockToggle');
        const cameraNames = ['camModel_DFRobot', 'camModel_PAJ7025R2', 'camModel_PAJ7025R3'];
        const cameraRadios = cameraNames.map((name, model) => track(radio('of-camera-model', text(name, 'text'), false, () => {
            state.setSetting(T.cameraModel, model);
            app.refresh();
        }), name));
        const camera = fieldset('cameraModelGroupbox', el('div', { class: 'row center' }, unlock),
            el('div', { class: 'row center wrap' }, ...cameraRadios));
        track(camera, 'cameraModelGroupbox');

        // ----- TinyUSB identifier ------------------------------------------------------------------
        let usbAdvanced = false;
        const usbToggle = checkbox(text('tinyUSBLayoutToggle', 'text'), false, (on) => { usbAdvanced = on; update(); });
        track(usbToggle, 'tinyUSBLayoutToggle');

        const usbRadios = [1, 2, 3, 4].map((n) => track(radio('of-usb-preset', `P${n}`, false, () => {
            state.setUsbPreset(n);
            app.refresh();
        }), 'tUSB_p1'));
        const usbSimple = el('div', { class: 'usb-view' },
            el('p', { class: 'hint', text: text('label_13', 'text') }),
            el('div', { class: 'row center wrap' }, ...usbRadios));

        const productId = spin({ hex: true, min: 0, max: 0xFFFF, prefix: '0x', onChange: (value) => { state.setUsbId(value); app.refresh(); } });
        track(productId.root, 'productIdInput');
        const productName = el('input', { type: 'text', class: 'text-input', attrs: { maxlength: OF.AppState.NAME_MAX, placeholder: text('productNameInput', 'placeholderText'), spellcheck: 'false' } });
        productName.addEventListener('input', () => {
            if (state.setUsbName(productName.value)) {
                productName.classList.remove('invalid');
            } else {
                productName.value = state.cur.tinyUSB.name;
                productName.classList.add('invalid');
            }
            app.refresh();
        });
        track(productName, 'productNameInput');
        const usbAdvancedView = el('div', { class: 'usb-view' },
            el('div', { class: 'hint rich', html: text('label_15', 'text') }),
            el('div', { class: 'grid-2' },
                el('label', { class: 'field-label', text: text('label_11', 'text') }), productId.root,
                el('label', { class: 'field-label', text: text('label_12', 'text') }), productName));

        const usb = fieldset('tinyUSBGroupbox',
            el('p', { class: 'hint', text: text('label_6', 'text') }),
            usbAdvancedView, usbSimple,
            el('div', { class: 'row end' }, usbToggle));

        const rootNode = el('div', { class: 'tab-body' },
            el('div', { class: 'tab-scroll' },
                el('div', { class: 'tab-content narrow settings' }, forceFeedback, lighting, input, uiux, i2c, camera, usb)),
            desc.root);

        function setChecked(name, index) {
            controls[name].input.checked = state.toggle(index);
        }

        function update() {
            if (!state.loaded) return;
            const settings = state.cur.settings;

            // Pin dependent groups (Qt pinBoxes_currentIndexChanged).
            forceFeedback.disabled = !(state.pinMapped(E.rumblePin) || state.pinMapped(E.solenoidPin));
            controls.rumbleFFBox.disabled = !state.pinMapped(E.rumblePin);
            controls.solenoidFFBox.disabled = !state.pinMapped(E.solenoidPin);
            solenoidTemp.disabled = !state.pinMapped(E.tempPin);
            solenoidTemp.hidden = !app.showUnsafe;
            neopixels.disabled = !state.pinMapped(E.neoPixel);
            OF.UI.setEnabled(controls.commonAnodeToggle, state.pinMapped(E.ledR) && state.pinMapped(E.ledG) && state.pinMapped(E.ledB));
            i2c.disabled = !(state.pinMapped(E.periphSDA) && state.pinMapped(E.periphSCL));

            for (const [name, index] of [['autofireToggle', B.autofire], ['solenoidToggle', B.solenoid], ['rumbleToggle', B.rumble],
                ['rumbleFFToggle', B.rumbleFF], ['commonAnodeToggle', B.commonAnode], ['invertStaticPixelsBox', B.invertStaticPixels],
                ['lowButtonsToggle', B.lowButtonsMode], ['simplePauseToggle', B.simplePause], ['holdToPauseToggle', B.holdToPause],
                ['i2cOLEDtoggle', B.i2cOLED], ['oledAltAddrsToggle', B.i2cOLEDaltAddr]])
                setChecked(name, index);

            rumbleSettings.disabled = !state.toggle(B.rumble);
            solenoidSettings.disabled = !state.toggle(B.solenoid);
            OF.UI.setEnabled(controls.autofireToggle, state.autofireAvailable);
            controls.holdToPauseLengthBox.disabled = !state.toggle(B.holdToPause);
            oled.disabled = !state.toggle(B.i2cOLED);

            for (const [name, index] of [['tempWarningBox', T.tempWarning], ['tempShutoffBox', T.tempShutdown],
                ['solenoidOnLengthBox', T.solenoidOnLength], ['solenoidOffLengthBox', T.solenoidOffLength],
                ['solenoidHoldLengthBox', T.solenoidHoldLength], ['rumbleIntensityBox', T.rumbleStrength],
                ['rumbleLengthBox', T.rumbleInterval], ['neopixelStrandLengthBox', T.customLEDcount],
                ['customLEDstaticSpinbox', T.customLEDstatic], ['holdToPauseLengthBox', T.holdToPauseLength]])
                controls[name].value = settings[index];

            controls.customLEDstaticSpinbox.prefix = t(state.toggle(B.invertStaticPixels) ? 'Last ' : 'First ');
            const staticCount = settings[T.customLEDstatic];
            colorButtons.forEach((button, i) => {
                button.disabled = i >= staticCount;
                button.swatch.style.background = OF.UI.toHex(settings[T.customLEDcolor1 + i]);
            });
            pixelNotice.hidden = !state.pixelsChanged();

            cameraRadios.forEach((node, model) => {
                node.input.checked = settings[T.cameraModel] === model;
                node.input.disabled = !cameraUnlocked;
            });
            unlock.input.checked = cameraUnlocked;

            usbToggle.input.checked = usbAdvanced;
            usbSimple.hidden = usbAdvanced;
            usbAdvancedView.hidden = !usbAdvanced;
            usbRadios.forEach((node, i) => { node.input.checked = state.usbPreset === i + 1; });
            productId.value = state.cur.tinyUSB.id;
            if (root.document.activeElement !== productName) productName.value = state.cur.tinyUSB.name;
        }

        return {
            root: rootNode,
            update,
            onShow() { desc.reset(); },
            /** New board data: camera lock and USB view restart like in the Qt App. */
            onLoad() {
                cameraUnlocked = false;
                usbAdvanced = state.usbPreset === 0;
                productName.classList.remove('invalid');
            },
            saveViewState() { return { cameraUnlocked, usbAdvanced }; },
            restoreViewState(view) { if (view) { cameraUnlocked = view.cameraUnlocked; usbAdvanced = view.usbAdvanced; } },
        };
    }

    OF.Tabs = OF.Tabs || {};
    OF.Tabs.settings = { id: 'settings', label: 'Gun Settings', icon: 'settings', build };
})(typeof globalThis !== 'undefined' ? globalThis : this);
