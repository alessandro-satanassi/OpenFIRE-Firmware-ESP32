const OpenFIREmaps = {
    mouseMap: [
        { name: "Left Click", val: 0b00000001 },
        { name: "Right Click", val: 0b00000010 },
        { name: "Middle Click", val: 0b00000100 },
        { name: "Side Button Back", val: 0b00001000 },
        { name: "Side Button Forward", val: 0b00010000 }
    ],
    gamepadMap: [
        { name: "A Button", val: 0 },
        { name: "B Button", val: 1 },
        { name: "X Button", val: 3 },
        { name: "Y Button", val: 4 },
        { name: "Left Shoulder", val: 6 },
        { name: "Right Shoulder", val: 7 },
        { name: "Left Trigger", val: 8 },
        { name: "Right Trigger", val: 9 },
        { name: "Select Button", val: 10 },
        { name: "Start Button", val: 11 },
        { name: "Left Stick Click", val: 13 },
        { name: "Right Stick Click", val: 14 },
        { name: "D-Pad Up", val: 15 },
        { name: "D-Pad Down", val: 16 },
        { name: "D-Pad Left", val: 17 },
        { name: "D-Pad Right", val: 18 }
    ],
    keyboardMap: [
        { name: "Player-relative Start Key", val: 0xFF },
        { name: "Player-relative Coin Key", val: 0xFE },
        { name: "Up Arrow", val: 0xDA },
        { name: "Down Arrow", val: 0xD9 },
        { name: "Left Arrow", val: 0xD8 },
        { name: "Right Arrow", val: 0xD7 },
        { name: "Enter/Return", val: 0xB0 },
        { name: "Backspace", val: 0xB2 },
        { name: "Escape", val: 0xB1 },
        { name: "Left Ctrl", val: 0x80 },
        { name: "Right Ctrl", val: 0x84 },
        { name: "Left Alt", val: 0x82 },
        { name: "Right Alt", val: 0x86 },
        { name: "Left Shift", val: 0x81 },
        { name: "Right Shift", val: 0x85 },
        { name: "Tab", val: 0xB3 },
        { name: "A", val: 97 }, { name: "B", val: 98 }, { name: "C", val: 99 }, { name: "D", val: 100 },
        { name: "E", val: 101 }, { name: "F", val: 102 }, { name: "G", val: 103 }, { name: "H", val: 104 },
        { name: "I", val: 105 }, { name: "J", val: 106 }, { name: "K", val: 107 }, { name: "L", val: 108 },
        { name: "M", val: 109 }, { name: "N", val: 110 }, { name: "O", val: 111 }, { name: "P", val: 112 },
        { name: "Q", val: 113 }, { name: "R", val: 114 }, { name: "S", val: 115 }, { name: "T", val: 116 },
        { name: "U", val: 117 }, { name: "V", val: 118 }, { name: "W", val: 119 }, { name: "X", val: 120 },
        { name: "Y", val: 121 }, { name: "Z", val: 122 },
        { name: "F1", val: 0xC2 }, { name: "F2", val: 0xC3 }, { name: "F3", val: 0xC4 }, { name: "F4", val: 0xC5 },
        { name: "F5", val: 0xC6 }, { name: "F6", val: 0xC7 }, { name: "F7", val: 0xC8 }, { name: "F8", val: 0xC9 },
        { name: "F9", val: 0xCA }, { name: "F10", val: 0xCB }, { name: "F11", val: 0xCC }, { name: "F12", val: 0xCD }
    ],
    funcTypes: ["Mouse", "Keyboard", "Gamepad"]
};


function drawBoardUI(boardName, gunConfig) {
    const boxPositions = OpenFIREshared.boardsBoxPositions[boardName];
    if (!boxPositions) return;

    const boardTitle = document.getElementById("board-title");
    boardTitle.innerHTML = `<span style="color:#aaa;">${gunConfig.currentProfile} | </span>` + (OpenFIREshared.boardNames ? (OpenFIREshared.boardNames[boardName] || boardName) : boardName);

    const pinsLeft = document.getElementById("pins-left");
    const pinsRight = document.getElementById("pins-right");
    const pinsMiddle = document.getElementById("pins-middle");
    pinsLeft.innerHTML = "";
    pinsRight.innerHTML = "";
    pinsMiddle.innerHTML = "";

    const posLeft = OpenFIREshared.boardBoxPositions_e.posLeft;
    const posRight = OpenFIREshared.boardBoxPositions_e.posRight;
    const posMiddle = OpenFIREshared.boardBoxPositions_e.posMiddle;
    const posCheck = OpenFIREshared.boardBoxPositions_e.posCheck;

    const allFunctions = Object.keys(gunConfig.pins).sort();
    const isCustomPins = gunConfig.toggles["CustomPins"] === true;
    const chkCustom = document.getElementById('chk-custom-pins');
    if (chkCustom) { chkCustom.checked = isCustomPins; }


    const elementsLeft = [];
    const elementsRight = [];
    const elementsMiddle = [];

    

    for (let gpio = 0; gpio < boxPositions.length; gpio++) {
        const val = boxPositions[gpio];
        if (val === 0) continue; 

        const group = val & posCheck;
        const order = val ^ group;

        let currentFunc = "-1";
        for (const funcName of allFunctions) {
            if (gunConfig.pins[funcName] === gpio) {
                currentFunc = funcName;
                break;
            }
        }

        let labelColor = "#aaaaaa";
        let capText = "GPIO" + gpio;
        if (capabilitiesMap) {
            const cap = capabilitiesMap[gpio];
            const OF_Const = OpenFIREshared.pinCapabilities_e;
            if (cap & OF_Const.pinAnyI2C) {
                labelColor = "#BE00B0";
                capText += " - I2C(*)";
            } else if (cap & OF_Const.pinCanI2C) {
                if (cap & OF_Const.pinIsI2C1) {
                    labelColor = "#FF8800";
                    capText += " - I2C1";
                } else {
                    labelColor = "#0099FF";
                    capText += " - I2C0";
                }
            }
        }

        const wrapper = document.createElement("div");
        wrapper.className = "pin-control";
        wrapper.dataset.order = order;

        let optionsHtml = `<option value="-1">(${i18n.t("Unmapped")})</option>`;
        for (const f of allFunctions) {
            const selected = (currentFunc === f) ? "selected" : "";
            optionsHtml += `<option value="${f}" ${selected}>${i18n.t(f)}</option>`;
        }

        wrapper.innerHTML = `
            <span class="pin-label" style="color: ${labelColor}" title="${capText}">«GPIO${gpio}»</span>
            <select class="pin-select" data-gpio="${gpio}" ${isCustomPins ? "" : "disabled"}>${optionsHtml}</select>
        `;

        if (group === posLeft) { wrapper.style.flexDirection = 'row-reverse'; elementsLeft.push(wrapper); }
        else if (group === posRight) {
            elementsRight.push(wrapper);
        }
        else if (group === posMiddle) {
            wrapper.style.flexDirection = "column";
            elementsMiddle.push(wrapper);
        }
    }

    elementsLeft.sort((a, b) => a.dataset.order - b.dataset.order).forEach(el => pinsLeft.appendChild(el));
    elementsRight.sort((a, b) => a.dataset.order - b.dataset.order).forEach(el => pinsRight.appendChild(el));
    elementsMiddle.sort((a, b) => a.dataset.order - b.dataset.order).forEach(el => pinsMiddle.appendChild(el));
}

function populateSettingsUI(gunConfig) {
    // Toggles (Checkbox)
    const togglesMap = {
        "Autofire": "tgl-Autofire",
        "Solenoid": "tgl-Solenoid",
        "Rumble": "tgl-Rumble",
        "RumbFFB": "tgl-RumbFFB",
        "LEDAnode": "tgl-LEDAnode",
        "LowButtons": "tgl-LowButtons",
        "SimplePause": "tgl-SimplePause",
        "HoldToPause": "tgl-HoldToPause",
        "i2cOLED": "tgl-i2cOLED",
        "i2cOLEDaltAddr": "tgl-i2cOLEDaltAddr"
    };

    for (const [key, id] of Object.entries(togglesMap)) {
        const el = document.getElementById(id);
        if (el && gunConfig.toggles.hasOwnProperty(key)) {
            el.checked = gunConfig.toggles[key];
        }
    }

    // Settings (Number inputs)
        const settingsMap = {
        "SolOn": "set-SolenoidTime",
        "SolOff": "set-SolenoidAutofireTime",
        "SolHold": "set-AutofireTriggerTime",
        "RumbPwr": "set-RumblePWM",
        "RumbTime": "set-RumbleTime",
        "CtmPixelsCount": "set-LEDPixelCount",
        "StaticPixels": "set-StaticPixelCount",
        "HoldToPauseLength": "set-HoldTime"
    };

    for (const [key, id] of Object.entries(settingsMap)) {
        const el = document.getElementById(id);
        if (el && gunConfig.settings.hasOwnProperty(key)) {
            el.value = gunConfig.settings[key];
        }
    }
}

function populateButtonsUI(gunConfig) {
    const tbody = document.getElementById("buttons-tbody");
    if (!tbody) return;
    tbody.innerHTML = "";

    const btnEntries = Object.entries(OpenFIREshared.boardInputs_Strings)
                             .filter(([name, val]) => val >= 0)
                             .sort((a, b) => a[1] - b[1]);

        const typeOptions = OpenFIREmaps.funcTypes.map((t, idx) => `<option value="${idx}">${i18n.t(t)}</option>`).join('');

    // Analog Stick Mode
    const analogSel = document.getElementById("sel-AnalogStickMode");
    if (analogSel) {
        analogSel.innerHTML = `
            <option value="0">${i18n.t("Gamepad Analog Stick (Left/Right)")}</option>
            <option value="1">${i18n.t("Gamepad D-Pad")}</option>
            <option value="2">${i18n.t("Keyboard Arrows")}</option>
        `;
        if (gunConfig.settings["Analog Mode"] !== undefined) {
            analogSel.value = gunConfig.settings["Analog Mode"];
        } else if (gunConfig.settings["Analog Stick Mode"] !== undefined) {
            analogSel.value = gunConfig.settings["Analog Stick Mode"];
        }
    }

    for (const [btnName, btnEnum] of btnEntries) {
        const btnData = gunConfig.buttons[btnName];
        if (!btnData || btnData.length < 6) continue;

        const onType = btnData[0];
        const onVal = btnData[1];
        const offType = btnData[2];
        const offVal = btnData[3];
        const gpType = btnData[4];
        const gpVal = btnData[5];

        const tr = document.createElement("tr");

        const tdName = document.createElement("td");
        tdName.innerText = i18n.t(btnName);
        tr.appendChild(tdName);

        const createCell = (type, val, isGamepadMode) => {
            const tdType = document.createElement("td");
            const tdVal = document.createElement("td");
            
            const selType = document.createElement("select");
            selType.innerHTML = typeOptions;
            selType.value = type;
            
            const selVal = document.createElement("select");
            const populateSelVal = (t) => {
                let mapArr = [];
                if (t == 0) mapArr = OpenFIREmaps.mouseMap;
                else if (t == 1) mapArr = OpenFIREmaps.keyboardMap;
                else if (t == 2) mapArr = OpenFIREmaps.gamepadMap;
                
                selVal.innerHTML = mapArr.map(item => `<option value="${item.val}">${i18n.t(item.name)}</option>`).join('');
            };
            
            selType.addEventListener('change', () => {
                populateSelVal(selType.value);
            });
            
            populateSelVal(type);
            selVal.value = val;

            if (isGamepadMode) {
                selType.disabled = true;
                selType.value = 2; // Forzato a Gamepad
                populateSelVal(2);
                selVal.value = gpVal;
            }

            tdType.appendChild(selType);
            tdVal.appendChild(selVal);
            tr.appendChild(tdType);
            tr.appendChild(tdVal);
        };

        createCell(onType, onVal, false);
        createCell(offType, offVal, false);
        createCell(gpType, gpVal, true);

        tbody.appendChild(tr);
    }
}

function populateProfilesUI(gunConfig) {
    const tbody = document.getElementById("profiles-tbody");
    if (!tbody) return;
    tbody.innerHTML = "";

    for (let i = 0; i < 4; i++) {
        const prof = gunConfig.profiles[i] || {};
        const isCurrent = (gunConfig.currentProfile === i);
        
        const tr = document.createElement("tr");

        // Nome
        const nameVal = prof.Name ? prof.Name.replace(/"/g, '&quot;') : `Profile ${i+1}`;
        
        // Offset
        const formatVal = (v) => v === undefined ? 0 : (Number.isInteger(v) ? v : Number(v).toFixed(2));
        const top = formatVal(prof.TopOffset);
        const btm = formatVal(prof.BtmOffset);
        const lft = formatVal(prof.LftOffset);
        const rht = formatVal(prof.RhtOffset);
        const tlled = formatVal(prof.TLLed);
        const trled = formatVal(prof.TRLed);

        // Sensibilità
        const optsSens = ["Default", "Higher", "Highest"].map((s, idx) => `<option value="${idx}" ${(prof.IrSens==idx)?'selected':''}>${i18n.t(s)}</option>`).join('');
        // Modalità
        const optsMode = ["Normal", "1-Frame Avg", "2-Frame Avg"].map((s, idx) => `<option value="${idx}" ${(prof.IrRunMode==idx)?'selected':''}>${i18n.t(s)}</option>`).join('');
        // Layout
        const optsLayout = ["Square", "Diamond"].map((s, idx) => `<option value="${idx}" ${(prof.IrLayout==idx)?'selected':''}>${i18n.t(s)}</option>`).join('');
        // Display
        const optsAr = ["16:9", "16:10", "3:2", "5:4", "4:3"].map((s, idx) => `<option value="${idx}" ${(prof.AspectRatio==idx)?'selected':''}>${s}</option>`).join(''); // AR labels don't need translation usually

        // Colore
        let colHex = "#000000";
        if (prof.Color !== undefined) {
            colHex = "#" + ("000000" + prof.Color.toString(16)).slice(-6);
        }

        tr.innerHTML = `
            <td><input type="radio" name="currentProf" value="${i}" ${isCurrent ? 'checked' : ''} style="transform: scale(1.5);"></td>
            <td><input type="text" value="${nameVal}" style="width: 100px; text-align: left;"></td>
            <td style="color:#aaa;">${top}</td>
            <td style="color:#aaa;">${btm}</td>
            <td style="color:#aaa;">${lft}</td>
            <td style="color:#aaa;">${rht}</td>
            <td style="color:#aaa;">${tlled}</td>
            <td style="color:#aaa;">${trled}</td>
            <td><select data-prof="${i}" data-field="IrSens">${optsSens}</select></td>
            <td><select data-prof="${i}" data-field="IrRunMode">${optsMode}</select></td>
            <td><select data-prof="${i}" data-field="IrLayout">${optsLayout}</select></td>
            <td><select data-prof="${i}" data-field="AspectRatio">${optsAr}</select></td>
            <td><input type="color" data-prof="${i}" data-field="Color" value="${colHex}" style="padding: 0; width: 25px; height: 25px; cursor: pointer; border: none; border-radius: 4px;"></td>
        `;

        tbody.appendChild(tr);
    }
}

function initGunTestsUI() {
    const grid = document.getElementById("tests-buttons-grid");
    if (!grid) return;
    grid.innerHTML = "";

    const btnEntries = Object.entries(OpenFIREshared.boardInputs_Strings)
                             .filter(([name, val]) => val >= 0)
                             .sort((a, b) => a[1] - b[1]);

    window.testBtnIndicators = [];
    
    for (const [btnName, btnIndex] of btnEntries) {
        const div = document.createElement("div");
        div.className = "test-btn-indicator";
        div.innerText = i18n.t(btnName);
        div.id = `test-btn-${btnIndex}`;
        grid.appendChild(div);
        window.testBtnIndicators[btnIndex] = div;
    }
}


function gatherGunConfigFromUI() {
    const gc = window.gunConfig;
    if (!gc) return;

    // Toggles
    const togglesMap = {
        "Autofire": "tgl-Autofire", "Solenoid": "tgl-Solenoid", "Rumble": "tgl-Rumble", "RumbFFB": "tgl-RumbFFB",
        "LEDAnode": "tgl-LEDAnode", "LowButtons": "tgl-LowButtons", "SimplePause": "tgl-SimplePause",
        "HoldToPause": "tgl-HoldToPause", "i2cOLED": "tgl-i2cOLED", "i2cOLEDaltAddr": "tgl-i2cOLEDaltAddr",
        "CustomPins": "tgl-CustomPins"
    };
    for (const [key, id] of Object.entries(togglesMap)) {
        const el = document.getElementById(id);
        if (el) gc.toggles[key] = el.checked;
    }

    // Pins
    for (const [key, val] of Object.entries(gc.pins)) {
        const sel = document.getElementById(`pin-${key}`);
        if (sel) gc.pins[key] = parseInt(sel.value);
    }

    // Settings
        const settingsMap = {
        "SolOn": "set-SolenoidTime",
        "SolOff": "set-SolenoidAutofireTime",
        "SolHold": "set-AutofireTriggerTime",
        "RumbPwr": "set-RumblePWM",
        "RumbTime": "set-RumbleTime",
        "CtmPixelsCount": "set-LEDPixelCount",
        "StaticPixels": "set-StaticPixelCount",
        "HoldToPauseLength": "set-HoldTime"
    };
    for (const [key, id] of Object.entries(settingsMap)) {
        const el = document.getElementById(id);
        if (el) gc.settings[key] = parseInt(el.value);
    }
    const selAnalog = document.getElementById("sel-AnalogStickMode");
    if (selAnalog) gc.settings["Analog Mode"] = parseInt(selAnalog.value);

    // Buttons
    const tbodyBtn = document.getElementById("buttons-tbody");
    if (tbodyBtn) {
        const rows = tbodyBtn.querySelectorAll("tr");
        rows.forEach(tr => {
            const btnNameStr = tr.cells[0].innerText;
            // Reverse translation is hard, let's just find the key that matches i18n.t(key)
            const btnKey = Object.keys(OpenFIREshared.boardInputs_Strings).find(k => i18n.t(k) === btnNameStr || k === btnNameStr);
            if (btnKey && gc.buttons[btnKey]) {
                gc.buttons[btnKey][0] = parseInt(tr.cells[1].querySelector("select").value);
                gc.buttons[btnKey][1] = parseInt(tr.cells[2].querySelector("select").value);
                gc.buttons[btnKey][2] = parseInt(tr.cells[3].querySelector("select").value);
                gc.buttons[btnKey][3] = parseInt(tr.cells[4].querySelector("select").value);
                gc.buttons[btnKey][4] = parseInt(tr.cells[5].querySelector("select").value);
                gc.buttons[btnKey][5] = parseInt(tr.cells[6].querySelector("select").value);
            }
        });
    }

    // Profiles
    const tbodyProf = document.getElementById("profiles-tbody");
    if (tbodyProf) {
        const rows = tbodyProf.querySelectorAll("tr");
        rows.forEach((tr, i) => {
            if (tr.querySelector("input[type=radio]").checked) gc.currentProfile = i;
            if (!gc.profiles[i]) gc.profiles[i] = {};
            
            gc.profiles[i].Name = tr.querySelector("input[type=text]").value;
            
            const selects = tr.querySelectorAll("select");
            selects.forEach(sel => {
                gc.profiles[i][sel.dataset.field] = parseInt(sel.value);
            });

            const colInp = tr.querySelector("input[type=color]");
            if (colInp) {
                const hex = colInp.value.replace("#", "");
                gc.profiles[i].Color = parseInt(hex, 16);
            }
        });
    }
}

    // --- Camera Tester Overlay ---
    let cameraTesterActive = false;
    const btnCam = document.getElementById("btn-test-camera");
    const overlayCam = document.getElementById("overlay-camera");
    const canvasCam = document.getElementById("canvas-camera");
    
    if (btnCam && overlayCam && canvasCam) {
        const ctx = canvasCam.getContext("2d");
        
        btnCam.addEventListener("click", async () => {
            cameraTesterActive = true;
            overlayCam.style.display = "flex";
            if(document.documentElement.requestFullscreen) {
                try {
                    await document.documentElement.requestFullscreen();
                    if (navigator.keyboard && navigator.keyboard.lock) await navigator.keyboard.lock(['Escape']);
                } catch (e) {}
            }
            
            canvasCam.width = window.innerWidth;
            canvasCam.height = window.innerHeight;
            
            window.ofProtocol.sendCommand(OpenFIREshared.serialCmdTypes_e.sIRTest, new Uint8Array([1]));
        });
        
        document.addEventListener("keydown", (e) => {
            if (e.key === "Escape" && cameraTesterActive) {
                e.preventDefault();
                cameraTesterActive = false;
                overlayCam.style.display = "none";
                if(document.fullscreenElement) document.exitFullscreen().catch(()=>{});
                if(navigator.keyboard && navigator.keyboard.unlock) navigator.keyboard.unlock();
                window.ofProtocol.sendCommand(OpenFIREshared.serialCmdTypes_e.sIRTest, new Uint8Array([0]));
            }
        });
        
        // Expose a draw function
        window.drawCameraTest = (coordsList) => {
            if (!cameraTesterActive) return;
            
            const w = canvasCam.width;
            const h = canvasCam.height;
            ctx.clearRect(0, 0, w, h);
            
            // The coordinates are presumably scaled to 1920x1080
            const scaleX = w / 1920.0;
            const scaleY = h / 1080.0;
            const scale = Math.min(scaleX, scaleY);
            const offsetX = (w - (1920.0 * scale)) / 2.0;
            const offsetY = (h - (1080.0 * scale)) / 2.0;
            
            ctx.save();
            ctx.translate(offsetX, offsetY);
            ctx.scale(scale, scale);
            
            const pointX = [0,0,0,0];
            const pointY = [0,0,0,0];
            const outsideFov = [false, false, false, false];
            
            for(let i=0; i<4; i++) {
                const encodedX = coordsList[i * 2];
                outsideFov[i] = (encodedX % 2) !== 0;
                pointX[i] = (encodedX - (outsideFov[i] ? 1 : 0)) / 2;
                pointY[i] = coordsList[(i * 2) + 1];
            }
            
            // Draw Box
            ctx.beginPath();
            ctx.moveTo(pointX[0], pointY[0]);
            ctx.lineTo(pointX[1], pointY[1]);
            ctx.lineTo(pointX[3], pointY[3]);
            ctx.lineTo(pointX[2], pointY[2]);
            ctx.closePath();
            ctx.strokeStyle = "gray";
            ctx.lineWidth = 2;
            ctx.stroke();
            
            // Draw Points
            const colors = ["lime", "lime", "cyan", "cyan"]; // TL, TR, BL, BR
            for(let i=0; i<4; i++) {
                ctx.beginPath();
                ctx.arc(pointX[i], pointY[i], 25, 0, 2*Math.PI);
                ctx.strokeStyle = colors[i];
                ctx.lineWidth = 3;
                if (outsideFov[i]) {
                    ctx.fillStyle = colors[i];
                    ctx.fill();
                } else {
                    ctx.stroke();
                }
            }
            
            // Point D (Red)
            ctx.beginPath();
            ctx.arc(coordsList[10], coordsList[11], 25, 0, 2*Math.PI);
            ctx.strokeStyle = "red";
            ctx.lineWidth = 3;
            ctx.stroke();
            
            ctx.restore();
            
            // Point Med (Gray) - scaled to full window
            ctx.save();
            ctx.scale(scaleX, scaleY);
            ctx.beginPath();
            ctx.arc(coordsList[8], coordsList[9], 25, 0, 2*Math.PI);
            ctx.strokeStyle = "gray";
            ctx.lineWidth = 3;
            ctx.stroke();
            ctx.restore();
        };
    }

function initBoardPreviewUI() {
    const sel = document.getElementById("preview-board-select");
    const container = document.getElementById("preview-board-container");
    if (!sel || !container) return;
    
    sel.innerHTML = "";
    Object.keys(OpenFIREshared.boardsBoxPositions).forEach(key => {
        if (key.includes("generic")) return; // Match Qt App: skip "generic" boards
        const name = (OpenFIREshared.boardNames && OpenFIREshared.boardNames[key]) ? OpenFIREshared.boardNames[key] : key;
        sel.innerHTML += `<option value="${key}">${name}</option>`;
    });

    const drawPreview = (boardName) => {
        const boxPositions = OpenFIREshared.boardsBoxPositions[boardName];
        if (!boxPositions) return;
        const presets = OpenFIREshared.boardsPresetsMap[boardName] || [];
        let capab = OpenFIREshared.mcuCapableMaps[boardName];
        if (!capab) {
            if (boardName.includes('esp32-s3')) capab = OpenFIREshared.mcuCapableMaps['esp32-s3'];
            else capab = OpenFIREshared.mcuCapableMaps['rp2040_235X'];
        }
        if (!capab) capab = [];

        // Reverse map from boardInputs_e value to name string
        const inputMapReverse = {};
        Object.entries(OpenFIREshared.boardInputs_Strings).forEach(([k, v]) => inputMapReverse[v] = k);

        let htmlLeft = "";
        let htmlRight = "";
        let htmlMiddle = "";

        const posLeft = OpenFIREshared.boardBoxPositions_e.posLeft;
        const posRight = OpenFIREshared.boardBoxPositions_e.posRight;
        const posMiddle = OpenFIREshared.boardBoxPositions_e.posMiddle;
        const posCheck = OpenFIREshared.boardBoxPositions_e.posCheck;

        const middleItems = [];
        const leftMap = {};
        const rightMap = {};
        let maxLeft = 0;
        let maxRight = 0;

        boxPositions.forEach((val, gpioPin) => {
            if (val === 0) return;
            const group = val & posCheck;
            const order = val ^ group;

            const funcVal = presets[gpioPin];
            let funcName = "Unmapped";
            if (funcVal !== undefined && inputMapReverse[funcVal]) funcName = inputMapReverse[funcVal];
            if (funcVal === OpenFIREshared.boardInputs_e.unavailable) funcName = "Unavailable";

            const cap = capab[gpioPin] || 0;
            const OF_Const = OpenFIREshared.pinCapabilities_e;
            const cStr = [];
            if (cap & OF_Const.pinHasADC) cStr.push('<span style="color:#FF0099; font-family:monospace; font-weight:bold;">ADC</span>');
            else cStr.push('<span style="color:#555555; font-family:monospace;">ADC</span>');
            
            let gpioColor = '#aaaaaa';
            if (cap & OF_Const.pinAnyI2C) {
                gpioColor = '#BE00B0';
                cStr.push('<span style="color:#BE00B0; font-family:monospace; font-weight:bold;">I2C(*)</span>');
            } else if (cap & OF_Const.pinCanI2C) {
                const isI2C1 = (cap & OF_Const.pinIsI2C1);
                const isSCL = (cap & OF_Const.pinIsI2CSCL);
                gpioColor = isI2C1 ? '#FF8800' : '#0099FF';
                cStr.push(`<span style="color:${gpioColor}; font-family:monospace; font-weight:bold;">I2C${isI2C1 ? '1' : '0'}${isSCL ? 'SCL' : 'SDA'}</span>`);
            } else {
                cStr.push('<span style="color:#555555; font-family:monospace;">I2C</span>');
            }
            
            if (cap & OF_Const.pinAnySPI) {
                cStr.push('<span style="color:#D1003D; font-family:monospace; font-weight:bold;">SPI(*)</span>');
            } else if (cap & OF_Const.pinCanSPI) {
                const isSPI1 = (cap & OF_Const.pinIsSPI1);
                let spiFunc = '';
                const st = (cap & OF_Const.pinCanSPI) >> 5;
                if (st === 1) spiFunc = 'RX';
                else if (st === 2) spiFunc = 'TX';
                else if (st === 3) spiFunc = 'SCK';
                else if (st === 4) spiFunc = 'CSn';
                cStr.push(`<span style="color:#009C3A; font-family:monospace; font-weight:bold;">SPI${isSPI1 ? '1' : '0'}${spiFunc}</span>`);
            } else {
                cStr.push('<span style="color:#555555; font-family:monospace;">SPI</span>');
            }
            
            const capHtml = `<span style="font-size:10px; margin:0 5px; white-space:nowrap;">${cStr.join(' ')}</span>`;
            const gpioHtml = `<span style="color:${gpioColor}; font-size:12px; white-space:nowrap;">«GPIO${gpioPin}»</span>`;
            const funcHtml = `<span id="func-gpio-${gpioPin}" style="color:#ddd; font-size:14px; transition: text-shadow 0.1s;">${i18n.t(funcName)}</span>`;
            
            // Hover logic added via inline events
            const hoverEvents = `onmouseenter="highlightBoardPin(${gpioPin}, true)" onmouseleave="highlightBoardPin(${gpioPin}, false)"`;
            
                        const tooltipText = `Pin GPIO N. ${gpioPin}.\n\nI pin con numero blu appartengono a I2C0.\nQuelli con numero arancione appartengono a I2C1.\nQuelli in viola possono selezionare automaticamente qualsiasi canale I2C.\nQuelli grigi non supportano I2C.\n\nADC indica la capacità di leggere input analogici.\nI2C e SPI indicano se il pin supporta tali dispositivi e su quale canale.\n(*) significa che il pin può usare la funzione tramite canali selezionabili dal software.`;
            const labelStr = (group === posLeft) ? 
                `<tr ${hoverEvents} style="background:transparent; cursor:default;" title="${tooltipText}"><td style="border:none; text-align:right; padding:2px 5px;">${funcHtml}</td><td style="border:none; text-align:center; padding:2px 5px;">${gpioHtml}</td><td style="border:none; text-align:left; padding:2px 5px;">${capHtml}</td></tr>` : 
                `<tr ${hoverEvents} style="background:transparent; cursor:default;" title="${tooltipText}"><td style="border:none; text-align:right; padding:2px 5px;">${capHtml}</td><td style="border:none; text-align:center; padding:2px 5px;">${gpioHtml}</td><td style="border:none; text-align:left; padding:2px 5px;">${funcHtml}</td></tr>`;
            
            if (group === posLeft) { leftMap[order] = labelStr; maxLeft = Math.max(maxLeft, order); }
            else if (group === posRight) { rightMap[order] = labelStr; maxRight = Math.max(maxRight, order); }
            else if (group === posMiddle) { middleItems.push({ order, html: `<div ${hoverEvents} style="display:flex; flex-direction:column; align-items:center; margin: 0 10px; cursor:default;" title="${tooltipText}"><div style="margin-bottom:2px;">${capHtml}</div><div style="margin-bottom:2px;">${gpioHtml}</div><div>${funcHtml}</div></div>` }); }
        });

        for(let i = 1; i <= maxLeft; i++) {
            if (leftMap[i]) htmlLeft += leftMap[i];
            else htmlLeft += `<tr style="background:transparent;"><td colspan="3" style="border:none; height:22px;"></td></tr>`;
        }
        for(let i = 1; i <= maxRight; i++) {
            if (rightMap[i]) htmlRight += rightMap[i];
            else htmlRight += `<tr style="background:transparent;"><td colspan="3" style="border:none; height:22px;"></td></tr>`;
        }
        middleItems.sort((a,b)=>a.order-b.order).forEach(x => htmlMiddle += x.html);

        container.innerHTML = `
<div style="width:100%; display:flex; flex-direction:column; align-items:center; justify-content:center;">
<div style="width:100%; display:flex; align-items:center; justify-content:center;">
                <div style="flex: 1 1 200px; min-width:150px; max-width:350px; text-align:right; padding-right:10px;"><table style="width:100%; border-collapse:collapse; background:transparent;">${htmlLeft}</table></div>
                <div style="flex: 1 1 150px; min-width:100px; max-width:350px; display:flex; flex-direction:column; align-items:center;">
                    <div id="board-svg-container" style="display:flex; align-items:center; justify-content:center; padding: 0 15px;">
                          ${OpenFIREshared.boardSVGsMap && OpenFIREshared.boardSVGsMap[boardName] ? OpenFIREshared.boardSVGsMap[boardName] : '<img src="boardPics/' + (OpenFIREshared.boardImagesMap[boardName] || 'generic.svg') + '" style="height: 100%; width: auto;">'}
                      </div>
                    <div style="width:100%; display:flex; justify-content:center; gap:20px; margin-top:10px;">${htmlMiddle}</div>
                </div>
                <div style="flex: 1 1 200px; min-width:150px; max-width:350px; text-align:left; padding-left:10px;"><table style="width:100%; border-collapse:collapse; background:transparent;">${htmlRight}</table></div>
            </div>
                          <div style="width: 100%; text-align:center; margin-top:20px; padding-top:10px; border-top:1px solid #444; font-size:14px; color:#ddd; line-height:1.4;">
                  ${boardName.includes('esp32-s3') ? 
                      `${i18n.t("Compatible with the <a href='https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32' target='_blank'><span style='text-decoration: underline; color:#8ab4f8;'>ESP-IDF fork of the OpenFIRE Firmware</span></a> by <i>Alessandro Satanassi.</i>")}<br>${i18n.t("Any issues should be reported <b><a href='https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/issues' target='_blank'><span style='text-decoration: underline; color:#8ab4f8;'>here!</span></a></b>")}` : 
                      i18n.t("Compatible with <a href='https://github.com/TeamOpenFIRE/OpenFIRE-Firmware' target='_blank'><span style='text-decoration: underline; color:#8ab4f8;'>upstream OpenFIRE Firmware</span></a> by <i>Team OpenFIRE.</i>")}
              </div>
</div>
`;
        
        // Fix SVG styling to match container with proper min/max bounds
        const svgEl = container.querySelector('#board-svg-container svg');
        if (svgEl) {
            if (svgEl.hasAttribute('width')) svgEl.removeAttribute('width');
            if (svgEl.hasAttribute('height')) svgEl.removeAttribute('height');
            
            svgEl.style.height = '55vh';
            svgEl.style.minHeight = '300px';
            svgEl.style.maxHeight = '600px';
            svgEl.style.width = '100%';
            svgEl.style.maxWidth = '320px';
            svgEl.style.display = 'block';
        }
    };

    sel.onchange = () => drawPreview(sel.value);
    
    // Select default or current
    if (window.gunConfig && window.gunConfig.boardName && OpenFIREshared.boardsBoxPositions[window.gunConfig.boardName]) {
        sel.value = window.gunConfig.boardName;
    }
        if (document.getElementById("menu-btn-preview").style.display !== "none") {
        drawPreview(sel.value);
    }
}

    // --- Calibration Overlay Logic ---
    let caliActive = false;
    const btnCaliStart = document.getElementById("btn-cali-start");
    const overlayCali = document.getElementById("overlay-cali");

    if (btnCaliStart && overlayCali) {
        // Create a crosshair element dynamically
        const crosshair = document.createElement("div");
        crosshair.innerHTML = `<svg width="50" height="50" viewBox="0 0 50 50">
            <circle cx="25" cy="25" r="20" stroke="red" stroke-width="4" fill="none"/>
            <line x1="25" y1="0" x2="25" y2="50" stroke="red" stroke-width="4"/>
            <line x1="0" y1="25" x2="50" y2="25" stroke="red" stroke-width="4"/>
        </svg>`;
        crosshair.style.position = "absolute";
        crosshair.style.transform = "translate(-50%, -50%)"; // Center the point exactly
        crosshair.style.transition = "top 0.3s, left 0.3s";
        crosshair.style.display = "none";
        overlayCali.appendChild(crosshair);

        btnCaliStart.addEventListener("click", async () => {
            caliActive = true;
            overlayCali.style.display = "flex";
            if(document.documentElement.requestFullscreen) {
                try {
                    await document.documentElement.requestFullscreen();
                    if (navigator.keyboard && navigator.keyboard.lock) await navigator.keyboard.lock(['Escape']);
                } catch (e) {}
            }
            
            // Start command: sCaliProfile (6) -> [sCaliStart (7), profileNum, (irSens) + (layoutType << 4)]
            const prof = window.gunConfig.currentProfile;
            const pObj = window.gunConfig.profiles[prof];
            const irSens = pObj ? pObj.IrSens : 0;
            const layout = pObj ? pObj.IrLayout : 0;

            const payload = new Uint8Array([7, prof, irSens + (layout << 4)]);
            window.ofProtocol.sendCommand(OpenFIREshared.serialCmdTypes_e.sCaliProfile, payload);
            
            // Set initial state
            crosshair.style.left = "50%";
            crosshair.style.top = "50%";
            crosshair.style.display = "block";
            
            const topText = overlayCali.querySelector('.top-text');
            if (topText) topText.innerHTML = i18n.t("Start calibration:<br>Shoot the target in the center to begin.");
        });

        document.addEventListener("keydown", (e) => {
            if (e.key === "Escape" && caliActive) {
                e.preventDefault();
                caliActive = false;
                overlayCali.style.display = "none";
                if(document.fullscreenElement) document.exitFullscreen().catch(()=>{});
                if(navigator.keyboard && navigator.keyboard.unlock) navigator.keyboard.unlock();
                // Abort with serialTerminator
                window.ofProtocol.sendCommand(OpenFIREshared.serialCmdTypes_e.serialTerminator);
            }
        });

        window.updateCaliStage = (stage) => {
            if (!caliActive) return;
            const topText = overlayCali.querySelector('.top-text');
            if (!topText) return;

            // Cali_Init = 0, Top = 1, Bottom = 2, Left = 3, Right = 4, Center = 5, Verify = 6, End = 7
            switch (stage) {
                case 0: // Init
                    crosshair.style.left = "50%"; crosshair.style.top = "50%";
                    topText.innerHTML = i18n.t("Start calibration:<br>Shoot the target in the center to begin.");
                    break;
                case 1: // Top
                    crosshair.style.left = "50%"; crosshair.style.top = "0%";
                    topText.innerHTML = i18n.t("Step 1:<br>Shoot the target at the TOP edge of the screen.");
                    break;
                case 2: // Bottom
                    crosshair.style.left = "50%"; crosshair.style.top = "100%";
                    topText.innerHTML = i18n.t("Step 2:<br>Shoot the target at the BOTTOM edge of the screen.");
                    break;
                case 3: // Left
                    crosshair.style.left = "0%"; crosshair.style.top = "50%";
                    topText.innerHTML = i18n.t("Step 3:<br>Shoot the target at the LEFT edge of the screen.");
                    break;
                case 4: // Right
                    crosshair.style.left = "100%"; crosshair.style.top = "50%";
                    topText.innerHTML = i18n.t("Step 4:<br>Shoot the target at the RIGHT edge of the screen.");
                    break;
                case 5: // Center
                    crosshair.style.left = "50%"; crosshair.style.top = "50%";
                    topText.innerHTML = i18n.t("Step 5:<br>Shoot the target at the CENTER of the screen.");
                    break;
                case 6: // Verify
                    crosshair.style.left = "50%"; crosshair.style.top = "50%";
                    topText.innerHTML = i18n.t("Verify Calibration:<br>Shoot off-screen to save.");
                    break;
                case 7: // End
                    caliActive = false;
                    overlayCali.style.display = "none";
                if(document.fullscreenElement) document.exitFullscreen();
                    break;
            }
        };
    }

// ============================================================================
document.addEventListener("DOMContentLoaded", () => {
    // Inizializza subito la UI delle preview in modo che sia disponibile prima della connessione
    initBoardPreviewUI();
    if (window.location.protocol === "http:" && window.location.hostname !== "localhost" && window.location.hostname !== "127.0.0.1") {
        document.getElementById('menu-btn-preview').style.display = 'none';
    }


    // Emitter Alignment Overlay
    const btnAlign = document.getElementById('menu-btn-align');
    const overlayAlign = document.getElementById('overlay-align');
    const canvasAlign = document.getElementById('canvas-align');
    if (btnAlign && overlayAlign && canvasAlign) {
        const drawAlign = () => {
            if (overlayAlign.style.display !== 'block') return;
            const ctx = canvasAlign.getContext('2d');
            const w = window.innerWidth;
            const h = window.innerHeight;
            canvasAlign.width = w;
            canvasAlign.height = h;
            
            ctx.clearRect(0, 0, w, h);
            
            const centerX = w / 2;
            const centerY = h / 2;
            const offset = (h * 0.711) / 2;
            const leftX = centerX - offset;
            const rightX = centerX + offset;
            
            const boxW = Math.max(30, w * 0.02);
            const boxH = Math.max(20, h * 0.02);
            
            ctx.strokeStyle = "firebrick";
            ctx.lineWidth = 2;
            ctx.strokeRect(leftX, 0, rightX - leftX, h);
            
            ctx.strokeStyle = "olivedrab";
            ctx.beginPath();
            ctx.moveTo(centerX, 0);
            ctx.lineTo(w, centerY);
            ctx.lineTo(centerX, h);
            ctx.lineTo(0, centerY);
            ctx.closePath();
            ctx.stroke();
            
            ctx.fillStyle = "firebrick";
            ctx.fillRect(leftX - boxW/2, 0, boxW, boxH);
            ctx.fillRect(rightX - boxW/2, 0, boxW, boxH);
            ctx.fillRect(leftX - boxW/2, h - boxH, boxW, boxH);
            ctx.fillRect(rightX - boxW/2, h - boxH, boxW, boxH);
            
            ctx.fillStyle = "olivedrab";
            ctx.fillRect(centerX - boxW/2, 0, boxW, boxH);
            ctx.fillRect(centerX - boxW/2, h - boxH, boxW, boxH);
            ctx.fillRect(0, centerY - boxW/2, boxH, boxW);
            ctx.fillRect(w - boxH, centerY - boxW/2, boxH, boxW);
        };

        btnAlign.addEventListener('click', async () => {
            overlayAlign.style.display = 'block';
            if(document.documentElement.requestFullscreen) {
                try {
                    await document.documentElement.requestFullscreen();
                    if (navigator.keyboard && navigator.keyboard.lock) await navigator.keyboard.lock(['Escape']);
                } catch (e) {}
            }
            drawAlign();
            setTimeout(drawAlign, 100);
            setTimeout(drawAlign, 500);
        });
        
        window.addEventListener('resize', drawAlign);
        
        document.addEventListener("keydown", (e) => {
            if (e.key === "Escape" && overlayAlign.style.display === "block") {
                e.preventDefault();
                overlayAlign.style.display = "none";
                if(document.fullscreenElement) document.exitFullscreen().catch(()=>{});
                if(navigator.keyboard && navigator.keyboard.unlock) navigator.keyboard.unlock();
            }
        });
    }

    // Modal Events


    document.getElementById('menu-btn-about').addEventListener('click', () => document.getElementById('modal-about').style.display = 'flex');
    document.getElementById('btn-close-about').addEventListener('click', () => document.getElementById('modal-about').style.display = 'none');
    
        // Gun Test Action Buttons
    const bindTest = (id, cmd) => {
        const btn = document.getElementById(id);
        if (btn) btn.addEventListener('mousedown', () => window.ofProtocol && window.ofProtocol.sendCommand(OpenFIREshared.serialCmdTypes_e[cmd]));
    };
    bindTest('btn-test-rumble', 'sTestRumble');
    bindTest('btn-test-solenoid', 'sTestSolenoid');
    bindTest('btn-test-led-r', 'sTestLEDR');
    bindTest('btn-test-led-g', 'sTestLEDG');
        bindTest('btn-test-led-b', 'sTestLEDB');
    bindTest('btn-restart-dfu', 'sRebootToBootloader');
    
    const btnFormat = document.getElementById('btn-format-mem');
    if (btnFormat) btnFormat.addEventListener('click', () => {
        if (confirm("Sei sicuro di voler cancellare tutta la memoria e riavviare la scheda?")) {
            window.ofProtocol && window.ofProtocol.sendCommand(OpenFIREshared.serialCmdTypes_e.sClearFlash);
        }
    });

    document.getElementById('menu-btn-preview').addEventListener('click', () => document.getElementById('modal-preview').style.display = 'flex');
    document.getElementById('btn-close-preview').addEventListener('click', () => document.getElementById('modal-preview').style.display = 'none');
    // Setup Tabs Navigation
    document.querySelectorAll('.tab-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(c => c.style.display = 'none');
            btn.classList.add('active');
            document.getElementById(btn.dataset.tab).style.display = 'block';
        });
    });

    // ============================================================================
    // Connection Logic
    // ============================================================================
    if (typeof OpenFIREConnection !== 'undefined') {
        window.ofProtocol = new OpenFIREConnection();
        window.ofProtocol.onEventReceived = (evt) => { console.log(evt); };

        async function doConnect() {
            const statusText = document.getElementById("status");
            statusText.innerText = i18n.t("Connecting...");
            const success = await window.ofProtocol.connect();
            
            if (success) {
                statusText.innerText = i18n.t("Connected! Starting Handshake...");
                try {
                    const boardInfo = await window.ofProtocol.beginDock();
                    
                    const imgBoard = document.getElementById("board-image");
                    if (imgBoard) {
                        if (window.ofProtocol.isWebSerial) {
                            imgBoard.src = `boardPics/${OpenFIREshared.boardImagesMap[boardInfo.boardName] || 'generic.svg'}`;
                        } else {
                            imgBoard.src = "board.svg";
                        }
                    }

                    statusText.innerText = i18n.t("Docked! Syncing Settings...");
                    
                    const gunConfig = await window.ofProtocol.syncSettings();
                    statusText.innerText = i18n.t("Sync complete! Active profile: ") + gunConfig.currentProfile;
                    
                    // Re-render things that depend on gunConfig
                    if (window.boardName) {
                        drawPreview(window.boardName, true);
                    }
                } catch (e) {
                    console.error(e);
                    statusText.innerText = i18n.t("Sync Error");
                }
            } else {
                statusText.innerText = i18n.t("Connection failed");
            }
        }

        if (window.location.protocol === 'file:') {
            const btn = document.createElement('button');
            btn.id = 'btn-connect-serial';
            btn.className = 'save-btn';
            btn.innerText = i18n.t('Connect (Web Serial)');
            btn.setAttribute('data-i18n', '');
            btn.style.marginLeft = '10px';
            btn.onclick = doConnect;
            document.querySelector('.status-bar').appendChild(btn);
            document.getElementById("status").innerText = i18n.t("Ready (File).");
        } else {
            // Auto-connect for WebSocket
            doConnect();
        }
    } else {
        document.getElementById("status").innerText = i18n.t("Ready (File).");
    }

});


// Global hover function for board preview
window.highlightBoardPin = function(gpioPin, isHover) {
    const pinObj = document.getElementById(`OF_pin${gpioPin}`);
    if (pinObj) {
        if (isHover) {
            // Check if we already have the original opacity stored
            if (!pinObj.dataset.origOpacity) {
                pinObj.dataset.origOpacity = pinObj.style.opacity || "0";
            }
            pinObj.style.opacity = "1";
        } else {
            pinObj.style.opacity = pinObj.dataset.origOpacity || "0";
        }
    }
    const funcSpan = document.getElementById(`func-gpio-${gpioPin}`);
    if (funcSpan) {
        funcSpan.style.webkitTextStroke = isHover ? "0.6px currentColor" : "0px";
        funcSpan.style.textShadow = isHover ? "0 0 1px rgba(255,255,255,0.3)" : "none";
    }
};
