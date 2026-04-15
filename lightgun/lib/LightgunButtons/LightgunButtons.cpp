/*!
 * @file LightgunButtons.cpp
 * @brief HID buttons originally intended for use with a light gun.
 *
 * @copyright Mike Lynch, 2021
 *
 *  LightgunButtons is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <Arduino.h>
#include "LightgunButtons.h"
#include <TinyUSB_Devices.h>
#include "../../src/boards/OpenFIREshared.h"  //696969 serve per gli indici di Pedal e Pedal2

// valutare 696969 tolti poichè già definiti in Arduino TinyUSB in hid.h // /// Standard Gamepad HAT/DPAD Buttons (from Linux input event codes)

#define GAMEPAD_HAT_CENTERED   hid_gamepad_hat_t::GAMEPAD_HAT_CENTERED   // 0
#define GAMEPAD_HAT_UP         hid_gamepad_hat_t::GAMEPAD_HAT_UP         // 1
#define GAMEPAD_HAT_UP_RIGHT   hid_gamepad_hat_t::GAMEPAD_HAT_UP_RIGHT   // 2
#define GAMEPAD_HAT_RIGHT      hid_gamepad_hat_t::GAMEPAD_HAT_RIGHT      // 3
#define GAMEPAD_HAT_DOWN_RIGHT hid_gamepad_hat_t::GAMEPAD_HAT_DOWN_RIGHT // 4
#define GAMEPAD_HAT_DOWN       hid_gamepad_hat_t::GAMEPAD_HAT_DOWN       // 5
#define GAMEPAD_HAT_DOWN_LEFT  hid_gamepad_hat_t::GAMEPAD_HAT_DOWN_LEFT  // 6
#define GAMEPAD_HAT_LEFT       hid_gamepad_hat_t::GAMEPAD_HAT_LEFT       // 7
#define GAMEPAD_HAT_UP_LEFT    hid_gamepad_hat_t::GAMEPAD_HAT_UP_LEFT    // 8


LightgunButtons::LightgunButtons(Data_t _data, unsigned int _count) :
    pressed(0),
    released(0),
    repeat(0),
    debounced(0),
    debouncing(0),
    pressedReleased(0),
    padMask(0),
    interval(33),
    report(0),
    lastMillis(0),
    lastRepeatMillis(0),
    pinState(0xFFFFFFFF),
    internalPressedReleased(0),
    reportedPressed(0),
    count(_count),
    stateFifo(_data.pArrFifo),
    debounceCount(_data.pArrDebounceCount)
{
}

void LightgunButtons::Begin()
{
    // set button pins to input with pullup
    for(unsigned int i = 0; i < count; ++i) {
        bool is_current_wireless_pedal = (TinyUSBDevices.is_pedal_wireless && (i == OF_Const::btnPedal || i == OF_Const::btnPedal2));  // 696969
        // do no setup if the pin is uninitialized
        if((ButtonDesc[i].pin >= 0) || is_current_wireless_pedal) {  //696969
            if (!is_current_wireless_pedal) pinMode(ButtonDesc[i].pin, INPUT_PULLUP);  //696969
            stateFifo[i] = 0xFFFFFFFF;
            debounceCount[i] = 0;
        }
    }
}

void LightgunButtons::Unset()
{
    // set button pins to normal input
    for(unsigned int i = 0; i < count; ++i) {
        bool is_current_wireless_pedal = (TinyUSBDevices.is_pedal_wireless && (i == OF_Const::btnPedal || i == OF_Const::btnPedal2));  // 696969
        // do no setup if the pin is uninitialized
        if((ButtonDesc[i].pin >= 0) || is_current_wireless_pedal) {  //696969
            if (!is_current_wireless_pedal) pinMode(ButtonDesc[i].pin, INPUT); // 696969
            debounceCount[i] = 0;
        }
    }
    pressed = 0;
    released = 0;
    debounced = 0;
    debouncing = 0;
    pressedReleased = 0;
    padMask = 0;
    lastMillis = 0;
    lastRepeatMillis = 0;
    internalPressedReleased = 0;
    reportedPressed = 0;
}

uint32_t LightgunButtons::Poll(unsigned long minTicks)
{
    unsigned long m = millis();
    unsigned long ticks = m - lastMillis;
    uint32_t bitMask;
    
    // reset pressed and released from last poll
    pressed = 0;
    released = 0;
    pressedReleased = 0;
    
    if(ticks < minTicks) {
        return 0;
    }
    lastMillis = m;    
    if(debouncing && ticks) {
        bitMask = 1;
        for(unsigned int i = 0; i < count; ++i, bitMask <<= 1) {
            const Desc_t& btn = ButtonDesc[i];
            // Condizione ottimizzata
            bool is_current_wireless_pedal = (TinyUSBDevices.is_pedal_wireless && (i == OF_Const::btnPedal || i == OF_Const::btnPedal2)); // 696969
            if((ButtonDesc[i].pin >= 0) || is_current_wireless_pedal) { // 696969
                if(debounceCount[i]) {
                    if(ticks < debounceCount[i]) {
                        debounceCount[i] -= ticks;
                    } else {
                        debounceCount[i] = 0;
                        debouncing &= ~bitMask;
                    }
                }
            }
        }
    }

    bitMask = 1;
    for(unsigned int i = 0; i < count; ++i, bitMask <<= 1) {
        const Desc_t& btn = ButtonDesc[i];
        
        bool is_current_wireless_pedal = (TinyUSBDevices.is_pedal_wireless && (i == OF_Const::btnPedal || i == OF_Const::btnPedal2));  // 696969
        // do no processing if the pin is uninitialized
        if((ButtonDesc[i].pin >= 0) || is_current_wireless_pedal) { // 696969
            // if not debouncing
            if(!debounceCount[i]) {
                uint32_t state; //696969
                
                // 696969
                if (is_current_wireless_pedal) {
                    // Mappa true (premuto) in 0, e false (rilasciato) in 1
                    //state = (i == OF_Const::btnPedal) ? (TinyUSBDevices.pedal_wireless ? 0 : 1) : (TinyUSBDevices.pedal2_wireless ? 0 : 1);
                    // Seleziona la maschera giusta in base al pedale che stiamo scansionando
                    uint8_t mask = (i == OF_Const::btnPedal) ? 0b00000001 : 0b00000010;
                    // Se il bit corrispondente è a 1 (premuto), state = 0, altrimenti state = 1
                    state = (TinyUSBDevices.pedals_wireless_state & mask) ? 0 : 1;
                } else {
                    // read the pin, expected to return 0 or 1
                    state = digitalRead(btn.pin);
                }
                // 696969

                // if a state fifo mask is defined
                if(btn.debounceFifoMask) {
                    // add the state to the fifo
                    stateFifo[i] <<= 1;
                    stateFifo[i] |= state;

                    // apply the mask and check the value
                    uint32_t m = stateFifo[i] & btn.debounceFifoMask;
                    if(!m) {
                        state = 0;
                    } else if(m == btn.debounceFifoMask) {
                        // use the bit mask for this button
                        state = bitMask;
                    } else {
                        // button is bouncing, continue to next button
                        continue;
                    }
                } else {
                    // no fifo, so change high state into bit mask
                    if(state) {
                        state = bitMask;
                    }
                }

                // if existing pin state does not match new state
                if((pinState & bitMask) != state) {
                    // update the pin state
                    pinState = (pinState & ~bitMask) | state;

                    // set the debounce counter and set the flag
                    debounceCount[i] = btn.debounceTicks;
                    debouncing |= bitMask;

                    if(!state) {
                        // state is low, button is pressed

                        // if reporting is enabled for the button (no trigger reporting)
                        if(i && report & bitMask) {
                            reportedPressed |= bitMask;
                            if(analogOutput) {
                                if(btn.reportType3 == ReportType_Mouse) {
                                    AbsMouse5.press(btn.reportCode3);
                                } else if(btn.reportType3 == ReportType_Keyboard) {
                                    Keyboard.press(btn.reportCode3);
                                } else if(btn.reportType3 == ReportType_Gamepad) {
                                    if(btn.reportCode3 < START_DPAD_KEY) {
                                        Gamepad16.press(btn.reportCode3);
                                    } else {
                                        bitSet(padMask, btn.reportCode3-START_DPAD_KEY);
                                        Gamepad16.padUpdate(PadMaskConvert(padMask));
                                    }
                                }
                            } else if(offScreen) {
                                bitSet(internalOffscreenMask, i);
                                if(btn.reportType2 == ReportType_Mouse) {
                                    AbsMouse5.press(btn.reportCode2);
                                } else if(btn.reportType2 == ReportType_Keyboard) {
                                    Keyboard.press(btn.reportCode2);
                                } else if(btn.reportType2 == ReportType_Gamepad) {
                                    if(btn.reportCode2 < START_DPAD_KEY) {
                                        Gamepad16.press(btn.reportCode2);
                                    } else {
                                        bitSet(padMask, btn.reportCode2-START_DPAD_KEY);
                                        Gamepad16.padUpdate(PadMaskConvert(padMask));
                                    }
                                }
                            } else {
                                if(btn.reportType == ReportType_Mouse) {
                                    AbsMouse5.press(btn.reportCode);
                                } else if(btn.reportType == ReportType_Keyboard) {
                                    Keyboard.press(btn.reportCode);
                                } else if(btn.reportType == ReportType_Gamepad) {
                                    if(btn.reportCode < START_DPAD_KEY) {
                                        Gamepad16.press(btn.reportCode);
                                    } else {
                                        bitSet(padMask, btn.reportCode-START_DPAD_KEY);
                                        Gamepad16.padUpdate(PadMaskConvert(padMask));
                                    }
                                }
                            }
                        }

                        // button is debounced pressed and add it to the pressed/released combo
                        debounced |= bitMask;
                        pressed |= bitMask;
                        internalPressedReleased |= bitMask;

                    } else {
                        // state high, button is not pressed

                        // if the button press was reported then report the release
                        // note that the report flag is ignored here to avoid stuck buttons
                        // in case the reporting is disabled while button(s) are pressed
                        if(i && reportedPressed & bitMask) {
                            reportedPressed &= ~bitMask;
                            if(analogOutput) {
                                if(btn.reportType3 == ReportType_Mouse) {
                                    AbsMouse5.release(btn.reportCode3);
                                } else if(btn.reportType3 == ReportType_Keyboard) {
                                    Keyboard.release(btn.reportCode3);
                                } else if(btn.reportType3 == ReportType_Gamepad) {
                                    if(btn.reportCode3 < START_DPAD_KEY) {
                                        Gamepad16.release(btn.reportCode3);
                                    } else {
                                        bitClear(padMask, btn.reportCode3-START_DPAD_KEY);
                                        Gamepad16.padUpdate(PadMaskConvert(padMask));
                                    }
                                }
                            } else if(bitRead(internalOffscreenMask, i)) {
                                bitClear(internalOffscreenMask, i);
                                if(btn.reportType2 == ReportType_Mouse) {
                                    AbsMouse5.release(btn.reportCode2);
                                } else if(btn.reportType2 == ReportType_Keyboard) {
                                    Keyboard.release(btn.reportCode2);
                                } else if(btn.reportType2 == ReportType_Gamepad) {
                                    if(btn.reportCode2 < START_DPAD_KEY) {
                                        Gamepad16.release(btn.reportCode2);
                                    } else {
                                        bitClear(padMask, btn.reportCode2-START_DPAD_KEY);
                                        Gamepad16.padUpdate(PadMaskConvert(padMask));
                                    }
                                }
                            } else {
                                if(btn.reportType == ReportType_Mouse) {
                                    AbsMouse5.release(btn.reportCode);
                                } else if(btn.reportType == ReportType_Keyboard) {
                                    Keyboard.release(btn.reportCode);
                                } else if(btn.reportType == ReportType_Gamepad) {
                                    if(btn.reportCode < START_DPAD_KEY) {
                                        Gamepad16.release(btn.reportCode);
                                    } else {
                                        bitClear(padMask, btn.reportCode-START_DPAD_KEY);
                                        Gamepad16.padUpdate(PadMaskConvert(padMask));
                                    }
                                }
                            }
                        }

                        // clear the debounced state and button is released
                        debounced &= ~bitMask;
                        released |= bitMask;

                        // if all buttons released
                        if(!debounced) {
                            // report the combination pressed/released state
                            pressedReleased = internalPressedReleased;
                            internalPressedReleased = 0;
                        }
                    }
                }
            }
        }
    }

    return pressed;
}

void LightgunButtons::SendReports(const bool &forceReportAll)
{
    
    
    if(TinyUSBDevices.newReport[ReportType_Mouse]) {
        AbsMouse5.report();
        if(!forceReportAll) return;
    }
    if(TinyUSBDevices.newReport[ReportType_Keyboard]) {
        Keyboard.report();
        if(!forceReportAll) return;
    }
    if(TinyUSBDevices.newReport[ReportType_Gamepad]) {
        Gamepad16.report();
        if(!forceReportAll) return;
    }

}

void LightgunButtons::ReleaseAll()
{
    AbsMouse5.releaseAll();
    Keyboard.releaseAll();
    Gamepad16.releaseAll();
    SendReports(true);
}

uint32_t LightgunButtons::Repeat()
{
    unsigned long m = millis();
    if(m - lastRepeatMillis >= interval) {
        lastRepeatMillis = m;
        repeat = debounced;
    } else {
        repeat = 0;
    }
    return repeat;
}

uint32_t LightgunButtons::PadMaskConvert(const uint32_t &mask)
{
    switch(mask) {
        case 1: // 0x00000001
            return GAMEPAD_HAT_UP;
            break;
        case 2: // 0x00000010
            return GAMEPAD_HAT_DOWN;
            break;
        case 4: // 0x00000100
            return GAMEPAD_HAT_LEFT;
            break;
        case 8: // 0x00001000
            return GAMEPAD_HAT_RIGHT;
            break;
        case 5: // 0x00000101
            return GAMEPAD_HAT_UP_LEFT;
            break;
        case 9: // 0x00001001
            return GAMEPAD_HAT_UP_RIGHT;
            break;
        case 6: // 0x00000110
            return GAMEPAD_HAT_DOWN_LEFT;
            break;
        case 10: // 0x00001010
            return GAMEPAD_HAT_DOWN_RIGHT;
            break;
        case 3: // 0x00000011
            return GAMEPAD_HAT_UP;
            break;
        case 12: // 0x00001100
            return GAMEPAD_HAT_LEFT;
            break;
        default: // 0x00000000
            return GAMEPAD_HAT_CENTERED;
            break;
    }
}

int LightgunButtons::MaskToIndex(uint32_t mask)
{
    uint32_t bitMask = 1;
    for(unsigned int i = 0; bitMask; ++i, bitMask <<= 1) {
        if(bitMask == mask) {
            return i;
        }
    }
    return -1;
}
