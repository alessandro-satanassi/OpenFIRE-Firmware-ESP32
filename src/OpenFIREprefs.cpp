/*!
 * @file OpenFIREprefs.cpp
 * @brief OpenFIRE file system loading/saving and presets access.
 *
 * @copyright Mike Lynch & That One Seong, 2021
 * @copyright GNU Lesser General Public License
 *
 * @author Mike Lynch
 * @author [That One Seong](SeongsSeongs@gmail.com)
 * @date 2025
 */

#include "OpenFIREprefs.h"

int OF_Prefs::InitFS()
{
    #ifdef ARDUINO_ARCH_ESP32  
    if(LittleFS.begin(true))
    #else
    if(LittleFS.begin())
    #endif
        return Error_Success;
    else return Error_NoData;
}

void OF_Prefs::Load()
{
    LoadToggles();
    if(toggles[OF_Const::customPins]) LoadPins();
    LoadSettings();
    LoadUSBID();
    LoadButtons();
    for(int i = 0; i < ButtonCount; ++i)
        memcpy(&LightgunButtons::ButtonDesc[i].reportType, OF_Prefs::backupButtonDesc[i], sizeof(OF_Prefs::backupButtonDesc[i]));
}

int OF_Prefs::LoadProfiles()
{
    File prefsFile = LittleFS.open("/profiles.conf", "r");
    if(prefsFile) {
        int profileNum = 0;
        char buf[32];
        size_t bWritten = 0;
        size_t readSize = 0;
        while(prefsFile.available()) {
            bWritten = prefsFile.readBytesUntil('\0', buf, 32);
            // readBytesUntil discards the terminator, so plop one at the end
            buf[bWritten++] = '\0';
            if(bWritten && OFPresets.profSettingTypes_Strings.count(buf)) {
                switch(OFPresets.profSettingTypes_Strings.at(buf)) {
                  case OF_Const::profCurrent:
                      currentProfile = prefsFile.read();
                      if(currentProfile >= PROFILE_COUNT) currentProfile = 0;
                      break;
                  default:
                      profileNum = prefsFile.read();
                      readSize = prefsFile.read();
                      profileNum < PROFILE_COUNT ? prefsFile.readBytes((char*)&profiles[profileNum] + (sizeof(uint32_t) * OFPresets.profSettingTypes_Strings.at(buf)), readSize) : prefsFile.seek(readSize, fs::SeekCur);
                      break;
                }
            } else {
                prefsFile.seek(1, fs::SeekCur);
                readSize = prefsFile.read();
                prefsFile.seek(readSize, fs::SeekCur);
            }
        }

        prefsFile.close();
        return Error_Success;
    } else return Error_Read;
}

int OF_Prefs::SaveProfiles()
{
    File prefsFile = LittleFS.open("/profiles.conf", "w");
    if(prefsFile) {
        bool currentProfLogged = false;
        for(size_t i = 0; i < PROFILE_COUNT; ++i) {
            for(auto &pair : OFPresets.profSettingTypes_Strings) {
                if(pair.second == OF_Const::profCurrent) {
                    if(!currentProfLogged) {
                        // only write string and profile num
                        prefsFile.write((const uint8_t*)pair.first.data(), pair.first.length()+1);
                        prefsFile.write((uint8_t)currentProfile);
                        currentProfLogged = true;
                    }
                } else {
                    // write data type:
                    prefsFile.write((const uint8_t*)pair.first.data(), pair.first.length()+1);
                    // Append profile number:
                    prefsFile.write((uint8_t*)&i, 1);

                    // data type:
                    switch(pair.second) {
                    // 16-bytes profile name
                    case OF_Const::profName:
                        prefsFile.write(sizeof(ProfileData_t::name));
                        prefsFile.write((uint8_t*)profiles[i].name, sizeof(ProfileData_t::name));
                        break;
                    // everything else is generic 32-bit data
                    default:
                        prefsFile.write(sizeof(int));
                        prefsFile.write((uint8_t*)&profiles[i] + (sizeof(int)*pair.second), sizeof(int));
                        break;
                    }
                }
            }
        }

        prefsFile.close();
        return Error_Success;
    } else return Error_Write;
}

int OF_Prefs::SaveToPtr(File prefsFile, void *dataPtr, const std::unordered_map<std::string_view, int> &mapPtr, const size_t &dataSize)
{
    if(prefsFile) {
        for(auto &pair : mapPtr) {
            if((pair.second >= 0 && dataPtr != backupButtonDesc) || dataPtr == backupButtonDesc && pair.second >= 0 && pair.second < ButtonCount) {
                prefsFile.write((const uint8_t*)pair.first.data(), pair.first.length()+1);
                prefsFile.write((uint8_t)dataSize);
                prefsFile.write((uint8_t*)dataPtr + (dataSize * pair.second), dataSize);
            }
        }
        
        prefsFile.close();
        return Error_Success;
    } else return Error_NoData;
}

int OF_Prefs::LoadToPtr(File prefsFile, void *dataPtr, const std::unordered_map<std::string_view, int> &mapPtr)
{
    if(prefsFile) {
        char buf[32];
        size_t bWritten = 0;
        size_t dataSize = 0;
        while(prefsFile.available()) {
            bWritten = prefsFile.readBytesUntil('\0', buf, 32);
            // readBytesUntil discards the terminator, so plop one at the end
            buf[bWritten++] = '\0';
            dataSize = prefsFile.read();
            if(bWritten && mapPtr.count(buf))
                prefsFile.readBytes((char*)dataPtr + (dataSize * mapPtr.at(buf)), dataSize);
            else prefsFile.seek(dataSize, fs::SeekCur);
        }

        prefsFile.close();
        return Error_Success;
    } else return Error_NoData;
}

int OF_Prefs::LoadUSBID()
{
    File idFile = LittleFS.open("/USB.conf", "r");
    if(idFile) {
        idFile.readBytes((char*)&usb, sizeof(USBMap_t));

        idFile.close();
        return Error_Success;
    } else return Error_NoData;
}

int OF_Prefs::SaveUSBID()
{
    File idFile = LittleFS.open("/USB.conf", "w");
    if(idFile) {
        idFile.write((uint8_t*)&usb, sizeof(USBMap_t));

        idFile.close();
        return Error_Success;
    } else return Error_NoData;
}

void OF_Prefs::ResetPreferences()
{
    LittleFS.format();
}

void OF_Prefs::LoadPresets()
{
    memset(pins, -1, sizeof(OF_Prefs::pins));

    if(OFPresets.boardsPresetsMap.count(OPENFIRE_BOARD)) {
        for(int i = 0; i < OFPresets.boardsPresetsMap.at(OPENFIRE_BOARD).size(); ++i)
            if(OFPresets.boardsPresetsMap.at(OPENFIRE_BOARD).at(i) > -1)
                pins[OFPresets.boardsPresetsMap.at(OPENFIRE_BOARD).at(i)] = i;
    }

    // save buttons map to backup descriptor (important for initializing backupDesc for saving)
    // couldn't find a better place for this tbh lol
    for(int i = 0; i < ButtonCount; ++i)
        memcpy(OF_Prefs::backupButtonDesc[i], &LightgunButtons::ButtonDesc[i].reportType, sizeof(OF_Prefs::backupButtonDesc[i]));
}

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)

int OF_Prefs::LoadWireless(uint8_t *channel, uint8_t *power)
{
    File wirelessFile = LittleFS.open("/wireless.conf", "r");
    if(wirelessFile) {
      //if (wirelessFile.available()) {
        int bWritten = wirelessFile.read(channel, sizeof(uint8_t));
        bWritten += wirelessFile.read(power, sizeof(uint8_t));
        wirelessFile.close();
        if (bWritten == 2) return Error_Success; else return Error_NoData;
      //}
    } else return Error_NoData;
}

int OF_Prefs::SaveWireless(uint8_t *channel, uint8_t *power)
{
    File wirelessFile = LittleFS.open("/wireless.conf", "w");
    if(wirelessFile) {
        int bWritten = wirelessFile.write(channel, sizeof(int8_t));
        bWritten += wirelessFile.write(power, sizeof(int8_t));
        wirelessFile.close();
        if (bWritten == 2) return Error_Success; else return Error_NoData;
    } else return Error_NoData;
}

int OF_Prefs::LoadLastDongleWireless(uint8_t *address, uint8_t *channel)
{
    File lastDongleFile = LittleFS.open("/lastDONGLE.conf", "r");
    if(lastDongleFile) {
        int bWritten = lastDongleFile.read(address, 6);
        bWritten += lastDongleFile.read(channel, 1);
        lastDongleFile.close();
        if (bWritten == (6 + 1)) return Error_Success; else return Error_NoData;
    } else return Error_NoData;
}

int OF_Prefs::SaveLastDongleWireless(uint8_t *address, uint8_t *channel)
{
    File lastDongleFile = LittleFS.open("/lastDONGLE.conf", "w");
    if(lastDongleFile) {
        int bWritten = lastDongleFile.write(address, 6);
        bWritten += lastDongleFile.write(channel, 1);
        lastDongleFile.close();
        if (bWritten == (6 + 1)) return Error_Success; else return Error_NoData;
    } else return Error_NoData;
}

#endif // defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)