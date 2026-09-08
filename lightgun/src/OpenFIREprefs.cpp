/*!
 * @file OpenFIREprefs.cpp
 * @brief OpenFIRE file system loading/saving and presets access.
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V2.0
 * @date 2026
 *
 * I thank you for producing the first original code:
 *
 * @copyright Mike Lynch & That One Seong, 2021
 * @copyright GNU Lesser General Public License
 *
 * @author Mike Lynch
 * @author [That One Seong](SeongsSeongs@gmail.com)
 * @date 2025
 */

#include "OpenFIREprefs.h"

static constexpr size_t PREF_NAME_SIZE = 32;

static bool ReadPreferenceName(File &prefsFile, char *name, const size_t nameSize)
{
    if(name == nullptr || nameSize == 0)
        return false;

    size_t length = 0;
    while(true) {
        const int value = prefsFile.read();
        if(value < 0)
            return false;

        if(value == '\0') {
            name[length] = '\0';
            return length > 0;
        }

        if(length >= nameSize - 1)
            return false;

        name[length++] = (char)value;
    }
}

static bool SkipPreferenceData(File &prefsFile, size_t length)
{
    while(length > 0) {
        if(prefsFile.read() < 0)
            return false;
        --length;
    }

    return true;
}

/*
void OF_Prefs::InitProfileDefaults(const CameraProfile& profile)
{
    const float centerX = (float)profile.mouseResX * 0.5f;
    const float centerY = (float)profile.mouseResY * 0.5f;
    for (int i = 0; i < PROFILE_COUNT; ++i) {
        profiles[i].adjX = centerX;
        profiles[i].adjY = centerY;
    }
}
*/

/*
void OF_Prefs::EnsureProfileDefaults(const CameraProfile& profile)
{
    const float centerX = (float)profile.mouseResX * 0.5f;
    const float centerY = (float)profile.mouseResY * 0.5f;

    for (int i = 0; i < PROFILE_COUNT; ++i) {
        if (profiles[i].adjX == 0.0f)
            profiles[i].adjX = centerX;

        if (profiles[i].adjY == 0.0f)
            profiles[i].adjY = centerY;
    }
}
*/

void OF_Prefs::InitProfileDefaults(const CameraProfile& profile)
{
    const float centerX = (float)profile.mouseResX * 0.5f;
    const float centerY = (float)profile.mouseResY * 0.5f;

    for (int i = 0; i < PROFILE_COUNT; ++i) {
        if (profiles[i].adjX == 0.0f &&
            profiles[i].adjY == 0.0f) {
            profiles[i].adjX = centerX;
            profiles[i].adjY = centerY;
        }
    }
}

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
    if(!prefsFile)
        return Error_Read;

    char name[PREF_NAME_SIZE];
    bool loaded = true;

    while(loaded && prefsFile.available()) {
        loaded = ReadPreferenceName(prefsFile, name, sizeof(name));
        if(!loaded)
            break;

        const auto field = OFPresets.profSettingTypes_Strings.find(name);
        if(field != OFPresets.profSettingTypes_Strings.end() &&
           field->second == OF_Const::profCurrent) {
            const int profileNum = prefsFile.read();
            if(profileNum < 0) {
                loaded = false;
            } else {
                currentProfile = profileNum < PROFILE_COUNT ? profileNum : 0;
            }
            continue;
        }

        const int profileNum = prefsFile.read();
        const int storedSize = prefsFile.read();
        if(profileNum < 0 || storedSize < 0) {
            loaded = false;
            break;
        }

        if(field == OFPresets.profSettingTypes_Strings.end() ||
           profileNum >= PROFILE_COUNT) {
            loaded = SkipPreferenceData(prefsFile, (uint8_t)storedSize);
            continue;
        }

        const int index = field->second;
        if(index < 0 || index >= OF_Const::profDataTypes) {
            loaded = false;
            break;
        }

        const size_t expectedSize = index == OF_Const::profName ?
                                    sizeof(ProfileData_t::name) :
                                    sizeof(uint32_t);

        uint8_t value[sizeof(ProfileData_t::name)];
        if((size_t)storedSize != expectedSize ||
           prefsFile.readBytes((char*)value, expectedSize) != expectedSize) {
            loaded = false;
        } else {
            memcpy((uint8_t*)&profiles[profileNum] +
                       (sizeof(uint32_t) * index),
                   value,
                   expectedSize);

            if(index == OF_Const::profName)
                profiles[profileNum].name[sizeof(ProfileData_t::name) - 1] = '\0';
        }
    }

    prefsFile.close();
    return loaded ? Error_Success : Error_Read;
}

int OF_Prefs::SaveProfiles()
{
    File prefsFile = LittleFS.open("/profiles.conf", "w");
    if(!prefsFile)
        return Error_Write;

    bool written = true;
    bool currentProfLogged = false;
    for(size_t i = 0; i < PROFILE_COUNT && written; ++i) {
        for(auto &pair : OFPresets.profSettingTypes_Strings) {
            if(pair.second == OF_Const::profCurrent) {
                if(!currentProfLogged) {
                    // Same bytes and order: field name, then current profile.
                    written = prefsFile.write((const uint8_t*)pair.first.data(), pair.first.length()+1) == pair.first.length()+1 &&
                              prefsFile.write((uint8_t)currentProfile) == 1;
                    currentProfLogged = true;
                }
            } else {
                // Field name, profile number, byte count, value (unchanged format).
                written = prefsFile.write((const uint8_t*)pair.first.data(), pair.first.length()+1) == pair.first.length()+1 &&
                          prefsFile.write((uint8_t*)&i, 1) == 1;
                if(written) {
                    switch(pair.second) {
                    case OF_Const::profName:
                        written = prefsFile.write((uint8_t)sizeof(ProfileData_t::name)) == 1 &&
                                  prefsFile.write((uint8_t*)profiles[i].name, sizeof(ProfileData_t::name)) == sizeof(ProfileData_t::name);
                        break;
                    default:
                        written = prefsFile.write((uint8_t)sizeof(int)) == 1 &&
                                  prefsFile.write((uint8_t*)&profiles[i] + (sizeof(int)*pair.second), sizeof(int)) == sizeof(int);
                        break;
                    }
                }
            }
            if(!written)
                break;
        }
    }
    prefsFile.close();
    return written ? Error_Success : Error_Write;
}

int OF_Prefs::SaveToPtr(
    File prefsFile,
    void *dataPtr,
    const std::unordered_map<std::string_view, int> &mapPtr,
    const size_t &dataSize)
{
    if(!prefsFile)
        return Error_Write;

    const bool buttonData = dataPtr == backupButtonDesc;
    bool written = true;

    for(auto &pair : mapPtr) {
        if(pair.second < 0 ||
           (buttonData && pair.second >= ButtonCount))
            continue;

        written =
            prefsFile.write((const uint8_t*)pair.first.data(),
                            pair.first.length() + 1) == pair.first.length() + 1 &&
            prefsFile.write((uint8_t)dataSize) == 1 &&
            prefsFile.write((uint8_t*)dataPtr + (dataSize * pair.second),
                            dataSize) == dataSize;

        if(!written)
            break;
    }

    prefsFile.close();
    return written ? Error_Success : Error_Write;
}

int OF_Prefs::LoadToPtr(File prefsFile, void *dataPtr, const std::unordered_map<std::string_view, int> &mapPtr, const size_t &dataSize)
{
    if(!prefsFile)
        return Error_NoData;

    char name[PREF_NAME_SIZE];
    bool loaded = true;

    while(loaded && prefsFile.available()) {
        loaded = ReadPreferenceName(prefsFile, name, sizeof(name));
        if(!loaded)
            break;

        const int storedSize = prefsFile.read();
        if(storedSize < 0) {
            loaded = false;
            break;
        }

        const auto field = mapPtr.find(name);
        if(field == mapPtr.end() || field->second < 0 ||
           (dataPtr == backupButtonDesc && field->second >= ButtonCount)) {
            loaded = SkipPreferenceData(prefsFile, (uint8_t)storedSize);
            continue;
        }

        if((size_t)storedSize != dataSize ||
           (size_t)prefsFile.available() < dataSize ||
           prefsFile.readBytes((char*)dataPtr + (dataSize * field->second),
                               dataSize) != dataSize) {
            loaded = false;
        }
    }

    prefsFile.close();
    return loaded ? Error_Success : Error_Read;
}

int OF_Prefs::LoadUSBID()
{
    File idFile = LittleFS.open("/USB.conf", "r");
    if(!idFile)
        return Error_NoData;

    USBMap_t loadedUSB;
    const bool loaded = idFile.readBytes((char*)&loadedUSB,
                                         sizeof(loadedUSB)) == sizeof(loadedUSB);

    idFile.close();

    if(!loaded)
        return Error_Read;

    loadedUSB.deviceName[sizeof(loadedUSB.deviceName) - 1] = '\0';
    usb = loadedUSB;
    return Error_Success;
}

int OF_Prefs::SaveUSBID()
{
    File idFile = LittleFS.open("/USB.conf", "w");
    if(idFile) {
        const bool written = idFile.write((uint8_t*)&usb, sizeof(USBMap_t)) == sizeof(USBMap_t);

        idFile.close();
        return written ? Error_Success : Error_Write;
    } else return Error_Write;
}

void OF_Prefs::ResetPreferences()
{
    LittleFS.format();
}

void OF_Prefs::LoadPresets()
{
    memset(pins, -1, sizeof(OF_Prefs::pins));

    const auto boardIt =
        OFPresets.boardsPresetsMap.find(OPENFIRE_BOARD);

    if(boardIt != OFPresets.boardsPresetsMap.end()) {
        const auto &boardPins = boardIt->second;

        for(int i = 0; i < boardPins.size(); ++i) {
            const int pinFunction = boardPins.at(i);

            if(pinFunction > -1)
                pins[pinFunction] = i;
        }
    }

    // Initialize the backup button descriptors used when saving.
    for(int i = 0; i < ButtonCount; ++i)
        memcpy(OF_Prefs::backupButtonDesc[i],
               &LightgunButtons::ButtonDesc[i].reportType,
               sizeof(OF_Prefs::backupButtonDesc[i]));
}

#if defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)

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

int OF_Prefs::LoadLastPedalWireless(uint8_t *address, uint8_t *channel)
{
    File lastPedalFile = LittleFS.open("/lastPEDAL.conf", "r");
    if(lastPedalFile) {
        int bWritten = lastPedalFile.read(address, 6);
        bWritten += lastPedalFile.read(channel, 1);
        lastPedalFile.close();
        if (bWritten == (6 + 1)) return Error_Success; else return Error_NoData;
    } else return Error_NoData;
}

int OF_Prefs::SaveLastPedalWireless(uint8_t *address, uint8_t *channel)
{
    File lastPedalFile = LittleFS.open("/lastPEDAL.conf", "w");
    if(lastPedalFile) {
        int bWritten = lastPedalFile.write(address, 6);
        bWritten += lastPedalFile.write(channel, 1);
        lastPedalFile.close();
        if (bWritten == (6 + 1)) return Error_Success; else return Error_NoData;
    } else return Error_NoData;
}

#endif // defined(OPENFIRE_WIRELESS_ENABLE) && defined(ARDUINO_ARCH_ESP32)
