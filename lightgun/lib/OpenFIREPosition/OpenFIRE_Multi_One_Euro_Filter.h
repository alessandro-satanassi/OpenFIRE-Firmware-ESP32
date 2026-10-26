#ifdef USE_MULTI_ONE_EURO_FILTER
/*!
 * @file OpenFIRE_Multi_One_Euro_Filter.h
 * @brief Library for One Euro Filter for 4 LED
 * @n CPP Library for One Euro Filter for 4 LED
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V2.0
 * @date 2026
 */
#ifndef OpenFIRE_Multi_One_Euro_Filter_h
#define OpenFIRE_Multi_One_Euro_Filter_h

#include <Arduino.h>
#include "OpenFIREConst.h"
#include <OpenFIRECameraProfile.h>

class OpenFIRE_One_Euro_Multi {
private:
    // EN: Each point retains its own history; all points share the final X/Y weights.
    // IT: Ogni punto ha la propria storia; i pesi finali X/Y sono comuni ai quattro punti.
    struct FilterState {
        float previousX, filteredX;
        float previousY, filteredY;
        float velocityX, velocityY;
    };
    FilterState states[4];
    unsigned long lastMicros;
    bool initialized = false;

    // EN/IT: Temporal tuning / taratura temporale, invariata.
    static constexpr float minCutoff = 0.1f;
    static constexpr float maxCutoff = 30.0f;
    static constexpr float derivativeBaseCutoff = 1.0f;
    static constexpr float derivativeSnapCutoff = 25.0f;
    static constexpr float microSnapThreshold = 0.5f;
    static constexpr float twoPi = 6.28318530718f;

    // EN: Velocities are normalized by the sensor's spatial/noise units.
    // IT: Le velocita sono normalizzate sulle unita spaziali/di rumore del sensore.
    static constexpr float snapGain = 31.25f;
    static constexpr float snapEdgeGain = 62.5f;
    static constexpr float betaGain = 0.2464f;
    float invNoiseX = 1.0f;
    float invNoiseY = 1.0f;

    // EN/IT: Geometry cached by configure() / geometria precalcolata in configure().
    float centerX = 0.0f, centerY = 0.0f;
    float invCenterX = 0.0f, invCenterY = 0.0f;
    int mouseMaxX = 0, mouseMaxY = 0;

    // EN: Optical tuning uses profile data only, without camera model checks.
    // IT: La taratura ottica usa solo i dati del profilo, senza riconoscere la camera.
    static constexpr float opticalReferenceFocal = 216.0629425048828125f;
    static constexpr float opticalNeutralTangentSq = 0.25f;
    static constexpr float opticalMotionStrength = 0.25f;
    static constexpr float opticalStabilityStrength = 0.25f;
    static constexpr float opticalLagRadius = 0.6f;
    float opticalMotionRange = 0.0f;
    float opticalSlow = 1.0f, opticalStability = 1.0f;
    bool opticalActive = false;
    float fallbackDt = 1.0f / 209.0f;

    template <bool Optical> inline void filterFrame(int* x, int* y, float* outX, float* outY);

public:
    OpenFIRE_One_Euro_Multi();
    // EN: Restart history after a profile/layout change.
    // IT: Riavvia la storia dopo un cambio di profilo/layout.
    void configure(const CameraProfile& profile);
    // EN: Four ordered points in mouse coordinates; subpixel output.
    // IT: Quattro punti ordinati in coordinate mouse; uscita subpixel.
    void process(int* x, int* y, float* outX, float* outY);
};

#endif // OpenFIRE_Multi_One_Euro_Filter_h
#endif // USE_MULTI_ONE_EURO_FILTER
