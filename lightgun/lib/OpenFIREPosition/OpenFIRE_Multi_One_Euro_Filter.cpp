#ifdef USE_MULTI_ONE_EURO_FILTER
/*!
 * @file OpenFIRE_Multi_One_Euro_Filter.cpp
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

#include "OpenFIRE_Multi_One_Euro_Filter.h"

// EN/IT: Paired weights, one shared reciprocal / pesi accoppiati, un solo reciproco.
// EN: Preserve operation order to retain the existing floating-point rounding.
// IT: L'ordine delle operazioni conserva gli arrotondamenti del filtro.
static inline void pairedAlpha(float fx, float fy, float dtTwoPi, float& ax, float& ay) {
    const float tx = fx * dtTwoPi, ty = fy * dtTwoPi;
    const float denX = 1.0f + tx, denY = 1.0f + ty;
    const float inverse = 1.0f / (denX * denY);
    ax = (tx * denY) * inverse;
    ay = (ty * denX) * inverse;
}

OpenFIRE_One_Euro_Multi::OpenFIRE_One_Euro_Multi() {
    lastMicros = 0;
    initialized = false;
    for (int i = 0; i < 4; i++) {
        states[i] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    }
}

void OpenFIRE_One_Euro_Multi::configure(const CameraProfile& profile) {
    // EN/IT: Profile-dependent values, calculated at setup / valori del profilo calcolati all'avvio.
    invNoiseX = 1.0f / (((float)profile.mouseResX / profile.sensorResX) * profile.noiseFactor);
    invNoiseY = 1.0f / (((float)profile.mouseResY / profile.sensorResY) * profile.noiseFactor);

    // EN: Ordinary optics remain neutral; wide-angle adaptation uses profile data only.
    // IT: Le ottiche ordinarie restano neutre; l'adattamento grandangolare usa solo il profilo.
    opticalMotionRange = 0.0f;
    opticalSlow = opticalStability = 1.0f;
    opticalActive = false;
    if (profile.lensFocalRatio > 0.0f && profile.sensorResX > 0) {
        float side = 1.0f + profile.lensRadialK1 + profile.lensRadialK2;
        if (side < profile.lensCorrectionMin)
            side = profile.lensCorrectionMin;
        if (side > profile.lensCorrectionMax)
            side = profile.lensCorrectionMax;
        const float tangent = side / profile.lensFocalRatio;

        float excess = tangent * tangent - opticalNeutralTangentSq;
        if (excess < 0.0f)
            excess = 0.0f;
        const float weight = excess / (1.0f + excess);
        const float angular = opticalReferenceFocal / (0.5f * profile.sensorResX * profile.lensFocalRatio);
        const float angularExcess = angular > 1.0f ? angular - 1.0f : 0.0f;
        const float motion = 1.0f + opticalMotionStrength * weight * angularExcess;
        opticalStability = 1.0f + opticalStabilityStrength * excess;
        opticalSlow = 1.0f / opticalStability;

        opticalMotionRange = motion - opticalSlow;
        opticalActive = excess > 0.0f;
    }

    mouseMaxX = profile.mouseMaxX;
    mouseMaxY = profile.mouseMaxY;
    centerX = (float)mouseMaxX * 0.5f;
    centerY = (float)mouseMaxY * 0.5f;
    invCenterX = 1.0f / centerX;
    invCenterY = 1.0f / centerY;
    fallbackDt = 1.0f / (float)profile.fps;

    initialized = false;
    lastMicros = 0;
}

template <bool Optical>
__attribute__((always_inline)) inline void OpenFIRE_One_Euro_Multi::filterFrame(
    int* x, int* y, float* outX, float* outY) {
    // 1. EN/IT: Elapsed time and restart / tempo trascorso e ripresa del tracking.
    // EN/IT: Unsigned subtraction handles clock wrap / sottrazione unsigned per l'overflow del clock.
    unsigned long currentMicros = micros();
    float dt = ((float)(currentMicros - lastMicros)) * 0.000001f;
    lastMicros = currentMicros;

    if (dt <= 0.0f)
        dt = fallbackDt;

    if (dt > 0.1f) {
        initialized = false;
    }

    if (!initialized) {
        #pragma GCC unroll 2
        for (int i = 0; i < 4; i++) {
            states[i].previousX = states[i].filteredX = (float)x[i];
            states[i].previousY = states[i].filteredY = (float)y[i];
            states[i].velocityX = states[i].velocityY = 0.0f;
            outX[i] = states[i].filteredX;
            outY[i] = states[i].filteredY;
        }
        initialized = true;
        return;
    }

    // 2. EN/IT: Edge response and derivative weights / risposta ai bordi e pesi delle derivate.
    const float dtTwoPi = dt * twoPi;

    float baseAlpha, snapAlpha;
    pairedAlpha(derivativeBaseCutoff, derivativeSnapCutoff, dtTwoPi, baseAlpha, snapAlpha);

    const float invDt = 1.0f / dt;

    // EN: The most central in-image point guides edge attenuation, with a 0.20 floor.
    // IT: Il punto piu centrale nell'immagine guida l'attenuazione, con minimo 0,20.
    float edge = 0.20f;

    #pragma GCC unroll 4
    for (int i = 0; i < 4; i++) {
        if ((uint32_t)x[i] <= (uint32_t)mouseMaxX && (uint32_t)y[i] <= (uint32_t)mouseMaxY) {
            float edgeX = ((float)x[i] - centerX) * invCenterX;
            float edgeY = ((float)y[i] - centerY) * invCenterY;

            float distanceSquared = (edgeX * edgeX) + (edgeY * edgeY);
            float attenuation = 1.0f - (distanceSquared * 0.8f);

            if (attenuation > edge) {
                edge = attenuation;
            }
        }
    }

    float adaptiveBeta = betaGain * edge;
    const float dynamicSnap = snapGain + snapEdgeGain * (1.0f - edge);
    const float lower = dynamicSnap * 0.8f, upper = dynamicSnap * 1.2f;
    const float invRange = 1.0f / (dynamicSnap * 0.4f);
    const float invDtX = invDt * invNoiseX;
    const float invDtY = invDt * invNoiseY;

    // 3. EN: Individual velocities and shared motion statistics.
    // IT: Velocita individuali e statistiche di movimento comuni.
    // EN: Explicit X/Y branches avoid extra stack pressure in the hot loop.
    // IT: I rami X/Y espliciti evitano stack aggiuntivo nel ciclo veloce.
    float maxVelocityX = 0.0f;
    float maxVelocityY = 0.0f;
    float energyX = 0.0f, energyY = 0.0f, alignmentX = 0.0f, alignmentY = 0.0f;
    float lagEnergy = 0.0f;

    #pragma GCC unroll 2
    for (int i = 0; i < 4; i++) {
        float dx = ((float)x[i] - states[i].previousX) * invDtX;
        float differenceX = fabsf(dx - states[i].velocityX);

        float derivativeAlphaX;
        if (differenceX <= lower) {
            derivativeAlphaX = baseAlpha;
        } else if (differenceX >= upper) {
            derivativeAlphaX = snapAlpha;
        } else {
            derivativeAlphaX = baseAlpha + ((differenceX - lower) * invRange) * (snapAlpha - baseAlpha);
        }

        states[i].velocityX += derivativeAlphaX * (dx - states[i].velocityX);

        float absVelocityX = fabsf(states[i].velocityX);
        if (absVelocityX > maxVelocityX)
            maxVelocityX = absVelocityX;

        float dy = ((float)y[i] - states[i].previousY) * invDtY;
        float differenceY = fabsf(dy - states[i].velocityY);

        float derivativeAlphaY;
        if (differenceY <= lower) {
            derivativeAlphaY = baseAlpha;
        } else if (differenceY >= upper) {
            derivativeAlphaY = snapAlpha;
        } else {
            derivativeAlphaY = baseAlpha + ((differenceY - lower) * invRange) * (snapAlpha - baseAlpha);
        }

        states[i].velocityY += derivativeAlphaY * (dy - states[i].velocityY);

        float absVelocityY = fabsf(states[i].velocityY);
        if (absVelocityY > maxVelocityY)
            maxVelocityY = absVelocityY;
        if (Optical) {
            const float ex = ((float)x[i] - states[i].filteredX) * invNoiseX;
            const float ey = ((float)y[i] - states[i].filteredY) * invNoiseY;
            lagEnergy += ex * ex + ey * ey;
        }
        energyX += states[i].velocityX * states[i].velocityX;
        energyY += states[i].velocityY * states[i].velocityY;
        alignmentX += dx * states[i].velocityX;
        alignmentY += dy * states[i].velocityY;
    }

    // 4. EN: Optical gain and braking share one reciprocal.
    // IT: Guadagno ottico e frenata condividono un solo reciproco.
    float opticalBrake = 0.0f;
    if (Optical) {
        const float energy = energyX + energyY, alignment = alignmentX + alignmentY;
        const float regularizer = 2.0f * dynamicSnap * dynamicSnap;
        const float radius = dynamicSnap * fallbackDt * opticalStability * opticalLagRadius;
        const float noise = 8.0f * radius * radius;
        const float brakeDen = energy + fabsf(alignment) + regularizer;
        const float lagDen = lagEnergy + noise;

        const float inverse = 1.0f / (brakeDen * lagDen);
        opticalBrake = ((energy > alignment ? energy - alignment : 0.0f) * lagDen) * inverse;
        const float activity = (lagEnergy * brakeDen) * inverse;
        adaptiveBeta *= opticalSlow + opticalMotionRange * activity;
    }

    // 5. EN/IT: Common position weights and braking / pesi di posizione comuni e frenata.
    // EN: Sum pointwise energies/dots so rotational velocities do not cancel each other.
    // IT: Sommare energie e prodotti dei singoli punti evita di cancellare le velocita nelle rotazioni.
    float cutoffX = minCutoff + adaptiveBeta * maxVelocityX;
    if (cutoffX > maxCutoff)
        cutoffX = maxCutoff;
    float alphaX, alphaY;

    float cutoffY = minCutoff + adaptiveBeta * maxVelocityY;
    if (cutoffY > maxCutoff)
        cutoffY = maxCutoff;
    pairedAlpha(cutoffX, cutoffY, dtTwoPi, alphaX, alphaY);
    if (alphaX != alphaY) {
        const float energy = energyX + energyY;
        const float alignment = alignmentX + alignmentY;
        const float regularizer = 2.0f * dynamicSnap * dynamicSnap;
        const float reduction = energy > alignment ? energy - alignment : 0.0f;
        float weight = Optical ? opticalBrake : reduction / (energy + fabsf(alignment) + regularizer);
        weight *= weight;
        const float commonAlpha = alphaX > alphaY ? alphaX : alphaY;
        alphaX += weight * (commonAlpha - alphaX);
        alphaY += weight * (commonAlpha - alphaY);
    }

    // 6. EN: Apply the shared weights to each point's history.
    // IT: Applichiamo i pesi comuni alla storia di ciascun punto.
    #pragma GCC unroll 2
    for (int i = 0; i < 4; i++) {
        // EN: Snap only a fractional residual when the raw coordinate is unchanged.
        // IT: Eliminiamo solo un residuo frazionario quando la coordinata raw e ferma.
        bool stillX = ((float)x[i] == states[i].previousX);
        bool stillY = ((float)y[i] == states[i].previousY);

        states[i].filteredX += alphaX * ((float)x[i] - states[i].filteredX);
        states[i].previousX = (float)x[i];

        states[i].filteredY += alphaY * ((float)y[i] - states[i].filteredY);
        states[i].previousY = (float)y[i];

        if (stillX && fabsf((float)x[i] - states[i].filteredX) < microSnapThreshold) {
            states[i].filteredX = (float)x[i];
        }
        if (stillY && fabsf((float)y[i] - states[i].filteredY) < microSnapThreshold) {
            states[i].filteredY = (float)y[i];
        }

        if (isnan(states[i].filteredX) || isnan(states[i].filteredY)) {
            states[i].filteredX = (float)x[i];
            states[i].filteredY = (float)y[i];
            states[i].velocityX = 0.0f;
            states[i].velocityY = 0.0f;
        }

        outX[i] = states[i].filteredX;
        outY[i] = states[i].filteredY;
    }
}

// EN: Keep the template above this wrapper so both optical paths compile inline.
// IT: Il template precede questo wrapper per espandere i due percorsi senza chiamate interne.
void OpenFIRE_One_Euro_Multi::process(int* x, int* y, float* outX, float* outY) {
    if (opticalActive)
        filterFrame<true>(x, y, outX, outY);
    else
        filterFrame<false>(x, y, outX, outY);
}

#endif
