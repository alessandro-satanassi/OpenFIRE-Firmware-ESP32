/*!
 * @file OpenFIRE_Anchor.h
 * @brief Light Gun library for 4 LED setup
 * @n Ricostruzione dei LED nascosti, comune a Square_Advanced e Diamond_Advanced
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#ifndef _OpenFIRE_Anchor_h_
#define _OpenFIRE_Anchor_h_

#include <stdint.h>
#include <math.h>

// Limite numerico delle coordinate (lo stesso di Square_Advanced), non il campo visivo.
static constexpr float ANCHOR_COORD_LIMIT = 16383.0f;

// Ricostruzione dei LED nascosti da una geometria di riferimento (es. l'ultima vista con 4 LED),
// che ne conserva la forma reale (prospettiva compresa): una similitudine (n = 2 LED visti) o
// un'affinita' (n = 3) porta i LED visti allora (px, py) su quelli visti ora (qx, qy) e si applica
// agli m nascosti (hx, hy: in ingresso nella geometria di riferimento, in uscita ora; validi solo
// se torna true). Per non deformare
// la forma con la prospettiva, ciascuna delle due viste si guarda da una camera virtuale ruotata
// (rotazione minima) verso la direzione media dei suoi LED visti: se l'arma ha solo ruotato le
// due viste differiscono allora per una rotazione attorno all'asse e la ricostruzione e' esatta.
// Con focale approssimata o punti ancora distorti dalla lente, e' una stima;
// con focal = 0 si lavora direttamente sull'immagine.
// false anche se un nascosto esce dal limite numerico delle coordinate.
struct AnchorView {
    float f, ox, oy;                                 // focale e centro dell'immagine
    float r[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};        // rotazione della camera virtuale (per righe)

    AnchorView(float focal, float cx, float cy) : f(focal), ox(cx), oy(cy) {}

    bool aim(const float* x, const float* y, uint8_t n) {
        if (f <= 0.0f) return true;
        float mx = 0.0f, my = 0.0f, mz = 0.0f;       // somma delle direzioni (versori)
        for (uint8_t i = 0; i < n; ++i) {
            const float dx = x[i] - ox, dy = y[i] - oy;
            const float k = 1.0f / sqrtf(dx * dx + dy * dy + f * f);
            mx += dx * k; my += dy * k; mz += f * k;
        }
        const float k = 1.0f / sqrtf(mx * mx + my * my + mz * mz);
        mx *= k; my *= k; mz *= k;
        if (!(mz > 0.0f)) return false;
        // rotazione minima che porta (mx, my, mz) sull'asse ottico: asse (my, -mx, 0)
        const float vx = my, vy = -mx, c = 1.0f / (1.0f + mz);
        r[0] = 1.0f - vy * vy * c; r[1] = vx * vy * c;         r[2] = vy;
        r[3] = vx * vy * c;         r[4] = 1.0f - vx * vx * c; r[5] = -vx;
        r[6] = -vy;                 r[7] = vx;                 r[8] = 1.0f - (vx * vx + vy * vy) * c;
        return true;
    }
    // dall'immagine al piano della camera virtuale (unita' mouse) e ritorno
    bool to(float x, float y, float& u, float& v) const {
        const float dx = x - ox, dy = y - oy;
        if (f <= 0.0f) { u = dx; v = dy; return true; }
        const float wz = r[6] * dx + r[7] * dy + r[8] * f;
        if (!(wz > 1e-3f * f)) return false;
        const float k = f / wz;
        u = (r[0] * dx + r[1] * dy + r[2] * f) * k;
        v = (r[3] * dx + r[4] * dy + r[5] * f) * k;
        return true;
    }
    bool from(float u, float v, float& x, float& y) const {
        if (f <= 0.0f) { x = u + ox; y = v + oy; return true; }
        const float wz = r[2] * u + r[5] * v + r[8] * f;
        if (!(wz > 1e-3f * f)) return false;
        const float k = f / wz;
        x = ox + (r[0] * u + r[3] * v + r[6] * f) * k;
        y = oy + (r[1] * u + r[4] * v + r[7] * f) * k;
        return true;
    }
};

// Cyclic order of four strictly convex points, without atan2. The hull is exact
// for the bounded integer coordinates used by Square. Diamond uses measured floats.
template <typename T>
inline bool OrderFourPoints(const T* x, const T* y, uint8_t (&order)[4]) {
    uint8_t sorted[4] = {0,1,2,3};
    auto exchange = [&](uint8_t a, uint8_t b) {
        const uint8_t i=sorted[a], j=sorted[b];
        if (x[i]>x[j] || (x[i]==x[j] && y[i]>y[j])) { sorted[a]=j; sorted[b]=i; }
    };
    exchange(0,1); exchange(2,3); exchange(0,2); exchange(1,3); exchange(1,2);
    auto turn = [&](uint8_t a, uint8_t b, uint8_t c) {
        return (x[b]-x[a])*(y[c]-y[a]) - (y[b]-y[a])*(x[c]-x[a]);
    };
    uint8_t hull[8], count=0;
    for (uint8_t i=0;i<4;++i) {
        while (count>=2 && turn(hull[count-2],hull[count-1],sorted[i])<=0) --count;
        hull[count++]=sorted[i];
    }
    const uint8_t lower=count;
    for (int8_t i=2;i>=0;--i) {
        while (count>lower && turn(hull[count-2],hull[count-1],sorted[i])<=0) --count;
        hull[count++]=sorted[i];
    }
    if (count!=5) return false;
    // Same cyclic origin as atan2 in [-pi, pi]: first point above the left ray.
    const float cy=((float)y[0]+y[1]+y[2]+y[3])*0.25f;
    uint8_t first=0;
    for (uint8_t i=0;i<4;++i)
        if (y[hull[i]]<cy && y[hull[(i+3)&3]]>=cy) {first=i;break;}
    for (uint8_t i=0;i<4;++i) order[i]=hull[(first+i)&3];
    return true;
}

// A physical parallelogram has one shared midpoint for its diagonals. Their
// image intersection yields relative depths, so lengths can be compared in 3D
// without a homography or square roots. Raw lens distortion makes this approximate.
struct CameraQuad {
    float x[4], y[4], z[4];

    template <typename T>
    bool project(float focal, float cx, float cy, const T* px, const T* py,
                 const uint8_t (&order)[4]) {
        if (!(focal>0.0f)) return false;
        const float ax=px[order[0]], ay=py[order[0]], bx=px[order[1]], by=py[order[1]];
        const float ux=px[order[2]]-ax, uy=py[order[2]]-ay;
        const float vx=px[order[3]]-bx, vy=py[order[3]]-by;
        const float det=ux*vy-uy*vx;
        if (!(fabsf(det)>1.0f)) return false;
        const float inv=1.0f/det, dx=bx-ax, dy=by-ay;
        const float t=(dx*vy-dy*vx)*inv, u=(dx*uy-dy*ux)*inv;
        if (!(t>0.0f && t<1.0f && u>0.0f && u<1.0f)) return false;
        const float depth[4]={1.0f-t,1.0f-u,t,u};
        for (uint8_t i=0;i<4;++i) {
            x[i]=(px[order[i]]-cx)*depth[i];
            y[i]=(py[order[i]]-cy)*depth[i];
            z[i]=focal*depth[i];
        }
        return true;
    }

    float distance2(uint8_t a, uint8_t b) const {
        const float dx=x[b]-x[a],dy=y[b]-y[a],dz=z[b]-z[a];
        return dx*dx+dy*dy+dz*dz;
    }
};


// Cache only the reference geometry. A repeated partial view needs to rotate/project
// the current rays, but need not recompute the unchanged complete reference.
// The key includes every reference and hidden point: no caller-side invalidation is needed.
class AnchorReconstruction {
public:
    void configure(float focal, float cx, float cy) {
        f = focal; ox = cx; oy = cy; ready = false;
    }

    bool reconstruct(uint8_t n, const float* px, const float* py,
                     const float* qx, const float* qy, uint8_t m, float* hx, float* hy) {
        if (!matches(n, px, py, m, hx, hy) && !prepare(n, px, py, m, hx, hy)) return false;
        AnchorView Q(f, ox, oy);
        if (!Q.aim(qx, qy, n)) return false;
        float qu[3] = {}, qv[3] = {};
        for (uint8_t i = 0; i < n; ++i)
            if (!Q.to(qx[i], qy[i], qu[i], qv[i])) return false;
        const float b1x = qu[1] - qu[0], b1y = qv[1] - qv[0];
        const float b2x = n == 2 ? -b1y : qu[2] - qu[0];
        const float b2y = n == 2 ?  b1x : qv[2] - qv[0];
        if (n == 3) {
            const float p = (b1x * a2y - b2x * a1y) * inv;
            const float q = (b2x * a1x - b1x * a2x) * inv;
            const float r = (b1y * a2y - b2y * a1y) * inv;
            const float t = (b2y * a1x - b1y * a2x) * inv;
            const float cs = (p+t)*(p+t) + (r-q)*(r-q);
            const float ac = (p-t)*(p-t) + (q+r)*(q+r);
            if (!(cs > 0.0f && ac * shearScale <= cs * shearLimit)) return false;
        }
        for (uint8_t j = 0; j < m; ++j) {
            if (!Q.from(qu[0] + hs[j]*b1x + ht[j]*b2x,
                        qv[0] + hs[j]*b1y + ht[j]*b2y, hx[j], hy[j])) return false;
            if (!(fabsf(hx[j]) <= ANCHOR_COORD_LIMIT && fabsf(hy[j]) <= ANCHOR_COORD_LIMIT)) return false;
        }
        return true;
    }

private:
    float f = 0.0f, ox = 0.0f, oy = 0.0f;
    bool ready = false;
    uint8_t count = 0, hiddenCount = 0;
    float refX[3], refY[3], hiddenX[2], hiddenY[2];
    float a1x, a1y, a2x, a2y, inv, shearScale, shearLimit;
    float hs[2], ht[2];

    bool matches(uint8_t n, const float* x, const float* y,
                 uint8_t m, const float* hx, const float* hy) const {
        if (!ready || count != n || hiddenCount != m) return false;
        for (uint8_t i = 0; i < n; ++i)
            if (refX[i] != x[i] || refY[i] != y[i]) return false;
        for (uint8_t i = 0; i < m; ++i)
            if (hiddenX[i] != hx[i] || hiddenY[i] != hy[i]) return false;
        return true;
    }

    bool prepare(uint8_t n, const float* x, const float* y,
                 uint8_t m, const float* hx, const float* hy) {
        ready = false;
        AnchorView P(f, ox, oy);
        if (!P.aim(x, y, n)) return false;
        float pu[3] = {}, pv[3] = {};
        for (uint8_t i = 0; i < n; ++i)
            if (!P.to(x[i], y[i], pu[i], pv[i])) return false;
        a1x = pu[1]-pu[0]; a1y = pv[1]-pv[0];
        a2x = n == 2 ? -a1y : pu[2]-pu[0];
        a2y = n == 2 ?  a1x : pv[2]-pv[0];
        const float l1 = a1x*a1x + a1y*a1y, l2 = a2x*a2x + a2y*a2y;
        const float det = a1x*a2y - a1y*a2x;
        if (!(fabsf(det) > 1e-3f*(l1+l2))) return false;
        inv = 1.0f/det;
        const float lmin = fminf(l1,l2), lmax = fmaxf(l1,l2);
        shearScale = 8.0f*lmax + 2.0f*lmin; shearLimit = lmin;
        for (uint8_t j = 0; j < m; ++j) {
            float u, v;
            if (!P.to(hx[j], hy[j], u, v)) return false;
            const float ex = u-pu[0], ey = v-pv[0];
            hs[j] = (ex*a2y-ey*a2x)*inv; ht[j] = (a1x*ey-a1y*ex)*inv;
            hiddenX[j] = hx[j]; hiddenY[j] = hy[j];
        }
        for (uint8_t i = 0; i < n; ++i) { refX[i]=x[i]; refY[i]=y[i]; }
        count = n; hiddenCount = m; ready = true;
        return true;
    }
};

// Stateless entry point retained for isolated geometry users.
inline bool LedsFromAnchor(float focal, float ox, float oy, uint8_t n,
                           const float* px, const float* py, const float* qx, const float* qy,
                           uint8_t m, float* hx, float* hy) {
    AnchorReconstruction reference;
    reference.configure(focal, ox, oy);
    return reference.reconstruct(n, px, py, qx, qy, m, hx, hy);
}

#endif // _OpenFIRE_Anchor_h_
