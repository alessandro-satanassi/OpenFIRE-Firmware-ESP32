#include "OpenFIRE_Calibration.h"
#include "OpenFIREConst.h"
#include "OpenFIRE_Perspective_Advanced.h"
#include <algorithm>
#include <cmath>

#ifdef USE_PERSPECTIVE_ADVANCED
namespace OpenFIRECalibration {
namespace {

// Calibration targets as fractions of the screen: centre, top, bottom, left, right, centre.
constexpr float TargetX[6] = {.5f, .5f, .5f, 0.f, 1.f, .5f};
constexpr float TargetY[6] = {.5f, 0.f, 1.f, .5f, .5f, .5f};
using Shots = std::array<Measurement, 6>;
constexpr unsigned ShotCount = std::tuple_size<Shots>::value;
// Screen parameters fitted from the shots: offset and scale on each axis.
constexpr unsigned FitParameters = 4;
// 2 ln 2: median of a squared two-dimensional normal error, in units of its per-axis variance.
constexpr float SquaredError2DMedian = 1.3862944f;

struct Point { float x = 0, y = 0; };
struct Line { float x = 0, y = 0, z = 0; };

float square(float x) { return x*x; }
Line cross(Line a, Line b) { return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x}; }

// Same radial model as OpenFIRE_Perspective (camera -> corrected) and its inverse.
float lensGain(const CameraProfile& p, float r2) {
    return std::clamp(1.f+p.lensRadialK1*r2+p.lensRadialK2*r2*r2, p.lensCorrectionMin, p.lensCorrectionMax);
}
Point corrected(const CameraProfile& p, Point q) {
    const float cx = p.mouseResX*.5f, cy = p.mouseResY*.5f;
    const float dx = q.x-cx, dy = q.y-cy, gain = lensGain(p, square(dx/cx)+square(dy/cy));
    return {cx+dx*gain, cy+dy*gain};
}
Point distorted(const CameraProfile& p, Point q) {
    const float cx = p.mouseResX*.5f, cy = p.mouseResY*.5f;
    const float dx = q.x-cx, dy = q.y-cy, radius = std::sqrt(square(dx/cx)+square(dy/cy));
    if(!(radius > 0) || (p.lensRadialK1 == 0 && p.lensRadialK2 == 0)) return q;
    // r*gain(r) grows with r inside the clamp limits: bisection to float resolution.
    float low = 0, high = radius/p.lensCorrectionMin;
    for(unsigned i = 0; i < 24; ++i) {
        const float r = (low+high)*.5f;
        if(r*lensGain(p, r*r) < radius) low = r; else high = r;
    }
    const float gain = (low+high)*.5f/radius;
    return {cx+dx*gain, cy+dy*gain};
}

// The two lines through the LED centre that the barrel crosses when aiming at the
// screen's vertical (axis 0) and horizontal (axis 1) centre lines.
void midlines(const CameraProfile& p, bool diamond, const Measurement& shot, Line lines[2]) {
    Line v[4];
    for(unsigned j = 0; j < 4; ++j) {
        const Point q = corrected(p, {shot.mean[2*j], shot.mean[2*j+1]});
        v[j] = {(q.x-p.mouseResX*.5f)/p.mouseResX, (q.y-p.mouseResY*.5f)/p.mouseResX, 1};
    }
    if(diamond) { lines[0] = cross(v[0],v[2]); lines[1] = cross(v[1],v[3]); }
    else {
        const Line center = cross(cross(v[0],v[3]), cross(v[1],v[2]));
        lines[0] = cross(center, cross(cross(v[0],v[2]), cross(v[1],v[3])));
        lines[1] = cross(center, cross(cross(v[0],v[1]), cross(v[2],v[3])));
    }
    for(unsigned axis = 0; axis < 2; ++axis) {
        const float norm = std::hypot(lines[axis].x, lines[axis].y);
        if(!(norm > 0)) { lines[axis] = {}; continue; }
        lines[axis] = {lines[axis].x/norm, lines[axis].y/norm, lines[axis].z/norm};
    }
}
bool onCentreLine(unsigned shot, unsigned axis) { return axis == 0 ? TargetX[shot] == .5f : TargetY[shot] == .5f; }

// Camera point of the barrel: least-squares point of the centre lines of all shots.
bool commonCenter(const CameraProfile& p, bool diamond, const Shots& shots, const float* weight, Point& source) {
    float a = 0, b = 0, c = 0, d = 0, e = 0;
    for(unsigned k = 0; k < 6; ++k) {
        Line lines[2];
        midlines(p, diamond, shots[k], lines);
        for(unsigned axis = 0; axis < 2; ++axis) if(onCentreLine(k, axis)) {
            const Line l = lines[axis]; const float w = weight[k];
            a += w*l.x*l.x; b += w*l.x*l.y; c += w*l.y*l.y; d += w*l.x*l.z; e += w*l.y*l.z;
        }
    }
    const float det = a*c-b*b;
    if(!(det > 0)) return false;
    const Point q = distorted(p, {p.mouseResX*(.5f+(b*e-c*d)/det), p.mouseResY*.5f+p.mouseResX*(b*d-a*e)/det});
    if(!std::isfinite(q.x) || !std::isfinite(q.y)) return false;
    source = q;
    return true;
}

// Perspective output of each shot, exactly as in Run mode.
bool aimPoints(const CameraProfile& p, bool diamond, const Shots& shots, const bool* used, Point source,
               float left, float right, Point z[6]) {
    for(unsigned k = 0; k < 6; ++k) {
        if(!used[k]) continue;
        OpenFIRE_Perspective per;   // a fresh instance reports (0,0) until a warp succeeds
        per.configure(p);
        per.source(source.x, source.y);
        const auto& q = shots[k].mean;
        if(diamond) per.warp(q[0],q[1],q[2],q[3],q[6],q[7],q[4],q[5],
                             res_x/2,0,0,res_y/2,res_x,res_y/2,res_x/2,res_y);
        else per.warp(q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],
                      left,0,right,0,left,res_y,right,res_y);
        z[k] = {float(per.getX()), float(per.getY())};
        if(z[k].x == 0 && z[k].y == 0) return false;
    }
    return true;
}

// Weighted straight line value = low + slope*target (target in screen units).
bool fitAxis(const float* value, const float* target, const float* w, float& low, float& slope) {
    float sw = 0, mt = 0, mv = 0;
    for(unsigned k = 0; k < 6; ++k) { sw += w[k]; mt += w[k]*target[k]; mv += w[k]*value[k]; }
    if(!(sw > 0)) return false;
    mt /= sw; mv /= sw;
    float stt = 0, stv = 0;
    for(unsigned k = 0; k < 6; ++k) { stt += w[k]*square(target[k]-mt); stv += w[k]*(target[k]-mt)*(value[k]-mv); }
    if(!(stt > 0)) return false;
    slope = stv/stt; low = mv-slope*mt;
    return slope > 0 && std::isfinite(slope) && std::isfinite(low);
}

} // namespace

bool Refine(const CameraProfile& p, bool diamond, const Shots& shots,
            const Parameters& original, Parameters& result) {
    result = original;
    if(p.mouseResX <= 0 || p.mouseResY <= 0 || !(p.lensCorrectionMin > 0) ||
       !(p.lensCorrectionMax >= p.lensCorrectionMin) || (!diamond && !(original.ledRight > original.ledLeft)))
        return false;
    // A shot taken while an LED was hidden or out of view is left out: the others still
    // fix both axes as long as each axis keeps two different targets.
    bool used[6];
    float tx[6], ty[6], weight[6], zx[6] = {}, zy[6] = {};
    for(unsigned k = 0; k < 6; ++k) {
        used[k] = shots[k].valid;
        for(unsigned a = 0; a < 8; ++a)
            if(!std::isfinite(shots[k].mean[a]) || shots[k].mean[a] < 0 || shots[k].mean[a] > (a&1 ? p.mouseMaxY : p.mouseMaxX)) used[k] = false;
        tx[k] = TargetX[k]*res_x; ty[k] = TargetY[k]*res_y; weight[k] = used[k];
    }
    Point source, z[6];
    float lowX = 0, slopeX = 0, lowY = 0, slopeY = 0;
    for(unsigned pass = 0; pass < ShotCount; ++pass) {
        if(!commonCenter(p, diamond, shots, weight, source)) return false;
        if(!aimPoints(p, diamond, shots, used, source, original.ledLeft, original.ledRight, z)) return false;
        for(unsigned k = 0; k < 6; ++k) { zx[k] = z[k].x; zy[k] = z[k].y; }
        if(!fitAxis(zx, tx, weight, lowX, slopeX) || !fitAxis(zy, ty, weight, lowY, slopeY)) return false;
        if(pass+1 == ShotCount) break;
        // Human aim errors are the same on every shot and in every direction: a shot far
        // from the line through the others is a mistake and counts less (Student-t weight,
        // tails set by the degrees of freedom left by the fit, spread from the median shot).
        float e[6], sorted[6];
        unsigned n = 0;
        for(unsigned k = 0; k < 6; ++k) if(used[k]) {
            e[k] = square((zx[k]-lowX)/slopeX-tx[k])+square((zy[k]-lowY)/slopeY-ty[k]);
            unsigned i = n++;
            for(; i > 0 && sorted[i-1] > e[k]; --i) sorted[i] = sorted[i-1];
            sorted[i] = e[k];
        }
        const float dof = 2.f*n-FitParameters;
        if(!(dof > 0)) break;
        const float variance = (sorted[(n-1)/2]+sorted[n/2])*.5f*(2.f*n/dof)/SquaredError2DMedian;
        if(!(variance > 0)) break;
        for(unsigned k = 0; k < 6; ++k) if(used[k]) weight[k] = 1/(dof+e[k]/variance);
    }
    Parameters n = original;
    n.sourceX = source.x; n.sourceY = source.y;
    const float offsets[4] = {lowY, res_y-lowY-slopeY*res_y, lowX, res_x-lowX-slopeX*res_x};
    for(unsigned k = 0; k < 4; ++k) {
        const float limit = k < 2 ? res_y : res_x;
        if(!std::isfinite(offsets[k]) || std::abs(offsets[k]) > limit) return false;
        n.offsets[k] = int32_t(std::lround(offsets[k]));
    }
    // MapCalibratedAxis keeps at least a quarter of each range for the screen.
    if(res_x-n.offsets[2]-n.offsets[3] < res_x/4 || res_y-n.offsets[0]-n.offsets[1] < res_y/4) return false;
    result = n;
    return true;
}

Session::Session(const CameraProfile& p, bool d):profile(p),diamond(d) {
    const uint32_t fps=p.fps?p.fps:209;
    // OLED/I/O work can briefly delay calibration reads. Keep those valid
    // samples, while timestamps still bound the total age of each window.
    maxGap=std::max<uint32_t>(50000,3000000u/fps);
    maxAge=std::max<uint32_t>(1,256000000u/fps)+maxGap;
}
void Session::invalidate() {count=0;next=0;}
void Session::reset() {invalidate();shots={};}
void Session::observe(const int32_t* xy, uint32_t now) {
    if(count && uint32_t(now-lastStamp)>maxGap) invalidate();
    for(unsigned a=0;a<8;++a) {
        if(xy[a]<0 || xy[a]>(a&1?profile.mouseMaxY:profile.mouseMaxX)) {invalidate();return;}
    }
    for(unsigned a=0;a<8;++a) history[next][a]=xy[a];
    stamps[next]=now;
    next=(next+1)&255;if(count<256) ++count;lastStamp=now;++frames;
}
bool Session::capture(uint8_t shot, uint32_t now) {
    if(shot>=shots.size()) return false;
    shots[shot]={};
    if(count<16 || uint32_t(now-lastStamp)>maxGap) {invalidate();return false;}
    constexpr unsigned sizes[5]={16,32,64,128,256};
    float mean[5][8]={},variance[8]={},chi[4]={},product[5]={1,0,0,0,0};
    unsigned windows=0;
    while(windows<5 && sizes[windows]<=count &&
          uint32_t(now-stamps[(next+256-sizes[windows])&255])<=maxAge) ++windows;
    if(!windows) {invalidate();return false;}
    for(unsigned a=0;a<8;++a) {
        // Integer sums stay exact in float (at most 256 values below 2^13).
        float sum=0; unsigned j=0;
        for(unsigned w=0;w<windows;++w) {
            for(;j<sizes[w];++j) sum+=history[(next+256-1-j)&255][a];
            mean[w][a]=sum/sizes[w];
        }
        for(j=0;j<sizes[0];++j) variance[a]+=square(history[(next+256-1-j)&255][a]-mean[0][a]);
        variance[a]/=sizes[0]-1;
        const float noise=variance[a]+square(profile.camToMouseMult)/12;
        for(unsigned w=1;w<windows;++w)
            chi[w-1]+=square(mean[w][a]-mean[w-1][a])/(noise*(1.0f/sizes[w-1]-1.0f/sizes[w]));
    }
    for(unsigned w=1;w<windows;++w) product[w]=product[w-1]/(1+square(chi[w-1]/24));
    float weights[5]={};
    for(unsigned w=0;w<windows;++w) weights[w]=product[w]-(w+1<windows?product[w+1]:0);
    for(unsigned a=0;a<8;++a)
        for(unsigned w=0;w<windows;++w) shots[shot].mean[a]+=weights[w]*mean[w][a];
    shots[shot].valid=true;
    invalidate(); // The next target starts its own pre-trigger history.
    return true;
}
bool Session::refine(const Parameters& original, Parameters& result) const {
    return Refine(profile,diamond,shots,original,result);
}

} // namespace OpenFIRECalibration
#endif // USE_PERSPECTIVE_ADVANCED
