#include "OpenFIRE_Calibration.h"
#include "OpenFIREConst.h"
#include "OpenFIRE_Perspective_Advanced.h"
#include <algorithm>
#include <cmath>

#ifdef USE_PERSPECTIVE_ADVANCED
namespace OpenFIRECalibration {
namespace {

// Calibration targets as fractions of the screen: centre, top, bottom, left, right, centre,
// then the corners top-right, bottom-right, bottom-left, top-left.
constexpr float TargetX[ShotSlots] = {.5f, .5f, .5f, 0.f, 1.f, .5f, 1.f, 1.f, 0.f, 0.f};
constexpr float TargetY[ShotSlots] = {.5f, 0.f, 1.f, .5f, .5f, .5f, 0.f, 1.f, 1.f, 0.f};
using Shots = std::array<Measurement, ShotSlots>;
constexpr unsigned ShotCount = ShotSlots;
// Screen parameters fitted from the shots: offset and scale on each axis.
constexpr unsigned FitParameters = 4;
// With the LED shape (see Shape) the screen model is projective: four parameters more.
constexpr unsigned ShapeParameters = 8;
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
    for(unsigned k = 0; k < ShotCount; ++k) {
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
               float left, float right, Point z[ShotCount]) {
    for(unsigned k = 0; k < ShotCount; ++k) {
        if(!used[k]) continue;
        OpenFIRE_Perspective per;
        per.configure(p);
        per.source(source.x, source.y);
        const auto& q = shots[k].mean;
        // (0,0) is a valid output (the top-left corner): failure comes from warp() itself.
        const bool done = diamond ? per.warp(q[0],q[1],q[2],q[3],q[6],q[7],q[4],q[5],
                                             res_x/2,0,0,res_y/2,res_x,res_y/2,res_x/2,res_y)
                                  : per.warp(q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],
                                             left,0,right,0,left,res_y,right,res_y);
        if(!done) return false;
        z[k] = {float(per.getX()), float(per.getY())};
    }
    return true;
}

// Weighted straight line value = low + slope*target (target in screen units).
bool fitAxis(const float* value, const float* target, const float* w, float& low, float& slope) {
    float sw = 0, mt = 0, mv = 0;
    for(unsigned k = 0; k < ShotCount; ++k) { sw += w[k]; mt += w[k]*target[k]; mv += w[k]*value[k]; }
    if(!(sw > 0)) return false;
    mt /= sw; mv /= sw;
    float stt = 0, stv = 0;
    for(unsigned k = 0; k < ShotCount; ++k) { stt += w[k]*square(target[k]-mt); stv += w[k]*(target[k]-mt)*(value[k]-mv); }
    if(!(stt > 0)) return false;
    slope = stv/stt; low = mv-slope*mt;
    return slope > 0 && std::isfinite(slope) && std::isfinite(low);
}

// Edge offsets of MapCalibratedAxis for value = low + slope*target on each axis.
bool edgeOffsets(float lowX, float slopeX, float lowY, float slopeY, int32_t result[4]) {
    const float offsets[4] = {lowY, res_y-lowY-slopeY*res_y, lowX, res_x-lowX-slopeX*res_x};
    int32_t rounded[4];
    for(unsigned k = 0; k < 4; ++k) {
        const float limit = k < 2 ? res_y : res_x;
        if(!std::isfinite(offsets[k]) || std::abs(offsets[k]) > limit) return false;
        rounded[k] = int32_t(std::lround(offsets[k]));
    }
    // MapCalibratedAxis keeps at least a quarter of each range for the screen.
    if(res_x-rounded[2]-rounded[3] < res_x/4 || res_y-rounded[0]-rounded[1] < res_y/4) return false;
    std::copy(rounded, rounded+4, result);
    return true;
}

// Projective screen model h (h[8] = 1) from the target t to the Perspective output z, both
// centred and in screen units: z = (h0 tx+h1 ty+h2, h3 tx+h4 ty+h5) / (h6 tx+h7 ty+1).
// Weighted sum of squares of its linear residuals.
float linearResiduals(const float h[8], const float* tx, const float* ty, const float* u, const float* v,
                      const float* w) {
    float sum = 0;
    for(unsigned k = 0; k < ShotCount; ++k) if(w[k] > 0) {
        const float d = h[6]*tx[k]+h[7]*ty[k]+1;
        sum += w[k]*(square(h[0]*tx[k]+h[1]*ty[k]+h[2]-u[k]*d)+square(h[3]*tx[k]+h[4]*ty[k]+h[5]-v[k]*d));
    }
    return sum;
}

// Weighted least squares of the linear residuals: normal equations, Gauss-Jordan elimination.
bool fitProjective(const float* tx, const float* ty, const float* u, const float* v, const float* w,
                   float h[8]) {
    float a[ShapeParameters][ShapeParameters+1] = {};
    for(unsigned k = 0; k < ShotCount; ++k) if(w[k] > 0) {
        const float rows[2][ShapeParameters+1] = {
            {tx[k], ty[k], 1, 0, 0, 0, -u[k]*tx[k], -u[k]*ty[k], u[k]},
            {0, 0, 0, tx[k], ty[k], 1, -v[k]*tx[k], -v[k]*ty[k], v[k]}};
        for(const auto& row : rows)
            for(unsigned i = 0; i < ShapeParameters; ++i)
                for(unsigned j = 0; j <= ShapeParameters; ++j) a[i][j] += w[k]*row[i]*row[j];
    }
    for(unsigned c = 0; c < ShapeParameters; ++c) {
        unsigned pivot = c;
        for(unsigned r = c+1; r < ShapeParameters; ++r) if(std::abs(a[r][c]) > std::abs(a[pivot][c])) pivot = r;
        if(!(std::abs(a[pivot][c]) > 0)) return false;
        std::swap(a[c], a[pivot]);
        for(unsigned r = 0; r < ShapeParameters; ++r) if(r != c) {
            const float f = a[r][c]/a[c][c];
            for(unsigned j = c; j <= ShapeParameters; ++j) a[r][j] -= f*a[c][j];
        }
    }
    for(unsigned i = 0; i < ShapeParameters; ++i) {
        h[i] = a[i][ShapeParameters]/a[i][i];
        if(!std::isfinite(h[i])) return false;
    }
    return true;
}

// LED shape from all the shots, with the weights of the per-axis fit. The four parameters that
// the projective model adds are shrunk toward the per-axis model as much as the shots fail to
// confirm them (positive-part James-Stein estimator on the linear residuals: no threshold, only
// the numbers of shots and parameters). The inverse model, from the Perspective output to the
// target, is split into the shape (applied first) and the per-axis part (the edge offsets):
// [a 0 b; 0 c d; 0 0 1] * [1 skewX 0; skewY 1 0; keystoneX keystoneY 1].
// On success low and slope describe the axes after the shape correction.
bool fitShape(const float* zx, const float* zy, const float* weight,
              float& lowX, float& slopeX, float& lowY, float& slopeY, Shape& shape) {
    float tx[ShotCount], ty[ShotCount], u[ShotCount], v[ShotCount];
    unsigned n = 0;
    for(unsigned k = 0; k < ShotCount; ++k) {
        tx[k] = TargetX[k]-.5f; ty[k] = TargetY[k]-.5f;
        u[k] = zx[k]/res_x-.5f; v[k] = zy[k]/res_y-.5f;
        n += weight[k] > 0;
    }
    const float axis[ShapeParameters] = {slopeX, 0, lowX/res_x+(slopeX-1)*.5f,
                                         0, slopeY, lowY/res_y+(slopeY-1)*.5f, 0, 0};
    float h[ShapeParameters];
    if(!fitProjective(tx, ty, u, v, weight, h)) return false;
    const float full = linearResiduals(h, tx, ty, u, v, weight);
    const float dof = 2.f*n-ShapeParameters, gain = linearResiduals(axis, tx, ty, u, v, weight)-full;
    if(!(dof > 0) || !(gain > 0)) return false;
    const float kept = 1-(ShapeParameters-FitParameters-2.f)*full/((dof+2)*gain);
    if(!(kept > 0)) return false;
    for(unsigned i = 0; i < ShapeParameters; ++i) h[i] = axis[i]+kept*(h[i]-axis[i]);
    // Inverse model from the adjugate, scaled to a unit last element.
    const float scale = h[0]*h[4]-h[1]*h[3];
    if(!(scale > 0)) return false;
    const float b = (h[1]*h[5]-h[2]*h[4])/scale, d = (h[2]*h[3]-h[0]*h[5])/scale;
    const float keystoneX = (h[3]*h[7]-h[4]*h[6])/scale, keystoneY = (h[1]*h[6]-h[0]*h[7])/scale;
    const float a = (h[4]-h[5]*h[7])/scale-b*keystoneX, c = (h[0]-h[2]*h[6])/scale-d*keystoneY;
    if(!(a > 0) || !(c > 0)) return false;
    // The correction must hold over the whole screen: positive denominator at every shot.
    for(unsigned k = 0; k < ShotCount; ++k)
        if(weight[k] > 0 && !(1+keystoneX*u[k]+keystoneY*v[k] > 0)) return false;
    shape.skewX = ((h[2]*h[7]-h[1])/scale-b*keystoneY)/a;
    shape.skewY = ((h[5]*h[6]-h[3])/scale-d*keystoneX)/c;
    shape.keystoneX = keystoneX; shape.keystoneY = keystoneY;
    slopeX = 1/a; lowX = res_x*(.5f-(.5f+b)/a);
    slopeY = 1/c; lowY = res_y*(.5f-(.5f+d)/c);
    return true;
}

// Mean of the four LEDs on the screen, through the screen model fitted with this source (per
// axis, or with the LED shape): where the shots place the LEDs.
bool ledMean(const CameraProfile& p, bool diamond, const Shots& shots, const bool* used, const float* tx,
             const float* ty, const float* weight, bool withShape, float left, float right, Point source,
             Point& mean, float fit[4], Shape& shape) {
    Point z[ShotCount];
    if(!aimPoints(p, diamond, shots, used, source, left, right, z)) return false;
    float zx[ShotCount], zy[ShotCount];
    for(unsigned k = 0; k < ShotCount; ++k) { zx[k] = z[k].x; zy[k] = z[k].y; }
    if(!fitAxis(zx, tx, weight, fit[0], fit[1]) || !fitAxis(zy, ty, weight, fit[2], fit[3])) return false;
    shape = {};
    if(withShape && !fitShape(zx, zy, weight, fit[0], fit[1], fit[2], fit[3], shape)) return false;
    // The LEDs in Perspective output, as in aimPoints.
    const Point led[4] = {{diamond ? res_x*.5f : left, 0}, {diamond ? 0 : right, diamond ? res_y*.5f : 0},
                          {diamond ? float(res_x) : left, diamond ? res_y*.5f : res_y},
                          {diamond ? res_x*.5f : right, float(res_y)}};
    mean = {};
    for(const Point& q : led) {
        // ApplyShape without rounding, then the edge offsets.
        const float u = q.x/res_x-.5f, v = q.y/res_y-.5f, w = 1+shape.keystoneX*u+shape.keystoneY*v;
        if(!(w > 0)) return false;
        mean.x += ((((u+shape.skewX*v)/w+.5f)*res_x)-fit[0])/fit[1]*.25f;
        mean.y += ((((shape.skewY*u+v)/w+.5f)*res_y)-fit[2])/fit[3]*.25f;
    }
    return std::isfinite(mean.x) && std::isfinite(mean.y);
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
    bool used[ShotCount];
    float tx[ShotCount], ty[ShotCount], weight[ShotCount], zx[ShotCount] = {}, zy[ShotCount] = {};
    unsigned taken = 0;
    for(unsigned k = 0; k < ShotCount; ++k) {
        used[k] = shots[k].valid;
        for(unsigned a = 0; a < 8; ++a)
            if(!std::isfinite(shots[k].mean[a]) || shots[k].mean[a] < 0 || shots[k].mean[a] > (a&1 ? p.mouseMaxY : p.mouseMaxX)) used[k] = false;
        tx[k] = TargetX[k]*res_x; ty[k] = TargetY[k]*res_y; weight[k] = used[k];
        taken += used[k];
    }
    Point source, z[ShotCount];
    float lowX = 0, slopeX = 0, lowY = 0, slopeY = 0;
    for(unsigned pass = 0; ; ++pass) {
        if(!commonCenter(p, diamond, shots, weight, source)) return false;
        if(!aimPoints(p, diamond, shots, used, source, original.ledLeft, original.ledRight, z)) return false;
        for(unsigned k = 0; k < ShotCount; ++k) { zx[k] = z[k].x; zy[k] = z[k].y; }
        if(!fitAxis(zx, tx, weight, lowX, slopeX) || !fitAxis(zy, ty, weight, lowY, slopeY)) return false;
        if(pass+1 >= taken) break;   // one pass per shot
        // Human aim errors are the same on every shot and in every direction: a shot far
        // from the line through the others is a mistake and counts less (Student-t weight,
        // tails set by the degrees of freedom left by the fit, spread from the median shot).
        float e[ShotCount], sorted[ShotCount];
        unsigned n = 0;
        for(unsigned k = 0; k < ShotCount; ++k) if(used[k]) {
            e[k] = square((zx[k]-lowX)/slopeX-tx[k])+square((zy[k]-lowY)/slopeY-ty[k]);
            unsigned i = n++;
            for(; i > 0 && sorted[i-1] > e[k]; --i) sorted[i] = sorted[i-1];
            sorted[i] = e[k];
        }
        const float dof = 2.f*n-FitParameters;
        if(!(dof > 0)) break;
        const float variance = (sorted[(n-1)/2]+sorted[n/2])*.5f*(2.f*n/dof)/SquaredError2DMedian;
        if(!(variance > 0)) break;
        for(unsigned k = 0; k < ShotCount; ++k) if(used[k]) weight[k] = 1/(dof+e[k]/variance);
    }
    Parameters n = original;
    n.sourceX = source.x; n.sourceY = source.y;
    n.shape = {};
    if(!edgeOffsets(lowX, slopeX, lowY, slopeY, n.offsets)) return false;
    // With every shot, the corners included, the LED shape is measured too. It is kept only
    // with valid edge offsets; otherwise the per-axis result above stays.
    Shape shape;
    const bool withShape = taken == ShotCount && fitShape(zx, zy, weight, lowX, slopeX, lowY, slopeY, shape) &&
                           edgeOffsets(lowX, slopeX, lowY, slopeY, n.offsets);
    if(withShape) n.shape = shape;
    // From one distance a shift of the source and a shift of the screen against the LEDs look the
    // same, and the centre lines above split them assuming a rectangle of LEDs: with uneven LEDs
    // (e.g. bars of different length) the split is wrong and the error grows away from the
    // calibration distance and with the gun tilted. The LEDs centred on the screen split them with
    // all the shots: the source moves until the fitted model puts the mean of the four LEDs at the
    // screen centre. The mean is almost linear in the source: Newton steps, derivatives over one
    // camera pixel; the first step lands within a fraction of a pixel, the second removes the rest.
    // On any failure the result above stays.
    float fit[4];
    Shape fitted;
    auto meanAt = [&](Point at, Point& m) {
        return ledMean(p, diamond, shots, used, tx, ty, weight, withShape, original.ledLeft, original.ledRight,
                       at, m, fit, fitted);
    };
    Point s = source, mean, movedX, movedY;
    const float h = float(p.camToMouseMult);
    bool centred = meanAt(s, mean);
    for(unsigned step = 0; centred && step < 2; ++step) {
        centred = meanAt({s.x+h, s.y}, movedX) && meanAt({s.x, s.y+h}, movedY);
        const float a = (movedX.x-mean.x)/h, b = (movedY.x-mean.x)/h, c = (movedX.y-mean.y)/h, d = (movedY.y-mean.y)/h;
        const float det = a*d-b*c, ex = mean.x-res_x*.5f, ey = mean.y-res_y*.5f;
        centred = centred && std::abs(det) > 0;
        if(centred) s = {s.x-(d*ex-b*ey)/det, s.y-(a*ey-c*ex)/det};
        // Last: the model fitted at the final source.
        centred = centred && meanAt(s, mean);
    }
    int32_t offsets[4];
    if(centred && std::isfinite(s.x) && std::isfinite(s.y) && edgeOffsets(fit[0], fit[1], fit[2], fit[3], offsets)) {
        n.sourceX = s.x; n.sourceY = s.y; n.shape = fitted;
        std::copy(offsets, offsets+4, n.offsets);
    }
    result = n;
    return true;
}

void ApplyShape(const Shape& shape, int32_t limit, int32_t& x, int32_t& y) {
    if(shape.skewX == 0 && shape.skewY == 0 && shape.keystoneX == 0 && shape.keystoneY == 0) return;
    const float u = float(x)/res_x-.5f, v = float(y)/res_y-.5f;
    const float w = 1+shape.keystoneX*u+shape.keystoneY*v;
    if(!(w > 0)) return;   // beyond the horizon of the correction, far off the screen
    const float range = float(limit);
    x = int32_t(std::lround(std::clamp(((u+shape.skewX*v)/w+.5f)*res_x, -range, range)));
    y = int32_t(std::lround(std::clamp(((shape.skewY*u+v)/w+.5f)*res_y, -range, range)));
}

void ApplyParallax(float offset, float calibrationArea, float area, float angle, int32_t span, int32_t limit,
                   int32_t& x, int32_t& y) {
    if(offset == 0 || !(calibrationArea > 0) || !(area > 0)) return;
    // The camera sits below the sight in the gun: with the gun tilted the shift turns with it.
    const float range = float(limit), shift = offset*span*(std::sqrt(calibrationArea/area)-1);
    x = int32_t(std::lround(std::clamp(x-shift*std::sin(angle), -range, range)));
    y = int32_t(std::lround(std::clamp(y+shift*std::cos(angle), -range, range)));
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
