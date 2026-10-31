/*  OpenFIRE Web App - lens refinement of the fitted camera ("Refine lens", Gun Settings > camera).

    The gun streams the four tracker vertices (sLensSample) while the user sweeps it over the
    whole screen. Seen by a real camera, every frame is a rigid view of the same flat LED
    frame: with the right lens values each view must be explained by a rotation and a position
    of the gun. The fit looks for the lens values (K1, K2 and the focal ratio) together with the
    unknown LED shape (aspect, skew, keystone) and the pose of every frame, minimising the
    reprojection error (bundle adjustment, Schur complement on the per-frame poses).

    Which values are refined depends on the data only, never on the camera model:
    record values -> K1 -> K1 + K2 -> K1 + K2 + focal ratio, each step taken only when it
    lowers the error significantly and its values are well determined; the curve must be
    physically valid (positive, radial order kept) over the whole field.

    The lens model is the firmware's (OpenFIRE_Perspective::applyLensCorrection):
        nx = (x - cx) / cx, ny = (y - cy) / cy, q = nx^2 + ny^2
        corrected = c + (p - c) * (1 + K1 q + K2 q^2)        (clamped to the record's limits)
    The firmware does not know the camera: it sends its coordinate space and the camera record
    values (sLensInfo), and only reads the refined values back from settings.

        const info = OF.LensFit.parseInfo(payload);              // sLensInfo
        const collector = new OF.LensFit.Collector(info);
        collector.add(OF.LensFit.parseSample(payload));          // sLensSample, every frame
        collector.progress();                                    // { total, screen, field, roll, ... }
        const result = OF.LensFit.solve(info, collector.frames); // { tier, k1, k2, focal, ... }
        OF.LensFit.toBits(value) / fromBits(bits)                // settings encoding (float bits)
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    // Firmware OF_Const::lensValueUnset: all ones (a NaN) = the camera record's value.
    const LENS_UNSET = 0xFFFFFFFF;

    // ----- Settings encoding --------------------------------------------------------------

    const f32 = new Float32Array(1);
    const u32 = new Uint32Array(f32.buffer);

    function toBits(value) {
        if (value === null || value === undefined || !Number.isFinite(value)) return LENS_UNSET;
        f32[0] = value;
        return u32[0] >>> 0;
    }

    function fromBits(bits) {
        bits = bits >>> 0;
        if (bits === LENS_UNSET) return null;
        u32[0] = bits;
        return Number.isFinite(f32[0]) ? f32[0] : null;
    }

    // Same check as the firmware (OpenFIRECamera::LensCurveValid): over q = 0..2 the factor
    // 1 + K1 q + K2 q^2 stays positive and the corrected radius grows with the radius,
    // 1 + 3 K1 q + 5 K2 q^2 > 0. Both are parabolas: their extremes are at the ends or at the
    // vertex -b/(2c) when inside, so the check is exact. Returns the extremes of the factor (new
    // limits). The values are checked as the firmware stores them (float32), with a small margin
    // so that the firmware never refuses a curve accepted here.
    const CURVE_MARGIN = 1e-3;

    function parabolaRange(b, c) {
        const end = 1 + 2 * b + 4 * c;
        let min = Math.min(1, end), max = Math.max(1, end);
        if (c !== 0) {
            const q = -b / (2 * c);
            if (q > 0 && q < 2) {
                const v = 1 + b * q + c * q * q;
                min = Math.min(min, v);
                max = Math.max(max, v);
            }
        }
        return { min, max };
    }

    function curveValid(k1, k2) {
        k1 = Math.fround(k1); k2 = Math.fround(k2);
        if (!(Math.abs(k1) <= 2) || !(Math.abs(k2) <= 2)) return { ok: false };
        const gain = parabolaRange(k1, k2), slope = parabolaRange(3 * k1, 5 * k2);
        if (!(gain.min > CURVE_MARGIN) || !(slope.min > CURVE_MARGIN)) return { ok: false };
        return { ok: true, min: gain.min, max: gain.max };
    }

    // ----- Serial payloads ----------------------------------------------------------------

    /** sLensInfo: version, mouseResX/Y (int32), record then active K1, K2, max, min, focal (float32). */
    function parseInfo(payload) {
        if (!payload || payload.length < 49 || payload[0] !== 1) return null;
        const view = new DataView(payload.buffer, payload.byteOffset, payload.byteLength);
        const group = (offset) => ({
            k1: view.getFloat32(offset, true), k2: view.getFloat32(offset + 4, true),
            max: view.getFloat32(offset + 8, true), min: view.getFloat32(offset + 12, true),
            focal: view.getFloat32(offset + 16, true),
        });
        const info = { resX: view.getInt32(1, true), resY: view.getInt32(5, true), record: group(9), active: group(29) };
        if (!(info.resX > 0) || !(info.resY > 0)) return null;
        return info;
    }

    /** sLensSample: version, flags (bit 0 Diamond, bit 1 complete), 4 x (x, y) int32. */
    function parseSample(payload) {
        if (!payload || payload.length < 34 || payload[0] !== 1) return null;
        const view = new DataView(payload.buffer, payload.byteOffset, payload.byteLength);
        const points = [];
        for (let i = 0; i < 4; ++i) points.push([view.getInt32(2 + i * 8, true), view.getInt32(6 + i * 8, true)]);
        return { diamond: !!(payload[1] & 1), complete: !!(payload[1] & 2), points };
    }

    // ----- Small linear algebra -------------------------------------------------------------

    function solveLinear(A, b) {
        // Gaussian elimination with partial pivoting; A is n x n (array of rows), returns x or null.
        const n = b.length;
        const M = A.map((row, i) => row.slice(0, n).concat([b[i]]));
        for (let c = 0; c < n; ++c) {
            let p = c;
            for (let r = c + 1; r < n; ++r) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
            if (!(Math.abs(M[p][c]) > 1e-300)) return null;
            if (p !== c) { const tmp = M[p]; M[p] = M[c]; M[c] = tmp; }
            for (let r = c + 1; r < n; ++r) {
                const f = M[r][c] / M[c][c];
                if (f === 0) continue;
                for (let k = c; k <= n; ++k) M[r][k] -= f * M[c][k];
            }
        }
        const x = new Array(n).fill(0);
        for (let r = n - 1; r >= 0; --r) {
            let s = M[r][n];
            for (let k = r + 1; k < n; ++k) s -= M[r][k] * x[k];
            x[r] = s / M[r][r];
        }
        return x.every(Number.isFinite) ? x : null;
    }

    function invert(A) {
        // Gauss-Jordan with partial pivoting, one pass for all columns.
        const n = A.length;
        const M = A.map((row) => row.slice(0, n));
        const I = A.map((_, r) => { const e = new Array(n).fill(0); e[r] = 1; return e; });
        for (let c = 0; c < n; ++c) {
            let p = c;
            for (let r = c + 1; r < n; ++r) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
            if (!(Math.abs(M[p][c]) > 1e-300)) return null;
            if (p !== c) { let t = M[p]; M[p] = M[c]; M[c] = t; t = I[p]; I[p] = I[c]; I[c] = t; }
            for (let r = c + 1; r < n; ++r) {
                const f = M[r][c] / M[c][c];
                if (f === 0) continue;
                for (let k = c; k < n; ++k) M[r][k] -= f * M[c][k];
                for (let k = 0; k < n; ++k) I[r][k] -= f * I[c][k];
            }
        }
        for (let r = n - 1; r >= 0; --r) {
            for (let k = r + 1; k < n; ++k) { const f = M[r][k]; if (f) for (let m = 0; m < n; ++m) I[r][m] -= f * I[k][m]; }
            const d = M[r][r];
            for (let m = 0; m < n; ++m) I[r][m] /= d;
        }
        return I.every((row) => row.every(Number.isFinite)) ? I : null;
    }

    function mul3(A, B) {
        const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        for (let i = 0; i < 3; ++i) for (let j = 0; j < 3; ++j)
            C[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
        return C;
    }

    function inv3(m) {
        const [a, b, c] = m[0], [d, e, f] = m[1], [g, h, i] = m[2];
        const A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g;
        const det = a * A + b * B + c * C;
        if (!(Math.abs(det) > 1e-300)) return null;
        const s = 1 / det;
        return [[A * s, -(b * i - c * h) * s, (b * f - c * e) * s],
            [B * s, (a * i - c * g) * s, -(a * f - c * d) * s],
            [C * s, -(a * h - b * g) * s, (a * e - b * d) * s]];
    }

    /** Unit square (0,0) (1,0) (1,1) (0,1) -> quad q[0..3] (perimeter order). Heckbert. */
    function squareToQuad(q) {
        const [x0, y0] = q[0], [x1, y1] = q[1], [x2, y2] = q[2], [x3, y3] = q[3];
        const sx = x0 - x1 + x2 - x3, sy = y0 - y1 + y2 - y3;
        const dx1 = x1 - x2, dx2 = x3 - x2, dy1 = y1 - y2, dy2 = y3 - y2;
        const den = dx1 * dy2 - dx2 * dy1;
        if (!(Math.abs(den) > 1e-300)) return null;
        const g = (sx * dy2 - dx2 * sy) / den, h = (dx1 * sy - sx * dy1) / den;
        return [[x1 - x0 + g * x1, x3 - x0 + h * x3, x0], [y1 - y0 + g * y1, y3 - y0 + h * y3, y0], [g, h, 1]];
    }

    function rotation(v) {
        // Rodrigues: rotation vector -> matrix.
        const th = Math.hypot(v[0], v[1], v[2]);
        if (th < 1e-12) return [[1, -v[2], v[1]], [v[2], 1, -v[0]], [-v[1], v[0], 1]];
        const k = [v[0] / th, v[1] / th, v[2] / th], c = Math.cos(th), s = Math.sin(th), t = 1 - c;
        return [[c + k[0] * k[0] * t, k[0] * k[1] * t - k[2] * s, k[0] * k[2] * t + k[1] * s],
            [k[1] * k[0] * t + k[2] * s, c + k[1] * k[1] * t, k[1] * k[2] * t - k[0] * s],
            [k[2] * k[0] * t - k[1] * s, k[2] * k[1] * t + k[0] * s, c + k[2] * k[2] * t]];
    }

    function rotationVector(R) {
        const c = Math.max(-1, Math.min(1, (R[0][0] + R[1][1] + R[2][2] - 1) / 2));
        const th = Math.acos(c);
        const w = [R[2][1] - R[1][2], R[0][2] - R[2][0], R[1][0] - R[0][1]];
        if (th < 1e-9) return [w[0] / 2, w[1] / 2, w[2] / 2];
        if (Math.PI - th < 1e-6) {
            // Near 180 degrees: axis from the diagonal.
            const ax = [Math.sqrt(Math.max(0, (R[0][0] + 1) / 2)), Math.sqrt(Math.max(0, (R[1][1] + 1) / 2)),
                Math.sqrt(Math.max(0, (R[2][2] + 1) / 2))];
            if (R[0][1] < 0) ax[1] = -ax[1];
            if (R[0][2] < 0) ax[2] = -ax[2];
            return ax.map((a) => a * th);
        }
        const s = th / (2 * Math.sin(th));
        return [w[0] * s, w[1] * s, w[2] * s];
    }

    function orthonormalise(R) {
        // Polar decomposition by averaging with the inverse transpose (converges in a few steps).
        let M = R.map((r) => r.slice());
        for (let it = 0; it < 30; ++it) {
            const inv = inv3(M);
            if (!inv) return null;
            const next = M.map((row, i) => row.map((v, j) => 0.5 * (v + inv[j][i])));
            let d = 0;
            for (let i = 0; i < 3; ++i) for (let j = 0; j < 3; ++j) d = Math.max(d, Math.abs(next[i][j] - M[i][j]));
            M = next;
            if (d < 1e-12) break;
        }
        return M;
    }

    // ----- Model --------------------------------------------------------------------------

    const SQUARE_PERIMETER = [0, 1, 3, 2];   // Square tracker order TL, TR, BL, BR
    const DIAMOND_PERIMETER = [0, 1, 2, 3];  // Diamond: top, right, bottom, left
    const UNIT = [[0, 0], [1, 0], [1, 1], [0, 1]];

    class Model {
        constructor(info) {
            this.cx = info.resX / 2;
            this.cy = info.resY / 2;
            this.f0 = info.resX / 2 * info.record.focal;   // nominal focal length, camera units
            this.record = info.record;
        }

        /** Lens-corrected point (unclamped when limits is null). */
        correct(p, k1, k2, limits) {
            const dx = p[0] - this.cx, dy = p[1] - this.cy;
            const nx = dx / this.cx, ny = dy / this.cy, q = nx * nx + ny * ny;
            let g = 1 + k1 * q + k2 * q * q;
            if (limits) g = Math.min(limits.max, Math.max(limits.min, g));
            return [this.cx + dx * g, this.cy + dy * g];
        }

        kinv(focalScale) {
            const f = this.f0 * focalScale;
            return [[1 / f, 0, -this.cx / f], [0, 1 / f, -this.cy / f], [0, 0, 1]];
        }
    }

    function shapeMatrix(s) {
        // LED frame: unit square -> plane, gauge fixed by vertex 0 at the origin and vertex 1 on +x.
        return [[1, s[1], 0], [0, s[0], 0], [s[2], s[3], 1]];
    }

    function shapePoints(s) {
        const A = shapeMatrix(s);
        return UNIT.map(([u, v]) => {
            const w = A[2][0] * u + A[2][1] * v + 1;
            return [(A[0][0] * u + A[0][1] * v) / w, (A[1][1] * v) / w];
        });
    }

    // Parameter vector of the globals: [k1, k2, focalScale, a, s, kx, ky] with a free mask.
    const G_K1 = 0, G_K2 = 1, G_F = 2;

    function expand(free, values, fixed) {
        const all = fixed.slice();
        let n = 0;
        for (let i = 0; i < all.length; ++i) if (free[i]) all[i] = values[n++];
        return all;
    }

    function geometricResiduals(model, frames, all) {
        const kinv = model.kinv(all[G_F]);
        const ainv = inv3(shapeMatrix(all.slice(3)));
        if (!ainv) return null;
        const out = new Float64Array(frames.length * 2);
        for (let i = 0; i < frames.length; ++i) {
            const q = frames[i].map((p) => model.correct(p, all[G_K1], all[G_K2], null));
            const G = squareToQuad(q);
            if (!G) return null;
            const M = mul3(mul3(kinv, G), ainv);
            const m1 = [M[0][0], M[1][0], M[2][0]], m2 = [M[0][1], M[1][1], M[2][1]];
            const n1 = Math.hypot(...m1), n2 = Math.hypot(...m2);
            out[2 * i] = (m1[0] * m2[0] + m1[1] * m2[1] + m1[2] * m2[2]) / (n1 * n2);
            out[2 * i + 1] = (n1 - n2) / (n1 + n2);
        }
        return out;
    }

    function sumSquares(r) { let s = 0; for (let i = 0; i < r.length; ++i) s += r[i] * r[i]; return s; }

    // The fits are generators: they yield after every iteration, so that solveAsync() can give
    // the page back to the browser (drawing, ESC) between two iterations; run() drives them in
    // one go (solve(), the tests).
    function run(steps) {
        let next = steps.next();
        while (!next.done) next = steps.next();
        return next.value;
    }

    /** Dense Levenberg-Marquardt with forward-difference Jacobian (few parameters). */
    function* levenbergSteps(fn, x0, steps, iterations) {
        let x = x0.slice();
        let r = fn(x);
        if (!r) return null;
        let cost = sumSquares(r);
        let lambda = 1e-3;
        const n = x.length;
        for (let it = 0; it < iterations; ++it) {
            yield;
            const J = [];
            for (let j = 0; j < n; ++j) {
                const xp = x.slice(); xp[j] += steps[j];
                const rp = fn(xp);
                if (!rp) return { x, cost };
                const col = new Float64Array(r.length);
                for (let k = 0; k < r.length; ++k) col[k] = (rp[k] - r[k]) / steps[j];
                J.push(col);
            }
            const H = [], g = [];
            for (let a = 0; a < n; ++a) {
                H.push(new Array(n).fill(0));
                let s = 0;
                for (let k = 0; k < r.length; ++k) s += J[a][k] * r[k];
                g.push(s);
                for (let b = 0; b <= a; ++b) {
                    let h = 0;
                    for (let k = 0; k < r.length; ++k) h += J[a][k] * J[b][k];
                    H[a][b] = h;
                }
            }
            for (let a = 0; a < n; ++a) for (let b = a + 1; b < n; ++b) H[a][b] = H[b][a];
            let improved = false;
            for (let tries = 0; tries < 10 && !improved; ++tries) {
                const D = H.map((row, a) => row.map((v, b) => a === b ? v * (1 + lambda) + 1e-30 : v));
                const dx = solveLinear(D, g.map((v) => -v));
                if (!dx) { lambda *= 10; continue; }
                const xn = x.map((v, j) => v + dx[j]);
                const rn = fn(xn);
                const cn = rn ? sumSquares(rn) : Infinity;
                if (cn < cost) {
                    const done = (cost - cn) < 1e-12 * cost;
                    x = xn; r = rn; cost = cn; lambda = Math.max(lambda / 10, 1e-9); improved = true;
                    if (done) return { x, cost };
                } else lambda *= 10;
            }
            if (!improved) break;
        }
        return { x, cost };
    }

    // ----- Initialisation -----------------------------------------------------------------

    /** LED shape from the most central frame, taken as seen face on. */
    function initialShape(model, frames, k1, k2) {
        let best = 0, bestD = Infinity;
        frames.forEach((f, i) => {
            const mx = (f[0][0] + f[1][0] + f[2][0] + f[3][0]) / 4 - model.cx;
            const my = (f[0][1] + f[1][1] + f[2][1] + f[3][1]) / 4 - model.cy;
            const d = mx * mx + my * my;
            if (d < bestD) { bestD = d; best = i; }
        });
        const kinv = model.kinv(1);
        const pts = frames[best].map((p) => {
            const c = model.correct(p, k1, k2, null);
            return [kinv[0][0] * c[0] + kinv[0][2], kinv[1][1] * c[1] + kinv[1][2]];
        });
        // Similarity: vertex 0 at the origin, vertex 1 on +x at distance 1.
        const ox = pts[0][0], oy = pts[0][1];
        const ang = Math.atan2(pts[1][1] - oy, pts[1][0] - ox);
        const len = Math.hypot(pts[1][0] - ox, pts[1][1] - oy);
        const c = Math.cos(-ang) / len, s = Math.sin(-ang) / len;
        const w = pts.map(([x, y]) => [(x - ox) * c - (y - oy) * s, (x - ox) * s + (y - oy) * c]);
        const A = squareToQuad(w);
        if (!A) return [1, 0, 0, 0];
        const sc = 1 / A[0][0];
        return [A[1][1] * sc, A[0][1] * sc, A[2][0], A[2][1]];
    }

    /** Pose of each frame from its homography (rotation vector, translation). */
    function initialPoses(model, frames, all) {
        const kinv = model.kinv(all[G_F]);
        const ainv = inv3(shapeMatrix(all.slice(3)));
        return frames.map((f) => {
            const q = f.map((p) => model.correct(p, all[G_K1], all[G_K2], null));
            const G = squareToQuad(q);
            let M = mul3(mul3(kinv, G), ainv);
            const lam = Math.hypot(M[0][0], M[1][0], M[2][0]);
            M = M.map((row) => row.map((v) => v / lam));
            if (M[2][2] < 0) M = M.map((row) => row.map((v) => -v));
            const r1 = [M[0][0], M[1][0], M[2][0]], r2 = [M[0][1], M[1][1], M[2][1]];
            const r3 = [r1[1] * r2[2] - r1[2] * r2[1], r1[2] * r2[0] - r1[0] * r2[2], r1[0] * r2[1] - r1[1] * r2[0]];
            const R = orthonormalise([[r1[0], r2[0], r3[0]], [r1[1], r2[1], r3[1]], [r1[2], r2[2], r3[2]]]) ||
                [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
            return rotationVector(R).concat([M[0][2], M[1][2], M[2][2]]);
        });
    }

    // ----- Bundle adjustment ----------------------------------------------------------------

    // Residual of each vertex: corrected point minus projection, brought back to raw camera
    // units through the inverse Jacobian of the correction (J = g I + d grad(g)^T), so that every
    // point weighs as its own measurement noise, at the corners as at the centre.
    function frameResiduals(model, frame, all, world, pose) {
        const R = rotation(pose);
        const f = model.f0 * all[G_F];
        const k1 = all[G_K1], k2 = all[G_K2];
        const out = new Float64Array(8);
        for (let j = 0; j < 4; ++j) {
            const [X, Y] = world[j];
            const pc0 = R[0][0] * X + R[0][1] * Y + pose[3];
            const pc1 = R[1][0] * X + R[1][1] * Y + pose[4];
            const pc2 = R[2][0] * X + R[2][1] * Y + pose[5];
            const c = model.correct(frame[j], k1, k2, null);
            const rx = c[0] - (model.cx + f * pc0 / pc2);
            const ry = c[1] - (model.cy + f * pc1 / pc2);
            const dx = frame[j][0] - model.cx, dy = frame[j][1] - model.cy;
            const q = (dx / model.cx) ** 2 + (dy / model.cy) ** 2;
            const g = 1 + k1 * q + k2 * q * q, dg = k1 + 2 * k2 * q;
            const gx = dg * 2 * dx / (model.cx * model.cx), gy = dg * 2 * dy / (model.cy * model.cy);
            const a = g + dx * gx, b = dx * gy, cc = dy * gx, d = g + dy * gy;
            const det = a * d - b * cc;
            if (det > 1e-6) {
                out[2 * j] = (d * rx - b * ry) / det;
                out[2 * j + 1] = (a * ry - cc * rx) / det;
            } else {
                out[2 * j] = rx;
                out[2 * j + 1] = ry;
            }
        }
        return out;
    }

    /**
     * Levenberg-Marquardt on the globals (free mask) and the 6 pose values of each frame, with
     * the Schur complement: the normal equations reduce to the few global parameters.
     * Returns { all, poses, cost, rms, covariance (free globals), frameRms }.
     */
    function bundle(model, frames, free, start, poses, iterations) {
        return run(bundleSteps(model, frames, free, start, poses, iterations));
    }

    function* bundleSteps(model, frames, free, start, poses, iterations) {
        const idx = free.map((on, i) => on ? i : -1).filter((i) => i >= 0);
        const P = idx.length, N = frames.length;
        let all = start.slice();
        let pose = poses.map((p) => p.slice());
        const gStep = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-6, 1e-6];
        const residuals = (a, ps) => {
            const world = shapePoints(a.slice(3));
            const out = [];
            for (let i = 0; i < N; ++i) out.push(frameResiduals(model, frames[i], a, world, ps[i]));
            return out;
        };
        const costOf = (rs) => rs.reduce((s, r) => s + sumSquares(r), 0);
        let rs = residuals(all, pose);
        let cost = costOf(rs);
        if (!Number.isFinite(cost)) return null;
        let lambda = 1e-3;
        let S = null;

        const linearise = function* () {
            // Per-frame Jacobians: globals (8 x P) and pose (8 x 6). Also yields every 64 frames,
            // so that one iteration over many frames does not hold the page for long.
            const world = shapePoints(all.slice(3));
            const Jg = [], Jp = [];
            const worldShift = idx.map((g) => {
                if (g < 3) return null;
                const a = all.slice(); a[g] += gStep[g];
                return { all: a, world: shapePoints(a.slice(3)) };
            });
            for (let i = 0; i < N; ++i) {
                if ((i & 63) === 63) yield;
                const r0 = rs[i];
                const jg = idx.map((g, n) => {
                    let rp;
                    if (g < 3) { const a = all.slice(); a[g] += gStep[g]; rp = frameResiduals(model, frames[i], a, world, pose[i]); }
                    else rp = frameResiduals(model, frames[i], worldShift[n].all, worldShift[n].world, pose[i]);
                    const col = new Float64Array(8);
                    for (let k = 0; k < 8; ++k) col[k] = (rp[k] - r0[k]) / gStep[g];
                    return col;
                });
                const jp = [];
                for (let k = 0; k < 6; ++k) {
                    const p = pose[i].slice();
                    const h = k < 3 ? 1e-6 : 1e-6 * Math.max(1, Math.abs(p[5]));
                    p[k] += h;
                    const rp = frameResiduals(model, frames[i], all, world, p);
                    const col = new Float64Array(8);
                    for (let m = 0; m < 8; ++m) col[m] = (rp[m] - r0[m]) / h;
                    jp.push(col);
                }
                Jg.push(jg); Jp.push(jp);
            }
            return { Jg, Jp };
        };

        const dot = (a, b) => { let s = 0; for (let k = 0; k < 8; ++k) s += a[k] * b[k]; return s; };

        for (let it = 0; it < iterations; ++it) {
            yield;
            const { Jg, Jp } = yield* linearise();
            // Blocks of the normal equations.
            const U = idx.map(() => new Array(P).fill(0)), bg = new Array(P).fill(0);
            const V = [], W = [], bp = [];
            for (let i = 0; i < N; ++i) {
                const r = rs[i];
                for (let a = 0; a < P; ++a) {
                    bg[a] += dot(Jg[i][a], r);
                    for (let b = 0; b <= a; ++b) U[a][b] += dot(Jg[i][a], Jg[i][b]);
                }
                const v = [], w = [], b6 = [];
                for (let a = 0; a < 6; ++a) {
                    v.push(new Array(6).fill(0));
                    b6.push(dot(Jp[i][a], r));
                    for (let b = 0; b <= a; ++b) v[a][b] = dot(Jp[i][a], Jp[i][b]);
                }
                for (let a = 0; a < 6; ++a) for (let b = a + 1; b < 6; ++b) v[a][b] = v[b][a];
                for (let a = 0; a < P; ++a) { w.push(new Array(6)); for (let b = 0; b < 6; ++b) w[a][b] = dot(Jg[i][a], Jp[i][b]); }
                V.push(v); W.push(w); bp.push(b6);
            }
            for (let a = 0; a < P; ++a) for (let b = a + 1; b < P; ++b) U[a][b] = U[b][a];

            let improved = false;
            for (let tries = 0; tries < 10 && !improved; ++tries) {
                if (tries) yield;
                const damp = (M) => M.map((row, a) => row.map((x, b) => a === b ? x * (1 + lambda) + 1e-20 : x));
                const Vinv = V.map((v) => invert(damp(v)));
                if (Vinv.some((v) => !v)) { lambda *= 10; continue; }
                const Sd = damp(U).map((row) => row.slice());
                const rhs = bg.map((v) => -v);
                for (let i = 0; i < N; ++i) {
                    // Y = W Vinv (P x 6)
                    const Y = W[i].map((row) => Vinv[i][0].map((_, c) => row.reduce((s, x, k) => s + x * Vinv[i][k][c], 0)));
                    for (let a = 0; a < P; ++a) {
                        for (let b = 0; b < P; ++b) Sd[a][b] -= Y[a].reduce((s, x, k) => s + x * W[i][b][k], 0);
                        rhs[a] += Y[a].reduce((s, x, k) => s + x * bp[i][k], 0);
                    }
                }
                const dg = P ? solveLinear(Sd, rhs) : [];
                if (!dg) { lambda *= 10; continue; }
                const allN = all.slice();
                idx.forEach((g, n) => { allN[g] += dg[n]; });
                const poseN = pose.map((p, i) => {
                    const t = bp[i].map((v, k) => -v - W[i].reduce((s, row, a) => s + row[k] * dg[a], 0));
                    const d = Vinv[i].map((row) => row.reduce((s, x, k) => s + x * t[k], 0));
                    return p.map((v, k) => v + d[k]);
                });
                const rsN = residuals(allN, poseN);
                const costN = costOf(rsN);
                if (costN < cost) {
                    const small = (cost - costN) < 1e-10 * cost;
                    all = allN; pose = poseN; rs = rsN; cost = costN;
                    lambda = Math.max(lambda / 10, 1e-12); improved = true;
                    if (small) it = iterations;
                } else lambda *= 10;
            }
            if (!improved) break;
        }

        // Marginal covariance of the free globals: sigma^2 * (U - W V^-1 W^T)^-1, undamped.
        const dof = N * 8 - (P + N * 6);
        const sigma2 = dof > 0 ? cost / dof : NaN;
        let covariance = null;
        if (P) {
            const { Jg, Jp } = yield* linearise();
            S = idx.map(() => new Array(P).fill(0));
            let singular = false;
            for (let i = 0; i < N; ++i) {
                const v = [];
                for (let a = 0; a < 6; ++a) { v.push(new Array(6)); for (let b = 0; b < 6; ++b) v[a][b] = dot(Jp[i][a], Jp[i][b]); }
                const vinv = invert(v);
                if (!vinv) singular = true;
                const w = idx.map((_, a) => Array.from({ length: 6 }, (__, b) => dot(Jg[i][a], Jp[i][b])));
                for (let a = 0; a < P; ++a) for (let b = 0; b < P; ++b) {
                    let s = dot(Jg[i][a], Jg[i][b]);
                    if (vinv) for (let k = 0; k < 6; ++k) for (let m = 0; m < 6; ++m) s -= w[a][k] * vinv[k][m] * w[b][m];
                    S[a][b] += s;
                }
            }
            // A pose that cannot be separated from the lens, or a covariance that is not a valid
            // one, gives no uncertainty at all: the tier is then refused (prudent fallback).
            const Sinv = singular ? null : invert(S);
            if (Sinv && Sinv.every((row, a) => row[a] > 0 && row.every(Number.isFinite)))
                covariance = Sinv.map((row) => row.map((x) => x * sigma2));
        }
        const frameRms = rs.map((r) => Math.sqrt(sumSquares(r) / 8));
        return { all, poses: pose, cost, rms: Math.sqrt(cost / (N * 8)), sigma2, covariance, idx, frameRms };
    }

    // ----- Solve ----------------------------------------------------------------------------

    const TIERS = [
        { name: 'record', free: [false, false, false] },
        { name: 'k1', free: [true, false, false] },
        { name: 'k1k2', free: [true, true, false] },
        { name: 'k1k2focal', free: [true, true, true] },
    ];

    // A tier is kept over the simpler one only with a clear gain (in units of the noise) and
    // values well determined over the part of the field the LEDs actually covered.
    const SIGNIFICANCE = 25;          // cost decrease / sigma^2 per added parameter
    const MAX_GAIN_STD = 0.004;       // 0.4% standard deviation of the correction factor
    const MAX_FOCAL_STD = 0.01;       // 1% standard deviation of the focal ratio

    function* fitTierSteps(model, frames, tier, shapeGuess, base) {
        // base: the lens values kept for what the tier does not refine (record by default).
        const b = base || { k1: model.record.k1, k2: model.record.k2, f: 1 };
        const free7 = tier.free.concat([true, true, true, true]);
        const fixed = [b.k1, b.k2, b.f].concat(shapeGuess);
        const x0 = fixed.filter((_, i) => free7[i]);
        const steps = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-6, 1e-6].filter((_, i) => free7[i]);
        const geo = yield* levenbergSteps((x) => geometricResiduals(model, frames, expand(free7, x, fixed)), x0, steps, 60);
        if (!geo) return null;
        const start = expand(free7, geo.x, fixed);
        const poses = initialPoses(model, frames, start);
        return yield* bundleSteps(model, frames, free7, start, poses, 40);
    }

    /** Cost of frames seen through fixed globals (lens and LED shape): only the poses are fitted. */
    function* heldOutCostSteps(model, frames, all) {
        const fit = yield* bundleSteps(model, frames, [false, false, false, false, false, false, false],
            all, initialPoses(model, frames, all), 20);
        return fit && Number.isFinite(fit.cost) ? fit.cost : Infinity;
    }

    function gainStd(fit, qObs) {
        // Largest standard deviation of 1 + K1 q + K2 q^2 over the covered radii (0..qObs):
        // K1 and K2 can compensate each other at one radius and not at another.
        if (!fit.covariance) return Infinity;
        const iK1 = fit.idx.indexOf(G_K1), iK2 = fit.idx.indexOf(G_K2);
        const n = fit.idx.length;
        let worst = 0;
        for (let s = 1; s <= 20; ++s) {
            const q = qObs * s / 20;
            const grad = new Array(n).fill(0);
            if (iK1 >= 0) grad[iK1] = q;
            if (iK2 >= 0) grad[iK2] = q * q;
            let v = 0;
            for (let a = 0; a < n; ++a) for (let b = 0; b < n; ++b) v += grad[a] * fit.covariance[a][b] * grad[b];
            if (!(v >= 0) || !Number.isFinite(v)) return Infinity;   // not a valid covariance
            worst = Math.max(worst, v);
        }
        return Math.sqrt(worst);
    }

    /** Lens values in use differ from the record (refined earlier). */
    function isCustom(info) {
        const a = info.active, r = info.record;
        if (!a) return false;
        const differs = (x, y) => Number.isFinite(x) && Number.isFinite(y) && Math.abs(x - y) > 1e-6 * Math.max(1, Math.abs(y));
        return differs(a.k1, r.k1) || differs(a.k2, r.k2) || differs(a.focal, r.focal);
    }

    /**
     * Fits the lens values. frames: arrays of 4 points in perimeter order, camera units.
     * Returns { ok, tier, k1, k2, focal (null = record value), rmsRecord, rms, rmsActive,
     * activeCustom, keepActive, validated, fits, frames, reason }.
     */
    function solve(info, framesIn) {
        return run(solveSteps(info, framesIn, {}));
    }

    /**
     * Same as solve(), giving the page back to the browser between iterations.
     * options.onProgress(fraction 0..1), options.isCancelled() -> true stops (result null).
     */
    async function solveAsync(info, framesIn, options = {}) {
        const ctx = { progress: 0 };
        const steps = solveSteps(info, framesIn, ctx);
        let slice = Date.now();
        for (;;) {
            const next = steps.next();
            if (next.done) return next.value;
            if (Date.now() - slice > 25) {
                if (options.onProgress) options.onProgress(ctx.progress);
                await new Promise((resolve) => setTimeout(resolve, 0));
                if (options.isCancelled && options.isCancelled()) return null;
                slice = Date.now();
            }
        }
    }

    function* solveSteps(info, framesIn, ctx = {}) {
        const model = new Model(info);
        const rec = info.record;
        if (!(rec.focal > 0)) return { ok: false, reason: 'focal' };
        if (!framesIn || framesIn.length < 40) return { ok: false, reason: 'frames' };
        let frames = framesIn;
        const custom = isCustom(info);
        const STAGES = 7 + (custom ? 1 : 0);
        let stage = 0;
        const done = () => { ctx.progress = Math.min(1, ++stage / STAGES); };

        // Outliers (a reflection taken for an LED): drop frames far from the first fit.
        const shape0 = initialShape(model, frames, rec.k1, rec.k2);
        const first = yield* fitTierSteps(model, frames, TIERS[2], shape0);
        done();
        if (first) {
            const sorted = first.frameRms.slice().sort((a, b) => a - b);
            const median = sorted[sorted.length >> 1];
            const keep = frames.filter((_, i) => first.frameRms[i] <= Math.max(4 * median, 1e-9));
            if (keep.length >= 40) frames = keep;
        }

        let qObs = 0;
        {
            const qs = [];
            for (const f of frames) for (const p of f) {
                const nx = (p[0] - model.cx) / model.cx, ny = (p[1] - model.cy) / model.cy;
                qs.push(nx * nx + ny * ny);
            }
            qs.sort((a, b) => a - b);
            qObs = qs[Math.floor(qs.length * 0.95)];
        }

        const shape = initialShape(model, frames, rec.k1, rec.k2);
        const fits = [];
        for (const tier of TIERS) { fits.push(yield* fitTierSteps(model, frames, tier, shape)); done(); }
        if (!fits[0]) return { ok: false, reason: 'fit' };

        // Noise level from the richest model that converged; every tier is compared with the
        // record values: gain in units of the noise, minus a price for each refined value.
        const noise = fits.slice().reverse().find((f) => f && Number.isFinite(f.sigma2) && f.sigma2 > 0);
        const sigma2 = noise ? noise.sigma2 : NaN;
        const candidates = [];
        const determined = [];   // valid and well determined, whether significant or not
        const checks = [];
        for (let t = 1; t < TIERS.length; ++t) {
            const fit = fits[t];
            if (!fit || !Number.isFinite(fit.cost) || !(sigma2 > 0)) { checks.push(null); continue; }
            const params = TIERS[t].free.filter(Boolean).length;
            const values = fit.all;
            const valid = curveValid(values[G_K1], values[G_K2]).ok;
            const gStd = gainStd(fit, qObs);
            const iF = fit.idx.indexOf(G_F);
            const fVar = iF >= 0 && fit.covariance ? fit.covariance[iF][iF] : 0;
            const fStd = fVar >= 0 && Number.isFinite(fVar) ? Math.sqrt(fVar) : Infinity;
            const focalOk = iF < 0 || (fStd < MAX_FOCAL_STD && values[G_F] > 0.5 && values[G_F] < 2);
            const score = (fits[0].cost - fit.cost) / sigma2 - SIGNIFICANCE * params;
            checks.push({ tier: TIERS[t].name, valid, gainStd: gStd, focalStd: fStd, score });
            if (valid && gStd < MAX_GAIN_STD && focalOk) determined.push({ t, score });
            if (valid && gStd < MAX_GAIN_STD && focalOk && score > 0) candidates.push({ t, score });
        }
        candidates.sort((a, b) => b.score - a.score);

        // Check on views not used for the fit: the frames are split in two by where the LEDs
        // fall on the camera (fine checkerboard, so both halves cover the whole field and the
        // nearly equal consecutive views stay together). A refined tier fitted on one half must
        // also see the other half better than the record values fitted on the same half.
        let chosen = 0, validated = null;
        if (candidates.length) {
            const half = [[], []];
            for (const f of frames) {
                const mx = (f[0][0] + f[1][0] + f[2][0] + f[3][0]) / 4, my = (f[0][1] + f[1][1] + f[2][1] + f[3][1]) / 4;
                const cell = Math.floor(mx / (2 * model.cx) * 8) + Math.floor(my / (2 * model.cy) * 8);
                half[cell & 1].push(f);
            }
            if (half[0].length >= 40 && half[1].length >= 40) {
                const [fitSet, testSet] = half[0].length >= half[1].length ? half : [half[1], half[0]];
                const shapeA = initialShape(model, fitSet, rec.k1, rec.k2);
                const recordA = yield* fitTierSteps(model, fitSet, TIERS[0], shapeA);
                const recordCost = recordA ? yield* heldOutCostSteps(model, testSet, recordA.all) : Infinity;
                done();
                validated = false;
                for (const c of candidates) {
                    const fitA = yield* fitTierSteps(model, fitSet, TIERS[c.t], shapeA);
                    const cost = fitA ? yield* heldOutCostSteps(model, testSet, fitA.all) : Infinity;
                    if (checks[c.t - 1]) checks[c.t - 1].heldOut = { record: recordCost, tier: cost };
                    if (cost < recordCost) { chosen = c.t; validated = true; break; }
                }
                done();
            } else chosen = candidates[0].t;   // too few views to split: statistical checks only
        }

        // Values already refined earlier: the new result must also beat them, otherwise they stay.
        const fit = fits[chosen];
        const tier = TIERS[chosen];
        let rmsActive = null, keepActive = false;
        if (custom) {
            const a = info.active;
            const fitActive = yield* fitTierSteps(model, frames, TIERS[0], shape,
                { k1: a.k1, k2: a.k2, f: a.focal > 0 ? a.focal / rec.focal : 1 });
            done();
            if (fitActive && Number.isFinite(fitActive.cost)) {
                rmsActive = fitActive.rms;
                // At least the gain asked of one parameter, also to go back to the record values:
                // views that cannot tell the two apart keep the values in use.
                const params = Math.max(1, tier.free.filter(Boolean).length);
                keepActive = !(sigma2 > 0) || (fitActive.cost - fit.cost) / sigma2 - SIGNIFICANCE * params <= 0;
            }
        }
        ctx.progress = 1;

        // The record values stay: the best tier that is valid and well determined but not
        // significant (or not confirmed on the other half) is still given, as values very close
        // to the record ones that the user may apply anyway.
        let near = null;
        if (chosen === 0 && determined.length) {
            const best = determined.reduce((a, b) => (b.score > a.score ? b : a));
            const nt = TIERS[best.t], nf = fits[best.t].all;
            near = { tier: nt.name, k1: nt.free[G_K1] ? nf[G_K1] : null, k2: nt.free[G_K2] ? nf[G_K2] : null,
                focal: nt.free[G_F] ? rec.focal * nf[G_F] : null };
        }

        return {
            ok: true,
            tier: tier.name,
            k1: tier.free[G_K1] ? fit.all[G_K1] : null,
            k2: tier.free[G_K2] ? fit.all[G_K2] : null,
            focal: tier.free[G_F] ? rec.focal * fit.all[G_F] : null,
            rmsRecord: fits[0].rms,
            rms: fit.rms,
            rmsActive,
            activeCustom: custom,
            keepActive,
            near,
            validated,
            qObs,
            frames: frames.length,
            dropped: framesIn.length - frames.length,
            checks,
            fits: fits.map((f, i) => f && { tier: TIERS[i].name, rms: f.rms, all: f.all, covariance: f.covariance }),
        };
    }

    // ----- Data collection ------------------------------------------------------------------

    const SCREEN_COLS = 6, SCREEN_ROWS = 4;      // cells to paint inside the screen
    const MIN_FRAMES = 250;
    // Frames per screen cell: a screen all green also holds the frames the fit needs, so that
    // "all green" never asks the user to go on without saying why.
    const SCREEN_TARGET = Math.ceil(MIN_FRAMES / (SCREEN_COLS * SCREEN_ROWS));
    // Camera field: the lens is radial, so what counts is how far from the centre the LEDs were
    // seen. Rings of the normalised radius r = sqrt(q) (1 = edge of the field at the sides).
    const FIELD_RINGS = [0, 0.3, 0.5, 0.65, 0.8];  // lower bound of each ring; the last one has no end
    const FIELD_TARGET = 12;                     // LED points per ring
    const ROLL_TARGET = 8 * Math.PI / 180;       // tilt reached on each side
    const ROLL_FRAMES = 15;                      // frames beyond the tilt, each side
    const AIM_INSET = 0.07;                      // guide targets at least 7% inside the screen edges (fully visible)
    const MAX_FRAMES = 900;
    const CELL_CAP = 12;                         // frames kept per screen cell and size band

    class Collector {
        constructor(info) {
            this.info = info;
            this.cx = info.resX / 2;
            this.cy = info.resY / 2;
            this.frames = [];
            this.screen = new Array(SCREEN_COLS * SCREEN_ROWS).fill(0);
            this.field = new Array(FIELD_RINGS.length).fill(0);
            this.rollLeft = 0;
            this.rollRight = 0;
            this.buckets = new Map();
            this.meta = [];           // what each kept frame added (parallel to frames)
            this.last = null;
            this.aim = null;          // last aim, as a screen fraction
            this.lastPoints = null;   // last vertices (camera units, perimeter order)
            this.lastComplete = false;
            this.live = null;         // last valid view: { pts, Hi, diamond, roll }, for the guide
            this.guideCell = -1;      // camera field cell the guide is pointing at (kept while useful)
        }

        /** Adds a sLensSample; returns true when the frame was kept for the fit. */
        add(sample) {
            if (!sample) return false;
            const order = sample.diamond ? DIAMOND_PERIMETER : SQUARE_PERIMETER;
            const pts = order.map((i) => sample.points[i]);
            this.lastPoints = pts;
            this.lastComplete = sample.complete;
            // An incomplete view (strict checks not met: the tracker's own vertices) still moves
            // the aim and the tilt on screen, so that they follow the gun smoothly; it is never
            // kept for the fit.
            const area = 0.5 * Math.abs((pts[2][0] - pts[0][0]) * (pts[3][1] - pts[1][1]) -
                (pts[3][0] - pts[1][0]) * (pts[2][1] - pts[0][1]));
            const fieldArea = this.info.resX * this.info.resY;
            if (!(area > fieldArea * 0.002)) return false;
            // Convex quad: the four turns have the same sign.
            let sign = 0;
            for (let i = 0; i < 4; ++i) {
                const a = pts[i], b = pts[(i + 1) & 3], c = pts[(i + 2) & 3];
                const z = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
                if (!z || (sign && Math.sign(z) !== sign)) return false;
                sign = Math.sign(z);
            }
            // Aim: the image centre in the LED frame, as a screen fraction (uncalibrated, only to paint).
            const H = squareToQuad(pts);
            const Hi = H && inv3(H);
            if (!Hi) return false;
            const w = Hi[2][0] * this.cx + Hi[2][1] * this.cy + Hi[2][2];
            const u = (Hi[0][0] * this.cx + Hi[0][1] * this.cy + Hi[0][2]) / w;
            const v = (Hi[1][0] * this.cx + Hi[1][1] * this.cy + Hi[1][2]) / w;
            const sx = sample.diamond ? 0.5 + 0.5 * (u - v) : u;
            const sy = sample.diamond ? 0.5 * (u + v) : v;
            this.aim = [sx, sy];
            {
                // Tilt of the gun (positive = clockwise seen from behind): the top edge of the LED
                // frame (Square), or left to right (Diamond). The camera image turns the other way,
                // unless the tracker mirrors it (the frame is then seen with the opposite winding).
                const a = sample.diamond ? pts[3] : pts[0], b = pts[1];
                const imageRoll = Math.atan2(b[1] - a[1], b[0] - a[0]);
                this.live = { pts, Hi, diamond: sample.diamond, roll: sign > 0 ? -imageRoll : imageRoll };
            }
            if (!sample.complete) return false;

            // A gun held still adds the same view again: wait for a move of 1% of the field.
            const scale = Math.hypot(this.info.resX, this.info.resY) * 0.01;
            if (this.last && pts.every((p, i) => Math.hypot(p[0] - this.last[i][0], p[1] - this.last[i][1]) < scale))
                return false;

            const size = Math.sqrt(area / fieldArea);
            const band = size < 0.25 ? 0 : size < 0.45 ? 1 : 2;
            const col = Math.floor(sx * SCREEN_COLS), row = Math.floor(sy * SCREEN_ROWS);
            const inside = col >= 0 && col < SCREEN_COLS && row >= 0 && row < SCREEN_ROWS;
            const key = (inside ? row * SCREEN_COLS + col : 'out' + Math.sign(Math.round(sx - 0.5)) + Math.sign(Math.round(sy - 0.5))) + ':' + band;
            const cells = [];
            for (const p of pts) {
                cells.push(this.ring(p));
            }
            const roll = this.live.roll;
            const meta = { key, cell: inside ? row * SCREEN_COLS + col : -1, cells,
                roll: roll > ROLL_TARGET ? 1 : roll < -ROLL_TARGET ? -1 : 0 };

            // Bounded collection that never gets stuck: when the screen cell (per size band) or
            // the total is full, a frame that adds something still missing replaces a redundant
            // one (whose removal leaves every count at its target), from the same cell when that
            // cell is full, otherwise from the most crowded one. A frame adding nothing is ignored.
            const bucketFull = (this.buckets.get(key) || 0) >= CELL_CAP;
            const totalFull = this.frames.length >= MAX_FRAMES;
            if (bucketFull || totalFull) {
                if (!this._adds(meta)) return false;
                let victim = this._redundant(bucketFull ? key : null);
                if (victim < 0 && bucketFull && totalFull) victim = this._redundant(null);
                if (victim >= 0) this._remove(victim);
                else if (totalFull) return false;
            }
            this.last = pts;
            this._keep(pts, meta);
            return true;
        }

        _keep(pts, meta) {
            this.frames.push(pts);
            this.meta.push(meta);
            this._count(meta, 1);
        }

        _count(meta, d) {
            this.buckets.set(meta.key, (this.buckets.get(meta.key) || 0) + d);
            if (meta.cell >= 0) this.screen[meta.cell] += d;
            for (const c of meta.cells) this.field[c] += d;
            if (meta.roll > 0) this.rollRight += d;
            if (meta.roll < 0) this.rollLeft += d;
        }

        /** true when the frame adds something that is still below its target. */
        _adds(meta) {
            return (meta.cell >= 0 && this.screen[meta.cell] < SCREEN_TARGET) ||
                meta.cells.some((c) => this.field[c] < FIELD_TARGET) ||
                (meta.roll > 0 && this.rollRight < ROLL_FRAMES) ||
                (meta.roll < 0 && this.rollLeft < ROLL_FRAMES);
        }

        /** Index of a kept frame whose removal leaves every count at its target, or -1. */
        _redundant(onlyKey) {
            let best = -1, bestCount = -1;
            for (let i = 0; i < this.meta.length; ++i) {
                const m = this.meta[i];
                if (onlyKey !== null && m.key !== onlyKey) continue;
                const count = this.buckets.get(m.key) || 0;
                if (count <= bestCount) continue;
                if (m.cell >= 0 && this.screen[m.cell] <= SCREEN_TARGET) continue;
                if (m.roll > 0 && this.rollRight <= ROLL_FRAMES) continue;
                if (m.roll < 0 && this.rollLeft <= ROLL_FRAMES) continue;
                if (m.cells.some((c) => this.field[c] - m.cells.filter((x) => x === c).length < FIELD_TARGET)) continue;
                best = i; bestCount = count;
            }
            return best;
        }

        _remove(i) {
            this._count(this.meta[i], -1);
            const last = this.frames.length - 1;
            this.frames[i] = this.frames[last]; this.meta[i] = this.meta[last];
            this.frames.pop(); this.meta.pop();
        }

        progress() {
            const screen = this.screen.reduce((s, n) => s + Math.min(1, n / SCREEN_TARGET), 0) / this.screen.length;
            const field = this.field.reduce((s, n) => s + Math.min(1, n / FIELD_TARGET), 0) / this.field.length;
            const rollLeft = Math.min(1, this.rollLeft / ROLL_FRAMES), rollRight = Math.min(1, this.rollRight / ROLL_FRAMES);
            const roll = (rollLeft + rollRight) / 2;
            const frames = Math.min(1, this.frames.length / MIN_FRAMES);
            const total = (screen + field + roll) / 3;
            return { total, screen, field, roll, rollLeft, rollRight, frames,
                done: screen >= 1 && field >= 1 && roll >= 1 && frames >= 1, count: this.frames.length };
        }

        /**
         * What the user has to do now, one step at a time:
         *   1 paint the screen (aim everywhere);
         *   2 bring the LEDs to the edges of the camera view: target = screen point to aim at,
         *     or closer = true when no target can be reached from this distance;
         *   3 tilt the gun left and right (roll = current tilt, radians; left/right = 0..1);
         *   4 done.
         */
        guide() {
            const p = this.progress();
            if (p.screen < 1 || p.frames < 1) return { step: 1, progress: p };
            if (p.field < 1) {
                // No reachable target: the outer rings need the gun closer, the inner ones farther.
                const target = this._fieldTarget();
                let missing = this.field.length - 1;
                while (missing >= 0 && this.field[missing] >= FIELD_TARGET) --missing;
                const outer = missing >= 2;
                return { step: 2, progress: p, target: target ? target.aim : null,
                    closer: !target && outer, farther: !target && !outer, cell: target ? target.cell : -1 };
            }
            if (p.roll < 1)
                return { step: 3, progress: p, roll: this.live ? this.live.roll : 0, left: p.rollLeft, right: p.rollRight,
                    limit: ROLL_TARGET };
            return { step: 4, progress: p };
        }

        /** Ring of the camera field of a point (camera units). */
        ring(p) {
            const r = Math.hypot((p[0] - this.cx) / this.cx, (p[1] - this.cy) / this.cy);
            let i = FIELD_RINGS.length - 1;
            while (i > 0 && r < FIELD_RINGS[i]) --i;
            return i;
        }

        /** Screen point to aim at so that one LED lands in a camera field ring still missing. */
        _fieldTarget() {
            const live = this.live;
            if (!live) return null;
            const { resX, resY } = this.info;
            const margin = 0.02 * Math.min(resX, resY);
            const toScreen = (x, y) => {
                const Hi = live.Hi;
                const w = Hi[2][0] * x + Hi[2][1] * y + Hi[2][2];
                const u = (Hi[0][0] * x + Hi[0][1] * y + Hi[0][2]) / w, v = (Hi[1][0] * x + Hi[1][1] * y + Hi[1][2]) / w;
                return live.diamond ? [0.5 + 0.5 * (u - v), 0.5 * (u + v)] : [u, v];
            };
            const aim = this.aim || [0.5, 0.5];
            // Moving the view by d moves every LED by d: the new aim is what now sits at centre - d.
            // Candidates: each LED brought to the middle of a missing ring, in 16 directions.
            const reach = (ring, dir) => {
                if (this.field[ring] >= FIELD_TARGET) return null;
                const r = ring + 1 < FIELD_RINGS.length ? (FIELD_RINGS[ring] + FIELD_RINGS[ring + 1]) / 2 : FIELD_RINGS[ring] + 0.06;
                const a = dir / 16 * 2 * Math.PI;
                const px = this.cx * (1 + r * Math.cos(a)), py = this.cy * (1 + r * Math.sin(a));
                let best = null;
                for (const led of live.pts) {
                    const dx = px - led[0], dy = py - led[1];
                    if (!live.pts.every((q) => q[0] + dx >= margin && q[0] + dx <= resX - margin &&
                        q[1] + dy >= margin && q[1] + dy <= resY - margin)) continue;
                    const s = toScreen(this.cx - dx, this.cy - dy);
                    if (!(s[0] >= AIM_INSET && s[0] <= 1 - AIM_INSET && s[1] >= AIM_INSET && s[1] <= 1 - AIM_INSET)) continue;
                    const dist = Math.hypot(s[0] - aim[0], s[1] - aim[1]);
                    if (!best || dist < best.dist) best = { aim: s, dist, cell: ring * 16 + dir };
                }
                return best;
            };
            // Keep the same target while it is still missing and reachable, so it does not jump around.
            if (this.guideCell >= 0) {
                const kept = reach(Math.floor(this.guideCell / 16), this.guideCell % 16);
                if (kept) return kept;
            }
            let best = null;
            for (let ring = this.field.length - 1; ring >= 0 && !best; --ring)   // the outermost missing ring first
                for (let dir = 0; dir < 16; ++dir) {
                    const c = reach(ring, dir);
                    if (c && (!best || c.dist < best.dist)) best = c;
                }
            this.guideCell = best ? best.cell : -1;
            return best;
        }
    }

    OF.LensFit = {
        LENS_UNSET, toBits, fromBits, curveValid, parseInfo, parseSample, solve, solveAsync, Collector,
        SCREEN_COLS, SCREEN_ROWS, SCREEN_TARGET, FIELD_RINGS, FIELD_TARGET, ROLL_TARGET,
        _internal: { Model, squareToQuad, rotation, rotationVector, bundle, initialShape, geometricResiduals },
    };
})(typeof globalThis !== 'undefined' ? globalThis : this);
