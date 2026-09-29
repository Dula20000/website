/* SRM internal ballistics — JavaScript port of the Python solver in code/srm/.
   Level-set (fast-marching) grain regression + lumped chamber ballistics.
   Mirrors geometry.py, propellant.py and ballistics.py line for line so that the
   browser tool reproduces the Python numbers (checked by code/tests/js_crosscheck.js).
   Works as a browser global (window.SRM) and as a CommonJS module (node). */
(function (root) {
  'use strict';
  const PI = Math.PI;
  const P_REF = 1.0e6;
  const G0 = 9.80665;

  // ------------------------------------------------------------------ geometry
  function linspace(a, b, n) {
    const out = new Float64Array(n);
    const step = (b - a) / (n - 1);
    for (let i = 0; i < n; i++) out[i] = a + i * step;
    if (n > 1) out[n - 1] = b;
    return out;
  }

  function makeGrid(R, N) {
    const h = 2.0 * R / (N - 5);
    const half = 0.5 * (N - 1) * h;
    const x = new Float64Array(N);
    for (let i = 0; i < N; i++) x[i] = -half + h * i;
    return { x, h, N };
  }

  function starVertices(points, rTip, rValley) {
    const v = [];
    for (let k = 0; k < points; k++) {
      const a = 2 * PI * k / points;
      v.push([rTip * Math.cos(a), rTip * Math.sin(a)]);
      const b = a + PI / points;
      v.push([rValley * Math.cos(b), rValley * Math.sin(b)]);
    }
    return v;
  }

  function sdPolygon(X, Y, verts) {
    let d2 = Infinity, inside = false;
    const n = verts.length;
    for (let i = 0; i < n; i++) {
      const ax = verts[i][0], ay = verts[i][1];
      const bx = verts[(i + 1) % n][0], by = verts[(i + 1) % n][1];
      const ex = bx - ax, ey = by - ay;
      const wx = X - ax, wy = Y - ay;
      let t = (wx * ex + wy * ey) / (ex * ex + ey * ey);
      t = Math.min(Math.max(t, 0), 1);
      const dx = wx - t * ex, dy = wy - t * ey;
      d2 = Math.min(d2, dx * dx + dy * dy);
      const cond = (ay > Y) !== (by > Y);
      const xint = ax + (Y - ay) * ex / (ey !== 0 ? ey : 1e-300);
      if (cond && X < xint) inside = !inside;
    }
    const d = Math.sqrt(d2);
    return inside ? -d : d;
  }

  function sdSlot(X, Y, theta, r0, r1, width) {
    const c = Math.cos(theta), s = Math.sin(theta);
    const u = X * c + Y * s - 0.5 * (r0 + r1);
    const v = -X * s + Y * c;
    const qx = Math.abs(u) - 0.5 * (r1 - r0);
    const qy = Math.abs(v) - 0.5 * width;
    const outside = Math.hypot(Math.max(qx, 0), Math.max(qy, 0));
    const inside = Math.min(Math.max(qx, qy), 0);
    return outside + inside;
  }

  /** Implicit port function phi0 (negative in the port) on the grid, row-major [iy*N+ix]. */
  function portSdf(geom, grid) {
    const { x, N } = grid;
    const phi = new Float64Array(N * N);
    let verts = null;
    if (geom.type === 'star') verts = starVertices(geom.points, geom.r_tip, geom.r_valley);
    for (let i = 0; i < N; i++) {
      const Y = x[i];
      for (let j = 0; j < N; j++) {
        const X = x[j];
        let v;
        if (geom.type === 'bates') v = Math.hypot(X, Y) - 0.5 * geom.core_d;
        else if (geom.type === 'moon') v = Math.hypot(X - geom.offset, Y) - 0.5 * geom.core_d;
        else if (geom.type === 'star') v = sdPolygon(X, Y, verts);
        else if (geom.type === 'finocyl') {
          v = Math.hypot(X, Y) - 0.5 * geom.core_d;
          for (let k = 0; k < geom.slots; k++) {
            const th = 2 * PI * k / geom.slots;
            v = Math.min(v, sdSlot(X, Y, th, 0.0, geom.slot_r, geom.slot_w));
          }
        } else throw new Error('unknown geometry ' + geom.type);
        phi[i * N + j] = v;
      }
    }
    return phi;
  }

  // binary min-heap on (T, k), same ordering as Python heapq with tuples
  function Heap() { this.T = []; this.K = []; }
  Heap.prototype.less = function (i, j) {
    const a = this.T[i], b = this.T[j];
    return a < b || (a === b && this.K[i] < this.K[j]);
  };
  Heap.prototype.swap = function (i, j) {
    let t = this.T[i]; this.T[i] = this.T[j]; this.T[j] = t;
    t = this.K[i]; this.K[i] = this.K[j]; this.K[j] = t;
  };
  Heap.prototype.push = function (t, k) {
    this.T.push(t); this.K.push(k);
    let i = this.T.length - 1;
    while (i > 0) {
      const p = (i - 1) >> 1;
      if (this.less(i, p)) { this.swap(i, p); i = p; } else break;
    }
  };
  Heap.prototype.pop = function () {
    const t = this.T[0], k = this.K[0];
    const lt = this.T.pop(), lk = this.K.pop();
    if (this.T.length) {
      this.T[0] = lt; this.K[0] = lk;
      let i = 0; const n = this.T.length;
      for (;;) {
        const l = 2 * i + 1, r = l + 1;
        let m = i;
        if (l < n && this.less(l, m)) m = l;
        if (r < n && this.less(r, m)) m = r;
        if (m === i) break;
        this.swap(i, m); i = m;
      }
    }
    return [t, k];
  };

  /** Second-order fast marching from the port surface (see geometry.fmm). */
  function fmm(phi0, N, h, bandCells) {
    if (bandCells === undefined) bandCells = 3.0;
    const n = N * N;
    const T = new Float64Array(n).fill(Infinity);
    const state = new Int8Array(n);
    const neg = new Uint8Array(n), band = new Uint8Array(n);
    for (let k = 0; k < n; k++) neg[k] = phi0[k] <= 0 ? 1 : 0;
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
      const k = i * N + j;
      if (i + 1 < N && neg[k] !== neg[k + N]) { band[k] = 1; band[k + N] = 1; }
      if (j + 1 < N && neg[k] !== neg[k + 1]) { band[k] = 1; band[k + 1] = 1; }
    }
    if (bandCells > 0) for (let k = 0; k < n; k++) if (phi0[k] > 0 && phi0[k] < bandCells * h) band[k] = 1;
    for (let k = 0; k < n; k++) {
      if (band[k]) { T[k] = phi0[k]; state[k] = 2; } else if (neg[k]) { T[k] = phi0[k]; state[k] = 3; }
    }
    const heap = new Heap();
    const hh = h * h;
    function axisTerm(i, j, di, dj) {
      let bc = 0, bm = 0, have = false;
      for (let s = -1; s <= 1; s += 2) {
        const i1 = i + s * di, j1 = j + s * dj;
        if (i1 < 0 || i1 >= N || j1 < 0 || j1 >= N) continue;
        const k1 = i1 * N + j1;
        if (state[k1] !== 2) continue;
        const t1 = T[k1];
        let c = 1.0, m = t1;
        const i2 = i1 + s * di, j2 = j1 + s * dj;
        if (i2 >= 0 && i2 < N && j2 >= 0 && j2 < N) {
          const k2 = i2 * N + j2;
          if (state[k2] === 2 && T[k2] <= t1) { c = 2.25; m = (4.0 * t1 - T[k2]) / 3.0; }
        }
        if (!have || m < bm) { bc = c; bm = m; have = true; }
      }
      return have ? [bc, bm] : null;
    }
    function solve(k) {
      const i = Math.floor(k / N), j = k - i * N;
      const ax = axisTerm(i, j, 0, 1), ay = axisTerm(i, j, 1, 0);
      const terms = [];
      if (ax) terms.push(ax);
      if (ay) terms.push(ay);
      if (terms.length === 2 && terms[1][1] < terms[0][1]) terms.reverse();
      const c1 = terms[0][0], m1 = terms[0][1];
      let tnew = m1 + h / Math.sqrt(c1);
      if (terms.length === 2 && tnew > terms[1][1]) {
        const c2 = terms[1][0], m2 = terms[1][1];
        const a = c1 + c2;
        const b = -2.0 * (c1 * m1 + c2 * m2);
        const c = c1 * m1 * m1 + c2 * m2 * m2 - hh;
        const disc = b * b - 4.0 * a * c;
        if (disc >= 0.0) tnew = (-b + Math.sqrt(disc)) / (2.0 * a);
      }
      return tnew;
    }
    const DI = [0, 0, 1, -1], DJ = [1, -1, 0, 0];
    function pushNeighbours(k) {
      const i = Math.floor(k / N), j = k - i * N;
      for (let q = 0; q < 4; q++) {
        const i1 = i + DI[q], j1 = j + DJ[q];
        if (i1 < 0 || i1 >= N || j1 < 0 || j1 >= N) continue;
        const k1 = i1 * N + j1;
        if ((state[k1] === 0 || state[k1] === 1) && !neg[k1]) {
          const tn = solve(k1);
          if (tn < T[k1]) { T[k1] = tn; state[k1] = 1; heap.push(tn, k1); }
        }
      }
    }
    for (let k = 0; k < n; k++) if (band[k] && !neg[k]) pushNeighbours(k);
    while (heap.T.length) {
      const [tk, k] = heap.pop();
      if (state[k] === 2 || tk > T[k]) continue;
      state[k] = 2;
      pushNeighbours(k);
    }
    return T;
  }

  /** Marching squares on field F at `level`; calls cb(x1,y1,x2,y2) with segments oriented so F<level is on the left. */
  function contour(F, grid, level, cb) {
    const { x, N, h } = grid;
    const ex = [0, 0, 0, 0], ey = [0, 0, 0, 0], has = [false, false, false, false];
    const tt = (va, vb) => va / (va - vb);
    for (let i = 0; i < N - 1; i++) {
      const Y0 = x[i];
      for (let j = 0; j < N - 1; j++) {
        const k = i * N + j;
        const v00 = F[k] - level, v10 = F[k + 1] - level, v11 = F[k + N + 1] - level, v01 = F[k + N] - level;
        const s00 = v00 < 0, s10 = v10 < 0, s11 = v11 < 0, s01 = v01 < 0;
        has[0] = s00 !== s10; has[1] = s10 !== s11; has[2] = s01 !== s11; has[3] = s00 !== s01;
        const cnt = has[0] + has[1] + has[2] + has[3];
        if (cnt === 0) continue;
        const X0 = x[j];
        if (has[0]) { ex[0] = X0 + h * tt(v00, v10); ey[0] = Y0; }
        if (has[1]) { ex[1] = X0 + h; ey[1] = Y0 + h * tt(v10, v11); }
        if (has[2]) { ex[2] = X0 + h * tt(v01, v11); ey[2] = Y0 + h; }
        if (has[3]) { ex[3] = X0; ey[3] = Y0 + h * tt(v00, v01); }
        const cxs = [X0, X0 + h, X0 + h, X0], cys = [Y0, Y0, Y0 + h, Y0 + h], vs = [v00, v10, v11, v01];
        // orientation: negative corners must lie left of the segment (iso: saddle, use isolated corner only)
        const emit = (a, b, iso) => {
          let x1 = ex[a], y1 = ey[a], x2 = ex[b], y2 = ey[b];
          const dx = x2 - x1, dy = y2 - y1;
          let score = 0;
          for (let ci = 0; ci < 4; ci++) {
            if (iso !== undefined && ci !== iso) continue;
            const wgt = vs[ci] < 0 ? 1.0 : -1.0;
            const cx = cxs[ci] - x1, cy = cys[ci] - y1;
            score += wgt * (dx * cy - dy * cx);
          }
          if (score < 0) { const tx = x1, ty = y1; x1 = x2; y1 = y2; x2 = tx; y2 = ty; }
          cb(x1, y1, x2, y2);
        };
        if (cnt === 2) {
          let a = -1, b = -1;
          for (let q = 0; q < 4; q++) if (has[q]) { if (a < 0) a = q; else b = q; }
          emit(a, b);
        } else if (cnt === 4) {
          const vc = 0.25 * (v00 + v10 + v11 + v01);
          if ((vc < 0) === s00) { emit(0, 1, 1); emit(2, 3, 3); } else { emit(0, 3, 0); emit(1, 2, 2); }
        }
      }
    }
  }

  function clippedSegLength(x1, y1, x2, y2, R) {
    const dx = x2 - x1, dy = y2 - y1;
    const a = dx * dx + dy * dy;
    const b = 2.0 * (x1 * dx + y1 * dy);
    const c = x1 * x1 + y1 * y1 - R * R;
    const disc = b * b - 4 * a * c;
    if (!(disc > 0 && a > 0)) return 0;
    const sq = Math.sqrt(disc);
    const s1 = (-b - sq) / (2 * a), s2 = (-b + sq) / (2 * a);
    const lo = Math.min(Math.max(s1, 0), 1), hi = Math.min(Math.max(s2, 0), 1);
    return Math.max(hi - lo, 0) * Math.sqrt(a);
  }

  function clippedLength(F, grid, level, R) {
    let s = 0;
    contour(F, grid, level, (x1, y1, x2, y2) => { s += clippedSegLength(x1, y1, x2, y2, R); });
    return s;
  }

  function enclosedArea(F, grid, level) {
    let s = 0;
    contour(F, grid, level, (x1, y1, x2, y2) => { s += x1 * y2 - x2 * y1; });
    return 0.5 * s;
  }

  function sampleMinOnCircle(F, grid, R, m) {
    m = m || 720;
    const { x, h, N } = grid;
    let mn = Infinity;
    for (let q = 0; q < m; q++) {
      const th = q * (2 * PI / m);
      const px = R * Math.cos(th), py = R * Math.sin(th);
      const fx = (px - x[0]) / h, fy = (py - x[0]) / h;
      const j = Math.floor(fx), i = Math.floor(fy);
      const tx = fx - j, ty = fy - i;
      const k = i * N + j;
      const v = F[k] * (1 - tx) * (1 - ty) + F[k + 1] * tx * (1 - ty) + F[k + N] * (1 - tx) * ty + F[k + N + 1] * tx * ty;
      if (v < mn) mn = v;
    }
    return mn;
  }

  /** Perimeter/port-area table versus web distance (see geometry.GrainTable). */
  function GrainTable(geom, R, N, nw) {
    const grid = makeGrid(R, N);
    const phi0 = portSdf(geom, grid);
    const T = fmm(phi0, N, grid.h);
    let wmax = -Infinity;
    const { x } = grid;
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
      if (Math.hypot(x[j], x[i]) <= R && T[i * N + j] > wmax) wmax = T[i * N + j];
    }
    const w = linspace(0, wmax, nw);
    const P = new Float64Array(nw), A = new Float64Array(nw);
    for (let q = 0; q < nw; q++) {
      let lev = w[q];
      if (q === nw - 1) lev -= 0.1 * (w[1] - w[0]);
      P[q] = clippedLength(T, grid, lev, R);
    }
    const A0 = enclosedArea(T, grid, 0.0);
    const disc = PI * R * R;
    let cum = 0;
    for (let q = 1; q < nw; q++) cum += 0.5 * (P[q] + P[q - 1]) * (w[q] - w[q - 1]);
    // mass-conserving correction (see geometry.GrainTable): A0 + int P dw fills the disc exactly
    const pScale = cum > 0 ? (disc - A0) / cum : 1.0;
    const Pc = new Float64Array(nw);
    for (let q = 0; q < nw; q++) Pc[q] = P[q] * pScale;
    let cumc = 0;
    A[0] = Math.min(A0, disc);
    for (let q = 1; q < nw; q++) {
      cumc += 0.5 * (Pc[q] + Pc[q - 1]) * (w[q] - w[q - 1]);
      A[q] = Math.min(A0 + cumc, disc);
    }
    this.geom = geom; this.R = R; this.grid = grid; this.T = T; this.phi0 = phi0;
    this.w = w; this.P_raw = P; this.P = Pc; this.A = A; this.A0 = A0; this.p_scale = pScale;
    this.w_max = wmax; this.w_contact = sampleMinOnCircle(T, grid, R);
    this.A_end_raw = A0 + cum; this.closure = this.A_end_raw / disc - 1;
  }

  function interp(xv, xp, fp) {
    const n = xp.length;
    if (xv <= xp[0]) return fp[0];
    if (xv >= xp[n - 1]) return fp[n - 1];
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xp[mid] <= xv) lo = mid; else hi = mid; }
    const slope = (fp[lo + 1] - fp[lo]) / (xp[lo + 1] - xp[lo]);
    return slope * (xv - xp[lo]) + fp[lo];
  }

  // ------------------------------------------------------------------ propellant + nozzle
  function Propellant(d) {
    this.name = d.name; this.rho = d.rho; this.laws = d.laws.map((l) => l.map(Number));
    this.cstar = d.cstar; this.gamma = d.gamma;
    const g = this.gamma;
    this.Gam = Math.sqrt(g) * Math.pow(2 / (g + 1), (g + 1) / (2 * (g - 1)));
    this.RT = Math.pow(this.cstar * this.Gam, 2);
  }
  Propellant.prototype.lawAt = function (P) {
    for (const l of this.laws) if (P < l[1]) return [l[2], l[3]];
    const L = this.laws[this.laws.length - 1];
    return [L[2], L[3]];
  };
  Propellant.prototype.rate = function (P) {
    if (P <= 0) return 0.0;
    const [a, n] = this.lawAt(P);
    return a * Math.pow(P / P_REF, n);
  };
  Propellant.prototype.equilibriumPc = function (Kn) {
    if (Kn <= 0) return 0.0;
    const K = this.rho * this.cstar * Kn;
    const nl = this.laws.length;
    let P = 0;
    for (let i = 0; i < nl; i++) {
      const [lo, hi, a, n] = this.laws[i];
      P = Math.pow(K * a / Math.pow(P_REF, n), 1.0 / (1.0 - n));
      const loE = i === 0 ? 0.0 : lo;
      const hiE = i === nl - 1 ? Infinity : hi;
      if (loE <= P && P < hiE) return P;
      if (i < nl - 1 && P >= hiE) {
        const a2 = this.laws[i + 1][2], n2 = this.laws[i + 1][3];
        if (K * a2 * Math.pow(hi / P_REF, n2) < hi) return hi;
      }
    }
    return P;
  };

  function areaRatio(M, g) {
    return (1.0 / M) * Math.pow((2 / (g + 1)) * (1 + 0.5 * (g - 1) * M * M), (g + 1) / (2 * (g - 1)));
  }
  function machFromArea(eps, g, supersonic) {
    let lo = supersonic ? 1.0 : 1e-9, hi = supersonic ? 100.0 : 1.0;
    for (let it = 0; it < 200; it++) {
      const mid = 0.5 * (lo + hi);
      const f = areaRatio(mid, g) - eps;
      if (supersonic) { if (f < 0) lo = mid; else hi = mid; } else { if (f < 0) hi = mid; else lo = mid; }
    }
    return 0.5 * (lo + hi);
  }
  const pRatio = (M, g) => Math.pow(1 + 0.5 * (g - 1) * M * M, -g / (g - 1));

  function Nozzle(eps, g) {
    this.eps = eps; this.g = g;
    this.pe_pc = pRatio(machFromArea(eps, g, true), g);
    this.psub = pRatio(machFromArea(eps, g, false), g);
    this.cf_mom = Math.sqrt(2 * g * g / (g - 1) * Math.pow(2 / (g + 1), (g + 1) / (g - 1)) *
      (1 - Math.pow(this.pe_pc, (g - 1) / g)));
  }
  Nozzle.prototype.cf = function (pc, pa) { return this.cf_mom + (this.pe_pc - pa / pc) * this.eps; };
  Nozzle.prototype.flow = function (pc, pa, At, cstar, RT) {
    if (pc <= pa) return [0.0, 0.0];
    const g = this.g;
    if (pa / pc <= this.psub) {
      const mdot = pc * At / cstar;
      return [mdot, Math.max(this.cf(pc, pa) * pc * At, 0.0)];
    }
    const x = pa / pc, Ae = At * this.eps;
    const mdot = Ae * pc * Math.sqrt(2 * g / ((g - 1) * RT)) *
      Math.sqrt(Math.max(Math.pow(x, 2 / g) - Math.pow(x, (g + 1) / g), 0.0));
    const ue = Math.sqrt(2 * g / (g - 1) * RT * Math.max(1 - Math.pow(x, (g - 1) / g), 0.0));
    return [mdot, mdot * ue];
  };

  function motorClass(I) {
    if (I <= 0) return '-';
    if (I <= 0.625) return '1/4A';
    if (I <= 1.25) return '1/2A';
    const k = Math.max(0, Math.ceil(Math.log2(I / 2.5) - 1e-12));
    return k < 26 ? String.fromCharCode(65 + k) : 'beyond Z';
  }

  // ------------------------------------------------------------------ motor
  const tableCache = new Map();
  function tableKey(geom, R, N, nw) { return JSON.stringify([geom, R, N, nw]); }
  function grainTable(geom, R, N, nw) {
    const key = tableKey(geom, R, N, nw);
    if (!tableCache.has(key)) {
      if (tableCache.size > 64) tableCache.clear();
      tableCache.set(key, new GrainTable(geom, R, N, nw));
    }
    return tableCache.get(key);
  }
  /** Distinct geometries of a spec, so a UI can build tables one per animation frame. */
  function pendingTables(spec) {
    const N = spec.grid.N, nw = spec.grid.nw;
    return spec.segments.map((s) => s.geom).filter((g, i, arr) =>
      !tableCache.has(tableKey(g, spec.R, N, nw)) && arr.findIndex((h) => JSON.stringify(h) === JSON.stringify(g)) === i);
  }

  function Motor(spec) {
    this.spec = spec;
    this.prop = new Propellant(spec.prop);
    this.nozzle = new Nozzle(spec.eps, this.prop.gamma);
    this.R = spec.R; this.pa = spec.pa !== undefined ? spec.pa : 101325.0; this.Dt = spec.Dt;
    const N = spec.grid.N, nw = spec.grid.nw;
    this.segs = [];
    for (const seg of spec.segments) {
      const tab = grainTable(seg.geom, this.R, N, nw);
      const cnt = seg.count || 1;
      for (let c = 0; c < cnt; c++) {
        const L0 = seg.L, ends = seg.ends || 0;
        const wEnd = Math.min(tab.w_max, ends ? L0 / ends : Infinity);
        this.segs.push({ tab, L0, ends, w_end: wEnd });
      }
    }
    this.K = this.segs.length;
    this.disc = PI * this.R * this.R;
    this.V_case = this.disc * this.segs.reduce((s, g) => s + g.L0, 0) + (spec.V_extra || 0);
    this.erosive = spec.erosive || { on: false };
    this.throat_erosion = spec.throat_erosion || 0.0;
    this.m_prop = this.prop.rho * this.segs.reduce((s, g) => s + this.vprop(g, 0.0), 0);
    this.w_web = Math.min(...this.segs.map((g) => Math.min(g.tab.w_contact, g.ends ? g.L0 / g.ends : Infinity)));
    this.dw_tab = Math.min(...this.segs.map((g) => g.tab.w[1]));
  }
  Motor.prototype.segState = function (sg, w) {
    const Lb = sg.L0 - sg.ends * w;
    if (w >= sg.w_end || Lb <= 0) return [0.0, this.disc, 0.0];
    const P = interp(w, sg.tab.w, sg.tab.P), A = interp(w, sg.tab.w, sg.tab.A);
    return [P * Lb + sg.ends * (this.disc - A), A, P];
  };
  Motor.prototype.vprop = function (sg, w) {
    const Lb = sg.L0 - sg.ends * w;
    if (w >= sg.w_end || Lb <= 0) return 0.0;
    const A = interp(w, sg.tab.w, sg.tab.A);
    return Math.max(this.disc - A, 0.0) * Lb;
  };
  Motor.prototype.abTotal = function (w) {
    let s = 0;
    for (const sg of this.segs) s += this.segState(sg, w)[0];
    return s;
  };
  Motor.prototype.At = function (t) { return 0.25 * PI * Math.pow(this.Dt + this.throat_erosion * (t || 0), 2); };

  Motor.prototype.quasiSteady = function (nPts) {
    nPts = nPts || 400;
    const wAll = Math.max(...this.segs.map((g) => g.w_end));
    const wf = linspace(0.0, wAll, nPts);
    const n = nPts - 1;
    const At = this.At(0.0), prop = this.prop, noz = this.nozzle;
    const w = [], t = [], Pc = [], F = [], Kn = [];
    for (let i = 0; i < n; i++) {
      const kn = this.abTotal(wf[i]) / At;
      const p = prop.equilibriumPc(kn);
      w.push(wf[i]); Kn.push(kn); Pc.push(p);
      F.push(noz.flow(p, this.pa, At, prop.cstar, prop.RT)[1]);
    }
    t.push(0);
    let I = 0;
    for (let i = 1; i < n; i++) {
      const r1 = prop.rate(Pc[i]), r0 = prop.rate(Pc[i - 1]);
      t.push(t[i - 1] + 0.5 * (1 / r1 + 1 / r0) * (w[i] - w[i - 1]));
      I += 0.5 * (F[i] + F[i - 1]) * (t[i] - t[i - 1]);
    }
    let pmax = -Infinity, pmin = Infinity, kmax = -Infinity, kmin = Infinity, first = null, last = null;
    for (let i = 0; i < n; i++) {
      if (w[i] <= 0.97 * this.w_web) {
        pmax = Math.max(pmax, Pc[i]); pmin = Math.min(pmin, Pc[i]);
        kmax = Math.max(kmax, Kn[i]); kmin = Math.min(kmin, Kn[i]);
        if (first === null) first = Kn[i];
        last = Kn[i];
      }
    }
    return { w, t, Pc, F, Kn, I, Pc_max: Math.max(...Pc), neutrality_pc: pmax / pmin, neutrality_kn: kmax / kmin,
      progressivity: last / first };
  };

  Motor.prototype.rhs = function (t, y) {
    const prop = this.prop, noz = this.nozzle;
    const Pc = y[0];
    const r0 = prop.rate(Pc);
    const At = this.At(t);
    const K = this.K;
    const states = new Array(K);
    for (let k = 0; k < K; k++) states[k] = this.segState(this.segs[k], y[1 + k]);
    const rates = new Array(K).fill(r0);
    if (this.erosive.on) {
      const alpha = this.erosive.alpha, beta = this.erosive.beta;
      let up = 0.0;
      for (let k = 0; k < K; k++) {
        const [Ab, A, P] = states[k];
        const mk = prop.rho * Ab * r0;
        if (Ab > 0 && P > 0 && A > 0) {
          const G = (up + 0.5 * mk) / A;
          const Dh = 4.0 * A / P;
          if (G > 0) rates[k] = r0 + alpha * Math.pow(G, 0.8) * Math.pow(Dh, -0.2) * Math.exp(-beta * r0 * prop.rho / G);
        }
        up += mk;
      }
    }
    let mgen = 0.0, dV = 0.0, vprop = 0.0, AbTot = 0.0;
    for (let k = 0; k < K; k++) {
      const Ab = states[k][0];
      mgen += prop.rho * Ab * rates[k];
      dV += Ab * rates[k];
      vprop += this.vprop(this.segs[k], y[1 + k]);
      AbTot += Ab;
    }
    const V = this.V_case - vprop;
    const [mout, F] = noz.flow(Pc, this.pa, At, prop.cstar, prop.RT);
    const rhoG = Pc / prop.RT;
    const dy = new Array(K + 1);
    dy[0] = prop.RT / V * (mgen - rhoG * dV - mout);
    for (let k = 0; k < K; k++) dy[1 + k] = states[k][0] > 0 ? rates[k] : 0.0;
    return [dy, { F, Kn: AbTot / At, V, mout, r: r0, At }];
  };

  Motor.prototype.transient = function (dtFrac, maxSteps) {
    dtFrac = dtFrac || 0.2; maxSteps = maxSteps || 400000;
    const prop = this.prop, K = this.K;
    let y = [this.pa].concat(new Array(K).fill(0.0));
    let t = 0.0, Pmax = this.pa;
    const out = { t: [], Pc: [], F: [], Kn: [], w: [], mdot: [] };
    const add = (a, b, s) => a.map((v, i) => v + s * b[i]);
    for (let step = 0; step < maxSteps; step++) {
      const [dy, aux] = this.rhs(t, y);
      out.t.push(t); out.Pc.push(y[0]); out.F.push(aux.F); out.Kn.push(aux.Kn); out.w.push(y[1]); out.mdot.push(aux.mout);
      Pmax = Math.max(Pmax, y[0]);
      let burned = true;
      for (let k = 0; k < K; k++) if (y[1 + k] < this.segs[k].w_end) { burned = false; break; }
      if (burned && y[0] - this.pa < 0.01 * (Pmax - this.pa)) break;
      const tau = aux.V * prop.cstar / (prop.RT * aux.At);
      let dt = dtFrac * tau;
      if (!burned) dt = Math.min(dt, 0.5 * this.dw_tab / Math.max(aux.r, 1e-6));
      const k1 = dy;
      const k2 = this.rhs(t + 0.5 * dt, add(y, k1, 0.5 * dt))[0];
      const k3 = this.rhs(t + 0.5 * dt, add(y, k2, 0.5 * dt))[0];
      const k4 = this.rhs(t + dt, add(y, k3, dt))[0];
      y = y.map((v, i) => v + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
      y[0] = Math.max(y[0], this.pa * 0.5);
      t += dt;
    }
    return out;
  };

  Motor.prototype.metrics = function (tr, qs) {
    const t = tr.t, F = tr.F, Pc = tr.Pc;
    let I = 0, Fmax = -Infinity;
    for (let i = 1; i < t.length; i++) I += 0.5 * (F[i] + F[i - 1]) * (t[i] - t[i - 1]);
    for (const f of F) Fmax = Math.max(Fmax, f);
    let i0 = -1, i1 = -1;
    for (let i = 0; i < F.length; i++) if (F[i] >= 0.1 * Fmax) { if (i0 < 0) i0 = i; i1 = i; }
    const t0 = t[i0], t1 = t[i1], tb = t1 - t0;
    const Favg = I / tb;
    const meanBetween = (a, b) => {
      let s = 0, ta = null, tl = null, pt = null, pf = null;
      for (let i = 0; i < t.length; i++) {
        if (t[i] >= a && t[i] <= b) {
          if (ta === null) ta = t[i];
          if (pt !== null) s += 0.5 * (F[i] + pf) * (t[i] - pt);
          pt = t[i]; pf = F[i]; tl = t[i];
        }
      }
      return s / (tl - ta);
    };
    const shape = meanBetween(t0 + 2 * tb / 3, t1) / meanBetween(t0, t0 + tb / 3);
    const cls = motorClass(I);
    const m = { I_total: I, F_max: Fmax, F_avg: Favg, t_burn: tb, t_start: t0, t_end: t1,
      Pc_max: Math.max(...Pc), m_prop: this.m_prop, Isp: I / (this.m_prop * G0), class: cls,
      designation: cls.length === 1 ? cls + Math.round(Favg) : cls,
      Kn_initial: tr.Kn[0], Kn_max: Math.max(...tr.Kn), At: this.At(0), w_web: this.w_web, shape_ratio: shape };
    if (qs) Object.assign(m, { neutrality_pc: qs.neutrality_pc, neutrality_kn: qs.neutrality_kn,
      progressivity: qs.progressivity, qs_I: qs.I, qs_Pc_max: qs.Pc_max });
    return m;
  };

  const SRM = { makeGrid, portSdf, fmm, contour, clippedLength, enclosedArea, GrainTable, grainTable, pendingTables,
    Propellant, Nozzle, Motor, motorClass, interp, starVertices, P_REF, G0 };
  if (typeof module !== 'undefined' && module.exports) module.exports = SRM;
  else root.SRM = SRM;
})(typeof window !== 'undefined' ? window : this);
