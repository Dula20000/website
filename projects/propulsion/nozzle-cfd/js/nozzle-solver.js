/* Quasi-1D Euler nozzle solver — JavaScript port of code/nozzlecfd (solver.py, exact.py, atmosphere.py).
   Same scheme: cell-centred FV, MUSCL (van Albada) on primitives, HLLC with Davis/Einfeldt bounds,
   SSP-RK3, reservoir inflow, back-pressure outflow with supersonic/subsonic switch.
   Non-dimensional: p0 = rho0 = 1 (R T0 = 1). Works in the browser (window.NozzleCFD) and in node (require). */
(function (root) {
  'use strict';
  const NG = 2;

  // ---------------------------------------------------------------- geometry
  function bellNozzle(eps, CR = 4.0, Lc = 0.3, Ld = 1.0) {
    return {
      kind: 'bell', eps, x_throat: Lc, L: Lc + Ld,
      A(x) {
        if (x <= Lc) { const s = Math.max((Lc - x) / Lc, 0); return 1 + (CR - 1) * s * s; }
        const t = Math.min(Math.max((x - Lc) / Ld, 0), 1);
        return 1 + (eps - 1) * (3 * t * t - 2 * t * t * t);
      },
    };
  }
  // Anderson family: A = 1 + 2.2 (x-1.5)^2 upstream, 1 + k (x-1.5)^2 downstream, k = (eps-1)/2.25.
  // eps = 5.95 gives Anderson's nozzle exactly.
  function andersonNozzle(eps = 5.95) {
    const k = eps === 5.95 ? 2.2 : (eps - 1) / 2.25;
    return {
      kind: 'anderson', eps, x_throat: 1.5, L: 3.0,
      A(x) { const d = x - 1.5; return x <= 1.5 ? 1 + 2.2 * d * d : 1 + k * d * d; },
    };
  }

  // ---------------------------------------------------------------- exact relations
  function areaMach(M, g) {
    return (1 / M) * Math.pow((2 / (g + 1)) * (1 + 0.5 * (g - 1) * M * M), (g + 1) / (2 * (g - 1)));
  }
  function bisect(f, a, b, n = 200) {
    let fa = f(a);
    for (let i = 0; i < n; i++) {
      const m = 0.5 * (a + b), fm = f(m);
      if (fm === 0 || (b - a) < 1e-15 * Math.max(1, Math.abs(m))) return m;
      if ((fa < 0) === (fm < 0)) { a = m; fa = fm; } else b = m;
    }
    return 0.5 * (a + b);
  }
  function machFromArea(ar, g, sup) {
    if (ar <= 1) return 1;
    const f = (M) => areaMach(M, g) - ar;
    if (sup) { let hi = 2; while (f(hi) < 0) hi *= 2; return bisect(f, 1, hi); }
    return bisect(f, 1e-12, 1);
  }
  const TT0 = (M, g) => 1 / (1 + 0.5 * (g - 1) * M * M);
  const pp0 = (M, g) => Math.pow(TT0(M, g), g / (g - 1));
  function shockPR(M, g) { return 1 + 2 * g / (g + 1) * (M * M - 1); }
  function p02p01(M1, g) {
    const s = M1 * M1;
    return Math.pow((g + 1) * s / (2 + (g - 1) * s), g / (g - 1)) *
           Math.pow(2 * g / (g + 1) * s - (g - 1) / (g + 1), -1 / (g - 1));
  }
  function critical(eps, g) {
    const Msub = machFromArea(eps, g, false), Msup = machFromArea(eps, g, true);
    const psup = pp0(Msup, g);
    return { p_sub: pp0(Msub, g), p_sup: psup, p_nse: psup * shockPR(Msup, g), Me_sup: Msup, Me_sub: Msub };
  }
  function shockAreaRatio(pb, eps, g) {
    const c = critical(eps, g);
    const gm = g - 1, k = Math.pow(2 / (g + 1), (g + 1) / gm), r = 1 / (pb * eps);
    const Me2 = -1 / gm + Math.sqrt(1 / (gm * gm) + (2 / gm) * k * r * r);
    const r0 = pb * Math.pow(1 + 0.5 * gm * Me2, g / gm);
    if (r0 >= 1) return { As: 1, r0: 1 };
    const M1 = bisect((M) => p02p01(M, g) - r0, 1 + 1e-12, c.Me_sup * (1 + 1e-9) + 1e-9);
    return { As: areaMach(M1, g), M1, r0 };
  }
  /** Exact quasi-1D solution on points xs for geometry geom and back pressure pb/p0. */
  function exactNozzle(xs, geom, pb, g) {
    const Ath = geom.A(geom.x_throat), Ae = geom.A(geom.L), eps = Ae / Ath;
    const c = critical(eps, g);
    const M = new Float64Array(xs.length), p = new Float64Array(xs.length);
    let regime, xShock = null, AsR = null;
    if (pb >= c.p_sub) {
      regime = 'subsonic';
      const Me = pb < 1 ? Math.sqrt(2 / (g - 1) * (Math.pow(pb, -(g - 1) / g) - 1)) : 0;
      const Astar = Me > 0 ? Ae / areaMach(Me, g) : Infinity;
      xs.forEach((x, i) => { M[i] = isFinite(Astar) ? machFromArea(geom.A(x) / Astar, g, false) : 0; p[i] = pp0(M[i], g); });
    } else {
      let r0 = 1, Astar2 = Ath;
      if (pb > c.p_nse) {
        regime = 'shock';
        const s = shockAreaRatio(pb, eps, g);
        r0 = s.r0; AsR = s.As; Astar2 = Ath / r0;
        xShock = s.As <= 1 + 1e-14 ? geom.x_throat : bisect((x) => geom.A(x) - s.As * Ath, geom.x_throat, geom.L);
      } else regime = 'supersonic';
      xs.forEach((x, i) => {
        const a = geom.A(x);
        if (x <= geom.x_throat) { M[i] = machFromArea(a / Ath, g, false); p[i] = pp0(M[i], g); }
        else if (xShock === null || x < xShock) { M[i] = machFromArea(a / Ath, g, true); p[i] = pp0(M[i], g); }
        else { M[i] = machFromArea(a / Astar2, g, false); p[i] = r0 * pp0(M[i], g); }
      });
    }
    return { M, p, regime, x_shock: xShock, As: AsR, crit: c };
  }

  // ---------------------------------------------------------------- atmosphere (US 1976, 0-51 km)
  const ATM = (function () {
    const G0 = 9.80665, R = 8314.32 / 28.9644, RE = 6356.766e3;
    const layers = [[0, -6.5e-3], [11000, 0], [20000, 1e-3], [32000, 2.8e-3], [47000, 0]];
    const Tb = [288.15], Pb = [101325];
    for (let i = 0; i < 4; i++) {
      const [h0, L] = layers[i], h1 = layers[i + 1][0], T0 = Tb[i], p0 = Pb[i], T1 = T0 + L * (h1 - h0);
      Pb.push(L === 0 ? p0 * Math.exp(-G0 * (h1 - h0) / (R * T0)) : p0 * Math.pow(T1 / T0, -G0 / (R * L)));
      Tb.push(T1);
    }
    return function (z) {
      const h = RE * z / (RE + z);
      let i = 0; while (i < 4 && h >= layers[i + 1][0]) i++;
      const [h0, L] = layers[i];
      if (L === 0) return { T: Tb[i], p: Pb[i] * Math.exp(-G0 * (h - h0) / (R * Tb[i])) };
      const T = Tb[i] + L * (h - h0);
      return { T, p: Pb[i] * Math.pow(T / Tb[i], -G0 / (R * L)) };
    };
  })();

  // ---------------------------------------------------------------- limiter / flux
  function vanAlbada(dm, dp) {
    if (dm * dp <= 0) return 0;
    const e = 1e-12;
    return (dm * (dp * dp + e) + dp * (dm * dm + e)) / (dm * dm + dp * dp + 2 * e);
  }
  function minmod(dm, dp) { return dm * dp <= 0 ? 0 : (Math.abs(dm) < Math.abs(dp) ? dm : dp); }

  // HLLC flux, writes into out[0..2]
  function hllc(rL, uL, pL, rR, uR, pR, g, out) {
    const EL = pL / (g - 1) + 0.5 * rL * uL * uL, ER = pR / (g - 1) + 0.5 * rR * uR * uR;
    const cL = Math.sqrt(g * pL / rL), cR = Math.sqrt(g * pR / rR);
    const sL = Math.sqrt(rL), sR = Math.sqrt(rR);
    const HL = (EL + pL) / rL, HR = (ER + pR) / rR;
    const ut = (sL * uL + sR * uR) / (sL + sR), Ht = (sL * HL + sR * HR) / (sL + sR);
    const ct = Math.sqrt(Math.max((g - 1) * (Ht - 0.5 * ut * ut), 1e-300));
    const SL = Math.min(Math.min(uL - cL, uR - cR), ut - ct);
    const SR = Math.max(Math.max(uR + cR, uL + cL), ut + ct);
    const Ss = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR));
    if (SL >= 0) { out[0] = rL * uL; out[1] = rL * uL * uL + pL; out[2] = (EL + pL) * uL; return; }
    if (SR <= 0) { out[0] = rR * uR; out[1] = rR * uR * uR + pR; out[2] = (ER + pR) * uR; return; }
    if (Ss >= 0) {
      const f = rL * (SL - uL) / (SL - Ss);
      out[0] = rL * uL + SL * (f - rL);
      out[1] = rL * uL * uL + pL + SL * (f * Ss - rL * uL);
      out[2] = (EL + pL) * uL + SL * (f * (EL / rL + (Ss - uL) * (Ss + pL / (rL * (SL - uL)))) - EL);
    } else {
      const f = rR * (SR - uR) / (SR - Ss);
      out[0] = rR * uR + SR * (f - rR);
      out[1] = rR * uR * uR + pR + SR * (f * Ss - rR * uR);
      out[2] = (ER + pR) * uR + SR * (f * (ER / rR + (Ss - uR) * (Ss + pR / (rR * (SR - uR)))) - ER);
    }
  }

  // ---------------------------------------------------------------- solver
  class Solver {
    /**
     * @param {object} geom  from bellNozzle / andersonNozzle
     * @param {number} N     interior cells
     * @param {number} g     gamma
     * @param {number} pb    back pressure / p0
     * @param {object} opt   {cfl, localDt, init: 'ramp'|'rest', limiter: 'vanalbada'|'minmod'|'none'}
     */
    constructor(geom, N, g, pb, opt = {}) {
      this.geom = geom; this.N = N; this.g = g; this.pb = pb;
      this.cfl = opt.cfl ?? 0.8; this.localDt = opt.localDt ?? true;
      this.limiter = opt.limiter ?? 'vanalbada';
      const L = geom.L; this.dx = L / N;
      this.xf = new Float64Array(N + 1); this.Af = new Float64Array(N + 1);
      this.x = new Float64Array(N); this.A = new Float64Array(N);
      for (let i = 0; i <= N; i++) { this.xf[i] = L * i / N; this.Af[i] = geom.A(this.xf[i]); }
      for (let i = 0; i < N; i++) { this.x[i] = 0.5 * (this.xf[i] + this.xf[i + 1]); this.A[i] = geom.A(this.x[i]); }
      const n = N + 2 * NG;
      this.W = [new Float64Array(n), new Float64Array(n), new Float64Array(n)];  // primitives with ghosts
      this.Q = [new Float64Array(N), new Float64Array(N), new Float64Array(N)];
      this.Q0 = this.Q.map((a) => new Float64Array(N)); this.Q1 = this.Q.map((a) => new Float64Array(N));
      this.K = this.Q.map(() => new Float64Array(N));
      this.FA = [new Float64Array(N + 1), new Float64Array(N + 1), new Float64Array(N + 1)];
      this.dt = new Float64Array(N);
      this.flux = new Float64Array(3);
      this.slope = [new Float64Array(n), new Float64Array(n), new Float64Array(n)];
      this.iter = 0; this.res0 = null; this.res = 1; this.exitSupersonic = false;
      this.init(opt.init ?? 'rest');
    }
    init(kind) {
      const { g, N, geom } = this; const xt = geom.x_throat, L = geom.L;
      for (let i = 0; i < N; i++) {
        let rho, u, p;
        if (kind === 'rest') { rho = 1; u = 0; p = 1; }
        else {
          const x = this.x[i];
          const M = x <= xt ? 0.2 + 0.8 * x / xt : 1 + (x - xt) / (L - xt);
          const T = 1 / (1 + 0.5 * (g - 1) * M * M);
          p = Math.pow(T, g / (g - 1)); rho = p / T; u = M * Math.sqrt(g * T);
        }
        this.Q[0][i] = rho; this.Q[1][i] = rho * u; this.Q[2][i] = p / (g - 1) + 0.5 * rho * u * u;
      }
      this.iter = 0; this.res0 = null; this.res = 1; this.history = [];
    }
    fillPrimitives(Q) {
      const { g, N, W } = this; const [R, U, P] = W;
      for (let i = 0; i < N; i++) {
        const r = Q[0][i], u = Q[1][i] / r;
        R[i + NG] = r; U[i + NG] = u; P[i + NG] = (g - 1) * (Q[2][i] - 0.5 * r * u * u);
      }
      // inflow ghosts: velocity extrapolated with a minmod-limited slope, kept subsonic; isentropic from reservoir
      const mm = (a, b) => (a * b <= 0 ? 0 : (Math.abs(a) < Math.abs(b) ? a : b));
      const su = mm(U[NG] - U[NG + 1], U[NG + 1] - U[NG + 2]);
      const umax = 0.95 * Math.sqrt(2 * g / (g + 1));
      for (let k = 1; k <= 2; k++) {
        const ug = Math.min(Math.max(U[NG] + k * su, 0), umax);
        const T = 1 - 0.5 * (g - 1) / g * ug * ug;
        const pg = Math.pow(T, g / (g - 1));
        R[NG - k] = pg / T; U[NG - k] = ug; P[NG - k] = pg;
      }
      // outflow ghosts
      const a = NG + N - 1;
      const sr = mm(R[a] - R[a - 1], R[a - 1] - R[a - 2]);
      const sv = mm(U[a] - U[a - 1], U[a - 1] - U[a - 2]);
      const sp = mm(P[a] - P[a - 1], P[a - 1] - P[a - 2]);
      const rf = Math.max(R[a] + 0.5 * sr, 1e-10), uf = U[a] + 0.5 * sv, pf = Math.max(P[a] + 0.5 * sp, 1e-10);
      const Mf = uf / Math.sqrt(g * pf / rf);
      let up = -1;
      if (Mf >= 1) this.exitSupersonic = this.pb <= pf * shockPR(Mf, g);
      else {
        this.exitSupersonic = false;
        // captured shock in the last cells that cannot stand at this pb: treat exit as supersonic
        for (const j of [a - 1, a - 2, a - 3]) {
          const Mj = U[j] / Math.sqrt(g * P[j] / R[j]);
          if (Mj >= 1) { if (this.pb <= P[j] * shockPR(Mj, g)) { this.exitSupersonic = true; up = j; } break; }
        }
      }
      for (let k = 1; k <= 2; k++) {
        const j = a + k;
        if (up >= 0) { R[j] = R[up]; U[j] = U[up]; P[j] = P[up]; continue; }
        R[j] = Math.max(R[a] + k * sr, 1e-10);
        U[j] = U[a] + k * sv;
        P[j] = this.exitSupersonic ? Math.max(P[a] + k * sp, 1e-10) : this.pb;
      }
    }
    rhs(Q, K) {
      const { g, N, W, FA, Af, flux, slope } = this;
      this.fillPrimitives(Q);
      const n = N + 2 * NG;
      const lim = this.limiter === 'minmod' ? minmod : this.limiter === 'none' ? () => 0 : vanAlbada;
      for (let v = 0; v < 3; v++) {
        const w = W[v], s = slope[v];
        for (let j = 1; j < n - 1; j++) s[j] = lim(w[j] - w[j - 1], w[j + 1] - w[j]);
      }
      for (let f = 0; f <= N; f++) {
        const l = f + NG - 1, r = f + NG;   // global cells either side of face f
        let rL = W[0][l] + 0.5 * slope[0][l], uL = W[1][l] + 0.5 * slope[1][l], pL = W[2][l] + 0.5 * slope[2][l];
        let rR = W[0][r] - 0.5 * slope[0][r], uR = W[1][r] - 0.5 * slope[1][r], pR = W[2][r] - 0.5 * slope[2][r];
        if (rL <= 0 || pL <= 0 || rR <= 0 || pR <= 0) {
          rL = W[0][l]; uL = W[1][l]; pL = W[2][l]; rR = W[0][r]; uR = W[1][r]; pR = W[2][r];
        }
        hllc(rL, uL, pL, rR, uR, pR, g, flux);
        FA[0][f] = flux[0] * Af[f]; FA[1][f] = flux[1] * Af[f]; FA[2][f] = flux[2] * Af[f];
      }
      for (let i = 0; i < N; i++) {
        const vol = this.A[i] * this.dx;
        K[0][i] = -(FA[0][i + 1] - FA[0][i]) / vol;
        K[1][i] = (-(FA[1][i + 1] - FA[1][i]) + W[2][i + NG] * (Af[i + 1] - Af[i])) / vol;
        K[2][i] = -(FA[2][i + 1] - FA[2][i]) / vol;
      }
    }
    computeDt() {
      const { g, N, Q, dt } = this; let mn = Infinity;
      for (let i = 0; i < N; i++) {
        const r = Q[0][i], u = Q[1][i] / r, p = (g - 1) * (Q[2][i] - 0.5 * r * u * u);
        dt[i] = this.cfl * this.dx / (Math.abs(u) + Math.sqrt(g * p / r));
        if (dt[i] < mn) mn = dt[i];
      }
      if (!this.localDt) dt.fill(mn);
    }
    /** One SSP-RK3 step; returns relative RMS density residual. */
    step() {
      const { N, Q, Q0, Q1, K, dt } = this;
      this.computeDt();
      for (let v = 0; v < 3; v++) Q0[v].set(Q[v]);
      this.rhs(Q0, K);
      let s = 0; for (let i = 0; i < N; i++) s += K[0][i] * K[0][i];
      const r = Math.sqrt(s / N);
      for (let v = 0; v < 3; v++) for (let i = 0; i < N; i++) Q1[v][i] = Q0[v][i] + dt[i] * K[v][i];
      this.rhs(Q1, K);
      for (let v = 0; v < 3; v++) for (let i = 0; i < N; i++) Q1[v][i] = 0.75 * Q0[v][i] + 0.25 * (Q1[v][i] + dt[i] * K[v][i]);
      this.rhs(Q1, K);
      for (let v = 0; v < 3; v++) for (let i = 0; i < N; i++) Q[v][i] = Q0[v][i] / 3 + 2 / 3 * (Q1[v][i] + dt[i] * K[v][i]);
      this.iter++;
      if (this.res0 === null) this.res0 = r;
      this.res = r / this.res0;
      if (!isFinite(r)) throw new Error('diverged');
      return this.res;
    }
    run(maxIter = 200000, tol = 1e-11) {
      while (this.iter < maxIter) { if (this.step() < tol) return true; }
      return false;
    }
    /** Primitive fields, Mach, exit state and face fluxes at the current iterate. */
    state() {
      const { g, N } = this;
      this.rhs(this.Q, this.K);
      const rho = new Float64Array(N), u = new Float64Array(N), p = new Float64Array(N), M = new Float64Array(N);
      for (let i = 0; i < N; i++) {
        rho[i] = this.Q[0][i]; u[i] = this.Q[1][i] / rho[i];
        p[i] = (g - 1) * (this.Q[2][i] - 0.5 * rho[i] * u[i] * u[i]); M[i] = u[i] / Math.sqrt(g * p[i] / rho[i]);
      }
      const ex = (q) => 1.5 * q[N - 1] - 0.5 * q[N - 2];
      const re = ex(rho), ue = ex(u), pe = ex(p);
      return { x: this.x, A: this.A, rho, u, p, M, exit: { rho: re, u: ue, p: pe, M: ue / Math.sqrt(g * pe / re) },
               mdotIn: this.FA[0][0], mdotOut: this.FA[0][N], momExit: this.FA[1][N], exitSupersonic: this.exitSupersonic };
    }
  }

  /** Shock position from sonic crossing (M >= 1 to M < 1) downstream of the throat. */
  function shockLocation(x, M, xt) {
    let found = null;
    for (let i = 0; i < x.length - 1; i++) {
      if (x[i] > xt && M[i] >= 1 && M[i + 1] < 1) found = x[i] + (1 - M[i]) * (x[i + 1] - x[i]) / (M[i + 1] - M[i]);
    }
    return found;
  }
  const schmucker = (Me) => Math.pow(1.88 * Me - 1, -0.64);

  const api = { bellNozzle, andersonNozzle, areaMach, machFromArea, pp0, critical, shockAreaRatio, exactNozzle,
                atmosphere: ATM, Solver, shockLocation, schmucker, mdotStar: (g) => Math.sqrt(g) * Math.pow(2 / (g + 1), (g + 1) / (2 * (g - 1))) };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.NozzleCFD = api;
})(typeof window !== 'undefined' ? window : this);
