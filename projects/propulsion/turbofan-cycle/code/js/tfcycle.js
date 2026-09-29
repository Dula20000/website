/* tfcycle.js -- JavaScript port of the Python `tfcycle` package (same equations,
   same algorithms): 1976 US Standard Atmosphere; variable-property gas model
   (NASA 7-coefficient polynomials for air and Jet-A combustion products) or
   Mattingly's two-gas / ideal model; parametric cycle analysis of a two-spool
   separate-exhaust turbofan with a simple turbine cooling-air bleed; closed-form
   ideal turbofan; TSFC-optimal fan pressure ratio; fixed-geometry off-design
   matching; T-s path through the station states.
   Works in the browser (window.TF) and in Node (module.exports). */
(function (root) {
  'use strict';
  const FT = 0.3048, G0 = 9.80665, R_AIR = 8314.32 / 28.9644, R_EARTH = 6356766.0;
  const TSFC_LB = 3600 * 9.80665;

  // ---------------- atmosphere ----------------
  const LH = [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852];
  const LL = [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002];
  const TB = [288.15], PB = [101325];
  for (let i = 0; i < LL.length; i++) {
    const dh = LH[i + 1] - LH[i], t1 = TB[i] + LL[i] * dh;
    const p1 = LL[i] === 0 ? PB[i] * Math.exp(-G0 * dh / (R_AIR * TB[i]))
                           : PB[i] * Math.pow(t1 / TB[i], -G0 / (R_AIR * LL[i]));
    TB.push(t1); PB.push(p1);
  }
  function atmosphere(h, geometric) {
    const H = geometric ? R_EARTH * h / (R_EARTH + h) : +h;
    if (H < -1000 || H > LH[LH.length - 1]) throw new Error('altitude out of range');
    let i = 0;
    while (i < LL.length - 1 && H > LH[i + 1]) i++;
    const dh = H - LH[i], T = TB[i] + LL[i] * dh;
    const P = LL[i] === 0 ? PB[i] * Math.exp(-G0 * dh / (R_AIR * TB[i]))
                          : PB[i] * Math.pow(T / TB[i], -G0 / (R_AIR * LL[i]));
    return { T, P, rho: P / (R_AIR * T), a: Math.sqrt(1.4 * R_AIR * T) };
  }

  // ---------------- thermodynamics (port of thermo.py) ----------------
  const RU = 8314.462618, T_REF = 298.15, T_MID = 1000.0;
  const NASA = {
    N2: [28.0134, [3.298677, 1.4082404e-3, -3.963222e-6, 5.641515e-9, -2.444854e-12, -1020.8999, 3.950372],
         [2.92664, 1.4879768e-3, -5.68476e-7, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528]],
    O2: [31.9988, [3.78245636, -2.99673416e-3, 9.84730201e-6, -9.68129509e-9, 3.24372837e-12, -1063.94356, 3.65767573],
         [3.28253784, 1.48308754e-3, -7.57966669e-7, 2.09470555e-10, -2.16717794e-14, -1088.45772, 5.45323129]],
    Ar: [39.948, [2.5, 0, 0, 0, 0, -745.375, 4.366], [2.5, 0, 0, 0, 0, -745.375, 4.366]],
    CO2: [44.0095, [2.35677352, 8.98459677e-3, -7.12356269e-6, 2.45919022e-9, -1.43699548e-13, -48371.9697, 9.90105222],
          [3.85746029, 4.41437026e-3, -2.21481404e-6, 5.23490188e-10, -4.72084164e-14, -48759.166, 2.27163806]],
    H2O: [18.01528, [4.19864056, -2.0364341e-3, 6.52040211e-6, -5.48797062e-9, 1.77197817e-12, -30293.7267, -0.849032208],
          [3.03399249, 2.17691804e-3, -1.64072518e-7, -9.7041987e-11, 1.68200992e-14, -30004.2971, 4.9667701]],
  };
  const SPECIES = ['N2', 'O2', 'Ar', 'CO2', 'H2O'];
  const AIR_X = { N2: 0.78084, O2: 0.20946, Ar: 0.00934, CO2: 0.00036 };
  let M_AIR = 0; for (const k in AIR_X) M_AIR += AIR_X[k] * NASA[k][0];
  const FUEL_Y = 23 / 12, M_FUEL = 12.011 + FUEL_Y * 1.00794;
  const hRT = (a, T) => a[0] + a[1] * T / 2 + a[2] * T ** 2 / 3 + a[3] * T ** 3 / 4 + a[4] * T ** 4 / 5 + a[5] / T;
  const sR = (a, T) => a[0] * Math.log(T) + a[1] * T + a[2] * T ** 2 / 2 + a[3] * T ** 3 / 3 + a[4] * T ** 4 / 4 + a[6];
  const cpR = (a, T) => a[0] + a[1] * T + a[2] * T ** 2 + a[3] * T ** 3 + a[4] * T ** 4;
  const COEF = {};
  for (const k of SPECIES) {
    const lo = NASA[k][1].slice(), hi = NASA[k][2].slice();
    hi[5] += (hRT(lo, T_MID) - hRT(hi, T_MID)) * T_MID;
    hi[6] += sR(lo, T_MID) - sR(hi, T_MID);
    COEF[k] = [lo, hi];
  }
  function composition(f) {
    const n = {};
    for (const k of SPECIES) n[k] = (AIR_X[k] || 0) / M_AIR;
    const nC = f / M_FUEL;
    n.CO2 += nC; n.H2O += nC * FUEL_Y / 2; n.O2 -= nC * (1 + FUEL_Y / 4);
    if (n.O2 < 0) throw new Error('fuel-air ratio above stoichiometric');
    for (const k of SPECIES) n[k] /= (1 + f);
    return n;
  }
  function Mixture(f) {
    f = f || 0;
    const n = composition(f), lo = [0, 0, 0, 0, 0, 0, 0], hi = [0, 0, 0, 0, 0, 0, 0];
    let nt = 0;
    for (const k of SPECIES) {
      for (let j = 0; j < 7; j++) { lo[j] += n[k] * RU * COEF[k][0][j]; hi[j] += n[k] * RU * COEF[k][1][j]; }
    }
    for (const k of SPECIES) nt += n[k];
    const R = RU * nt, href = hRT(lo, T_REF) * T_REF;
    const A = (T) => (T <= T_MID ? lo : hi);
    const g = {
      kind: 'varcp', f, R,
      h: (T) => hRT(A(T), T) * T - href,
      cp: (T) => cpR(A(T), T),
      s0: (T) => sR(A(T), T),
    };
    g.gamma = (T) => { const c = g.cp(T); return c / (c - R); };
    g.T_from_h = (h) => {
      let T = Math.min(Math.max(T_REF + h / 1150, 150), 4000);
      for (let k = 0; k < 60; k++) { const dT = (g.h(T) - h) / g.cp(T); T -= dT; if (Math.abs(dT) < 1e-12 * T) break; }
      return T;
    };
    g.T_from_s0 = (s, guess) => {
      let T = guess == null ? 800 : guess;
      for (let k = 0; k < 60; k++) { const d = (g.s0(T) - s) / g.cp(T); T *= Math.exp(-d); if (Math.abs(d) < 1e-13) break; }
      return T;
    };
    g.choke_T = (Tt) => {
      const ht = g.h(Tt);
      let T = 2 * Tt / (g.gamma(Tt) + 1);
      for (let k = 0; k < 60; k++) {
        const gm = g.gamma(T), F = 2 * (ht - g.h(T)) - gm * R * T, dT = F / (2 * g.cp(T) + gm * R);
        T += dT; if (Math.abs(dT) < 1e-12 * T) break;
      }
      return T;
    };
    return g;
  }
  function CPG(gamma, cp) {
    const R = (gamma - 1) / gamma * cp;
    return { kind: 'cpg', g0: gamma, cp0: cp, R, h: (T) => cp * T, cp: () => cp, gamma: () => gamma, s0: (T) => cp * Math.log(T),
             T_from_h: (h) => h / cp, T_from_s0: (s) => Math.exp(s / cp), choke_T: (Tt) => 2 * Tt / (gamma + 1) };
  }
  function makeThermo(i) {
    if (i.thermo === 'cpg') { const c = CPG(i.gamma_c, i.cp_c), t = CPG(i.gamma_t, i.cp_t); return { air: c, gas: () => t }; }
    const air = Mixture(0);
    return { air, gas: (f) => (f > 0 ? Mixture(f) : air) };
  }
  const T_isen = (g, T1, PR) => g.kind === 'varcp' ? g.T_from_s0(g.s0(T1) + g.R * Math.log(PR)) : T1 * PR ** (g.R / g.cp0);
  const pr_isen = (g, T1, T2) => g.kind === 'cpg' ? (T2 / T1) ** (g.cp0 / g.R) : Math.exp((g.s0(T2) - g.s0(T1)) / g.R);
  const compressPoly = (g, T1, pi, e) => g.kind === 'cpg' ? T1 * pi ** (g.R / (g.cp0 * e))
    : g.T_from_s0(g.s0(T1) + g.R * Math.log(pi) / e, T1 * pi ** 0.3);
  const turbPiPoly = (g, T1, T2, e) => g.kind === 'cpg' ? (T2 / T1) ** (g.cp0 / (g.R * e)) : Math.exp((g.s0(T2) - g.s0(T1)) / (e * g.R));
  const compEta = (g, T1, T2, pi) => T2 === T1 ? 1 : (g.h(T_isen(g, T1, pi)) - g.h(T1)) / (g.h(T2) - g.h(T1));
  const turbEta = (g, T1, T2, pi) => T2 === T1 ? 1 : (g.h(T1) - g.h(T2)) / (g.h(T1) - g.h(T_isen(g, T1, pi)));
  const compPiEta = (g, T1, dh, eta) => pr_isen(g, T1, g.T_from_h(g.h(T1) + eta * dh));
  const turbPiEta = (g, T1, dh, eta) => pr_isen(g, T1, g.T_from_h(g.h(T1) - dh / eta));
  function machState(g, Tt, T) { const V = Math.sqrt(Math.max(0, 2 * (g.h(Tt) - g.h(T)))); return [V / Math.sqrt(g.gamma(T) * g.R * T), V]; }
  function nozzle(g, Tt, Pt, P0, kind) {
    if (Pt <= P0 * (1 + 1e-12)) return null;
    if (kind === 'convergent') {
      const Tc = g.choke_T(Tt), Pc = Pt * pr_isen(g, Tt, Tc);
      if (Pc > P0) return { P: Pc, T: Tc, M: 1, V: Math.sqrt(2 * (g.h(Tt) - g.h(Tc))) };
    }
    const T = T_isen(g, Tt, P0 / Pt), mv = machState(g, Tt, T);
    return { P: P0, T, M: mv[0], V: mv[1] };
  }
  function massFlux(g, Tt, Pt, P) {
    const Tc = g.choke_T(Tt), Pc = Pt * pr_isen(g, Tt, Tc);
    let T;
    if (P <= Pc) { T = Tc; P = Pc; } else T = T_isen(g, Tt, P / Pt);
    return P / (g.R * T) * Math.sqrt(Math.max(0, 2 * (g.h(Tt) - g.h(T))));
  }
  const vFE = (g, Tt, Pt, P0) => { const T = T_isen(g, Tt, P0 / Pt); return Math.sqrt(Math.max(0, 2 * (g.h(Tt) - g.h(T)))); };

  const DEFAULTS = {
    M0: 0.78, alt: 10668, T0: null, P0: null, alpha: 5.3, pi_f: 1.65, pi_cL: 2.7, pi_cH: 11.0,
    Tt4: 1500, eps_cool: 0.12, thermo: 'varcp', gamma_c: 1.4, cp_c: 1004, gamma_t: 1.33, cp_t: 1156, hPR: 42.8e6,
    pi_d_max: 0.995, pi_b: 0.96, eta_b: 0.995, e_f: 0.89, e_cL: 0.89, e_cH: 0.90, e_tH: 0.89,
    e_tL: 0.90, eta_mH: 0.99, eta_mL: 0.99, pi_n: 0.99, pi_fn: 0.99,
    core_nozzle: 'convergent', fan_nozzle: 'convergent', neglect_fuel_mass: false, mil_spec_recovery: true,
  };
  function ramRecovery(M0) {
    if (M0 <= 1) return 1;
    if (M0 < 5) return 1 - 0.075 * Math.pow(M0 - 1, 1.35);
    return 800 / (Math.pow(M0, 4) + 935);
  }
  function freestream(i, air) {
    const atm = atmosphere(i.alt);
    const T0 = i.T0 != null ? i.T0 : atm.T, P0 = i.P0 != null ? i.P0 : atm.P;
    const a0 = Math.sqrt(air.gamma(T0) * air.R * T0), V0 = a0 * i.M0;
    const Tt0 = air.kind === 'varcp' ? air.T_from_h(air.h(T0) + 0.5 * V0 * V0) : T0 * (1 + 0.5 * (air.g0 - 1) * i.M0 ** 2);
    return { T0, P0, a0, V0, Tt0, Pt0: P0 * pr_isen(air, T0, Tt0) };
  }
  function burnerFar(air, hot, Tt3, Tt4, eta_b, hPR, nf) {
    const h3 = air.h(Tt3);
    if (nf) return (hot(0.02).h(Tt4) - h3) / (eta_b * hPR);
    let fb = 0.02;
    for (let k = 0; k < 100; k++) {
      const h4 = hot(fb).h(Tt4), fn = (h4 - h3) / (eta_b * hPR - h4);
      if (Math.abs(fn - fb) < 1e-15) { fb = fn; break; }
      fb = fn;
    }
    return fb;
  }

  // ---------------- on-design ----------------
  function design(user) {
    const i = Object.assign({}, DEFAULTS, user);
    const th = makeThermo(i), air = th.air;
    const fs = freestream(i, air), T0 = fs.T0, P0 = fs.P0, V0 = fs.V0;
    const pi_d = i.pi_d_max * (i.mil_spec_recovery ? ramRecovery(i.M0) : 1);
    const nf = !!i.neglect_fuel_mass, eps = i.eps_cool;
    const Tt2 = fs.Tt0, Pt2 = fs.Pt0 * pi_d;
    const Tt13 = compressPoly(air, Tt2, i.pi_f, i.e_f), Pt13 = Pt2 * i.pi_f;
    const Tt25 = compressPoly(air, Tt2, i.pi_cL, i.e_cL), Pt25 = Pt2 * i.pi_cL;
    const Tt3 = compressPoly(air, Tt25, i.pi_cH, i.e_cH), Pt3 = Pt25 * i.pi_cH;
    const Tt4 = i.Tt4, Pt4 = Pt3 * i.pi_b;
    const out = { inputs: i, valid: false, T0, P0, a0: fs.a0, V0, pi_d, OPR: i.pi_cL * i.pi_cH, f: NaN };
    if (Tt4 <= Tt3) { out.reason = 'Tt4 below compressor exit temperature'; return out; }
    const fb = burnerFar(air, th.gas, Tt3, Tt4, i.eta_b, i.hPR, nf);
    const f = (1 - eps) * fb, mb = (1 - eps) * (nf ? 1 : 1 + fb), m45 = nf ? 1 : 1 + f;
    const g4 = th.gas(fb), g45 = eps > 0 ? th.gas(f) : g4;
    out.f = f; out.f_b = fb;
    const h4 = g4.h(Tt4), W_HPC = air.h(Tt3) - air.h(Tt25), h44 = h4 - W_HPC / (i.eta_mH * mb);
    const W_LP_load = (air.h(Tt25) - air.h(Tt2)) + i.alpha * (air.h(Tt13) - air.h(Tt2));
    if (!(fb > 0) || h44 <= g4.h(150)) { out.reason = 'HPT cannot drive the HPC'; return out; }
    const Tt44 = g4.T_from_h(h44), pi_tH = turbPiPoly(g4, Tt4, Tt44, i.e_tH), Pt44 = Pt4 * pi_tH;
    let h45 = h44, Tt45 = Tt44;
    if (eps > 0) { h45 = (mb * h44 + eps * air.h(Tt3)) / m45; Tt45 = g45.T_from_h(h45); }
    const Pt45 = Pt44, h5 = h45 - W_LP_load / (i.eta_mL * m45);
    if (h5 <= g45.h(150)) { out.reason = 'LPT cannot drive the fan and LPC'; return out; }
    const Tt5 = g45.T_from_h(h5), pi_tL = turbPiPoly(g45, Tt45, Tt5, i.e_tL), Pt5 = Pt45 * pi_tL;
    const Tt9 = Tt5, Pt9 = Pt5 * i.pi_n, Tt19 = Tt13, Pt19 = Pt13 * i.pi_fn;
    const n9 = nozzle(g45, Tt9, Pt9, P0, i.core_nozzle);
    const n19 = i.alpha > 0 ? nozzle(air, Tt19, Pt19, P0, i.fan_nozzle) : { P: P0, T: T0, M: i.M0, V: V0 };
    if (!n9 || !n19) { out.reason = 'nozzle total pressure below ambient'; return out; }
    const Fc = m45 * n9.V - V0 + m45 * g45.R * n9.T * (1 - P0 / n9.P) / n9.V;
    const Fb = n19.V - V0 + (i.alpha > 0 ? air.R * n19.T * (1 - P0 / n19.P) / n19.V : 0);
    const Fs = (Fc + i.alpha * Fb) / (1 + i.alpha), S = Fs > 0 ? f / ((1 + i.alpha) * Fs) : NaN;
    const Vfe9 = vFE(g45, Tt9, Pt9, P0), Vfe19 = i.alpha > 0 ? vFE(air, Tt19, Pt19, P0) : V0;
    const dKE = 0.5 * (m45 * Vfe9 ** 2 + i.alpha * Vfe19 ** 2 - (1 + i.alpha) * V0 ** 2);
    const q = f * i.hPR, eta_th = dKE / q, eta_o = (1 + i.alpha) * Fs * V0 / q;
    Object.assign(out, {
      valid: Fs > 0, Fs, S, S_mg: S * 1e6, S_lb: S * TSFC_LB, eta_th, eta_o, eta_p: eta_th > 0 ? eta_o / eta_th : NaN,
      V9: n9.V, V19: n19.V, Vfe9, Vfe19, M9: n9.M, M19: n19.M, P9: n9.P, P19: n19.P, T9: n9.T, T19: n19.T,
      core_choked: n9.P > P0 * (1 + 1e-12), fan_choked: n19.P > P0 * (1 + 1e-12),
      tau_tH: Tt44 / Tt4, tau_tL: Tt5 / Tt45, pi_tH, pi_tL, mb, m45, m4: m45,
      eta_ad: { f: compEta(air, Tt2, Tt13, i.pi_f), cL: compEta(air, Tt2, Tt25, i.pi_cL), cH: compEta(air, Tt25, Tt3, i.pi_cH),
                tH: turbEta(g4, Tt4, Tt44, pi_tH), tL: turbEta(g45, Tt45, Tt5, pi_tL) },
      W_HPT: mb * (h4 - h44), W_LPT: m45 * (h45 - h5), W_HPC, W_LP_load,
      dh_f: air.h(Tt13) - air.h(Tt2), dh_cL: air.h(Tt25) - air.h(Tt2),
      stations: {
        '0': { Tt: fs.Tt0, Pt: fs.Pt0 }, '2': { Tt: Tt2, Pt: Pt2 }, '13': { Tt: Tt13, Pt: Pt13 }, '19': { Tt: Tt19, Pt: Pt19 },
        '2.5': { Tt: Tt25, Pt: Pt25 }, '3': { Tt: Tt3, Pt: Pt3 }, '4': { Tt: Tt4, Pt: Pt4 }, '4.4': { Tt: Tt44, Pt: Pt44 },
        '4.5': { Tt: Tt45, Pt: Pt45 }, '5': { Tt: Tt5, Pt: Pt5 }, '9': { Tt: Tt9, Pt: Pt9 },
      },
    });
    return out;
  }

  function idealTurbofan(M0, T0, gamma, cp, hPR, Tt4, pi_c, pi_f, alpha) {
    const R = (gamma - 1) / gamma * cp, a0 = Math.sqrt(gamma * R * T0);
    const tr = 1 + 0.5 * (gamma - 1) * M0 * M0, tl = Tt4 / T0;
    const tc = Math.pow(pi_c, (gamma - 1) / gamma), tf = Math.pow(pi_f, (gamma - 1) / gamma);
    const V9a = Math.sqrt(2 / (gamma - 1) * (tl - tr * (tc - 1 + alpha * (tf - 1)) - tl / (tr * tc)));
    const V19a = Math.sqrt(2 / (gamma - 1) * (tr * tf - 1));
    const Fs = a0 / (1 + alpha) * (V9a - M0 + alpha * (V19a - M0));
    const f = cp * T0 / hPR * (tl - tr * tc), S = f / ((1 + alpha) * Fs);
    const eta_th = 1 - 1 / (tr * tc);
    const eta_p = 2 * M0 * (V9a + alpha * V19a - (1 + alpha) * M0) / (V9a * V9a + alpha * V19a * V19a - (1 + alpha) * M0 * M0);
    const tfo = (tl - tr * (tc - 1) - tl / (tr * tc) + alpha * tr + 1) / (tr * (1 + alpha));
    return { Fs, S, S_lb: S * TSFC_LB, f, eta_th, eta_p, eta_o: eta_th * eta_p, V9: V9a * a0, V19: V19a * a0,
             tau_f_opt: tfo, pi_f_opt: Math.pow(tfo, gamma / (gamma - 1)) };
  }

  function piFMax(inp) {
    let a = 1.01, b = 8;
    if (design(Object.assign({}, inp, { pi_f: b })).valid) return b;
    for (let k = 0; k < 50; k++) {
      const m = 0.5 * (a + b);
      if (design(Object.assign({}, inp, { pi_f: m })).valid) a = m; else b = m;
    }
    return a;
  }
  function golden(fn, a, b, tol) {
    const g = (Math.sqrt(5) - 1) / 2;
    let x1 = b - g * (b - a), x2 = a + g * (b - a), f1 = fn(x1), f2 = fn(x2);
    while (b - a > tol) {
      if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = fn(x1); }
      else { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = fn(x2); }
    }
    return 0.5 * (a + b);
  }
  function optimumFanPR(inp, tol) {
    const hi = piFMax(inp);
    const fn = (pf) => { const r = design(Object.assign({}, inp, { pi_f: pf })); return r.valid && r.S > 0 ? r.S : 1e9; };
    const pf = golden(fn, 1.02, hi, tol || 1e-8);
    return { pi_f: pf, result: design(Object.assign({}, inp, { pi_f: pf })) };
  }

  // ---------------- Brent root finder (port of scipy's brentq) ----------------
  function brentq(f, xa, xb, xtol, rtol, maxiter) {
    xtol = xtol || 1e-14; rtol = rtol || 1e-15; maxiter = maxiter || 300;
    let xpre = xa, xcur = xb, xblk = 0, fpre = f(xpre), fcur = f(xcur), fblk = 0, spre = 0, scur = 0;
    if (fpre === 0) return xpre;
    if (fcur === 0) return xcur;
    if (Math.sign(fpre) === Math.sign(fcur)) throw new Error('brentq: no sign change');
    for (let it = 0; it < maxiter; it++) {
      if (fpre !== 0 && fcur !== 0 && Math.sign(fpre) !== Math.sign(fcur)) { xblk = xpre; fblk = fpre; spre = scur = xcur - xpre; }
      if (Math.abs(fblk) < Math.abs(fcur)) { xpre = xcur; xcur = xblk; xblk = xpre; fpre = fcur; fcur = fblk; fblk = fpre; }
      const delta = (xtol + rtol * Math.abs(xcur)) / 2, sbis = (xblk - xcur) / 2;
      if (fcur === 0 || Math.abs(sbis) < delta) return xcur;
      if (Math.abs(spre) > delta && Math.abs(fcur) < Math.abs(fpre)) {
        let stry;
        if (xpre === xblk) stry = -fcur * (xcur - xpre) / (fcur - fpre);
        else {
          const dpre = (fpre - fcur) / (xpre - xcur), dblk = (fblk - fcur) / (xblk - xcur);
          stry = -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre));
        }
        if (2 * Math.abs(stry) < Math.min(Math.abs(spre), 3 * Math.abs(sbis) - delta)) { spre = scur; scur = stry; }
        else { spre = sbis; scur = sbis; }
      } else { spre = sbis; scur = sbis; }
      xpre = xcur; fpre = fcur;
      xcur += Math.abs(scur) > delta ? scur : (sbis > 0 ? delta : -delta);
      fcur = f(xcur);
    }
    return xcur;
  }


  // ---------------- off-design (port of offdesign.py) ----------------
  function Engine(user, F_design) {
    const inp = Object.assign({}, DEFAULTS, user);
    if (inp.core_nozzle !== 'convergent' || inp.fan_nozzle !== 'convergent') throw new Error('off-design needs convergent nozzles');
    if (inp.neglect_fuel_mass) throw new Error('off-design needs fuel-mass terms');
    const r = design(inp);
    if (!r.valid) throw new Error('design not valid');
    const th = makeThermo(inp), air = th.air, e = r.eta_ad, eps = inp.eps_cool;
    const mdot0 = F_design != null ? F_design / r.Fs : 100;
    const st = r.stations, mc = mdot0 / (1 + inp.alpha);
    const g4d = th.gas(r.f_b), g45d = eps > 0 ? th.gas(r.f) : g4d;
    const A4 = mc * r.mb / massFlux(g4d, st['4'].Tt, st['4'].Pt, 0);
    const A45 = mc * r.m45 / massFlux(g45d, st['4.5'].Tt, st['4.5'].Pt, 0);
    const A8 = mc * r.m45 / massFlux(g45d, st['9'].Tt, st['9'].Pt, r.P0);
    const A18 = inp.alpha * mc / massFlux(air, st['19'].Tt, st['19'].Pt, r.P0);
    const K = r.dh_cL / r.dh_f;
    const ref = { dh_f: r.dh_f, dh_cH: r.W_HPC, Tt2: st['2'].Tt, Tt25: st['2.5'].Tt, Pt2: st['2'].Pt, Tt3: st['3'].Tt };

    function hpSpool(Tt25, Tt4, s) {
      let Tt3 = s.Tt3 != null ? s.Tt3 : ref.Tt3 * Tt25 / ref.Tt25, Tt44 = NaN;
      const h25 = air.h(Tt25);
      for (let it = 0; it < 200; it++) {
        if (Tt3 >= Tt4) return null;
        const fb = burnerFar(air, th.gas, Tt3, Tt4, inp.eta_b, inp.hPR, false);
        const f = (1 - eps) * fb, mb = (1 - eps) * (1 + fb), m45 = 1 + f;
        const g4 = th.gas(fb), g45 = eps > 0 ? th.gas(f) : g4, h4 = g4.h(Tt4), h3 = air.h(Tt3);
        const target = A4 * massFlux(g4, Tt4, 1, 0) * m45 / mb;
        const mix = (T44) => eps === 0 ? T44 : g45.T_from_h((mb * g4.h(T44) + eps * h3) / m45);
        const G = (T44) => A45 * turbPiEta(g4, Tt4, h4 - g4.h(T44), e.tH) * massFlux(g45, mix(T44), 1, 0) - target;
        const lo = 0.35 * Tt4, hi = Tt4 * (1 - 1e-12);
        if (G(hi) < 0 || G(lo) > 0) return null;
        Tt44 = brentq(G, lo, hi, 1e-14, 1e-15, 200);
        const Tt3n = air.T_from_h(h25 + inp.eta_mH * mb * (h4 - g4.h(Tt44)));
        const done = Math.abs(Tt3n - Tt3) < 1e-11 * Tt3;
        Tt3 = Tt3n;
        if (done) break;
      }
      const fb = burnerFar(air, th.gas, Tt3, Tt4, inp.eta_b, inp.hPR, false);
      const f = (1 - eps) * fb, mb = (1 - eps) * (1 + fb), m45 = 1 + f;
      const g4 = th.gas(fb), g45 = eps > 0 ? th.gas(f) : g4;
      const pi_tH = turbPiEta(g4, Tt4, g4.h(Tt4) - g4.h(Tt44), e.tH);
      const Tt45 = eps === 0 ? Tt44 : g45.T_from_h((mb * g4.h(Tt44) + eps * air.h(Tt3)) / m45);
      s.Tt3 = Tt3;
      return { Tt3, fb, f, mb, m45, g4, g45, Tt44, Tt45, pi_tH, W_HPC: air.h(Tt3) - h25 };
    }

    function evaluate(dh_f, s) {
      const h2 = air.h(s.Tt2), dh_cL = K * dh_f, P0 = s.P0, Tt4 = s.Tt4;
      const Tt13 = air.T_from_h(h2 + dh_f), Tt25 = air.T_from_h(h2 + dh_cL);
      const pi_f = compPiEta(air, s.Tt2, dh_f, e.f), pi_cL = compPiEta(air, s.Tt2, dh_cL, e.cL);
      const hp = hpSpool(Tt25, Tt4, s);
      if (!hp) return null;
      const pi_cH = compPiEta(air, Tt25, hp.W_HPC, e.cH);
      const Pt4 = s.Pt2 * pi_cL * pi_cH * inp.pi_b;
      const g4 = hp.g4, g45 = hp.g45, m45 = hp.m45;
      const mdot4 = A4 * massFlux(g4, Tt4, Pt4, 0), mcore = mdot4 / hp.mb;
      const Pt45 = Pt4 * hp.pi_tH, Tt45 = hp.Tt45, h45 = g45.h(Tt45), mflow = mcore * m45;
      const g = (T5) => {
        const Pt9 = Pt45 * turbPiEta(g45, Tt45, h45 - g45.h(T5), e.tL) * inp.pi_n;
        return Pt9 <= P0 ? -mflow : A8 * massFlux(g45, T5, Pt9, P0) - mflow;
      };
      let lo = g45.T_from_h(h45 - e.tL * (h45 - g45.h(0.3 * Tt45)) * 0.999);
      lo = Math.max(lo, 0.3 * Tt45);
      const hi = Tt45 * (1 - 1e-12);
      if (g(hi) < 0 || g(lo) > 0) return null;
      const Tt5 = brentq(g, lo, hi, 1e-14, 1e-15, 200);
      const pi_tL = turbPiEta(g45, Tt45, h45 - g45.h(Tt5), e.tL);
      const Pt13 = s.Pt2 * pi_f, Pt19 = Pt13 * inp.pi_fn;
      const mbyp = Pt19 > P0 ? A18 * massFlux(air, Tt13, Pt19, P0) : 0;
      const alpha = mbyp / mcore;
      const supply = inp.eta_mL * m45 * (h45 - g45.h(Tt5)), demand = dh_cL + alpha * dh_f;
      return { dh_f, pi_f, pi_cL, pi_cH, f: hp.f, m45, alpha, mcore, g45, Tt25, Tt45, Pt45, Tt5, Pt5: Pt45 * pi_tL, Tt13, Pt13,
               W_HPC: hp.W_HPC, res: supply - demand, demand };
    }

    function operate(M0, alt, Tt4) {
      const fs = freestream(Object.assign({}, inp, { M0, alt, T0: null, P0: null }), air), P0 = fs.P0, V0 = fs.V0;
      const pi_d = inp.pi_d_max * (inp.mil_spec_recovery ? ramRecovery(M0) : 1);
      const s = { Tt2: fs.Tt0, Pt2: fs.Pt0 * pi_d, Tt4, P0 };
      const bad = { valid: false, M0, alt, Tt4 };
      const R = (x) => { const ev = evaluate(x * ref.dh_f, s); return ev === null ? -1e9 : ev.res / ref.dh_f; };
      const lo = 1e-7;
      if (R(lo) <= 0) return bad;
      let hi = 1.5, k = 0;
      while (R(hi) > 0) { hi *= 1.5; if (++k > 40) return bad; }
      const x = brentq(R, lo, hi, 1e-14, 1e-15, 300);
      const ev = evaluate(x * ref.dh_f, s);
      if (!ev) return bad;
      const g45 = ev.g45, m45 = ev.m45, alpha = ev.alpha, f = ev.f;
      const Pt9 = ev.Pt5 * inp.pi_n, Pt19 = ev.Pt13 * inp.pi_fn;
      const n9 = nozzle(g45, ev.Tt5, Pt9, P0, 'convergent'), n19 = nozzle(air, ev.Tt13, Pt19, P0, 'convergent');
      if (!n9) return bad;
      const Fc = m45 * n9.V - V0 + m45 * g45.R * n9.T * (1 - P0 / n9.P) / n9.V;
      const Fb = n19 ? n19.V - V0 + air.R * n19.T * (1 - P0 / n19.P) / n19.V : 0;
      const Fs = (Fc + alpha * Fb) / (1 + alpha), mdot0 = ev.mcore * (1 + alpha), F = mdot0 * Fs;
      const S = Fs > 0 ? f / ((1 + alpha) * Fs) : NaN;
      const Vfe9 = vFE(g45, ev.Tt5, Pt9, P0), Vfe19 = n19 ? vFE(air, ev.Tt13, Pt19, P0) : 0;
      const dKE = 0.5 * (m45 * Vfe9 ** 2 + alpha * Vfe19 ** 2 - (1 + alpha) * V0 ** 2);
      const eta_th = dKE / (f * inp.hPR), eta_o = (1 + alpha) * Fs * V0 / (f * inp.hPR);
      return { valid: Fs > 0, M0, alt, T0: fs.T0, P0, Tt4, F, Fs, mdot0, S, S_lb: S * TSFC_LB, f, alpha,
               pi_f: ev.pi_f, OPR: ev.pi_cL * ev.pi_cH, eta_th, eta_o, eta_p: eta_o / eta_th,
               V9: n9.V, V19: n19 ? n19.V : 0, core_choked: n9.P > P0 * (1 + 1e-12), fan_choked: !!n19 && n19.P > P0 * (1 + 1e-12),
               NL_corr: Math.sqrt(ev.dh_f / s.Tt2 / (ref.dh_f / ref.Tt2)), lp_residual: ev.res / ev.demand };
    }
    return { inp, des: r, F_des: mdot0 * r.Fs, mdot0_des: mdot0, operate };
  }

  function schedule(inp, M0, alt, dTmax) {
    const thd = atmosphere(inp.alt).T * (1 + 0.2 * inp.M0 * inp.M0) / 288.15;
    const th = atmosphere(alt).T * (1 + 0.2 * M0 * M0) / 288.15;
    return Math.min(inp.Tt4 + dTmax, inp.Tt4 * th / thd);
  }


  // ---------------- T-s path (port of ondesign.ts_path) ----------------
  // core indices: 0:0 static, 1:0, 2:2, 14:2.5, 26:3, 38:4, 50:4.4, 51:4.5, 63:5, 64:9 static
  function tsPath(r, n) {
    n = n || 12;
    const i = r.inputs, th = makeThermo(i), air = th.air;
    const g4 = th.gas(r.f_b), g45 = i.eps_cool > 0 ? th.gas(r.f) : g4;
    const st = r.stations, T0 = r.T0, P0 = r.P0;
    function leg(pts, s1, T1, P1, T2, P2, g, curved) {
      const s2 = s1 + g.s0(T2) - g.s0(T1) - g.R * Math.log(P2 / P1), k = curved === false ? 1 : n, d = g.s0(T2) - g.s0(T1);
      for (let j = 1; j <= k; j++) {
        const T = T1 * (T2 / T1) ** (j / k), x = d !== 0 ? (g.s0(T) - g.s0(T1)) / d : j / k;
        pts.push([s1 + x * (s2 - s1), T]);
      }
      return s2;
    }
    const S = (k) => [st[k].Tt, st[k].Pt];
    const core = [[0, T0]];
    let s = leg(core, 0, T0, P0, ...S('0'), air, false);
    s = leg(core, s, ...S('0'), ...S('2'), air, false);
    const s2 = s;
    s = leg(core, s, ...S('2'), ...S('2.5'), air);
    s = leg(core, s, ...S('2.5'), ...S('3'), air);
    s = leg(core, s, ...S('3'), ...S('4'), g4);
    s = leg(core, s, ...S('4'), ...S('4.4'), g4);
    s = leg(core, s, ...S('4.4'), ...S('4.5'), g45, false);
    s = leg(core, s, ...S('4.5'), ...S('5'), g45);
    s = leg(core, s, ...S('5'), r.T9, r.P9, g45, false);
    const byp = [[s2, st['2'].Tt]];
    let sb = leg(byp, s2, ...S('2'), ...S('13'), air);
    sb = leg(byp, sb, ...S('13'), r.T19, r.P19, air, false);
    const smax = Math.max(s, sb) * 1.05 + 50, amb = [];
    for (let j = 0; j <= 20; j++) {
      const x = smax * j / 20;
      amb.push([x, air.kind === 'varcp' ? air.T_from_s0(air.s0(T0) + x) : T0 * Math.exp(x / air.cp0)]);
    }
    return { core, bypass: byp, ambient: amb };
  }

  const TF = { FT, TSFC_LB, DEFAULTS, atmosphere, Mixture, CPG, makeThermo, design, idealTurbofan, optimumFanPR, piFMax, golden,
               brentq, Engine, schedule, tsPath };
  if (typeof module !== 'undefined' && module.exports) module.exports = TF;
  else root.TF = TF;
})(typeof self !== 'undefined' ? self : this);
