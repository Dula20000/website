/* eqsolver.js — browser/Node port of the Python `eqsolver` package.
 *
 * Gibbs free-energy minimisation with the NASA RP-1311 reduced Newton
 * iteration (element potentials π_i, Δln n, Δln T), HP / SP / TP problems,
 * equilibrium derivatives (γ_s, cp_eq), and ideal-rocket performance
 * (infinite-area combustor, throat by u = a, exit by area ratio) for frozen
 * and shifting expansion.  Same algorithms, same step control and the same
 * thermo/propellant data (loaded from data/thermo.json, which run_study.py
 * writes from the Python tables) — so the two should agree to round-off.
 *
 * Usage:  const S = EqSolver.create(thermoJson);
 *         S.performance('LCH4', 'LOX', 3.6, 30e6, 40, {frozen:false})
 */
(function (root) {
  'use strict';
  const R = 8.314462618, PREF = 1e5, G0 = 9.80665, PATM = 101325;
  const LN_TRACE = -18.420681, LN_FLOOR = -80.0;

  function solveLinear(A, b) {            // Gaussian elimination, partial pivoting
    const n = b.length, M = A.map((r, i) => r.concat([b[i]]));
    for (let k = 0; k < n; k++) {
      let p = k, mx = Math.abs(M[k][k]);
      for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > mx) { mx = Math.abs(M[i][k]); p = i; }
      if (mx === 0) throw new Error('singular matrix');
      if (p !== k) { const t = M[k]; M[k] = M[p]; M[p] = t; }
      for (let i = k + 1; i < n; i++) {
        const f = M[i][k] / M[k][k];
        if (f === 0) continue;
        for (let j = k; j <= n; j++) M[i][j] -= f * M[k][j];
      }
    }
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      let s = M[i][n];
      for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
      x[i] = s / M[i][i];
    }
    return x;
  }
  const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
  const sum = (a) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s; };

  function create(data) {
    const SP = data.species, AM = data.atomic_mass;

    function GasSystem(names, elements) {
      const sp = names.filter((n) => Object.keys(SP[n].comp).every((e) => elements.includes(e)));
      const A = elements.map((e) => sp.map((n) => SP[n].comp[e] || 0));
      const low = sp.map((n) => SP[n].low), high = sp.map((n) => SP[n].high), tmid = sp.map((n) => SP[n].tmid);
      function thermo(T) {
        const ns = sp.length, cp = new Float64Array(ns), h = new Float64Array(ns), s = new Float64Array(ns);
        const lT = Math.log(T);
        for (let j = 0; j < ns; j++) {
          const a = T >= tmid[j] ? high[j] : low[j];
          cp[j] = a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4])));
          h[j] = a[0] + T * (a[1] / 2 + T * (a[2] / 3 + T * (a[3] / 4 + T * a[4] / 5))) + a[5] / T;
          s[j] = a[0] * lT + T * (a[1] + T * (a[2] / 2 + T * (a[3] / 3 + T * a[4] / 4))) + a[6];
        }
        return { cp, h, s };
      }
      return { names: sp, elements, A, thermo };
    }

    function defaultSystem(elements) {
      let names = data.cho_species.slice();
      if (elements.includes('N')) names = names.concat(['N2', 'NO']);
      if (elements.includes('Ar')) names.push('Ar');
      return GasSystem(names, elements);
    }

    function makeState(sys, T, P, nj, extra) {
      const st = Object.assign({ sys, T, P, nj, n: sum(nj), frozen: false }, extra || {});
      return st;
    }
    function stH(st) { const t = st.sys.thermo(st.T); return dot(st.nj, t.h) * R * st.T; }
    function stS(st) {
      const t = st.sys.thermo(st.T), ntot = sum(st.nj), lp = Math.log(st.P / PREF);
      let s = 0;
      for (let j = 0; j < st.nj.length; j++) if (st.nj[j] > 0) s += st.nj[j] * (t.s[j] - Math.log(st.nj[j] / ntot) - lp);
      return s * R;
    }
    function stRho(st) { return st.P / (sum(st.nj) * R * st.T); }
    function stGamma(st) { return st.frozen ? st.gamma_fr : st.gamma_s; }
    function stA(st) { return Math.sqrt(sum(st.nj) * R * st.T * stGamma(st)); }

    function equilibrate(sys, b0, P, opt) {
      opt = opt || {};
      const mode = opt.T != null ? 'TP' : (opt.h0 != null ? 'HP' : 'SP');
      const tol = opt.tol || 1e-11, maxIter = opt.maxIter || 400;
      const A = sys.A, nel = A.length, ns = sys.names.length;
      let T = opt.T != null ? opt.T : 3800.0, nj, n;
      if (opt.init) {
        nj = Float64Array.from(opt.init.nj); n = sum(nj);
        if (mode !== 'TP') T = opt.init.T;
      } else {
        n = 0.1 * Math.max(1.0, sum(b0));
        nj = new Float64Array(ns).fill(n / ns);
      }
      for (let j = 0; j < ns; j++) nj[j] = Math.max(nj[j], n * Math.exp(LN_FLOOR));
      const lnP = Math.log(P / PREF);
      const size = nel + 1 + (mode === 'TP' ? 0 : 1), iT = nel + 1;
      const history = [];
      let it, err = Infinity;
      const mu = new Float64Array(ns), lnx = new Float64Array(ns), dl = new Float64Array(ns), sj = new Float64Array(ns);
      for (it = 1; it <= maxIter; it++) {
        const th = sys.thermo(T), cpR = th.cp, hRT = th.h, sR = th.s;
        const nsum = sum(nj);
        for (let j = 0; j < ns; j++) { lnx[j] = Math.log(nj[j] / n); mu[j] = hRT[j] - sR[j] + lnx[j] + lnP; }
        const G = []; for (let i = 0; i < size; i++) G.push(new Array(size).fill(0));
        const rhs = new Array(size).fill(0);
        const An = A.map((row) => row.map((a, j) => a * nj[j]));
        for (let k = 0; k < nel; k++) {
          for (let i = 0; i < nel; i++) G[k][i] = dot(An[k], A[i]);
          const s = sum(An[k]);
          G[k][nel] = s; G[nel][k] = s;
          let bk = 0; for (let j = 0; j < ns; j++) bk += A[k][j] * nj[j];
          rhs[k] = b0[k] - bk + dot(An[k], mu);
        }
        G[nel][nel] = nsum - n;
        rhs[nel] = n - nsum + dot(nj, mu);
        if (mode === 'HP') {
          const njh = dot(nj, hRT);
          for (let k = 0; k < nel; k++) { const v = dot(An[k], hRT); G[k][iT] = v; G[iT][k] = v; }
          G[nel][iT] = njh; G[iT][nel] = njh;
          let a = 0, b = 0, c = 0;
          for (let j = 0; j < ns; j++) { a += nj[j] * cpR[j]; b += nj[j] * hRT[j] * hRT[j]; c += nj[j] * hRT[j] * mu[j]; }
          G[iT][iT] = a + b;
          rhs[iT] = opt.h0 / (R * T) - njh + c;
        } else if (mode === 'SP') {
          for (let j = 0; j < ns; j++) sj[j] = sR[j] - lnx[j] - lnP;
          for (let k = 0; k < nel; k++) { G[k][iT] = dot(An[k], hRT); G[iT][k] = dot(An[k], sj); }
          G[nel][iT] = dot(nj, hRT); G[iT][nel] = dot(nj, sj);
          let a = 0, b = 0, c = 0;
          for (let j = 0; j < ns; j++) { a += nj[j] * cpR[j]; b += nj[j] * hRT[j] * sj[j]; c += nj[j] * sj[j] * mu[j]; }
          G[iT][iT] = a + b;
          rhs[iT] = opt.s0 / R - dot(nj, sj) + n - nsum + c;
        }
        const sol = solveLinear(G, rhs);
        const dlnn = sol[nel], dlnT = mode === 'TP' ? 0 : sol[iT];
        for (let j = 0; j < ns; j++) {
          let s = 0; for (let i = 0; i < nel; i++) s += A[i][j] * sol[i];
          dl[j] = -mu[j] + s + dlnn + hRT[j] * dlnT;
        }
        // RP-1311 step control
        let m = Math.max(5 * Math.abs(dlnT), 5 * Math.abs(dlnn));
        let lam2 = 1.0;
        for (let j = 0; j < ns; j++) {
          const big = lnx[j] > LN_TRACE;
          if (big && dl[j] > 0) m = Math.max(m, Math.abs(dl[j]));
          if (!big && dl[j] >= 0) {
            const den = dl[j] - dlnn;
            if (den > 0) lam2 = Math.min(lam2, Math.abs((-lnx[j] - 9.2103404) / den));
          }
        }
        // (Python takes min over candidates only when any exist; lam2 starts at 1 in both)
        const lam1 = m > 0 ? 2.0 / m : 1.0;
        const lam = Math.min(1.0, lam1, lam2);
        let e1 = 0; for (let j = 0; j < ns; j++) e1 = Math.max(e1, nj[j] * Math.abs(dl[j]));
        err = Math.max(e1 / nsum, n * Math.abs(dlnn) / nsum, Math.abs(dlnT));
        history.push(err);
        n = Math.exp(Math.log(n) + lam * dlnn);
        const floor = Math.log(n) + LN_FLOOR;
        for (let j = 0; j < ns; j++) nj[j] = Math.exp(Math.max(Math.log(nj[j]) + lam * dl[j], floor));
        if (mode !== 'TP') T = Math.min(Math.max(T * Math.exp(lam * dlnT), 150.0), 7000.0);
        if (err < tol && lam === 1.0) break;
      }
      if (it > maxIter) throw new Error(mode + ' equilibrium did not converge (err=' + err.toExponential(3) + ')');
      const st = makeState(sys, T, P, nj, { iterations: it, history });
      derivatives(st);
      return st;
    }

    function derivatives(st) {
      const A = st.sys.A, nel = A.length, nj = st.nj, th = st.sys.thermo(st.T);
      const An = A.map((row) => row.map((a, j) => a * nj[j]));
      const G = []; for (let i = 0; i <= nel; i++) G.push(new Array(nel + 1).fill(0));
      for (let k = 0; k < nel; k++) {
        for (let i = 0; i < nel; i++) G[k][i] = dot(An[k], A[i]);
        const s = sum(An[k]); G[k][nel] = s; G[nel][k] = s;
      }
      const ntot = sum(nj);
      const rT = An.map((r) => -dot(r, th.h)).concat([-dot(nj, th.h)]);
      const rP = An.map((r) => sum(r)).concat([ntot]);
      const xT = solveLinear(G, rT), xP = solveLinear(G, rP);
      st.dlnV_dlnT = 1 + xT[nel];
      st.dlnV_dlnP = -1 + xP[nel];
      let cp = xT[nel] * dot(nj, th.h) + dot(nj, th.cp);
      for (let k = 0; k < nel; k++) cp += dot(An[k], th.h) * xT[k];
      for (let j = 0; j < nj.length; j++) cp += nj[j] * th.h[j] * th.h[j];
      st.cp_eq = cp * R;
      st.cp_fr = dot(nj, th.cp) * R;
      const cv = st.cp_eq + ntot * R * st.dlnV_dlnT * st.dlnV_dlnT / st.dlnV_dlnP;
      st.gamma_s = -(st.cp_eq / cv) / st.dlnV_dlnP;
      st.gamma_fr = st.cp_fr / (st.cp_fr - ntot * R);
      return st;
    }

    function frozenState(ref, P, s0) {
      const nj = Float64Array.from(ref.nj), ntot = sum(nj), lnP = Math.log(P / PREF);
      const lnx = Array.from(nj, (v) => Math.log(v / ntot));
      let T = ref.T, k;
      for (k = 0; k < 100; k++) {
        const th = ref.sys.thermo(T);
        let s = 0, cp = 0;
        for (let j = 0; j < nj.length; j++) { s += nj[j] * (th.s[j] - lnx[j] - lnP); cp += nj[j] * th.cp[j]; }
        const d = (s0 - s * R) / (cp * R);
        T *= Math.exp(d);
        if (Math.abs(d) < 1e-12) break;
      }
      if (k === 100) throw new Error('frozen isentrope did not converge');
      const st = makeState(ref.sys, T, P, nj, { frozen: true });
      const th = ref.sys.thermo(T);
      st.cp_fr = st.cp_eq = dot(nj, th.cp) * R;
      st.gamma_fr = st.gamma_s = st.cp_fr / (st.cp_fr - ntot * R);
      return st;
    }

    // ---- propellants -----------------------------------------------------
    const PROP = data.propellants;
    const molarMass = (comp) => Object.entries(comp).reduce((s, [e, k]) => s + AM[e] * k, 0) * 1e-3;
    function mix(list) {        // [[reactant, massFraction], ...]
      const elements = Array.from(new Set(list.flatMap(([r]) => Object.keys(r.comp)))).sort();
      const b0 = new Array(elements.length).fill(0);
      let h0 = 0;
      for (const [r, w] of list) {
        const nm = w / molarMass(r.comp);
        h0 += nm * r.h;
        for (const [e, k] of Object.entries(r.comp)) b0[elements.indexOf(e)] += nm * k;
      }
      return { elements, b0, h0 };
    }
    function bipropellant(fuel, ox, of) {
      return mix([[PROP[fuel], 1 / (1 + of)], [PROP[ox], of / (1 + of)]]);
    }
    function bulkDensity(fuel, ox, of) {
      return (1 + of) / (of / PROP[ox].density + 1 / PROP[fuel].density);
    }

    function chamber(fuel, ox, of, Pc, init) {
      const m = bipropellant(fuel, ox, of);
      const sys = defaultSystem(m.elements);
      let ch;
      try { ch = equilibrate(sys, m.b0, Pc, { h0: m.h0, init: init || null }); }
      catch (e) { if (!init) throw e; ch = equilibrate(sys, m.b0, Pc, { h0: m.h0 }); }
      ch.b0 = m.b0; ch.h0 = m.h0; ch.s0 = stS(ch);
      return ch;
    }

    function expand(ch, P, frozen, prev) {
      if (frozen) return frozenState(ch, P, ch.s0);
      return equilibrate(ch.sys, ch.b0, P, { s0: ch.s0, init: prev || ch });
    }

    function throat(ch, frozen) {
      const g = frozen ? ch.gamma_fr : ch.gamma_s;
      let x = g / (g - 1) * Math.log((g + 1) / 2), prev = null, st, u2, k;
      for (k = 0; k < 60; k++) {
        st = expand(ch, ch.P * Math.exp(-x), frozen, prev); prev = st;
        u2 = 2 * (ch.h0 - stH(st));
        const a2 = Math.pow(stA(st), 2), f = u2 - a2, gm = stGamma(st);
        x += -f * gm / ((gm + 1) * a2);
        if (Math.abs(f) / a2 < 1e-11) break;
      }
      if (k === 60) throw new Error('throat iteration failed');
      return { st, u: Math.sqrt(Math.max(u2, 0)) };
    }

    function exitState(ch, eps, frozen, th, mflux) {
      const xt = Math.log(ch.P / th.P), g = stGamma(th);
      let M = 2 + Math.log(eps), t = 1;
      for (let k = 0; k < 50; k++) {
        t = 1 + 0.5 * (g - 1) * M * M;
        const e = (1 / M) * Math.pow((2 / (g + 1)) * t, (g + 1) / (2 * (g - 1)));
        const dM = (Math.log(e) - Math.log(eps)) / (-1 / M + M / t);
        M = Math.max(M - dM, 1.0001);
        if (Math.abs(dM) < 1e-12) break;
      }
      let x = xt + g / (g - 1) * Math.log(t / (1 + 0.5 * (g - 1)));
      x = Math.max(x, xt + 1e-3);
      let lo = xt, hi = null, prev = th, st, u, k;
      for (k = 0; k < 80; k++) {
        st = expand(ch, ch.P * Math.exp(-x), frozen, prev); prev = st;
        const u2 = 2 * (ch.h0 - stH(st));
        u = Math.sqrt(u2);
        const e = mflux / (stRho(st) * u), f = Math.log(e) - Math.log(eps);
        if (Math.abs(f) < 1e-11) break;
        if (f > 0) hi = x; else lo = x;
        const a2 = Math.pow(stA(st), 2), slope = (1 - a2 / u2) / stGamma(st);
        let xn = slope > 1e-6 ? x - f / slope : NaN;
        if (!(lo < xn && (hi === null || xn < hi)) || !isFinite(xn)) xn = hi !== null ? 0.5 * (lo + hi) : x + 1;
        x = xn;
      }
      if (k === 80) throw new Error('exit area-ratio iteration failed');
      return { st, u };
    }

    function nozzle(ch, eps, opt) {
      opt = opt || {};
      const frozen = !!opt.frozen, Pa = opt.Pa != null ? opt.Pa : PATM;
      const t = throat(ch, frozen);
      const mflux = stRho(t.st) * t.u, cstar = ch.P / mflux;
      const ex = exitState(ch, eps, frozen, t.st, mflux);
      const ivac = ex.u + ex.st.P * eps / mflux, isl = ivac - Pa * eps / mflux;
      return {
        mode: frozen ? 'frozen' : 'shifting', Pc_MPa: ch.P / 1e6, eps,
        Tc: ch.T, M_c: 1e3 / sum(ch.nj), gamma_s_c: ch.gamma_s, gamma_fr_c: ch.gamma_fr, cp_eq_c: ch.cp_eq,
        Tt: t.st.T, Pt_MPa: t.st.P / 1e6, Te: ex.st.T, Pe_kPa: ex.st.P / 1e3, M_e: 1e3 / sum(ex.st.nj),
        gamma_e: stGamma(ex.st),
        cstar, CF_vac: ivac / cstar, CF_sl: isl / cstar, Isp_vac: ivac / G0, Isp_sl: isl / G0, u_e: ex.u,
        chamber: ch, exit: ex.st,
      };
    }

    function performance(fuel, ox, of, Pc, eps, opt) {
      opt = opt || {};
      const ch = opt.chamber || chamber(fuel, ox, of, Pc, opt.init);
      const r = nozzle(ch, eps, opt);
      r.of = of;
      return r;
    }

    function moleFractions(st, thr) {
      const ntot = sum(st.nj), out = {};
      st.sys.names.forEach((n, j) => { const x = st.nj[j] / ntot; if (x > (thr || 0)) out[n] = x; });
      return out;
    }

    return { R, G0, PATM, equilibrate, defaultSystem, chamber, nozzle, performance, bipropellant,
             bulkDensity, moleFractions, h: stH, s: stS, frozenState };
  }

  const api = { create };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.EqSolver = api;
})(typeof self !== 'undefined' ? self : this);
