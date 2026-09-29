/* nozzle-cfd page: fills numbers/tables/figures from data/*.json and drives the live solver. */
(function () {
  'use strict';
  const $ = (id) => document.getElementById(id);
  const ENG_COL = { merlin: 0, rs25: 1, raptor: 2 };

  // ------------------------------------------------------------ formatting
  const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  function sci(v, sig = 2) {
    if (v === null || v === undefined || !isFinite(v)) return '–';
    if (v === 0) return '0';
    const a = Math.abs(v);
    if (a >= 1e-3 && a < 1e5) return Number(v.toPrecision(sig)).toString();
    let e = Math.floor(Math.log10(a)), m = v / Math.pow(10, e);
    if (Math.abs(Number(m.toFixed(Math.max(sig - 1, 0)))) >= 10) { e += 1; m /= 10; }
    return m.toFixed(Math.max(sig - 1, 0)) + '×10' + String(e).split('').map((c) => SUP[c]).join('');
  }
  function fx(v, d) { return v === null || v === undefined || !isFinite(v) ? '–' : Number(v).toFixed(d); }
  function get(root, path) {
    return path.split('.').reduce((o, k) => (o == null ? undefined : o[k]), root);
  }
  function fillSpans(root) {
    document.querySelectorAll('[data-v]').forEach((el) => {
      let v = get(root, el.dataset.v);
      if (v === undefined) { console.warn('missing', el.dataset.v); return; }
      if (typeof v === 'number' && el.dataset.x) v *= Number(el.dataset.x);
      el.textContent = el.dataset.s ? sci(v, Number(el.dataset.s)) : el.dataset.d ? fx(v, Number(el.dataset.d)) : String(v);
    });
  }
  function fillTable(id, rows) {
    const tb = $(id).querySelector('tbody');
    tb.innerHTML = rows.map((r) => '<tr>' + r.map((c, i) => {
      const cls = typeof c === 'object' && c !== null ? c.cls : (i === 0 ? '' : 'num');
      const txt = typeof c === 'object' && c !== null ? c.t : c;
      return `<td class="${cls}">${txt}</td>`;
    }).join('') + '</tr>').join('');
  }
  const flag = (risk) => (risk ? { t: 'risk', cls: 'c fail' } : { t: 'clear', cls: 'c pass' });

  // ------------------------------------------------------------ static figures
  function hline(y, color, dash = 'dash', xref = 'paper') {
    return { type: 'line', xref, x0: 0, x1: 1, yref: 'y', y0: y, y1: y, line: { color, width: 1.2, dash } };
  }

  function figures(V, Eng, DS) {
    const s = PR.series;
    const iso = V.profile_isentropic;
    PR.plot('fig-iso-M', [
      { x: iso.x_fine, y: iso.M_fine, name: 'exact', mode: 'lines', line: { color: s(1), width: 2 } },
      { x: iso.x, y: iso.M, name: 'CFD N=100', mode: 'markers', marker: { color: s(0), size: 5 } },
    ], { xaxis: { title: 'x' }, yaxis: { title: 'Mach number' } });
    PR.plot('fig-iso-p', [
      { x: iso.x_fine, y: iso.p_fine, name: 'exact', mode: 'lines', line: { color: s(1), width: 2 } },
      { x: iso.x, y: iso.p, name: 'CFD N=100', mode: 'markers', marker: { color: s(0), size: 5 } },
    ], { xaxis: { title: 'x' }, yaxis: { title: 'p / p<sub>0</sub>' } });

    const G = V.grid_isentropic, dx = G.dx;
    const mu = G.vanalbada.L1_M, fo = G.none.L1_M, n = dx.length - 1;
    PR.plot('fig-grid', [
      { x: dx, y: mu, name: 'MUSCL + HLLC (van Albada)', mode: 'lines+markers', line: { color: s(0) } },
      { x: dx, y: fo, name: 'first order', mode: 'lines+markers', line: { color: s(1) } },
      { x: dx, y: dx.map((h) => mu[n] * Math.pow(h / dx[n], 2)), name: 'slope 2', mode: 'lines', line: { color: s(0), dash: 'dot', width: 1 } },
      { x: dx, y: dx.map((h) => fo[n] * (h / dx[n])), name: 'slope 1', mode: 'lines', line: { color: s(1), dash: 'dot', width: 1 } },
    ], { xaxis: { title: 'Δx', type: 'log', tickvals: dx, ticktext: dx.map((h) => h.toPrecision(2)) },
         yaxis: { title: 'L1 error in Mach', type: 'log', exponentformat: 'power', dtick: 1 },
         legend: { orientation: 'v', x: 0.99, xanchor: 'right', y: 0.03, yanchor: 'bottom' } });

    const sh = V.profile_shock;
    PR.plot('fig-sh-M', [
      { x: sh.x_fine, y: sh.M_fine, name: 'exact', mode: 'lines', line: { color: s(1), width: 2 } },
      { x: sh.x, y: sh.M, name: 'CFD N=200', mode: 'markers', marker: { color: s(0), size: 4 } },
    ], { xaxis: { title: 'x' }, yaxis: { title: 'Mach number' } });
    PR.plot('fig-sh-p', [
      { x: sh.x_fine, y: sh.p_fine, name: 'exact', mode: 'lines', line: { color: s(1), width: 2 } },
      { x: sh.x, y: sh.p, name: 'CFD N=200', mode: 'markers', marker: { color: s(0), size: 4 } },
    ], { xaxis: { title: 'x' }, yaxis: { title: 'p / p<sub>0</sub>' } });

    const sw = V.shock_sweep, cr = V.shock_case.crit;
    PR.plot('fig-sweep', [
      { x: sw.pb_fine, y: sw.x_exact_fine, name: 'exact', mode: 'lines', line: { color: s(1), width: 2 } },
      { x: sw.pb, y: sw.x_cfd, name: 'CFD N=200', mode: 'markers', marker: { color: s(0), size: 8, symbol: 'circle-open', line: { width: 2 } } },
    ], {
      xaxis: { title: 'back pressure p<sub>b</sub> / p<sub>0</sub>' }, yaxis: { title: 'shock station x<sub>s</sub>', range: [1.45, 3.05] },
      shapes: [cr.p_nse, cr.p_sub].map((x) => ({ type: 'line', x0: x, x1: x, yref: 'paper', y0: 0, y1: 1, line: { color: PR.css('--ink-3'), dash: 'dot', width: 1 } })),
      annotations: [{ x: cr.p_nse, y: 1.5, text: 'p<sub>NSE</sub>', showarrow: false, xanchor: 'left', font: { size: 11, color: PR.css('--ink-2') } },
                    { x: cr.p_sub, y: 3.0, text: 'p<sub>sub</sub>', showarrow: false, xanchor: 'right', font: { size: 11, color: PR.css('--ink-2') } }],
    });

    const H = V.residual_history;
    const lbl = { isentropic_anderson_N200: 'Anderson isentropic, N=200', shock_N200: 'shock case, N=200', merlin_N400: 'Merlin-class, γ=1.2, N=400' };
    PR.plot('fig-res', Object.keys(H).map((k, i) => ({ x: H[k].iter, y: H[k].res, name: lbl[k], mode: 'lines', line: { color: s(i), width: 1.6 } })),
      { xaxis: { title: 'iteration' }, yaxis: { title: 'RMS residual (rel.)', type: 'log', exponentformat: 'power', dtick: 2 } });

    // engines
    const E = Eng.engines;
    const tr1 = [], shapes6 = [hline(1, PR.css('--ink-3'), 'dot'), hline(Eng.summerfield, PR.css('--bad'), 'dash')];
    E.forEach((e) => {
      const c = s(ENG_COL[e.key]);
      tr1.push({ x: e.sweep.z_km, y: e.sweep.pe_pa, name: e.name, mode: 'lines', line: { color: c, width: e.key === 'merlin' ? 4.5 : 2, dash: e.key === 'raptor' ? 'dash' : 'solid' } });
      shapes6.push({ type: 'line', xref: 'x', x0: 0, x1: 4, yref: 'y', y0: e.schmucker_ratio, y1: e.schmucker_ratio, line: { color: c, width: 2, dash: 'dot' } });
    });
    PR.plot('fig-pepa', tr1, {
      xaxis: { title: 'altitude (km)', range: [0, 50] }, yaxis: { title: 'p<sub>e</sub> / p<sub>a</sub>', type: 'log', exponentformat: 'power', dtick: 1 }, shapes: shapes6,
      annotations: [
        { xref: 'paper', x: 1, y: Math.log10(1), text: 'p<sub>e</sub> = p<sub>a</sub>', showarrow: false, xanchor: 'right', yanchor: 'bottom', font: { size: 11, color: PR.css('--ink-2') } },
        { xref: 'paper', x: 1, y: Math.log10(Eng.summerfield), text: 'Summerfield 0.4', showarrow: false, xanchor: 'right', yanchor: 'top', font: { size: 11, color: PR.css('--bad') } },
      ],
    });

    const tr2 = [];
    E.forEach((e) => {
      const c = s(ENG_COL[e.key]);
      tr2.push({ x: e.sweep.z_km, y: e.sweep.cf, name: e.name, mode: 'lines', line: { color: c, width: 2 } });
      const zr = [], cr2 = [];
      e.sweep.summerfield_risk.forEach((r, i) => { if (r) { zr.push(e.sweep.z_km[i]); cr2.push(e.sweep.cf[i]); } });
      if (zr.length) tr2.push({ x: zr, y: cr2, name: e.name + ': Summerfield flag', mode: 'lines', line: { color: PR.css('--bad'), width: 6 }, opacity: 0.35, showlegend: true });
      const k = e.sweep.z_km.map((_, i) => i).filter((i) => i % 10 === 0);
      tr2.push({ x: k.map((i) => e.sweep.z_km[i]), y: k.map((i) => e.sweep.cf_ideal[i]), name: 'ideal', mode: 'markers', showlegend: false, marker: { color: c, size: 6, symbol: 'circle-open' } });
    });
    PR.plot('fig-cf', tr2, { xaxis: { title: 'altitude (km)', range: [0, 50] }, yaxis: { title: 'thrust coefficient C<sub>F</sub>' } });

    const tr3 = [];
    E.forEach((e) => {
      const c = s(ENG_COL[e.key]);
      tr3.push({ x: e.shock_curve.pb_over_pa_sl, y: e.shock_curve.As_over_eps, name: e.name, mode: 'lines', line: { color: c, width: 2 } });
      tr3.push({ x: [e.p_sep_summerfield_kPa * 1e3 / Eng.p_sl], y: [1], mode: 'markers', showlegend: false, name: e.name + ' Summerfield onset',
                 marker: { symbol: 'diamond', size: 11, color: c, line: { color: PR.css('--ink'), width: 1 } } });
    });
    PR.plot('fig-shockeng', tr3, {
      xaxis: { title: 'back pressure / sea-level ambient', type: 'log' }, yaxis: { title: 'A<sub>s</sub> / ε  (1 = exit plane)', range: [0, 1.08] },
      shapes: [{ type: 'line', x0: 1, x1: 1, yref: 'paper', y0: 0, y1: 1, line: { color: PR.css('--ink-3'), dash: 'dot' } }],
      annotations: [{ x: 0, y: 0.08, xref: 'x', text: 'sea level', showarrow: false, xanchor: 'left', font: { size: 11, color: PR.css('--ink-2') } }],
    });

    const tr4 = [
      { x: DS.Pc_MPa, y: DS.eps_summerfield, name: 'Summerfield ε<sub>max</sub> (γ 1.2)', mode: 'lines', line: { color: s(0), width: 2.5 } },
      { x: DS.Pc_MPa, y: DS.eps_summerfield_g115, name: 'Summerfield, γ 1.15 / 1.25', mode: 'lines', line: { color: s(0), width: 1, dash: 'dot' } },
      { x: DS.Pc_MPa, y: DS.eps_summerfield_g125, name: 'γ 1.25', showlegend: false, mode: 'lines', line: { color: s(0), width: 1, dash: 'dot' } },
      { x: DS.Pc_MPa, y: DS.eps_schmucker, name: 'Schmucker ε<sub>max</sub>', mode: 'lines', line: { color: s(3), width: 2, dash: 'dash' } },
      { x: DS.Pc_MPa, y: DS.eps_optimal, name: 'p<sub>e</sub> = p<sub>a</sub> at sea level', mode: 'lines', line: { color: PR.css('--ink-3'), width: 1.2 } },
    ];
    DS.engines.forEach((e) => tr4.push({ x: [e.Pc_MPa], y: [e.eps], mode: 'markers+text', name: e.name, showlegend: false,
      text: [e.name.replace(/ \(.*\)/, '').replace('-class', '')], textposition: 'top left', textfont: { size: 11, color: PR.css('--ink') },
      marker: { size: 12, symbol: 'star', color: s(ENG_COL[e.key]), line: { color: PR.css('--ink'), width: 1 } } }));
    PR.plot('fig-design', tr4, { xaxis: { title: 'chamber pressure p<sub>c</sub> (MPa)', type: 'log' }, yaxis: { title: 'expansion ratio ε', type: 'log' } });
  }

  function tables(V, Eng) {
    const A = V.anderson;
    const rel = (a, b) => sci(Math.abs(a / b - 1), 2);
    fillTable('tab-anderson', [
      ['exit Mach M<sub>e</sub>', fx(A.M_exit_cfd100, 4), fx(A.M_exit_exact, 4), rel(A.M_exit_cfd100, A.M_exit_exact)],
      ['exit pressure p<sub>e</sub>/p<sub>0</sub>', fx(A.p_exit_cfd100, 5), fx(A.p_exit_exact, 5), rel(A.p_exit_cfd100, A.p_exit_exact)],
      ['mass flow ṁ/(p<sub>0</sub>A*/√(RT<sub>0</sub>))', fx(A.mdot_cfd100, 5), fx(A.mdot_exact, 5), rel(A.mdot_cfd100, A.mdot_exact)],
      ['ṁ spread across all faces', sci(A.mdot_spread100, 2), '0', '–'],
    ]);
    const G = V.grid_isentropic, a = G.vanalbada, f = G.none, o = (v) => (v === null ? '–' : fx(v, 2));
    fillTable('tab-grid', G.N.map((N, i) => [{ t: N, cls: 'num' }, sci(a.L1_M[i], 3), o(a.order_L1_M[i]), sci(a.L1_p[i], 3), o(a.order_L1_p[i]),
      sci(a.mdot_err[i], 3), o(a.order_mdot[i]), sci(f.L1_M[i], 3), o(f.order_L1_M[i]), a.iters[i]]));
    const S = V.shock_case;
    fillTable('tab-shock', S.N.map((N, i) => [{ t: N, cls: 'num' }, fx(S.x_shock_cfd[i], 5), sci(S.x_shock_err[i], 2), fx(S.x_shock_err_cells[i], 3),
      o(S.order_x_shock[i]), sci(S.mdot_err[i], 2), sci(S.mdot_spread[i], 2), sci(S.p_exit_err[i], 2), S.iters[i]]));

    const E = Eng.engines;
    fillTable('tab-inputs', E.map((e) => [e.name, { t: e.prop, cls: '' }, fx(e.Pc / 1e6, 1), fx(e.eps, 0), fx(e.gamma, 2), { t: e.source + '. ' + e.note, cls: 'wrapcell' }]));
    fillTable('tab-engcfd', E.map((e) => [e.name, fx(e.cfd.M_exit, 4), sci(e.err.M_exit, 2), sci(e.err.pe, 2), fx(e.cf_vac, 4), sci(e.err.cf_vac, 2),
      sci(e.err.mdot, 2), sci(e.shock_check.As_err_rel, 2), e.cfd.iters]));
    const zc = (v) => (v === null ? '> 50' : v === 0 ? 'at SL' : fx(v, 1));
    fillTable('tab-eng', E.map((e) => [e.name, fx(e.pe_kPa, 1), fx(e.pe_pa_sl, 3), flag(e.summerfield_risk_sl), fx(e.schmucker_ratio, 3),
      flag(e.schmucker_risk_sl), fx(e.z_opt_km, 1), zc(e.z_summerfield_clear_km), fx(e.cf_sl, 4), fx(e.cf_vac, 4)]));
    fillTable('tab-crit', E.map((e) => [e.name, fx(e.p_nse_MPa, 3), fx(e.p_nse_over_pa_sl, 1), fx(e.p_mid_MPa, 3), fx(e.p_sub_MPa, 2), fx(e.p_sep_summerfield_kPa, 1)]));
    fillTable('tab-margin', E.map((e) => {
      const r = e.eps_over_eps_max, cls = r > 1 ? 'num fail' : r > 0.85 ? 'num marg' : 'num pass';
      const t = e.throttle_margin, cls2 = t > 1 ? 'num fail' : t > 0.85 ? 'num marg' : 'num pass';
      return [e.name, fx(e.eps, 0), fx(e.eps_max_summerfield, 1), { t: fx(r, 2), cls }, fx(e.eps_max_schmucker, 1), fx(e.eps_opt_sl, 1),
        fx(e.pc_min_summerfield_MPa, 1), { t: fx(t, 2), cls: cls2 }];
    }));
    const rows = [];
    E.forEach((e) => e.gamma_sensitivity.forEach((g, i) => rows.push([i === 0 ? e.name : '', fx(g.gamma, 2), fx(g.Me, 3), fx(g.pe_kPa, 1), fx(g.pe_pa_sl, 3),
      flag(g.pe_pa_sl < Eng.summerfield), fx(g.schmucker, 3)])));
    fillTable('tab-gamma', rows);
  }

  function crossText(XC) {
    const el = $('xc-text');
    if (!XC.available || !XC.cases.length) { el.textContent = 'The Node cross-check against the Python solver was not available when the data were generated.'; return; }
    const same = XC.cases.every((c) => c.iters_py === c.iters_js);
    el.innerHTML = `I checked the port against the Python solver by running it under Node (${XC.node}) on ${XC.cases.length} cases (${XC.cases.map((c) => c.name).join('; ')}). ` +
      `${same ? 'Both codes take the same number of iterations in every case, and' : 'Across the cases'} the converged Mach fields agree to max |ΔM| = ${sci(XC.max_dM, 2)}, i.e. round-off.`;
  }

  async function loadAll() {
    try {
      const [V, Eng, DS, XC, sum] = await Promise.all(['verification', 'engines', 'design_space', 'crosscheck', 'summary']
        .map((n) => PR.loadJSON('data/' + n + '.json')));
      const E = {}; Eng.engines.forEach((e) => { E[e.key] = e; });
      const root = {
        V, Eng, DS, XC, sum, E,
        maxInv: Math.max(...Eng.engines.map((e) => e.cfd.invariance_rel)),
        maxAsErr: Math.max(...Eng.engines.map((e) => e.shock_check.As_err_rel)),
        minPnse: Math.min(...Eng.engines.map((e) => e.p_nse_over_pa_sl)),
        maxPnse: Math.max(...Eng.engines.map((e) => e.p_nse_over_pa_sl)),
      };
      fillSpans(root);
      tables(V, Eng);
      crossText(XC);
      figures(V, Eng, DS);
    } catch (err) {
      console.error(err);
      document.querySelectorAll('[data-v]').forEach((el) => { if (!el.textContent) el.textContent = '–'; });
      $('xc-text').textContent = 'Could not load the study data (' + err.message + ').';
    }
  }

  // ============================================================ interactive tool
  const S = window.NozzleCFD;
  const PRESETS = {
    merlin: { geom: 'bell', eps: 16, gamma: 1.2, bmode: 'alt', pc: 9.7, alt: 0, N: 160 },
    rs25: { geom: 'bell', eps: 69, gamma: 1.2, bmode: 'alt', pc: 20.6, alt: 0, N: 160 },
    raptor: { geom: 'bell', eps: 40, gamma: 1.2, bmode: 'alt', pc: 30, alt: 0, N: 160 },
    anderson: { geom: 'anderson', eps: 5.95, gamma: 1.4, bmode: 'ratio', pb: 0.01, N: 100 },
    andshock: { geom: 'anderson', eps: 1.5002, gamma: 1.4, bmode: 'ratio', pb: 0.6784, N: 120 },
  };
  const st = { geom: 'bell', eps: 16, gamma: 1.2, bmode: 'alt', pc: 9.7, alt: 0, pb: 0.01, N: 160, speed: 12, dt: 'global' };
  let solver = null, geom = null, ex = null, running = true, converged = false, jumping = false;
  let hist = { it: [], r: [] }, stride = 1, lastPlot = 0, xsFine = [];

  function pbNow() { return st.bmode === 'alt' ? S.atmosphere(st.alt * 1000).p / (st.pc * 1e6) : st.pb; }

  function syncControls() {
    $('c-geom').value = st.geom; $('c-eps').value = st.eps; $('c-gamma').value = st.gamma; $('c-bmode').value = st.bmode;
    $('c-pc').value = st.pc; $('c-alt').value = st.alt; $('c-pb').value = Math.log10(st.pb); $('c-N').value = st.N;
    $('c-speed').value = st.speed; $('c-dt').value = st.dt;
    $('w-pc').classList.toggle('hidden', st.bmode !== 'alt'); $('w-alt').classList.toggle('hidden', st.bmode !== 'alt');
    $('w-pb').classList.toggle('hidden', st.bmode !== 'ratio');
    labels();
  }
  function labels() {
    $('o-eps').textContent = st.eps.toFixed(st.eps < 10 ? 2 : 1); $('o-gamma').textContent = st.gamma.toFixed(2);
    $('o-pc').textContent = st.pc.toFixed(1) + ' MPa'; $('o-alt').textContent = st.alt.toFixed(1) + ' km';
    $('o-pb').textContent = sci(st.pb, 3); $('o-N').textContent = st.N; $('o-speed').textContent = st.speed;
    $('r-pepa-l').innerHTML = st.bmode === 'alt' ? 'p<sub>e</sub>/p<sub>a</sub>' : 'p<sub>e</sub>/p<sub>b</sub>';
  }

  function rebuild() {
    geom = st.geom === 'anderson' ? S.andersonNozzle(st.eps) : S.bellNozzle(st.eps);
    solver = new S.Solver(geom, st.N, st.gamma, pbNow(), { localDt: st.dt === 'local', init: 'rest', cfl: 0.8 });
    xsFine = Array.from({ length: 301 }, (_, i) => geom.L * i / 300);
    hist = { it: [], r: [] }; stride = 1; converged = false;
    exactUpdate();
    render(true);
  }
  function exactUpdate() { ex = S.exactNozzle(xsFine, geom, solver.pb, st.gamma); }
  function newPb() { solver.pb = pbNow(); converged = false; exactUpdate(); render(true); }

  function record(it, r) {
    if (it % stride) return;
    hist.it.push(it); hist.r.push(r);
    if (hist.it.length > 1500) { hist.it = hist.it.filter((_, i) => i % 2 === 0); hist.r = hist.r.filter((_, i) => i % 2 === 0); stride *= 2; }
  }

  function march(budgetMs, maxSteps) {
    const t0 = performance.now();
    let n = 0;
    try {
      while (n < maxSteps && performance.now() - t0 < budgetMs) {
        const r = solver.step(); n++;
        record(solver.iter, r);
        if (r < 1e-10 && solver.iter > 50) { converged = true; break; }
        if (solver.iter >= 300000) { converged = true; break; }
      }
    } catch (e) {
      running = false; $('b-run').textContent = 'Run';
      $('r-regime').textContent = 'The solver diverged for these inputs; try a smaller CFL-limited change or restart.';
    }
  }

  let visible = true;
  function loop() {
    if (visible && solver && (running || jumping) && !converged) {
      march(jumping ? 40 : 12, jumping ? 1e9 : st.speed);
      if (converged) jumping = false;
    }
    if (visible) render(false);
    requestAnimationFrame(loop);
  }

  // colour map (viridis stops)
  const VIR = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]];
  function cmap(t) {
    t = Math.min(Math.max(t, 0), 1) * (VIR.length - 1);
    const i = Math.min(Math.floor(t), VIR.length - 2), f = t - i, a = VIR[i], b = VIR[i + 1];
    return `rgb(${a.map((v, k) => Math.round(v + f * (b[k] - v))).join(',')})`;
  }

  function drawNozzle(sol) {
    const cv = $('nz-canvas'), dpr = window.devicePixelRatio || 1;
    const W = cv.clientWidth, H = cv.clientHeight;
    if (!W) return;
    if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(H * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); }
    const g = cv.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, W, H);
    const padL = 10, padR = 70, padY = 14, cy = H / 2;
    const w = W - padL - padR;
    let Amax = 0; for (let i = 0; i <= solver.N; i++) Amax = Math.max(Amax, solver.Af[i]);
    const rs = (A) => Math.sqrt(A / Amax) * (H / 2 - padY);
    let Mmax = 1; for (const m of ex.M) Mmax = Math.max(Mmax, m); for (const m of sol.M) if (isFinite(m)) Mmax = Math.max(Mmax, m);
    const N = solver.N, L = geom.L;
    for (let i = 0; i < N; i++) {
      const x0 = padL + w * solver.xf[i] / L, x1 = padL + w * solver.xf[i + 1] / L + 0.6;
      const r0 = rs(solver.Af[i]), r1 = rs(solver.Af[i + 1]);
      g.fillStyle = cmap(sol.M[i] / Mmax);
      g.beginPath(); g.moveTo(x0, cy - r0); g.lineTo(x1, cy - r1); g.lineTo(x1, cy + r1); g.lineTo(x0, cy + r0); g.closePath(); g.fill();
    }
    g.strokeStyle = PR.css('--ink'); g.lineWidth = 2;
    for (const sgn of [-1, 1]) {
      g.beginPath();
      for (let i = 0; i <= N; i++) { const x = padL + w * solver.xf[i] / L, y = cy + sgn * rs(solver.Af[i]); i ? g.lineTo(x, y) : g.moveTo(x, y); }
      g.stroke();
    }
    g.setLineDash([4, 4]); g.strokeStyle = PR.css('--line-2'); g.lineWidth = 1;
    g.beginPath(); g.moveTo(padL, cy); g.lineTo(padL + w, cy); g.stroke(); g.setLineDash([]);
    const xs = S.shockLocation(sol.x, sol.M, geom.x_throat);
    if (xs !== null) {
      const X = padL + w * xs / L, r = rs(geom.A(xs));
      g.strokeStyle = '#ffffff'; g.lineWidth = 2.5; g.beginPath(); g.moveTo(X, cy - r); g.lineTo(X, cy + r); g.stroke();
    }
    // colour bar
    const bx = W - padR + 22, by = padY, bh = H - 2 * padY;
    for (let k = 0; k < bh; k++) { g.fillStyle = cmap(1 - k / bh); g.fillRect(bx, by + k, 12, 1.2); }
    g.fillStyle = PR.css('--ink-2'); g.font = '11px "JetBrains Mono", monospace'; g.textBaseline = 'middle';
    g.fillText(Mmax.toFixed(1), bx + 16, by + 6);
    g.fillText('0', bx + 16, by + bh - 4);
    g.save(); g.translate(bx - 8, by + bh / 2); g.rotate(-Math.PI / 2); g.textAlign = 'center'; g.fillText('Mach', 0, 0); g.restore();
  }

  function render(force) {
    if (!solver) return;
    const now = performance.now();
    if (!force && now - lastPlot < 110) return;
    lastPlot = now;
    const sol = solver.state();
    const pb = solver.pb, g = st.gamma, eps = geom.A(geom.L) / geom.A(geom.x_throat);
    drawNozzle(sol);
    // readouts
    $('r-iter').textContent = solver.iter;
    $('r-res').textContent = solver.iter ? sci(solver.res, 2) : '–';
    const MeEx = ex.M[ex.M.length - 1];
    $('r-Me').textContent = `${fx(sol.exit.M, 3)} (${fx(MeEx, 3)})`;
    const pepa = sol.exit.p / pb;
    $('r-pepa').textContent = fx(pepa, 3);
    $('r-cf').textContent = fx(sol.momExit - pb * eps, 4);
    const xs = S.shockLocation(sol.x, sol.M, geom.x_throat);
    const AsC = xs === null ? null : geom.A(xs);
    $('r-shock').textContent = (AsC === null ? 'none' : fx(AsC, 2)) + ' (' + (ex.As === null ? 'none' : fx(ex.As, 2)) + ')';
    $('r-mdot').textContent = sci(sol.mdotOut / S.mdotStar(g) - 1, 2);
    const fl = $('r-flag');
    if (!sol.exitSupersonic || xs !== null) { fl.textContent = 'n/a'; fl.className = 'v flag-na'; }
    else if (pepa < 0.4) { fl.textContent = 'RISK'; fl.className = 'v flag-risk'; }
    else { fl.textContent = 'clear'; fl.className = 'v flag-ok'; }
    const c = ex.crit;
    const reg = ex.regime === 'supersonic' ? (pb > c.p_sup ? 'supersonic exit, over-expanded' : 'supersonic exit, under-expanded')
      : ex.regime === 'shock' ? 'normal shock in the divergent section' : 'subsonic throughout (not choked)';
    $('r-regime').textContent = `Exact solution: ${reg}. p_b/p₀ = ${sci(pb, 3)}; critical values p_sup/p₀ = ${sci(c.p_sup, 3)}, ` +
      `p_NSE/p₀ = ${sci(c.p_nse, 3)}, p_sub/p₀ = ${sci(c.p_sub, 3)}. ${ex.regime === 'supersonic' ? `Schmucker limit at this M_e: ${fx(S.schmucker(MeEx), 2)}. ` : ''}` +
      (converged ? 'Converged.' : running || jumping ? 'Marching…' : 'Paused.');
    // plots
    const xL = Array.from(sol.x, (x) => x / geom.L), xf = xsFine.map((x) => x / geom.L);
    const s = PR.series;
    PR.plot('t-M', [
      { x: xf, y: Array.from(ex.M), name: 'exact steady', mode: 'lines', line: { color: s(1), dash: 'dash', width: 1.6 } },
      { x: xL, y: Array.from(sol.M), name: 'CFD (current iterate)', mode: 'lines', line: { color: s(0), width: 2.2 } },
    ], { xaxis: { title: 'x / L', range: [0, 1] }, yaxis: { title: 'Mach' }, margin: { t: 34 } }, { displayModeBar: false });
    PR.plot('t-p', [
      { x: xf, y: Array.from(ex.p), name: 'exact steady', mode: 'lines', line: { color: s(1), dash: 'dash', width: 1.6 } },
      { x: xL, y: Array.from(sol.p), name: 'CFD', mode: 'lines', line: { color: s(0), width: 2.2 } },
      { x: [0, 1], y: [pb, pb], name: 'p<sub>b</sub>', mode: 'lines', line: { color: s(5), width: 1, dash: 'dot' } },
    ], { xaxis: { title: 'x / L', range: [0, 1] }, yaxis: { title: 'p / p<sub>0</sub>', type: 'log', exponentformat: 'power', dtick: 1 }, margin: { t: 34 } }, { displayModeBar: false });
    PR.plot('t-res', [{ x: hist.it.slice(), y: hist.r.slice(), name: 'residual', mode: 'lines', line: { color: s(3), width: 1.5 } }],
      { xaxis: { title: 'iteration' }, yaxis: { title: 'residual', type: 'log', exponentformat: 'power', dtick: 1 }, showlegend: false, margin: { t: 12 } }, { displayModeBar: false });
  }

  function initTool() {
    if (!S || !window.Plotly) return;
    const on = (id, ev, fn) => $(id).addEventListener(ev, fn);
    on('c-preset', 'change', (e) => { Object.assign(st, PRESETS[e.target.value]); syncControls(); running = true; $('b-run').textContent = 'Pause'; rebuild(); });
    on('c-geom', 'change', (e) => { st.geom = e.target.value; rebuild(); });
    on('c-eps', 'input', (e) => { st.eps = +e.target.value; labels(); });
    on('c-eps', 'change', () => rebuild());
    on('c-gamma', 'input', (e) => { st.gamma = +e.target.value; labels(); });
    on('c-gamma', 'change', () => rebuild());
    on('c-N', 'input', (e) => { st.N = +e.target.value; labels(); });
    on('c-N', 'change', () => rebuild());
    on('c-bmode', 'change', (e) => {
      st.bmode = e.target.value;
      if (st.bmode === 'ratio') st.pb = Math.min(Math.max(S.atmosphere(st.alt * 1000).p / (st.pc * 1e6), 3.2e-4), 0.99);
      syncControls(); newPb();
    });
    on('c-pc', 'input', (e) => { st.pc = +e.target.value; labels(); newPb(); });
    on('c-alt', 'input', (e) => { st.alt = +e.target.value; labels(); newPb(); });
    on('c-pb', 'input', (e) => { st.pb = Math.pow(10, +e.target.value); labels(); newPb(); });
    on('c-speed', 'input', (e) => { st.speed = +e.target.value; labels(); });
    on('c-dt', 'change', (e) => { st.dt = e.target.value; solver.localDt = st.dt === 'local'; converged = false; });
    on('b-run', 'click', () => { running = !running; jumping = false; $('b-run').textContent = running ? 'Pause' : 'Run'; if (running) converged = false; });
    on('b-restart', 'click', () => { running = true; $('b-run').textContent = 'Pause'; rebuild(); });
    on('b-converge', 'click', () => { jumping = true; converged = false; });
    window.addEventListener('resize', () => render(true));
    window.matchMedia('(prefers-color-scheme: light)').addEventListener('change', () => render(true));
    Object.assign(st, PRESETS.merlin);
    syncControls();
    rebuild();
    // pause the animation while the tool is off-screen
    if ('IntersectionObserver' in window) {
      new IntersectionObserver((en) => { visible = en[0].isIntersecting; }).observe($('nz-tool'));
    }
    requestAnimationFrame(loop);
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', () => { loadAll(); initTool(); });
  else { loadAll(); initTool(); }
})();
