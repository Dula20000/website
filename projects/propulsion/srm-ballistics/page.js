/* Research page: fills every number from data/*.json, draws the figures, and runs the interactive tool
   on the JavaScript solver (srm.js). */
(async function () {
  'use strict';
  const $ = (id) => document.getElementById(id);
  const D = {};
  try {
    const files = ['verification', 'results', 'montecarlo', 'cases', 'js_crosscheck'];
    const got = await Promise.all(files.map((f) => PR.loadJSON('data/' + f + '.json')));
    [D.v, D.r, D.mc, D.c, D.js] = got;
  } catch (e) {
    const s = $('status');
    if (s) { s.textContent = 'could not load data: ' + e.message; s.className = 'status warn'; }
    return;
  }

  // ------------------------------------------------------------------ formatting + binding
  const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  function sci(v, d) {
    if (v === 0) return '0';
    let e = Math.floor(Math.log10(Math.abs(v)));
    let m = v / Math.pow(10, e);
    if (Math.abs(+m.toFixed(d)) >= 10) { e += 1; m /= 10; }
    return m.toFixed(d) + ' × 10' + String(e).split('').map((c) => SUP[c]).join('');
  }
  function bigNum(v) {
    if (Math.abs(v) >= 1e7) return sci(v, 2);
    return Math.round(v).toLocaleString('en-US');
  }
  const FMT = {
    s: (v) => String(v),
    f0: (v) => v.toFixed(0), f1: (v) => v.toFixed(1), f2: (v) => v.toFixed(2), f3: (v) => v.toFixed(3),
    int: (v) => Math.round(v).toLocaleString('en-US'),
    big: bigNum,
    pct0: (v) => (100 * v).toFixed(0) + ' %', pct1: (v) => (100 * v).toFixed(1) + ' %', pct2: (v) => (100 * v).toFixed(2) + ' %',
    mpa2: (v) => (v / 1e6).toFixed(2), mpa3: (v) => (v / 1e6).toFixed(3),
    e1: (v) => sci(v, 1),
  };
  const get = (path) => path.split('.').reduce((o, k) => (o == null ? undefined : o[k]), D);

  // ------------------------------------------------------------------ derived quantities (from the JSON)
  const mc = D.mc, R = D.r, V = D.v;
  const src = mc.src_lnPc, srcA = mc.anchor.src_lnPc;
  const params = ['a', 'n', 'Dt', 'rho'];
  const varTot = params.reduce((s, p) => s + src[p] * src[p], 0);
  const varTotA = params.reduce((s, p) => s + srcA[p] * srcA[p], 0);
  const n0 = mc.n0;
  const lever = Math.log(mc.nominal.Pc_max / 1e6);
  const knsbSpec = D.c.presets.knsb;
  const b = R.bates;
  const ruleNearest = b.rule_nseg.map((x) => Math.round(x));
  const nMatch = ruleNearest.filter((x, i) => x === b.best_nseg_per_d[i]).length;
  // SRB bucket: first peak, minimum, second peak of the ODE thrust
  const so = R.cases.srb.ode, sm = R.cases.srb.metrics;
  const tb = sm.t_end;
  let i1 = 0, i2 = 0;
  so.t.forEach((t, i) => { if (t < 0.45 * tb && so.F[i] > so.F[i1]) i1 = i; });
  so.t.forEach((t, i) => { if (t > 0.45 * tb && t < 0.85 * tb && (i2 === 0 || so.F[i] > so.F[i2])) i2 = i; });
  let ib = i1;
  for (let i = i1; i <= i2; i++) if (so.F[i] < so.F[ib]) ib = i;
  D.x = {
    share_n_pct: 100 * src.n * src.n / varTot,
    star_closure_161: Math.abs(V.closure.star[1]),
    dI_knsb: Math.abs(R.cases.knsb.qs_vs_ode.dI), dI_apcp: Math.abs(R.cases.apcp.qs_vs_ode.dI),
    knsb_Lrule: 1e3 * (3 * 2 * knsbSpec.R + knsbSpec.segments[0].geom.core_d) / 2,
    throat_rate_mm: 1e3 * R.cases.apcp_throat.rate,
    srb_mass_t: sm.m_prop / 1e3, srb_Fmax_MN: sm.F_max / 1e6,
    srb_bucket_t: so.t[ib], srb_bucket_ratio: so.F[ib] / so.F[i1],
    geo_pt: R.geometry.Pc_target / 1e6,
    best_Lseg_mm: 1e3 * b.best.L_seg, best_Lrule_mm: 1e3 * b.best.L_rule,
    rule_lo: Math.min(...b.rule_nseg), rule_hi: Math.max(...b.rule_nseg),
    rule_match: nMatch + ' of ' + b.dD.length,
    mass_lo: b.mass[b.mass.length - 1][0], mass_hi: b.mass[0][0],
    sig_a: mc.sig.a, sig_Dt: mc.sig.Dt, sig_rho: mc.sig.rho,
    cvI: mc.I.std / mc.I.mean,
    amp: 1 / (1 - n0),
    sweep_lo: mc.sweep.MEOP_factor[0], sweep_hi: mc.sweep.MEOP_factor[mc.sweep.MEOP_factor.length - 1],
    lever, n_equiv: mc.sig.n * lever,
    cv: mc.Pc_max.std / mc.Pc_max.mean, cv_anchor: mc.anchor.Pc_max.std / mc.anchor.Pc_max.mean,
    share_n_anchor: 100 * srcA.n * srcA.n / varTotA,
  };
  document.querySelectorAll('[data-v]').forEach((el) => {
    const [p, f] = el.dataset.v.split('|');
    const v = get(p);
    el.textContent = v === undefined || v === null ? '—' : (FMT[f || 'f2'] || FMT.f2)(v);
  });
  if (window.renderMathInElement) {
    renderMathInElement(document.body, { delimiters: [{ left: '$$', right: '$$', display: true }, { left: '\\(', right: '\\)', display: false }], throwOnError: false });
  }

  // ------------------------------------------------------------------ tables
  function fillTable(id, rows) {
    const tb = document.querySelector('#' + id + ' tbody');
    tb.innerHTML = rows.map((r) => '<tr>' + r.map((c) => {
      const cls = typeof c === 'object' && c !== null ? c.cls : '';
      const txt = typeof c === 'object' && c !== null ? c.t : c;
      return '<td' + (cls ? ' class="' + cls + '"' : '') + '>' + txt + '</td>';
    }).join('') + '</tr>').join('');
  }
  const num = (t) => ({ t, cls: 'num' });
  const ok = (b) => (b ? { t: 'pass', cls: 'pass' } : { t: 'fail', cls: 'fail' });

  const cons = V.conservation;
  const massMax = Math.max(...Object.values(cons).map((c) => Math.abs(c.mass_balance)));
  const predCv = Math.sqrt(mc.sig.a ** 2 + mc.sig.rho ** 2 + (2 * mc.sig.Dt) ** 2 + (mc.sig.n * lever) ** 2) / (1 - n0);
  const cvRatio = D.x.cv / predCv;
  const avgOrder = (a) => a.reduce((s, x) => s + x, 0) / a.length;
  fillTable('tab-ver', [
    ['Circular port perimeter vs 2π(r₀+w), N = 161', num(sci(V.circle.err_max[2], 1) + ' (order ' + V.circle.order[3].toFixed(2) + ')'), 'max rel. error &lt; 2×10⁻⁴, order 1.8–2.2', ok(V.circle.err_max[2] < 2e-4 && V.circle.order[3] > 1.8 && V.circle.order[3] < 2.2)],
    ['Circular port initial area, N = 161', num(sci(V.circle.err_area0[2], 1)), 'rel. error &lt; 10⁻³', ok(V.circle.err_area0[2] < 1e-3)],
    ['BATES A<sub>b</sub>(w) vs closed form (4 grains, N = 161)', num(sci(V.bates.err_max, 1)), 'max rel. error &lt; 10⁻³', ok(V.bates.err_max < 1e-3)],
    ['Star perimeter vs exact offset, FMM', num(sci(V.star.err_fmm[2], 1) + ' @201, ' + sci(V.star.err_fmm[4], 1) + ' @801'), 'mean rel. error &lt; 3×10⁻³ at N = 201, decreasing (mean order ' + avgOrder(V.star.order_fmm).toFixed(2) + ')', ok(V.star.err_fmm[2] < 3e-3 && V.star.err_fmm[4] < V.star.err_fmm[2])],
    ['FMM field vs exact distance, star (L1 / R)', num(sci(V.star.field_l1[4] / 0.045, 1) + ' @801'), 'converging (order ' + V.star.order_field_l1[3].toFixed(2) + ' on the finest pair)', ok(V.star.order_field_l1[3] > 0.5)],
    ['Raw co-area closure, circle / 8-pt star (N = 161)', num(sci(Math.abs(V.closure.bates[1]), 1) + ' / ' + (100 * Math.abs(V.closure.star[1])).toFixed(2) + ' %'), '&lt; 1.2 %, first order for corners', ok(Math.abs(V.closure.star[1]) < 0.012 && Math.abs(V.closure.star[2]) < 0.6 * Math.abs(V.closure.star[1]))],
    ['Transient mass balance ∫ṁ dt / m<sub>p</sub> − 1 (all cases)', num(sci(massMax, 1)), '&lt; 2×10⁻³', ok(massMax < 2e-3)],
    ['ODE vs equilibrium P<sub>c</sub> at mid-web (APCP / KNSB)', num((100 * R.cases.apcp.qs_vs_ode.mid_burn_dPc).toFixed(2) + ' % / ' + (100 * R.cases.knsb.qs_vs_ode.mid_burn_dPc).toFixed(2) + ' %'), 'APCP &lt; 1.5 %, KNSB &lt; 2 % (ρ<sub>g</sub>/ρ<sub>p</sub> term)', ok(Math.abs(R.cases.apcp.qs_vs_ode.mid_burn_dPc) < 0.015 && Math.abs(R.cases.knsb.qs_vs_ode.mid_burn_dPc) < 0.02)],
    ['Monte Carlo σ/μ of P<sub>c,max</sub> vs linearised prediction', num(D.x.cv.toFixed(4) + ' / ' + predCv.toFixed(4)), 'ratio within 8 % (skew is second order)', ok(Math.abs(cvRatio - 1) < 0.08)],
    ['Browser JS vs Python (KNSB + APCP, all metrics)', num(sci(D.js.max_rel_diff, 1)), 'max rel. difference &lt; 10⁻⁶', ok(D.js.max_rel_diff < 1e-6)],
  ]);

  const props = D.c.props;
  const lawTxt = (p) => p.laws.length > 1 ? 'piecewise, ' + p.laws.length + ' ranges (table below)'
    : 'a = ' + (p.laws[0][2] * 1e3).toFixed(2) + ' mm/s, n = ' + p.laws[0][3].toFixed(2);
  fillTable('tab-props', ['knsb', 'apcp', 'pban'].map((k) => [props[k].name, num(props[k].rho.toFixed(0)), lawTxt(props[k]),
    num(props[k].cstar.toFixed(0)), num(props[k].gamma.toFixed(3)), { t: D.c.prop_sources[k], cls: 'wrap-cell' }]));
  fillTable('tab-knsb', props.knsb.laws.map((l, i) => [
    i === 0 ? '&lt; ' + l[1] / 1e6 + ' (fit from ' + l[0] / 1e6 + ')' : i === props.knsb.laws.length - 1 ? '≥ ' + l[0] / 1e6 + ' (fit to ' + l[1] / 1e6 + ')' : l[0] / 1e6 + ' – ' + l[1] / 1e6,
    num((l[2] * 1e3).toFixed(3)), num(l[3].toFixed(3))]));
  function grainText(s) {
    const segs = s.segments;
    const g = segs[0].geom, mm = (v) => (v * 1e3).toFixed(v < 0.1 ? 1 : 0);
    if (segs.length > 1) return '11-pt star (' + segs[0].L + ' m) + ' + (segs.length - 1) + ' stepped-bore slices, bore ' + segs[1].geom.core_d + '–' + segs[segs.length - 1].geom.core_d + ' m';
    const n = segs[0].count || 1, ends = segs[0].ends ? segs[0].ends + ' end(s) burning' : 'ends inhibited';
    let d;
    if (g.type === 'bates') d = 'BATES, core ' + mm(g.core_d) + ' mm';
    else if (g.type === 'star') d = g.points + '-pt star, tip r ' + mm(g.r_tip) + ', valley r ' + mm(g.r_valley) + ' mm';
    else if (g.type === 'finocyl') d = 'core ' + mm(g.core_d) + ' mm + ' + g.slots + ' slots ' + mm(g.slot_w) + ' mm wide to r ' + mm(g.slot_r) + ' mm';
    else d = 'moon, core ' + mm(g.core_d) + ' mm offset ' + mm(g.offset) + ' mm';
    return n + ' × ' + d + ', L ' + mm(segs[0].L) + ' mm, ' + ends;
  }
  const motorRows = [['knsb', R.cases.knsb.metrics], ['apcp', R.cases.apcp.metrics], ['srb', R.cases.srb.metrics]]
    .concat(['bates', 'star', 'finocyl', 'moon'].map((k) => ['geo_' + k, R.geometry.items[k].metrics]));
  fillTable('tab-motors', motorRows.map(([k, m]) => {
    const s = D.c.presets[k];
    return [s.name, { t: grainText(s), cls: 'wrap-cell' }, num((2 * s.R * 1e3).toFixed(0)), num((s.Dt * 1e3).toFixed(2)), num(s.eps.toFixed(1)), num(m.m_prop < 1000 ? m.m_prop.toFixed(3) : bigNum(m.m_prop))];
  }));
  fillTable('tab-results', ['knsb', 'apcp', 'srb'].map((k) => {
    const m = R.cases[k].metrics;
    return [D.c.presets[k].name, m.class.length === 1 ? m.designation : '(off the scale)', num(bigNum(m.I_total)), num(bigNum(m.F_avg)),
      num(m.t_burn.toFixed(2)), num((m.Pc_max / 1e6).toFixed(2)), num(m.Isp.toFixed(1)),
      num(m.Kn_initial.toFixed(0) + ' / ' + m.Kn_max.toFixed(0)), num(m.neutrality_pc.toFixed(3))];
  }));
  const geoNames = { bates: 'Circular core (tube)', star: '8-point star', finocyl: 'Slotted core (finocyl)', moon: 'Moon burner' };
  fillTable('tab-geo', ['bates', 'star', 'finocyl', 'moon'].map((k) => {
    const it = R.geometry.items[k], m = it.metrics;
    return [geoNames[k], it.shape, num(m.shape_ratio.toFixed(2)), num(m.neutrality_pc.toFixed(2)), num((it.Dt * 1e3).toFixed(2)),
      num(m.t_burn.toFixed(2)), num(Math.round(m.I_total).toLocaleString('en-US')), num(m.Isp.toFixed(1))];
  }));
  const P = mc.Pc_max;
  fillTable('tab-mc', [
    ['Nominal peak P<sub>c</sub> (equilibrium)', num((mc.nominal.Pc_max / 1e6).toFixed(3) + ' MPa')],
    ['Mean / σ of peak P<sub>c</sub>', num((P.mean / 1e6).toFixed(3) + ' / ' + (P.std / 1e6).toFixed(3) + ' MPa (σ/μ ' + (100 * D.x.cv).toFixed(2) + ' %)')],
    ['99.865th percentile (3σ-equivalent MEOP)', num((P.p99865 / 1e6).toFixed(3) + ' MPa [' + (P.p99865_ci[0] / 1e6).toFixed(2) + ', ' + (P.p99865_ci[1] / 1e6).toFixed(2) + '] 95 %')],
    ['Normal-theory μ + 3σ', num((P.mu3s / 1e6).toFixed(3) + ' MPa')],
    ['MEOP / nominal peak', num(mc.MEOP_factor.toFixed(3))],
    ['Required burst pressure (FS ' + mc.FS_ultimate.toFixed(1) + ' × MEOP)', num((mc.burst_required / 1e6).toFixed(2) + ' MPa  (' + mc.design_factor.toFixed(2) + ' × nominal)')],
    ['Total impulse mean ± σ', num(Math.round(mc.I.mean).toLocaleString('en-US') + ' ± ' + Math.round(mc.I.std).toLocaleString('en-US') + ' N·s')],
    ['Total impulse 0.135 / 99.865 percentiles', num(Math.round(mc.I.p00135).toLocaleString('en-US') + ' / ' + Math.round(mc.I.p99865).toLocaleString('en-US') + ' N·s')],
    ['Tolerances referenced to operating pressure: MEOP', num((mc.anchor.Pc_max.p99865 / 1e6).toFixed(3) + ' MPa (σ/μ ' + (100 * D.x.cv_anchor).toFixed(2) + ' %)')],
  ]);

  // ------------------------------------------------------------------ figures
  const S = (i) => PR.series(i);
  const ink3 = () => PR.css('--ink-3');
  const lin = (x, y, name, color, extra) => Object.assign({ x, y, name, type: 'scatter', mode: 'lines', line: { color, width: 2 } }, extra || {});
  const logAx = (title) => ({ title, type: 'log', exponentformat: 'power' });

  // Fig 1 convergence
  {
    const hc = V.circle.h.map((h) => h * 1e3), hs = V.star.h.map((h) => h * 1e3);
    const ref = (hArr, e0, p) => hArr.map((h) => e0 * Math.pow(h / hArr[hArr.length - 1], p));
    PR.plot('fig-conv', [
      lin(hc, V.circle.err_max, 'circle perimeter (max)', S(0), { mode: 'lines+markers' }),
      lin(hs, V.star.err_fmm, 'star perimeter, FMM (mean)', S(1), { mode: 'lines+markers' }),
      lin(hs, V.star.err_exact, 'star, exact-distance contour', S(2), { mode: 'lines+markers', line: { color: S(2), width: 2, dash: 'dot' } }),
      lin(hc, ref(hc, V.circle.err_max[V.circle.err_max.length - 1], 2), 'slope 2', ink3(), { line: { color: ink3(), width: 1, dash: 'dash' } }),
      lin(hs, ref(hs, V.star.err_fmm[V.star.err_fmm.length - 1] * 0.6, 1), 'slope 1', ink3(), { line: { color: ink3(), width: 1, dash: 'dashdot' } }),
    ], { xaxis: logAx('grid spacing h [mm]'), yaxis: logAx('relative perimeter error') });
  }
  // Fig 2
  {
    const bt = V.bates, st = V.star.curve;
    const every = (a, k) => a.filter((_, i) => i % k === 0);
    PR.plot('fig-bates', [
      lin(bt.w, bt.Ab_exact.map((v) => v * 1e4), 'closed form', S(1)),
      { x: every(bt.w, 4), y: every(bt.Ab_grid, 4).map((v) => v * 1e4), name: 'level set', type: 'scatter', mode: 'markers', marker: { color: S(0), size: 6 } },
    ], { xaxis: { title: 'web burned w [mm]' }, yaxis: { title: 'A<sub>b</sub> [cm²]' } });
    PR.plot('fig-starP', [
      lin(st.w, st.P_exact, 'exact offset', S(1)),
      { x: st.w, y: st.P_fmm, name: 'FMM, N = 201', type: 'scatter', mode: 'markers', marker: { color: S(0), size: 6 } },
    ], { xaxis: { title: 'web burned w [mm]' }, yaxis: { title: 'perimeter [mm]' } });
  }
  // Fig 3
  function curvePlot(id, rec, fs, funit) {
    const o = rec.ode, q = rec.qs;
    PR.plot(id, [
      lin(o.t, o.F.map((v) => v / fs), 'thrust (ODE)', S(0)),
      lin(q.t, q.F.map((v) => v / fs), 'thrust (equilibrium)', S(0), { line: { color: S(0), width: 1.5, dash: 'dash' } }),
      lin(o.t, o.Pc.map((v) => v / 1e6), 'P<sub>c</sub> (ODE)', S(1), { yaxis: 'y2' }),
      lin(q.t, q.Pc.map((v) => v / 1e6), 'P<sub>c</sub> (equilibrium)', S(1), { yaxis: 'y2', line: { color: S(1), width: 1.5, dash: 'dash' } }),
    ], { xaxis: { title: 'time [s]' }, yaxis: { title: 'thrust [' + funit + ']', rangemode: 'tozero' },
      yaxis2: { title: 'P<sub>c</sub> [MPa]', overlaying: 'y', side: 'right', showgrid: false, rangemode: 'tozero' },
      margin: { r: 60 } });
  }
  curvePlot('fig-knsb', R.cases.knsb, 1, 'N');
  curvePlot('fig-apcp', R.cases.apcp, 1e3, 'kN');
  curvePlot('fig-srb', R.cases.srb, 1e6, 'MN');
  // Fig 4
  {
    const pc = (o) => o.Pc.map((v) => v / 1e6);
    PR.plot('fig-eros', [lin(R.cases.knsb.ode.t, pc(R.cases.knsb.ode), 'baseline', S(1)),
      lin(R.cases.knsb_erosive.ode.t, pc(R.cases.knsb_erosive.ode), 'erosive burning on', S(0))],
    { xaxis: { title: 'time [s]' }, yaxis: { title: 'P<sub>c</sub> [MPa]', rangemode: 'tozero' } });
    PR.plot('fig-throat', [lin(R.cases.apcp.ode.t, pc(R.cases.apcp.ode), 'baseline', S(1)),
      lin(R.cases.apcp_throat.ode.t, pc(R.cases.apcp_throat.ode), 'throat erosion on', S(0))],
    { xaxis: { title: 'time [s]' }, yaxis: { title: 'P<sub>c</sub> [MPa]', rangemode: 'tozero' } });
  }
  // contour helpers
  function circleTrace(r, cx, color) {
    const x = [], y = [];
    for (let i = 0; i <= 180; i++) { const a = 2 * Math.PI * i / 180; x.push(cx + r * Math.cos(a)); y.push(r * Math.sin(a)); }
    return { x, y, type: 'scatter', mode: 'lines', line: { color, width: 1.5 }, hoverinfo: 'skip', showlegend: false };
  }
  function contourTraces(levels, scale, dx) {
    return levels.map((lv, i) => ({
      x: lv.x.map((v) => (v === null ? null : v * scale + dx)), y: lv.y.map((v) => (v === null ? null : v * scale)),
      type: 'scatter', mode: 'lines', hoverinfo: 'skip', showlegend: false,
      line: { color: i === 0 ? S(0) : S(1), width: i === 0 ? 2.2 : 1.1 }, opacity: i === 0 ? 1 : 0.35 + 0.65 * (1 - i / levels.length),
    }));
  }
  const eqAx = { xaxis: { showgrid: false, zeroline: false }, yaxis: { showgrid: false, zeroline: false, scaleanchor: 'x', scaleratio: 1 } };
  {
    const Rm = R.cases.srb.R_mm / 1e3, off = 2.3 * Rm;
    PR.plot('fig-srbx', [circleTrace(Rm, 0, ink3()), circleTrace(Rm, off, ink3())]
      .concat(contourTraces(R.cases.srb.contours, 1e-3, 0)).concat(contourTraces(R.cases.srb.contours_aft, 1e-3, off)),
    Object.assign({ margin: { l: 40, r: 10, t: 10, b: 36 } }, eqAx, { xaxis: { showgrid: false, zeroline: false, title: 'x [m]' } }));
  }
  // Fig 6/7 geometry study
  {
    const host = $('geo-xsec');
    const keys = ['bates', 'star', 'finocyl', 'moon'];
    keys.forEach((k) => {
      const div = document.createElement('div');
      div.className = 'plot short'; div.id = 'geo-' + k;
      host.appendChild(div);
    });
    keys.forEach((k) => {
      const it = R.geometry.items[k];
      PR.plot('geo-' + k, [circleTrace(R.geometry.R_mm, 0, ink3())].concat(contourTraces(it.contours, 1, 0)),
        Object.assign({ margin: { l: 40, r: 10, t: 34, b: 30 }, title: { text: geoNames[k] + ' · ' + it.shape, font: { size: 13 }, x: 0.02, y: 0.97 } }, eqAx));
    });
    PR.plot('fig-geoF', keys.map((k, i) => lin(R.geometry.items[k].ode.t, R.geometry.items[k].ode.F.map((v) => v / 1e3), geoNames[k], S(i))),
      { xaxis: { title: 'time [s]' }, yaxis: { title: 'thrust [kN]', rangemode: 'tozero' } });
  }
  // Fig 8 BATES heat map
  {
    const z = b.neutrality.map((row) => row.map((v) => Math.log10(v)));
    const ticks = [1.2, 1.5, 2, 3, 5, 10, 30, 100];
    PR.plot('fig-bheat', [
      { type: 'heatmap', x: b.nseg, y: b.dD, z, colorscale: 'Viridis', reversescale: true, zmin: Math.log10(1.1), zmax: 2,
        customdata: b.neutrality, hovertemplate: 'N = %{x}<br>d/D = %{y}<br>P<sub>c</sub> max/min = %{customdata:.3f}<extra></extra>',
        colorbar: { title: { text: 'max/min', side: 'right' }, tickvals: ticks.map(Math.log10), ticktext: ticks.map(String), thickness: 12 } },
      lin(b.rule_nseg, b.dD, 'rule L = (3D+d)/2', '#ffffff', { line: { color: '#ffffff', width: 2, dash: 'dash' } }),
      { x: [b.best.nseg], y: [b.best.dD], type: 'scatter', mode: 'markers', name: 'most neutral', marker: { symbol: 'star', size: 16, color: S(0), line: { color: '#000', width: 1 } } },
    ], { xaxis: { title: 'number of segments N', dtick: 1 }, yaxis: { title: 'core / OD, d/D' } });
  }
  // Fig 9 histograms
  function histPlot(id, h, scale, xt, lines) {
    const c = h.edges.slice(0, -1).map((e, i) => 0.5 * (e + h.edges[i + 1]) * scale);
    const shapes = lines.map((l) => ({ type: 'line', x0: l.x, x1: l.x, yref: 'paper', y0: 0, y1: 1, line: { color: l.c, width: 2, dash: l.d || 'solid' } }));
    const ann = lines.map((l, i) => ({ x: l.x, yref: 'paper', y: 0.98 - 0.1 * i, text: l.t, showarrow: false, xanchor: 'left', font: { size: 11, color: l.c } }));
    PR.plot(id, [{ x: c, y: h.counts, type: 'bar', marker: { color: S(1), opacity: 0.75 }, name: 'samples', hovertemplate: '%{x:.3g}: %{y}<extra></extra>' }],
      { xaxis: { title: xt }, yaxis: { title: 'count' }, shapes, annotations: ann, bargap: 0.02, showlegend: false });
  }
  histPlot('fig-hP', mc.hist_Pc, 1, 'peak P<sub>c</sub> [MPa]', [
    { x: mc.nominal.Pc_max / 1e6, t: ' nominal', c: PR.css('--ink-2') },
    { x: P.mu3s / 1e6, t: ' μ+3σ', c: S(2), d: 'dot' },
    { x: P.p99865 / 1e6, t: ' 99.865 %', c: S(0) }]);
  histPlot('fig-hI', mc.hist_I, 1, 'total impulse [N·s]', [
    { x: mc.nominal.I, t: ' nominal', c: PR.css('--ink-2') }]);
  // Fig 10 sensitivities
  {
    const lab = { a: 'a (coefficient)', n: 'n (exponent)', Dt: 'D<sub>t</sub>', rho: 'ρ<sub>p</sub>' };
    const o = mc.oat.params;
    const ys = params.map((p) => lab[p]);
    const xs = ['a', 'n', 'D<sub>t</sub>', 'ρ<sub>p</sub>'];
    PR.plot('fig-sens', [
      { type: 'bar', orientation: 'h', y: ys, x: params.map((p) => 100 * o[p].dPc_plus), name: '+1σ', marker: { color: S(0) } },
      { type: 'bar', orientation: 'h', y: ys, x: params.map((p) => 100 * o[p].dPc_minus), name: '−1σ', marker: { color: S(1) } },
      { type: 'scatter', mode: 'markers', y: ys.concat(ys), x: params.map((p) => 100 * o[p].linear_dlnPc).concat(params.map((p) => -100 * o[p].linear_dlnPc)),
        name: 'linearised', marker: { symbol: 'diamond', size: 9, color: PR.css('--ink'), line: { width: 0 } } },
    ], { barmode: 'overlay', xaxis: { title: 'change in peak P<sub>c</sub> [%]' }, yaxis: { automargin: true }, margin: { l: 110 } });
    PR.plot('fig-share', [
      { type: 'bar', x: xs, y: params.map((p) => 100 * src[p] * src[p] / varTot), name: 'a referenced to 1 MPa', marker: { color: S(0) } },
      { type: 'bar', x: xs, y: params.map((p) => 100 * srcA[p] * srcA[p] / varTotA), name: 'referenced to operating P<sub>c</sub>', marker: { color: S(1) } },
    ], { barmode: 'group', xaxis: { tickangle: 0 }, yaxis: { title: 'share of variance [%]', rangemode: 'tozero' } });
  }
  // Fig 11 sweep
  {
    const sw = mc.sweep;
    PR.plot('fig-sweep', [
      lin(sw.n, sw.MEOP_factor, 'MEOP / nominal', S(0), { mode: 'lines+markers' }),
      lin(sw.n, sw.cv.map((v) => 100 * v), 'σ/μ, a ref. 1 MPa', S(1), { yaxis: 'y2', mode: 'lines+markers', line: { color: S(1), width: 1.5, dash: 'dot' } }),
      lin(sw.n, sw.cv_anchor.map((v) => 100 * v), 'σ/μ, ref. operating P<sub>c</sub>', S(2), { yaxis: 'y2', mode: 'lines+markers', line: { color: S(2), width: 1.5, dash: 'dot' } }),
    ], { xaxis: { title: 'nominal burn-rate exponent n' }, yaxis: { title: 'MEOP / nominal peak' },
      yaxis2: { title: 'σ/μ of peak P<sub>c</sub> [%]', overlaying: 'y', side: 'right', showgrid: false, rangemode: 'tozero' }, margin: { r: 60 } });
  }

  // ------------------------------------------------------------------ interactive tool
  initTool();

  function initTool() {
    const presets = D.c.presets;
    const order = ['knsb', 'apcp', 'srb', 'geo_bates', 'geo_star', 'geo_finocyl', 'geo_moon'];
    const sel = $('c-preset');
    sel.innerHTML = order.map((k) => '<option value="' + k + '">' + presets[k].name + '</option>').join('') + '<option value="custom" disabled>custom</option>';
    const clone = (o) => JSON.parse(JSON.stringify(o));
    const GRID = D.c.js_grid;
    const ODmin = 20, ODmax = 4000;
    const odFromSlider = (v) => ODmin * Math.pow(ODmax / ODmin, v / 1000);
    const sliderFromOd = (od) => Math.round(1000 * Math.log(od / ODmin) / Math.log(ODmax / ODmin));

    const GEO = {
      bates: [{ k: 'core', l: 'Core diameter d/D', min: 0.08, max: 0.8, step: 0.005, def: 0.33 }],
      star: [{ k: 'points', l: 'Star points', min: 3, max: 14, step: 1, def: 6, int: true },
        { k: 'tip', l: 'Tip radius / R', min: 0.25, max: 0.92, step: 0.005, def: 0.6 },
        { k: 'valley', l: 'Valley radius / R', min: 0.05, max: 0.85, step: 0.005, def: 0.25 }],
      finocyl: [{ k: 'core', l: 'Core diameter d/D', min: 0.06, max: 0.6, step: 0.005, def: 0.18 },
        { k: 'slots', l: 'Slots', min: 2, max: 12, step: 1, def: 6, int: true },
        { k: 'slotw', l: 'Slot width / D', min: 0.01, max: 0.2, step: 0.0025, def: 0.06 },
        { k: 'slotr', l: 'Slot outer radius / R', min: 0.1, max: 0.92, step: 0.005, def: 0.52 }],
      moon: [{ k: 'core', l: 'Core diameter d/D', min: 0.08, max: 0.7, step: 0.005, def: 0.33 },
        { k: 'off', l: 'Core offset / R', min: 0.0, max: 0.85, step: 0.005, def: 0.45 }],
    };
    let current = null;          // spec being simulated
    let geoVals = {};
    function geoFromSpec(g, Rr) {
      const D2 = 2 * Rr;
      if (g.type === 'bates') return { core: g.core_d / D2 };
      if (g.type === 'star') return { points: g.points, tip: g.r_tip / Rr, valley: g.r_valley / Rr };
      if (g.type === 'finocyl') return { core: g.core_d / D2, slots: g.slots, slotw: g.slot_w / D2, slotr: g.slot_r / Rr };
      return { core: g.core_d / D2, off: g.offset / Rr };
    }
    function clampGeo(type, v) {
      if (type === 'star') v.valley = Math.min(v.valley, v.tip - 0.04);
      if (type === 'finocyl') v.slotr = Math.max(v.slotr, v.core + 0.03);
      if (type === 'moon') v.off = Math.min(v.off, Math.max(0, 1 - v.core - 0.06));
      return v;
    }
    function geomFromVals(type, v, Rr) {
      const D2 = 2 * Rr;
      if (type === 'bates') return { type, core_d: v.core * D2 };
      if (type === 'star') return { type, points: v.points, r_tip: v.tip * Rr, r_valley: v.valley * Rr };
      if (type === 'finocyl') return { type, core_d: v.core * D2, slots: v.slots, slot_w: v.slotw * D2, slot_r: v.slotr * Rr };
      return { type, core_d: v.core * D2, offset: v.off * Rr };
    }
    function buildGeoControls(type) {
      const host = $('geo-ctls');
      host.innerHTML = GEO[type].map((c) => '<label class="ctl"><span class="row"><span>' + c.l + '</span><output id="o-g-' + c.k + '"></output></span>' +
        '<input id="c-g-' + c.k + '" type="range" min="' + c.min + '" max="' + c.max + '" step="' + c.step + '"></label>').join('');
      GEO[type].forEach((c) => {
        const el = $('c-g-' + c.k);
        el.value = geoVals[c.k];
        el.addEventListener('input', () => { geoVals[c.k] = +el.value; clampGeo(type, geoVals); shapeChanged(); });
      });
    }
    function showOutputs() {
      const Rr = odFromSlider(+$('c-od').value) / 2e3;
      $('o-od').textContent = (2 * Rr * 1e3).toFixed(2 * Rr < 0.2 ? 1 : 0) + ' mm';
      $('o-ld').textContent = (+$('c-ld').value).toFixed(2) + '  (' + (2 * Rr * $('c-ld').value * 1e3).toFixed(0) + ' mm)';
      $('o-ns').textContent = $('c-ns').value;
      $('o-dt').textContent = (+$('c-dt').value).toFixed(3) + '  (' + (current ? current.Dt * 1e3 : 0).toFixed(2) + ' mm)';
      $('o-eps').textContent = (+$('c-eps').value).toFixed(1);
      const type = $('c-geom').value;
      GEO[type].forEach((c) => {
        const o = $('o-g-' + c.k), el = $('c-g-' + c.k);
        if (!o) return;
        el.value = geoVals[c.k];
        o.textContent = c.int ? String(geoVals[c.k]) : (+geoVals[c.k]).toFixed(3);
      });
    }
    function setControlsFromSpec(s) {
      const seg = s.segments[0];
      const Rr = s.R;
      $('c-prop').value = s.prop_key;
      $('c-geom').value = seg.geom.type;
      $('c-od').value = sliderFromOd(2 * Rr * 1e3);
      $('c-ld').value = seg.L / (2 * Rr);
      $('c-ns').value = seg.count || 1;
      $('c-ends').value = String(seg.ends || 0);
      $('c-dt').value = s.Dt / (2 * Rr);
      $('c-eps').value = s.eps;
      $('c-eros').checked = !!(s.erosive && s.erosive.on);
      $('c-thr').checked = !!s.throat_erosion;
      geoVals = geoFromSpec(seg.geom, Rr);
      buildGeoControls(seg.geom.type);
    }
    function loadPreset(k) {
      current = clone(presets[k]);
      current.grid = GRID;
      setControlsFromSpec(current);
      showOutputs();
      recompute();
    }
    // global controls patch the current spec (keeps multi-segment presets intact)
    function globalChanged() {
      const Rr = current.R;
      current.Dt = +$('c-dt').value * 2 * Rr;
      current.eps = +$('c-eps').value;
      current.prop_key = $('c-prop').value;
      current.prop = D.c.props[current.prop_key];
      current.erosive = { on: $('c-eros').checked, alpha: 5e-6, beta: 53.0 };
      current.throat_erosion = $('c-thr').checked ? 1e-4 : 0.0;
      $('c-preset').value = 'custom';
      showOutputs();
      schedule();
    }
    // shape controls rebuild a uniform motor from all controls
    function shapeChanged() {
      const Rr = odFromSlider(+$('c-od').value) / 2e3;
      const type = $('c-geom').value;
      clampGeo(type, geoVals);
      current = {
        name: 'custom', R: Rr, prop_key: $('c-prop').value, prop: D.c.props[$('c-prop').value],
        segments: [{ geom: geomFromVals(type, geoVals, Rr), L: +$('c-ld').value * 2 * Rr, ends: +$('c-ends').value, count: +$('c-ns').value }],
        Dt: +$('c-dt').value * 2 * Rr, eps: +$('c-eps').value, pa: 101325.0, grid: GRID,
        erosive: { on: $('c-eros').checked, alpha: 5e-6, beta: 53.0 }, throat_erosion: $('c-thr').checked ? 1e-4 : 0.0,
      };
      $('c-preset').value = 'custom';
      showOutputs();
      schedule();
    }
    sel.addEventListener('change', () => loadPreset(sel.value));
    ['c-dt', 'c-eps'].forEach((id) => $(id).addEventListener('input', globalChanged));
    ['c-prop', 'c-eros', 'c-thr'].forEach((id) => $(id).addEventListener('change', globalChanged));
    ['c-od', 'c-ld', 'c-ns'].forEach((id) => $(id).addEventListener('input', shapeChanged));
    $('c-ends').addEventListener('change', shapeChanged);
    $('c-geom').addEventListener('change', () => {
      const type = $('c-geom').value;
      geoVals = {};
      GEO[type].forEach((c) => { geoVals[c.k] = c.def; });
      buildGeoControls(type);
      shapeChanged();
    });

    // ---------- compute (tables one per frame, then ballistics)
    let timer = null, job = 0;
    const nextFrame = () => new Promise((res) => setTimeout(res, 0));  // yields to the UI; also runs in hidden tabs
    function schedule() { clearTimeout(timer); timer = setTimeout(recompute, 120); }
    function setStatus(t, warn) { const s = $('status'); s.textContent = t; s.className = 'status' + (warn ? ' warn' : ''); }

    let sim = null;
    async function recompute() {
      const my = ++job;
      const spec = current;
      setStatus('computing level set…');
      const pend = SRM.pendingTables(spec);
      for (let i = 0; i < pend.length; i++) {
        await nextFrame();
        if (my !== job) return;
        setStatus('fast marching ' + (i + 1) + '/' + pend.length + '…');
        SRM.grainTable(pend[i], spec.R, spec.grid.N, spec.grid.nw);
      }
      await nextFrame();
      if (my !== job) return;
      const t0 = performance.now();
      const m = new SRM.Motor(spec);
      const qs = m.quasiSteady();
      const warns = [];
      if (!(qs.Pc_max > 0.3e6)) {
        setStatus('equilibrium pressure below 0.3 MPa: enlarge the burning area or shrink the throat', true);
        sim = null; drawEmpty(); return;
      }
      const tr = m.transient(0.2, 200000);
      const met = m.metrics(tr, qs);
      const laws = m.prop.laws;
      if (met.Pc_max > laws[laws.length - 1][1] && laws.length > 1) warns.push('above the burn-rate fit range');
      if (met.Pc_max > 25e6) warns.push('very high pressure');
      if (tr.t.length >= 200000) warns.push('trace truncated');
      const ms = performance.now() - t0;
      setStatus('solved in ' + ms.toFixed(0) + ' ms' + (warns.length ? ' · ' + warns.join(', ') : ''), warns.length > 0);
      sim = prepareSim(m, tr, qs, met);
      renderReadout(met, m);
      startAnim();
    }

    // ---------- canvas
    const cv = $('xsec'), ctx = cv.getContext('2d');
    let pix = null;
    function resolveColor(c) {
      const t = document.createElement('canvas'); t.width = t.height = 1;
      const x = t.getContext('2d'); x.fillStyle = c; x.fillRect(0, 0, 1, 1);
      const d = x.getImageData(0, 0, 1, 1).data; return [d[0], d[1], d[2], d[3]];
    }
    let COL = null;
    function colors() {
      const bg = resolveColor(PR.css('--bg'));
      const acc = resolveColor(PR.css('--accent'));
      const prop = resolveColor(PR.css('--s5'));
      const mix = (a, c, f) => [0, 1, 2].map((i) => Math.round(a[i] * (1 - f) + c[i] * f)).concat([255]);
      const ink3c = resolveColor(PR.css('--ink-3'));
      COL = { port: mix(bg, bg, 0), burnt: mix(bg, acc, 0.22), prop: mix(bg, prop, 0.72), cas: mix(bg, ink3c, ink3c[3] / 255),
        front: PR.css('--accent'), init: PR.css('--ink-2') };
    }
    function prepareSim(m, tr, qs, met) {
      const tab = m.segs[0].tab;
      // down-sample the trace for plotting (~400 points uniform in time)
      const tEnd = tr.t[tr.t.length - 1];
      const idx = [];
      let j = 0;
      for (let q = 0; q <= 400; q++) {
        const tq = tEnd * q / 400;
        while (j < tr.t.length - 1 && tr.t[j] < tq) j++;
        if (!idx.length || idx[idx.length - 1] !== j) idx.push(j);
      }
      const pick = (a) => idx.map((i) => a[i]);
      const Fmax = met.F_max;
      const fs = Fmax > 5e5 ? 1e6 : Fmax > 5e3 ? 1e3 : 1;
      const fu = fs === 1e6 ? 'MN' : fs === 1e3 ? 'kN' : 'N';
      return { m, tab, tr, qs, met, tEnd, pt: pick(tr.t), pF: pick(tr.F).map((v) => v / fs), pP: pick(tr.Pc).map((v) => v / 1e6), fu,
        K: m.K };
    }
    function buildPixels() {
      if (!sim) return;
      const rect = cv.getBoundingClientRect();
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const Sz = Math.max(120, Math.min(720, Math.round(rect.width * dpr)));
      cv.width = Sz; cv.height = Sz;
      const tab = sim.tab, g = tab.grid, N = g.N, Rr = tab.R, Rv = Rr * 1.06;
      const field = new Float32Array(Sz * Sz), where = new Uint8Array(Sz * Sz); // 0 outside, 1 case ring, 2 grain
      for (let py = 0; py < Sz; py++) {
        const y = Rv - (py + 0.5) / Sz * 2 * Rv;
        for (let px = 0; px < Sz; px++) {
          const x = -Rv + (px + 0.5) / Sz * 2 * Rv;
          const r = Math.hypot(x, y), k = py * Sz + px;
          if (r > Rv) { where[k] = 0; continue; }
          if (r > Rr) { where[k] = 1; continue; }
          where[k] = 2;
          const fx = (x - g.x[0]) / g.h, fy = (y - g.x[0]) / g.h;
          const jx = Math.min(N - 2, Math.max(0, Math.floor(fx))), iy = Math.min(N - 2, Math.max(0, Math.floor(fy)));
          const tx = fx - jx, ty = fy - iy, q = iy * N + jx, T = tab.T;
          field[k] = T[q] * (1 - tx) * (1 - ty) + T[q + 1] * tx * (1 - ty) + T[q + N] * (1 - tx) * ty + T[q + N + 1] * tx * ty;
        }
      }
      pix = { Sz, field, where, Rv, img: ctx.createImageData(Sz, Sz) };
    }
    function drawEmpty() { ctx.clearRect(0, 0, cv.width, cv.height); $('readout').innerHTML = ''; }
    function drawSection(w) {
      if (!sim || !pix) return;
      if (!COL) colors();
      const { Sz, field, where, img, Rv } = pix;
      const d = img.data;
      for (let k = 0; k < Sz * Sz; k++) {
        const o = 4 * k, wh = where[k];
        let c;
        if (wh === 0) { d[o + 3] = 0; continue; }
        if (wh === 1) c = COL.cas;
        else { const T = field[k]; c = T < 0 ? COL.port : T < w ? COL.burnt : COL.prop; }
        d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2]; d[o + 3] = 255;
      }
      ctx.putImageData(img, 0, 0);
      const s = Sz / (2 * Rv), X = (x) => (x + Rv) * s, Y = (y) => (Rv - y) * s, Rr = sim.tab.R;
      const stroke = (level, color, width, dash) => {
        ctx.beginPath();
        SRM.contour(sim.tab.T, sim.tab.grid, level, (x1, y1, x2, y2) => {
          if (Math.hypot(0.5 * (x1 + x2), 0.5 * (y1 + y2)) > Rr) return;
          ctx.moveTo(X(x1), Y(y1)); ctx.lineTo(X(x2), Y(y2));
        });
        ctx.setLineDash(dash || []); ctx.lineWidth = width; ctx.strokeStyle = color; ctx.stroke(); ctx.setLineDash([]);
      };
      const dpr = Sz / Math.max(1, cv.getBoundingClientRect().width);
      stroke(0, COL.init, 1 * dpr, [4 * dpr, 3 * dpr]);
      if (w > 0 && w < sim.tab.w_max) stroke(w, COL.front, 2.2 * dpr);
    }

    // ---------- plot + animation
    let playing = true, animT = 0, lastTs = null, raf = null, frame = 0;
    const DUR = 5.0;
    function plotBase() {
      const s = sim;
      PR.plot('tool-plot', [
        lin([], [], 'thrust', S(0)),
        lin([], [], 'P<sub>c</sub>', S(1), { yaxis: 'y2' }),
        lin(s.pt, s.pF, 'thrust (full)', S(0), { opacity: 0.22, showlegend: false, hoverinfo: 'skip' }),
        lin(s.pt, s.pP, 'P<sub>c</sub> (full)', S(1), { yaxis: 'y2', opacity: 0.22, showlegend: false, hoverinfo: 'skip' }),
      ], { xaxis: { title: 'time [s]', range: [0, s.tEnd] },
        yaxis: { title: 'thrust [' + s.fu + ']', range: [0, 1.08 * Math.max(...s.pF)] },
        yaxis2: { title: 'P<sub>c</sub> [MPa]', overlaying: 'y', side: 'right', showgrid: false, range: [0, 1.08 * Math.max(...s.pP)] },
        margin: { r: 56, t: 30 } });
    }
    function wAt(t) {
      const tr = sim.tr;
      let lo = 0, hi = tr.t.length - 1;
      if (t >= tr.t[hi]) return tr.w[hi];
      while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (tr.t[mid] <= t) lo = mid; else hi = mid; }
      const f = (t - tr.t[lo]) / (tr.t[hi] - tr.t[lo]);
      return tr.w[lo] + f * (tr.w[hi] - tr.w[lo]);
    }
    function showTime(tn) {
      if (!sim) return;
      const t = tn * sim.tEnd;
      const w = wAt(t);
      drawSection(w);
      let n = 0;
      while (n < sim.pt.length && sim.pt[n] <= t) n++;
      if (window.Plotly && (frame++ % 2 === 0 || tn >= 1)) {
        Plotly.restyle('tool-plot', { x: [sim.pt.slice(0, n), sim.pt.slice(0, n)], y: [sim.pF.slice(0, n), sim.pP.slice(0, n)] }, [0, 1]);
      }
      $('c-scrub').value = Math.round(1000 * tn);
      $('xsec-cap').textContent = 't = ' + t.toFixed(t < 10 ? 2 : 1) + ' s · web burned ' + (w * 1e3).toFixed(w < 0.1 ? 1 : 0) + ' mm' +
        (sim.K > 1 ? ' · segment 1 of ' + sim.K + ' shown' : '');
    }
    function tick(ts) {
      raf = null;
      if (!playing || !sim) return;
      if (lastTs === null) lastTs = ts;
      animT = Math.min(1, animT + (ts - lastTs) / 1000 / DUR);
      lastTs = ts;
      showTime(animT);
      if (animT < 1) raf = requestAnimationFrame(tick);
      else { playing = false; $('b-play').textContent = 'Play'; }
    }
    function startAnim() {
      buildPixels();
      plotBase();
      animT = 0; lastTs = null; playing = true; frame = 0;
      $('b-play').textContent = 'Pause';
      showTime(0);
      if (raf) cancelAnimationFrame(raf);
      raf = requestAnimationFrame(tick);
    }
    $('b-play').addEventListener('click', () => {
      if (!sim) return;
      playing = !playing;
      if (playing && animT >= 1) animT = 0;
      $('b-play').textContent = playing ? 'Pause' : 'Play';
      lastTs = null;
      if (playing && !raf) raf = requestAnimationFrame(tick);
    });
    $('b-replay').addEventListener('click', () => { if (sim) { animT = 0; lastTs = null; playing = true; $('b-play').textContent = 'Pause'; if (!raf) raf = requestAnimationFrame(tick); } });
    $('c-scrub').addEventListener('input', () => {
      if (!sim) return;
      playing = false; $('b-play').textContent = 'Play';
      animT = +$('c-scrub').value / 1000; frame = 0; showTime(animT);
    });
    let rz = null;
    window.addEventListener('resize', () => { clearTimeout(rz); rz = setTimeout(() => { buildPixels(); showTime(animT); }, 150); });
    window.matchMedia('(prefers-color-scheme: light)').addEventListener('change', () => { COL = null; showTime(animT); });

    function renderReadout(met, m) {
      const st = (v, l) => '<div class="stat"><div class="v">' + v + '</div><div class="l">' + l + '</div></div>';
      const I = met.I_total;
      const Itxt = I >= 1e7 ? sci(I, 2) : Math.round(I).toLocaleString('en-US');
      $('readout').innerHTML = [
        st(Itxt + '<small>N·s</small>', 'total impulse'),
        st(met.class.length === 1 ? met.designation : '—', 'motor class (letter + avg. thrust)'),
        st((met.Pc_max / 1e6).toFixed(2) + '<small>MPa</small>', 'peak chamber pressure'),
        st(met.Kn_initial.toFixed(0) + ' / ' + met.Kn_max.toFixed(0), 'K<sub>n</sub> initial / max'),
        st(met.t_burn.toFixed(met.t_burn < 10 ? 2 : 1) + '<small>s</small>', 'burn time (10 % of peak)'),
        st(met.Isp.toFixed(1) + '<small>s</small>', 'I<sub>sp</sub> (ideal, sea level)'),
        st(met.m_prop < 1000 ? met.m_prop.toFixed(met.m_prop < 10 ? 3 : 1) + '<small>kg</small>' : (met.m_prop / 1e3).toFixed(0) + '<small>t</small>', 'propellant mass'),
        st(met.neutrality_pc.toFixed(2), 'P<sub>c</sub> max/min over web'),
      ].join('');
    }

    // live cross-check against the Python reference (KNSB preset, same grid)
    {
      const spec = clone(presets.knsb); spec.grid = GRID;
      const m = new SRM.Motor(spec);
      const met = m.metrics(m.transient(), m.quasiSteady());
      const ref = D.v.js_ref.knsb;
      const dI = Math.abs(met.I_total / ref.I_total - 1), dP = Math.abs(met.Pc_max / ref.Pc_max - 1);
      $('live-check').innerHTML = 'Live check in this browser, KNSB preset: total impulse ' + met.I_total.toFixed(2) + ' N·s vs Python ' +
        ref.I_total.toFixed(2) + ' N·s (rel. diff ' + sci(Math.max(dI, 1e-16), 1) + '); peak P<sub>c</sub> rel. diff ' + sci(Math.max(dP, 1e-16), 1) + '.';
    }
    loadPreset('knsb');
  }
})().catch((e) => {
  console.error(e);
  const s = document.getElementById('status');
  if (s) { s.textContent = 'error: ' + e.message; s.className = 'status warn'; }
});
