// Cross-check: run the browser solver (../../srm.js) under node on the KNSB and APCP presets at the
// tool's grid and compare with the Python reference written by run_study.py (verification.json -> js_ref).
// Writes ../../data/js_crosscheck.json.  Usage: node code/tests/js_crosscheck.js
'use strict';
const fs = require('fs');
const path = require('path');
const SRM = require(path.join(__dirname, '..', '..', 'srm.js'));

const DATA = path.join(__dirname, '..', '..', 'data');
const cases = JSON.parse(fs.readFileSync(path.join(DATA, 'cases.json'), 'utf8'));
const ver = JSON.parse(fs.readFileSync(path.join(DATA, 'verification.json'), 'utf8'));

const keys = ['I_total', 'Pc_max', 'F_max', 't_burn', 'Isp', 'Kn_initial', 'neutrality_pc', 'qs_I', 'qs_Pc_max'];
const out = { grid: cases.js_grid, cases: {}, max_rel_diff: 0 };
for (const name of ['knsb', 'apcp']) {
  const spec = JSON.parse(JSON.stringify(cases.presets[name]));
  spec.grid = cases.js_grid;
  const t0 = Date.now();
  const m = new SRM.Motor(spec);
  const qs = m.quasiSteady();
  const tr = m.transient();
  const met = m.metrics(tr, qs);
  const ms = Date.now() - t0;
  const ref = ver.js_ref[name];
  const row = { ms, js: {}, py: {}, rel: {} };
  const tab = m.segs[0].tab;
  const extra = { P0: tab.P[0], P_mid: tab.P[Math.floor(cases.js_grid.nw / 2)], w_max: tab.w_max };
  for (const k of keys.concat(Object.keys(extra))) {
    const js = k in extra ? extra[k] : met[k];
    const py = ref[k];
    row.js[k] = js; row.py[k] = py;
    row.rel[k] = Math.abs(js / py - 1);
    out.max_rel_diff = Math.max(out.max_rel_diff, row.rel[k]);
  }
  out.cases[name] = row;
  console.log(`  JS ${name}: ${ms} ms, I=${met.I_total.toFixed(2)} (py ${ref.I_total.toFixed(2)}), ` +
    `Pc_max=${(met.Pc_max / 1e6).toFixed(4)} MPa (py ${(ref.Pc_max / 1e6).toFixed(4)})`);
}
fs.writeFileSync(path.join(DATA, 'js_crosscheck.json'), JSON.stringify(out));
console.log('  JS max relative difference:', out.max_rel_diff.toExponential(2));
if (out.max_rel_diff > 1e-6) process.exit(1);
