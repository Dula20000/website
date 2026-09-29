// Cross-check: run the JavaScript port against the Python results in ../data.
// Usage: node check.mjs <data-dir>   (prints a JSON summary on stdout)
import { createRequire } from 'module';
import { readFileSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';

const require = createRequire(import.meta.url);
const TF = require(join(dirname(fileURLToPath(import.meta.url)), 'tfcycle.js'));
const dataDir = process.argv[2] || join(dirname(fileURLToPath(import.meta.url)), '..', '..', 'data');
const cases = JSON.parse(readFileSync(join(dataDir, 'cases.json'), 'utf8')).cases;
const od = JSON.parse(readFileSync(join(dataDir, 'offdesign.json'), 'utf8'));

const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);
let worstDesign = 0, worstStation = 0, worstOff = 0, nOff = 0;
const details = [];
for (const c of cases) {
  const r = TF.design(c.inputs), p = c.result;
  let w = 0;
  for (const k of ['Fs', 'S', 'f', 'eta_th', 'eta_p', 'V9', 'V19']) w = Math.max(w, rel(r[k], p[k]));
  let ws = 0;
  for (const s of Object.keys(p.stations)) ws = Math.max(ws, rel(r.stations[s].Tt, p.stations[s].Tt), rel(r.stations[s].Pt, p.stations[s].Pt));
  worstDesign = Math.max(worstDesign, w); worstStation = Math.max(worstStation, ws);
  details.push({ key: c.key, design: w, stations: ws, S_lb_js: r.S_lb, S_lb_py: p.S_lb });
}
for (const eng of od.engines) {
  const c = cases.find((x) => x.key === eng.key);
  const E = TF.Engine(c.inputs, eng.F_des);
  const sls = E.operate(0, 0, TF.schedule(E.inp, 0, 0, od.DTmax));
  let w = 0;
  for (const L of eng.lapse_alt.filter((x) => x.M0 === 0.78 || x.M0 === 0)) {
    L.alt.forEach((h, j) => {
      const o = E.operate(L.M0, h, L.Tt4[j]);
      w = Math.max(w, rel(o.F / sls.F, L.F_rel[j]), rel(o.S_lb, L.S_lb[j])); nOff++;
    });
  }
  const H = eng.throttle;
  H.Tt4.forEach((T4, j) => {
    const o = E.operate(c.inputs.M0, c.inputs.alt, T4);
    w = Math.max(w, rel(o.F / E.F_des, H.F_rel[j]), rel(o.S_lb, H.S_lb[j])); nOff++;
  });
  worstOff = Math.max(worstOff, w);
  details.push({ key: eng.key, offdesign: w });
}
console.log(JSON.stringify({ worst_design: worstDesign, worst_station: worstStation, worst_offdesign: worstOff,
  n_offdesign_points: nOff, n_cases: cases.length, details }));
