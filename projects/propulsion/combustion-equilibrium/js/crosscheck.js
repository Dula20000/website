// Node cross-check: runs the JS port on the cases listed in argv[2] (JSON file)
// and prints the results as JSON.  Called by code/run_study.py.
'use strict';
const fs = require('fs');
const path = require('path');
const { create } = require('./eqsolver.js');
const thermo = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'data', 'thermo.json'), 'utf8'));
const S = create(thermo);
const cases = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const t0 = Date.now();
const out = cases.map((c) => {
  const r = S.performance(c.fuel, c.ox, c.of, c.Pc, c.eps, { frozen: c.frozen });
  const keep = {};
  for (const k of ['Tc', 'M_c', 'gamma_s_c', 'cstar', 'CF_vac', 'Isp_vac', 'Isp_sl', 'Te', 'Pe_kPa']) keep[k] = r[k];
  return Object.assign({ id: c.id }, keep);
});
process.stdout.write(JSON.stringify({ results: out, ms: Date.now() - t0 }));
