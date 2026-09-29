// Runs the browser solver (../js/nozzle-solver.js) under node and prints JSON for run_study.py to compare.
// Usage: node crosscheck_js.cjs '<json list of cases>'
const path = require('path');
const S = require(path.join(__dirname, '..', 'js', 'nozzle-solver.js'));
const cases = JSON.parse(process.argv[2]);
const out = cases.map((c) => {
  const geom = c.geom === 'anderson' ? S.andersonNozzle(c.eps) : S.bellNozzle(c.eps);
  const s = new S.Solver(geom, c.N, c.gamma, c.pb, { cfl: c.cfl, localDt: true, init: 'rest' });
  const t0 = Date.now();
  const conv = s.run(c.max_iter, c.tol);
  const st = s.state();
  return { name: c.name, iters: s.iter, converged: conv, ms: Date.now() - t0,
           M: Array.from(st.M), p: Array.from(st.p), mdot: st.mdotOut, mom_exit: st.momExit,
           x_shock: S.shockLocation(st.x, st.M, geom.x_throat) };
});
process.stdout.write(JSON.stringify(out));
