# nozzlecfd: quasi-1D Euler CFD of rocket nozzles

Propulsion Research Series, No. 1. This is the tool behind `../index.html`.

## What it does

- `nozzlecfd/solver.py`: a finite-volume solver for the quasi-1D Euler equations in a convergent-divergent nozzle.
  - The scheme is HLLC with Davis/Einfeldt wave-speed bounds, MUSCL reconstruction of (ρ, u, p) with a van Albada limiter (minmod, van Leer and first order are also available), and three-stage SSP-RK3 time marching. Local time stepping is optional.
  - The area source term is written in a well-balanced form: gas at rest is an exact discrete steady state.
  - Inflow is a subsonic reservoir with given p0 and T0.
  - Outflow is a back-pressure boundary that switches between supersonic extrapolation and imposed p_b, using the normal-shock pressure at the exit Mach number.
  - The gas is calorically perfect with user-set γ. Everything is non-dimensional (p0 = ρ0 = 1).
- `nozzlecfd/exact.py`: the exact quasi-1D solution. It covers the isentropic area-Mach relation on both branches, the normal-shock jump relations, the critical back pressures (p_sub, p_NSE, p_sup), the shock position for a given back pressure and its inverse, the choked mass flux, and the ideal C_F.
- `nozzlecfd/atmosphere.py`: the 1976 US Standard Atmosphere from 0 to 51 km (five layers, geometric to geopotential conversion).
- `nozzlecfd/performance.py`: C_F from the CFD exit flux, the optimal-expansion altitude, the Summerfield (p_e/p_a < 0.4) and Schmucker separation-risk flags, and the maximum sea-level-safe ε and minimum p_c.
- `nozzlecfd/geometry.py`: the Anderson CD nozzle A = 1 + 2.2(x-1.5)², its shock-test variant, and a smooth bell-like area law for a given ε.
- `../js/nozzle-solver.js`: a line-by-line JavaScript port of the solver, exact solution and atmosphere, used by the page. `crosscheck_js.cjs` runs it under Node so `run_study.py` can compare it with Python.

## Run

```
/opt/homebrew/bin/python3 run_study.py            # about 2.5 min; writes ../data/*.json and CSVs
for t in tests/test_*.py; do /opt/homebrew/bin/python3 $t; done   # or: python3 -m pytest tests
```

The tests also work under pytest. Each test file can be run on its own with plain `python3`.

## Results summary

These numbers come from `../data/summary.json`, `verification.json`, `engines.json` and `crosscheck.json`. All engine inputs are approximate public values, and γ = 1.2.

| Check / case | Result |
|---|---|
| Isentropic Anderson nozzle, L1 Mach error, N = 50 to 800 | observed order 1.88, 1.94, 1.97, 1.99 (first-order variant: 0.97 to 0.99) |
| Richardson extrapolation of M_e, ε = 16, γ = 1.2 | order 1.86, 1.95; extrapolated M_e matches exact to 1e-6 |
| Shock in divergent section (exact x_s = 2.25) | error 0.13–0.17 Δx on every grid, first-order convergence |
| Mass conservation | face-to-face ṁ spread ≤ 2e-12; ṁ/ṁ* error 5e-5 at N = 200 (second order) |
| JS port vs Python (3 cases incl. a shock) | identical iteration counts, max \|ΔM\| ≈ 2e-14 |

| Engine (p_c, ε) | p_e (kPa) | p_e/p_a SL | Summerfield at SL | p_e = p_a at | C_F SL / vac | p_NSE / p_SL |
|---|---|---|---|---|---|---|
| Merlin 1D-class (9.7 MPa, 16) | 65.7 | 0.65 | clear | 3.5 km | 1.630 / 1.797 | 9.1 |
| RS-25-class (20.6 MPa, 69) | 21.6 | 0.21 | risk (clears 5.0 km) | 11.3 km | 1.587 / 1.927 | 4.9 |
| Raptor-class (30 MPa, ≈40) | 62.6 | 0.62 | clear | 3.9 km | 1.749 / 1.884 | 12.1 |

## Limitations

- The model is inviscid and quasi-1D, so it **cannot model flow separation**. Separation risk is only flagged, using empirical criteria (Summerfield, Schmucker). Real free-shock and restricted-shock separation is 2D, viscous and unsteady.
- The exhaust has a single frozen γ. The RS-25-class verdict is sensitive to γ (p_e/p_a at sea level ranges from 0.27 to 0.14 for γ = 1.15 to 1.3).
- There are no divergence, boundary-layer or chamber losses. The engine inputs are approximate, and the Raptor figures are particularly uncertain.
- Numerics:
  - The shock position is first-order accurate, as for any shock-capturing scheme.
  - Within about 1 % of p_NSE, the choice between "shock just inside the exit" and "supersonic exit" is resolved only to discretisation accuracy.
  - From a supersonic initial guess (`init="ramp"`), a shock that belongs just inside the exit can lock onto the outflow boundary. The default start from rest avoids this.
