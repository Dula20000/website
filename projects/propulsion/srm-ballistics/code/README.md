# SRM internal ballistics with level-set grain regression

Propulsion Research Series, No. 3. This is a solid rocket motor design tool with two parts:

- **Grain regression.** A 2-D second-order fast-marching (level-set) method regresses any port cross-section: BATES, N-point star, finocyl-style slotted core, or moon burner. The burning perimeter and port area come from marching-squares contours clipped to the case. Segments can have 0, 1 or 2 burning ends.
- **Lumped chamber ballistics.** Burn rate follows Saint-Robert/Vieille, piecewise if needed. The tool solves both the quasi-steady equilibrium and the transient mass-balance ODE (RK4 from ambient, including filling and blow-down), with an ideal nozzle `C_F` for any expansion ratio and ambient pressure. Throat erosion and Lenoir–Robillard erosive burning are optional switches, off by default.
- **Metrics.** Total impulse, average and maximum thrust, burn time (10 % of peak), Isp, NAR/TRA class letter, `K_n(t)`, neutrality (`P_c` max/min over the web burn) and thrust-shape ratio.
- **Monte Carlo.** Vectorised equilibrium model sampling `a`, `n`, throat diameter and density. It reports the 99.865th-percentile (3σ-equivalent) MEOP, bootstrap intervals, one-at-a-time and regression sensitivities, and an exponent sweep.

`../srm.js` is a line-by-line JavaScript port used by the interactive page. `tests/js_crosscheck.js` runs it under node and compares it with the Python results.

## Run

```
python3 run_study.py                 # regenerates ../data/*.json and thrust_*.csv, runs the node cross-check
python3 tests/test_geometry.py       # or: python3 -m pytest tests
python3 tests/test_ballistics.py
python3 tests/test_montecarlo.py
node tests/js_crosscheck.js          # JS vs Python (needs data from run_study.py)
```

Requires numpy only. Node is optional.

## Layout

| file | contents |
|---|---|
| `srm/geometry.py` | port primitives (signed distance), FMM, marching squares, `GrainTable`, analytic references |
| `srm/propellant.py` | burn-rate laws, equilibrium pressure, ideal nozzle |
| `srm/ballistics.py` | `Motor`: burning area, quasi-steady and transient solutions, metrics |
| `srm/montecarlo.py` | vectorised equilibrium ensemble, Monte Carlo, sensitivities |
| `srm/cases.py` | KNSB / APCP / PBAN propellants and the case-study motors (approximate public values) |
| `run_study.py` | every result used by the page |

## Summary of results (from `../data/`)

| item | value |
|---|---|
| Circular-port perimeter vs `2π(r0+w)` | max rel. error 9.0e-5 at N = 161; observed order 2.00 |
| Star perimeter vs exact offset (FMM) | mean rel. error 1.6e-3 (N = 201), 2.6e-4 (N = 801); first order (corners) |
| BATES `A_b` vs closed form | max rel. error 2.8e-4 |
| Transient mass balance | ≤ 2.4e-4 |
| JS vs Python | max rel. difference 4e-12 |
| KNSB 54 mm BATES (4 grains) | J518, 970 N·s, 4.35 MPa peak, Isp 142.9 s (ideal) |
| APCP 98 mm BATES (4 grains) | N2685, 11 859 N·s, 5.09 MPa peak, Isp 235.0 s (ideal) |
| SRB-class star + stepped bore | qualitative two-hump thrust bucket; ~1.1e9 N·s |
| BATES trade (L_total = 0.52 m) | most neutral: 3 segments, d/D = 0.55 (P_c max/min 1.127); rule (3D+d)/2 agrees for 13 of 15 core ratios |
| Monte Carlo (20 000 samples) | MEOP (99.865 %) 6.41 MPa = 1.255 × nominal; n contributes ~49 % of the variance |

## Limitations

- 0-D chamber: no port flow or axial pressure drop. Erosive burning is a per-segment proxy.
- 2-D regression per segment with a uniform burn rate. Tapers are approximated by steps.
- Perimeter accuracy is only first order for cornered ports. `P(w)` is rescaled so that `A_p(0) + ∫P dw = πR²` (mass conservation). The raw closure error is reported.
- Ideal, frozen thermochemistry. There are no nozzle or two-phase losses and no separation model: negative over-expanded thrust is clipped to zero.
- Propellant inputs are approximate public values (Nakka's KNSB fit is recalled). Monte Carlo tolerances are assumed and treated as independent.
