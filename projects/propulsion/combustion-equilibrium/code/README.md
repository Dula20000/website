# eqsolver: chemical equilibrium and rocket performance (Propulsion Research No. 2)

This is a small CEA-style code, written from scratch, for the C–H–O(–N–Ar) ideal-gas system.

- `eqsolver/thermo.py`: NASA 7-coefficient polynomials for H2, O2, H, O, OH, H2O, HO2, H2O2, CO, CO2, CH4, CH3, HCO and CH2O, plus N2, NO and Ar. The coefficients were typed in from GRI-Mech 3.0. OH uses the later Burcat/ATcT fit. Handbook reference values at 298 K are included for testing.
- `eqsolver/equilibrium.py`: Gibbs minimisation with the NASA RP-1311 reduced Newton iteration (element potentials, Δln n, Δln T). It handles TP, HP and SP problems, uses RP-1311 step control, and computes equilibrium and frozen derivatives (cp, γ_s, sound speed).
- `eqsolver/rocket.py`: propellants (LOX, LH2, LCH4, RP-1 surrogate CH1.9423) with assigned enthalpies, and infinite-area-combustor rocket performance. The throat is found where u = a, and the exit from the area ratio. Frozen and shifting expansion are both supported, giving c*, C_F, and vacuum and sea-level Isp.
- `run_study.py`: regenerates everything in `../data/` (`thermo.json`, `cases.json`, `sweeps.json`, `verification.json`, `sweeps_of.csv`). If `node` is available it also runs the JS port (`../js/eqsolver.js`) and records the cross-check.
- `tests/`: verification tests. They run with `python3 -m pytest tests` or directly with `python3 tests/test_thermo.py`, `python3 tests/test_equilibrium.py` and `python3 tests/test_rocket.py`.

## Run

```
python3 run_study.py          # ~2 s
python3 tests/test_thermo.py && python3 tests/test_equilibrium.py && python3 tests/test_rocket.py
```

## Verification summary (see `data/verification.json`)

| Check | Result |
|---|---|
| Δh_f and s° at 298 K, all 17 species | within ±0.05 kJ/mol (major) / ±1.5 kJ/mol (minor) |
| cp, h, s continuity at 1000 K | < 2e-6 relative |
| Element conservation, HP energy balance | ~1e-14 relative |
| Law of mass action, 6 reactions | ~1e-14 in ln K |
| Random element-preserving perturbations | G always increases, ΔG ∝ h² |
| SciPy SLSQP direct minimisation | same minimum, major species within 5e-6 |
| Argon nozzle vs closed-form γ = 5/3 ideal rocket | c*, Pe/Pc, C_F to ~5e-13 |
| H2–O2 / CH4–air stoich. T_ad, 1 atm | 3074 K / 2224 K (textbook ≈3080 / ≈2226 K) |
| JS port vs Python, 8 engine cases | ~3e-13 relative |

## Case-study results (ideal, shifting; `data/cases.json`)

| Engine class (public inputs, approx.) | Tc (K) | c* (m/s) | Isp vac shift / frozen (s) | delivered vac ≈ | η_vac |
|---|---|---|---|---|---|
| Merlin 1D (LOX/RP-1, 9.7 MPa, ε 16, O/F 2.36) | 3644 | 1811 | 336.5 / 321.2 | 311 | 0.92 |
| RS-25 (LOX/LH2, 20.6 MPa, ε 69, O/F 6.0) | 3598 | 2323 | 462.9 / 445.5 | 452 | 0.98 |
| Raptor SL (LOX/CH4, 30 MPa, ε 40, O/F 3.6) | 3763 | 1851 | 371.1 / 348.4 | 347 | 0.94 |
| RD-180 (LOX/RP-1, 26.7 MPa, ε 36.9, O/F 2.72) | 3890 | 1810 | 358.6 / 336.5 | 338 | 0.94 |

Peak shifting vacuum Isp at 10 MPa and ε = 40 falls at O/F ≈ 4.76 (LH2), 3.44 (CH4) and 2.78 (RP-1).

## Limitations

- The thermo fits are nominally valid to 3500 K. Several chambers reach 3600–3900 K, so the polynomials are extrapolated there.
- The model is gas-phase only: no condensed carbon, no ions. At the O/F ranges studied, O/C ≥ 1.4.
- RP-1 is a single-formula surrogate. The propellant enthalpies are approximate recollections of the CEA `thermo.inp` values.
- The nozzle is ideal: infinite-area combustor, 1-D isentropic flow, no kinetics, no divergence or boundary-layer losses, and no sea-level flow separation.
- Engine inputs and delivered Isp values are rounded public figures.
