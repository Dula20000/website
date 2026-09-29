# tfcycle: turbofan parametric and off-design cycle analysis

Propulsion Research Series, No. 4: *Why turbofans keep getting bigger*.

This is a small Python package for the cycle analysis of a **two-spool, separate-exhaust turbofan**. It follows the parametric
and off-design analysis in Mattingly's *Elements of Propulsion*, uses a variable-property gas model, and includes a simple
turbine cooling-air bleed. A line-by-line JavaScript port runs the interactive tool on the research page.

## What it does

| Module | Contents |
|---|---|
| `tfcycle/atmosphere.py` | 1976 US Standard Atmosphere, all layers to 84.852 km geopotential |
| `tfcycle/thermo.py` | gas models. `varcp` (default): thermally perfect air and Jet-A (CH1.92) combustion products from NASA 7-coefficient polynomials (GRI-Mech 3.0 data for N2, O2, Ar, CO2, H2O), sensible enthalpy from 298.15 K, entropy function s°(T), T(h), T(s°), choking, mass flux. `cpg`: calorically perfect gas (Mattingly two-gas / ideal fallback) |
| `tfcycle/ondesign.py` | parametric cycle: polytropic fan/LPC/HPC/HPT/LPT via s°; inlet (MIL-E-5007), burner and nozzle pressure ratios; burner enthalpy balance `h_air(Tt3) + f_b·ηb·hPR = (1+f_b)·h_prod(Tt4, f_b)`; cooling bleed ε taken at station 3 and remixed after the HPT (4.4 → 4.5); convergent or ideally expanded nozzles; specific thrust, TSFC, efficiencies; closed-form ideal turbofan; TSFC-optimal fan pressure ratio; exact marginal optimum condition; T-s path |
| `tfcycle/offdesign.py` | fixed-geometry off-design matching: choked HPT and LPT nozzle throats (LPT throat after coolant mixing), fixed convergent core and fan nozzles, constant adiabatic efficiencies and cooling fraction, LPC work tied to fan work, Brent solves |
| `tfcycle/installation.py` | fan diameter from fan-face Mach number and hub/tip ratio; nacelle skin-friction drag plus weight drag scaled with D²; installed TSFC |
| `js/tfcycle.js` | JavaScript port of all of the above (browser and Node) |
| `js/check.mjs` | runs the JS port against the Python results in `../data` |

## Run

```bash
python3 run_study.py            # regenerates ../data/*.json and CSV (about 40 s); runs the JS cross-check if node is installed
python3 tests/test_cycle.py     # 13 verification tests (or: python3 -m pytest tests)
node js/check.mjs ../data       # JS vs Python comparison on its own
```

Requires Python 3 with scipy (`brentq`). No other dependencies.

## Verification (see `tests/test_cycle.py`, `../data/verification.json`)

| Check | Worst relative error |
|---|---|
| `thermo='cpg'` with all losses off vs closed-form ideal turbofan (5 input sets, M 0–2, BPR 0–12) | 1.0e-15 |
| Ideal turbojet limit vs independent closed form | 1.3e-16 |
| Whole-engine first-law balance (sensible enthalpies); HP and LP spool power balances | 1.3e-15 |
| Species cp, s°, ΔH(1000 K) vs JANAF tables; formation enthalpies of CO2, H2O | 0.20 %; < 0.01 % |
| T(h), T(s°) round trip; continuity at 1000 K | 1.8e-15 |
| Air Δh across each burner vs an independent cubic cp(T) fit (Çengel) | 0.09 % |
| Numerical optimum π_f vs analytic ideal optimum τ_f* | 1.0e-8 (golden-section tolerance on a flat minimum) |
| Numerical real-cycle optimum vs exact marginal condition V19/V9 = η_mL[1−x19(1−e_f)]/[1+x9(1/e_tL−1)] (variable cp, with cooling) | 1.8e-7 |
| Mattingly's approximate condition V19/V9 ≈ η_f η_tL | 2.3–6.2 % (it is an approximation) |
| Off-design solver at the design point vs on-design | 6.8e-14 |
| Off-design matching residuals (80 points) | 3.0e-13 |
| 1976 atmosphere layer-base pressures | 1.8e-7 |
| JS port vs Python (4 cases at design; 405 off-design points) | 1.6e-15; 5.4e-14 |

## Results summary (representative public inputs; Tt4, efficiencies and ε = 0.12 assumed)

| Engine class (cruise point) | BPR | π_f | OPR | Tt4 [K] | F/ṁ0 [N·s/kg] | TSFC [lb/(lbf·h)] (Tt4 ∓ 100 K) | v1 TSFC (two-gas, uncooled) |
|---|---|---|---|---|---|---|---|
| CFM56-7B-class, M 0.78 / 35 kft | 5.3 | 1.70 | 31.9 | 1500 | 186 | 0.610 (0.583–0.646) | 0.655 |
| LEAP-1A-class, M 0.78 / 35 kft | 11 | 1.45 | 40.0 | 1600 | 126 | 0.522 (0.498–0.553) | 0.553 |
| GE90-115B-class, M 0.83 / 35 kft | 8 | 1.55 | 42.1 | 1550 | 140 | 0.563 (0.546–0.592) | 0.600 |
| F100-class, M 0.9 / 35 kft | 0.36 | 3.8 | 31.9 | 1550 | 555 | 0.978 (0.928–1.031) | 1.053 |

* **CFM56 class:** now about 3 % below the ≈0.63–0.67 public cruise-TSFC region; only the hot end of the Tt4 band reaches it.
  The inputs were not retuned. The one input change is that Tt4 was raised by 100 K for every case, because Tt4 is now the
  combustor exit temperature, with cooling modelled explicitly, instead of an equivalent uncooled value.
* **LEAP vs CFM56 class:** 14.4 % lower TSFC (5–23 % across the Tt4 bands), against a ≈15 % manufacturer claim. About 67 % of
  the gain is propulsive.
* **Trade baseline** (M 0.78 / 35 kft, OPR 40, Tt4 1600 K, ε 0.12): the uninstalled optimum is at BPR ≈ 47. The installed
  optimum is at BPR ≈ 14, flat within 0.5 % from 12 to 17. It moves to ≈ 9 or ≈ 21 when the penalty is doubled or halved, and
  rises by ≈ 2.3 per 100 K of Tt4.
* **Old two-gas burner balance:** it needed 5.6–7.6 % more fuel than the enthalpy balance at the case engines' burner conditions.
  Over the whole cycle, switching the gas model lowers CFM56-class TSFC by 9.4 % at unchanged inputs.

## Limitations

* Simple cooling model: one coolant stream (fixed fraction) bypasses the burner and HPT and is remixed losslessly ahead of the
  LPT. There is no NGV/rotor split, customer bleed or power extraction.
* Complete lean combustion without dissociation, a CH1.92 fuel surrogate, LHV and fuel temperature at 298.15 K.
* Off-design uses constant component efficiencies and choked turbine nozzles. There are no maps, surge margin or speed limits.
* The installation model is crude: skin friction plus weight ∝ D², with no spillage, core-cowl scrubbing or ground clearance.
* The low-BPR case is modelled as separate-flow and dry. The real F100 is mixed-flow with an afterburner.
