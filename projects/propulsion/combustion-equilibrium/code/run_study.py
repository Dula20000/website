"""Regenerate every number used by the research page.

Writes (relative to this file):
    ../data/thermo.json        coefficient + propellant tables (input to the JS port)
    ../data/verification.json  V&V results
    ../data/cases.json         engine case studies
    ../data/sweeps.json        O/F and Pc sweeps, exhaust composition
    ../data/sweeps_*.csv       the O/F sweeps as CSV

Run:  python3 run_study.py
"""
from __future__ import annotations

import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "tests"))
DATA = os.path.join(HERE, "..", "data")

from eqsolver import thermo  # noqa: E402
from eqsolver.thermo import CHO_SPECIES, ATOMIC_MASS  # noqa: E402
from eqsolver.equilibrium import equilibrate, default_system  # noqa: E402
from eqsolver.rocket import (PROPELLANTS, PAIRS, chamber, nozzle, performance,  # noqa: E402
                             bulk_density, bipropellant, G0, P_ATM)
import test_thermo  # noqa: E402
import test_equilibrium as teq  # noqa: E402
import test_rocket as troc  # noqa: E402


def rnd(x, n=6):
    if isinstance(x, dict):
        return {k: rnd(v, n) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, n) for v in x]
    if isinstance(x, (float, np.floating)):
        x = float(x)
        if x == 0 or not math.isfinite(x):
            return x
        return float(f"{x:.{n}g}")
    if isinstance(x, np.integer):
        return int(x)
    return x


def dump(name, obj):
    with open(os.path.join(DATA, name), "w") as f:
        json.dump(obj, f, indent=1)
    print("wrote", name)


# --------------------------------------------------------------------------
PC_SWEEP = 10e6          # Pa, common chamber pressure for the O/F sweeps
EPS_VAC = 40.0           # vacuum nozzle for the sweeps
EPS_SL = 12.0            # sea-level nozzle for the sweeps (near-optimal at 10 MPa)

OF_GRID = {
    "LOX/LH2": np.round(np.arange(2.5, 10.001, 0.1), 3),
    "LOX/LCH4": np.round(np.arange(2.0, 5.001, 0.05), 3),
    "LOX/RP-1": np.round(np.arange(1.6, 3.801, 0.05), 3),
}

ENGINES = [
    # representative public values (approximate); delivered Isp as publicly quoted (≈)
    {"id": "merlin1d", "name": "Merlin 1D-class", "pair": "LOX/RP-1", "cycle": "gas generator",
     "Pc": 9.7e6, "eps": 16.0, "of": 2.36, "isp_sl_pub": 282.0, "isp_vac_pub": 311.0,
     "source": "SpaceX public figures / press material (approx.)"},
    {"id": "rs25", "name": "RS-25-class", "pair": "LOX/LH2", "cycle": "fuel-rich staged combustion",
     "Pc": 20.6e6, "eps": 69.0, "of": 6.0, "isp_sl_pub": 366.0, "isp_vac_pub": 452.0,
     "source": "NASA / Aerojet Rocketdyne public data (approx.)"},
    {"id": "raptor", "name": "Raptor-class (sea-level)", "pair": "LOX/LCH4",
     "cycle": "full-flow staged combustion",
     "Pc": 30.0e6, "eps": 40.0, "of": 3.6, "isp_sl_pub": 327.0, "isp_vac_pub": 347.0,
     "source": "SpaceX public statements (approx.; varies by version)"},
    {"id": "rd180", "name": "RD-180-class", "pair": "LOX/RP-1", "cycle": "ox-rich staged combustion",
     "Pc": 26.7e6, "eps": 36.9, "of": 2.72, "isp_sl_pub": 311.0, "isp_vac_pub": 338.0,
     "source": "NPO Energomash / ULA public data (approx.)"},
]

PC_LIST = [0.5e6, 1e6, 2e6, 3e6, 5e6, 7e6, 10e6, 15e6, 20e6, 25e6, 30e6]
PC_OF = {"LOX/LH2": 6.0, "LOX/LCH4": 3.6, "LOX/RP-1": 2.6}


def thermo_json():
    obj = {"species": thermo.to_json_dict(), "cho_species": CHO_SPECIES,
           "atomic_mass": ATOMIC_MASS, "R": thermo.R, "P_ref": thermo.P_REF,
           "propellants": {k: {"comp": v.comp, "h": v.h, "T": v.T,
                               "density": None if math.isnan(v.density) else v.density,
                               "note": v.note}
                           for k, v in PROPELLANTS.items()},
           "pairs": {k: {"fuel": f, "ox": o} for k, (f, o) in PAIRS.items()}}
    dump("thermo.json", obj)


def composition(st, thr=1e-7):
    return {k: v for k, v in st.mole_fractions().items() if v > thr}


# --------------------------------------------------------------------------
def sweeps():
    out = {"conditions": {"Pc_MPa": PC_SWEEP / 1e6, "eps_vac": EPS_VAC, "eps_sl": EPS_SL,
                          "Pa_kPa": P_ATM / 1e3}, "of": {}, "pc": {}, "optimum": {}}
    csv_rows = []
    for pair, (fuel, ox) in PAIRS.items():
        grid = OF_GRID[pair]
        rec = {k: [] for k in ("of", "Tc", "M_c", "gamma_s_c", "gamma_fr_c", "cstar_shift", "cstar_froz",
                                "isp_vac_shift", "isp_vac_froz", "isp_sl_shift", "isp_sl_froz",
                                "rho_bulk", "rho_isp_vac", "rho_isp_sl")}
        comp = {}
        for of in grid:
            ch = chamber(fuel, ox, float(of), PC_SWEEP)
            sv = nozzle(ch, EPS_VAC)
            fv = nozzle(ch, EPS_VAC, frozen=True)
            ss = nozzle(ch, EPS_SL)
            fs = nozzle(ch, EPS_SL, frozen=True)
            rho = bulk_density(fuel, ox, float(of))
            vals = dict(of=float(of), Tc=ch.T, M_c=ch.molar_mass * 1e3, gamma_s_c=ch.gamma_s,
                        gamma_fr_c=ch.gamma_fr, cstar_shift=sv.cstar, cstar_froz=fv.cstar,
                        isp_vac_shift=sv.Isp_vac, isp_vac_froz=fv.Isp_vac,
                        isp_sl_shift=ss.Isp_sl, isp_sl_froz=fs.Isp_sl, rho_bulk=rho,
                        rho_isp_vac=rho * sv.Isp_vac / 1e3, rho_isp_sl=rho * ss.Isp_sl / 1e3)
            for k, v in vals.items():
                rec[k].append(v)
            for sp, x in ch.mole_fractions().items():
                comp.setdefault(sp, []).append(x)
            csv_rows.append({"pair": pair, **vals})
        # keep species that ever exceed 1e-6
        comp = {k: v for k, v in comp.items() if max(v) > 1e-6}
        rec["composition"] = comp
        out["of"][pair] = rnd(rec)

        # optimum O/F by golden-section on shifting vac Isp and on Tc / density-Isp
        def golden(fn, a, b, tol=1e-4):
            gr = (math.sqrt(5) - 1) / 2
            c, d = b - gr * (b - a), a + gr * (b - a)
            fc, fd = fn(c), fn(d)
            while b - a > tol:
                if fc > fd:
                    b, d, fd = d, c, fc
                    c = b - gr * (b - a)
                    fc = fn(c)
                else:
                    a, c, fc = c, d, fd
                    d = a + gr * (b - a)
                    fd = fn(d)
            x = 0.5 * (a + b)
            return x, fn(x)

        g = np.array(rec["of"])
        def bracket(key):
            i = int(np.argmax(rec[key]))
            return float(g[max(i - 2, 0)]), float(g[min(i + 2, len(g) - 1)])
        opt = {}
        a, b = bracket("isp_vac_shift")
        opt["isp_vac_shift"] = golden(lambda of: performance(fuel, ox, of, PC_SWEEP, EPS_VAC).Isp_vac, a, b)
        a, b = bracket("isp_vac_froz")
        opt["isp_vac_froz"] = golden(lambda of: performance(fuel, ox, of, PC_SWEEP, EPS_VAC, frozen=True).Isp_vac, a, b)
        a, b = bracket("isp_sl_shift")
        opt["isp_sl_shift"] = golden(lambda of: performance(fuel, ox, of, PC_SWEEP, EPS_SL).Isp_sl, a, b)
        a, b = bracket("Tc")
        opt["Tc"] = golden(lambda of: chamber(fuel, ox, of, PC_SWEEP).T, a, b)
        a, b = bracket("rho_isp_vac")
        if b >= g[-1] - 1e-9:
            opt["rho_isp_vac"] = (float(g[-1]), float(rec["rho_isp_vac"][-1]))
            opt["rho_isp_vac_at_edge"] = True
        else:
            opt["rho_isp_vac"] = golden(lambda of: bulk_density(fuel, ox, of) / 1e3 *
                                        performance(fuel, ox, of, PC_SWEEP, EPS_VAC).Isp_vac, a, b)
        el, b0, _ = bipropellant(fuel, ox, float(grid[0]))
        opt["oc_min"] = float(b0[el.index("O")] / b0[el.index("C")]) if "C" in el else None
        out["optimum"][pair] = rnd({k: ({"of": v[0], "value": v[1]} if isinstance(v, tuple) else v)
                                    for k, v in opt.items()})

    # Pc sweep at fixed O/F, fixed ε (vacuum), shifting & frozen
    for pair, (fuel, ox) in PAIRS.items():
        of = PC_OF[pair]
        rec = {k: [] for k in ("Pc_MPa", "Tc", "M_c", "isp_vac_shift", "isp_vac_froz",
                                "x_dissoc", "cstar_shift")}
        for Pc in PC_LIST:
            ch = chamber(fuel, ox, of, Pc)
            sv = nozzle(ch, EPS_VAC)
            fv = nozzle(ch, EPS_VAC, frozen=True)
            x = ch.mole_fractions()
            rec["Pc_MPa"].append(Pc / 1e6)
            rec["Tc"].append(ch.T)
            rec["M_c"].append(ch.molar_mass * 1e3)
            rec["isp_vac_shift"].append(sv.Isp_vac)
            rec["isp_vac_froz"].append(fv.Isp_vac)
            rec["cstar_shift"].append(sv.cstar)
            rec["x_dissoc"].append(sum(x.get(k, 0) for k in ("H", "O", "OH", "HO2")))
        out["pc"][pair] = rnd({"of": of, **rec})

    dump("sweeps.json", out)
    with open(os.path.join(DATA, "sweeps_of.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
        w.writeheader()
        for r in csv_rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote sweeps_of.csv")
    return out


# --------------------------------------------------------------------------
def cases():
    out = []
    for e in ENGINES:
        fuel, ox = PAIRS[e["pair"]]
        ch = chamber(fuel, ox, e["of"], e["Pc"])
        s = nozzle(ch, e["eps"], of=e["of"])
        fz = nozzle(ch, e["eps"], of=e["of"], frozen=True)
        rec = {k: v for k, v in e.items()}
        rec["Pc_MPa"] = e["Pc"] / 1e6
        rec["shifting"] = s.summary()
        rec["frozen"] = fz.summary()
        rec["rho_bulk"] = bulk_density(fuel, ox, e["of"])
        rec["eta_vac_shift"] = e["isp_vac_pub"] / s.Isp_vac
        rec["eta_sl_shift"] = e["isp_sl_pub"] / s.Isp_sl
        rec["eta_vac_froz"] = e["isp_vac_pub"] / fz.Isp_vac
        rec["comp_chamber"] = composition(ch, 1e-5)
        rec["comp_throat"] = composition(s.throat, 1e-5)
        rec["comp_exit_shift"] = composition(s.exit, 1e-5)
        out.append(rnd(rec))
    dump("cases.json", {"engines": out})
    return out


# --------------------------------------------------------------------------
def js_crosscheck(case_list):
    node = shutil.which("node")
    if not node:
        return {"available": False}
    specs = []
    for e in case_list:
        fuel, ox = PAIRS[e["pair"]]
        for mode in ("shifting", "frozen"):
            specs.append({"id": f'{e["id"]}-{mode}', "fuel": fuel, "ox": ox, "of": e["of"],
                          "Pc": e["Pc_MPa"] * 1e6, "eps": e["eps"], "frozen": mode == "frozen"})
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(specs, f)
        tmp = f.name
    try:
        res = subprocess.run([node, os.path.join(HERE, "..", "js", "crosscheck.js"), tmp],
                             capture_output=True, text=True, check=True, timeout=120)
    finally:
        os.unlink(tmp)
    js = json.loads(res.stdout)
    # compare against unrounded python
    rows, worst = [], 0.0
    for spec, r in zip(specs, js["results"]):
        p = performance(spec["fuel"], spec["ox"], spec["of"], spec["Pc"], spec["eps"],
                        frozen=spec["frozen"]).summary()
        diffs = {k: abs(r[k] / p[k] - 1) for k in ("Tc", "M_c", "cstar", "CF_vac", "Isp_vac", "Isp_sl")}
        worst = max(worst, max(diffs.values()))
        rows.append({"id": spec["id"], "py_isp_vac": p["Isp_vac"], "js_isp_vac": r["Isp_vac"],
                     "py_Tc": p["Tc"], "js_Tc": r["Tc"], "max_rel_diff": max(diffs.values())})
    return {"available": True, "max_rel_diff": worst, "rows": rows, "node_ms": js["ms"]}


def verification(case_list):
    v = {}
    v["thermo"] = test_thermo.thermo_check_table()
    # element conservation, mass action, Gibbs minimum over a set of chamber states
    cons = []
    for fuel, of, Pc in [("LH2", 6.0, 20.6e6), ("LCH4", 3.6, 30e6), ("RP-1", 2.36, 9.7e6),
                         ("LH2", 3.0, 1e6), ("RP-1", 3.5, 1e6)]:
        ch = chamber(fuel, "LOX", of, Pc)
        cons.append({"case": f"LOX/{fuel} O/F={of} Pc={Pc/1e6:g} MPa",
                     "elem_rel": float(np.max(np.abs(ch.element_residual(ch._b0)) / ch._b0)),
                     "mass_action": float(max(teq.mass_action_residuals(ch))),
                     "energy_rel": abs(ch.h - ch._h0) / abs(ch._h0),
                     "iterations": ch.iterations})
    v["conservation"] = cons
    ch = chamber("LCH4", "LOX", 3.6, 10e6)
    st = equilibrate(ch.system, ch._b0, ch.P, T=ch.T)
    mn, order = teq.perturbation_check(st)
    v["gibbs_min"] = {"min_dG_over_RT": mn, "order": order, "samples": 200,
                      "case": "LOX/LCH4 O/F 3.6, 10 MPa, T = T_ad"}
    v["scipy"] = teq.scipy_crosscheck()
    v["flame"] = teq.flame_temperatures()
    el, b0, h0 = bipropellant("LH2", "LOX", 6.0)
    st = equilibrate(default_system(el), b0, 20e6, h0=h0)
    v["newton_history"] = st.history
    v["argon"] = troc.argon_exact()
    v["cstar"] = troc.cstar_consistency()
    # published ballpark checks
    p = performance("LH2", "LOX", 6.0, 20.6e6, 69)
    rp = [(of, chamber("RP-1", "LOX", of, 6.895e6).T) for of in np.arange(2.2, 3.41, 0.02)]
    of_pk, T_pk = max(rp, key=lambda t: t[1])
    v["published"] = {"lh2_of6_eps69_isp_vac": p.Isp_vac, "lh2_of6_eps69_Tc": p.chamber.T,
                      "rp1_Tc_peak_1000psia": T_pk, "rp1_Tc_peak_of": float(of_pk)}
    # frozen vs shifting
    fs = []
    for fuel, of in [("LH2", 6.0), ("LCH4", 3.6), ("RP-1", 2.36)]:
        s = performance(fuel, "LOX", of, 10e6, 40)
        f = performance(fuel, "LOX", of, 10e6, 40, frozen=True)
        fs.append({"pair": f"LOX/{fuel}", "of": of, "shift": s.Isp_vac, "froz": f.Isp_vac,
                   "diff_pct": 100 * (s.Isp_vac / f.Isp_vac - 1)})
    v["frozen_vs_shifting"] = fs
    # timing
    t = time.perf_counter()
    for _ in range(5):
        performance("LCH4", "LOX", 3.6, 30e6, 40)
    v["py_ms_per_case"] = (time.perf_counter() - t) / 5 * 1e3
    v["js"] = js_crosscheck(case_list)
    dump("verification.json", rnd(v, 8))
    return v


def main():
    os.makedirs(DATA, exist_ok=True)
    thermo_json()
    c = cases()
    s = sweeps()
    v = verification(c)
    # console summary
    print("\nCase studies (shifting):")
    for e in c:
        sh, fz = e["shifting"], e["frozen"]
        print(f'  {e["name"]:26s} Tc={sh["Tc"]:.0f} K  c*={sh["cstar"]:.0f}  Isp vac {sh["Isp_vac"]:.1f} '
              f'(frozen {fz["Isp_vac"]:.1f})  SL {sh["Isp_sl"]:.1f}  eta_vac={e["eta_vac_shift"]:.3f}')
    print("Optimum O/F:", json.dumps(s["optimum"], indent=None))
    print("Flame T:", v["flame"])
    print("JS cross-check max rel diff:", v["js"].get("max_rel_diff"))


if __name__ == "__main__":
    main()
