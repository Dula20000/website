"""Regenerate every result used by the research page into ../data/*.json (+ CSV thrust curves).

Run:  python3 run_study.py          (about a minute on a laptop)
"""
from __future__ import annotations

import copy
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from srm import cases  # noqa: E402
from srm.ballistics import Motor, downsample, motor_class  # noqa: E402
from srm.geometry import (GrainTable, contour_segments, enclosed_area, fmm, make_grid,  # noqa: E402
                          port_sdf, clipped_length, star_offset_perimeter, bates_burning_area)
from srm.montecarlo import QSEnsemble, monte_carlo, oat_sensitivity  # noqa: E402
from srm.propellant import P_REF  # noqa: E402

DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
os.makedirs(DATA, exist_ok=True)

JS_GRID = {"N": 101, "nw": 120}  # grid used by the in-browser tool


def rnd(x, s=6):
    """Round floats (recursively) to s significant digits for compact JSON."""
    if isinstance(x, dict):
        return {k: rnd(v, s) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, s) for v in x]
    if isinstance(x, (np.floating, float)):
        if not math.isfinite(x):
            return None
        return float(f"{x:.{s}g}")
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.ndarray):
        return rnd(x.tolist(), s)
    return x


def dump(name, obj, sig=6):
    with open(os.path.join(DATA, name), "w") as f:
        json.dump(rnd(obj, sig), f, separators=(",", ":"))
    print("  wrote", name)


def order(errs, hs):
    e, h = np.log(np.asarray(errs)), np.log(np.asarray(hs))
    return [float((e[i + 1] - e[i]) / (h[i + 1] - h[i])) for i in range(len(e) - 1)]


def chain_segments(x1, y1, x2, y2):
    """Join oriented marching-squares segments into polylines (lists of points)."""
    key = lambda x, y: (round(x * 1e9), round(y * 1e9))  # noqa: E731
    start = {}
    for k in range(len(x1)):
        start.setdefault(key(x1[k], y1[k]), []).append(k)
    used = np.zeros(len(x1), dtype=bool)
    lines = []
    for k0 in range(len(x1)):
        if used[k0]:
            continue
        line = [(x1[k0], y1[k0])]
        k = k0
        while k is not None and not used[k]:
            used[k] = True
            line.append((x2[k], y2[k]))
            nxt = [q for q in start.get(key(x2[k], y2[k]), []) if not used[q]]
            k = nxt[0] if nxt else None
        lines.append(line)
    return lines


def contour_xy(tab, levels, scale=1e3):
    """Null-separated polylines (mm) of the burning front at the given web distances (clipped to the case)."""
    out = []
    for lv in levels:
        xs, ys = [], []
        for line in chain_segments(*tab.contour(lv)):
            xs += [round(float(p[0]) * scale, 3) for p in line] + [None]
            ys += [round(float(p[1]) * scale, 3) for p in line] + [None]
        out.append({"w": lv * scale, "x": xs, "y": ys})
    return out


def run_case(spec, n_pts=500):
    m = Motor(spec)
    qs = m.quasi_steady()
    tr = m.transient()
    met = m.metrics(tr, qs)
    mass_out = float(np.sum(0.5 * (tr["mdot"][1:] + tr["mdot"][:-1]) * np.diff(tr["t"])))
    met["mass_balance_err"] = mass_out / m.m_prop - 1  # gas left in the chamber at cut-off is ~1e-4
    ds = downsample(tr, n_pts)
    qsd = {k: qs[k][:: max(1, len(qs["w"]) // 300)].tolist() for k in ("w", "t", "Pc", "F", "Kn")}
    return m, qs, tr, met, {"ode": ds, "qs": qsd, "metrics": met}


def write_csv(name, tr):
    with open(os.path.join(DATA, name), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["t_s", "thrust_N", "Pc_Pa", "Kn"])
        for i in range(len(tr["t"])):
            wr.writerow([f"{tr['t'][i]:.5f}", f"{tr['F'][i]:.3f}", f"{tr['Pc'][i]:.1f}", f"{tr['Kn'][i]:.2f}"])
    print("  wrote", name)


# --------------------------------------------------------------------------
# 1. verification
# --------------------------------------------------------------------------

def verification():
    print("verification ...")
    V = {}
    R = 0.045
    # (i) circular port: perimeter 2 pi (r0 + w)
    circ = {"N": [], "h": [], "err_max": [], "err_area0": []}
    for N in (41, 81, 161, 321, 641):
        g = GrainTable({"type": "bates", "core_d": 0.030}, R, N=N, nw=40)
        m = g.w < g.w_contact - g.w[1]
        ex = 2 * math.pi * (0.015 + g.w)
        circ["N"].append(N); circ["h"].append(g.h)
        circ["err_max"].append(float(np.max(np.abs(g.P_raw[m] - ex[m]) / ex[m])))
        circ["err_area0"].append(abs(g.A0 / (math.pi * 0.015 ** 2) - 1))
    circ["order"] = order(circ["err_max"], circ["h"])
    V["circle"] = circ

    # (ii) star: analytic offset perimeter, FMM field vs exact-distance contouring
    geo = {"type": "star", "points": 6, "r_tip": 0.030, "r_valley": 0.012}
    star = {"geom": geo, "N": [], "h": [], "err_fmm": [], "err_exact": [], "err_fmm_max": [],
            "field_linf": [], "field_l1": []}
    for N in (51, 101, 201, 401, 801):
        x, X, Y, h = make_grid(R, N)
        phi = port_sdf(geo, X, Y)
        T = fmm(phi, h)
        _, wv = star_offset_perimeter(6, 0.030, 0.012, [0.0])
        w_hi = min(wv, R - 0.030) - 0.0005
        w = np.linspace(0.001, w_hi, 25)
        Pex, _ = star_offset_perimeter(6, 0.030, 0.012, w)
        Pf = np.array([clipped_length(*contour_segments(T, x, l), R) for l in w])
        Pe = np.array([clipped_length(*contour_segments(phi, x, l), R) for l in w])
        pos = (phi > 0) & (np.hypot(X, Y) <= R)
        star["N"].append(N); star["h"].append(h)
        star["err_fmm"].append(float(np.mean(np.abs(Pf - Pex) / Pex)))
        star["err_fmm_max"].append(float(np.max(np.abs(Pf - Pex) / Pex)))
        star["err_exact"].append(float(np.mean(np.abs(Pe - Pex) / Pex)))
        star["field_linf"].append(float(np.max(np.abs(T - phi)[pos])))
        star["field_l1"].append(float(np.mean(np.abs(T - phi)[pos])))
        if N == 201:
            star["curve"] = {"w": (w * 1e3).tolist(), "P_exact": (Pex * 1e3).tolist(), "P_fmm": (Pf * 1e3).tolist()}
    star["w_valid"] = wv
    star["order_fmm"] = order(star["err_fmm"], star["h"])
    star["order_exact"] = order(star["err_exact"], star["h"])
    star["order_field_l1"] = order(star["field_l1"], star["h"])
    V["star"] = star

    # (iii) BATES burning area of the APCP motor vs closed form
    spec = cases.apcp_bates()
    m = Motor(spec)
    w = np.linspace(0, m.w_web * 0.999, 120)
    Ab = np.array([m.ab_total(x) for x in w])
    D, d, L = 2 * spec["R"], spec["segments"][0]["geom"]["core_d"], spec["segments"][0]["L"]
    Abx = 4 * bates_burning_area(D, d, L, w, ends=2)
    V["bates"] = {"w": (w * 1e3).tolist(), "Ab_grid": Ab.tolist(), "Ab_exact": Abx.tolist(),
                  "err_max": float(np.max(np.abs(Ab - Abx) / Abx)), "N": spec["grid"]["N"],
                  "vprop_err": m.m_prop / (spec["prop"]["rho"] * 4 * 0.25 * math.pi * (D * D - d * d) * L) - 1}
    # (iv) raw area closure (co-area consistency) vs grid for the four cross-sections of study (a)
    geos = cases.equal_mass_geometries(R)
    geos["finocyl"]["slot_r"] = 0.023
    clo = {"N": [81, 161, 321]}
    for k, g in geos.items():
        clo[k] = [GrainTable(g, R, N=N, nw=160).closure for N in clo["N"]]
    V["closure"] = clo
    return V


# --------------------------------------------------------------------------
# 2. case studies
# --------------------------------------------------------------------------

def case_studies(V):
    print("case studies ...")
    out = {}
    specs = {"knsb": cases.knsb_bates(), "apcp": cases.apcp_bates(), "srb": cases.srb_class()}
    conserv = {}
    for key, spec in specs.items():
        t0 = time.time()
        m, qs, tr, met, rec = run_case(spec)
        tab = m.segs[0]["tab"]
        levels = list(np.linspace(0, tab.w_max * 0.98, 7))
        rec["contours"] = contour_xy(tab, levels)
        rec["R_mm"] = spec["R"] * 1e3
        rec["name"] = spec["name"]
        if key == "srb":
            tab2 = m.segs[1]["tab"]
            rec["contours_aft"] = contour_xy(tab2, list(np.linspace(0, tab2.w_max * 0.98, 5)))
        rec["qs_vs_ode"] = {"dI": met["qs_I"] / met["I_total"] - 1, "dPc": met["qs_Pc_max"] / met["Pc_max"] - 1}
        # mid-burn agreement (ODE vs equilibrium at the same web distance)
        wm = 0.5 * m.w_web
        i = int(np.searchsorted(tr["w"], wm))
        rec["qs_vs_ode"]["mid_burn_dPc"] = float(np.interp(wm, qs["w"], qs["Pc"]) / tr["Pc"][i] - 1)
        conserv[key] = {"area_closure": float(max(abs(sg["tab"].closure) for sg in m.segs)),
                        "mass_balance": met["mass_balance_err"]}
        out[key] = rec
        write_csv(f"thrust_{key}.csv", tr)
        print(f"  {key}: {time.time() - t0:.1f}s  I={met['I_total']:.4g}  Pc={met['Pc_max'] / 1e6:.3f} MPa  {met['designation']}")
    V["conservation"] = conserv

    # switches: erosive burning on the KNSB motor, throat erosion on the APCP motor
    s = cases.knsb_bates(); s["erosive"] = {"on": True, "alpha": 5.0e-6, "beta": 53.0}
    m, qs, tr, met, rec = run_case(s, 400)
    out["knsb_erosive"] = {"ode": rec["ode"], "metrics": met, "alpha": 5.0e-6, "beta": 53.0,
                           "port_to_throat": math.pi * 0.0075 ** 2 / m.At()}
    s = cases.apcp_bates(); s["throat_erosion"] = 1.0e-4
    m, qs, tr, met, rec = run_case(s, 400)
    out["apcp_throat"] = {"ode": rec["ode"], "metrics": met, "rate": 1.0e-4}
    return out, specs


# --------------------------------------------------------------------------
# 3. study (a): geometry vs thrust-curve shape at equal propellant mass
# --------------------------------------------------------------------------

def size_throat_for_pc(spec, Pc_target):
    """Throat diameter giving a quasi-steady peak pressure Pc_target (single-law propellant)."""
    m = Motor(spec)
    a, n = m.prop.laws[0][2], m.prop.laws[0][3]
    w_all = max(sg["w_end"] for sg in m.segs)
    Abmax = max(m.ab_total(x) for x in np.linspace(0, w_all, 600)[:-1])
    Kn = Pc_target ** (1 - n) * P_REF ** n / (m.prop.rho * a * m.prop.cstar)
    return math.sqrt(4 * Abmax / Kn / math.pi)


def study_geometry():
    print("study (a): geometry ...")
    R, L, Pt = 0.045, 0.52, 6.0e6
    geoms = cases.equal_mass_geometries(R)
    A_target = math.pi * 0.015 ** 2
    # finocyl slot length to match the port area (area from the contour of phi0 on a fine grid)
    x, X, Y, h = make_grid(R, 801)
    lo, hi = 0.012, 0.040
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        g = dict(geoms["finocyl"], slot_r=mid)
        A = enclosed_area(*contour_segments(port_sdf(g, X, Y), x, 0.0))
        lo, hi = (mid, hi) if A < A_target else (lo, mid)
    geoms["finocyl"]["slot_r"] = round(0.5 * (lo + hi), 6)
    res = {"R_mm": R * 1e3, "L": L, "Pc_target": Pt, "A_port0_target": A_target, "items": {}}
    for key, geo in geoms.items():
        spec = cases.base("apcp", R, [{"geom": geo, "L": L, "ends": 0, "count": 1}], Dt=0.02, eps=8.0,
                          name=f"equal-mass {key}")
        spec["Dt"] = size_throat_for_pc(spec, Pt)
        m, qs, tr, met, rec = run_case(spec, 400)
        tab = m.segs[0]["tab"]
        rec["contours"] = contour_xy(tab, list(np.linspace(0, tab.w_max * 0.98, 8)))
        rec["geom"] = geo
        rec["Dt"] = spec["Dt"]
        rec["A_port0"] = tab.A0
        prog = met["shape_ratio"]
        rec["shape"] = "progressive" if prog > 1.15 else ("regressive" if prog < 0.87 else "near-neutral")
        # perimeter history (normalised) for the geometry plot
        rec["perimeter"] = {"w": (tab.w * 1e3).tolist()[::2], "P": (tab.P * 1e3).tolist()[::2]}
        res["items"][key] = rec
        print(f"  {key}: Dt={spec['Dt'] * 1e3:.2f} mm  prog={prog:.3f}  I={met['I_total']:.0f}  tb={met['t_burn']:.2f}")
    return res, geoms


# --------------------------------------------------------------------------
# 4. study (b): BATES segment count x core diameter
# --------------------------------------------------------------------------

def study_bates():
    print("study (b): BATES trade ...")
    R, Ltot, Pt = 0.045, 0.52, 6.0e6
    D = 2 * R
    nseg = list(range(1, 9))
    dD = [round(x, 3) for x in np.arange(0.20, 0.551, 0.025)]
    neut = np.zeros((len(dD), len(nseg))); imp = np.zeros_like(neut); mass = np.zeros_like(neut)
    avgpc = np.zeros_like(neut); tb = np.zeros_like(neut)
    for i, r in enumerate(dD):
        for j, N in enumerate(nseg):
            spec = cases.base("apcp", R, [{"geom": {"type": "bates", "core_d": r * D}, "L": Ltot / N,
                                           "ends": 2, "count": N}], Dt=0.02, eps=8.0)
            spec["Dt"] = size_throat_for_pc(spec, Pt)
            m = Motor(spec)
            qs = m.quasi_steady(300)
            neut[i, j] = qs["neutrality_pc"]
            imp[i, j] = qs["I"]; mass[i, j] = m.m_prop
            web = qs["w"] <= m.w_web
            avgpc[i, j] = float(np.trapezoid(qs["Pc"][web], qs["t"][web]) / qs["t"][web][-1])
            tb[i, j] = qs["t"][-1]
    k = np.unravel_index(np.argmin(neut), neut.shape)
    # Nakka-style rule: L_opt = (3 D + d) / 2  ->  optimum segment count for each d
    rule = [Ltot / ((3 * D + r * D) / 2) for r in dD]
    best_per_d = [nseg[int(np.argmin(neut[i]))] for i in range(len(dD))]
    return {"nseg": nseg, "dD": dD, "neutrality": neut.tolist(), "impulse": imp.tolist(), "mass": mass.tolist(),
            "avg_pc": avgpc.tolist(), "t_burn": tb.tolist(), "rule_nseg": rule, "best_nseg_per_d": best_per_d,
            "best": {"dD": dD[k[0]], "nseg": nseg[k[1]], "neutrality": float(neut[k]),
                     "L_seg": Ltot / nseg[k[1]], "L_rule": (3 * D + dD[k[0]] * D) / 2},
            "Ltot": Ltot, "D": D, "Pc_target": Pt}


# --------------------------------------------------------------------------
# 5. study (c): Monte Carlo
# --------------------------------------------------------------------------

def study_mc():
    print("study (c): Monte Carlo ...")
    spec = cases.apcp_bates()
    m = Motor(spec)
    sig = {"a": 0.03, "n": 0.02, "Dt": 0.005, "rho": 0.01}
    mc = monte_carlo(m, sig, n_samples=20000)
    nomP = mc["nominal"]["Pc_max"]
    mc_anchor = monte_carlo(m, sig, n_samples=20000, anchor_P=nomP)
    oat = oat_sensitivity(m, sig)

    def hist(v, bins=60):
        c, e = np.histogram(v, bins=bins)
        return {"counts": c.tolist(), "edges": e.tolist()}

    FS = 1.5
    res = {"sig": sig, "n_samples": mc["n_samples"], "nominal": mc["nominal"],
           "Pc_max": mc["Pc_max"], "I": mc["I"], "src_lnPc": mc["src_lnPc"], "src_lnI": mc["src_lnI"],
           "r2_lnPc": mc["r2_lnPc"],
           "hist_Pc": hist(mc["samples"]["Pc_max"] / 1e6), "hist_I": hist(mc["samples"]["I"]),
           "anchor": {"P": nomP, "Pc_max": mc_anchor["Pc_max"], "src_lnPc": mc_anchor["src_lnPc"]},
           "oat": oat, "n0": m.prop.laws[0][3],
           "MEOP": mc["Pc_max"]["p99865"], "MEOP_factor": mc["Pc_max"]["p99865"] / nomP,
           "FS_ultimate": FS, "burst_required": FS * mc["Pc_max"]["p99865"],
           "design_factor": FS * mc["Pc_max"]["p99865"] / nomP}
    # sweep of the nominal exponent at fixed nominal peak pressure
    sweep = {"n": [], "MEOP_factor": [], "cv": [], "share_n": [], "cv_anchor": []}
    ens0 = QSEnsemble(m)
    for n0 in np.round(np.arange(0.20, 0.701, 0.05), 3):
        s2 = copy.deepcopy(spec)
        KnMax = float(ens0.Ab.max() / m.At())
        a_new = nomP ** (1 - n0) * P_REF ** n0 / (m.prop.rho * m.prop.cstar * KnMax)
        s2["prop"]["laws"] = [[0.0, 1e12, a_new, float(n0)]]
        m2 = Motor(s2)
        r = monte_carlo(m2, sig, n_samples=8000, seed=7)
        ra = monte_carlo(m2, sig, n_samples=8000, seed=7, anchor_P=nomP)
        v = np.array([r["src_lnPc"][p] for p in ("a", "n", "Dt", "rho")]) ** 2
        sweep["n"].append(float(n0)); sweep["MEOP_factor"].append(r["Pc_max"]["p99865"] / r["nominal"]["Pc_max"])
        sweep["cv"].append(r["Pc_max"]["std"] / r["Pc_max"]["mean"])
        sweep["cv_anchor"].append(ra["Pc_max"]["std"] / ra["Pc_max"]["mean"])
        sweep["share_n"].append(float(v[1] / v.sum()))
    res["sweep"] = sweep
    print(f"  nominal Pc_max {nomP / 1e6:.3f} MPa, MEOP {res['MEOP'] / 1e6:.3f} MPa (x{res['MEOP_factor']:.3f})")
    return res


# --------------------------------------------------------------------------
# 6. JS cross-check reference
# --------------------------------------------------------------------------

def js_reference():
    ref = {}
    for key, fn in (("knsb", cases.knsb_bates), ("apcp", cases.apcp_bates)):
        spec = fn(N=JS_GRID["N"], nw=JS_GRID["nw"])
        m = Motor(spec)
        qs = m.quasi_steady()
        tr = m.transient()
        met = m.metrics(tr, qs)
        ref[key] = {k: met[k] for k in ("I_total", "Pc_max", "F_max", "t_burn", "Isp", "Kn_initial",
                                         "neutrality_pc", "qs_I", "qs_Pc_max")}
        ref[key]["P0"] = float(m.segs[0]["tab"].P[0]); ref[key]["P_mid"] = float(m.segs[0]["tab"].P[JS_GRID["nw"] // 2])
        ref[key]["w_max"] = float(m.segs[0]["tab"].w_max)
    return ref


def run_node_check():
    node = shutil.which("node")
    if not node:
        print("  node not found: skipping JS cross-check")
        return None
    script = os.path.join(HERE, "tests", "js_crosscheck.js")
    r = subprocess.run([node, script], capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        print(r.stderr)
        return None
    with open(os.path.join(DATA, "js_crosscheck.json")) as f:
        return json.load(f)


def presets_for_js(specs, geo_study):
    """Case inputs exported for the page (tool presets and inputs table)."""
    P = {}
    for k, s in specs.items():
        P[k] = s
    for k, it in geo_study["items"].items():
        s = cases.base("apcp", geo_study["R_mm"] / 1e3, [{"geom": it["geom"], "L": geo_study["L"], "ends": 0,
                                                         "count": 1}], Dt=it["Dt"], eps=8.0,
                       name=f"APCP equal-mass {k}")
        P["geo_" + k] = s
    props = {k: f().to_dict() for k, f in cases.PROPS.items()}
    return {"presets": P, "props": props, "prop_sources": cases.PROP_SOURCES, "js_grid": JS_GRID,
            "knsb_thermo": {"T0": 1600.0, "M": 39.9, "gamma": 1.1361, "rho_ideal": 1841.0, "rho_frac": 0.95}}


def main():
    t0 = time.time()
    V = verification()
    out, specs = case_studies(V)
    geo, geoms = study_geometry()
    bates = study_bates()
    mc = study_mc()
    V["js_ref"] = js_reference()
    dump("verification.json", V, 12)
    dump("cases.json", presets_for_js(specs, geo), 12)
    dump("results.json", {"cases": out, "geometry": geo, "bates": bates})
    dump("montecarlo.json", mc)
    js = run_node_check()
    if js:
        print("  JS cross-check max rel diff:", js["max_rel_diff"])
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
