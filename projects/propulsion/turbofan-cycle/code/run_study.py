"""Regenerate every result used by the research page into ../data/*.json (+ CSV).

    python3 run_study.py

Writes: cases.json, verification.json, trades.json, offdesign.json, carpet.csv,
cases_summary.csv.  If ``node`` is on the PATH the JavaScript port
(js/tfcycle.js) is run against the Python results and the comparison is stored in
verification.json["js_check"].
"""
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.join(HERE, "..", "data")
os.makedirs(DATA, exist_ok=True)

from tfcycle import atmosphere, DesignInputs, design, ideal_turbofan, optimum_fan_pr, Engine  # noqa: E402
from tfcycle.ondesign import energy_residual, ts_path, TSFC_LB, marginal_optimum_ratio, dh_air_poly, burner_far
from tfcycle.thermo import Mixture, CPG  # noqa: E402
from tfcycle.installation import fan_diameter, penalty_constant, installed  # noqa: E402
sys.path.insert(0, os.path.join(HERE, "tests"))
import test_cycle as T  # noqa: E402

FT = 0.3048
G = 9.80665


def clean(x, sig=None):
    """Make an object JSON-safe (NaN/inf -> None), optionally rounding floats."""
    if isinstance(x, dict):
        return {k: clean(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v, sig) for v in x]
    if isinstance(x, float):
        if not math.isfinite(x):
            return None
        if sig:
            return float(f"{x:.{sig}g}")
        return x
    return x


def dump(name, obj, sig=None):
    with open(os.path.join(DATA, name), "w") as fh:
        json.dump(clean(obj, sig), fh, separators=(",", ":"))
    print("wrote", name)


def summary(r):
    keys = ["valid", "Fs", "S", "S_mg", "S_lb", "f", "f_b", "eta_th", "eta_p", "eta_o", "V9", "V19",
            "Vfe9", "Vfe19", "M9", "M19", "core_choked", "fan_choked", "T0", "P0", "V0", "OPR",
            "tau_tH", "tau_tL", "pi_tH", "pi_tL", "T9", "P9", "T19", "P19"]
    out = {k: r.get(k) for k in keys}
    if r.get("valid"):
        out["eta_ad"] = r["eta_ad"]
        out["stations"] = r["stations"]
    return out


# ---------------------------------------------------------------------------
# Case-study engines: representative values from public sources + assumptions
# ---------------------------------------------------------------------------
CRUISE = dict(M0=0.78, alt=35000 * FT)
CASES = [
    dict(key="cfm56", name="CFM56-7B-class",
         inputs=DesignInputs(M0=0.78, alt=35000 * FT, alpha=5.3, pi_f=1.70, pi_cL=2.9, pi_cH=11.0,
                             Tt4=1500.0, e_f=0.90, e_cL=0.90, e_cH=0.905, e_tH=0.89, e_tL=0.91,
                             pi_b=0.95, pi_n=0.99, pi_fn=0.99),
         aircraft=dict(label="737-800-class", mass=65000.0, LoD=17.0, n_eng=2),
         public=dict(BPR="≈5.1–5.5", OPR="≈32–33 (take-off)", fanPR="≈1.6–1.7",
                     fan_D_m=61 * 0.0254, fan_D_label="61 in",
                     tsfc_range=[0.63, 0.67], tsfc_note="≈ cruise TSFC region quoted publicly")),
    dict(key="leap1a", name="LEAP-1A-class",
         inputs=DesignInputs(M0=0.78, alt=35000 * FT, alpha=11.0, pi_f=1.45, pi_cL=1.82, pi_cH=22.0,
                             Tt4=1600.0, e_f=0.92, e_cL=0.91, e_cH=0.92, e_tH=0.90, e_tL=0.92,
                             pi_b=0.96, pi_n=0.995, pi_fn=0.995),
         aircraft=dict(label="A320neo-class", mass=64000.0, LoD=17.5, n_eng=2),
         public=dict(BPR="≈11", OPR="≈40 (take-off)", fanPR="≈1.4–1.5 (estimate)",
                     fan_D_m=78 * 0.0254, fan_D_label="78 in",
                     improvement=0.15, improvement_note="≈15 % lower fuel burn than CFM56 (manufacturer claim)")),
    dict(key="ge90", name="GE90-115B-class",
         inputs=DesignInputs(M0=0.83, alt=35000 * FT, alpha=8.0, pi_f=1.55, pi_cL=1.83, pi_cH=23.0,
                             Tt4=1550.0, e_f=0.91, e_cL=0.905, e_cH=0.91, e_tH=0.895, e_tL=0.915,
                             pi_b=0.955, pi_n=0.99, pi_fn=0.99),
         aircraft=dict(label="777-300ER-class", mass=290000.0, LoD=19.0, n_eng=2),
         public=dict(BPR="≈8–9 (≈7 quoted for some ratings)", OPR="≈42", fanPR="≈1.5–1.6 (estimate)",
                     fan_D_m=128 * 0.0254, fan_D_label="128 in")),
    dict(key="f100", name="F100-class (low BPR)",
         inputs=DesignInputs(M0=0.90, alt=35000 * FT, alpha=0.36, pi_f=3.8, pi_cL=3.8, pi_cH=8.4,
                             Tt4=1550.0, e_f=0.88, e_cL=0.88, e_cH=0.90, e_tH=0.89, e_tL=0.90,
                             pi_b=0.95, pi_n=0.99, pi_fn=0.99),
         aircraft=None,
         public=dict(BPR="≈0.36–0.7", OPR="≈25–32", fanPR="≈3–3.8 (3-stage fan)",
                     note="real engine is mixed-flow with afterburner; modelled here dry, separate streams")),
]
TT4_BAND = 100.0
EPS_COOL = DesignInputs().eps_cool   # assumed turbine cooling-air fraction (all cases, default 0.12)
TT4_PREV = {"cfm56": 1400.0, "leap1a": 1500.0, "ge90": 1450.0, "f100": 1450.0}   # v1 uncooled-equivalent values


def ideal_counterpart(inp, T0):
    return ideal_turbofan(inp.M0, T0, 1.4, 1004.0, inp.hPR, inp.Tt4, inp.OPR, inp.pi_f, inp.alpha)


def loss_audit(inp):
    """TSFC as non-idealities are switched on one group at a time."""
    steps = []
    ideal = T.ideal_inputs(M0=inp.M0, alt=inp.alt, alpha=inp.alpha, pi_f=inp.pi_f,
                           pi_cL=inp.pi_cL, pi_cH=inp.pi_cH, Tt4=inp.Tt4, hPR=inp.hPR)
    cur = ideal
    seq = [
        ("Ideal cycle", {}),
        ("+ variable cp(T, f), air & products", dict(thermo="varcp")),
        ("+ fuel mass, burner efficiency", dict(neglect_fuel_mass=False, eta_b=inp.eta_b)),
        ("+ inlet, burner, nozzle Δpt", dict(pi_d_max=inp.pi_d_max, pi_b=inp.pi_b, pi_n=inp.pi_n, pi_fn=inp.pi_fn,
                                             mil_spec_recovery=True)),
        ("+ fan & compressor losses", dict(e_f=inp.e_f, e_cL=inp.e_cL, e_cH=inp.e_cH)),
        ("+ turbine & shaft losses", dict(e_tH=inp.e_tH, e_tL=inp.e_tL, eta_mH=inp.eta_mH, eta_mL=inp.eta_mL)),
        ("+ turbine cooling air", dict(eps_cool=inp.eps_cool)),
        ("+ convergent nozzles", dict(core_nozzle=inp.core_nozzle, fan_nozzle=inp.fan_nozzle)),
    ]
    for label, kw in seq:
        cur = cur.with_(**kw)
        r = design(cur)
        steps.append(dict(label=label, S_lb=r["S_lb"], Fs=r["Fs"], eta_th=r["eta_th"], eta_p=r["eta_p"]))
    return steps


def run_cases():
    out = []
    for c in CASES:
        inp = c["inputs"]
        r = design(inp)
        assert r["valid"], c["key"]
        ic = ideal_counterpart(inp, r["T0"])
        tj = design(inp.with_(alpha=0.0))
        band = [design(inp.with_(Tt4=inp.Tt4 + d)) for d in (-TT4_BAND, TT4_BAND)]
        pf_opt, ropt = optimum_fan_pr(inp)
        T_old = TT4_PREV[c["key"]]
        steps = []
        for label, kw in [("v1: two-gas, uncooled, Tt4 %.0f K" % T_old, dict(thermo="cpg", eps_cool=0.0, Tt4=T_old)),
                          ("variable cp, uncooled, Tt4 %.0f K" % T_old, dict(eps_cool=0.0, Tt4=T_old)),
                          ("variable cp, %.0f %% cooling, Tt4 %.0f K" % (100 * inp.eps_cool, T_old), dict(Tt4=T_old)),
                          ("v2: variable cp, %.0f %% cooling, Tt4 %.0f K" % (100 * inp.eps_cool, inp.Tt4), {})]:
            rr = design(inp.with_(**kw))
            steps.append(dict(label=label, S_lb=rr["S_lb"], Fs=rr["Fs"], V19_V9=rr["V19"] / rr["V9"]))
        item = dict(model_steps=steps, key=c["key"], name=c["name"], inputs=inp.to_dict(), public=c["public"],
                    result=summary(r),
                    ideal=dict(Fs=ic["Fs"], S_lb=ic["S"] * TSFC_LB, eta_th=ic["eta_th"], eta_p=ic["eta_p"],
                               eta_o=ic["eta_o"], pi_f_opt=ic["pi_f_opt"]),
                    turbojet=dict(Fs=tj["Fs"], S_lb=tj["S_lb"], eta_th=tj["eta_th"], eta_p=tj["eta_p"],
                                  eta_o=tj["eta_o"]),
                    band=dict(dT=TT4_BAND, S_lb=[band[0]["S_lb"], band[1]["S_lb"]],
                              Fs=[band[0]["Fs"], band[1]["Fs"]]),
                    opt=dict(pi_f=pf_opt, S_lb=ropt["S_lb"], Fs=ropt["Fs"], V19_V9=ropt["V19"] / ropt["V9"]),
                    V19_V9=r["V19"] / r["V9"],
                    ts=ts_path(r), audit=loss_audit(inp))
        if c["aircraft"]:
            a = c["aircraft"]
            Freq = a["mass"] * G / a["LoD"] / a["n_eng"]
            mdot = Freq / r["Fs"]
            st2 = r["stations"]["2"]
            D = fan_diameter(mdot, st2["Tt"], st2["Pt"])
            item["sizing"] = dict(aircraft=a, F_cruise=Freq, mdot0=mdot, D_fan=D,
                                  D_public=c["public"]["fan_D_m"], ratio=D / c["public"]["fan_D_m"],
                                  M2=0.60, hub_tip=0.30)
        out.append(item)
    by = {o["key"]: o for o in out}
    cfm, leap = by["cfm56"], by["leap1a"]
    imp = 1 - leap["result"]["S_lb"] / cfm["result"]["S_lb"]
    # widest spread from the two Tt4 bands (independent +/-100 K on each engine)
    combos = [1 - l / c_ for l in [leap["result"]["S_lb"]] + leap["band"]["S_lb"]
              for c_ in [cfm["result"]["S_lb"]] + cfm["band"]["S_lb"]]
    # gain split: ln(S ratio) = -ln(eta_th ratio) - ln(eta_p ratio)
    lo = math.log(leap["result"]["eta_o"] / cfm["result"]["eta_o"])
    split = dict(thermal=math.log(leap["result"]["eta_th"] / cfm["result"]["eta_th"]) / lo,
                 propulsive=math.log(leap["result"]["eta_p"] / cfm["result"]["eta_p"]) / lo)
    comp = dict(leap_vs_cfm=imp, leap_vs_cfm_range=[min(combos), max(combos)], split=split,
                cfm_in_public_band=cfm["public"]["tsfc_range"][0] <= cfm["result"]["S_lb"] <= cfm["public"]["tsfc_range"][1])
    return out, comp


# ---------------------------------------------------------------------------
# Trade studies (common baseline: LEAP-level technology, M 0.78 / 35 kft)
# ---------------------------------------------------------------------------
BASE = DesignInputs(M0=0.78, alt=35000 * FT, alpha=10.0, pi_f=1.5, pi_cL=1.8, pi_cH=40 / 1.8, Tt4=1600.0,
                    e_f=0.92, e_cL=0.91, e_cH=0.92, e_tH=0.90, e_tL=0.92, pi_b=0.96,
                    pi_n=0.995, pi_fn=0.995)


def golden(fun, a, b, tol=1e-6):
    g = (math.sqrt(5) - 1) / 2
    x1, x2 = b - g * (b - a), a + g * (b - a)
    f1, f2 = fun(x1), fun(x2)
    while b - a > tol:
        if f1 < f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - g * (b - a)
            f1 = fun(x1)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + g * (b - a)
            f2 = fun(x2)
    return 0.5 * (a + b)


def pf_max(inp):
    a, b = 1.01, 8.0
    if design(inp.with_(pi_f=b))["valid"]:
        return b
    for _ in range(50):
        m = 0.5 * (a + b)
        if design(inp.with_(pi_f=m))["valid"]:
            a = m
        else:
            b = m
    return a


def run_trades():
    t0 = time.time()
    out = dict(base=BASE.to_dict())
    # carpet: TSFC vs specific thrust over (BPR, pi_f)
    bprs = [0, 1, 2, 4, 6, 8, 10, 12, 15]
    lines = []
    rows = []
    for a in bprs:
        inp = BASE.with_(alpha=float(a))
        if a == 0:
            r = design(inp.with_(pi_f=1.0))
            lines.append(dict(alpha=a, pi_f=[1.0], Fs=[r["Fs"]], S_lb=[r["S_lb"]]))
            rows.append([a, 1.0, r["Fs"], r["S_lb"]])
            continue
        hi = pf_max(inp)
        pfs = [1.15 + k * 0.025 for k in range(400) if 1.15 + k * 0.025 < hi]
        L = dict(alpha=a, pi_f=[], Fs=[], S_lb=[])
        for pf in pfs:
            r = design(inp.with_(pi_f=pf))
            if r["valid"]:
                L["pi_f"].append(pf); L["Fs"].append(r["Fs"]); L["S_lb"].append(r["S_lb"])
                rows.append([a, pf, r["Fs"], r["S_lb"]])
        pfo, ro = optimum_fan_pr(inp, tol=1e-7)
        L["opt"] = dict(pi_f=pfo, Fs=ro["Fs"], S_lb=ro["S_lb"])
        lines.append(L)
    cross = []
    for pf in [1.3, 1.5, 1.7, 2.0, 2.5, 3.0]:
        C = dict(pi_f=pf, alpha=[], Fs=[], S_lb=[])
        for a in [x * 0.25 for x in range(0, 61)]:
            r = design(BASE.with_(alpha=a, pi_f=pf))
            if r["valid"]:
                C["alpha"].append(a); C["Fs"].append(r["Fs"]); C["S_lb"].append(r["S_lb"])
        cross.append(C)
    out["carpet"] = dict(lines=lines, cross=cross)
    with open(os.path.join(DATA, "carpet.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["BPR", "pi_f", "F_over_mdot0_N_s_per_kg", "TSFC_lb_per_lbf_h"])
        w.writerows([[x if isinstance(x, int) else f"{x:.6g}" for x in row] for row in rows])

    # efficiency split vs BPR at the TSFC-optimal fan pressure ratio
    sp = dict(alpha=[], pi_f=[], Fs=[], S_lb=[], eta_th=[], eta_p=[], eta_o=[], V19_V9=[])
    for k in range(0, 61):
        a = 0.5 * k
        inp = BASE.with_(alpha=a)
        if a == 0:
            r, pfo = design(inp.with_(pi_f=1.0)), 1.0
        else:
            pfo, r = optimum_fan_pr(inp, tol=1e-7)
        for key, val in [("alpha", a), ("pi_f", pfo), ("Fs", r["Fs"]), ("S_lb", r["S_lb"]),
                         ("eta_th", r["eta_th"]), ("eta_p", r["eta_p"]), ("eta_o", r["eta_o"]),
                         ("V19_V9", r["V19"] / r["V9"])]:
            sp[key].append(val)
    a_u = golden(lambda a: optimum_fan_pr(BASE.with_(alpha=a), tol=1e-8)[1]["S"], 15.0, 100.0, 1e-3)
    pf_u, r_u = optimum_fan_pr(BASE.with_(alpha=a_u), tol=1e-8)
    sp["uninstalled_opt"] = dict(alpha=a_u, pi_f=pf_u, S_lb=r_u["S_lb"], Fs=r_u["Fs"])
    out["split"] = sp

    # Tt4 and OPR: TSFC and specific thrust at BPR = 11 with optimal pi_f
    tt = dict(alpha=11.0, Tt4=[1400, 1500, 1600, 1700, 1800], OPR=[15 + 2.5 * k for k in range(19)], S_lb=[], Fs=[], pi_f=[])
    for T4 in tt["Tt4"]:
        rs, rf, rp = [], [], []
        for opr in tt["OPR"]:
            inp = BASE.with_(alpha=11.0, Tt4=float(T4), pi_cH=opr / BASE.pi_cL)
            if not design(inp.with_(pi_f=1.05))["valid"]:
                rs.append(None); rf.append(None); rp.append(None); continue
            pfo, r = optimum_fan_pr(inp, tol=1e-7)
            rs.append(r["S_lb"]); rf.append(r["Fs"]); rp.append(pfo)
        tt["S_lb"].append(rs); tt["Fs"].append(rf); tt["pi_f"].append(rp)
    out["tt4_opr"] = tt

    # Installation: fan diameter and installed TSFC vs BPR
    ref = design(BASE)
    mults = [0.5, 1.0, 2.0]
    pens = {m: penalty_constant(ref, k_mult=m) for m in mults}
    Freq = 64000.0 * G / 17.5 / 2          # A320neo-class cruise thrust per engine
    inst = dict(F_req=Freq, penalty={str(m): pens[m] for m in mults}, curves=[],
                assumptions=dict(M2=0.60, hub_tip=0.30, kD=1.20, kL=1.50, FQ=1.40, m_ref=2400.0,
                                 D_ref=1.55, LoD=17.0))
    alphas = [2 + 0.5 * k for k in range(0, 57)]
    unin = dict(alpha=[], S_lb=[], pi_f=[], D=[], Fs=[])
    for a in alphas:
        inp = BASE.with_(alpha=a)
        pfo, r = optimum_fan_pr(inp, tol=1e-7)
        st2 = r["stations"]["2"]
        pen = pens[1.0]
        phi = pen["c_pen"] / r["Fs"]
        mdot = Freq / (1 - phi) / r["Fs"]
        unin["alpha"].append(a); unin["S_lb"].append(r["S_lb"]); unin["pi_f"].append(pfo)
        unin["Fs"].append(r["Fs"]); unin["D"].append(fan_diameter(mdot, st2["Tt"], st2["Pt"]))
    inst["uninstalled"] = unin
    for m in mults:
        pen = pens[m]
        cur = dict(mult=m, alpha=[], pi_f=[], S_inst=[], S_unin=[], phi=[], D=[], phi_nac=[], phi_wt=[])
        for a in alphas:
            inp = BASE.with_(alpha=a)
            hi = pf_max(inp)

            def Sinst(pf):
                r = design(inp.with_(pi_f=pf))
                if not r["valid"]:
                    return 1e9
                v = installed(r, pen)["S_inst"]
                return v if math.isfinite(v) else 1e9
            pfo = golden(Sinst, 1.05, hi, 1e-7)
            r = design(inp.with_(pi_f=pfo))
            ins = installed(r, pen)
            st2 = r["stations"]["2"]
            mdot = Freq / (1 - ins["phi"]) / r["Fs"]
            cur["alpha"].append(a); cur["pi_f"].append(pfo); cur["S_inst"].append(ins["S_inst"] * TSFC_LB)
            cur["S_unin"].append(r["S_lb"]); cur["phi"].append(ins["phi"])
            cur["phi_nac"].append(pen["c_nac"] / r["Fs"]); cur["phi_wt"].append(pen["c_wt"] / r["Fs"])
            cur["D"].append(fan_diameter(mdot, st2["Tt"], st2["Pt"]))
        # refine the optimum BPR with a golden search on the continuous alpha
        def best_at(a):
            inp = BASE.with_(alpha=a)
            hi = pf_max(inp)

            def Sinst(pf):
                r = design(inp.with_(pi_f=pf))
                if not r["valid"]:
                    return 1e9
                v = installed(r, pen)["S_inst"]
                return v if math.isfinite(v) else 1e9
            pfo = golden(Sinst, 1.05, hi, 1e-7)
            return Sinst(pfo), pfo
        j = min(range(len(alphas)), key=lambda k: cur["S_inst"][k])
        lo_a, hi_a = alphas[max(j - 1, 0)], alphas[min(j + 1, len(alphas) - 1)]
        a_opt = golden(lambda a: best_at(a)[0], lo_a, hi_a, 1e-4)
        S_opt, pf_opt = best_at(a_opt)
        r = design(BASE.with_(alpha=a_opt, pi_f=pf_opt))
        ins = installed(r, pen)
        st2 = r["stations"]["2"]
        cur["opt"] = dict(alpha=a_opt, pi_f=pf_opt, S_inst=S_opt * TSFC_LB, S_unin=r["S_lb"], phi=ins["phi"],
                          Fs=r["Fs"], D=fan_diameter(Freq / (1 - ins["phi"]) / r["Fs"], st2["Tt"], st2["Pt"]))
        # how flat is the optimum: BPR range within +0.5 % of min installed TSFC
        ok = [alphas[k] for k in range(len(alphas)) if cur["S_inst"][k] <= 1.005 * cur["opt"]["S_inst"]]
        cur["flat_range"] = [min(ok), max(ok)]
        inst["curves"].append(cur)
    out["install"] = inst

    # optimum BPR vs Tt4 and vs OPR (installed with baseline penalty; uninstalled too)
    pen = pens[1.0]

    def opt_bpr(b):
        def best(a, installed_flag=True):
            inp = b.with_(alpha=a)
            hi = pf_max(inp)

            def fn(pf):
                r = design(inp.with_(pi_f=pf))
                if not r["valid"]:
                    return 1e9
                v = installed(r, pen)["S_inst"] if installed_flag else r["S"]
                return v if math.isfinite(v) else 1e9
            pfo = golden(fn, 1.02, hi, 1e-7)
            return fn(pfo)
        grid = [1 + 0.5 * k for k in range(0, 70)]
        vals = [best(a) for a in grid]
        j = min(range(len(grid)), key=lambda k: vals[k])
        a_opt = golden(best, grid[max(j - 1, 0)], grid[min(j + 1, len(grid) - 1)], 1e-3)
        a_un = golden(lambda a: best(a, False), 5.0, 150.0, 1e-2)
        return a_opt, best(a_opt) * TSFC_LB, a_un, best(a_un, False) * TSFC_LB

    ob = dict(Tt4=[1400, 1500, 1600, 1700, 1800], alpha_inst=[], alpha_unin=[], S_inst=[], S_unin=[])
    for T4 in ob["Tt4"]:
        a1, s1, a2, s2 = opt_bpr(BASE.with_(Tt4=float(T4)))
        ob["alpha_inst"].append(a1); ob["S_inst"].append(s1); ob["alpha_unin"].append(a2); ob["S_unin"].append(s2)
    oo = dict(OPR=[20, 30, 40, 50, 60], alpha_inst=[], alpha_unin=[], S_inst=[], S_unin=[])
    for opr in oo["OPR"]:
        a1, s1, a2, s2 = opt_bpr(BASE.with_(pi_cH=opr / BASE.pi_cL))
        oo["alpha_inst"].append(a1); oo["S_inst"].append(s1); oo["alpha_unin"].append(a2); oo["S_unin"].append(s2)
    out["opt_bpr_vs_opr"] = oo
    out["opt_bpr_vs_tt4"] = ob
    print(f"trades: {time.time() - t0:.1f} s")
    return out


# ---------------------------------------------------------------------------
# Off-design
# ---------------------------------------------------------------------------
DTMAX = 150.0


def schedule(E, M0, alt):
    """Max-climb-like schedule: Tt4/theta0 held at design, capped at Tt4_des + 150 K."""
    inp = E.inp
    T0d, _, _, _ = atmosphere(inp.alt)
    th_d = T0d * (1 + 0.2 * inp.M0 ** 2) / 288.15
    T0, _, _, _ = atmosphere(alt)
    th = T0 * (1 + 0.2 * M0 ** 2) / 288.15
    return min(inp.Tt4 + DTMAX, inp.Tt4 * th / th_d)


def run_offdesign(cases):
    out = dict(DTmax=DTMAX, engines=[])
    for c in CASES[:3]:
        cs = next(x for x in cases if x["key"] == c["key"])
        E = Engine(c["inputs"], F_design=cs["sizing"]["F_cruise"])
        sls = E.operate(0.0, 0.0, schedule(E, 0.0, 0.0))
        eng = dict(key=c["key"], name=c["name"], F_des=E.F_des, F_sls=sls["F"], S_sls=sls["S_lb"],
                   Tt4_des=c["inputs"].Tt4, Tt4_max=c["inputs"].Tt4 + DTMAX, lapse_alt=[], lapse_mach=[])
        alts = [k * 500.0 for k in range(0, 27)]
        for M0 in [0.0, 0.3, 0.5, 0.78, 0.9]:
            L = dict(M0=M0, alt=[], F_rel=[], S_lb=[], Tt4=[], NL=[], alpha=[])
            for h in alts:
                T4 = schedule(E, M0, h)
                o = E.operate(M0, h, T4)
                if not o["valid"]:
                    continue
                L["alt"].append(h); L["F_rel"].append(o["F"] / sls["F"]); L["S_lb"].append(o["S_lb"])
                L["Tt4"].append(T4); L["NL"].append(o["NL_corr"]); L["alpha"].append(o["alpha"])
            eng["lapse_alt"].append(L)
        for h in [0.0, 5000.0, 35000 * FT]:
            L = dict(alt=h, M0=[], F_rel=[], S_lb=[], Tt4=[])
            for k in range(0, 19):
                M0 = 0.05 * k
                T4 = schedule(E, M0, h)
                o = E.operate(M0, h, T4)
                if not o["valid"]:
                    continue
                L["M0"].append(M0); L["F_rel"].append(o["F"] / sls["F"]); L["S_lb"].append(o["S_lb"]); L["Tt4"].append(T4)
            eng["lapse_mach"].append(L)
        # constant-Tt4 lapse at cruise Mach for contrast
        L = dict(alt=[], F_rel=[], S_lb=[])
        s0 = E.operate(c["inputs"].M0, 0.0, c["inputs"].Tt4)
        for h in alts:
            o = E.operate(c["inputs"].M0, h, c["inputs"].Tt4)
            if o["valid"]:
                L["alt"].append(h); L["F_rel"].append(o["F"] / s0["F"]); L["S_lb"].append(o["S_lb"])
        eng["lapse_constTt4"] = L
        # throttle hook at the design flight condition
        H = dict(Tt4=[], F_rel=[], S_lb=[], alpha=[], OPR=[], NL=[], eta_th=[], eta_p=[])
        T4 = 900.0
        while T4 <= c["inputs"].Tt4 + DTMAX + 1e-9:
            o = E.operate(c["inputs"].M0, c["inputs"].alt, T4)
            if o["valid"]:
                H["Tt4"].append(T4); H["F_rel"].append(o["F"] / E.F_des); H["S_lb"].append(o["S_lb"])
                H["alpha"].append(o["alpha"]); H["OPR"].append(o["OPR"]); H["NL"].append(o["NL_corr"])
                H["eta_th"].append(o["eta_th"]); H["eta_p"].append(o["eta_p"])
            T4 += 10.0
        eng["throttle"] = H
        jm = min(range(len(H["S_lb"])), key=lambda k: H["S_lb"][k])
        eng["throttle_min"] = dict(Tt4=H["Tt4"][jm], F_rel=H["F_rel"][jm], S_lb=H["S_lb"][jm])
        # cruise/SLS ratio at the schedule
        cr = E.operate(c["inputs"].M0, c["inputs"].alt, schedule(E, c["inputs"].M0, c["inputs"].alt))
        eng["cruise_over_sls"] = cr["F"] / sls["F"]
        eng["sls"] = dict(alpha=sls["alpha"], OPR=sls["OPR"], NL=sls["NL_corr"], Tt4=schedule(E, 0.0, 0.0),
                          mc2_rel=sls["mc2_rel"])
        out["engines"].append(eng)
    return out


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
def run_verification():
    v = {}
    # atmosphere
    ref = {0: 101325.0, 11000: 22632.06, 20000: 5474.889, 32000: 868.0187, 47000: 110.9063,
           51000: 66.93887, 71000: 3.956420}
    v["atmosphere"] = [dict(H=H, P=atmosphere(H)[1], P_ref=P, rel=abs(atmosphere(H)[1] - P) / P) for H, P in ref.items()]
    v["atm_35kft"] = dict(zip(["T", "P", "rho", "a"], atmosphere(35000 * FT)))
    v["thermo"] = dict(janaf_worst=T.test_thermo_data_vs_janaf(), inversion_worst=T.test_thermo_inversions_and_continuity(),
                       eps_cool=EPS_COOL)
    # ideal-cycle
    rows = []
    for s in T.IDEAL_SETS:
        inp = T.ideal_inputs(**s)
        r = design(inp)
        ic = ideal_turbofan(inp.M0, r["T0"], 1.4, 1004.0, inp.hPR, inp.Tt4, inp.OPR, inp.pi_f, inp.alpha)
        errs = {}
        for k in ["Fs", "S", "f", "V9", "V19", "eta_th", "eta_p"]:
            if k == "V19" and inp.alpha == 0:
                continue
            if k == "eta_p" and inp.M0 == 0:
                continue
            errs[k] = abs(r[k] - ic[k]) / abs(ic[k])
        rows.append(dict(set=s, Fs=r["Fs"], S_lb=r["S_lb"], errs=errs, worst=max(errs.values())))
    v["ideal"] = rows
    inp = T.ideal_inputs(M0=0.9, alt=10000.0, alpha=0.0, pi_f=1.0, pi_cL=4.0, pi_cH=6.0, Tt4=1600.0)
    r = design(inp)
    Fs, S = T.ideal_turbojet(inp.M0, r["T0"], 1.4, 1004.0, inp.hPR, inp.Tt4, inp.OPR)
    v["turbojet"] = dict(Fs=r["Fs"], Fs_ref=Fs, err=max(abs(r["Fs"] - Fs) / Fs, abs(r["S"] - S) / S))
    # energy & spool balance
    eb = []
    for c in CASES:
        r = design(c["inputs"])
        res, ref_ = energy_residual(r)
        i = c["inputs"]
        eb.append(dict(key=c["key"], energy=abs(res) / ref_,
                       hp=abs(i.eta_mH * r["W_HPT"] - r["W_HPC"]) / r["W_HPC"],
                       lp=abs(i.eta_mL * r["W_LPT"] - r["W_LP_load"]) / r["W_LP_load"],
                       W_HPT=r["W_HPT"], W_HPC=r["W_HPC"], W_LPT=r["W_LPT"], W_LP=r["W_LP_load"]))
    v["balances"] = eb
    # optimum fan pressure ratio: ideal (analytic) and real (approximate condition)
    opt_ideal = []
    for s in T.IDEAL_SETS[:3]:
        inp = T.ideal_inputs(**s)
        pf, r = optimum_fan_pr(inp, tol=1e-11)
        ic = ideal_turbofan(inp.M0, r["T0"], 1.4, 1004.0, inp.hPR, inp.Tt4, inp.OPR, 1.5, inp.alpha)
        opt_ideal.append(dict(set=s, pf_num=pf, pf_an=ic["pi_f_opt"], rel=abs(pf - ic["pi_f_opt"]) / ic["pi_f_opt"],
                              V19_V9=r["V19"] / r["V9"]))
    v["opt_ideal"] = opt_ideal
    # curves: TSFC vs pi_f for ideal and real at a few BPR, with optimum markers
    curves = []
    s0 = T.IDEAL_SETS[0]
    for a in [4.0, 8.0, 12.0]:
        ii = T.ideal_inputs(**dict(s0, alpha=a))
        rr = DesignInputs(M0=s0["M0"], alt=s0["alt"], alpha=a, pi_cL=s0["pi_cL"], pi_cH=s0["pi_cH"], Tt4=s0["Tt4"],
                          core_nozzle="expanded", fan_nozzle="expanded")
        C = dict(alpha=a, pf=[], S_ideal=[], S_real=[])
        for k in range(0, 181):
            pf = 1.1 + 0.025 * k
            a1, a2 = design(ii.with_(pi_f=pf)), design(rr.with_(pi_f=pf))
            C["pf"].append(pf)
            C["S_ideal"].append(a1["S_lb"] if a1["valid"] else None)
            C["S_real"].append(a2["S_lb"] if a2["valid"] else None)
        pfi, ri = optimum_fan_pr(ii, tol=1e-10)
        ic = ideal_turbofan(ii.M0, ri["T0"], 1.4, 1004.0, ii.hPR, ii.Tt4, ii.OPR, 1.5, a)
        pfr, rr_ = optimum_fan_pr(rr, tol=1e-10)
        pred = rr_["eta_ad"]["f"] * rr_["eta_ad"]["tL"] * rr.eta_mL
        exact = marginal_optimum_ratio(rr_)
        C.update(pred_exact=exact, err_exact=abs(rr_["V19"] / rr_["V9"] / exact - 1),
                 err_approx=abs(rr_["V19"] / rr_["V9"] / pred - 1), pf_opt_ideal_num=pfi, pf_opt_ideal_an=ic["pi_f_opt"], S_opt_ideal=ri["S_lb"],
                 pf_opt_real=pfr, S_opt_real=rr_["S_lb"], V19_V9_real=rr_["V19"] / rr_["V9"], pred_real=pred,
                 eta_f=rr_["eta_ad"]["f"], eta_tL=rr_["eta_ad"]["tL"])
        curves.append(C)
    v["opt_curves"] = curves
    v["opt_real_worst"] = max(C["err_approx"] for C in curves)
    v["opt_real_exact_worst"] = max(C["err_exact"] for C in curves)
    # burner: two-gas (Mattingly modified cycle) vs variable-property enthalpy balance
    bb = []
    air = Mixture(0.0)
    for c in CASES:
        i = c["inputs"]
        r = design(i)
        Tt3, Tt4 = r["stations"]["3"]["Tt"], i.Tt4
        fb_var = burner_far(air, Mixture, Tt3, Tt4, i.eta_b, i.hPR)
        fb_cpg = (i.cp_t * Tt4 - i.cp_c * Tt3) / (i.eta_b * i.hPR - i.cp_t * Tt4)
        dh_nasa_air = air.h(Tt4) - air.h(Tt3)
        dh_poly_air = dh_air_poly(Tt3, Tt4)
        bb.append(dict(key=c["key"], Tt3=Tt3, Tt4=Tt4, fb_varcp=fb_var, fb_cpg=fb_cpg, excess=fb_cpg / fb_var - 1,
                       dh_two_gas=i.cp_t * Tt4 - i.cp_c * Tt3, dh_air_nasa=dh_nasa_air, dh_air_poly=dh_poly_air,
                       air_poly_vs_nasa=dh_poly_air / dh_nasa_air - 1))
    v["burner_bias"] = bb
    # off-design reproduction and matching residuals
    rep = []
    for c in CASES[:3]:
        E = Engine(c["inputs"], F_design=20000.0)
        o = E.operate(c["inputs"].M0, c["inputs"].alt, c["inputs"].Tt4)
        rep.append(dict(key=c["key"], err=max(abs(o["Fs"] / E.des["Fs"] - 1), abs(o["S"] / E.des["S"] - 1),
                                              abs(o["alpha"] / c["inputs"].alpha - 1), abs(o["OPR"] / c["inputs"].OPR - 1))))
    v["od_reproduce"] = rep
    E = Engine(CASES[0]["inputs"], F_design=20000.0)
    worst, n = 0.0, 0
    for h in [0.0, 3000.0, 6000.0, 9000.0, 12000.0]:
        for M0 in [0.0, 0.3, 0.6, 0.85]:
            for T4 in [1250.0, 1400.0, 1550.0, 1650.0]:
                o = E.operate(M0, h, T4)
                if o["valid"]:
                    n += 1
                    worst = max(worst, abs(o["lp_residual"]), abs(o["mass_err_45"]), abs(o["mass_err_8"]))
    v["od_residual"] = dict(worst=worst, n=n)
    # 'Mattingly-style' real-cycle input set (textbook-like inputs, English-unit example style)
    ms = DesignInputs(M0=0.8, T0=216.67, P0=atmosphere(35000 * FT)[1], alpha=8.0, pi_f=1.7, pi_cL=1.7,
                      pi_cH=36 / 1.7, Tt4=1666.7, pi_d_max=0.99, pi_b=0.96, eta_b=0.99, e_f=0.89, e_cL=0.90,
                      e_cH=0.90, e_tH=0.89, e_tL=0.89, eta_mH=0.99, eta_mL=0.99, pi_n=0.99, pi_fn=0.99,
                      core_nozzle="expanded", fan_nozzle="expanded", thermo="cpg", eps_cool=0.0)
    r = design(ms)
    rv = design(ms.with_(thermo="varcp"))
    ii = T.ideal_inputs(M0=0.8, T0=216.67, alpha=8.0, pi_f=1.7, pi_cL=1.7, pi_cH=36 / 1.7, Tt4=1666.7)
    ri = design(ii)
    v["mattingly_style"] = dict(inputs=ms.to_dict(), Fs=r["Fs"], S_lb=r["S_lb"], eta_th=r["eta_th"],
                                eta_p=r["eta_p"], f=r["f"], Fs_ideal=ri["Fs"], S_ideal=ri["S_lb"],
                                Fs_lbf=r["Fs"] / 9.80665, Fs_ideal_lbf=ri["Fs"] / 9.80665,
                                Fs_varcp=rv["Fs"], S_varcp=rv["S_lb"])
    return v


def js_check(cases, od):
    node = shutil.which("node")
    if not node:
        return dict(ran=False, reason="node not found")
    p = subprocess.run([node, os.path.join(HERE, "js", "check.mjs"), DATA], capture_output=True, text=True, timeout=300)
    if p.returncode != 0:
        return dict(ran=False, reason=p.stderr[-2000:])
    res = json.loads(p.stdout)
    res["ran"] = True
    return res


def main():
    t0 = time.time()
    cases, comp = run_cases()
    dump("cases.json", dict(cases=cases, comparison=comp, TT4_BAND=TT4_BAND))
    with open(os.path.join(DATA, "cases_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "BPR", "pi_f", "OPR", "Tt4_K", "F_over_mdot0", "TSFC_lb_lbf_h",
                    "TSFC_lo_Tt4-100", "TSFC_hi_Tt4+100", "eta_th", "eta_p", "eta_o"])
        for c in cases:
            i, r = c["inputs"], c["result"]
            w.writerow([c["name"], i["alpha"], i["pi_f"], f'{i["pi_cL"] * i["pi_cH"]:.4g}', i["Tt4"],
                        f'{r["Fs"]:.5g}', f'{r["S_lb"]:.5g}', f'{c["band"]["S_lb"][0]:.5g}',
                        f'{c["band"]["S_lb"][1]:.5g}', f'{r["eta_th"]:.4g}', f'{r["eta_p"]:.4g}', f'{r["eta_o"]:.4g}'])
    od = run_offdesign(cases)
    dump("offdesign.json", od)
    tr = run_trades()
    dump("trades.json", tr, sig=8)
    v = run_verification()
    v["js_check"] = js_check(cases, od)
    dump("verification.json", v)
    print(f"done in {time.time() - t0:.1f} s")
    print("CFM56-class TSFC %.4f, LEAP-class %.4f, improvement %.1f %%" %
          (cases[0]["result"]["S_lb"], cases[1]["result"]["S_lb"], 100 * comp["leap_vs_cfm"]))
    print("installed optimum BPR (x0.5, x1, x2):", [round(c["opt"]["alpha"], 2) for c in tr["install"]["curves"]])
    print("JS check:", {k: v["js_check"][k] for k in v["js_check"] if k != "details"})


if __name__ == "__main__":
    main()
