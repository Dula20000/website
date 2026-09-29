"""Regenerate every result used by the nozzle-cfd page into ../data/*.json.

    /opt/homebrew/bin/python3 run_study.py

Writes: verification.json, engines.json, design_space.json, crosscheck.json, summary.json
(and grid_convergence.csv, engines_altitude.csv).
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from nozzlecfd import exact, geometry, performance as perf, solver  # noqa: E402
from nozzlecfd.atmosphere import atmosphere, P_SL  # noqa: E402

DATA = os.path.join(HERE, "..", "data")
os.makedirs(DATA, exist_ok=True)
TOL = 1e-11


def rnd(v, sig=7):
    """Round floats (recursively) to `sig` significant digits for compact JSON."""
    if isinstance(v, dict):
        return {k: rnd(x, sig) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [rnd(x, sig) for x in v]
    if isinstance(v, np.ndarray):
        return rnd(v.tolist(), sig)
    if isinstance(v, (np.floating, float)):
        v = float(v)
        if not np.isfinite(v) or v == 0:
            return v if np.isfinite(v) else None
        return float(f"{v:.{sig}g}")
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def dump(name, obj):
    with open(os.path.join(DATA, name), "w") as f:
        json.dump(rnd(obj), f, separators=(",", ":"))
    print("wrote", name)


def solve(geom, N, g, pb, limiter="vanalbada"):
    t = time.time()
    r = solver.NozzleSolver(geom, N, g, pb, limiter=limiter).run(max_iter=200000, tol=TOL, record_every=10)
    r.wall = time.time() - t
    return r


def order(e):
    return [None] + [float(np.log2(e[i - 1] / e[i])) for i in range(1, len(e))]


# =====================================================================
# 1. Verification
# =====================================================================
def verification():
    out = {}
    # ---- 1a isentropic grid convergence, Anderson nozzle, gamma 1.4 --------
    geom = geometry.AndersonNozzle()
    g, pb = 1.4, 0.01
    Ns = [50, 100, 200, 400, 800]
    rows = {"N": Ns, "dx": [geom.L / n for n in Ns]}
    for lim in ("vanalbada", "none"):
        eM, ep, einf, emd, iters, Mexit = [], [], [], [], [], []
        for N in Ns:
            r = solve(geom, N, g, pb, limiter=lim)
            ex = exact.exact_nozzle(r.x, geom, pb, g)
            eM.append(np.mean(np.abs(r.M - ex["M"])))
            ep.append(np.mean(np.abs(r.p - ex["p"])))
            einf.append(np.max(np.abs(r.M - ex["M"])))
            emd.append(abs(r.mdot_faces.mean() / exact.mdot_star(g) - 1))
            iters.append(r.iters)
            Mexit.append(r.exit_state()["M"])
        rows[lim] = dict(L1_M=eM, L1_p=ep, Linf_M=einf, mdot_err=emd, iters=iters, M_exit=Mexit,
                         order_L1_M=order(eM), order_L1_p=order(ep), order_mdot=order(emd))
    out["grid_isentropic"] = rows
    Me_ex = exact.mach_from_area(geom.eps, g, True)
    out["anderson"] = dict(eps=geom.eps, gamma=g, pb=pb, M_exit_exact=Me_ex,
                           p_exit_exact=float(exact.p_p0(Me_ex, g)), T_exit_exact=float(exact.T_T0(Me_ex, g)),
                           rho_exit_exact=float(exact.rho_rho0(Me_ex, g)), mdot_exact=float(exact.mdot_star(g)))
    # profile at N = 100 with table at a few stations
    r = solve(geom, 100, g, pb)
    ex = exact.exact_nozzle(r.x, geom, pb, g)
    out["profile_isentropic"] = dict(x=r.x, A=r.A, M=r.M, p=r.p, rho=r.rho, M_exact=ex["M"], p_exact=ex["p"],
                                     mdot=r.mdot_faces)
    xf = np.linspace(0, geom.L, 301)
    exf = exact.exact_nozzle(xf, geom, pb, g)
    out["profile_isentropic"]["x_fine"] = xf
    out["profile_isentropic"]["M_fine"] = exf["M"]
    out["profile_isentropic"]["p_fine"] = exf["p"]
    ex_state = r.exit_state()
    out["anderson"].update(M_exit_cfd100=ex_state["M"], p_exit_cfd100=ex_state["p"],
                           mdot_cfd100=float(r.mdot_faces.mean()),
                           mdot_spread100=float((r.mdot_faces.max() - r.mdot_faces.min()) / r.mdot_faces.mean()))

    # Richardson (exact-free) observed order on exit Mach for a rocket nozzle, gamma 1.2
    bell = geometry.BellLikeNozzle(16.0)
    NsR = [100, 200, 400, 800]
    MeR = [solve(bell, N, 1.2, 101325 / 9.7e6).exit_state()["M"] for N in NsR]
    pR = [float(np.log2((MeR[i] - MeR[i + 1]) / (MeR[i + 1] - MeR[i + 2]))) for i in range(len(NsR) - 2)]
    MeRich = MeR[-1] + (MeR[-1] - MeR[-2]) / (2 ** pR[-1] - 1)
    out["richardson_bell16"] = dict(N=NsR, M_exit=MeR, observed_order=pR, M_exit_extrapolated=MeRich,
                                    M_exit_exact=exact.mach_from_area(16.0, 1.2, True))

    # ---- 1b shock in divergent section ------------------------------------
    sgeom = geometry.AndersonShockNozzle()
    As_mid = float(sgeom.A(2.25))
    pb_mid = exact.back_pressure_for_shock_at(As_mid, sgeom.eps, g)
    crit = exact.critical_back_pressures(sgeom.eps, g)
    Ns2 = [50, 100, 200, 400, 800]
    xs_err, md_err, md_spread, xs_cfd, pe_err, it2 = [], [], [], [], [], []
    for N in Ns2:
        r = solve(sgeom, N, g, pb_mid)
        exs = exact.exact_nozzle(r.x, sgeom, pb_mid, g)
        xs = solver.shock_location(r, sgeom.x_throat)
        xs_cfd.append(xs)
        xs_err.append(abs(xs - exs["x_shock"]))
        md_err.append(abs(r.mdot_faces.mean() / exact.mdot_star(g) - 1))
        md_spread.append((r.mdot_faces.max() - r.mdot_faces.min()) / r.mdot_faces.mean())
        pe_err.append(abs(r.exit_state()["p"] / pb_mid - 1))
        it2.append(r.iters)
        if N == 200:
            out["profile_shock"] = dict(x=r.x, M=r.M, p=r.p, M_exact=exs["M"], p_exact=exs["p"])
            xf = np.linspace(0, sgeom.L, 601)
            exf = exact.exact_nozzle(xf, sgeom, pb_mid, g)
            out["profile_shock"].update(x_fine=xf, M_fine=exf["M"], p_fine=exf["p"])
    out["shock_case"] = dict(eps=sgeom.eps, gamma=g, pb=pb_mid, x_shock_exact=2.25, A_shock=As_mid, crit=crit,
                             N=Ns2, dx=[sgeom.L / n for n in Ns2], x_shock_cfd=xs_cfd, x_shock_err=xs_err,
                             x_shock_err_cells=[e / (sgeom.L / n) for e, n in zip(xs_err, Ns2)],
                             order_x_shock=order(xs_err), mdot_err=md_err, mdot_spread=md_spread,
                             p_exit_err=pe_err, iters=it2)
    # Anderson's pb/p0 = 0.6784 back pressure on the same nozzle
    r = solve(sgeom, 200, g, 0.6784)
    exa = exact.exact_nozzle(r.x, sgeom, 0.6784, g)
    out["anderson_shock"] = dict(pb=0.6784, N=200, x_shock_exact=exa["x_shock"],
                                 x_shock_cfd=solver.shock_location(r, sgeom.x_throat),
                                 A_shock_exact=exa["A_shock"])
    # sweep of back pressure: shock station vs pb
    pbs = np.linspace(crit["p_nse"] + 0.01 * (crit["p_sub"] - crit["p_nse"]),
                      crit["p_sub"] - 0.02 * (crit["p_sub"] - crit["p_nse"]), 9)
    sw = dict(pb=pbs, x_cfd=[], x_exact=[])
    for p in pbs:
        r = solve(sgeom, 200, g, float(p))
        sw["x_cfd"].append(solver.shock_location(r, sgeom.x_throat))
        sw["x_exact"].append(exact.exact_nozzle(r.x[:2], sgeom, float(p), g)["x_shock"])
    pf = np.linspace(crit["p_nse"], crit["p_sub"], 200)
    sw["pb_fine"] = pf
    sw["x_exact_fine"] = [exact.exact_nozzle(np.array([0.0]), sgeom, float(p), g)["x_shock"] for p in pf]
    out["shock_sweep"] = sw

    # ---- 1c convergence histories ------------------------------------------
    hist = {}
    for key, (gm, N, gg, p) in {
        "isentropic_anderson_N200": (geom, 200, 1.4, 0.01),
        "shock_N200": (sgeom, 200, 1.4, pb_mid),
        "merlin_N400": (geometry.BellLikeNozzle(16.0), 400, 1.2, 101325 / 9.7e6),
    }.items():
        r = solve(gm, N, gg, p)
        h = np.array(r.residual)
        hist[key] = dict(iter=h[:, 0].astype(int).tolist(), res=h[:, 1], iters=r.iters, wall_s=r.wall)
    out["residual_history"] = hist
    return out


# =====================================================================
# 2. Engines along an ascent
# =====================================================================
ENGINES = [
    dict(key="merlin", name="Merlin 1D-class (sea level)", Pc=9.7e6, eps=16.0, gamma=1.2,
         prop="LOX / RP-1", source="Manufacturer public specifications and press material (approx.)",
         note="Pc ≈ 9.7 MPa, ε ≈ 16"),
    dict(key="rs25", name="RS-25-class (SSME)", Pc=20.6e6, eps=69.0, gamma=1.2,
         prop="LOX / LH2", source="NASA public fact sheets; Sutton & Biblarz engine tables (approx.)",
         note="Pc ≈ 20.6 MPa (109 % power level), ε ≈ 69"),
    dict(key="raptor", name="Raptor-class (sea level)", Pc=30.0e6, eps=40.0, gamma=1.2,
         prop="LOX / CH4", source="Widely cited public figures; not officially published in detail",
         note="Pc ≈ 30 MPa, ε ≈ 40 (uncertain; values from ≈ 34 to 40 are quoted)"),
]
N_ENGINE = 400


def engines():
    z = np.linspace(0.0, 50e3, 201)
    out = dict(z_note="geometric altitude, US Standard Atmosphere 1976", N=N_ENGINE,
               summerfield=perf.SUMMERFIELD, p_sl=P_SL, engines=[])
    csv_rows = []
    for e in ENGINES:
        g, eps, Pc = e["gamma"], e["eps"], e["Pc"]
        geom = geometry.BellLikeNozzle(eps)
        pa_sl = P_SL
        r = solve(geom, N_ENGINE, g, pa_sl / Pc)
        st = r.exit_state()
        Me_ex = exact.mach_from_area(eps, g, True)
        pe_ex = float(exact.p_p0(Me_ex, g))
        mom = float(r.mom_faces[-1])
        # exit state invariance: rerun at 10 km and 30 km ambient
        inv = []
        for zz in (10e3, 30e3):
            r2 = solve(geom, N_ENGINE, g, atmosphere(zz)[1] / Pc)
            inv.append(abs(r2.mom_faces[-1] - mom) / mom)
        pe = st["p"] * Pc
        sweep = perf.ascent_sweep(st["p"], st["M"], mom, eps, Pc, g, z)
        z_opt = perf.optimal_altitude(lambda zz: pe / atmosphere(zz)[1])
        z_sf = perf.altitude_for_ratio(pe, perf.SUMMERFIELD)
        sch = float(perf.schmucker_ratio(st["M"]))
        z_sch = perf.altitude_for_ratio(pe, sch)
        c = exact.critical_back_pressures(eps, g)
        As_mid = 0.5 * (1 + eps)
        p_mid = exact.back_pressure_for_shock_at(As_mid, eps, g)
        cf_sl = mom - pa_sl / Pc * eps
        cf_vac = mom
        # CFD check of a shock inside this nozzle (hypothetical back pressure placing it at A = (1+eps)/2)
        rs = solve(geom, N_ENGINE, g, p_mid)
        xs_cfd = solver.shock_location(rs, geom.x_throat)
        As_cfd = float(geom.A(xs_cfd)) if xs_cfd is not None else None
        # Mach profile (downsampled) for plotting
        k = max(1, N_ENGINE // 100)
        prof = dict(x=(r.x[::k] / geom.L), r=np.sqrt(r.A[::k]), M=r.M[::k], p=r.p[::k])
        # exact shock-location curve vs back pressure
        pbs = np.geomspace(c["p_nse"], c["p_sub"], 120)
        As_curve = [exact.shock_area_ratio(float(p), eps, g)[0] for p in pbs]
        # gamma sensitivity of the sea-level expansion state
        gsens = []
        for gg in (1.15, 1.2, 1.25, 1.3):
            Mg = exact.mach_from_area(eps, gg, True)
            peg = float(exact.p_p0(Mg, gg)) * Pc
            gsens.append(dict(gamma=gg, Me=Mg, pe_kPa=peg / 1e3, pe_pa_sl=peg / pa_sl,
                              schmucker=float(perf.schmucker_ratio(Mg))))
        rec = dict(
            key=e["key"], name=e["name"], prop=e["prop"], source=e["source"], note=e["note"],
            Pc=Pc, eps=eps, gamma=g,
            cfd=dict(M_exit=st["M"], pe_p0=st["p"], mdot=float(r.mdot_faces[-1]), mom_exit=mom,
                     iters=r.iters, converged=r.converged, exit_supersonic=r.exit_supersonic,
                     invariance_rel=max(inv)),
            exact=dict(M_exit=Me_ex, pe_p0=pe_ex, cf_vac=float(exact.cf_ideal(pe_ex, 0.0, eps, g)),
                       cf_sl=float(exact.cf_ideal(pe_ex, pa_sl / Pc, eps, g))),
            err=dict(M_exit=abs(st["M"] / Me_ex - 1), pe=abs(st["p"] / pe_ex - 1),
                     cf_vac=abs(cf_vac / float(exact.cf_ideal(pe_ex, 0.0, eps, g)) - 1),
                     mdot=abs(float(r.mdot_faces[-1]) / exact.mdot_star(g) - 1)),
            pe_kPa=pe / 1e3, pe_pa_sl=pe / pa_sl, cf_sl=cf_sl, cf_vac=cf_vac,
            z_opt_km=None if z_opt is None else z_opt / 1e3,
            z_summerfield_clear_km=None if z_sf is None else z_sf / 1e3,
            schmucker_ratio=sch, z_schmucker_clear_km=None if z_sch is None else z_sch / 1e3,
            summerfield_risk_sl=pe / pa_sl < perf.SUMMERFIELD, schmucker_risk_sl=pe / pa_sl < sch,
            p_nse_MPa=c["p_nse"] * Pc / 1e6, p_mid_MPa=p_mid * Pc / 1e6, p_sub_MPa=c["p_sub"] * Pc / 1e6,
            p_nse_over_pa_sl=c["p_nse"] * Pc / pa_sl, p_sub_over_pa_sl=c["p_sub"] * Pc / pa_sl,
            p_sep_summerfield_kPa=pe / perf.SUMMERFIELD / 1e3,
            shock_check=dict(pb_p0=p_mid, As_exact=As_mid, As_cfd=As_cfd,
                             As_err_rel=None if As_cfd is None else abs(As_cfd / As_mid - 1)),
            pc_min_summerfield_MPa=perf.pc_min_summerfield(eps, pa_sl, g) / 1e6,
            eps_max_summerfield=perf.eps_max_summerfield(Pc, pa_sl, g),
            eps_max_schmucker=perf.eps_max_schmucker(Pc, pa_sl, g),
            eps_opt_sl=perf.eps_optimal(Pc, pa_sl, g),
            sweep=sweep, profile=prof,
            shock_curve=dict(pb_over_pa_sl=(pbs * Pc / pa_sl), As_over_eps=np.array(As_curve) / eps),
            gamma_sensitivity=gsens,
        )
        rec["throttle_margin"] = rec["pc_min_summerfield_MPa"] * 1e6 / Pc
        rec["eps_over_eps_max"] = eps / rec["eps_max_summerfield"]
        # same nozzle at the Merlin-class chamber pressure (isolates the effect of Pc)
        rec["pe_pa_sl_at_ref_pc"] = st["p"] * ENGINES[0]["Pc"] / pa_sl
        out["engines"].append(rec)
        for i, zz in enumerate(sweep["z_km"]):
            csv_rows.append([e["key"], zz, sweep["pa"][i], sweep["pe_pa"][i], sweep["cf"][i],
                             int(sweep["summerfield_risk"][i])])
    with open(os.path.join(DATA, "engines_altitude.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["engine", "z_km", "pa_Pa", "pe_over_pa", "CF", "summerfield_risk"])
        w.writerows(csv_rows)
    return out


def design_space():
    g = 1.2
    Pc = np.geomspace(2e6, 40e6, 120)
    out = dict(gamma=g, pa=P_SL, Pc_MPa=Pc / 1e6,
               eps_summerfield=[perf.eps_max_summerfield(p, P_SL, g) for p in Pc],
               eps_schmucker=[perf.eps_max_schmucker(p, P_SL, g) for p in Pc],
               eps_optimal=[perf.eps_optimal(p, P_SL, g) for p in Pc],
               engines=[dict(key=e["key"], name=e["name"], Pc_MPa=e["Pc"] / 1e6, eps=e["eps"]) for e in ENGINES])
    # gamma band for the Summerfield line
    out["eps_summerfield_g115"] = [perf.eps_max_summerfield(p, P_SL, 1.15) for p in Pc]
    out["eps_summerfield_g125"] = [perf.eps_max_summerfield(p, P_SL, 1.25) for p in Pc]
    return out


# =====================================================================
# 3. JS port cross-check
# =====================================================================
def crosscheck():
    node = shutil.which("node")
    cases = [
        dict(name="Anderson isentropic, N=100", geom="anderson", eps=5.95, N=100, gamma=1.4, pb=0.01),
        dict(name="Merlin-class, N=200, sea level", geom="bell", eps=16.0, N=200, gamma=1.2, pb=P_SL / 9.7e6),
        dict(name="Merlin-class nozzle, shock at A=8.5", geom="bell", eps=16.0, N=200, gamma=1.2,
             pb=exact.back_pressure_for_shock_at(8.5, 16.0, 1.2)),
    ]
    for c in cases:
        c.update(cfl=0.8, max_iter=200000, tol=TOL)
    if node is None:
        return dict(available=False, cases=[])
    js = json.loads(subprocess.run([node, os.path.join(HERE, "crosscheck_js.cjs"), json.dumps(cases)],
                                   capture_output=True, text=True, check=True).stdout)
    res = []
    for c, j in zip(cases, js):
        geom = geometry.AndersonNozzle() if c["geom"] == "anderson" else geometry.BellLikeNozzle(c["eps"])
        r = solve(geom, c["N"], c["gamma"], c["pb"])
        dM = float(np.max(np.abs(np.array(j["M"]) - r.M)))
        dp = float(np.max(np.abs(np.array(j["p"]) - r.p) / r.p))
        xs_py = solver.shock_location(r, geom.x_throat)
        res.append(dict(name=c["name"], N=c["N"], iters_py=r.iters, iters_js=j["iters"], max_dM=dM, max_rel_dp=dp,
                        mdot_py=float(r.mdot_faces[-1]), mdot_js=j["mdot"],
                        cf_vac_py=float(r.mom_faces[-1]), cf_vac_js=j["mom_exit"],
                        x_shock_py=xs_py, x_shock_js=j["x_shock"], js_ms=j["ms"], py_s=r.wall))
    return dict(available=True, node=subprocess.run([node, "--version"], capture_output=True, text=True).stdout.strip(),
                cases=res, max_dM=max(c["max_dM"] for c in res))


def main():
    t0 = time.time()
    ver = verification()
    dump("verification.json", ver)
    with open(os.path.join(DATA, "grid_convergence.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["N", "L1_M_muscl", "order", "L1_M_first_order", "order_first", "mdot_err_muscl"])
        gi = ver["grid_isentropic"]
        for i, N in enumerate(gi["N"]):
            w.writerow([N, gi["vanalbada"]["L1_M"][i], gi["vanalbada"]["order_L1_M"][i],
                        gi["none"]["L1_M"][i], gi["none"]["order_L1_M"][i], gi["vanalbada"]["mdot_err"][i]])
    eng = engines()
    dump("engines.json", eng)
    ds = design_space()
    dump("design_space.json", ds)
    xc = crosscheck()
    dump("crosscheck.json", xc)
    gi = ver["grid_isentropic"]["vanalbada"]
    sc = ver["shock_case"]
    E = {e["key"]: e for e in eng["engines"]}
    summary = dict(
        order_L1_M_finest=gi["order_L1_M"][-1],
        order_L1_M_mean=float(np.mean(gi["order_L1_M"][1:])),
        L1_M_finest=gi["L1_M"][-1], N_finest=ver["grid_isentropic"]["N"][-1],
        shock_err_cells_N200=sc["x_shock_err_cells"][2], shock_err_N200=sc["x_shock_err"][2],
        mdot_spread_max=max(sc["mdot_spread"]),
        rs25_pe_pa_sl=E["rs25"]["pe_pa_sl"], merlin_pe_pa_sl=E["merlin"]["pe_pa_sl"],
        raptor_pe_pa_sl=E["raptor"]["pe_pa_sl"],
        rs25_pc_min_MPa=E["rs25"]["pc_min_summerfield_MPa"],
        merlin_z_opt_km=E["merlin"]["z_opt_km"], rs25_z_opt_km=E["rs25"]["z_opt_km"],
        raptor_z_opt_km=E["raptor"]["z_opt_km"],
        min_pnse_over_pa=min(e["p_nse_over_pa_sl"] for e in eng["engines"]),
        js_max_dM=xc.get("max_dM"),
        runtime_s=time.time() - t0,
    )
    dump("summary.json", summary)
    print("done in %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
