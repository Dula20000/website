"""Verification tests for tfcycle.  Run with ``python3 -m pytest`` or directly with
``python3 tests/test_cycle.py`` (no pytest needed)."""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tfcycle import atmosphere, DesignInputs, design, ideal_turbofan, optimum_fan_pr, Engine  # noqa: E402
from tfcycle.ondesign import energy_residual  # noqa: E402

FT = 0.3048


def ideal_inputs(**kw):
    """Real-cycle inputs collapsed to the ideal cycle: all efficiencies and pressure
    ratios = 1, one gas, ideally expanded nozzles, fuel mass neglected."""
    base = dict(gamma_c=1.4, cp_c=1004.0, gamma_t=1.4, cp_t=1004.0, pi_d_max=1.0, pi_b=1.0,
                eta_b=1.0, e_f=1.0, e_cL=1.0, e_cH=1.0, e_tH=1.0, e_tL=1.0, eta_mH=1.0,
                eta_mL=1.0, pi_n=1.0, pi_fn=1.0, core_nozzle="expanded", fan_nozzle="expanded",
                neglect_fuel_mass=True, mil_spec_recovery=False, thermo="cpg", eps_cool=0.0)
    base.update(kw)
    return DesignInputs(**base)


IDEAL_SETS = [
    dict(M0=0.8, alt=35000 * FT, alpha=8.0, pi_f=1.7, pi_cL=1.7, pi_cH=36 / 1.7, Tt4=1666.7),
    dict(M0=0.0, alt=0.0, alpha=5.0, pi_f=1.6, pi_cL=2.0, pi_cH=15.0, Tt4=1500.0),
    dict(M0=0.85, alt=11000.0, alpha=12.0, pi_f=1.4, pi_cL=1.5, pi_cH=30.0, Tt4=1700.0),
    dict(M0=2.0, alt=15000.0, alpha=0.5, pi_f=3.0, pi_cL=3.0, pi_cH=5.0, Tt4=1900.0),
    dict(M0=0.5, alt=3000.0, alpha=0.0, pi_f=1.0, pi_cL=1.0, pi_cH=20.0, Tt4=1400.0),
]


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def test_atmosphere_layer_bases():
    # 1976 US Standard Atmosphere tabulated base pressures [Pa] (geopotential km)
    ref = {0: 101325.0, 11000: 22632.06, 20000: 5474.889, 32000: 868.0187,
           47000: 110.9063, 51000: 66.93887, 71000: 3.956420}
    for H, P in ref.items():
        assert rel(atmosphere(H)[1], P) < 2e-6, (H, atmosphere(H)[1], P)
    assert abs(atmosphere(11000)[0] - 216.65) < 1e-9
    assert abs(atmosphere(0)[3] - 340.294) < 1e-3


def test_ideal_cycle_matches_closed_form():
    worst = 0.0
    for s in IDEAL_SETS:
        inp = ideal_inputs(**s)
        r = design(inp)
        assert r["valid"]
        T0 = r["T0"]
        ic = ideal_turbofan(inp.M0, T0, 1.4, 1004.0, inp.hPR, inp.Tt4, inp.OPR, inp.pi_f, inp.alpha)
        keys = ["Fs", "S", "f", "V9"] + (["V19"] if inp.alpha > 0 else [])
        if inp.M0 > 0:
            keys += ["eta_th", "eta_p"]
        else:
            assert rel(r["eta_th"], ic["eta_th"]) < 1e-12
        for k in keys:
            e = rel(r[k], ic[k])
            worst = max(worst, e)
            assert e < 1e-12, (s, k, r[k], ic[k])
    return worst


def ideal_turbojet(M0, T0, gamma, cp, hPR, Tt4, pi_c):
    """Independent closed form for the ideal turbojet (Mattingly ch. 5)."""
    R = (gamma - 1) / gamma * cp
    a0 = math.sqrt(gamma * R * T0)
    tr = 1 + 0.5 * (gamma - 1) * M0 ** 2
    tl = Tt4 / T0
    tc = pi_c ** ((gamma - 1) / gamma)
    tt = 1 - tr / tl * (tc - 1)
    V9a = math.sqrt(2 / (gamma - 1) * tl / (tr * tc) * (tr * tc * tt - 1))
    Fs = a0 * (V9a - M0)
    f = cp * T0 / hPR * (tl - tr * tc)
    return Fs, f / Fs


def test_ideal_turbojet_limit():
    inp = ideal_inputs(M0=0.9, alt=10000.0, alpha=0.0, pi_f=1.0, pi_cL=4.0, pi_cH=6.0, Tt4=1600.0)
    r = design(inp)
    Fs, S = ideal_turbojet(inp.M0, r["T0"], 1.4, 1004.0, inp.hPR, inp.Tt4, inp.OPR)
    assert rel(r["Fs"], Fs) < 1e-12 and rel(r["S"], S) < 1e-12


def real_cases():
    return [
        DesignInputs(),
        DesignInputs(alpha=11, pi_f=1.45, pi_cL=1.82, pi_cH=22, Tt4=1600),
        DesignInputs(M0=0.9, alpha=0.36, pi_f=3.8, pi_cL=3.8, pi_cH=8.4, Tt4=1550),
        DesignInputs(M0=0.0, alt=0.0, alpha=8, pi_f=1.55, pi_cL=1.83, pi_cH=23, Tt4=1700,
                     core_nozzle="expanded", fan_nozzle="expanded"),
        DesignInputs(thermo="cpg", eps_cool=0.0, Tt4=1450),
        DesignInputs(thermo="cpg", eps_cool=0.15, Tt4=1550),
        DesignInputs(eps_cool=0.0, Tt4=1450),
    ]


# JANAF-table reference values (J/(mol K) and kJ/mol), rounded as tabulated
JANAF = {  # species: (cp 298.15, s 298.15, cp 1000, H(1000)-H(298.15))
    "N2": (29.124, 191.609, 32.698, 21.460),
    "O2": (29.376, 205.147, 34.870, 22.707),
    "Ar": (20.786, 154.845, 20.786, 14.577),
    "CO2": (37.129, 213.795, 54.308, 33.397),
    "H2O": (33.588, 188.834, 41.268, 25.978),
}


def test_thermo_data_vs_janaf():
    """NASA-polynomial species data reproduce JANAF cp, s and sensible enthalpy
    within 0.3 %; formation enthalpies of CO2 and H2O within 0.01 %."""
    from tfcycle.thermo import species_props, NASA, _h_RT, RU
    worst = 0.0
    for k, (cp298, s298, cp1000, dh1000) in JANAF.items():
        c1, _, s1 = species_props(k, 298.15)
        c2, h2, _ = species_props(k, 1000.0)
        for a, b in [(c1 / 1e3, cp298), (s1 / 1e3, s298), (c2 / 1e3, cp1000), (h2 / 1e6, dh1000)]:
            worst = max(worst, rel(a, b))
    assert worst < 3e-3, worst
    for k, hf in [("CO2", -393.522), ("H2O", -241.826)]:
        v = _h_RT(NASA[k][1], 298.15) * 298.15 * RU / 1e6
        assert rel(v, hf) < 1e-4, (k, v)
    return worst


def test_thermo_inversions_and_continuity():
    """T(h), T(s0) invert h, s0 to round-off, and h, s0 are continuous at 1000 K."""
    from tfcycle.thermo import Mixture
    worst = 0.0
    for f in [0.0, 0.02, 0.05]:
        g = Mixture(f)
        for T in [220.0, 600.0, 999.9, 1000.0, 1000.1, 1600.0, 2200.0]:
            worst = max(worst, rel(g.T_from_h(g.h(T)), T), rel(g.T_from_s0(g.s0(T)), T))
        assert abs(g.h(1000 + 1e-9) - g.h(1000 - 1e-9)) < 1e-3
        assert abs(g.s0(1000 + 1e-9) - g.s0(1000 - 1e-9)) < 1e-6
    assert worst < 1e-12, worst
    return worst


def test_burner_balance_varcp():
    """Burner energy balance h_air(Tt3) + f_b eta_b hPR = (1+f_b) h_prod(Tt4, f_b)
    closes, and the result differs from the two-gas formula by the known bias."""
    from tfcycle.thermo import Mixture
    from tfcycle.ondesign import burner_far
    air = Mixture(0.0)
    fb = burner_far(air, Mixture, 750.0, 1600.0, 0.995, 42.8e6)
    res = air.h(750.0) + fb * 0.995 * 42.8e6 - (1 + fb) * Mixture(fb).h(1600.0)
    assert abs(res) < 1e-8 * air.h(750.0)
    f_cpg = (1156 * 1600 - 1004 * 750) / (0.995 * 42.8e6 - 1156 * 1600)
    assert f_cpg > 1.05 * fb      # two-gas model overestimates fuel flow


def test_energy_and_spool_balances():
    worst = 0.0
    for inp in real_cases():
        r = design(inp)
        assert r["valid"]
        res, ref = energy_residual(r)
        e1 = abs(res) / ref
        e2 = abs(inp.eta_mH * r["W_HPT"] - r["W_HPC"]) / r["W_HPC"]
        e3 = abs(inp.eta_mL * r["W_LPT"] - r["W_LP_load"]) / r["W_LP_load"]
        worst = max(worst, e1, e2, e3)
        assert e1 < 1e-13 and e2 < 1e-13 and e3 < 1e-13, (e1, e2, e3)
    return worst


def test_ideal_optimum_fan_pressure_ratio():
    worst = 0.0
    for s in IDEAL_SETS[:3]:
        inp = ideal_inputs(**s)
        pf, r = optimum_fan_pr(inp, tol=1e-11)
        ic = ideal_turbofan(inp.M0, r["T0"], 1.4, 1004.0, inp.hPR, inp.Tt4, inp.OPR, 1.5, inp.alpha)
        e = rel(pf, ic["pi_f_opt"])
        worst = max(worst, e)
        assert e < 1e-6, (pf, ic["pi_f_opt"])
        assert rel(r["V19"], r["V9"]) < 1e-5     # equal-jet-velocity condition
    return worst


def test_real_optimum_velocity_ratio():
    """Real cycle, expanded nozzles: at the TSFC-optimal fan pressure ratio the jet
    velocity ratio V19/V9 is close to eta_f * eta_tL * eta_mL (Mattingly's
    approximate optimum condition).  Tolerance 5 %: the condition is approximate."""
    inp = DesignInputs(alpha=8, pi_cL=1.8, pi_cH=22, Tt4=1600,
                       core_nozzle="expanded", fan_nozzle="expanded")
    pf, r = optimum_fan_pr(inp)
    pred = r["eta_ad"]["f"] * r["eta_ad"]["tL"] * inp.eta_mL
    assert rel(r["V19"] / r["V9"], pred) < 0.05


def test_real_optimum_exact_marginal_condition():
    """The TSFC-optimal pi_f found numerically satisfies the exact first-order
    condition derived for the real cycle (ondesign.marginal_optimum_ratio)."""
    from tfcycle.ondesign import marginal_optimum_ratio
    worst = 0.0
    for a, th in [(4.0, "varcp"), (8.0, "varcp"), (12.0, "varcp"), (8.0, "cpg")]:
        inp = DesignInputs(alpha=a, pi_cL=1.8, pi_cH=22, Tt4=1600, thermo=th,
                           core_nozzle="expanded", fan_nozzle="expanded")
        pf, r = optimum_fan_pr(inp, tol=1e-11)
        e = rel(r["V19"] / r["V9"], marginal_optimum_ratio(r))
        worst = max(worst, e)
        assert e < 1e-5, (a, e)
    return worst


def test_offdesign_reproduces_design():
    worst = 0.0
    for inp in real_cases()[:3] + real_cases()[4:]:
        E = Engine(inp, F_design=20000.0)
        o = E.operate(inp.M0, inp.alt, inp.Tt4)
        r = dict(E.des, alpha=inp.alpha, pi_f=inp.pi_f)
        for k in ["Fs", "S", "alpha", "pi_f", "f", "V9", "V19"]:
            e = rel(o[k], r[k])
            worst = max(worst, e)
            assert e < 1e-10, (k, o[k], r[k])
        assert rel(o["OPR"], inp.OPR) < 1e-9 and rel(o["F"], 20000.0) < 1e-9
    return worst


def test_offdesign_matching_residuals():
    E = Engine(DesignInputs(), F_design=20000.0)
    worst = 0.0
    for alt in [0.0, 6000.0, 12000.0]:
        for M0 in [0.0, 0.4, 0.8]:
            for Tt4 in [1300.0, 1500.0, 1650.0]:
                o = E.operate(M0, alt, Tt4)
                assert o["valid"], (alt, M0, Tt4)
                w = max(abs(o["lp_residual"]), abs(o["mass_err_45"]), abs(o["mass_err_8"]))
                worst = max(worst, w)
                assert w < 1e-10, (alt, M0, Tt4, o["lp_residual"], o["mass_err_45"], o["mass_err_8"])
    return worst


def test_offdesign_trends():
    """Physical sanity: at fixed Tt4 and Mach, thrust falls with altitude below the
    tropopause; at fixed flight condition thrust rises with Tt4."""
    E = Engine(DesignInputs(), F_design=20000.0)
    Fs = [E.operate(0.5, h, 1500.0)["F"] for h in [0, 2000, 4000, 6000, 8000, 10000]]
    assert all(b < a for a, b in zip(Fs, Fs[1:]))
    Ft = [E.operate(0.78, 35000 * FT, T)["F"] for T in [1300, 1400, 1500, 1600]]
    assert all(b > a for a, b in zip(Ft, Ft[1:]))


if __name__ == "__main__":
    n = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            out = fn()
            n += 1
            print(f"PASS {name}" + (f"  (worst rel. error {out:.2e})" if isinstance(out, float) else ""))
    print(f"{n} tests passed")
