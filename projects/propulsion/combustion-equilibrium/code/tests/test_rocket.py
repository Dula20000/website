"""Verification of the nozzle / performance layer.

* Exact solution: a calorically perfect monatomic gas (argon, cp = 5/2 R
  exactly in the NASA fit) must reproduce the closed-form ideal-rocket
  relations (c*, P_e/P_c, C_F) with γ = 5/3 to solver tolerance.
* Ideal-rocket consistency for real propellants: c* ≈ sqrt(γ R T_c / M)/Γ(γ).
* Frozen Isp < shifting Isp, by a few per cent.
* Published CEA-type theoretical figures (only coarse, well-known ones).
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from eqsolver.equilibrium import GasSystem  # noqa: E402
from eqsolver.rocket import (Reactant, mix, chamber_from_mix, nozzle, performance,  # noqa: E402
                             chamber, G0)
from eqsolver.thermo import SPECIES, R  # noqa: E402


def big_gamma(g):
    return math.sqrt(g) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))


def ideal_rocket(g, T0, M, eps):
    """Closed-form ideal (calorically perfect) rocket: returns c*, Pe/Pc, CF_vac."""
    cstar = math.sqrt(R / M * T0) / big_gamma(g)
    Me = 3.0
    for _ in range(100):
        t = 1 + 0.5 * (g - 1) * Me ** 2
        e = (1 / Me) * ((2 / (g + 1)) * t) ** ((g + 1) / (2 * (g - 1)))
        dM = (math.log(e) - math.log(eps)) / (-1 / Me + Me / t)
        Me -= dM
        if abs(dM) < 1e-15:
            break
    pr = (1 + 0.5 * (g - 1) * Me ** 2) ** (-g / (g - 1))
    cf = big_gamma(g) * math.sqrt(2 * g / (g - 1) * (1 - pr ** ((g - 1) / g))) + pr * eps
    return cstar, pr, cf


def argon_exact(T0=3000.0, Pc=5e6, eps=10.0):
    ar = Reactant("Ar-hot", {"Ar": 1}, SPECIES["Ar"].h(T0), T0)
    el, b0, h0 = mix([(ar, 1.0)])
    ch = chamber_from_mix(el, b0, h0, Pc, system=GasSystem(["Ar"], el))
    out = {}
    for frozen in (False, True):
        p = nozzle(ch, eps, frozen=frozen)
        cs, pr, cf = ideal_rocket(5 / 3, T0, SPECIES["Ar"].molar_mass * 1e-3, eps)
        out["frozen" if frozen else "shifting"] = {
            "Tc": ch.T, "cstar": p.cstar, "cstar_exact": cs,
            "pr": p.Pe / Pc, "pr_exact": pr, "CF": p.CF_vac, "CF_exact": cf,
            "err_cstar": abs(p.cstar / cs - 1), "err_pr": abs(p.Pe / Pc / pr - 1),
            "err_CF": abs(p.CF_vac / cf - 1), "err_Tc": abs(ch.T / T0 - 1)}
    return out


def test_argon_exact():
    for mode, r in argon_exact().items():
        for k in ("err_cstar", "err_pr", "err_CF", "err_Tc"):
            assert r[k] < 1e-8, (mode, k, r)


def cstar_consistency(fuel="LCH4", ox="LOX", of=3.6, Pc=10e6):
    ch = chamber(fuel, ox, of, Pc)
    out = {}
    for frozen in (False, True):
        p = nozzle(ch, 40, frozen=frozen)
        g = ch.gamma_fr if frozen else ch.gamma_s
        # c* = sqrt(γ R T_c / M) / Γ(γ),  Γ(γ) = γ (2/(γ+1))^((γ+1)/(2(γ-1)))
        Gam = g * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
        c_ideal = math.sqrt(g * R * ch.T / ch.molar_mass) / Gam
        out["frozen" if frozen else "shifting"] = {"cstar": p.cstar, "cstar_ideal": c_ideal,
                                                   "gamma": g, "rel": p.cstar / c_ideal - 1}
    return out


def test_cstar_consistency():
    r = cstar_consistency()
    assert abs(r["frozen"]["rel"]) < 0.01, r
    assert abs(r["shifting"]["rel"]) < 0.01, r


def test_frozen_below_shifting():
    for fuel, of in [("LH2", 6.0), ("LCH4", 3.6), ("RP-1", 2.36)]:
        s = performance(fuel, "LOX", of, 10e6, 40)
        f = performance(fuel, "LOX", of, 10e6, 40, frozen=True)
        d = s.Isp_vac / f.Isp_vac - 1
        assert 0.005 < d < 0.10, (fuel, d)
        assert f.cstar < s.cstar


def test_published_ballpark():
    p = performance("LH2", "LOX", 6.0, 20.6e6, 69)
    assert p.Isp_vac > 450, p.Isp_vac                 # ≈ 450+ s theoretical
    Tmax = max(chamber("RP-1", "LOX", of, 6.895e6).T for of in np.arange(2.2, 3.4, 0.05))
    assert 3550 < Tmax < 3750, Tmax                   # ≈ 3600-3700 K peak at 1000 psia


def test_momentum_energy_consistency():
    """Exit velocity must equal sqrt(2(h_c - h_e)) and mass flux ratio must equal ε."""
    p = performance("LCH4", "LOX", 3.6, 30e6, 40)
    ue = math.sqrt(2 * (p.chamber._h0 - p.exit.h))
    assert abs(ue / p.u_e - 1) < 1e-12
    eps = (p.throat.rho * p.u_t) / (p.exit.rho * p.u_e)
    assert abs(eps / 40 - 1) < 1e-9
    assert abs(p.u_t / p.throat.sound_speed - 1) < 1e-9


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
