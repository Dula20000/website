"""Ballistics verification: equilibrium law, ODE -> QS, mass conservation, nozzle, class letters."""
import math

import numpy as np

import _runner  # noqa: F401
from srm import cases
from srm.ballistics import Motor, motor_class
from srm.propellant import Nozzle, P_REF, Propellant, cstar_from_thermo


def test_equilibrium_single_law_closed_form():
    p = Propellant("t", 1750.0, [(0, 1e12, 4e-3, 0.35)], 1500.0, 1.2)
    for Kn in (100, 250, 400):
        P = p.equilibrium_pc(Kn)
        exact = (1750 * 4e-3 * 1500 * Kn / P_REF ** 0.35) ** (1 / 0.65)
        assert abs(P / exact - 1) < 1e-12
        assert abs(P - 1750 * 1500 * Kn * p.rate(P)) / P < 1e-12


def test_equilibrium_piecewise_is_self_consistent():
    p = cases.knsb()
    for Kn in np.linspace(60, 600, 40):
        P = p.equilibrium_pc(Kn)
        g = p.rho * p.cstar * Kn * p.rate(P)
        # either a true root or the jump pressure of a fit discontinuity
        on_boundary = any(abs(P - l[1]) < 1e-6 for l in p.laws)
        assert on_boundary or abs(g / P - 1) < 1e-9, (Kn, P, g)


def test_knsb_cstar_derived():
    c = cstar_from_thermo(1600.0, 39.9, 1.1361)
    assert 850 < c < 950


def test_ode_matches_qs_mid_burn_and_mass_conserved():
    m = Motor(cases.apcp_bates(N=121, nw=120))
    qs = m.quasi_steady()
    tr = m.transient()
    wm = 0.5 * m.w_web
    i = int(np.searchsorted(tr["w"], wm))
    rel = np.interp(wm, qs["w"], qs["Pc"]) / tr["Pc"][i] - 1
    # QS neglects the gas-density term rho_g/rho_p (~1%) scaled by 1/(1-n); agreement to ~1 %
    assert abs(rel) < 0.015, rel
    mass = np.sum(0.5 * (tr["mdot"][1:] + tr["mdot"][:-1]) * np.diff(tr["t"]))
    assert abs(mass / m.m_prop - 1) < 2e-3


def test_ode_time_step_converged():
    s = cases.knsb_bates(N=101, nw=120)
    a = Motor(s).transient(dt_frac=0.2)
    b = Motor(s).transient(dt_frac=0.05)
    Ia = np.sum(0.5 * (a["F"][1:] + a["F"][:-1]) * np.diff(a["t"]))
    Ib = np.sum(0.5 * (b["F"][1:] + b["F"][:-1]) * np.diff(b["t"]))
    assert abs(Ia / Ib - 1) < 1e-3
    assert abs(a["Pc"].max() / b["Pc"].max() - 1) < 1e-3


def test_nozzle_limits():
    g = 1.2
    n = Nozzle(10.0, g)
    # vacuum CF below the infinite-expansion limit, and CF(pa=0) consistent with momentum + pressure terms
    cf_inf = math.sqrt(2 * g * g / (g - 1) * (2 / (g + 1)) ** ((g + 1) / (g - 1)))
    assert n.cf(5e6, 0.0) < cf_inf
    # mass flow continuous at the choking boundary
    pc = 101325.0 / n.psub
    RT = (1500 * math.sqrt(g) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))) ** 2
    m1, _ = n.flow(pc * (1 + 1e-7), 101325.0, 1e-4, 1500.0, RT)
    m2, _ = n.flow(pc * (1 - 1e-7), 101325.0, 1e-4, 1500.0, RT)
    assert abs(m1 / m2 - 1) < 1e-4


def test_motor_class_letters():
    assert motor_class(2.5) == "A" and motor_class(2.6) == "B"
    assert motor_class(640) == "I" and motor_class(641) == "J"
    assert motor_class(10240) == "M" and motor_class(11859) == "N"


def test_bates_rule_neutral_length():
    """d Ab/dw = 0 at mid-web when L = (3D + d)/2 (Ab(0) = Ab(burnout) for the closed form)."""
    from srm.geometry import bates_burning_area
    D, d = 0.09, 0.03
    L = (3 * D + d) / 2
    web = (D - d) / 2
    a0 = bates_burning_area(D, d, L, 0.0)
    a1 = bates_burning_area(D, d, L, web * (1 - 1e-12))
    assert abs(a0 / a1 - 1) < 1e-9


if __name__ == "__main__":
    _runner.run(globals())
