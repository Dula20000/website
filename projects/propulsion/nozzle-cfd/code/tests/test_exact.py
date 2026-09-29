"""Checks of the exact quasi-1D relations against textbook values and self-consistency."""
import _path  # noqa: F401
import numpy as np
from nozzlecfd import exact


def test_area_mach_roundtrip():
    for g in (1.2, 1.4):
        for M in (0.1, 0.5, 0.9, 1.5, 3.0, 5.0):
            ar = float(exact.area_mach(M, g))
            assert abs(exact.mach_from_area(ar, g, M > 1) - M) < 1e-10


def test_normal_shock_textbook_M2():
    # gamma = 1.4, M1 = 2: M2 = 0.5774, p2/p1 = 4.5, p02/p01 = 0.7209 (standard normal-shock tables)
    M2, pr, rr, tr = exact.normal_shock(2.0, 1.4)
    assert abs(M2 - 0.57735) < 1e-4
    assert abs(pr - 4.5) < 1e-12
    assert abs(rr - 2.6667) < 1e-4
    assert abs(exact.p02_p01(2.0, 1.4) - 0.7209) < 1e-4


def test_area_mach_textbook():
    # isentropic tables, gamma = 1.4: A/A* = 1.6875 at M = 2, p/p0 = 0.1278
    assert abs(float(exact.area_mach(2.0, 1.4)) - 1.6875) < 1e-4
    assert abs(float(exact.p_p0(2.0, 1.4)) - 0.12780) < 1e-5


def test_shock_position_inverse():
    for g, eps in ((1.4, 1.5), (1.4, 5.95), (1.2, 16.0), (1.2, 69.0)):
        c = exact.critical_back_pressures(eps, g)
        for f in (0.1, 0.4, 0.8):
            As = 1 + f * (eps - 1)
            pb = exact.back_pressure_for_shock_at(As, eps, g)
            assert c["p_nse"] <= pb <= c["p_sub"]
            As2 = exact.shock_area_ratio(pb, eps, g)[0]
            assert abs(As2 - As) / As < 1e-9


def test_limits_of_shock_range():
    g, eps = 1.4, 2.0
    c = exact.critical_back_pressures(eps, g)
    assert abs(exact.shock_area_ratio(c["p_nse"], eps, g)[0] - eps) < 1e-6
    assert abs(exact.shock_area_ratio(c["p_sub"], eps, g)[0] - 1.0) < 1e-6


def test_cf_ideal_vacuum_limit():
    # CF at pe = pa must equal the pure momentum term; large eps approaches CF_vac,max
    g = 1.2
    Me = exact.mach_from_area(16.0, g, True)
    pe = float(exact.p_p0(Me, g))
    cf = exact.cf_ideal(pe, pe, 16.0, g)
    mdot = exact.mdot_star(g)
    ue = Me * np.sqrt(g * float(exact.T_T0(Me, g)))
    assert abs(cf - mdot * ue) < 1e-10


if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v(); print("ok", k)
