"""Geometry verification: exact solutions, convergence order, area closure."""
import math

import numpy as np

import _runner  # noqa: F401  (sets sys.path)
from srm.geometry import (GrainTable, bates_burning_area, contour_segments, clipped_length, fmm, make_grid,
                          port_sdf, star_offset_perimeter)

R = 0.045


def test_circle_perimeter_second_order():
    errs, hs = [], []
    for N in (41, 81, 161):
        g = GrainTable({"type": "bates", "core_d": 0.030}, R, N=N, nw=30)
        m = g.w < g.w_contact - g.w[1]
        ex = 2 * math.pi * (0.015 + g.w)
        errs.append(np.max(np.abs(g.P_raw[m] - ex[m]) / ex[m])); hs.append(g.h)
    p = math.log(errs[1] / errs[2]) / math.log(hs[1] / hs[2])
    assert errs[-1] < 2e-4, errs
    assert 1.8 < p < 2.2, p


def test_fmm_reproduces_distance_for_circle():
    x, X, Y, h = make_grid(R, 121)
    phi = port_sdf({"type": "bates", "core_d": 0.030}, X, Y)
    T = fmm(phi, h)
    pos = phi > 0
    assert np.max(np.abs(T - phi)[pos]) < 0.05 * h


def test_bates_area_closed_form():
    g = GrainTable({"type": "bates", "core_d": 0.030}, R, N=161, nw=160)
    L = 0.13
    w = np.linspace(0, 0.0299, 50)
    Ab = g.perimeter(w) * (L - 2 * w) + 2 * (math.pi * R * R - g.port_area(w))
    Abx = bates_burning_area(2 * R, 0.030, L, w, ends=2)
    assert np.max(np.abs(Ab - Abx) / Abx) < 1e-3


def test_star_perimeter_converges():
    geo = {"type": "star", "points": 6, "r_tip": 0.030, "r_valley": 0.012}
    w = np.linspace(0.001, 0.0145, 15)
    Pex, wv = star_offset_perimeter(6, 0.030, 0.012, w)
    assert wv > 0.0145
    errs = []
    for N in (101, 201):
        x, X, Y, h = make_grid(R, N)
        T = fmm(port_sdf(geo, X, Y), h)
        P = np.array([clipped_length(*contour_segments(T, x, l), R) for l in w])
        errs.append(np.mean(np.abs(P - Pex) / Pex))
    assert errs[1] < 0.6 * errs[0] and errs[1] < 3e-3, errs


def test_area_closure_all_geometries():
    """Raw port area A0 + integral of P dw must fill the disc at burnout (co-area formula);
    first-order convergence for cornered ports, and exact closure after the correction."""
    geos = [{"type": "bates", "core_d": 0.030},
            {"type": "star", "points": 8, "r_tip": 0.030, "r_valley": 0.0077},
            {"type": "finocyl", "core_d": 0.016, "slots": 6, "slot_w": 0.0055, "slot_r": 0.023},
            {"type": "moon", "core_d": 0.030, "offset": 0.020}]
    for geo in geos:
        g1 = GrainTable(geo, R, N=161, nw=160)
        g2 = GrainTable(geo, R, N=321, nw=160)
        assert abs(g1.closure) < 0.012, (geo, g1.closure)
        assert abs(g2.closure) < 0.6 * abs(g1.closure) or abs(g2.closure) < 1e-3, (geo, g1.closure, g2.closure)
        assert abs(g1.A[-1] / (math.pi * R * R) - 1) < 1e-12


def test_initial_port_area_star_exact():
    n, rt, rv = 8, 0.030, 0.0077
    g = GrainTable({"type": "star", "points": n, "r_tip": rt, "r_valley": rv}, R, N=241, nw=10)
    exact = n * rt * rv * math.sin(math.pi / n)
    assert abs(g.A0 / exact - 1) < 2e-3


if __name__ == "__main__":
    _runner.run(globals())
