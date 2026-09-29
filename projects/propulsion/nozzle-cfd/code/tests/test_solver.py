"""Verification of the finite-volume solver against exact quasi-1D solutions."""
import _path  # noqa: F401
import numpy as np
from nozzlecfd import exact, geometry, solver


def test_quiescent_well_balanced():
    g = geometry.BellLikeNozzle(16.0)
    s = solver.NozzleSolver(g, 80, 1.2, pb=1.0, init="rest")
    for _ in range(50):
        s.step()
    rho, u, p = s._prim(s.Q)
    assert np.max(np.abs(u)) < 1e-13 and np.max(np.abs(p - 1)) < 1e-13


def _l1_mach_error(N, geom, g, pb):
    r = solver.NozzleSolver(geom, N, g, pb).run(max_iter=100000, tol=1e-12)
    ex = exact.exact_nozzle(r.x, geom, pb, g)
    return np.mean(np.abs(r.M - ex["M"])), r


def test_isentropic_second_order():
    geom = geometry.AndersonNozzle()
    e1, r1 = _l1_mach_error(100, geom, 1.4, 0.01)
    e2, r2 = _l1_mach_error(200, geom, 1.4, 0.01)
    assert r1.converged and r2.converged and r2.exit_supersonic
    order = np.log2(e1 / e2)
    assert e2 < 1e-3
    assert order > 1.7, order


def test_mass_conservation_and_choked_flow():
    geom = geometry.AndersonNozzle()
    r = solver.NozzleSolver(geom, 200, 1.4, 0.01).run(max_iter=100000, tol=1e-12)
    m = r.mdot_faces
    assert (m.max() - m.min()) / m.mean() < 1e-9
    assert abs(m.mean() / exact.mdot_star(1.4) - 1) < 2e-4


def test_shock_location():
    geom = geometry.AndersonShockNozzle()
    pb = exact.back_pressure_for_shock_at(float(geom.A(2.25)), geom.eps, 1.4)
    N = 200
    r = solver.NozzleSolver(geom, N, 1.4, pb).run(max_iter=100000, tol=1e-10)
    assert r.converged and not r.exit_supersonic
    xs = solver.shock_location(r, geom.x_throat)
    ex = exact.exact_nozzle(r.x, geom, pb, 1.4)
    assert abs(xs - ex["x_shock"]) < 0.5 * geom.L / N
    assert abs(r.p[-1] - pb) / pb < 2e-2


def test_exit_state_independent_of_ambient_when_supersonic():
    geom = geometry.BellLikeNozzle(16.0)
    c = exact.critical_back_pressures(16.0, 1.2)
    res = []
    for pb in (0.2 * c["p_sup"], c["p_sup"], 0.9 * c["p_nse"]):
        r = solver.NozzleSolver(geom, 120, 1.2, pb).run(max_iter=100000, tol=1e-11)
        assert r.exit_supersonic
        res.append(r.mom_faces[-1])
    assert max(res) - min(res) < 1e-9


def test_thrust_coefficient_vs_ideal():
    g, eps = 1.2, 16.0
    geom = geometry.BellLikeNozzle(eps)
    pa = 101325.0 / 9.7e6
    r = solver.NozzleSolver(geom, 400, g, pa).run(max_iter=100000, tol=1e-11)
    cf = r.mom_faces[-1] - pa * eps
    Me = exact.mach_from_area(eps, g, True)
    cf_i = exact.cf_ideal(float(exact.p_p0(Me, g)), pa, eps, g)
    assert abs(cf / cf_i - 1) < 2e-3


if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v(); print("ok", k)
