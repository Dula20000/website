"""Verification of the Gibbs minimiser.

* element conservation to machine precision,
* law of mass action for independent reactions (analytic K_p from the same
  thermo data) -- an exact check of the stationarity conditions,
* the converged composition is a true constrained minimum of G
  (random element-preserving perturbations always raise G, quadratically),
* independent cross-check against a general-purpose constrained optimiser
  (scipy SLSQP minimising G directly),
* HP energy balance and textbook adiabatic flame temperatures,
* quadratic Newton convergence.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from eqsolver.equilibrium import equilibrate, gibbs, default_system, GasSystem  # noqa: E402
from eqsolver.rocket import mix, PROPELLANTS, bipropellant, chamber  # noqa: E402
from eqsolver.thermo import SPECIES, P_REF, R  # noqa: E402


def h2o2_stoich():
    H2, O2 = PROPELLANTS["H2g"], PROPELLANTS["O2g"]
    m = 2 * H2.molar_mass + O2.molar_mass
    return mix([(H2, 2 * H2.molar_mass / m), (O2, O2.molar_mass / m)])


def ch4_air_stoich():
    C, O, N = PROPELLANTS["CH4g"], PROPELLANTS["O2g"], PROPELLANTS["N2g"]
    m = C.molar_mass + 2 * O.molar_mass + 2 * 3.76 * N.molar_mass
    return mix([(C, C.molar_mass / m), (O, 2 * O.molar_mass / m), (N, 7.52 * N.molar_mass / m)])


def flame_temperatures():
    out = {}
    el, b0, h0 = h2o2_stoich()
    out["H2_O2_stoich_1atm"] = equilibrate(default_system(el), b0, 101325.0, h0=h0).T
    el, b0, h0 = ch4_air_stoich()
    out["CH4_air_stoich_1atm"] = equilibrate(default_system(el), b0, 101325.0, h0=h0).T
    return out


def mass_action_residuals(st):
    """|ln K_p(thermo) - ln K_p(composition)| for independent reactions."""
    names = st.names
    _, hRT, sR = st.system.thermo(st.T)
    g = dict(zip(names, hRT - sR))            # g0/RT
    x = dict(zip(names, st.x))
    lp = math.log(st.P / P_REF)
    rx = [({"H2O": 1}, {"H2": 1, "O2": 0.5}),
          ({"H2O": 1}, {"OH": 1, "H": 1}),
          ({"H2": 1}, {"H": 2}),
          ({"O2": 1}, {"O": 2}),
          ({"CO2": 1}, {"CO": 1, "O": 1}),
          ({"HO2": 1}, {"H": 1, "O2": 1})]
    res = []
    for reac, prod in rx:
        if not all(k in g for k in list(reac) + list(prod)):
            continue
        dg = sum(v * g[k] for k, v in prod.items()) - sum(v * g[k] for k, v in reac.items())
        lnK = -dg
        dnu = sum(prod.values()) - sum(reac.values())
        lnQ = (sum(v * math.log(x[k]) for k, v in prod.items())
               - sum(v * math.log(x[k]) for k, v in reac.items()) + dnu * lp)
        res.append(abs(lnK - lnQ))
    return res


def test_element_conservation():
    for fuel, ox, of in [("LH2", "LOX", 6.0), ("LCH4", "LOX", 3.6), ("RP-1", "LOX", 2.36),
                         ("LH2", "LOX", 2.0), ("RP-1", "LOX", 6.0)]:
        ch = chamber(fuel, ox, of, 10e6)
        r = np.abs(ch.element_residual(ch._b0)) / ch._b0
        assert r.max() < 1e-13, (fuel, of, r)


def test_energy_balance_hp():
    ch = chamber("LCH4", "LOX", 3.6, 30e6)
    assert abs(ch.h - ch._h0) < 1e-9 * abs(ch._h0) + 1e-6


def test_mass_action():
    for fuel, ox, of, P in [("LH2", "LOX", 6.0, 20.6e6), ("RP-1", "LOX", 2.36, 9.7e6),
                            ("LCH4", "LOX", 3.6, 1e5)]:
        ch = chamber(fuel, ox, of, P)
        assert max(mass_action_residuals(ch)) < 1e-9


def perturbation_check(st, nsamples=200, delta=1e-4, seed=1):
    """Return (min ΔG/RT over random null-space perturbations, observed order)."""
    rng = np.random.default_rng(seed)
    A = st.system.A
    # null space of A (element-preserving directions)
    _, sv, vt = np.linalg.svd(A)
    rank = int((sv > 1e-12).sum())
    Nsp = vt[rank:].T
    g0 = gibbs(st)
    mins, orders = [], []
    for _ in range(nsamples):
        d = Nsp @ rng.standard_normal(Nsp.shape[1])
        d /= np.linalg.norm(d)
        # scale by sqrt(nj) so trace species are not driven negative
        w = np.sqrt(st.nj)
        step = d * w
        step /= np.max(np.abs(step) / st.nj)          # |Δn_j| <= n_j
        step = Nsp @ (Nsp.T @ step)                   # re-project
        dg = []
        for h in (delta, delta / 2):
            n1 = st.nj + h * step
            if (n1 <= 0).any():
                break
            dg.append(gibbs(st, n1) - g0)
        if len(dg) < 2:
            continue
        mins.append(dg[0])
        orders.append(math.log(dg[0] / dg[1], 2) if dg[1] > 0 else float("nan"))
    return float(min(mins)), float(np.nanmedian(orders))


def test_gibbs_true_minimum():
    ch = chamber("LCH4", "LOX", 3.6, 10e6)
    # evaluate at the converged T, P (a TP minimum)
    st = equilibrate(ch.system, ch._b0, ch.P, T=ch.T)
    mn, order = perturbation_check(st)
    assert mn > 0
    assert abs(order - 2.0) < 0.15          # ΔG ~ h^2: a smooth minimum


def scipy_crosscheck(T=3500.0, P=5e6):
    """Minimise G directly with SLSQP (independent of the RP-1311 algorithm)."""
    from scipy.optimize import minimize
    el, b0, _ = bipropellant("LCH4", "LOX", 3.2)
    sysm = default_system(el)
    ref = equilibrate(sysm, b0, P, T=T)
    _, hRT, sR = sysm.thermo(T)
    g0 = hRT - sR + math.log(P / P_REF)
    scale = ref.nj.sum()

    def G(y):
        n = np.exp(y) * scale
        return float(n @ (g0 + np.log(n / n.sum()))) / scale

    def dG(y):
        n = np.exp(y) * scale
        return (g0 + np.log(n / n.sum())) * n / scale

    cons = {"type": "eq", "fun": lambda y: (sysm.A @ (np.exp(y) * scale) - b0) / b0,
            "jac": lambda y: (sysm.A * (np.exp(y) * scale)) / b0[:, None]}
    y0 = np.log(np.full(len(sysm.names), 1.0 / len(sysm.names)))
    r = minimize(G, y0, jac=dG, constraints=[cons], method="SLSQP",
                 options={"ftol": 1e-15, "maxiter": 2000})
    n = np.exp(r.x) * scale
    x_sc = n / n.sum()
    major = ref.x > 1e-4
    rel = np.abs(x_sc[major] - ref.x[major]) / ref.x[major]
    return {"max_rel_diff_major": float(rel.max()),
            "G_rpl311": gibbs(ref), "G_slsqp": float(n @ (g0 + np.log(n / n.sum())))}


def test_scipy_crosscheck():
    r = scipy_crosscheck()
    assert r["max_rel_diff_major"] < 1e-4, r
    assert r["G_rpl311"] <= r["G_slsqp"] + 1e-8 * abs(r["G_slsqp"])


def test_flame_temperatures():
    t = flame_temperatures()
    assert abs(t["H2_O2_stoich_1atm"] - 3080) < 15, t     # textbook ≈ 3080 K
    assert abs(t["CH4_air_stoich_1atm"] - 2226) < 15, t   # textbook ≈ 2226 K


def test_quadratic_convergence():
    el, b0, h0 = bipropellant("LH2", "LOX", 6.0)
    st = equilibrate(default_system(el), b0, 20e6, h0=h0)
    h = st.history
    assert st.iterations < 40
    # final full Newton steps: e_{k+1} <~ C e_k^2
    tail = [e for e in h if e < 1e-2 and e > 1e-14]
    ratios = [math.log(tail[i + 1]) / math.log(tail[i]) for i in range(len(tail) - 1)]
    assert max(ratios[-2:]) > 1.7, (tail, ratios)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
