"""Monte Carlo helper: vectorised QS matches the scalar solver; linearised 1/(1-n) sensitivities."""
import math

import numpy as np

import _runner  # noqa: F401
from srm import cases
from srm.ballistics import Motor
from srm.montecarlo import QSEnsemble, monte_carlo, oat_sensitivity


def _motor():
    return Motor(cases.apcp_bates(N=121, nw=120))


def test_vectorised_matches_scalar_qs():
    m = _motor()
    ens = QSEnsemble(m)
    nom = ens.nominal()
    qs = m.quasi_steady(500)
    assert abs(nom["Pc_max"] / qs["Pc_max"] - 1) < 2e-3
    assert abs(nom["I"] / qs["I"] - 1) < 5e-3


def test_linearised_sensitivities():
    m = _motor()
    sig = {"a": 0.03, "n": 0.02, "Dt": 0.005, "rho": 0.01}
    o = oat_sensitivity(m, sig)["params"]
    for p in ("a", "rho", "Dt", "n"):
        central = 0.5 * (math.log1p(o[p]["dPc_plus"]) - math.log1p(o[p]["dPc_minus"]))
        assert abs(central / o[p]["linear_dlnPc"] - 1) < 0.02, (p, central, o[p]["linear_dlnPc"])


def test_mc_statistics_and_anchor():
    m = _motor()
    sig = {"a": 0.03, "n": 0.02, "Dt": 0.005, "rho": 0.01}
    r = monte_carlo(m, sig, n_samples=4000, seed=1)
    n = m.prop.laws[0][3]
    pred = math.sqrt(sig["a"] ** 2 + sig["rho"] ** 2 + (2 * sig["Dt"]) ** 2
                     + (sig["n"] * math.log(r["nominal"]["Pc_max"] / 1e6)) ** 2) / (1 - n)
    cv = r["Pc_max"]["std"] / r["Pc_max"]["mean"]
    assert abs(cv / pred - 1) < 0.08, (cv, pred)
    ra = monte_carlo(m, sig, n_samples=4000, seed=1, anchor_P=r["nominal"]["Pc_max"])
    assert abs(ra["src_lnPc"]["n"]) < 0.1 * abs(r["src_lnPc"]["n"])


if __name__ == "__main__":
    _runner.run(globals())
