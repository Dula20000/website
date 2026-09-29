"""Monte Carlo pressure dispersion with the vectorised quasi-steady solution.

For a single-law propellant the equilibrium pressure at web distance w is
    Pc(w) = (rho a c* Ab(w) / (At Pref^n))^(1/(1-n)),
so every multiplicative uncertainty (a, rho, At) is raised to 1/(1-n), and the
exponent itself enters as d ln Pc / dn = ln(Pc/Pref)/(1-n) (at fixed a).
The geometry Ab(w) is independent of the sampled parameters, so all samples
share one table and the whole ensemble is a few array operations.
"""
from __future__ import annotations

import math

import numpy as np

from .ballistics import Motor
from .propellant import G0, P_REF

PARAMS = ("a", "n", "Dt", "rho")


class QSEnsemble:
    def __init__(self, motor: Motor, n_w=500):
        if len(motor.prop.laws) != 1:
            raise ValueError("Monte Carlo helper assumes a single burn-rate law")
        self.m = motor
        w_all = max(sg["w_end"] for sg in motor.segs)
        self.w = np.linspace(0.0, w_all, n_w)[:-1]
        self.Ab = np.array([motor.ab_total(wi) for wi in self.w])
        self.V0 = motor.m_prop / motor.prop.rho
        _, a0, self.n0 = motor.prop.laws[0][1], motor.prop.laws[0][2], motor.prop.laws[0][3]
        self.a0 = a0
        noz = motor.nozzle
        self.cf_mom, self.pe_pc, self.eps = noz.cf_mom, noz.pe_pc, noz.eps

    def run(self, a, n, Dt, rho):
        """Vectorised QS solution for arrays of samples; returns dict of per-sample metrics."""
        a = np.atleast_1d(a)[:, None]; n = np.atleast_1d(n)[:, None]
        Dt = np.atleast_1d(Dt)[:, None]; rho = np.atleast_1d(rho)[:, None]
        m, cs = self.m, self.m.prop.cstar
        At = 0.25 * math.pi * Dt ** 2
        Kn = self.Ab[None, :] / At
        Pc = (rho * a * cs * Kn / P_REF ** n) ** (1.0 / (1.0 - n))
        r = a * (Pc / P_REF) ** n
        cf = self.cf_mom + (self.pe_pc - m.pa / Pc) * self.eps
        F = np.maximum(cf * Pc * At, 0.0)
        dw = np.diff(self.w)[None, :]
        I = np.sum(0.5 * (F[:, 1:] / r[:, 1:] + F[:, :-1] / r[:, :-1]) * dw, axis=1)
        tb = np.sum(0.5 * (1 / r[:, 1:] + 1 / r[:, :-1]) * dw, axis=1)
        mp = rho[:, 0] * self.V0
        return {"Pc_max": Pc.max(axis=1), "I": I, "t_b": tb, "Isp": I / (mp * G0)}

    def nominal(self):
        return {k: float(v[0]) for k, v in self.run(self.a0, self.n0, self.m.Dt, self.m.prop.rho).items()}


def monte_carlo(motor: Motor, sig: dict, n_samples=20000, seed=20260929, anchor_P=None):
    """Sample a (rel), n (abs), Dt (rel), rho (rel) as independent normals.

    If anchor_P is given, the burn-rate uncertainty is defined at that pressure:
    r(anchor_P) is perturbed by sig['a'] and n by sig['n'] about it, i.e. a and n
    become correlated through a = r_anchor / (anchor_P/Pref)^n.
    """
    ens = QSEnsemble(motor)
    rng = np.random.default_rng(seed)
    xi = rng.standard_normal((n_samples, 4))
    a0, n0, D0, r0 = ens.a0, ens.n0, motor.Dt, motor.prop.rho
    n = n0 + sig["n"] * xi[:, 1]
    if anchor_P is None:
        a = a0 * (1 + sig["a"] * xi[:, 0])
    else:
        r_anchor = a0 * (anchor_P / P_REF) ** n0 * (1 + sig["a"] * xi[:, 0])
        a = r_anchor / (anchor_P / P_REF) ** n
    Dt = D0 * (1 + sig["Dt"] * xi[:, 2])
    rho = r0 * (1 + sig["rho"] * xi[:, 3])
    res = ens.run(a, n, Dt, rho)
    nom = ens.nominal()
    P = res["Pc_max"]; I = res["I"]
    q = 99.865  # one-sided 3-sigma-equivalent percentile
    boot = np.random.default_rng(seed + 1)
    bs = [np.percentile(P[boot.integers(0, n_samples, n_samples)], q) for _ in range(400)]
    # standardised regression coefficients on ln(Pc_max) -> variance shares
    Xs = (xi - xi.mean(0)) / xi.std(0)
    y = np.log(P)
    coef, *_ = np.linalg.lstsq(np.column_stack([np.ones(n_samples), Xs]), y, rcond=None)
    src = coef[1:] * 1.0
    yI = np.log(I)
    coefI, *_ = np.linalg.lstsq(np.column_stack([np.ones(n_samples), Xs]), yI, rcond=None)
    return {
        "n_samples": n_samples, "sig": sig, "anchor_P": anchor_P, "nominal": nom,
        "Pc_max": {"mean": float(P.mean()), "std": float(P.std()), "p99865": float(np.percentile(P, q)),
                   "p99865_ci": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                   "mu3s": float(P.mean() + 3 * P.std()), "p50": float(np.median(P)), "min": float(P.min()),
                   "max": float(P.max())},
        "I": {"mean": float(I.mean()), "std": float(I.std()), "p00135": float(np.percentile(I, 100 - q)),
              "p99865": float(np.percentile(I, q))},
        "src_lnPc": dict(zip(PARAMS, map(float, src))),
        "src_lnI": dict(zip(PARAMS, map(float, coefI[1:]))),
        "r2_lnPc": float(1 - np.var(y - np.column_stack([np.ones(n_samples), Xs]) @ coef) / np.var(y)),
        "samples": {"Pc_max": P, "I": I},
    }


def oat_sensitivity(motor: Motor, sig: dict):
    """One-at-a-time +/-1 sigma effect on peak Pc and total impulse, with the linearised prediction."""
    ens = QSEnsemble(motor)
    a0, n0, D0, r0 = ens.a0, ens.n0, motor.Dt, motor.prop.rho
    nom = ens.nominal()
    base = dict(a=a0, n=n0, Dt=D0, rho=r0)
    out = {}
    Pn = nom["Pc_max"]
    lin = {"a": sig["a"] / (1 - n0), "rho": sig["rho"] / (1 - n0), "Dt": 2 * sig["Dt"] / (1 - n0),
           "n": sig["n"] * math.log(Pn / P_REF) / (1 - n0)}
    for p in PARAMS:
        vals = {}
        for s in (-1, 1):
            kw = dict(base)
            kw[p] = base[p] + s * sig[p] if p == "n" else base[p] * (1 + s * sig[p])
            r = ens.run(**kw)
            vals[s] = {k: float(v[0]) for k, v in r.items()}
        out[p] = {"dPc_plus": vals[1]["Pc_max"] / Pn - 1, "dPc_minus": vals[-1]["Pc_max"] / Pn - 1,
                  "dI_plus": vals[1]["I"] / nom["I"] - 1, "dI_minus": vals[-1]["I"] / nom["I"] - 1,
                  "linear_dlnPc": lin[p] * (-1 if p == "Dt" else 1)}
    return {"nominal": nom, "params": out}
