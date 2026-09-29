"""Exact quasi-one-dimensional nozzle relations (calorically perfect gas).

Everything here is analytic or a scalar root solve of an analytic relation:
isentropic area-Mach relation, normal-shock jump relations, and the classical
normal-shock-in-nozzle solution for a prescribed back pressure (Anderson,
*Modern Compressible Flow*, ch. 5). These are the reference solutions the
finite-volume solver is verified against.

Pressures are returned as ratios to the reservoir stagnation pressure p0,
temperatures as ratios to T0, densities as ratios to rho0 = p0/(R T0).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

__all__ = [
    "area_mach", "mach_from_area", "p_p0", "T_T0", "rho_rho0",
    "normal_shock", "p02_p01", "critical_back_pressures", "shock_area_ratio",
    "back_pressure_for_shock_at", "exact_nozzle", "mdot_star", "cf_ideal",
]


def area_mach(M, g):
    """A/A* for Mach M (isentropic)."""
    M = np.asarray(M, dtype=float)
    return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** ((g + 1.0) / (2.0 * (g - 1.0)))


def mach_from_area(ar, g, supersonic):
    """Invert the area-Mach relation for A/A* = ar >= 1 on the chosen branch."""
    if ar < 1.0:
        if ar > 1.0 - 1e-12:
            return 1.0
        raise ValueError("A/A* < 1")
    if abs(ar - 1.0) < 1e-14:
        return 1.0
    f = lambda M: area_mach(M, g) - ar
    if supersonic:
        hi = 2.0
        while f(hi) < 0:
            hi *= 2.0
        return brentq(f, 1.0, hi, xtol=1e-15, rtol=1e-15, maxiter=500)
    return brentq(f, 1e-12, 1.0, xtol=1e-15, rtol=1e-15, maxiter=500)


def T_T0(M, g):
    return 1.0 / (1.0 + 0.5 * (g - 1.0) * np.asarray(M, float) ** 2)


def p_p0(M, g):
    return T_T0(M, g) ** (g / (g - 1.0))


def rho_rho0(M, g):
    return T_T0(M, g) ** (1.0 / (g - 1.0))


def normal_shock(M1, g):
    """Return (M2, p2/p1, rho2/rho1, T2/T1) across a normal shock."""
    M1s = M1 * M1
    M2 = np.sqrt((1.0 + 0.5 * (g - 1.0) * M1s) / (g * M1s - 0.5 * (g - 1.0)))
    pr = 1.0 + 2.0 * g / (g + 1.0) * (M1s - 1.0)
    rr = (g + 1.0) * M1s / (2.0 + (g - 1.0) * M1s)
    return M2, pr, rr, pr / rr


def p02_p01(M1, g):
    """Stagnation-pressure ratio across a normal shock."""
    M1s = M1 * M1
    a = ((g + 1.0) * M1s / (2.0 + (g - 1.0) * M1s)) ** (g / (g - 1.0))
    b = (2.0 * g / (g + 1.0) * M1s - (g - 1.0) / (g + 1.0)) ** (-1.0 / (g - 1.0))
    return a * b


def critical_back_pressures(eps, g):
    """Characteristic back-pressure ratios pb/p0 for a nozzle of exit area ratio eps.

    Returns dict with
      p_sub  : isentropic subsonic exit (throat just choked, no shock),
      p_nse  : normal shock standing exactly at the exit plane,
      p_sup  : isentropic supersonic exit (design / perfectly expanded),
      Me_sup, Me_sub : the corresponding exit Mach numbers.
    """
    Msub = mach_from_area(eps, g, False)
    Msup = mach_from_area(eps, g, True)
    psub = float(p_p0(Msub, g))
    psup = float(p_p0(Msup, g))
    pnse = psup * normal_shock(Msup, g)[1]
    return dict(p_sub=psub, p_nse=float(pnse), p_sup=psup, Me_sup=Msup, Me_sub=Msub)


def shock_area_ratio(pb, eps, g):
    """Area ratio A_s/A* at which a normal shock stands for back pressure pb/p0.

    Valid for p_nse <= pb <= p_sub. Uses the closed-form exit Mach number from
    pe*Ae/(p01*A*) and then inverts the stagnation-pressure loss for M1.
    Returns (As/A*, M1, Me, p02/p01).
    """
    c = critical_back_pressures(eps, g)
    if not (c["p_nse"] - 1e-14 <= pb <= c["p_sub"] + 1e-14):
        raise ValueError("back pressure outside shock-in-nozzle range")
    k = (2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))
    r = 1.0 / (pb * eps)                           # p01 A* / (pe Ae)
    gm = g - 1.0
    Me2 = -1.0 / gm + np.sqrt(1.0 / gm ** 2 + (2.0 / gm) * k * r * r)
    Me = np.sqrt(Me2)
    p0e = pb * (1.0 + 0.5 * gm * Me2) ** (g / gm)   # = p02/p01
    if p0e >= 1.0 - 1e-13:
        return 1.0, 1.0, Me, 1.0
    M1 = brentq(lambda M: p02_p01(M, g) - p0e, 1.0 + 1e-12, c["Me_sup"] * (1 + 1e-9) + 1e-9,
                xtol=1e-15, rtol=1e-15, maxiter=500)
    return float(area_mach(M1, g)), float(M1), float(Me), float(p0e)


def back_pressure_for_shock_at(As, eps, g):
    """Inverse of shock_area_ratio: pb/p0 that places the shock at area ratio As (1<=As<=eps)."""
    M1 = mach_from_area(As, g, True)
    M2, pr, _, _ = normal_shock(M1, g)
    r0 = p02_p01(M1, g)
    Me = mach_from_area(eps / (1.0 / r0), g, False)  # downstream A2* = A*/r0
    return float(r0 * p_p0(Me, g))


def exact_nozzle(x, geom, pb, g):
    """Exact quasi-1D solution at points x for a nozzle geometry ``geom``.

    ``geom`` provides ``A(x)``, ``x_throat`` (location of the minimum area) and ``L``.
    pb is the back pressure ratio pb/p0. Returns dict with M, p, T, rho (ratios to
    the reservoir stagnation values), the regime ('subsonic', 'shock',
    'supersonic'), the shock position x_shock and area A_shock (if any), and the
    critical back pressures of the nozzle.
    """
    x = np.asarray(x, float)
    A = geom.A(x)
    Ath = float(geom.A(geom.x_throat))
    Ae = float(geom.A(geom.L))
    eps = Ae / Ath
    c = critical_back_pressures(eps, g)
    M = np.empty_like(x)
    p0loc = np.ones_like(x)
    xs = None
    As = None
    if pb >= c["p_sub"]:
        regime = "subsonic"
        Me = np.sqrt(2.0 / (g - 1.0) * (pb ** (-(g - 1.0) / g) - 1.0)) if pb < 1 else 0.0
        Astar = Ae / area_mach(Me, g) if Me > 0 else np.inf
        for i, a in enumerate(A):
            M[i] = mach_from_area(max(a / Astar, 1.0), g, False) if np.isfinite(Astar) else 0.0
    else:
        conv = x <= geom.x_throat
        for i in np.where(conv)[0]:
            M[i] = mach_from_area(max(A[i] / Ath, 1.0), g, False)
        div = ~conv
        if pb > c["p_nse"]:
            regime = "shock"
            As_ratio, M1, Me, r0 = shock_area_ratio(pb, eps, g)
            As = As_ratio * Ath
            # A is monotone increasing on the divergent side: root-find the shock station
            if As_ratio <= 1.0 + 1e-14:
                xs = float(geom.x_throat)
            else:
                xs = float(brentq(lambda s: float(geom.A(s)) - As, geom.x_throat, geom.L, xtol=1e-14))
            Astar2 = Ath / r0
            for i in np.where(div)[0]:
                if x[i] < xs:
                    M[i] = mach_from_area(max(A[i] / Ath, 1.0), g, True)
                else:
                    M[i] = mach_from_area(max(A[i] / Astar2, 1.0), g, False)
                    p0loc[i] = r0
        else:
            regime = "supersonic"
            for i in np.where(div)[0]:
                M[i] = mach_from_area(max(A[i] / Ath, 1.0), g, True)
    T = T_T0(M, g)
    p = p0loc * p_p0(M, g)
    rho = p / T
    return dict(M=M, p=p, T=T, rho=rho, regime=regime, x_shock=xs, A_shock=As, crit=c)


def mdot_star(g):
    """Choked mass flux per unit throat area, non-dimensional (p0 = rho0 = 1, R T0 = 1)."""
    return np.sqrt(g) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))


def cf_ideal(pe_p0, pa_p0, eps, g):
    """Ideal thrust coefficient (Sutton & Biblarz form)."""
    mom = np.sqrt(2.0 * g * g / (g - 1.0) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))
                  * (1.0 - pe_p0 ** ((g - 1.0) / g)))
    return mom + (pe_p0 - pa_p0) * eps
