"""Nozzle performance along an ascent: thrust coefficient, expansion state and separation-risk flags.

Separation criteria are empirical rules of thumb, not physics resolved by the
quasi-1D solver (which is inviscid and cannot separate):

* Summerfield criterion: separation likely when p_e/p_a < ~0.4 (often quoted as
  0.25-0.4; 0.4 is used here as the conservative end).
* Schmucker criterion: p_sep/p_a = (1.88 M_e - 1)^(-0.64), evaluated with the
  exit Mach number as the incipient-separation Mach number.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from . import exact
from .atmosphere import atmosphere, P_SL

SUMMERFIELD = 0.4


def schmucker_ratio(Me):
    """Schmucker separation pressure ratio p_sep/p_a for an exit Mach number Me."""
    return (1.88 * np.asarray(Me, float) - 1.0) ** (-0.64)


def cf_from_exit(mom_exit, pa_p0, eps):
    """Thrust coefficient from the non-dimensional exit momentum+pressure flux (rho u^2 + p) A_e / (p0 A*)."""
    return mom_exit - pa_p0 * eps


def optimal_altitude(pe_pa_fn, zmax=50e3):
    """Geometric altitude where p_e = p_a (None if outside 0..zmax)."""
    f = lambda z: np.log(pe_pa_fn(z))
    if f(0.0) >= 0:
        return 0.0 if abs(f(0.0)) < 1e-12 else None  # already under-expanded at sea level
    if f(zmax) < 0:
        return None
    return float(brentq(f, 0.0, zmax, xtol=0.1))


def altitude_for_ratio(pe, ratio, zmax=50e3):
    """Altitude where p_e/p_a equals ``ratio`` (e.g. 0.4); 0 if already above at sea level."""
    g = lambda z: pe / atmosphere(z)[1] - ratio
    if g(0.0) >= 0:
        return 0.0
    if g(zmax) < 0:
        return None
    return float(brentq(g, 0.0, zmax, xtol=0.1))


def eps_for_pe_p0(target_pe_p0, g):
    """Expansion ratio whose isentropic supersonic exit pressure is target (pe/p0)."""
    Me = np.sqrt(2.0 / (g - 1.0) * (target_pe_p0 ** (-(g - 1.0) / g) - 1.0))
    return float(exact.area_mach(Me, g)), float(Me)


def eps_max_summerfield(Pc, pa, g, k=SUMMERFIELD):
    """Largest expansion ratio with pe/pa >= k at chamber pressure Pc and ambient pa."""
    return eps_for_pe_p0(k * pa / Pc, g)[0]


def eps_max_schmucker(Pc, pa, g):
    """Largest expansion ratio with pe/pa >= (1.88 Me - 1)^-0.64 (Schmucker)."""
    def h(eps):
        Me = exact.mach_from_area(eps, g, True)
        return exact.p_p0(Me, g) * Pc / pa - schmucker_ratio(Me)
    lo, hi = 1.05, 2.0
    if h(lo) < 0:
        return None
    while h(hi) > 0:
        hi *= 2.0
        if hi > 1e5:
            return None
    return float(brentq(h, lo, hi, xtol=1e-8))


def eps_optimal(Pc, pa, g):
    """Expansion ratio for perfect expansion p_e = p_a."""
    return eps_for_pe_p0(pa / Pc, g)[0]


def pc_min_summerfield(eps, pa, g, k=SUMMERFIELD):
    """Lowest chamber pressure at which a nozzle of ratio eps keeps pe/pa >= k."""
    Me = exact.mach_from_area(eps, g, True)
    return float(k * pa / exact.p_p0(Me, g))


def ascent_sweep(pe_p0, Me, mom_exit, eps, Pc, g, z):
    """Performance vs geometric altitude z (m) for a fixed (supersonic) exit state."""
    _, pa, _ = atmosphere(z)
    pe = pe_p0 * Pc
    cf = cf_from_exit(mom_exit, pa / Pc, eps)
    cf_ideal = exact.cf_ideal(pe_p0, pa / Pc, eps, g)
    return dict(z_km=(np.asarray(z) / 1e3).tolist(), pa=pa.tolist(), pe_pa=(pe / pa).tolist(),
                cf=cf.tolist(), cf_ideal=np.asarray(cf_ideal).tolist(),
                summerfield_risk=(pe / pa < SUMMERFIELD).tolist(),
                schmucker_risk=(pe / pa < schmucker_ratio(Me)).tolist())
