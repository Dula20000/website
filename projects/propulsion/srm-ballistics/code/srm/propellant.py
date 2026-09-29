"""Propellant burn-rate laws and ideal-nozzle relations.

Burn rate follows the Saint-Robert/Vieille law r = a (Pc / 1 MPa)^n, optionally
piecewise in pressure (as in Nakka's KNSB fits).  a is stored in m/s so that
r [m/s] = a * (P[Pa] / 1e6) ** n.
"""
from __future__ import annotations

import math

import numpy as np

P_REF = 1.0e6  # Pa, reference pressure of the a coefficients
G0 = 9.80665


class Propellant:
    """Propellant with (possibly piecewise) burn-rate law and ideal thermochemistry.

    laws: list of (p_min [Pa], p_max [Pa], a [m/s], n).  Below the first range
    the first law is extrapolated; above the last the last law is used.
    rho: delivered density [kg/m^3]; cstar [m/s]; gamma (ratio of specific heats).
    """

    def __init__(self, name, rho, laws, cstar, gamma):
        self.name, self.rho, self.laws = name, rho, [tuple(map(float, l)) for l in laws]
        self.cstar, self.gamma = cstar, gamma

    @property
    def Gam(self):
        g = self.gamma
        return math.sqrt(g) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))

    @property
    def RT(self):
        """R*T0 of the combustion gas from c* = sqrt(R T0) / Gamma."""
        return (self.cstar * self.Gam) ** 2

    def law_at(self, P):
        for lo, hi, a, n in self.laws:
            if P < hi:
                return a, n
        return self.laws[-1][2], self.laws[-1][3]

    def rate(self, P):
        if P <= 0:
            return 0.0
        a, n = self.law_at(P)
        return a * (P / P_REF) ** n

    def equilibrium_pc(self, Kn, rho_frac=1.0):
        """Lowest stable root of P = rho c* Kn r(P) (quasi-steady chamber pressure).

        For a single law: P = (rho a c* Kn / Pref^n)^(1/(1-n)).  For piecewise
        laws the lowest self-consistent root is returned; if the fit has a gap
        (discontinuity) with no root, the jump pressure is returned.
        """
        if Kn <= 0:
            return 0.0
        K = self.rho * rho_frac * self.cstar * Kn
        nl = len(self.laws)
        for i, (lo, hi, a, n) in enumerate(self.laws):
            P = (K * a / P_REF ** n) ** (1.0 / (1.0 - n))
            lo_eff = 0.0 if i == 0 else lo
            hi_eff = math.inf if i == nl - 1 else hi
            if lo_eff <= P < hi_eff:
                return P
            # discontinuity at hi: g(hi-) > 0 and g(hi+) < 0
            if i < nl - 1 and P >= hi_eff:
                a2, n2 = self.laws[i + 1][2], self.laws[i + 1][3]
                if K * a2 * (hi / P_REF) ** n2 < hi:
                    return hi
        return P

    def to_dict(self):
        return {"name": self.name, "rho": self.rho, "laws": [list(l) for l in self.laws],
                "cstar": self.cstar, "gamma": self.gamma}


def cstar_from_thermo(T0, M, gamma):
    """Ideal characteristic velocity from chamber temperature [K], molar mass [kg/kmol], gamma."""
    R = 8314.462618 / M
    Gam = math.sqrt(gamma) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
    return math.sqrt(R * T0) / Gam


# --------------------------------------------------------------------------
# nozzle
# --------------------------------------------------------------------------

def area_ratio(M, g):
    return (1.0 / M) * ((2 / (g + 1)) * (1 + 0.5 * (g - 1) * M * M)) ** ((g + 1) / (2 * (g - 1)))


def mach_from_area(eps, g, supersonic=True):
    """Invert the isentropic area-Mach relation by bisection."""
    lo, hi = (1.0, 100.0) if supersonic else (1e-9, 1.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f = area_ratio(mid, g) - eps
        if supersonic:
            lo, hi = (mid, hi) if f < 0 else (lo, mid)
        else:
            lo, hi = (lo, mid) if f < 0 else (mid, hi)
    return 0.5 * (lo + hi)


def p_ratio(M, g):
    return (1 + 0.5 * (g - 1) * M * M) ** (-g / (g - 1))


class Nozzle:
    """Ideal nozzle: throat diameter Dt [m], exit/throat area ratio eps.

    Ideal thrust coefficient (Sutton & Biblarz, ch. 3):
      CF = sqrt(2 g^2/(g-1) (2/(g+1))^((g+1)/(g-1)) (1 - (pe/pc)^((g-1)/g)))
           + (pe - pa)/pc * eps
    When pa/pc exceeds the subsonic-exit pressure ratio the throat unchokes and
    mass flow and thrust are computed from isentropic flow expanding to pa.
    """

    def __init__(self, eps, gamma):
        self.eps, self.g = eps, gamma
        self.pe_pc = p_ratio(mach_from_area(eps, gamma, True), gamma)
        self.psub = p_ratio(mach_from_area(eps, gamma, False), gamma)
        g = gamma
        self.cf_mom = math.sqrt(2 * g * g / (g - 1) * (2 / (g + 1)) ** ((g + 1) / (g - 1))
                                * (1 - self.pe_pc ** ((g - 1) / g)))

    def cf(self, pc, pa):
        return self.cf_mom + (self.pe_pc - pa / pc) * self.eps

    def flow(self, pc, pa, At, cstar, RT):
        """Returns (mdot, thrust)."""
        if pc <= pa:
            return 0.0, 0.0
        g = self.g
        if pa / pc <= self.psub:  # choked throat
            mdot = pc * At / cstar
            F = max(self.cf(pc, pa) * pc * At, 0.0)
            return mdot, F
        x = pa / pc
        Ae = At * self.eps
        mdot = Ae * pc * math.sqrt(2 * g / ((g - 1) * RT)) * math.sqrt(max(x ** (2 / g) - x ** ((g + 1) / g), 0.0))
        ue = math.sqrt(2 * g / (g - 1) * RT * max(1 - x ** ((g - 1) / g), 0.0))
        return mdot, mdot * ue
