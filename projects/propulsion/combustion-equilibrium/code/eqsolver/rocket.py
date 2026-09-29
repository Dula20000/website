"""Propellant definitions and ideal-rocket performance (infinite-area combustor).

Procedure (after RP-1311 ch. 6):

1. Chamber: HP equilibrium at P_c with the assigned enthalpy of the reactants.
2. Throat: isentropic expansion (s = s_c) to the pressure at which the flow
   velocity u = sqrt(2 (h_c - h)) equals the local sound speed.  Solved by
   Newton on ln(P_c/P) with the analytic slope d(u^2 - a^2)/d ln(P_c/P)
   ~= (γ+1) a^2 / γ.
3. Exit: the supersonic pressure giving the requested area ratio
   ε = (ρu)_t / (ρu)_e, solved by Newton on ln(P_c/P_e) using
   d ln ε / d ln(P_c/P) = (1 - 1/M^2)/γ, with bisection safeguarding.

"Shifting" means every expanded state is an SP equilibrium; "frozen" keeps the
chamber composition fixed (frozen at the combustor, as in CEA's default frozen
option).  Isp_vac = u_e + P_e ε / (ρu)_t;  Isp_SL subtracts P_a ε / (ρu)_t.
Flow separation at sea level is not modelled.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .thermo import R, SPECIES, ATOMIC_MASS
from .equilibrium import (GasSystem, equilibrate, frozen_state, default_system,
                          ConvergenceError)

G0 = 9.80665
P_ATM = 101325.0


@dataclass(frozen=True)
class Reactant:
    """A propellant with an assigned molar enthalpy.

    h is the absolute (formation-based) molar enthalpy in J/mol at the storage
    temperature T, consistent with the NASA-polynomial reference state.
    """
    name: str
    comp: dict
    h: float                 # J/mol
    T: float                 # K (informational)
    density: float = float("nan")   # kg/m^3 (liquid, at storage T)
    note: str = ""

    @property
    def molar_mass(self):     # kg/mol
        return sum(ATOMIC_MASS[e] * k for e, k in self.comp.items()) * 1e-3


# Representative values from public sources, all approximate.  Enthalpies are
# the RP-1311 / CEA thermo.inp "reactant" entries as I recall them (rounded):
#   O2(L)  at 90.17 K  : -12.979 kJ/mol
#   H2(L)  at 20.27 K  :  -9.012 kJ/mol
#   CH4(L) at 111.6 K  : ~ -89.0 kJ/mol   (gas Hf -74.6, minus vaporisation and
#                                           sensible enthalpy 298 -> 111 K)
#   RP-1  surrogate CH1.9423, Hf ~ -24.72 kJ per mol of CH1.9423 (CEA 'RP-1')
# Densities (normal boiling point / typical loading): LOX 1141, LH2 70.8,
#   LCH4 422.6, RP-1 ~ 810 kg/m^3.
PROPELLANTS = {
    "LOX":  Reactant("LOX", {"O": 2}, -12979.0, 90.17, 1141.0,
                     "O2(L) at NBP; h from CEA thermo.inp (approx.)"),
    "LH2":  Reactant("LH2", {"H": 2}, -9012.0, 20.27, 70.8,
                     "H2(L) at NBP; h from CEA thermo.inp (approx.)"),
    "LCH4": Reactant("LCH4", {"C": 1, "H": 4}, -89000.0, 111.6, 422.6,
                     "CH4(L) at NBP; h approx."),
    "RP-1": Reactant("RP-1", {"C": 1, "H": 1.9423}, -24717.0, 298.15, 810.0,
                     "surrogate CH1.9423 (CEA RP-1 entry, approx.)"),
    # gaseous reactants at 298.15 K for flame-temperature validation
    "H2g":  Reactant("H2g", {"H": 2}, SPECIES["H2"].h(298.15), 298.15),
    "O2g":  Reactant("O2g", {"O": 2}, SPECIES["O2"].h(298.15), 298.15),
    "CH4g": Reactant("CH4g", {"C": 1, "H": 4}, SPECIES["CH4"].h(298.15), 298.15),
    "N2g":  Reactant("N2g", {"N": 2}, SPECIES["N2"].h(298.15), 298.15),
}

PAIRS = {
    "LOX/LH2":  ("LH2", "LOX"),
    "LOX/LCH4": ("LCH4", "LOX"),
    "LOX/RP-1": ("RP-1", "LOX"),
}


def mix(reactants_mass):
    """Element vector b0 (mol/kg), element list and h0 (J/kg) for a list of
    (Reactant, mass_fraction) pairs."""
    elements = sorted({e for r, _ in reactants_mass for e in r.comp})
    b0 = np.zeros(len(elements))
    h0 = 0.0
    for r, w in reactants_mass:
        nm = w / r.molar_mass          # mol of reactant per kg of mixture
        h0 += nm * r.h
        for e, k in r.comp.items():
            b0[elements.index(e)] += nm * k
    return elements, b0, h0


def bipropellant(fuel, ox, of):
    f = PROPELLANTS[fuel] if isinstance(fuel, str) else fuel
    o = PROPELLANTS[ox] if isinstance(ox, str) else ox
    return mix([(f, 1.0 / (1.0 + of)), (o, of / (1.0 + of))])


def bulk_density(fuel, ox, of):
    f, o = PROPELLANTS[fuel], PROPELLANTS[ox]
    return (1.0 + of) / (of / o.density + 1.0 / f.density)


# --------------------------------------------------------------------------
def _expand(chamber, P, frozen, prev):
    if frozen:
        return frozen_state(chamber, P, s0=chamber._s0)
    return equilibrate(chamber.system, chamber._b0, P, s0=chamber._s0, init=prev or chamber)


def _throat(ch, frozen, tol=1e-11):
    g = ch.gamma_fr if frozen else ch.gamma_s
    x = g / (g - 1.0) * math.log((g + 1.0) / 2.0)     # ln(Pc/Pt) initial guess
    prev = None
    for _ in range(60):
        st = _expand(ch, ch.P * math.exp(-x), frozen, prev)
        prev = st
        u2 = 2.0 * (ch._h0 - st.h)
        a2 = st.sound_speed ** 2
        f = u2 - a2
        gm = st.gamma
        dx = -f * gm / ((gm + 1.0) * a2)
        x += dx
        if abs(f) / a2 < tol:
            break
    else:
        raise ConvergenceError("throat iteration failed")
    return st, math.sqrt(max(u2, 0.0))


def _exit(ch, eps, frozen, throat, mflux_t, tol=1e-11):
    xt = math.log(ch.P / throat.P)
    g = throat.gamma
    # ideal-gas estimate of supersonic Mach for the requested area ratio
    M = 2.0 + math.log(eps)
    for _ in range(50):
        t = 1 + 0.5 * (g - 1) * M * M
        e = (1 / M) * ((2 / (g + 1)) * t) ** ((g + 1) / (2 * (g - 1)))
        # Newton on ln eps with d ln e / dM = -1/M + M/t
        dM = (math.log(e) - math.log(eps)) / (-1 / M + M / t)
        M = max(M - dM, 1.0001)
        if abs(dM) < 1e-12:
            break
    x = math.log(ch.P / throat.P) + g / (g - 1) * math.log(t / (1 + 0.5 * (g - 1)))
    x = max(x, xt + 1e-3)
    lo, hi = xt, None
    prev = throat
    for _ in range(80):
        st = _expand(ch, ch.P * math.exp(-x), frozen, prev)
        prev = st
        u2 = 2.0 * (ch._h0 - st.h)
        u = math.sqrt(u2)
        e = mflux_t / (st.rho * u)
        f = math.log(e) - math.log(eps)
        if abs(f) < tol:
            break
        if f > 0:
            hi = x
        else:
            lo = x
        a2 = st.sound_speed ** 2
        slope = (1.0 - a2 / u2) / st.gamma
        xn = x - f / slope if slope > 1e-6 else float("nan")
        if not (lo < xn and (hi is None or xn < hi)) or not math.isfinite(xn):
            xn = 0.5 * (lo + hi) if hi is not None else x + 1.0
        x = xn
    else:
        raise ConvergenceError("exit area-ratio iteration failed")
    return st, u


@dataclass
class Performance:
    of: float
    Pc: float
    eps: float
    frozen: bool
    chamber: object
    throat: object
    exit: object
    cstar: float
    u_t: float
    u_e: float
    Isp_vac: float     # s
    Isp_sl: float      # s
    CF_vac: float
    CF_sl: float
    Pe: float

    def summary(self):
        ch = self.chamber
        return {k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in {
            "of": self.of, "Pc_MPa": self.Pc / 1e6, "eps": self.eps,
            "mode": "frozen" if self.frozen else "shifting",
            "Tc": ch.T, "M_c": ch.molar_mass * 1e3, "gamma_s_c": ch.gamma_s,
            "gamma_fr_c": ch.gamma_fr, "cp_eq_c": ch.cp_eq,
            "Tt": self.throat.T, "Pt_MPa": self.throat.P / 1e6,
            "Te": self.exit.T, "Pe_kPa": self.Pe / 1e3, "M_e": self.exit.molar_mass * 1e3,
            "gamma_e": self.exit.gamma,
            "cstar": self.cstar, "CF_vac": self.CF_vac, "CF_sl": self.CF_sl,
            "Isp_vac": self.Isp_vac, "Isp_sl": self.Isp_sl, "u_e": self.u_e,
        }.items()}


def chamber(fuel, ox, of, Pc):
    """HP (adiabatic, isobaric) combustion of a bipropellant at chamber pressure Pc."""
    elements, b0, h0 = bipropellant(fuel, ox, of)
    return chamber_from_mix(elements, b0, h0, Pc)


def chamber_from_mix(elements, b0, h0, Pc, system=None):
    sysm = system or default_system(elements)
    ch = equilibrate(sysm, b0, Pc, h0=h0)
    ch._b0, ch._h0, ch._s0 = b0, h0, ch.s
    return ch


def performance(fuel, ox, of, Pc, eps, *, frozen=False, Pa=P_ATM, ch=None):
    """Theoretical rocket performance for a bipropellant pair."""
    if ch is None:
        ch = chamber(fuel, ox, of, Pc)
    return nozzle(ch, eps, frozen=frozen, Pa=Pa, of=of)


def nozzle(ch, eps, *, frozen=False, Pa=P_ATM, of=float("nan")):
    th, ut = _throat(ch, frozen)
    mflux_t = th.rho * ut
    cstar = ch.P / mflux_t
    ex, ue = _exit(ch, eps, frozen, th, mflux_t)
    ivac = ue + ex.P * eps / mflux_t
    isl = ivac - Pa * eps / mflux_t
    return Performance(of, ch.P, eps, frozen, ch, th, ex, cstar, ut, ue,
                       ivac / G0, isl / G0, ivac / cstar, isl / cstar, ex.P)
