"""Gibbs free-energy minimisation for ideal-gas mixtures (NASA RP-1311 formulation).

The unknowns are the log mole numbers ln n_j of every gas species (per kg of
reactants), the log of the total mole number ln n, and for HP/SP problems ln T.
Following Gordon & McBride (1994), the Newton corrections Δln n_j are
eliminated analytically through the stationarity condition

    Δln n_j = -μ_j/RT + Σ_i a_ij π_i + Δln n + (H_j/RT) Δln T,

leaving a small linear system in the element potentials π_i (one per element),
Δln n and (for HP/SP) Δln T.  The resulting "reduced" Newton matrix is of size
(n_elements + 1) or (n_elements + 2) regardless of the number of species.

Step control follows RP-1311 §3.3 (λ1 limits large corrections of the major
species / T / n, λ2 keeps trace species from jumping past ~1e-4 mole fraction
in one step).

Only gaseous species are considered; condensed carbon is outside scope because
no condensed phase is thermodynamically stable for the O/F ranges studied
here (all of them are well on the oxygen-rich side of the soot limit).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .thermo import R, P_REF, SPECIES, CHO_SPECIES


class ConvergenceError(RuntimeError):
    pass


class GasSystem:
    """A set of gaseous species and the element-composition matrix."""

    def __init__(self, species_names, elements=None):
        sp = [SPECIES[n] for n in species_names]
        if elements is None:
            elements = sorted({e for s in sp for e in s.comp})
        self.elements = list(elements)
        # drop species that contain elements not present
        sp = [s for s in sp if all(e in self.elements for e in s.comp)]
        self.species = sp
        self.names = [s.name for s in sp]
        self.A = np.array([[s.comp.get(e, 0) for s in sp] for e in self.elements], float)
        self.M = np.array([s.molar_mass for s in sp]) * 1e-3          # kg/mol
        self.low = np.array([s.low for s in sp])
        self.high = np.array([s.high for s in sp])
        self.tmid = np.array([s.t_mid for s in sp])

    # vectorised thermo --------------------------------------------------
    def thermo(self, T):
        """Return (cp/R, h/RT, s0/R) arrays at temperature T."""
        c = np.where((T >= self.tmid)[:, None], self.high, self.low)
        a1, a2, a3, a4, a5, a6, a7 = c.T
        cp = a1 + T * (a2 + T * (a3 + T * (a4 + T * a5)))
        h = a1 + T * (a2 / 2 + T * (a3 / 3 + T * (a4 / 4 + T * a5 / 5))) + a6 / T
        s = a1 * math.log(T) + T * (a2 + T * (a3 / 2 + T * (a4 / 3 + T * a5 / 4))) + a7
        return cp, h, s


@dataclass
class State:
    """Converged equilibrium (or frozen) state, per kg of mixture."""
    system: GasSystem
    T: float
    P: float
    nj: np.ndarray               # mol/kg
    n: float                     # total mol/kg (= Σ nj at convergence)
    iterations: int = 0
    history: list = field(default_factory=list)   # convergence residual per iteration
    frozen: bool = False
    # filled in by derivatives()
    cp_eq: float = float("nan")          # J/(kg K)
    cp_fr: float = float("nan")
    dlnV_dlnT: float = 1.0
    dlnV_dlnP: float = -1.0
    gamma_s: float = float("nan")        # isentropic exponent (RP-1311 γ_s)
    gamma_fr: float = float("nan")

    # ---- mixture properties (per kg) -----------------------------------
    @property
    def names(self):
        return self.system.names

    @property
    def x(self):
        return self.nj / self.nj.sum()

    @property
    def molar_mass(self):
        """Mixture molecular weight M = 1/n, kg/mol  (RP-1311 'MW')."""
        return 1.0 / self.nj.sum()

    def _thermo(self):
        return self.system.thermo(self.T)

    @property
    def h(self):
        """Specific enthalpy, J/kg."""
        _, hRT, _ = self._thermo()
        return float(np.dot(self.nj, hRT) * R * self.T)

    def sj_R(self):
        _, _, sR = self._thermo()
        nsum = self.nj.sum()
        with np.errstate(divide="ignore"):
            lx = np.log(self.nj / nsum)
        return sR - lx - math.log(self.P / P_REF)

    @property
    def s(self):
        """Specific entropy, J/(kg K)."""
        sj = self.sj_R()
        mask = self.nj > 0
        return float(np.dot(self.nj[mask], sj[mask]) * R)

    @property
    def g(self):
        """Specific Gibbs energy, J/kg."""
        return self.h - self.T * self.s

    @property
    def rho(self):
        return self.P / (self.nj.sum() * R * self.T)

    @property
    def gamma(self):
        return self.gamma_fr if self.frozen else self.gamma_s

    @property
    def sound_speed(self):
        return math.sqrt(self.nj.sum() * R * self.T * self.gamma)

    def element_residual(self, b0):
        return self.system.A @ self.nj - b0

    def mole_fractions(self, threshold=0.0):
        x = self.x
        return {n: float(v) for n, v in zip(self.names, x) if v > threshold}


# --------------------------------------------------------------------------
#  Newton iteration
# --------------------------------------------------------------------------
LN_TRACE = -18.420681   # ln(1e-8): RP-1311 trace-species threshold
LN_FLOOR = -80.0        # composition floor, ln(n_j/n)


def _initial(system, b0):
    ns = len(system.names)
    n = 0.1 * max(1.0, float(b0.sum()))
    return np.full(ns, n / ns), n


def equilibrate(system, b0, P, *, T=None, h0=None, s0=None, init=None,
                tol=1e-11, max_iter=400):
    """Solve for chemical equilibrium.

    Exactly one of ``T`` (TP problem), ``h0`` (HP, J/kg) or ``s0`` (SP, J/kg/K)
    fixes the energy condition.  ``b0`` is the element-abundance vector in
    mol/kg (order = system.elements).  ``init`` may be a previous State used as
    a warm start.
    """
    mode = "TP" if T is not None else ("HP" if h0 is not None else "SP")
    if mode != "TP" and init is None and T is None:
        T = 3800.0
    if init is not None:
        nj = init.nj.copy()
        n = float(init.nj.sum())
        if mode != "TP":
            T = init.T
    else:
        nj, n = _initial(system, b0)
        if T is None:
            T = 3800.0
    nj = np.maximum(nj, n * math.exp(LN_FLOOR))
    A = system.A
    nel = A.shape[0]
    lnP = math.log(P / P_REF)
    size = nel + 1 + (0 if mode == "TP" else 1)
    iT = nel + 1
    history = []

    for it in range(1, max_iter + 1):
        cpR, hRT, sR = system.thermo(T)
        lnx = np.log(nj / n)
        mu = hRT - sR + lnx + lnP                   # μ_j / RT
        G = np.zeros((size, size))
        rhs = np.zeros(size)
        An = A * nj                                 # a_ij n_j
        G[:nel, :nel] = An @ A.T
        G[:nel, nel] = An.sum(axis=1)
        G[nel, :nel] = An.sum(axis=1)
        G[nel, nel] = nj.sum() - n
        rhs[:nel] = b0 - A @ nj + An @ mu
        rhs[nel] = n - nj.sum() + nj @ mu
        if mode == "HP":
            v = An @ hRT
            G[:nel, iT] = v
            G[iT, :nel] = v
            G[nel, iT] = nj @ hRT
            G[iT, nel] = nj @ hRT
            G[iT, iT] = nj @ cpR + nj @ (hRT * hRT)
            rhs[iT] = h0 / (R * T) - nj @ hRT + nj @ (hRT * mu)
        elif mode == "SP":
            sj = sR - lnx - lnP
            G[:nel, iT] = An @ hRT
            G[nel, iT] = nj @ hRT
            G[iT, :nel] = An @ sj
            G[iT, nel] = nj @ sj
            G[iT, iT] = nj @ cpR + nj @ (hRT * sj)
            rhs[iT] = s0 / R - nj @ sj + n - nj.sum() + nj @ (sj * mu)
        sol = np.linalg.solve(G, rhs)
        pi = sol[:nel]
        dlnn = sol[nel]
        dlnT = sol[iT] if mode != "TP" else 0.0
        dlnnj = -mu + A.T @ pi + dlnn + hRT * dlnT

        # ---- RP-1311 step-size control ----
        big = lnx > LN_TRACE
        cand = [5 * abs(dlnT), 5 * abs(dlnn)]
        pos = big & (dlnnj > 0)
        if pos.any():
            cand.append(np.max(np.abs(dlnnj[pos])))
        m = max(cand)
        lam1 = 2.0 / m if m > 0 else 1.0
        lam2 = 1.0
        small = (~big) & (dlnnj >= 0)
        if small.any():
            denom = dlnnj[small] - dlnn
            ok = denom > 0
            if ok.any():
                lam2 = float(np.min(np.abs((-lnx[small][ok] - 9.2103404) / denom[ok])))
        lam = min(1.0, lam1, lam2)

        # convergence measure (RP-1311 criteria, tightened)
        nsum = nj.sum()
        err = max(float(np.max(nj * np.abs(dlnnj)) / nsum), n * abs(dlnn) / nsum, abs(dlnT))
        history.append(err)

        lnnj = np.log(nj) + lam * dlnnj
        n = math.exp(math.log(n) + lam * dlnn)
        lnnj = np.maximum(lnnj, math.log(n) + LN_FLOOR)
        nj = np.exp(lnnj)
        if mode != "TP":
            T = min(max(T * math.exp(lam * dlnT), 150.0), 7000.0)
        if err < tol and lam == 1.0:
            break
    else:
        raise ConvergenceError(f"{mode} equilibrium did not converge (err={err:.3e})")

    st = State(system, T, P, nj, float(nj.sum()), it, history)
    derivatives(st)
    return st


def derivatives(st: State):
    """Equilibrium and frozen thermodynamic derivatives (RP-1311 §2.5)."""
    sysm = st.system
    A = sysm.A
    nel = A.shape[0]
    nj = st.nj
    cpR, hRT, _ = sysm.thermo(st.T)
    An = A * nj
    size = nel + 1
    G = np.zeros((size, size))
    G[:nel, :nel] = An @ A.T
    G[:nel, nel] = An.sum(axis=1)
    G[nel, :nel] = An.sum(axis=1)
    # ∂/∂lnT at constant P
    rT = np.concatenate([-(An @ hRT), [-(nj @ hRT)]])
    xT = np.linalg.solve(G, rT)
    # ∂/∂lnP at constant T
    rP = np.concatenate([An.sum(axis=1), [nj.sum()]])
    xP = np.linalg.solve(G, rP)
    dlnn_dlnT = xT[nel]
    dlnn_dlnP = xP[nel]
    st.dlnV_dlnT = 1.0 + dlnn_dlnT
    st.dlnV_dlnP = -1.0 + dlnn_dlnP
    cp_eq_R = ((An @ hRT) @ xT[:nel] + (nj @ hRT) * dlnn_dlnT
               + nj @ cpR + nj @ (hRT * hRT))
    st.cp_fr = float(nj @ cpR) * R
    st.cp_eq = float(cp_eq_R) * R
    ntot = nj.sum()
    cv_eq = st.cp_eq + ntot * R * st.dlnV_dlnT ** 2 / st.dlnV_dlnP
    st.gamma_s = -(st.cp_eq / cv_eq) / st.dlnV_dlnP
    st.gamma_fr = st.cp_fr / (st.cp_fr - ntot * R)
    return st


def frozen_state(ref: State, P, *, s0=None, T=None, tol=1e-12):
    """State with the composition of ``ref`` frozen, at pressure P and either
    temperature T or entropy s0 (Newton on ln T)."""
    nj = ref.nj.copy()
    sysm = ref.system
    ntot = nj.sum()
    lnx = np.log(nj / ntot)
    lnP = math.log(P / P_REF)
    if T is None:
        T = ref.T
        for _ in range(100):
            cpR, _, sR = sysm.thermo(T)
            s = R * float(nj @ (sR - lnx - lnP))
            cp = R * float(nj @ cpR)
            d = (s0 - s) / cp
            T *= math.exp(d)
            if abs(d) < tol:
                break
        else:
            raise ConvergenceError("frozen isentrope did not converge")
    st = State(sysm, T, P, nj, float(ntot), frozen=True)
    cpR, _, _ = sysm.thermo(T)
    st.cp_fr = st.cp_eq = float(nj @ cpR) * R
    st.gamma_fr = st.cp_fr / (st.cp_fr - ntot * R)
    st.gamma_s = st.gamma_fr
    return st


def gibbs(st: State, nj=None):
    """Total Gibbs energy G/RT (per kg) of composition nj at the state's T, P."""
    nj = st.nj if nj is None else nj
    _, hRT, sR = st.system.thermo(st.T)
    n = nj.sum()
    return float(nj @ (hRT - sR + np.log(nj / n) + math.log(st.P / P_REF)))


def default_system(elements):
    names = list(CHO_SPECIES)
    if "N" in elements:
        names += ["N2", "NO"]
    if "Ar" in elements:
        names += ["Ar"]
    return GasSystem(names, elements)
