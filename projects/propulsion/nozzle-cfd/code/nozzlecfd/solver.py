"""Quasi-one-dimensional Euler finite-volume solver for convergent-divergent nozzles.

Governing equations (per unit area state Q = [rho, rho u, E]):

    d(Q A)/dt + d(F A)/dx = S,   F = [rho u, rho u^2 + p, (E + p) u],   S = [0, p dA/dx, 0]

Discretisation
--------------
* Cell-centred finite volumes on a uniform grid, cell volume V_i = A(x_i) dx,
  face areas A_{i+1/2} = A(x_{i+1/2}).
* Pressure source written as p_i (A_{i+1/2} - A_{i-1/2}); with this form a fluid at
  rest in any duct is an exact discrete steady state (well-balanced for u = 0).
* MUSCL reconstruction of primitive variables (rho, u, p) with a slope limiter
  (van Albada by default; minmod and van Leer available).
* HLLC approximate Riemann solver (Toro) with Einfeldt/Roe-average wave-speed
  bounds, which keeps the transonic rarefaction at the throat entropy-satisfying
  without an explicit fix.
* Three-stage, third-order SSP Runge-Kutta (Shu-Osher). Local (per-cell) time
  stepping may be used to accelerate convergence to the steady state.

Boundary conditions (two ghost cells each end)
----------------------------------------------
* Inflow: subsonic reservoir. Velocity is extrapolated from the interior; static
  T, p, rho follow from the prescribed p0, T0 isentropically.
* Outflow: if the (extrapolated) exit Mach number is >= 1 all primitives are
  extrapolated; otherwise rho and u are extrapolated and p is set to the back
  pressure. The switch is re-evaluated every stage.

Non-dimensionalisation: p0 = 1, rho0 = 1 (so R T0 = 1, reference speed sqrt(R T0)).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

NG = 2  # ghost cells per side


@dataclass
class Result:
    x: np.ndarray            # cell centres
    xf: np.ndarray           # faces
    A: np.ndarray            # area at cell centres
    Af: np.ndarray           # area at faces
    rho: np.ndarray
    u: np.ndarray
    p: np.ndarray
    M: np.ndarray
    mdot_faces: np.ndarray   # mass flow through each face (non-dim)
    mom_faces: np.ndarray    # momentum + pressure flux (rho u^2 + p) A through each face (non-dim)
    residual: list = field(default_factory=list)
    iters: int = 0
    converged: bool = False
    exit_supersonic: bool = True
    gamma: float = 1.4

    # face-extrapolated exit state (second-order extrapolation to x = L)
    def exit_state(self):
        def ex(q):
            return 1.5 * q[-1] - 0.5 * q[-2]
        rho, u, p = ex(self.rho), ex(self.u), ex(self.p)
        c = np.sqrt(self.gamma * p / rho)
        return dict(rho=rho, u=u, p=p, M=u / c, A=self.Af[-1], mdot=float(self.mdot_faces[-1]))


# ---------------------------------------------------------------- limiters
def _limit(dm, dp, kind):
    if kind == "minmod":
        return np.where(dm * dp > 0, np.sign(dm) * np.minimum(np.abs(dm), np.abs(dp)), 0.0)
    if kind == "vanleer":
        s = dm + dp
        return np.where(dm * dp > 0, 2 * dm * dp / np.where(s == 0, 1, s), 0.0)
    if kind == "vanalbada":
        e = 1e-12
        num = dm * (dp * dp + e) + dp * (dm * dm + e)
        return np.where(dm * dp > 0, num / (dm * dm + dp * dp + 2 * e), 0.0)
    if kind == "none":  # first order
        return np.zeros_like(dm)
    raise ValueError(kind)


# ---------------------------------------------------------------- HLLC flux
def hllc(rL, uL, pL, rR, uR, pR, g):
    """HLLC numerical flux for the 1-D Euler equations (vectorised)."""
    EL = pL / (g - 1) + 0.5 * rL * uL * uL
    ER = pR / (g - 1) + 0.5 * rR * uR * uR
    cL = np.sqrt(g * pL / rL)
    cR = np.sqrt(g * pR / rR)
    # Roe averages for Einfeldt bounds
    sL, sR = np.sqrt(rL), np.sqrt(rR)
    HL = (EL + pL) / rL
    HR = (ER + pR) / rR
    ut = (sL * uL + sR * uR) / (sL + sR)
    Ht = (sL * HL + sR * HR) / (sL + sR)
    ct = np.sqrt(np.maximum((g - 1) * (Ht - 0.5 * ut * ut), 1e-300))
    SL = np.minimum(np.minimum(uL - cL, uR - cR), ut - ct)
    SR = np.maximum(np.maximum(uR + cR, uL + cL), ut + ct)
    Ss = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR))

    FL = np.stack([rL * uL, rL * uL * uL + pL, (EL + pL) * uL])
    FR = np.stack([rR * uR, rR * uR * uR + pR, (ER + pR) * uR])
    UL = np.stack([rL, rL * uL, EL])
    UR = np.stack([rR, rR * uR, ER])

    def star(r, u, p, E, S):
        f = r * (S - u) / (S - Ss)
        return np.stack([f, f * Ss, f * (E / r + (Ss - u) * (Ss + p / (r * (S - u))))])

    UsL = star(rL, uL, pL, EL, SL)
    UsR = star(rR, uR, pR, ER, SR)
    F = np.where(SL >= 0, FL,
        np.where(Ss >= 0, FL + SL * (UsL - UL),
        np.where(SR > 0, FR + SR * (UsR - UR), FR)))
    return F


# ---------------------------------------------------------------- solver
class NozzleSolver:
    """Steady quasi-1D Euler solver.

    Parameters
    ----------
    geom : object with ``A(x)``, ``L``, ``x_throat``
    N : number of interior cells
    gamma : ratio of specific heats
    pb : back pressure ratio pb/p0
    limiter : 'vanalbada' | 'vanleer' | 'minmod' | 'none'
    cfl : Courant number
    local_dt : use local time stepping (steady-state acceleration)
    init : 'rest' (default) - reservoir gas at rest fills the nozzle and the back pressure is
           applied at t = 0, i.e. the physical start-up; 'ramp' - a generic Mach ramp
           (0.2 -> 1 -> 2). Starting from the supersonic ramp, a shock that belongs just
           inside the exit can lock onto the outflow boundary as a spurious discrete steady
           state; the start-up from rest reaches the correct position.
    """

    def __init__(self, geom, N, gamma, pb, limiter="vanalbada", cfl=0.8, local_dt=True, init="rest"):
        self.geom, self.N, self.g, self.pb = geom, int(N), float(gamma), float(pb)
        self.limiter, self.cfl, self.local_dt = limiter, float(cfl), bool(local_dt)
        L = geom.L
        self.dx = L / N
        self.xf = np.linspace(0.0, L, N + 1)
        self.x = 0.5 * (self.xf[:-1] + self.xf[1:])
        self.Af = geom.A(self.xf)
        self.A = geom.A(self.x)
        self.vol = self.A * self.dx
        self.dA = self.Af[1:] - self.Af[:-1]
        self.exit_supersonic = False
        self._init(init)

    # -- initial condition -------------------------------------------------
    def _init(self, kind):
        g, x, xt, L = self.g, self.x, self.geom.x_throat, self.geom.L
        if kind == "rest":
            # reservoir gas at rest filling the nozzle, back pressure applied at the exit
            rho = np.ones(self.N)
            u = np.zeros(self.N)
            p = np.ones(self.N)
        else:
            # generic Mach ramp (0.2 at inlet, 1 at throat, 2 at exit); not the exact solution
            M = np.where(x <= xt, 0.2 + 0.8 * x / xt, 1.0 + 1.0 * (x - xt) / (L - xt))
            T = 1.0 / (1.0 + 0.5 * (g - 1) * M * M)
            p = T ** (g / (g - 1))
            rho = p / T
            u = M * np.sqrt(g * T)
        E = p / (g - 1) + 0.5 * rho * u * u
        self.Q = np.stack([rho, rho * u, E])

    # -- primitives / BCs ----------------------------------------------------
    def _prim(self, Q):
        g = self.g
        rho = Q[0]
        u = Q[1] / rho
        p = (g - 1) * (Q[2] - 0.5 * rho * u * u)
        return rho, u, p

    def _with_ghosts(self, rho, u, p):
        """Fill two ghost cells per side.

        Extrapolation uses a minmod-limited slope of the last three cells, so it is
        linear in smooth flow but falls back towards zeroth order when a wave sits
        at the boundary (avoids spurious states when a shock passes the exit).
        """
        g, N = self.g, self.N
        R = np.empty(N + 2 * NG); U = np.empty_like(R); P = np.empty_like(R)
        R[NG:-NG], U[NG:-NG], P[NG:-NG] = rho, u, p

        def mm(a, b):
            return 0.0 if a * b <= 0 else (a if abs(a) < abs(b) else b)

        # inflow: velocity extrapolated (limited, kept subsonic), static state isentropic from reservoir
        su = mm(u[0] - u[1], u[1] - u[2])
        umax = 0.95 * np.sqrt(2.0 * g / (g + 1.0))          # ~ Mach 0.95 at stagnation enthalpy R T0 = 1
        for k in (1, 2):
            ug = min(max(u[0] + k * su, 0.0), umax)
            T = 1.0 - 0.5 * (g - 1) / g * ug * ug   # c_p T + u^2/2 = c_p T0 with R T0 = 1
            pg = T ** (g / (g - 1))
            R[NG - k], U[NG - k], P[NG - k] = pg / T, ug, pg
        # outflow
        sr = mm(rho[-1] - rho[-2], rho[-2] - rho[-3])
        sv = mm(u[-1] - u[-2], u[-2] - u[-3])
        sp = mm(p[-1] - p[-2], p[-2] - p[-3])
        # exit-face state (half-cell extrapolation to x = L) decides the regime
        rf = max(rho[-1] + 0.5 * sr, 1e-10)
        uf = u[-1] + 0.5 * sv
        pf = max(p[-1] + 0.5 * sp, 1e-10)
        Mf = uf / np.sqrt(g * pf / rf)
        prs = lambda M: 1.0 + 2.0 * g / (g + 1.0) * (M * M - 1.0)   # normal-shock pressure ratio
        up = None
        if Mf >= 1.0:
            # supersonic exit stays supersonic unless pb exceeds the pressure behind a normal shock at Mf
            self.exit_supersonic = self.pb <= pf * prs(Mf)
        else:
            self.exit_supersonic = False
            # A captured shock sitting in the last cells: if pb is below the pressure behind a normal
            # shock at the upstream (supersonic) Mach number, no steady shock can stand there and the
            # exit is treated as supersonic using that upstream state (lets a start-up shock leave).
            Mc = u / np.sqrt(g * p / rho)
            for j in (N - 2, N - 3, N - 4):
                if Mc[j] >= 1.0:
                    if self.pb <= p[j] * prs(Mc[j]):
                        self.exit_supersonic = True
                        up = j
                    break
        for k in (1, 2):
            j = -NG - 1 + k
            if up is not None:
                R[j], U[j], P[j] = rho[up], u[up], p[up]
            else:
                R[j] = max(rho[-1] + k * sr, 1e-10)
                U[j] = u[-1] + k * sv
                P[j] = max(p[-1] + k * sp, 1e-10) if self.exit_supersonic else self.pb
        return R, U, P

    # -- residual -------------------------------------------------------------
    def rhs(self, Q):
        """Return dQ/dt for interior cells (per unit volume) and the face fluxes F*A (3 x N+1)."""
        g = self.g
        rho, u, p = self._prim(Q)
        R, U, P = self._with_ghosts(rho, u, p)
        W = np.stack([R, U, P])
        d = np.diff(W, axis=1)                       # N+2NG-1 differences
        s = _limit(d[:, :-1], d[:, 1:], self.limiter)  # slopes for cells 1..N+2NG-2
        # cells with slopes: indices 1 .. N+2NG-2 ; faces between cell NG-1..NG+N-1 and next
        Wc = W[:, 1:-1]
        WL = Wc[:, :-1] + 0.5 * s[:, :-1]            # left state at face between c and c+1
        WR = Wc[:, 1:] - 0.5 * s[:, 1:]
        # keep faces NG-1+.. : interior faces are between cells NG-1..NG+N-1 (global)
        # Wc index j corresponds to global cell j+1; face j between global j+1 and j+2
        i0 = NG - 2
        WL = WL[:, i0:i0 + self.N + 1]
        WR = WR[:, i0:i0 + self.N + 1]
        # positivity fallback to first order where reconstruction goes bad
        bad = (WL[0] <= 0) | (WL[2] <= 0) | (WR[0] <= 0) | (WR[2] <= 0)
        if np.any(bad):
            WL[:, bad] = W[:, NG - 1:NG + self.N][:, bad]
            WR[:, bad] = W[:, NG:NG + self.N + 1][:, bad]
        F = hllc(WL[0], WL[1], WL[2], WR[0], WR[1], WR[2], g)
        FA = F * self.Af
        dQ = -(FA[:, 1:] - FA[:, :-1])
        dQ[1] += p * self.dA
        return dQ / self.vol, FA

    def _dt(self, Q):
        rho, u, p = self._prim(Q)
        c = np.sqrt(self.g * p / rho)
        dt = self.cfl * self.dx / (np.abs(u) + c)
        return dt if self.local_dt else np.full_like(dt, dt.min())

    def step(self):
        """One SSP-RK3 step. Returns the RMS density residual of the first stage."""
        Q0 = self.Q
        dt = self._dt(Q0)
        k1, _ = self.rhs(Q0)
        Q1 = Q0 + dt * k1
        k2, _ = self.rhs(Q1)
        Q2 = 0.75 * Q0 + 0.25 * (Q1 + dt * k2)
        k3, _ = self.rhs(Q2)
        self.Q = Q0 / 3.0 + 2.0 / 3.0 * (Q2 + dt * k3)
        return float(np.sqrt(np.mean(k1[0] ** 2)))

    def run(self, max_iter=200000, tol=1e-11, verbose=False, record_every=1):
        res_hist = []
        r0 = None
        it = 0
        conv = False
        for it in range(1, max_iter + 1):
            r = self.step()
            if r0 is None:
                r0 = r if r > 0 else 1.0
            if it % record_every == 0 or it == 1:
                res_hist.append((it, r / r0, r))
            if not np.isfinite(r):
                raise FloatingPointError("solver diverged")
            if r / r0 < tol:
                conv = True
                res_hist.append((it, r / r0, r))
                break
        return self.result(res_hist, it, conv)

    def result(self, res_hist=None, it=0, conv=False):
        rho, u, p = self._prim(self.Q)
        _, FA = self.rhs(self.Q)
        M = u / np.sqrt(self.g * p / rho)
        return Result(x=self.x, xf=self.xf, A=self.A, Af=self.Af, rho=rho, u=u, p=p, M=M,
                      mdot_faces=FA[0].copy(), mom_faces=FA[1].copy(), residual=res_hist or [], iters=it, converged=conv,
                      exit_supersonic=bool(self.exit_supersonic), gamma=self.g)


def shock_location(res, x_throat):
    """Shock position from a captured solution: sonic crossing (M from >1 to <1) downstream of the throat,
    linearly interpolated between cell centres. Returns None if no such crossing."""
    x, M = res.x, res.M
    idx = np.where((x[:-1] > x_throat) & (M[:-1] >= 1.0) & (M[1:] < 1.0))[0]
    if len(idx) == 0:
        return None
    i = idx[-1]
    return float(x[i] + (1.0 - M[i]) * (x[i + 1] - x[i]) / (M[i + 1] - M[i]))
