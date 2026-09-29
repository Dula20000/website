"""Lumped (0-D) internal ballistics of a segmented solid rocket motor.

Two solutions are provided for the same motor:

* quasi-steady (QS): at every web distance w the chamber sits at the
  equilibrium pressure Pc = (rho_p a c* Kn / Pref^n)^(1/(1-n)), Kn = Ab/At, and
  time follows from dt = dw / r(Pc);
* transient ODE: mass balance of the chamber gas with a choked (or unchoked)
  nozzle, integrated with fixed-per-step RK4 from ambient pressure:
      V dPc/dt = R T0 [ (rho_p - rho_g) Ab r - mdot_noz ],   dw/dt = r,
  which captures the filling transient and the blow-down tail-off.

Optional, off by default: throat erosion (constant recession rate of the
throat diameter) and Lenoir-Robillard erosive burning evaluated per axial
segment from the upstream mass flux.
"""
from __future__ import annotations

import copy
import json
import math

import numpy as np

from .geometry import GrainTable
from .propellant import G0, Nozzle, Propellant

_TABLE_CACHE: dict = {}


def grain_table(geom, R, N, nw):
    key = json.dumps([geom, R, N, nw], sort_keys=True)
    if key not in _TABLE_CACHE:
        _TABLE_CACHE[key] = GrainTable(geom, R, N=N, nw=nw)
    return _TABLE_CACHE[key]


def motor_class(I):
    """NAR/TRA-style impulse class letter (A: 1.26-2.5 N s, each letter doubles)."""
    if I <= 0:
        return "-"
    if I <= 0.625:
        return "1/4A"
    if I <= 1.25:
        return "1/2A"
    k = max(0, math.ceil(math.log2(I / 2.5) - 1e-12))
    return chr(ord("A") + k) if k < 26 else "beyond Z"


class Motor:
    """Motor built from a JSON-style spec (see cases.py for examples).

    spec keys: R (grain outer radius, m), segments [{geom, L, ends, count}],
    prop (Propellant.to_dict()), Dt (m), eps, pa (Pa), grid {N, nw},
    optional erosive {on, alpha, beta}, throat_erosion (m/s of diameter).
    Segments are ordered from the head end to the nozzle.
    """

    def __init__(self, spec: dict):
        self.spec = copy.deepcopy(spec)
        s = self.spec
        p = s["prop"]
        self.prop = Propellant(p["name"], p["rho"], p["laws"], p["cstar"], p["gamma"])
        self.nozzle = Nozzle(s["eps"], self.prop.gamma)
        self.R = s["R"]
        self.pa = s.get("pa", 101325.0)
        self.Dt = s["Dt"]
        N, nw = s["grid"]["N"], s["grid"]["nw"]
        self.segs = []
        for seg in s["segments"]:
            tab = grain_table(seg["geom"], self.R, N, nw)
            for _ in range(seg.get("count", 1)):
                L0, ends = seg["L"], seg.get("ends", 0)
                w_end = min(tab.w_max, L0 / ends if ends else math.inf)
                self.segs.append({"tab": tab, "L0": L0, "ends": ends, "w_end": w_end})
        self.K = len(self.segs)
        self.disc = math.pi * self.R ** 2
        self.V_case = self.disc * sum(sg["L0"] for sg in self.segs) + s.get("V_extra", 0.0)
        self.erosive = s.get("erosive", {"on": False})
        self.throat_erosion = s.get("throat_erosion", 0.0)
        self.m_prop = self.prop.rho * sum(self.vprop(sg, 0.0) for sg in self.segs)
        self.w_web = min(min(sg["tab"].w_contact, sg["L0"] / sg["ends"] if sg["ends"] else math.inf)
                         for sg in self.segs)
        self.dw_tab = min(sg["tab"].w[1] for sg in self.segs)

    # ---------------- geometry of one segment ----------------
    def seg_state(self, sg, w):
        """(Ab, Aport, perimeter) of one segment at web distance w."""
        Lb = sg["L0"] - sg["ends"] * w
        if w >= sg["w_end"] or Lb <= 0:
            return 0.0, self.disc, 0.0
        tab = sg["tab"]
        P = float(np.interp(w, tab.w, tab.P))
        A = float(np.interp(w, tab.w, tab.A))
        return P * Lb + sg["ends"] * (self.disc - A), A, P

    def vprop(self, sg, w):
        Lb = sg["L0"] - sg["ends"] * w
        if w >= sg["w_end"] or Lb <= 0:
            return 0.0
        A = float(np.interp(w, sg["tab"].w, sg["tab"].A))
        return max(self.disc - A, 0.0) * Lb

    def ab_total(self, w):
        return sum(self.seg_state(sg, w)[0] for sg in self.segs)

    def At(self, t=0.0):
        return 0.25 * math.pi * (self.Dt + self.throat_erosion * t) ** 2

    # ---------------- quasi-steady solution ----------------
    def quasi_steady(self, n_pts=400):
        """Equilibrium-pressure solution marched in web distance."""
        w_all = max(sg["w_end"] for sg in self.segs)
        w = np.linspace(0.0, w_all, n_pts)[:-1]  # last point is burnout (Ab = 0)
        At = self.At(0.0)
        prop, noz = self.prop, self.nozzle
        Ab = np.array([self.ab_total(wi) for wi in w])
        Kn = Ab / At
        Pc = np.array([prop.equilibrium_pc(k) for k in Kn])
        r = np.array([prop.rate(p) for p in Pc])
        F = np.array([noz.flow(p, self.pa, At, prop.cstar, prop.RT)[1] for p in Pc])
        inv = 1.0 / r
        t = np.concatenate([[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(w))])
        I = float(np.sum(0.5 * (F[1:] + F[:-1]) * np.diff(t)))
        web = w <= 0.97 * self.w_web  # web-burn window (excludes the burn-out transition)
        return {"w": w, "t": t, "Pc": Pc, "F": F, "Kn": Kn, "r": r, "I": I,
                "Pc_max": float(Pc.max()),
                "neutrality_pc": float(Pc[web].max() / Pc[web].min()),
                "neutrality_kn": float(Kn[web].max() / Kn[web].min()),
                "progressivity": float(Kn[web][-1] / Kn[web][0])}

    # ---------------- transient ODE ----------------
    def _rhs(self, t, y):
        prop, noz = self.prop, self.nozzle
        Pc = y[0]
        r0 = prop.rate(Pc)
        At = self.At(t)
        states = [self.seg_state(sg, y[1 + k]) for k, sg in enumerate(self.segs)]
        rates = [r0] * self.K
        if self.erosive.get("on"):
            alpha, beta = self.erosive["alpha"], self.erosive["beta"]
            up = 0.0
            for k, (Ab, A, P) in enumerate(states):
                mk = prop.rho * Ab * r0
                if Ab > 0 and P > 0 and A > 0:
                    G = (up + 0.5 * mk) / A
                    Dh = 4.0 * A / P
                    if G > 0:
                        rates[k] = r0 + alpha * G ** 0.8 * Dh ** -0.2 * math.exp(-beta * r0 * prop.rho / G)
                up += mk
        mgen = 0.0; dV = 0.0; vprop = 0.0
        for k, sg in enumerate(self.segs):
            Ab = states[k][0]
            mgen += prop.rho * Ab * rates[k]
            dV += Ab * rates[k]
            vprop += self.vprop(sg, y[1 + k])
        V = self.V_case - vprop
        mout, F = noz.flow(Pc, self.pa, At, prop.cstar, prop.RT)
        rho_g = Pc / prop.RT
        dP = prop.RT / V * (mgen - rho_g * dV - mout)
        dy = [dP] + [rates[k] if states[k][0] > 0 else 0.0 for k in range(self.K)]
        Ab_tot = sum(st[0] for st in states)
        return dy, {"F": F, "Kn": Ab_tot / At, "V": V, "mout": mout, "r": r0, "At": At}

    def transient(self, dt_frac=0.2, max_steps=400000):
        """RK4 integration from ambient pressure through burnout and blow-down."""
        prop = self.prop
        y = [self.pa] + [0.0] * self.K
        t = 0.0
        out = {"t": [], "Pc": [], "F": [], "Kn": [], "w": [], "mdot": []}
        Pmax = self.pa
        for _ in range(max_steps):
            dy, aux = self._rhs(t, y)
            out["t"].append(t); out["Pc"].append(y[0]); out["F"].append(aux["F"])
            out["Kn"].append(aux["Kn"]); out["w"].append(y[1]); out["mdot"].append(aux["mout"])
            Pmax = max(Pmax, y[0])
            burned = all(y[1 + k] >= sg["w_end"] for k, sg in enumerate(self.segs))
            if burned and y[0] - self.pa < 0.01 * (Pmax - self.pa):
                break
            tau = aux["V"] * prop.cstar / (prop.RT * aux["At"])
            dt = dt_frac * tau
            if not burned:
                dt = min(dt, 0.5 * self.dw_tab / max(aux["r"], 1e-6))
            k1 = dy
            k2, _ = self._rhs(t + 0.5 * dt, [a + 0.5 * dt * b for a, b in zip(y, k1)])
            k3, _ = self._rhs(t + 0.5 * dt, [a + 0.5 * dt * b for a, b in zip(y, k2)])
            k4, _ = self._rhs(t + dt, [a + dt * b for a, b in zip(y, k3)])
            y = [a + dt / 6.0 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(y, k1, k2, k3, k4)]
            y[0] = max(y[0], self.pa * 0.5)
            t += dt
        return {k: np.array(v) for k, v in out.items()}

    # ---------------- metrics ----------------
    def metrics(self, tr, qs=None):
        t, F, Pc = tr["t"], tr["F"], tr["Pc"]
        I = float(np.sum(0.5 * (F[1:] + F[:-1]) * np.diff(t)))
        Fmax = float(F.max())
        on = np.where(F >= 0.1 * Fmax)[0]
        t0, t1 = float(t[on[0]]), float(t[on[-1]])
        tb = t1 - t0
        Favg = I / tb
        cls = motor_class(I)
        # thrust-shape ratio: mean thrust over the last third of the action time / first third
        def mean_between(a, b):
            k = (t >= a) & (t <= b)
            return float(np.sum(0.5 * (F[k][1:] + F[k][:-1]) * np.diff(t[k])) / (t[k][-1] - t[k][0]))
        shape = mean_between(t0 + 2 * tb / 3, t1) / mean_between(t0, t0 + tb / 3)
        m = {"I_total": I, "F_max": Fmax, "F_avg": Favg, "t_burn": tb, "t_start": t0, "t_end": t1,
             "Pc_max": float(Pc.max()), "m_prop": self.m_prop, "Isp": I / (self.m_prop * G0),
             "class": cls, "designation": cls + str(int(round(Favg))) if len(cls) == 1 else cls,
             "Kn_initial": float(tr["Kn"][0]), "Kn_max": float(tr["Kn"].max()),
             "At": self.At(0.0), "w_web": self.w_web, "shape_ratio": shape}
        if qs is not None:
            m.update({"neutrality_pc": qs["neutrality_pc"], "neutrality_kn": qs["neutrality_kn"],
                      "progressivity": qs["progressivity"], "qs_I": qs["I"], "qs_Pc_max": qs["Pc_max"]})
        return m


def downsample(tr, n=500):
    """Keep ~n points uniform in time (plus start-up and the extremes) for JSON output."""
    t = np.asarray(tr["t"])
    L = len(t)
    grid = np.searchsorted(t, np.linspace(0.0, t[-1], n)).clip(0, L - 1)
    idx = np.unique(np.concatenate([grid, np.arange(min(L, 40)),
                                    [int(np.argmax(tr["Pc"])), int(np.argmax(tr["F"]))]]))
    return {k: [float(x) for x in np.asarray(v)[idx]] for k, v in tr.items()}
