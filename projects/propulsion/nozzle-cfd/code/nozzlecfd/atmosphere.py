"""U.S. Standard Atmosphere, 1976 — lowest five layers (0 to 51 km geopotential).

Layer bases (geopotential km), base temperature (K) and lapse rate (K/km) are the
standard's defining values; pressure follows from hydrostatics with
g0 = 9.80665 m/s^2, M0 = 28.9644 kg/kmol, R* = 8314.32 J/(kmol K).
"""
from __future__ import annotations

import numpy as np

R_EARTH = 6356.766e3          # m, effective earth radius used by the standard
G0 = 9.80665
R_AIR = 8314.32 / 28.9644     # J/(kg K)
P_SL = 101325.0
T_SL = 288.15

# (base geopotential altitude m, lapse K/m)
_LAYERS = [
    (0.0, -6.5e-3),
    (11000.0, 0.0),
    (20000.0, 1.0e-3),
    (32000.0, 2.8e-3),
    (47000.0, 0.0),
    (51000.0, None),
]


def _bases():
    Tb, pb = [T_SL], [P_SL]
    for i in range(len(_LAYERS) - 2):
        h0, L = _LAYERS[i]
        h1 = _LAYERS[i + 1][0]
        T0, p0 = Tb[-1], pb[-1]
        T1 = T0 + L * (h1 - h0)
        if L == 0.0:
            p1 = p0 * np.exp(-G0 * (h1 - h0) / (R_AIR * T0))
        else:
            p1 = p0 * (T1 / T0) ** (-G0 / (R_AIR * L))
        Tb.append(T1)
        pb.append(p1)
    return Tb, pb


_TB, _PB = _bases()


def geopotential(z):
    """Geometric altitude z (m) to geopotential altitude h (m)."""
    return R_EARTH * np.asarray(z, float) / (R_EARTH + np.asarray(z, float))


def atmosphere(z):
    """Return (T [K], p [Pa], rho [kg/m^3]) at geometric altitude z (m), 0 <= z <~ 51.4 km."""
    h = np.atleast_1d(geopotential(z))
    if np.any(h > 51000.0 + 1e-6) or np.any(h < -1e-6):
        raise ValueError("altitude outside 0-51 km geopotential")
    T = np.empty_like(h)
    p = np.empty_like(h)
    for k, hk in enumerate(h):
        i = 0
        while i < 4 and hk >= _LAYERS[i + 1][0]:
            i += 1
        h0, L = _LAYERS[i]
        if L == 0.0:
            T[k] = _TB[i]
            p[k] = _PB[i] * np.exp(-G0 * (hk - h0) / (R_AIR * _TB[i]))
        else:
            T[k] = _TB[i] + L * (hk - h0)
            p[k] = _PB[i] * (T[k] / _TB[i]) ** (-G0 / (R_AIR * L))
    rho = p / (R_AIR * T)
    if np.ndim(z) == 0:
        return float(T[0]), float(p[0]), float(rho[0])
    return T, p, rho
