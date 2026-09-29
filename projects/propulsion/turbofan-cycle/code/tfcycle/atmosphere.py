"""US Standard Atmosphere, 1976 (lower 86 km).

Temperature is piecewise linear in *geopotential* altitude H; pressure follows the
hydrostatic equation layer by layer.  Aviation flight levels ("35,000 ft ISA") are
pressure altitudes, i.e. geopotential altitudes in this model, so ``atmosphere``
takes geopotential altitude by default.

Constants are those of the 1976 standard: g0 = 9.80665 m/s^2,
R* = 8314.32 J/(kmol K), M0 = 28.9644 kg/kmol, r0 = 6 356 766 m.
"""
import math

FT = 0.3048                       # m per ft
G0 = 9.80665
R_AIR = 8314.32 / 28.9644         # 287.0531... J/(kg K)
GAMMA = 1.4
R_EARTH = 6356766.0
# layer base geopotential altitude [m] and lapse rate [K/m]
_H = [0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0]
_L = [-0.0065, 0.0, 0.0010, 0.0028, 0.0, -0.0028, -0.0020]
_TB = [288.15]
_PB = [101325.0]
for _i in range(len(_L)):
    _dh = _H[_i + 1] - _H[_i]
    _t1 = _TB[_i] + _L[_i] * _dh
    if _L[_i] == 0.0:
        _p1 = _PB[_i] * math.exp(-G0 * _dh / (R_AIR * _TB[_i]))
    else:
        _p1 = _PB[_i] * (_t1 / _TB[_i]) ** (-G0 / (R_AIR * _L[_i]))
    _TB.append(_t1)
    _PB.append(_p1)


def geometric_to_geopotential(z):
    return R_EARTH * z / (R_EARTH + z)


def atmosphere(h, geometric=False):
    """Return (T [K], P [Pa], rho [kg/m^3], a [m/s]) at altitude ``h`` [m].

    ``h`` is geopotential (pressure) altitude unless ``geometric=True``.
    Valid for -1 km <= H <= 84.852 km.
    """
    H = geometric_to_geopotential(h) if geometric else float(h)
    if H < -1000.0 or H > _H[-1]:
        raise ValueError("altitude outside 1976 US Standard Atmosphere range")
    i = 0
    while i < len(_L) - 1 and H > _H[i + 1]:
        i += 1
    dh = H - _H[i]
    T = _TB[i] + _L[i] * dh
    if _L[i] == 0.0:
        P = _PB[i] * math.exp(-G0 * dh / (R_AIR * _TB[i]))
    else:
        P = _PB[i] * (T / _TB[i]) ** (-G0 / (R_AIR * _L[i]))
    rho = P / (R_AIR * T)
    a = math.sqrt(GAMMA * R_AIR * T)
    return T, P, rho, a
