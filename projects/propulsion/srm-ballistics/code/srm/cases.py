"""Case-study inputs.  Representative values from public sources, all approximate.

KNSB (KNO3/sorbitol 65/35): piecewise burn-rate fit, ideal density, chamber
temperature, molar mass and gamma as published by Richard Nakka on his amateur
rocketry web pages (recalled; treat as approximate).  c* is derived here from
T0, M and gamma rather than quoted.

APCP: AP/HTPB/Al composite with a, n, c*, gamma and density picked inside the
typical ranges given for composite propellants in Sutton & Biblarz.

PBAN (SRB-class): approximate widely quoted figures for the Space Shuttle SRB
propellant (burn rate ~0.37 in/s at ~625 psia, n ~0.35); used only for a
qualitative thrust-shape reproduction.
"""
from __future__ import annotations

import math

from .propellant import Propellant, cstar_from_thermo

MPA = 1.0e6
MM = 1.0e-3


def knsb():
    rho_ideal = 1841.0
    T0, M, gamma = 1600.0, 39.9, 1.1361
    laws = [  # (Pmin, Pmax, a [m/s @ MPa], n) -- Nakka's KNSB fit, r in mm/s with P in MPa
        (0.103 * MPA, 0.807 * MPA, 10.708 * MM, 0.625),
        (0.807 * MPA, 1.503 * MPA, 8.763 * MM, -0.314),
        (1.503 * MPA, 3.792 * MPA, 7.852 * MM, -0.013),
        (3.792 * MPA, 7.033 * MPA, 3.907 * MM, 0.535),
        (7.033 * MPA, 10.67 * MPA, 9.653 * MM, 0.064),
    ]
    return Propellant("KNSB (Nakka)", 0.95 * rho_ideal, laws, cstar_from_thermo(T0, M, gamma), gamma)


def apcp():
    return Propellant("APCP (AP/HTPB/Al, typical)", 1750.0,
                      [(0.0, 1e12, 4.0 * MM, 0.35)], 1500.0, 1.20)


def pban():
    # r = 0.368 in/s at 625 psia (approx.), n = 0.35  ->  a at 1 MPa
    r_ref = 0.368 * 0.0254
    p_ref = 625 * 6894.757
    a = r_ref / (p_ref / MPA) ** 0.35
    return Propellant("PBAN (SRB-class, approx.)", 1750.0, [(0.0, 1e12, a, 0.35)], 1530.0, 1.18)


PROPS = {"knsb": knsb, "apcp": apcp, "pban": pban}

PROP_SOURCES = {
    "knsb": "R. Nakka, KNSB propellant data (amateur rocketry web pages); c* derived from T0, M, gamma",
    "apcp": "Typical AP/HTPB composite ranges, Sutton & Biblarz, Rocket Propulsion Elements",
    "pban": "Widely quoted approximate Shuttle SRB propellant figures; qualitative use only",
}


def base(prop_key, R, segments, Dt, eps, N=161, nw=160, **kw):
    spec = {"R": R, "segments": segments, "prop": PROPS[prop_key]().to_dict(), "prop_key": prop_key,
            "Dt": Dt, "eps": eps, "pa": 101325.0, "grid": {"N": N, "nw": nw},
            "erosive": {"on": False, "alpha": 2.0e-6, "beta": 53.0}, "throat_erosion": 0.0}
    spec.update(kw)
    return spec


def knsb_bates(N=161, nw=160):
    """54 mm-class amateur motor: four 45 x 15 x 70 mm BATES grains, all ends burning."""
    return base("knsb", 0.0225, [{"geom": {"type": "bates", "core_d": 0.015}, "L": 0.070, "ends": 2, "count": 4}],
                Dt=0.0105, eps=7.0, N=N, nw=nw, name="KNSB 54 mm BATES (4 grains)")


def apcp_bates(N=161, nw=160):
    """98 mm-class high-power motor: four 90 x 30 x 130 mm BATES grains."""
    return base("apcp", 0.045, [{"geom": {"type": "bates", "core_d": 0.030}, "L": 0.130, "ends": 2, "count": 4}],
                Dt=0.022, eps=8.0, N=N, nw=nw, name="APCP 98 mm BATES (4 grains)")


def srb_class(N=161, nw=160):
    """Reduced-fidelity SRB-class booster: forward 11-point star segment plus a stepped-bore aft grain.

    The aft grain is represented by 16 cylindrical-bore slices whose bore grows
    linearly from 1.4 m to 2.0 m (a crude stand-in for a tapered bore), ends
    inhibited.  Qualitative shape reproduction only -- not a model of the real motor.
    """
    R = 1.80
    star = {"type": "star", "points": 11, "r_tip": 1.50, "r_valley": 0.60}
    segs = [{"geom": star, "L": 6.0, "ends": 0, "count": 1}]
    nb = 16
    for k in range(nb):
        d = round(1.4 + 0.6 * k / (nb - 1), 4)
        segs.append({"geom": {"type": "bates", "core_d": d}, "L": 1.8, "ends": 0, "count": 1})
    return base("pban", R, segs, Dt=1.37, eps=7.5, N=N, nw=nw, name="SRB-class star + stepped bore (qualitative)")


def equal_mass_geometries(R=0.045, core_d=0.030):
    """Four 2-D cross-sections with the same initial port area (hence same propellant mass)."""
    A = math.pi * (0.5 * core_d) ** 2
    n = 8; r_tip = 0.030
    r_valley = A / (n * r_tip * math.sin(math.pi / n))
    return {
        "bates": {"type": "bates", "core_d": core_d},
        "star": {"type": "star", "points": n, "r_tip": r_tip, "r_valley": r_valley},
        "finocyl": {"type": "finocyl", "core_d": 0.016, "slots": 6, "slot_w": 0.0055, "slot_r": None},
        "moon": {"type": "moon", "core_d": core_d, "offset": 0.020},
    }
