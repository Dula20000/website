"""Fan-diameter estimate and a deliberately simple installation penalty.

Fan diameter
    The fan face (station 2) passes the total air flow at an assumed axial Mach
    number M2 through an annulus of hub/tip ratio h:
        A2 = mdot0 sqrt(Tt2) / (Pt2 MFP(M2)),   D = sqrt(4 A2 / (pi (1 - h^2))).

Installation penalty (per engine, at the cruise point)
    * nacelle skin-friction drag  D_nac = q0 Cf FQ S_wet,
      S_wet = pi (k_D D)(k_L D)  (max nacelle diameter k_D*D, length k_L*D),
      Cf from the Schlichting turbulent flat-plate law at Re based on nacelle length,
      FQ = form factor x interference factor;
    * weight drag  D_W = m_eng g / (L/D)_aircraft, with propulsion-system mass scaled
      as m_eng = m_ref (D / D_ref)^2.
    Both scale with D^2, i.e. with mdot0, so for a fixed thrust the penalty is
        phi = (D_nac + D_W) / F_gross = c_pen / (F/mdot0),
    with c_pen [m/s] a constant for the flight condition.  The installed TSFC is
        S_inst = mdot_f / (F_gross - D_nac - D_W) = S / (1 - phi).
"""
import math
from .gas import Gas


def fan_diameter(mdot0, Tt2, Pt2, M2=0.60, hub_tip=0.30, gamma=1.4, cp=1004.0):
    g = Gas(gamma, cp)
    A = mdot0 * math.sqrt(Tt2) / (Pt2 * g.MFP(M2))
    return math.sqrt(4.0 * A / (math.pi * (1.0 - hub_tip ** 2)))


def schlichting_cf(Re):
    return 0.455 / math.log10(Re) ** 2.58


def mu_sutherland(T):
    return 1.458e-6 * T ** 1.5 / (T + 110.4)


def penalty_constant(r, M2=0.60, hub_tip=0.30, kD=1.20, kL=1.50, FQ=1.40,
                     m_ref=2400.0, D_ref=1.55, LoD=17.0, k_mult=1.0):
    """Installation penalty constant c_pen [m/s] for a design-point result ``r``
    (from ``ondesign.design``) and the breakdown of its two parts.

    Because D^2 is proportional to mdot0, c_pen does not depend on engine size:
    it is evaluated here for mdot0 = 1 kg/s, with the nacelle Reynolds number taken
    at a representative 1.6 m fan (weak, logarithmic dependence)."""
    T0, P0 = r["T0"], r["P0"]
    i = r["inputs"]
    gamma = i["gamma_c"]
    R = (gamma - 1) / gamma * i["cp_c"]
    M0 = i["M0"]
    q0 = 0.5 * gamma * P0 * M0 ** 2
    st = r["stations"]["2"]
    D1 = fan_diameter(1.0, st["Tt"], st["Pt"], M2, hub_tip, gamma, i["cp_c"])  # m per sqrt(kg/s)
    D2_per_mdot = D1 ** 2
    rho0 = P0 / (R * T0)
    V0 = r["V0"]
    Re = rho0 * V0 * (kL * 1.6) / mu_sutherland(T0)
    Cf = schlichting_cf(Re)
    Swet_per_D2 = math.pi * kD * kL
    c_nac = q0 * Cf * FQ * Swet_per_D2 * D2_per_mdot
    c_wt = m_ref / D_ref ** 2 * 9.80665 / LoD * D2_per_mdot
    return dict(c_pen=k_mult * (c_nac + c_wt), c_nac=k_mult * c_nac, c_wt=k_mult * c_wt,
                Cf=Cf, Re=Re, q0=q0, D2_per_mdot=D2_per_mdot)


def installed(r, pen):
    """Installed TSFC and penalty fraction for a design-point result and a penalty
    dict from ``penalty_constant``."""
    phi = pen["c_pen"] / r["Fs"]
    if phi >= 1:
        return dict(phi=phi, S_inst=float("inf"))
    return dict(phi=phi, S_inst=r["S"] / (1 - phi))
