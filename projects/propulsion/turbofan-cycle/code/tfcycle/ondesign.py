"""Parametric (on-design) cycle analysis of a two-spool, separate-exhaust turbofan.

Formulation follows Mattingly, *Elements of Propulsion* (parametric "real cycle"
analysis) with two selectable gas models (see ``thermo.py``):

* ``thermo='varcp'`` (default): thermally perfect air and combustion products with
  cp(T, f) from NASA 7-coefficient polynomials.  The fuel-air ratio comes from an
  enthalpy balance on a common 298.15 K reference,
  ``h_air(Tt3) + f_b eta_b hPR = (1 + f_b) h_prod(Tt4, f_b)``, and the entropy
  function s0(T) is used through compressors, turbines and nozzles;
* ``thermo='cpg'``: Mattingly's two-gas "modified" cycle, (gamma_c, cp_c) cold and
  (gamma_t, cp_t) hot, which also contains the ideal cycle as a limit.

Components: polytropic efficiencies for fan, LPC, HPC, HPT, LPT (for any gas model
s0(T2) - s0(T1) = R ln(pi)/e for compression and e R ln(pi) for expansion);
total-pressure ratios for inlet (MIL-E-5007 ram recovery above M = 1), burner and
nozzles; burner efficiency; one mechanical efficiency per spool; convergent
(choked/unchoked) or ideally expanded nozzles.

Turbine cooling (simple bleed model): a fraction ``eps_cool`` of the core air is
taken off after the HPC (station 3), bypasses the burner and HPT, and is remixed at
constant total pressure after the HPT (station 4.4 -> 4.5, ahead of the LPT).  The
coolant does no HPT work, so the model is on the conservative side for a given Tt4.

Station numbering (Mattingly / ARP 755): 0 free stream, 2 fan face, 13 fan exit,
18/19 fan-nozzle throat/exit, 2.5 LPC (fan root + booster) exit, 3 HPC exit,
4 burner exit, 4.4 HPT exit, 4.5 LPT inlet after coolant mixing, 5 LPT exit,
8/9 core-nozzle throat/exit.

Spool arrangement: the LPT drives the fan (bypass, pi_f) and LPC (core, pi_cL);
the HPT drives the HPC (pi_cH).  OPR = pi_cL * pi_cH.  Customer bleed and shaft
power extraction are not modelled.
"""
from dataclasses import dataclass, asdict, replace
import math

from .atmosphere import atmosphere
from .thermo import (Thermo, nozzle, v_fully_expanded, compress_polytropic, turbine_pi_polytropic,
                     compressor_eta, turbine_eta, pr_isentropic)

TSFC_LB = 3600.0 * 9.80665          # kg/(N s)  ->  lbm/(lbf h)


@dataclass
class DesignInputs:
    # flight condition
    M0: float = 0.78
    alt: float = 10668.0            # geopotential (pressure) altitude [m]; 35 000 ft
    T0: float = None                # optional overrides of the atmosphere
    P0: float = None
    # cycle
    alpha: float = 5.3              # bypass ratio
    pi_f: float = 1.65              # fan (bypass stream) total-pressure ratio
    pi_cL: float = 2.7              # core-stream LPC ratio (fan root + booster)
    pi_cH: float = 11.0             # HPC ratio
    Tt4: float = 1500.0             # burner exit total temperature [K]
    eps_cool: float = 0.12          # turbine cooling air, fraction of core air (assumed)
    # gas model and fuel
    thermo: str = "varcp"           # 'varcp' (NASA polynomials) or 'cpg' (two-gas / ideal)
    gamma_c: float = 1.40           # used only when thermo == 'cpg'
    cp_c: float = 1004.0
    gamma_t: float = 1.33
    cp_t: float = 1156.0
    hPR: float = 42.8e6             # Jet-A lower heating value at 298.15 K [J/kg]
    # component figures of merit
    pi_d_max: float = 0.995
    pi_b: float = 0.96
    eta_b: float = 0.995
    e_f: float = 0.89
    e_cL: float = 0.89
    e_cH: float = 0.90
    e_tH: float = 0.89
    e_tL: float = 0.90
    eta_mH: float = 0.99
    eta_mL: float = 0.99
    pi_n: float = 0.99
    pi_fn: float = 0.99
    core_nozzle: str = "convergent"  # 'convergent' or 'expanded'
    fan_nozzle: str = "convergent"
    neglect_fuel_mass: bool = False  # Mattingly ideal-cycle assumption (1+f -> 1)
    mil_spec_recovery: bool = True   # apply MIL-E-5007 ram recovery above M0 = 1

    @property
    def OPR(self):
        return self.pi_cL * self.pi_cH

    def with_(self, **kw):
        return replace(self, **kw)

    def to_dict(self):
        return asdict(self)


def ram_recovery(M0):
    """MIL-E-5007 inlet ram recovery eta_r (=1 subsonic)."""
    if M0 <= 1.0:
        return 1.0
    if M0 < 5.0:
        return 1.0 - 0.075 * (M0 - 1.0) ** 1.35
    return 800.0 / (M0 ** 4 + 935.0)


def make_thermo(i):
    return Thermo(i.thermo, i.gamma_c, i.cp_c, i.gamma_t, i.cp_t)


def freestream(i, air):
    """Static and total free-stream state."""
    T0a, P0a, _, _ = atmosphere(i.alt)
    T0 = i.T0 if i.T0 is not None else T0a
    P0 = i.P0 if i.P0 is not None else P0a
    a0 = math.sqrt(air.gamma(T0) * air.R * T0)
    V0 = a0 * i.M0
    Tt0 = air.T_from_h(air.h(T0) + 0.5 * V0 * V0) if air.kind == "varcp" \
        else T0 * (1.0 + 0.5 * (air.g0 - 1.0) * i.M0 ** 2)
    Pt0 = P0 * pr_isentropic(air, T0, Tt0)
    return T0, P0, a0, V0, Tt0, Pt0


def burner_far(air, hot_of_f, Tt3, Tt4, eta_b, hPR, neglect_fuel_mass=False):
    """Burner fuel-air ratio f_b (per kg of burner air) from
    h_air(Tt3) + f_b eta_b hPR = (1 + f_b) h_prod(Tt4, f_b)   (fixed point)."""
    h3 = air.h(Tt3)
    if neglect_fuel_mass:
        return (hot_of_f(0.02).h(Tt4) - h3) / (eta_b * hPR)
    fb = 0.02
    for _ in range(100):
        h4 = hot_of_f(fb).h(Tt4)
        fn = (h4 - h3) / (eta_b * hPR - h4)
        if abs(fn - fb) < 1e-15:
            fb = fn
            break
        fb = fn
    return fb


def design(inp: DesignInputs):
    """Parametric cycle analysis.  Returns a dict of performance, stations and
    internals.  ``valid`` is False if the cycle cannot close (turbine cannot supply
    the compressor work or a nozzle total pressure is below ambient)."""
    th = make_thermo(inp)
    air = th.air()
    T0, P0, a0, V0, Tt0, Pt0 = freestream(inp, air)
    pi_d = inp.pi_d_max * (ram_recovery(inp.M0) if inp.mil_spec_recovery else 1.0)
    nf = inp.neglect_fuel_mass
    eps = inp.eps_cool

    Tt2, Pt2 = Tt0, Pt0 * pi_d
    Tt13, Pt13 = compress_polytropic(air, Tt2, inp.pi_f, inp.e_f), Pt2 * inp.pi_f
    Tt25, Pt25 = compress_polytropic(air, Tt2, inp.pi_cL, inp.e_cL), Pt2 * inp.pi_cL
    Tt3, Pt3 = compress_polytropic(air, Tt25, inp.pi_cH, inp.e_cH), Pt25 * inp.pi_cH
    Tt4, Pt4 = inp.Tt4, Pt3 * inp.pi_b

    out = dict(inputs=inp.to_dict(), valid=False, T0=T0, P0=P0, a0=a0, V0=V0,
               tau_r=Tt0 / T0, pi_r=Pt0 / P0, pi_d=pi_d, OPR=inp.OPR)
    if Tt4 <= Tt3:
        out["reason"] = "Tt4 below compressor exit temperature"
        out["f"] = float("nan")
        return out
    fb = burner_far(air, th.gas, Tt3, Tt4, inp.eta_b, inp.hPR, nf)
    f = (1.0 - eps) * fb                      # fuel per unit core air
    mb = (1.0 - eps) * (1.0 if nf else 1.0 + fb)   # HPT flow per unit core air
    m45 = 1.0 if nf else 1.0 + f                 # LPT / nozzle flow per unit core air
    g4 = th.gas(fb)
    g45 = th.gas(f) if eps > 0 else g4
    out["f"], out["f_b"] = f, fb

    h4 = g4.h(Tt4)
    W_HPC = air.h(Tt3) - air.h(Tt25)
    h44 = h4 - W_HPC / (inp.eta_mH * mb)
    W_LP_load = (air.h(Tt25) - air.h(Tt2)) + inp.alpha * (air.h(Tt13) - air.h(Tt2))
    if fb <= 0 or h44 <= g4.h(150.0):
        out["reason"] = "HPT cannot drive the HPC"
        return out
    Tt44 = g4.T_from_h(h44)
    pi_tH = turbine_pi_polytropic(g4, Tt4, Tt44, inp.e_tH)
    Pt44 = Pt4 * pi_tH
    if eps > 0:
        h45 = (mb * h44 + eps * air.h(Tt3)) / m45
        Tt45 = g45.T_from_h(h45)
    else:
        h45, Tt45 = h44, Tt44
    Pt45 = Pt44
    h5 = h45 - W_LP_load / (inp.eta_mL * m45)
    if h5 <= g45.h(150.0):
        out["reason"] = "LPT cannot drive the fan and LPC"
        return out
    Tt5 = g45.T_from_h(h5)
    pi_tL = turbine_pi_polytropic(g45, Tt45, Tt5, inp.e_tL)
    Pt5 = Pt45 * pi_tL
    Tt9, Pt9 = Tt5, Pt5 * inp.pi_n
    Tt19, Pt19 = Tt13, Pt13 * inp.pi_fn

    n9 = nozzle(g45, Tt9, Pt9, P0, inp.core_nozzle)
    n19 = nozzle(air, Tt19, Pt19, P0, inp.fan_nozzle) if inp.alpha > 0 else dict(P=P0, T=T0, M=inp.M0, V=V0)
    if n9 is None or n19 is None:
        out["reason"] = "nozzle total pressure below ambient"
        return out

    Fc = m45 * n9["V"] - V0 + m45 * g45.R * n9["T"] * (1.0 - P0 / n9["P"]) / n9["V"]
    Fb = n19["V"] - V0 + (air.R * n19["T"] * (1.0 - P0 / n19["P"]) / n19["V"] if inp.alpha > 0 else 0.0)
    Fs = (Fc + inp.alpha * Fb) / (1.0 + inp.alpha)
    S = f / ((1.0 + inp.alpha) * Fs) if Fs > 0 else float("nan")

    Vfe9 = v_fully_expanded(g45, Tt9, Pt9, P0)
    Vfe19 = v_fully_expanded(air, Tt19, Pt19, P0) if inp.alpha > 0 else V0
    dKE = 0.5 * (m45 * Vfe9 ** 2 + inp.alpha * Vfe19 ** 2 - (1.0 + inp.alpha) * V0 ** 2)
    q_in = f * inp.hPR
    eta_th = dKE / q_in
    eta_o = (1.0 + inp.alpha) * Fs * V0 / q_in
    eta_p = eta_o / eta_th if eta_th > 0 else float("nan")

    eta = dict(
        f=compressor_eta(air, Tt2, Tt13, inp.pi_f),
        cL=compressor_eta(air, Tt2, Tt25, inp.pi_cL),
        cH=compressor_eta(air, Tt25, Tt3, inp.pi_cH),
        tH=turbine_eta(g4, Tt4, Tt44, pi_tH),
        tL=turbine_eta(g45, Tt45, Tt5, pi_tL),
    )
    st = {
        "0":   dict(Tt=Tt0, Pt=Pt0, g="air"),
        "2":   dict(Tt=Tt2, Pt=Pt2, g="air"),
        "13":  dict(Tt=Tt13, Pt=Pt13, g="air"),
        "19":  dict(Tt=Tt19, Pt=Pt19, g="air"),
        "2.5": dict(Tt=Tt25, Pt=Pt25, g="air"),
        "3":   dict(Tt=Tt3, Pt=Pt3, g="air"),
        "4":   dict(Tt=Tt4, Pt=Pt4, g="burner"),
        "4.4": dict(Tt=Tt44, Pt=Pt44, g="burner"),
        "4.5": dict(Tt=Tt45, Pt=Pt45, g="mixed"),
        "5":   dict(Tt=Tt5, Pt=Pt5, g="mixed"),
        "9":   dict(Tt=Tt9, Pt=Pt9, g="mixed"),
    }
    out.update(
        valid=bool(Fs > 0),
        Fs=Fs, S=S, S_mg=S * 1e6, S_lb=S * TSFC_LB,
        eta_th=eta_th, eta_p=eta_p, eta_o=eta_o,
        V9=n9["V"], V19=n19["V"], Vfe9=Vfe9, Vfe19=Vfe19, M9=n9["M"], M19=n19["M"],
        P9=n9["P"], P19=n19["P"], T9=n9["T"], T19=n19["T"],
        core_choked=bool(n9["P"] > P0 * (1 + 1e-12)), fan_choked=bool(n19["P"] > P0 * (1 + 1e-12)),
        Fs_core=Fc, Fs_bypass=Fb,
        tau_tH=Tt44 / Tt4, tau_tL=Tt5 / Tt45, pi_tH=pi_tH, pi_tL=pi_tL, eta_ad=eta,
        mb=mb, m45=m45, m4=m45,
        W_HPT=mb * (h4 - h44), W_LPT=m45 * (h45 - h5), W_HPC=W_HPC, W_LP_load=W_LP_load,
        dh_f=air.h(Tt13) - air.h(Tt2), dh_cL=air.h(Tt25) - air.h(Tt2),
        stations=st,
    )
    return out


def energy_residual(r):
    """First-law closure of the whole engine per unit core air flow [J/kg]:
    in  = (1+alpha) h_air(Tt0) + f eta_b hPR
    out = m45 (h_gas(T9) + V9^2/2) + alpha (h_air(T19) + V19^2/2) + shaft losses
    (sensible enthalpies, common 298.15 K reference for 'varcp'; for 'cpg' the
    two-gas model's own datum).  Returns (residual, reference magnitude)."""
    i = DesignInputs(**r["inputs"])
    th = make_thermo(i)
    air = th.air()
    g45 = th.gas(r["f"]) if i.eps_cool > 0 else th.gas(r["f_b"])
    a = i.alpha
    Tt0 = r["stations"]["0"]["Tt"]
    e_in = (1 + a) * air.h(Tt0) + r["f"] * i.eta_b * i.hPR
    loss = (1 - i.eta_mH) * r["W_HPT"] + (1 - i.eta_mL) * r["W_LPT"]
    e_out = r["m45"] * (g45.h(r["T9"]) + 0.5 * r["V9"] ** 2) + loss
    if a > 0:
        e_out += a * (air.h(r["T19"]) + 0.5 * r["V19"] ** 2)
    return e_in - e_out, e_in


def ts_path(r, n=12):
    """Temperature-entropy path through the station states, for plotting.

    Entropy relative to the free-stream static state, built leg by leg with the gas
    of that leg: ds = s0(T2) - s0(T1) - R ln(P2/P1); the burner leg uses the
    products' properties and the coolant-mixing leg (4.4 -> 4.5) shows only the
    constant-pressure temperature drop of the core gas.  Along a leg, s is
    interpolated linearly in s0(T) (exact for polytropic processes).
    Point indices on the core path: 0:0 static, 1:0, 2:2, 14:2.5, 26:3, 38:4,
    50:4.4, 51:4.5, 63:5, 64:9 static; bypass: 0:2, 12:13, 13:19 static."""
    i = DesignInputs(**r["inputs"])
    th = make_thermo(i)
    air = th.air()
    g4 = th.gas(r["f_b"])
    g45 = th.gas(r["f"]) if i.eps_cool > 0 else g4
    st = r["stations"]
    T0, P0 = r["T0"], r["P0"]

    def leg(pts, s1, T1, P1, T2, P2, gas, curved=True):
        s2 = s1 + gas.s0(T2) - gas.s0(T1) - gas.R * math.log(P2 / P1)
        k = n if curved else 1
        d = gas.s0(T2) - gas.s0(T1)
        for j in range(1, k + 1):
            T = T1 * (T2 / T1) ** (j / k)
            x = (gas.s0(T) - gas.s0(T1)) / d if d != 0 else j / k
            pts.append((s1 + x * (s2 - s1), T))
        return s2

    S = lambda k: (st[k]["Tt"], st[k]["Pt"])  # noqa: E731
    core = [(0.0, T0)]
    s = leg(core, 0.0, T0, P0, *S("0"), air, False)
    s = leg(core, s, *S("0"), *S("2"), air, False)
    s2 = s
    s = leg(core, s, *S("2"), *S("2.5"), air)
    s = leg(core, s, *S("2.5"), *S("3"), air)
    s = leg(core, s, *S("3"), *S("4"), g4)
    s = leg(core, s, *S("4"), *S("4.4"), g4)
    s = leg(core, s, *S("4.4"), *S("4.5"), g45, False)
    s = leg(core, s, *S("4.5"), *S("5"), g45)
    s = leg(core, s, *S("5"), r["T9"], r["P9"], g45, False)
    byp = [(s2, st["2"]["Tt"])]
    sb = leg(byp, s2, *S("2"), *S("13"), air)
    sb = leg(byp, sb, *S("13"), r["T19"], r["P19"], air, False)
    smax = max(s, sb) * 1.05 + 50
    amb = []
    for j in range(21):
        x = smax * j / 20
        amb.append((x, air.T_from_s0(air.s0(T0) + x) if air.kind == "varcp" else T0 * math.exp(x / air.cp0)))
    return dict(core=core, bypass=byp, ambient=amb)


def ideal_turbofan(M0, T0, gamma, cp, hPR, Tt4, pi_c, pi_f, alpha):
    """Closed-form ideal separate-exhaust turbofan (Mattingly, ideal cycle):
    isentropic components, P9 = P19 = P0, constant (gamma, cp), f << 1.
    Returns dict(Fs, S, f, eta_th, eta_p, eta_o, V9, V19, tau_f_opt)."""
    R = (gamma - 1) / gamma * cp
    a0 = math.sqrt(gamma * R * T0)
    tau_r = 1 + 0.5 * (gamma - 1) * M0 ** 2
    tau_lam = Tt4 / T0
    tau_c = pi_c ** ((gamma - 1) / gamma)
    tau_f = pi_f ** ((gamma - 1) / gamma)
    V9a = math.sqrt(2 / (gamma - 1) * (tau_lam - tau_r * (tau_c - 1 + alpha * (tau_f - 1))
                                          - tau_lam / (tau_r * tau_c)))
    V19a = math.sqrt(2 / (gamma - 1) * (tau_r * tau_f - 1))
    Fs = a0 / (1 + alpha) * (V9a - M0 + alpha * (V19a - M0))
    f = cp * T0 / hPR * (tau_lam - tau_r * tau_c)
    S = f / ((1 + alpha) * Fs)
    eta_th = 1 - 1 / (tau_r * tau_c)
    eta_p = 2 * M0 * (V9a + alpha * V19a - (1 + alpha) * M0) / (V9a ** 2 + alpha * V19a ** 2 - (1 + alpha) * M0 ** 2)
    tau_f_opt = (tau_lam - tau_r * (tau_c - 1) - tau_lam / (tau_r * tau_c) + alpha * tau_r + 1) / (tau_r * (1 + alpha))
    return dict(Fs=Fs, S=S, f=f, eta_th=eta_th, eta_p=eta_p, eta_o=eta_th * eta_p,
                V9=V9a * a0, V19=V19a * a0, tau_f_opt=tau_f_opt,
                pi_f_opt=tau_f_opt ** (gamma / (gamma - 1)))


def optimum_fan_pr(inp: DesignInputs, lo=1.02, hi=None, tol=1e-10):
    """Fan pressure ratio that minimises TSFC for the given bypass ratio (golden
    section on the feasible interval).  Returns (pi_f_opt, result dict)."""
    def S_of(pf):
        r = design(inp.with_(pi_f=pf))
        return r["S"] if r["valid"] and r["S"] > 0 else 1e9

    if hi is None:
        # largest feasible pi_f (turbine exhausted) by bisection
        a, b = lo, 8.0
        if design(inp.with_(pi_f=b))["valid"]:
            hi = b
        else:
            for _ in range(80):
                m = 0.5 * (a + b)
                if design(inp.with_(pi_f=m))["valid"]:
                    a = m
                else:
                    b = m
            hi = a
    g = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    x1, x2 = b - g * (b - a), a + g * (b - a)
    f1, f2 = S_of(x1), S_of(x2)
    while b - a > tol:
        if f1 < f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - g * (b - a)
            f1 = S_of(x1)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + g * (b - a)
            f2 = S_of(x2)
    pf = 0.5 * (a + b)
    return pf, design(inp.with_(pi_f=pf))



def marginal_optimum_ratio(r):
    """Exact first-order optimality condition for the fan pressure ratio of the
    real cycle with ideally expanded nozzles (derived in the write-up):

        V19/V9 = eta_mL * [1 - x19 (1 - e_f)] / [1 + x9 (1/e_tL - 1)],
        x19 = T19/Tt19,  x9 = T9/Tt9.

    Moving dw of work per unit bypass air from the core to the bypass stream raises
    the bypass-jet kinetic energy by dw[1 - x19(1 - e_f)] (the marginal fan
    efficiency) and lowers the core-jet kinetic energy by (dw/eta_mL)[1 + x9(1/e_tL
    - 1)] (reheat penalty of a polytropic turbine); fuel flow is unchanged.  Setting
    dF = 0 gives the ratio above; with e_f = e_tL = eta_mL = 1 it reduces to the
    ideal equal-velocity condition V19 = V9.  Mattingly's simpler approximation is
    V19/V9 ~ eta_f * eta_tL (adiabatic efficiencies)."""
    i = r["inputs"]
    x19 = r["T19"] / r["stations"]["19"]["Tt"]
    x9 = r["T9"] / r["stations"]["9"]["Tt"]
    return i["eta_mL"] * (1 - x19 * (1 - i["e_f"])) / (1 + x9 * (1 / i["e_tL"] - 1))


def cp_air_poly(T):
    """cp of air [J/(kg K)] from the cubic fit a + bT + cT^2 + dT^3 (kJ/(kmol K)),
    a = 28.11, b = 0.1967e-2, c = 0.4802e-5, d = -1.966e-9, 273-1800 K
    (Cengel & Boles, ideal-gas specific heat table; quoted max error < 1 %)."""
    return (28.11 + 0.1967e-2 * T + 0.4802e-5 * T ** 2 - 1.966e-9 * T ** 3) / 28.97 * 1000.0


def dh_air_poly(T1, T2):
    """Enthalpy rise of air between T1 and T2 from the cubic cp fit [J/kg]."""
    def h(T):
        return (28.11 * T + 0.1967e-2 * T ** 2 / 2 + 0.4802e-5 * T ** 3 / 3 - 1.966e-9 * T ** 4 / 4) / 28.97 * 1000.0
    return h(T2) - h(T1)
