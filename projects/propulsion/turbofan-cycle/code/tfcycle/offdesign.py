"""Off-design performance of a fixed-geometry two-spool separate-exhaust turbofan.

Method (Mattingly, *Elements of Propulsion*, off-design / "performance" analysis),
written for either gas model of ``thermo.py``:

* component adiabatic efficiencies (eta_f, eta_cL, eta_cH, eta_tH, eta_tL, computed
  from the design-point polytropic values), total-pressure ratios, burner and
  mechanical efficiencies and the cooling-air fraction are held at design values;
* the HPT nozzle (station 4) and the LPT nozzle (station 4.5, after coolant mixing)
  throats are choked and of fixed area; the core (8) and fan (18) nozzles are
  convergent with fixed throat areas, choked or unchoked;
* the LPC work tracks the fan work on the shared shaft, dh_cL = K dh_f with K fixed
  at design (Mattingly's tau_cL = tau_f is K = 1 for a calorically perfect gas);
* the unknown fan work dh_f closes the LP-spool power balance (Brent).  Inside each
  evaluation the HP spool, burner and HPT are converged by fixed-point iteration on
  Tt3 with a Brent solve for the HPT exit temperature that passes the flow through
  the choked 4.5 throat; the LPT exit temperature comes from a Brent solve of the
  core-nozzle flow match; the bypass ratio from the fan-nozzle flow.

With a calorically perfect gas and no cooling this reduces exactly to Mattingly's
fixed tau_tH, pi_tH.  Throttle is set by Tt4.  Absolute size enters only through the
design mass flow; ratios (thrust lapse, TSFC) are size independent.
"""
import math
from scipy.optimize import brentq

from .ondesign import design, DesignInputs, ram_recovery, TSFC_LB, make_thermo, freestream, burner_far
from .thermo import (nozzle, mass_flux, v_fully_expanded, compressor_pi_from_eta, turbine_pi_from_eta)

XTOL = 1e-14
RTOL = 1e-15


class Engine:
    def __init__(self, inp: DesignInputs, F_design=None, mdot0=None):
        if inp.core_nozzle != "convergent" or inp.fan_nozzle != "convergent":
            raise ValueError("off-design model assumes fixed convergent nozzles")
        if inp.neglect_fuel_mass:
            raise ValueError("off-design model needs the real fuel-mass terms")
        r = design(inp)
        if not r["valid"]:
            raise ValueError("design point not valid: " + r.get("reason", ""))
        self.inp, self.des = inp, r
        self.th = make_thermo(inp)
        self.air = self.th.air()
        if F_design is not None:
            mdot0 = F_design / r["Fs"]
        self.mdot0_des = mdot0 if mdot0 is not None else 100.0
        self.F_des = self.mdot0_des * r["Fs"]
        st = r["stations"]
        mc = self.mdot0_des / (1 + inp.alpha)
        g4 = self.th.gas(r["f_b"])
        g45 = self.th.gas(r["f"]) if inp.eps_cool > 0 else g4
        P0 = r["P0"]
        self.A4 = mc * r["mb"] / mass_flux(g4, st["4"]["Tt"], st["4"]["Pt"], 0.0)
        self.A45 = mc * r["m45"] / mass_flux(g45, st["4.5"]["Tt"], st["4.5"]["Pt"], 0.0)
        self.A8 = mc * r["m45"] / mass_flux(g45, st["9"]["Tt"], st["9"]["Pt"], P0)
        self.A18 = inp.alpha * mc / mass_flux(self.air, st["19"]["Tt"], st["19"]["Pt"], P0)
        self.eta = r["eta_ad"]
        self.K = r["dh_cL"] / r["dh_f"]
        self.ref = dict(dh_f=r["dh_f"], dh_cH=r["W_HPC"], Tt2=st["2"]["Tt"], Tt25=st["2.5"]["Tt"],
                        Pt2=st["2"]["Pt"], Tt3=st["3"]["Tt"])

    # ------------------------------------------------------------------
    def _hp_spool(self, Tt25, Tt4, s):
        """Converge HPC work, burner and HPT for given LPC exit and Tt4."""
        i, air, e = self.inp, self.air, self.eta
        eps = i.eps_cool
        Tt3 = s.get("Tt3", self.ref["Tt3"] * Tt25 / self.ref["Tt25"])
        h25 = air.h(Tt25)
        for _ in range(200):
            if Tt3 >= Tt4:
                return None
            fb = burner_far(air, self.th.gas, Tt3, Tt4, i.eta_b, i.hPR)
            f = (1 - eps) * fb
            mb, m45 = (1 - eps) * (1 + fb), 1 + f
            g4 = self.th.gas(fb)
            g45 = self.th.gas(f) if eps > 0 else g4
            h4 = g4.h(Tt4)
            h3 = air.h(Tt3)
            phi4 = mass_flux(g4, Tt4, 1.0, 0.0)
            target = self.A4 * phi4 * m45 / mb

            def mix(Tt44):
                if eps == 0:
                    return Tt44
                return g45.T_from_h((mb * g4.h(Tt44) + eps * h3) / m45)

            def G(Tt44):
                pi = turbine_pi_from_eta(g4, Tt4, h4 - g4.h(Tt44), e["tH"])
                return self.A45 * pi * mass_flux(g45, mix(Tt44), 1.0, 0.0) - target

            lo, hi = 0.35 * Tt4, Tt4 * (1 - 1e-12)
            if G(hi) < 0 or G(lo) > 0:
                return None
            Tt44 = brentq(G, lo, hi, xtol=XTOL, rtol=RTOL, maxiter=200)
            h3n = h25 + i.eta_mH * mb * (h4 - g4.h(Tt44))
            Tt3n = air.T_from_h(h3n)
            done = abs(Tt3n - Tt3) < 1e-11 * Tt3
            Tt3 = Tt3n
            if done:
                break
        # final consistent pass at the converged Tt3
        fb = burner_far(air, self.th.gas, Tt3, Tt4, i.eta_b, i.hPR)
        f = (1 - eps) * fb
        mb, m45 = (1 - eps) * (1 + fb), 1 + f
        g4 = self.th.gas(fb)
        g45 = self.th.gas(f) if eps > 0 else g4
        pi_tH = turbine_pi_from_eta(g4, Tt4, g4.h(Tt4) - g4.h(Tt44), e["tH"])
        Tt45 = Tt44 if eps == 0 else g45.T_from_h((mb * g4.h(Tt44) + eps * air.h(Tt3)) / m45)
        s["Tt3"] = Tt3
        return dict(Tt3=Tt3, fb=fb, f=f, mb=mb, m45=m45, g4=g4, g45=g45, Tt44=Tt44, Tt45=Tt45,
                    pi_tH=pi_tH, W_HPC=air.h(Tt3) - h25)

    def _evaluate(self, dh_f, s):
        """March through the engine for fan work dh_f; returns the state with the LP
        power-balance residual ``res`` (supply - demand, J/kg core air)."""
        i, air, e = self.inp, self.air, self.eta
        Tt2, Pt2, P0, Tt4 = s["Tt2"], s["Pt2"], s["P0"], s["Tt4"]
        h2 = air.h(Tt2)
        dh_cL = self.K * dh_f
        Tt13 = air.T_from_h(h2 + dh_f)
        Tt25 = air.T_from_h(h2 + dh_cL)
        pi_f = compressor_pi_from_eta(air, Tt2, dh_f, e["f"])
        pi_cL = compressor_pi_from_eta(air, Tt2, dh_cL, e["cL"])
        hp = self._hp_spool(Tt25, Tt4, s)
        if hp is None:
            return None
        pi_cH = compressor_pi_from_eta(air, Tt25, hp["W_HPC"], e["cH"])
        Pt25 = Pt2 * pi_cL
        Pt3 = Pt25 * pi_cH
        Pt4 = Pt3 * i.pi_b
        g4, g45, m45 = hp["g4"], hp["g45"], hp["m45"]
        mdot4 = self.A4 * mass_flux(g4, Tt4, Pt4, 0.0)
        mcore = mdot4 / hp["mb"]
        Pt45 = Pt4 * hp["pi_tH"]
        Tt45 = hp["Tt45"]
        h45 = g45.h(Tt45)
        mflow = mcore * m45

        def g(Tt5):
            pi = turbine_pi_from_eta(g45, Tt45, h45 - g45.h(Tt5), e["tL"])
            Pt9 = Pt45 * pi * i.pi_n
            if Pt9 <= P0:
                return -mflow
            return self.A8 * mass_flux(g45, Tt5, Pt9, P0) - mflow

        lo = g45.T_from_h(h45 - e["tL"] * (h45 - g45.h(0.3 * Tt45)) * 0.999)
        lo = max(lo, 0.3 * Tt45)
        hi = Tt45 * (1 - 1e-12)
        if g(hi) < 0 or g(lo) > 0:
            return None
        Tt5 = brentq(g, lo, hi, xtol=XTOL, rtol=RTOL, maxiter=200)
        pi_tL = turbine_pi_from_eta(g45, Tt45, h45 - g45.h(Tt5), e["tL"])
        Pt5 = Pt45 * pi_tL
        Pt13 = Pt2 * pi_f
        Pt19 = Pt13 * i.pi_fn
        mb_ = self.A18 * mass_flux(air, Tt13, Pt19, P0) if Pt19 > P0 else 0.0
        alpha = mb_ / mcore
        supply = i.eta_mL * m45 * (h45 - g45.h(Tt5))
        demand = dh_cL + alpha * dh_f
        return dict(dh_f=dh_f, pi_f=pi_f, pi_cL=pi_cL, pi_cH=pi_cH, pi_tL=pi_tL, f=hp["f"], fb=hp["fb"],
                    m45=m45, mb=hp["mb"], alpha=alpha, mcore=mcore, mdot4=mdot4, mbyp=mb_, g45=g45,
                    Tt25=Tt25, Tt3=hp["Tt3"], Pt3=Pt3, Pt4=Pt4, Tt44=hp["Tt44"], Tt45=Tt45, Pt45=Pt45,
                    Tt5=Tt5, Pt5=Pt5, Tt13=Tt13, Pt13=Pt13, W_HPC=hp["W_HPC"],
                    res=supply - demand, supply=supply, demand=demand)

    def operate(self, M0, alt, Tt4, T0=None, P0=None):
        """Operating point at flight Mach ``M0``, geopotential altitude ``alt`` [m] and
        burner exit temperature ``Tt4`` [K].  Returns a dict (``valid`` False if no
        matched solution exists, e.g. at very low throttle)."""
        i, air = self.inp, self.air
        fs = i.with_(M0=M0, alt=alt, T0=T0, P0=P0)
        T0, P0, a0, V0, Tt0, Pt0 = freestream(fs, air)
        pi_d = i.pi_d_max * (ram_recovery(M0) if i.mil_spec_recovery else 1.0)
        s = dict(Tt2=Tt0, Pt2=Pt0 * pi_d, Tt4=Tt4, P0=P0)
        bad = dict(valid=False, M0=M0, alt=alt, Tt4=Tt4)
        ref = self.ref["dh_f"]

        def R(x):
            ev = self._evaluate(x * ref, s)
            return -1e9 if ev is None else ev["res"] / ref

        lo = 1e-7
        if R(lo) <= 0:
            return bad
        hi = 1.5
        k = 0
        while R(hi) > 0:
            hi *= 1.5
            k += 1
            if k > 40:
                return bad
        x = brentq(R, lo, hi, xtol=XTOL, rtol=RTOL, maxiter=300)
        ev = self._evaluate(x * ref, s)
        if ev is None:
            return bad
        g45 = ev["g45"]
        m45, alpha, f = ev["m45"], ev["alpha"], ev["f"]
        Pt9 = ev["Pt5"] * i.pi_n
        Pt19 = ev["Pt13"] * i.pi_fn
        n9 = nozzle(g45, ev["Tt5"], Pt9, P0, "convergent")
        n19 = nozzle(air, ev["Tt13"], Pt19, P0, "convergent")
        if n9 is None:
            return bad
        Fc = m45 * n9["V"] - V0 + m45 * g45.R * n9["T"] * (1 - P0 / n9["P"]) / n9["V"]
        Fb = 0.0 if n19 is None else n19["V"] - V0 + air.R * n19["T"] * (1 - P0 / n19["P"]) / n19["V"]
        Fs = (Fc + alpha * Fb) / (1 + alpha)
        mdot0 = ev["mcore"] * (1 + alpha)
        F = mdot0 * Fs
        S = f / ((1 + alpha) * Fs) if Fs > 0 else float("nan")
        V19 = 0.0 if n19 is None else n19["V"]
        Vfe9 = v_fully_expanded(g45, ev["Tt5"], Pt9, P0)
        Vfe19 = 0.0 if n19 is None else v_fully_expanded(air, ev["Tt13"], Pt19, P0)
        dKE = 0.5 * (m45 * Vfe9 ** 2 + alpha * Vfe19 ** 2 - (1 + alpha) * V0 ** 2)
        eta_th = dKE / (f * i.hPR)
        eta_o = (1 + alpha) * Fs * V0 / (f * i.hPR)
        rf = self.ref
        Tt2, Pt2 = s["Tt2"], s["Pt2"]
        m45_chk = self.A45 * mass_flux(g45, ev["Tt45"], ev["Pt45"], 0.0)
        m8_chk = self.A8 * mass_flux(g45, ev["Tt5"], Pt9, P0)
        return dict(
            valid=bool(Fs > 0), M0=M0, alt=alt, T0=T0, P0=P0, V0=V0, Tt4=Tt4,
            F=F, Fs=Fs, mdot0=mdot0, S=S, S_mg=S * 1e6, S_lb=S * TSFC_LB, f=f,
            alpha=alpha, pi_f=ev["pi_f"], pi_cL=ev["pi_cL"], pi_cH=ev["pi_cH"],
            OPR=ev["pi_cL"] * ev["pi_cH"], tau_tL=ev["Tt5"] / ev["Tt45"],
            eta_th=eta_th, eta_o=eta_o, eta_p=eta_o / eta_th if eta_th > 0 else float("nan"),
            V9=n9["V"], V19=V19, Vfe9=Vfe9, Vfe19=Vfe19, M9=n9["M"], M19=0.0 if n19 is None else n19["M"],
            core_choked=bool(n9["P"] > P0 * (1 + 1e-12)),
            fan_choked=bool(n19 is not None and n19["P"] > P0 * (1 + 1e-12)),
            NL_corr=math.sqrt(ev["dh_f"] / Tt2 / (rf["dh_f"] / rf["Tt2"])),
            NL_phys=math.sqrt(ev["dh_f"] / rf["dh_f"]),
            NH_corr=math.sqrt(ev["W_HPC"] / ev["Tt25"] / (rf["dh_cH"] / rf["Tt25"])),
            NH_phys=math.sqrt(ev["W_HPC"] / rf["dh_cH"]),
            mc2_rel=(mdot0 * math.sqrt(Tt2 / 288.15) / (Pt2 / 101325.0)) /
                    (self.mdot0_des * math.sqrt(rf["Tt2"] / 288.15) / (rf["Pt2"] / 101325.0)),
            lp_residual=ev["res"] / ev["demand"],
            mass_err_45=m45_chk / (ev["mcore"] * m45) - 1, mass_err_8=m8_chk / (ev["mcore"] * m45) - 1,
        )
