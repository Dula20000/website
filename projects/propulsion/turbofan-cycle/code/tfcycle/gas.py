"""Calorically perfect gas with fixed (gamma, cp): the building block of Mattingly's
"modified" cycle analysis, which uses one (gamma_c, cp_c) pair upstream of the
burner and another (gamma_t, cp_t) pair downstream of it."""
import math


class Gas:
    def __init__(self, gamma, cp):
        self.gamma = float(gamma)
        self.cp = float(cp)
        self.R = (self.gamma - 1.0) / self.gamma * self.cp

    # --- isentropic helpers -------------------------------------------------
    def T_ratio(self, M):
        """Tt/T at Mach M."""
        return 1.0 + 0.5 * (self.gamma - 1.0) * M * M

    def P_ratio(self, M):
        """Pt/P at Mach M."""
        return self.T_ratio(M) ** (self.gamma / (self.gamma - 1.0))

    def M_from_PtP(self, PtP):
        """Mach number that gives total-to-static pressure ratio PtP (>=1)."""
        g = self.gamma
        return math.sqrt(max(0.0, 2.0 / (g - 1.0) * (PtP ** ((g - 1.0) / g) - 1.0)))

    def MFP(self, M):
        """Mass-flow parameter  mdot*sqrt(Tt)/(Pt*A)  [kg K^0.5 /(N s)]."""
        g = self.gamma
        return M * math.sqrt(g / self.R) * self.T_ratio(M) ** (-(g + 1.0) / (2.0 * (g - 1.0)))

    @property
    def MFP_star(self):
        return self.MFP(1.0)

    @property
    def crit_PR(self):
        """Pt/P* for a choked (M = 1) throat."""
        return self.P_ratio(1.0)
