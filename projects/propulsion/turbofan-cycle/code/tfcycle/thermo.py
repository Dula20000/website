"""Gas property models.

Two interchangeable gas models share one interface (``h``, ``cp``, ``s0``, ``R``,
``gamma``, ``T_from_h``, ``T_from_s0``):

``CPG``      calorically perfect gas, h = cp T, s0 = cp ln T.  Used for the ideal
             cycle and for Mattingly's two-gas "modified" cycle (cold gas for
             f = 0, hot gas for f > 0).
``Mixture``  thermally perfect ideal-gas mixture with temperature-dependent cp from
             NASA 7-coefficient polynomials (GRI-Mech 3.0 thermodynamic data set,
             200-1000 K / 1000-3500 K) for N2, O2, Ar, CO2 and H2O.

Enthalpies are *sensible* enthalpies measured from 298.15 K for every species, so the
burner balance with a lower heating value hPR referred to 298.15 K is consistent:

    h_air(Tt3) + f_b (h_fuel + eta_b hPR) = (1 + f_b) h_prod(Tt4, f_b),

with the fuel entering at 298.15 K (h_fuel = 0).  s0(T) is the standard-state
entropy function (J/(kg K)); along an isentrope of fixed composition
s0(T2) - s0(T1) = R ln(P2/P1).  Mixing entropy is constant for a fixed composition
and drops out.

Combustion products: complete lean combustion of a CH_y fuel (y = 23/12, a C12H23
Jet-A surrogate) in dry air (N2 0.78084, O2 0.20946, Ar 0.00934, CO2 0.00036 by mole).

To keep h and s0 continuous (and T(h), T(s0) single valued) the high-temperature
branch of each species is shifted by a constant so that h and s0 match the
low-temperature branch at 1000 K; the shifts are below 0.1 % of h.
"""
import math

RU = 8314.462618          # J/(kmol K)
T_REF = 298.15
T_MID = 1000.0

# species: molar mass [kg/kmol], low (200-1000 K) and high (1000-3500 K) a1..a7
NASA = {
    "N2": (28.0134,
           [3.298677, 1.4082404e-3, -3.963222e-6, 5.641515e-9, -2.444854e-12, -1020.8999, 3.950372],
           [2.92664, 1.4879768e-3, -5.68476e-7, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528]),
    "O2": (31.9988,
           [3.78245636, -2.99673416e-3, 9.84730201e-6, -9.68129509e-9, 3.24372837e-12, -1063.94356, 3.65767573],
           [3.28253784, 1.48308754e-3, -7.57966669e-7, 2.09470555e-10, -2.16717794e-14, -1088.45772, 5.45323129]),
    "Ar": (39.948,
           [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366],
           [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366]),
    "CO2": (44.0095,
            [2.35677352, 8.98459677e-3, -7.12356269e-6, 2.45919022e-9, -1.43699548e-13, -48371.9697, 9.90105222],
            [3.85746029, 4.41437026e-3, -2.21481404e-6, 5.23490188e-10, -4.72084164e-14, -48759.166, 2.27163806]),
    "H2O": (18.01528,
            [4.19864056, -2.0364341e-3, 6.52040211e-6, -5.48797062e-9, 1.77197817e-12, -30293.7267, -0.849032208],
            [3.03399249, 2.17691804e-3, -1.64072518e-7, -9.7041987e-11, 1.68200992e-14, -30004.2971, 4.96677010]),
}
SPECIES = ["N2", "O2", "Ar", "CO2", "H2O"]
AIR_X = {"N2": 0.78084, "O2": 0.20946, "Ar": 0.00934, "CO2": 0.00036}
M_AIR = sum(AIR_X[k] * NASA[k][0] for k in AIR_X)
FUEL_Y = 23.0 / 12.0                       # H/C atom ratio of CH_y
M_FUEL = 12.011 + FUEL_Y * 1.00794         # kg per kmol of CH_y


def _h_RT(a, T):
    return a[0] + a[1] * T / 2 + a[2] * T ** 2 / 3 + a[3] * T ** 3 / 4 + a[4] * T ** 4 / 5 + a[5] / T


def _s_R(a, T):
    return a[0] * math.log(T) + a[1] * T + a[2] * T ** 2 / 2 + a[3] * T ** 3 / 3 + a[4] * T ** 4 / 4 + a[6]


def _cp_R(a, T):
    return a[0] + a[1] * T + a[2] * T ** 2 + a[3] * T ** 3 + a[4] * T ** 4


def _continuous_high(lo, hi):
    """High branch with a6, a7 shifted so h and s0 are continuous at 1000 K."""
    hi = list(hi)
    hi[5] += (_h_RT(lo, T_MID) - _h_RT(hi, T_MID)) * T_MID
    hi[6] += _s_R(lo, T_MID) - _s_R(hi, T_MID)
    return hi


COEF = {k: (NASA[k][0], list(NASA[k][1]), _continuous_high(NASA[k][1], NASA[k][2])) for k in SPECIES}


def species_props(name, T):
    """(cp, h - h(298.15), s0) of one species in J/(kmol K), J/kmol, J/(kmol K);
    uses the *unshifted* published coefficients (for validation against tables)."""
    M, lo, hi = NASA[name]
    a = lo if T <= T_MID else hi
    h298 = _h_RT(lo, T_REF) * T_REF
    return _cp_R(a, T) * RU, (_h_RT(a, T) * T - h298) * RU, _s_R(a, T) * RU


def composition(f):
    """kmol of each species per kg of mixture for fuel-air ratio f (lean)."""
    n = {k: AIR_X.get(k, 0.0) / M_AIR for k in SPECIES}
    nC = f / M_FUEL
    n["CO2"] += nC
    n["H2O"] += nC * FUEL_Y / 2
    n["O2"] -= nC * (1 + FUEL_Y / 4)
    if n["O2"] < 0:
        raise ValueError("fuel-air ratio above stoichiometric")
    return {k: v / (1 + f) for k, v in n.items()}


F_STOICH = (AIR_X["O2"] / M_AIR) / (1 + FUEL_Y / 4) * M_FUEL


class CPG:
    """Calorically perfect gas: h = cp T (datum 0 K), s0 = cp ln T."""
    kind = "cpg"

    def __init__(self, gamma, cp):
        self.cp0 = float(cp)
        self.g0 = float(gamma)
        self.R = (self.g0 - 1) / self.g0 * self.cp0

    def h(self, T):
        return self.cp0 * T

    def cp(self, T):
        return self.cp0

    def gamma(self, T):
        return self.g0

    def s0(self, T):
        return self.cp0 * math.log(T)

    def T_from_h(self, h):
        return h / self.cp0

    def T_from_s0(self, s):
        return math.exp(s / self.cp0)

    def choke_T(self, Tt):
        return 2 * Tt / (self.g0 + 1)


class Mixture:
    """Thermally perfect mixture, NASA polynomials, per kg of mixture."""
    kind = "varcp"

    def __init__(self, f=0.0):
        self.f = f
        n = composition(f)
        self.lo = [0.0] * 7
        self.hi = [0.0] * 7
        for k in SPECIES:
            _, lo, hi = COEF[k]
            for j in range(7):
                self.lo[j] += n[k] * RU * lo[j]
                self.hi[j] += n[k] * RU * hi[j]
        self.R = RU * sum(n.values())
        self.h_ref = _h_RT(self.lo, T_REF) * T_REF

    def _a(self, T):
        return self.lo if T <= T_MID else self.hi

    def h(self, T):
        return _h_RT(self._a(T), T) * T - self.h_ref

    def cp(self, T):
        return _cp_R(self._a(T), T)

    def gamma(self, T):
        c = self.cp(T)
        return c / (c - self.R)

    def s0(self, T):
        return _s_R(self._a(T), T)

    def T_from_h(self, h):
        T = T_REF + h / 1150.0
        T = min(max(T, 150.0), 4000.0)
        for _ in range(60):
            dT = (self.h(T) - h) / self.cp(T)
            T -= dT
            if abs(dT) < 1e-12 * T:
                break
        return T

    def T_from_s0(self, s, guess=800.0):
        T = guess
        for _ in range(60):
            d = (self.s0(T) - s) / self.cp(T)      # Newton step in ln T
            T *= math.exp(-d)
            if abs(d) < 1e-13:
                break
        return T

    def choke_T(self, Tt):
        """Static temperature at M = 1 for an isentropic expansion from Tt:
        2 (h(Tt) - h(T)) = gamma(T) R T  (Newton)."""
        ht = self.h(Tt)
        T = 2 * Tt / (self.gamma(Tt) + 1)
        for _ in range(60):
            g = self.gamma(T)
            F = 2 * (ht - self.h(T)) - g * self.R * T
            dT = F / (2 * self.cp(T) + g * self.R)
            T += dT
            if abs(dT) < 1e-12 * T:
                break
        return T


# ---------------------------------------------------------------------------
# helpers built on the common interface
# ---------------------------------------------------------------------------
def T_isentropic(gas, T1, PR):
    """Temperature after an isentropic change with pressure ratio PR = P2/P1."""
    return gas.T_from_s0(gas.s0(T1) + gas.R * math.log(PR)) if gas.kind == "varcp" \
        else T1 * PR ** (gas.R / gas.cp0)


def pr_isentropic(gas, T1, T2):
    """P2/P1 of an isentropic change from T1 to T2."""
    if gas.kind == "cpg":
        return (T2 / T1) ** (gas.cp0 / gas.R)
    return math.exp((gas.s0(T2) - gas.s0(T1)) / gas.R)


def compress_polytropic(gas, T1, pi, e):
    """Exit temperature of a polytropic compression (efficiency e)."""
    if gas.kind == "cpg":
        return T1 * pi ** (gas.R / (gas.cp0 * e))
    return gas.T_from_s0(gas.s0(T1) + gas.R * math.log(pi) / e, guess=T1 * pi ** 0.3)


def turbine_pi_polytropic(gas, T1, T2, e):
    """Pressure ratio of a polytropic expansion T1 -> T2 (efficiency e)."""
    if gas.kind == "cpg":
        return (T2 / T1) ** (gas.cp0 / (gas.R * e))
    return math.exp((gas.s0(T2) - gas.s0(T1)) / (e * gas.R))


def compressor_eta(gas, T1, T2, pi):
    """Adiabatic efficiency of a compression T1 -> T2 with pressure ratio pi."""
    if T2 == T1:
        return 1.0
    T2s = T_isentropic(gas, T1, pi)
    return (gas.h(T2s) - gas.h(T1)) / (gas.h(T2) - gas.h(T1))


def turbine_eta(gas, T1, T2, pi):
    if T2 == T1:
        return 1.0
    T2s = T_isentropic(gas, T1, pi)
    return (gas.h(T1) - gas.h(T2)) / (gas.h(T1) - gas.h(T2s))


def compressor_pi_from_eta(gas, T1, dh, eta):
    """Pressure ratio for specific work dh at adiabatic efficiency eta."""
    T2s = gas.T_from_h(gas.h(T1) + eta * dh)
    return pr_isentropic(gas, T1, T2s)


def turbine_pi_from_eta(gas, T1, dh, eta):
    """Pressure ratio (<1) for extracted work dh at adiabatic efficiency eta."""
    T2s = gas.T_from_h(gas.h(T1) - dh / eta)
    return pr_isentropic(gas, T1, T2s)


def mach_from_state(gas, Tt, T):
    V = math.sqrt(max(0.0, 2 * (gas.h(Tt) - gas.h(T))))
    return V / math.sqrt(gas.gamma(T) * gas.R * T), V


def nozzle(gas, Tt, Pt, P0, kind):
    """Nozzle exit state (P, T, M, V) from total state; convergent nozzles choke
    when the M = 1 static pressure exceeds P0.  None if Pt <= P0."""
    if Pt <= P0 * (1 + 1e-12):
        return None
    P = P0
    if kind == "convergent":
        Tc = gas.choke_T(Tt)
        Pc = Pt * pr_isentropic(gas, Tt, Tc)
        if Pc > P0:
            P = Pc
            V = math.sqrt(2 * (gas.h(Tt) - gas.h(Tc)))
            return dict(P=P, T=Tc, M=1.0, V=V)
    T = T_isentropic(gas, Tt, P / Pt)
    M, V = mach_from_state(gas, Tt, T)
    return dict(P=P, T=T, M=M, V=V)


def mass_flux(gas, Tt, Pt, P):
    """Mass flow per unit area rho V at static pressure P (limited to choking)."""
    Tc = gas.choke_T(Tt)
    Pc = Pt * pr_isentropic(gas, Tt, Tc)
    if P <= Pc:
        T = Tc
        P = Pc
    else:
        T = T_isentropic(gas, Tt, P / Pt)
    V = math.sqrt(max(0.0, 2 * (gas.h(Tt) - gas.h(T))))
    return P / (gas.R * T) * V


def v_fully_expanded(gas, Tt, Pt, P0):
    T = T_isentropic(gas, Tt, P0 / Pt)
    return math.sqrt(max(0.0, 2 * (gas.h(Tt) - gas.h(T))))


class Thermo:
    """Gas-model selector.  mode 'varcp': NASA-polynomial air/products;
    mode 'cpg': Mattingly two-gas model (cold gas for f == 0, hot gas for f > 0)."""

    def __init__(self, mode="varcp", gamma_c=1.4, cp_c=1004.0, gamma_t=1.33, cp_t=1156.0):
        self.mode = mode
        if mode == "cpg":
            self.cold = CPG(gamma_c, cp_c)
            self.hot = CPG(gamma_t, cp_t)
        elif mode == "varcp":
            self.cold = Mixture(0.0)
        else:
            raise ValueError(mode)

    def air(self):
        return self.cold

    def gas(self, f):
        if self.mode == "cpg":
            return self.hot
        return Mixture(f) if f > 0 else self.cold
