"""NASA 7-coefficient thermodynamic data for the C-H-O(-N-Ar) gas-phase system.

Each species carries two coefficient sets (a1..a7) valid below and above a
common breakpoint T_mid = 1000 K:

    cp/R  = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
    h/RT  = a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T
    s0/R  = a1 ln T + a2 T + a3 T^2/2 + a4 T^3/3 + a5 T^4/4 + a7

The coefficients were typed in from the GRI-Mech 3.0 ``thermo30.dat`` file
(which itself draws on the Burcat / TPIS / JANAF compilations).  Two
deliberate deviations, both documented and tested:

* OH uses the later Burcat (ATcT-based, Ruscic et al.) fit instead of the
  GRI-Mech 3.0 RUS78 fit, because the older fit carries a heat of formation of
  about +39 kJ/mol whereas the currently accepted value is about +37.3 kJ/mol.
* N2 uses the classic Chemkin 121286 fit (300-5000 K), which is what GRI-Mech
  ships as well.

The published fits are nominally valid to 3500 K (GRI) or 5000-6000 K (N2, NO).
Rocket chamber temperatures reach ~3700 K, so the high-range polynomials are
used with mild extrapolation above 3500 K; this is a stated limitation.

Every entry is checked in ``tests/test_thermo.py`` against the standard heat of
formation and standard entropy at 298.15 K and for continuity of cp, h and s at
the 1000 K breakpoint.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

R = 8.314462618          # J/(mol K), CODATA 2018 molar gas constant
T_REF = 298.15           # K
P_REF = 1.0e5            # Pa (1 bar standard state, as in CEA)

ATOMIC_MASS = {"C": 12.0107, "H": 1.00794, "O": 15.9994, "N": 14.0067, "Ar": 39.948}


@dataclass
class Species:
    name: str
    comp: dict                    # element -> atoms per molecule
    low: tuple                    # a1..a7 for T < T_mid
    high: tuple                   # a1..a7 for T >= T_mid
    t_mid: float = 1000.0
    t_range: tuple = (200.0, 3500.0)
    source: str = "GRI-Mech 3.0"
    molar_mass: float = field(init=False)  # g/mol

    def __post_init__(self):
        self.molar_mass = sum(ATOMIC_MASS[e] * k for e, k in self.comp.items())

    def coeffs(self, T):
        return self.high if T >= self.t_mid else self.low

    def cp_R(self, T):
        a = self.coeffs(T)
        return a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4])))

    def h_RT(self, T):
        a = self.coeffs(T)
        return (a[0] + T * (a[1] / 2 + T * (a[2] / 3 + T * (a[3] / 4 + T * a[4] / 5)))
                + a[5] / T)

    def s_R(self, T):
        a = self.coeffs(T)
        return (a[0] * math.log(T) + T * (a[1] + T * (a[2] / 2 + T * (a[3] / 3 + T * a[4] / 4)))
                + a[6])

    # convenience in SI
    def h(self, T):
        """Molar enthalpy, J/mol (includes the heat of formation)."""
        return self.h_RT(T) * R * T

    def s0(self, T):
        """Standard-state molar entropy at 1 bar, J/(mol K)."""
        return self.s_R(T) * R

    def cp(self, T):
        return self.cp_R(T) * R


def _sp(name, comp, low, high, **kw):
    return Species(name, comp, tuple(low), tuple(high), **kw)


# --------------------------------------------------------------------------
#  Coefficient table.  low = 200-1000 K, high = 1000-3500 K (GRI-Mech 3.0)
# --------------------------------------------------------------------------
SPECIES_LIST = [
    _sp("H2", {"H": 2},
        [2.34433112E+00, 7.98052075E-03, -1.94781510E-05, 2.01572094E-08, -7.37611761E-12,
         -9.17935173E+02, 6.83010238E-01],
        [3.33727920E+00, -4.94024731E-05, 4.99456778E-07, -1.79566394E-10, 2.00255376E-14,
         -9.50158922E+02, -3.20502331E+00]),
    _sp("O2", {"O": 2},
        [3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12,
         -1.06394356E+03, 3.65767573E+00],
        [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14,
         -1.08845772E+03, 5.45323129E+00]),
    _sp("H", {"H": 1},
        [2.50000000E+00, 7.05332819E-13, -1.99591964E-15, 2.30081632E-18, -9.27732332E-22,
         2.54736599E+04, -4.46682853E-01],
        [2.50000001E+00, -2.30842973E-11, 1.61561948E-14, -4.73515235E-18, 4.98197357E-22,
         2.54736599E+04, -4.46682914E-01]),
    _sp("O", {"O": 1},
        [3.16826710E+00, -3.27931884E-03, 6.64306396E-06, -6.12806624E-09, 2.11265971E-12,
         2.91222592E+04, 2.05193346E+00],
        [2.56942078E+00, -8.59741137E-05, 4.19484589E-08, -1.00177799E-11, 1.22833691E-15,
         2.92175791E+04, 4.78433864E+00]),
    _sp("OH", {"O": 1, "H": 1},
        [3.99198424E+00, -2.40106655E-03, 4.61664033E-06, -3.87916306E-09, 1.36319502E-12,
         3.36889836E+03, -1.03998477E-01],
        [2.83853033E+00, 1.10741289E-03, -2.94000209E-07, 4.20698729E-11, -2.42289890E-15,
         3.69780808E+03, 5.84494652E+00],
        t_range=(200.0, 6000.0), source="Burcat (ATcT-based) update of GRI-Mech OH"),
    _sp("H2O", {"H": 2, "O": 1},
        [4.19864056E+00, -2.03643410E-03, 6.52040211E-06, -5.48797062E-09, 1.77197817E-12,
         -3.02937267E+04, -8.49032208E-01],
        [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14,
         -3.00042971E+04, 4.96677010E+00]),
    _sp("HO2", {"H": 1, "O": 2},
        [4.30179801E+00, -4.74912051E-03, 2.11582891E-05, -2.42763894E-08, 9.29225124E-12,
         2.94808040E+02, 3.71666245E+00],
        [4.01721090E+00, 2.23982013E-03, -6.33658150E-07, 1.14246370E-10, -1.07908535E-14,
         1.11856713E+02, 3.78510215E+00]),
    _sp("H2O2", {"H": 2, "O": 2},
        [4.27611269E+00, -5.42822417E-04, 1.67335701E-05, -2.15770813E-08, 8.62454363E-12,
         -1.77025821E+04, 3.43505074E+00],
        [4.16500285E+00, 4.90831694E-03, -1.90139225E-06, 3.71185986E-10, -2.87908305E-14,
         -1.78617877E+04, 2.91615662E+00]),
    _sp("CO", {"C": 1, "O": 1},
        [3.57953347E+00, -6.10353680E-04, 1.01681433E-06, 9.07005884E-10, -9.04424499E-13,
         -1.43440860E+04, 3.50840928E+00],
        [2.71518561E+00, 2.06252743E-03, -9.98825771E-07, 2.30053008E-10, -2.03647716E-14,
         -1.41518724E+04, 7.81868772E+00]),
    _sp("CO2", {"C": 1, "O": 2},
        [2.35677352E+00, 8.98459677E-03, -7.12356269E-06, 2.45919022E-09, -1.43699548E-13,
         -4.83719697E+04, 9.90105222E+00],
        [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14,
         -4.87591660E+04, 2.27163806E+00]),
    _sp("CH4", {"C": 1, "H": 4},
        [5.14987613E+00, -1.36709788E-02, 4.91800599E-05, -4.84743026E-08, 1.66693956E-11,
         -1.02466476E+04, -4.64130376E+00],
        [7.48514950E-02, 1.33909467E-02, -5.73285809E-06, 1.22292535E-09, -1.01815230E-13,
         -9.46834459E+03, 1.84373180E+01]),
    _sp("CH3", {"C": 1, "H": 3},
        [3.67359040E+00, 2.01095175E-03, 5.73021856E-06, -6.87117425E-09, 2.54385734E-12,
         1.64449988E+04, 1.60456433E+00],
        [2.28571772E+00, 7.23990037E-03, -2.98714348E-06, 5.95684644E-10, -4.67154394E-14,
         1.67755843E+04, 8.48007179E+00]),
    _sp("HCO", {"H": 1, "C": 1, "O": 1},
        [4.22118584E+00, -3.24392532E-03, 1.37799446E-05, -1.33144093E-08, 4.33768865E-12,
         3.83956496E+03, 3.39437243E+00],
        [2.77217438E+00, 4.95695526E-03, -2.48445613E-06, 5.89161778E-10, -5.33508711E-14,
         4.01191815E+03, 9.79834492E+00]),
    _sp("CH2O", {"H": 2, "C": 1, "O": 1},
        [4.79372315E+00, -9.90833369E-03, 3.73220008E-05, -3.79285261E-08, 1.31772652E-11,
         -1.43089567E+04, 6.02812900E-01],
        [1.76069008E+00, 9.20000082E-03, -4.42258813E-06, 1.00641212E-09, -8.83855640E-14,
         -1.39958323E+04, 1.36563230E+01]),
    # ---- nitrogen / argon, used only for the fuel-air validation cases ----
    _sp("N2", {"N": 2},
        [3.29867700E+00, 1.40824040E-03, -3.96322200E-06, 5.64151500E-09, -2.44485400E-12,
         -1.02089990E+03, 3.95037200E+00],
        [2.92664000E+00, 1.48797680E-03, -5.68476000E-07, 1.00970380E-10, -6.75335100E-15,
         -9.22797700E+02, 5.98052800E+00],
        t_range=(300.0, 5000.0), source="Chemkin 121286 (as shipped in GRI-Mech 3.0)"),
    _sp("NO", {"N": 1, "O": 1},
        [4.21847630E+00, -4.63897600E-03, 1.10410220E-05, -9.33613540E-09, 2.80357700E-12,
         9.84462300E+03, 2.28084640E+00],
        [3.26060560E+00, 1.19110430E-03, -4.29170480E-07, 6.94576690E-11, -4.03360990E-15,
         9.92097460E+03, 6.36930270E+00],
        t_range=(200.0, 6000.0)),
    _sp("Ar", {"Ar": 1},
        [2.5, 0.0, 0.0, 0.0, 0.0, -7.45375000E+02, 4.36600100E+00],
        [2.5, 0.0, 0.0, 0.0, 0.0, -7.45375000E+02, 4.36600100E+00],
        t_range=(300.0, 5000.0)),
]

SPECIES = {sp.name: sp for sp in SPECIES_LIST}

# Species considered for C-H-O rocket problems.  N2/NO/Ar are added only when
# the reactants contain N or Ar.
CHO_SPECIES = ["H2", "O2", "H", "O", "OH", "H2O", "HO2", "H2O2",
               "CO", "CO2", "CH4", "CH3", "HCO", "CH2O"]

# Reference values at 298.15 K used by the tests (NIST-JANAF / ATcT tables,
# rounded; kJ/mol and J/(mol K)).  These are well-established handbook values.
REFERENCE_298 = {
    #        dHf      S0
    "H2":   (0.0,     130.68),
    "O2":   (0.0,     205.15),
    "H":    (217.998, 114.72),
    "O":    (249.18,  161.06),
    "OH":   (37.3,    183.7),
    "H2O":  (-241.826, 188.83),
    "HO2":  (12.0,    229.0),
    "H2O2": (-136.1,  232.9),
    "CO":   (-110.53, 197.66),
    "CO2":  (-393.51, 213.79),
    "CH4":  (-74.6,   186.3),
    "CH3":  (146.0,   194.0),
    "HCO":  (42.0,    224.6),
    "CH2O": (-108.6,  218.8),
    "N2":   (0.0,     191.61),
    "NO":   (90.3,    210.76),
    "Ar":   (0.0,     154.85),
}


def species_table(names):
    """Return list of Species objects for the given names."""
    return [SPECIES[n] for n in names]


def to_json_dict():
    """Serialisable copy of the coefficient table (consumed by the JS port)."""
    return {sp.name: {"comp": sp.comp, "low": list(sp.low), "high": list(sp.high),
                      "tmid": sp.t_mid, "M": sp.molar_mass, "source": sp.source}
            for sp in SPECIES_LIST}
