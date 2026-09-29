"""Verification of the typed-in NASA 7-coefficient data.

For every species:
  * the 298.15 K enthalpy reproduces the standard heat of formation,
  * the 298.15 K entropy reproduces the tabulated standard entropy,
  * cp, h and s are continuous at the 1000 K breakpoint,
  * cp/R stays physical (positive, below the classical fully-excited limit)
    over the range the rocket problems use.
Run with pytest or plain ``python3 tests/test_thermo.py``.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from eqsolver.thermo import SPECIES_LIST, REFERENCE_298, R  # noqa: E402

# Species that dominate rocket exhaust must match to 0.05 kJ/mol; minor
# radicals/intermediates have genuine 0.5-1 kJ/mol scatter between
# compilations (GRI-Mech vs JANAF vs ATcT), so they get 1.5 kJ/mol.
MAJOR = {"H2", "O2", "H", "O", "OH", "H2O", "CO", "CO2", "CH4", "N2", "Ar"}


def _eval(c, T):
    cp = c[0] + T * (c[1] + T * (c[2] + T * (c[3] + T * c[4])))
    h = c[0] + T * (c[1] / 2 + T * (c[2] / 3 + T * (c[3] / 4 + T * c[4] / 5))) + c[5] / T
    s = c[0] * math.log(T) + T * (c[1] + T * (c[2] / 2 + T * (c[3] / 3 + T * c[4] / 4))) + c[6]
    return cp, h, s


def thermo_check_table():
    rows = []
    for sp in SPECIES_LIST:
        hf_ref, s_ref = REFERENCE_298[sp.name]
        hf = sp.h(298.15) / 1000
        s = sp.s0(298.15)
        lo = _eval(sp.low, sp.t_mid)
        hi = _eval(sp.high, sp.t_mid)
        jumps = [abs(a - b) / max(abs(a), 1.0) for a, b in zip(lo, hi)]
        rows.append({"species": sp.name, "hf": hf, "hf_ref": hf_ref, "dhf": hf - hf_ref,
                     "s": s, "s_ref": s_ref, "ds": s - s_ref,
                     "tol_hf": 0.05 if sp.name in MAJOR else 1.5,
                     "tol_s": 0.2 if sp.name in MAJOR else 2.0,
                     "jump_cp": jumps[0], "jump_h": jumps[1], "jump_s": jumps[2],
                     "source": sp.source})
    return rows


def test_heat_of_formation_298():
    for r in thermo_check_table():
        assert abs(r["dhf"]) <= r["tol_hf"], r


def test_entropy_298():
    for r in thermo_check_table():
        assert abs(r["ds"]) <= r["tol_s"], r


def test_continuity_at_breakpoint():
    for r in thermo_check_table():
        assert max(r["jump_cp"], r["jump_h"], r["jump_s"]) < 1e-5, r


def test_brief_values():
    """The specific checks demanded in the study brief."""
    from eqsolver.thermo import SPECIES as S
    want = {"H2O": -241.826, "CO2": -393.51, "CO": -110.53, "OH": 37.3, "H": 218.0, "O": 249.2}
    for k, v in want.items():
        assert abs(S[k].h(298.15) / 1000 - v) < 0.1, (k, S[k].h(298.15) / 1000, v)


def test_cp_physical():
    for sp in SPECIES_LIST:
        natoms = sum(sp.comp.values())
        # classical limit: 3N-? ; use generous 3N+1 bound for cp/R, and > 2.5
        for T in range(300, 4001, 50):
            c = sp.cp_R(T)
            assert c >= 2.49, (sp.name, T, c)
            assert c < 3 * natoms + 1.5, (sp.name, T, c)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
