"""tfcycle -- parametric (on-design) and off-design cycle analysis of a two-spool,
separate-exhaust turbofan, following the formulation of J. D. Mattingly,
*Elements of Propulsion: Gas Turbines and Rockets* (AIAA), with a 1976 US
Standard Atmosphere.

Modules
-------
atmosphere   US Standard Atmosphere 1976 (0-86 km)
gas          calorically perfect gas helpers (mass-flow parameter, nozzle relations)
ondesign     parametric cycle analysis (real "modified" cycle + ideal closed form)
offdesign    fixed-geometry off-design matching (choked HPT/LPT throats)
installation fan-diameter estimate and a simple nacelle-drag / weight penalty
"""
from .atmosphere import atmosphere, FT
from .gas import Gas
from .ondesign import DesignInputs, design, ideal_turbofan, optimum_fan_pr, TSFC_LB
from .offdesign import Engine

__all__ = ["atmosphere", "FT", "Gas", "DesignInputs", "design", "ideal_turbofan",
           "optimum_fan_pr", "Engine", "TSFC_LB"]
