"""Nozzle area distributions A(x), normalised so the throat area is 1 unless stated.

Each geometry is a small object with attributes ``L`` (length), ``x_throat`` and
``eps`` (exit/throat area ratio) and a vectorised ``A(x)``.
"""
from __future__ import annotations

import numpy as np


class AndersonNozzle:
    """Anderson's textbook CD nozzle, A(x) = 1 + 2.2 (x - 1.5)^2 on 0 <= x <= 3.

    Throat at x = 1.5, A_e/A* = 5.95. (Anderson, *Computational Fluid Dynamics:
    The Basics with Applications*, ch. 7.)
    """
    L = 3.0
    x_throat = 1.5
    name = "Anderson CD nozzle"

    def A(self, x):
        return 1.0 + 2.2 * (np.asarray(x, float) - 1.5) ** 2

    @property
    def eps(self):
        return float(self.A(self.L))


class AndersonShockNozzle:
    """Same convergent section; gentler divergence A = 1 + 0.2223 (x - 1.5)^2 for x > 1.5.

    A_e/A* = 1.5002. With pb/p0 = 0.6784 this is the shock-capturing example
    in Anderson's CFD text; here it is simply a convenient shock-in-nozzle test.
    """
    L = 3.0
    x_throat = 1.5
    name = "Anderson shock nozzle"

    def A(self, x):
        x = np.asarray(x, float)
        return np.where(x <= 1.5, 1.0 + 2.2 * (x - 1.5) ** 2, 1.0 + 0.2223 * (x - 1.5) ** 2)

    @property
    def eps(self):
        return float(self.A(self.L))


class BellLikeNozzle:
    """Smooth rocket-nozzle area law for a prescribed expansion ratio.

    Convergent: A = 1 + (CR - 1) * s^2, s = (x_t - x)/x_t   (parabolic, CR = contraction ratio)
    Divergent:  A = 1 + (eps - 1) * (3 t^2 - 2 t^3), t = (x - x_t)/L_d
    The divergent law has zero slope at the throat and at the exit, a crude
    stand-in for a bell. In quasi-1D inviscid flow the exit state depends only
    on eps, so the contour details matter only for the grid, not the answer.
    Length unit: the divergent length L_d = 1.
    """
    name = "bell-like"

    def __init__(self, eps, CR=4.0, Lc=0.3, Ld=1.0):
        self.eps_ = float(eps)
        self.CR = float(CR)
        self.x_throat = float(Lc)
        self.Ld = float(Ld)
        self.L = self.x_throat + self.Ld

    def A(self, x):
        x = np.asarray(x, float)
        s = np.clip((self.x_throat - x) / self.x_throat, 0.0, None)
        t = np.clip((x - self.x_throat) / self.Ld, 0.0, 1.0)
        conv = 1.0 + (self.CR - 1.0) * s * s
        div = 1.0 + (self.eps_ - 1.0) * (3 * t * t - 2 * t ** 3)
        return np.where(x <= self.x_throat, conv, div)

    @property
    def eps(self):
        return self.eps_
