"""2-D grain cross-section regression by the fast-marching level-set method.

The port (gas region) of one axial grain segment is described by an implicit
function phi0(x, y): negative inside the port, positive in the propellant, and
equal to the signed distance near the port boundary.  With a burn rate that is
uniform over the cross-section, the burning surface after a web distance w has
been consumed is the level set {T = w} of the arrival-time field T, which obeys
the eikonal equation |grad T| = 1 with T = phi0 on the initial surface.  T is
computed once per geometry with a second-order fast-marching method (FMM); the
burning perimeter P(w) and port area A(w) are then read off the field by
marching squares, clipped to the (inhibited) outer case radius R.

Primitives: BATES circular core, N-point star polygon, finocyl-style slotted
core ("wagon wheel"), and moon-burner (offset circular core).  Union of
primitives uses min(), which is exact on the propellant side where it matters.
"""
from __future__ import annotations

import heapq
import math

import numpy as np

# --------------------------------------------------------------------------
# signed-distance primitives (negative inside the port)
# --------------------------------------------------------------------------

def sd_circle(X, Y, cx, cy, r):
    return np.hypot(X - cx, Y - cy) - r


def sd_slot(X, Y, theta, r0, r1, width):
    """Signed distance to a rectangle of given width running radially from r0 to r1 at angle theta."""
    c, s = math.cos(theta), math.sin(theta)
    u = X * c + Y * s - 0.5 * (r0 + r1)
    v = -X * s + Y * c
    qx = np.abs(u) - 0.5 * (r1 - r0)
    qy = np.abs(v) - 0.5 * width
    outside = np.hypot(np.maximum(qx, 0.0), np.maximum(qy, 0.0))
    inside = np.minimum(np.maximum(qx, qy), 0.0)
    return outside + inside


def star_vertices(points, r_tip, r_valley):
    """Vertices of an N-point star: tips at 2*pi*k/N, valleys half-way between."""
    v = []
    for k in range(points):
        a = 2 * math.pi * k / points
        v.append((r_tip * math.cos(a), r_tip * math.sin(a)))
        b = a + math.pi / points
        v.append((r_valley * math.cos(b), r_valley * math.sin(b)))
    return np.array(v)


def sd_polygon(X, Y, verts):
    """Exact signed distance to a simple polygon (negative inside)."""
    d2 = np.full(X.shape, np.inf)
    inside = np.zeros(X.shape, dtype=bool)
    n = len(verts)
    for i in range(n):
        ax, ay = verts[i]
        bx, by = verts[(i + 1) % n]
        ex, ey = bx - ax, by - ay
        wx, wy = X - ax, Y - ay
        t = np.clip((wx * ex + wy * ey) / (ex * ex + ey * ey), 0.0, 1.0)
        dx, dy = wx - t * ex, wy - t * ey
        d2 = np.minimum(d2, dx * dx + dy * dy)
        # even-odd crossing test
        cond = (ay > Y) != (by > Y)
        xint = ax + (Y - ay) * ex / (ey if ey != 0 else 1e-300)
        inside ^= cond & (X < xint)
    d = np.sqrt(d2)
    return np.where(inside, -d, d)


def port_sdf(geom: dict, X, Y):
    """Implicit port function phi0 for a geometry spec (SI units, metres)."""
    t = geom["type"]
    if t == "bates":
        return sd_circle(X, Y, 0.0, 0.0, 0.5 * geom["core_d"])
    if t == "moon":
        return sd_circle(X, Y, geom["offset"], 0.0, 0.5 * geom["core_d"])
    if t == "star":
        return sd_polygon(X, Y, star_vertices(geom["points"], geom["r_tip"], geom["r_valley"]))
    if t == "finocyl":
        rc = 0.5 * geom["core_d"]
        phi = sd_circle(X, Y, 0.0, 0.0, rc)
        for k in range(geom["slots"]):
            th = 2 * math.pi * k / geom["slots"]
            phi = np.minimum(phi, sd_slot(X, Y, th, 0.0, geom["slot_r"], geom["slot_w"]))
        return phi
    raise ValueError("unknown geometry type " + t)


# --------------------------------------------------------------------------
# grid
# --------------------------------------------------------------------------

def make_grid(R, N):
    """N x N node grid covering the grain disc of radius R plus a 2-cell margin."""
    h = 2.0 * R / (N - 5)
    half = 0.5 * (N - 1) * h
    x = -half + h * np.arange(N)
    X, Y = np.meshgrid(x, x)  # X[iy, ix]
    return x, X, Y, h


# --------------------------------------------------------------------------
# fast marching (second-order upwind, Sethian 1999)
# --------------------------------------------------------------------------

def fmm(phi0: np.ndarray, h: float, band_cells: float = 3.0) -> np.ndarray:
    """Arrival-time field T (|grad T| = 1) marched outward from the port.

    Nodes with a 4-neighbour of opposite sign, plus propellant nodes within
    `band_cells` cells of the surface, form the frozen initial band and keep
    T = phi0 (the implicit function is a signed distance there; this narrow-band
    initialisation removes most of the first-order error the march otherwise
    picks up near sharp star tips and slot corners).  Port
    interior nodes keep phi0 but are never used in stencils.  Propellant nodes
    are accepted in increasing T with a second-order one-sided difference
    where two accepted upwind nodes exist, first order otherwise.
    Deterministic tie-break on (T, flat index) so the JS port reproduces it.
    """
    ny, nx = phi0.shape
    P = phi0.ravel().tolist()
    n = nx * ny
    T = [math.inf] * n
    state = [0] * n  # 0 far, 1 trial, 2 accepted, 3 port interior (unusable)
    neg = phi0 <= 0.0
    band = np.zeros_like(neg)
    band[:-1, :] |= neg[:-1, :] != neg[1:, :]
    band[1:, :] |= neg[1:, :] != neg[:-1, :]
    band[:, :-1] |= neg[:, :-1] != neg[:, 1:]
    band[:, 1:] |= neg[:, 1:] != neg[:, :-1]
    if band_cells > 0:
        band |= (phi0 > 0) & (phi0 < band_cells * h)
    bandf = band.ravel().tolist()
    negf = neg.ravel().tolist()
    for k in range(n):
        if bandf[k]:
            T[k] = P[k]
            state[k] = 2
        elif negf[k]:
            T[k] = P[k]
            state[k] = 3
    heap = []
    hh = h * h

    def axis_term(k, i, j, di, dj):
        """Best upwind (coef, centre) along one axis, or None."""
        best = None
        for s in (-1, 1):
            i1, j1 = i + s * di, j + s * dj
            if not (0 <= i1 < ny and 0 <= j1 < nx):
                continue
            k1 = i1 * nx + j1
            if state[k1] != 2:
                continue
            t1 = T[k1]
            term = (1.0, t1)
            i2, j2 = i1 + s * di, j1 + s * dj
            if 0 <= i2 < ny and 0 <= j2 < nx:
                k2 = i2 * nx + j2
                if state[k2] == 2 and T[k2] <= t1:
                    term = (2.25, (4.0 * t1 - T[k2]) / 3.0)
            if best is None or term[1] < best[1]:
                best = term
        return best

    def solve(k):
        i, j = divmod(k, nx)
        ax = axis_term(k, i, j, 0, 1)
        ay = axis_term(k, i, j, 1, 0)
        terms = [t for t in (ax, ay) if t is not None]
        terms.sort(key=lambda t: t[1])
        c1, m1 = terms[0]
        tnew = m1 + h / math.sqrt(c1)
        if len(terms) == 2 and tnew > terms[1][1]:
            c2, m2 = terms[1]
            a = c1 + c2
            b = -2.0 * (c1 * m1 + c2 * m2)
            c = c1 * m1 * m1 + c2 * m2 * m2 - hh
            disc = b * b - 4.0 * a * c
            if disc >= 0.0:
                tnew = (-b + math.sqrt(disc)) / (2.0 * a)
        return tnew

    def push_neighbours(k):
        i, j = divmod(k, nx)
        for di, dj in ((0, 1), (0, -1), (1, 0), (-1, 0)):
            i1, j1 = i + di, j + dj
            if 0 <= i1 < ny and 0 <= j1 < nx:
                k1 = i1 * nx + j1
                if state[k1] in (0, 1) and not negf[k1]:
                    tn = solve(k1)
                    if tn < T[k1]:
                        T[k1] = tn
                        state[k1] = 1
                        heapq.heappush(heap, (tn, k1))

    for k in range(n):
        if bandf[k] and not negf[k]:
            push_neighbours(k)
    while heap:
        tk, k = heapq.heappop(heap)
        if state[k] == 2 or tk > T[k]:
            continue
        state[k] = 2
        push_neighbours(k)
    return np.array(T).reshape(ny, nx)


# --------------------------------------------------------------------------
# marching squares
# --------------------------------------------------------------------------

def _cell_values(F):
    v00 = F[:-1, :-1]; v10 = F[:-1, 1:]; v11 = F[1:, 1:]; v01 = F[1:, :-1]
    return v00, v10, v11, v01


def contour_segments(F, x, level):
    """Marching-squares segments of {F = level}; returns (x1, y1, x2, y2) arrays.

    Segments are oriented so that the region F < level lies on their left.
    Saddle cells are resolved with the cell-centre average.
    """
    h = x[1] - x[0]
    v00, v10, v11, v01 = [v - level for v in _cell_values(F)]
    X0 = np.broadcast_to(x[:-1][None, :], v00.shape)
    Y0 = np.broadcast_to(x[:-1][:, None], v00.shape)

    def cross(va, vb):
        with np.errstate(divide="ignore", invalid="ignore"):
            return va / (va - vb)

    # crossing points on edges: 0 bottom (00-10), 1 right (10-11), 2 top (01-11), 3 left (00-01)
    ex = [X0 + h * cross(v00, v10), X0 + h, X0 + h * cross(v01, v11), X0]
    ey = [Y0, Y0 + h * cross(v10, v11), Y0 + h, Y0 + h * cross(v00, v01)]
    has = [(v00 < 0) != (v10 < 0), (v10 < 0) != (v11 < 0), (v01 < 0) != (v11 < 0), (v00 < 0) != (v01 < 0)]
    cnt = sum(hv.astype(int) for hv in has)
    corners = [(0.0, 0.0, v00), (h, 0.0, v10), (h, h, v11), (0.0, h, v01)]
    parts = []

    def add(a, b, m, iso=None):
        """Segment between crossings on edges a and b, oriented with F<level on the left.

        Orientation: sum over the cell corners of sign-weighted cross products
        (negative corners must lie left of the segment).  In saddle cells only the
        corner isolated by this segment is used.
        """
        x1, y1, x2, y2 = ex[a][m], ey[a][m], ex[b][m], ey[b][m]
        dx, dy = x2 - x1, y2 - y1
        score = np.zeros_like(x1)
        for ci, (ox, oy, vv) in enumerate(corners):
            if iso is not None and ci != iso:
                continue
            wgt = np.where(vv[m] < 0, 1.0, -1.0)
            cx = X0[m] + ox - x1; cy = Y0[m] + oy - y1
            score += wgt * (dx * cy - dy * cx)
        flip = score < 0
        parts.append((np.where(flip, x2, x1), np.where(flip, y2, y1), np.where(flip, x1, x2), np.where(flip, y1, y2)))

    # two-crossing cells: all 6 edge pairs
    for a in range(4):
        for b in range(a + 1, 4):
            m = (cnt == 2) & has[a] & has[b]
            if m.any():
                add(a, b, m)
    # saddle cells, resolved with the cell-centre average
    sad = cnt == 4
    if sad.any():
        vc = 0.25 * (v00 + v10 + v11 + v01)
        same = (vc < 0) == (v00 < 0)
        m1 = sad & same      # isolate corners 10 and 01: pairs (0,1) and (2,3)
        m2 = sad & ~same     # isolate corners 00 and 11: pairs (0,3) and (1,2)
        for (a, b), m, iso in (((0, 1), m1, 1), ((2, 3), m1, 3), ((0, 3), m2, 0), ((1, 2), m2, 2)):
            if m.any():
                add(a, b, m, iso)
    if not parts:
        z = np.zeros(0)
        return z, z, z, z
    return tuple(np.concatenate([p[k] for p in parts]) for k in range(4))


def clipped_length(x1, y1, x2, y2, R):
    """Total length of segments lying inside the disc |p| <= R."""
    dx, dy = x2 - x1, y2 - y1
    a = dx * dx + dy * dy
    b = 2.0 * (x1 * dx + y1 * dy)
    c = x1 * x1 + y1 * y1 - R * R
    disc = b * b - 4 * a * c
    ok = (disc > 0) & (a > 0)
    sq = np.sqrt(np.where(ok, disc, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        s1 = np.where(ok, (-b - sq) / (2 * a), 0.0)
        s2 = np.where(ok, (-b + sq) / (2 * a), 0.0)
    lo = np.clip(s1, 0.0, 1.0); hi = np.clip(s2, 0.0, 1.0)
    frac = np.where(ok, np.maximum(hi - lo, 0.0), 0.0)
    return float(np.sum(frac * np.sqrt(a)))


def enclosed_area(x1, y1, x2, y2):
    """Area of {F < level} from consistently oriented closed contour segments (shoelace)."""
    return float(0.5 * np.sum(x1 * y2 - x2 * y1))


def sample_on_circle(F, x, R, m=720):
    """Bilinear samples of F on the circle of radius R (used for case-contact web)."""
    h = x[1] - x[0]
    th = np.linspace(0, 2 * np.pi, m, endpoint=False)
    px, py = R * np.cos(th), R * np.sin(th)
    fx = (px - x[0]) / h; fy = (py - x[0]) / h
    j = np.floor(fx).astype(int); i = np.floor(fy).astype(int)
    tx = fx - j; ty = fy - i
    return (F[i, j] * (1 - tx) * (1 - ty) + F[i, j + 1] * tx * (1 - ty)
            + F[i + 1, j] * (1 - tx) * ty + F[i + 1, j + 1] * tx * ty)


# --------------------------------------------------------------------------
# geometry table
# --------------------------------------------------------------------------

class GrainTable:
    """Burning perimeter and port area versus web distance for one cross-section.

    Attributes: w (uniform grid, 0..w_max), P, A (port area, clipped to disc),
    w_contact (first contact of the front with the case), w_max (burnout),
    R, h, T (arrival-time field) and x (grid coordinates).
    """

    def __init__(self, geom: dict, R: float, N: int = 161, nw: int = 160, exact_field=False):
        self.geom, self.R, self.N, self.nw = geom, R, N, nw
        x, X, Y, h = make_grid(R, N)
        phi0 = port_sdf(geom, X, Y)
        self.x, self.h = x, h
        self.T = phi0 if exact_field else fmm(phi0, h)
        inside = np.hypot(X, Y) <= R
        self.w_max = float(np.max(self.T[inside]))
        self.w_contact = float(np.min(sample_on_circle(self.T, x, R)))
        self.A0 = enclosed_area(*contour_segments(self.T, x, 0.0))
        self.w = np.linspace(0.0, self.w_max, nw)
        # the last level coincides with the largest nodal value, where the contour degenerates
        # onto the clip circle; evaluate it a tenth of a level-spacing earlier (left limit).
        lev = self.w.copy()
        lev[-1] -= 0.1 * (self.w[1] - self.w[0])
        self.P_raw = np.array([clipped_length(*contour_segments(self.T, x, wl), R) for wl in lev])
        dw = np.diff(self.w)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (self.P_raw[1:] + self.P_raw[:-1]) * dw)])
        disc = math.pi * R * R
        self.A_end_raw = float(self.A0 + cum[-1])
        self.closure = self.A_end_raw / disc - 1.0
        # Mass-conserving correction: scale P so that A0 + int P dw fills the disc exactly.
        # The raw closure error (first order for sharp stars/slots, ~1e-6 for circles) is reported.
        self.p_scale = (disc - self.A0) / cum[-1] if cum[-1] > 0 else 1.0
        self.P = self.P_raw * self.p_scale
        cumc = np.concatenate([[0.0], np.cumsum(0.5 * (self.P[1:] + self.P[:-1]) * dw)])
        self.A = np.minimum(self.A0 + cumc, disc)

    def perimeter(self, w):
        return np.interp(w, self.w, self.P, right=0.0)

    def port_area(self, w):
        return np.interp(w, self.w, self.A, right=math.pi * self.R ** 2)

    def contour(self, level):
        """Contour segments at a web distance, clipped for display (NaN separated polylines)."""
        x1, y1, x2, y2 = contour_segments(self.T, self.x, level)
        keep = (np.hypot(x1, y1) <= self.R * 1.0005) & (np.hypot(x2, y2) <= self.R * 1.0005)
        return x1[keep], y1[keep], x2[keep], y2[keep]


# --------------------------------------------------------------------------
# analytic references
# --------------------------------------------------------------------------

def star_offset_perimeter(points, r_tip, r_valley, w):
    """Exact perimeter of the outward offset of a star polygon (valid until an edge vanishes).

    P(w) = P0 + w * (sum of turning angles at convex tips)
               - 2 w * sum tan(|turning angle| / 2) at reflex valleys.
    Returns (P(w), w_valid) where w_valid is where the first edge is consumed.
    """
    v = star_vertices(points, r_tip, r_valley)
    n = len(v)
    P0 = sum(math.dist(v[i], v[(i + 1) % n]) for i in range(n))
    conv = 0.0; refl = 0.0; wv = math.inf
    for i in range(n):
        a, b, c = v[i - 1], v[i], v[(i + 1) % n]
        e1 = b - a; e2 = c - b
        turn = math.atan2(e1[0] * e2[1] - e1[1] * e2[0], e1 @ e2)
        if turn > 0:
            conv += turn
        else:
            tt = math.tan(-turn / 2)
            refl += tt
            wv = min(wv, math.dist(a, b) / tt)
    return P0 + np.asarray(w) * (conv - 2 * refl), wv


def bates_burning_area(D, d, L, w, ends=2):
    """Closed-form BATES burning area for one segment (outer surface inhibited)."""
    w = np.asarray(w, dtype=float)
    dc = d + 2 * w
    Lb = L - ends * w
    alive = (dc < D) & (Lb > 0)
    Ab = math.pi * dc * Lb + ends * 0.25 * math.pi * (D * D - dc * dc)
    return np.where(alive, Ab, 0.0)
