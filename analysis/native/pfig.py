"""Common tools for the pore-scale figures (Physics of Fluids style).

Reads ANSYS Fluent binary .cas/.dat files directly (fluentbin.py), builds a triangulation of the
polygonal cells, interpolates cell values to nodes, and computes the stream function from the
solver's own face mass fluxes. All figures: white background, vector PDF, text kept as text
(TrueType, editable in Illustrator). Colour fields default to vector triangles (FULL_VECTOR=1).
The common save function enforces Arial 12 pt and black text; photo/SHG panels remain raster.
"""
import os, sys, collections, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.tri as mtri
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fluentbin import read_cas, read_dat

RHO = 1012.0
FULL_VECTOR = os.environ.get("FULL_VECTOR", "1") == "1"
DPI_RASTER = 600

from figure_style import STYLE, WIDTH_IN, save_figure
plt.rcParams.update(STYLE)
WALL, INL, OUT, GREY = "#2b2b28", "#eb6834", "#2a78d6", "#8a8a85"
COL_W1, COL_W2 = 8.5/2.54, WIDTH_IN  # final single / double column width (in)


class Case:
    """One Fluent case + data. Coordinates in micrometres."""

    def __init__(self, cas, dat=None):
        P, F, names = read_cas(cas); self.D = read_dat(dat or cas[:-4]+".dat")
        self.P = P*1e6; self.F = F; self.names = names
        self.zid = {n: z for z, (t, n) in names.items()}
        self.fluid = [z for z, (t, n) in names.items() if t == "fluid"][0]
        self._faces(); self._tris(); self._cellgeom()

    # ---------- mesh ----------
    def _faces(self):
        n0, n1, c0, c1, fl, zz = [], [], [], [], [], []
        for z, (bc, f0, fa) in self.F.items():
            n0.append(fa[:, 0]-1); n1.append(fa[:, 1]-1); c0.append(fa[:, 2]-1); c1.append(fa[:, 3]-1)
            fl.append(self.D[(18, z)][1] if (18, z) in self.D else np.zeros(len(fa))); zz.append(np.full(len(fa), z))
        self.n0, self.n1, self.c0, self.c1 = map(np.concatenate, (n0, n1, c0, c1))
        self.flux = np.concatenate(fl); self.fzone = np.concatenate(zz)
        self.nc = int(max(self.c0.max(), self.c1.max())+1)

    def _tris(self):
        nxt = collections.defaultdict(dict)
        for a, b, p, q in zip(self.n0, self.n1, self.c0, self.c1):
            nxt[p][a] = b
            if q >= 0: nxt[q][b] = a
        tris, owner = [], []
        for c, d in nxt.items():
            s = next(iter(d)); ring = [s]
            while True:
                w = d[ring[-1]]
                if w == s: break
                ring.append(w)
            for k in range(1, len(ring)-1): tris.append((ring[0], ring[k], ring[k+1])); owner.append(c)
        self.tris = np.array(tris); self.tri_cell = np.array(owner)
        self.T = mtri.Triangulation(self.P[:, 0], self.P[:, 1], self.tris)

    def _cellgeom(self):
        x, y = self.P[self.tris, 0], self.P[self.tris, 1]
        a = 0.5*np.abs((x[:, 1]-x[:, 0])*(y[:, 2]-y[:, 0])-(x[:, 2]-x[:, 0])*(y[:, 1]-y[:, 0]))
        self.cell_area = np.bincount(self.tri_cell, a, self.nc)
        # node weights: area of attached triangles
        self._nw = np.zeros(len(self.P)); [np.add.at(self._nw, self.tris[:, k], a/3) for k in range(3)]
        self._ta = a

    # ---------- fields ----------
    def cell(self, fid):
        return self.D[(fid, self.fluid)][1]

    def to_nodes(self, cv):
        """Area-weighted average of cell values onto nodes."""
        v = np.zeros(len(self.P)); tv = cv[self.tri_cell]*self._ta/3
        for k in range(3): np.add.at(v, self.tris[:, k], tv)
        return v/np.maximum(self._nw, 1e-300)

    def speed(self):
        return np.hypot(self.cell(111), self.cell(112))

    def q_inlet(self):
        """Volumetric flow per unit depth through the inlet (m^2/s)."""
        return -self.D[(18, self.zid["inlet"])][1].sum()/RHO

    def zone_segments(self, name):
        a = self.F[self.zid[name]][2]
        return np.stack([self.P[a[:, 0]-1], self.P[a[:, 1]-1]], 1)

    def stream_function(self):
        """psi at nodes from the face mass fluxes: psi(n1) - psi(n0) = flux(face, c0 -> c1).
        Normalised so that psi = 0 on left_wall and psi = 1 on right_wall."""
        adj = collections.defaultdict(list)
        for a, b, f in zip(self.n0, self.n1, self.flux): adj[a].append((b, f)); adj[b].append((a, -f))
        start = self.F[self.zid["left_wall"]][2][0, 0]-1
        psi = np.full(len(self.P), np.nan); psi[start] = 0.0; dq = collections.deque([start])
        while dq:
            a = dq.popleft()
            for b, f in adj[a]:
                if np.isnan(psi[b]): psi[b] = psi[a]+f; dq.append(b)
        self.psi_mismatch = float(np.abs(psi[self.n1]-psi[self.n0]-self.flux).max()/np.abs(self.flux).max())
        rw = psi[self.F[self.zid["right_wall"]][2][:, :2]-1].mean()
        return psi/rw


def draw_field(ax, case, nodal, cmap, vmin, vmax, levels=None):
    """Smooth colour field (Gouraud on nodal values). Rasterised unless FULL_VECTOR."""
    if FULL_VECTOR:
        cv = np.zeros(case.nc); cnt = np.zeros(case.nc)
        tv = nodal[case.tris].mean(1); np.add.at(cv, case.tri_cell, tv); np.add.at(cnt, case.tri_cell, 1)
        m = ax.tripcolor(case.T, facecolors=(cv/cnt)[case.tri_cell], cmap=cmap, vmin=vmin, vmax=vmax,
                         edgecolors="face", linewidth=0.05)
    else:
        m = ax.tripcolor(case.T, nodal, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
    return m


def draw_walls(ax, case, lw=0.7):
    from matplotlib.collections import LineCollection
    for nm, c, w in (("left_wall", WALL, lw), ("right_wall", WALL, lw), ("inlet", INL, lw*1.6), ("outlet", OUT, lw*1.6)):
        ax.add_collection(LineCollection(case.zone_segments(nm), colors=c, linewidths=w, capstyle="round", zorder=5))


def save(fig, stem, outdir):
    save_figure(fig, stem, outdir)
