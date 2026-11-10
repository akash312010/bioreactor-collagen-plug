"""Convert a 2D ANSYS Fluent binary case/data pair to a VTK legacy file for VisIt / ParaView.

Cells are written as the original Fluent polygons (VTK_POLYGON), coordinates in micrometres.
Cell data:  p (Pa), u, v (m/s), speed (m/s), and every user scalar present (uds0 ...).
Point data: psi, the stream function from the solver's own face mass fluxes,
            normalised to 0 on left_wall and 1 on right_wall (streamlines = contours of psi).
Field data: q (m^2/s through the inlet), Lc (um), A (um^2), p_in (length-weighted inlet pressure, Pa).

Usage: python fluent2vtk.py CASE.cas DATA.dat OUT.vtk
"""
import sys, collections, numpy as np
from fluentbin import read_cas, read_dat
RHO = 1012.0

def convert(cas, dat, out):
    P, F, names = read_cas(cas); D = read_dat(dat); P = P*1e6
    zid = {n: z for z, (t, n) in names.items()}; fluid = [z for z, (t, n) in names.items() if t == "fluid"][0]
    n0, n1, c0, c1, fl = [], [], [], [], []
    for z, (bc, f0, fa) in F.items():
        n0.append(fa[:, 0]-1); n1.append(fa[:, 1]-1); c0.append(fa[:, 2]-1); c1.append(fa[:, 3]-1)
        fl.append(D[(18, z)][1] if (18, z) in D else np.zeros(len(fa)))
    n0, n1, c0, c1, fl = map(np.concatenate, (n0, n1, c0, c1, fl))
    nc = int(max(c0.max(), c1.max())+1)
    # polygons: c0 lies to the left of n0 -> n1, so walking n0 -> n1 keeps c0 counter-clockwise
    nxt = collections.defaultdict(dict)
    for a, b, p, q in zip(n0, n1, c0, c1):
        nxt[p][a] = b
        if q >= 0: nxt[q][b] = a
    rings = []
    for c in range(nc):
        d = nxt[c]; s = next(iter(d)); ring = [s]
        while True:
            w = d[ring[-1]]
            if w == s: break
            ring.append(w)
        rings.append(ring)
    # stream function at nodes
    adj = collections.defaultdict(list)
    for a, b, f in zip(n0, n1, fl): adj[a].append((b, f)); adj[b].append((a, -f))
    start = F[zid["left_wall"]][2][0, 0]-1
    psi = np.full(len(P), np.nan); psi[start] = 0.0; dq = collections.deque([start])
    while dq:
        a = dq.popleft()
        for b, f in adj[a]:
            if np.isnan(psi[b]): psi[b] = psi[a]+f; dq.append(b)
    psi /= psi[F[zid["right_wall"]][2][:, :2]-1].mean()
    mism = float(np.abs((psi[n1]-psi[n0])-fl/np.abs(fl).max()*0).max())
    # fields
    cell = {"p": D[(1, fluid)][1], "u": D[(111, fluid)][1], "v": D[(112, fluid)][1]}
    cell["speed"] = np.hypot(cell["u"], cell["v"])
    for k in sorted(k for k in D if k[1] == fluid and 700 <= k[0] < 750):
        cell[f"uds{k[0]-700}"] = D[k][1]
    seg = P[F[zid["inlet"]][2][:, :2]-1]; L = np.linalg.norm(seg[:, 1]-seg[:, 0], axis=1)
    pin = float(np.sum(D[(1, zid["inlet"])][1]*L)/L.sum()) if (1, zid["inlet"]) in D else float("nan")
    x, y = P[:, 0], P[:, 1]
    area = sum(0.5*abs(np.dot(x[r], np.roll(y[r], -1))-np.dot(y[r], np.roll(x[r], -1))) for r in rings)
    q = float(-D[(18, zid["inlet"])][1].sum()/RHO)
    with open(out, "w") as f:
        f.write("# vtk DataFile Version 3.0\nFluent 2D case converted by fluent2vtk.py (um)\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write("FIELD FieldData 4\n")
        for nm, val in (("q_m2s", q), ("Lc_um", float(x.max())), ("A_um2", area), ("p_in_Pa", pin)):
            f.write(f"{nm} 1 1 double\n{val:.10g}\n")
        f.write(f"POINTS {len(P)} double\n"); np.savetxt(f, np.c_[x, y, np.zeros(len(P))], fmt="%.9g")
        f.write(f"CELLS {nc} {sum(len(r)+1 for r in rings)}\n")
        for r in rings: f.write(f"{len(r)} " + " ".join(map(str, r)) + "\n")
        f.write(f"CELL_TYPES {nc}\n"); f.write("\n".join(["7"]*nc)+"\n")
        f.write(f"CELL_DATA {nc}\n")
        for nm, val in cell.items():
            f.write(f"SCALARS {nm} double 1\nLOOKUP_TABLE default\n"); np.savetxt(f, val, fmt="%.9g")
        f.write(f"VECTORS vel double\n"); np.savetxt(f, np.c_[cell["u"], cell["v"], np.zeros(nc)], fmt="%.9g")
        f.write(f"POINT_DATA {len(P)}\nSCALARS psi double 1\nLOOKUP_TABLE default\n"); np.savetxt(f, psi, fmt="%.9g")
    return dict(cells=nc, q=q, Lc=float(x.max()), A=area, p_in=pin, fields=list(cell))

if __name__ == "__main__":
    print(convert(*sys.argv[1:4]))
