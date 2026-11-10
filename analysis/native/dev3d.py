"""Read the 3D tetra device mesh (Fluent binary .cas) and cut it with a plane through the axis."""
import re, numpy as np
HDR = re.compile(rb"\((3010|3013) \(([0-9a-fA-F ]+)\)\s*\(")
def read3d(fn):
    b = open(fn, "rb").read(); nodes = None; faces = {}
    for m in HDR.finditer(b):
        sid = int(m.group(1)); h = [int(x, 16) for x in m.group(2).split()]; st = m.end()
        if sid == 3010:
            zone, f, l, typ, nd = h; nodes = np.frombuffer(b, "<f8", (l-f+1)*nd, st).reshape(-1, nd)
        else:
            zone, f, l, bc, typ = h; assert typ == 3, typ
            faces[zone] = np.frombuffer(b, "<i4", (l-f+1)*5, st).reshape(-1, 5).astype(np.int64)
    names = {int(z): (t.decode(), n.decode()) for z, t, n in re.findall(rb"\(39 \((\d+) ([a-z-]+) ([A-Za-z0-9_\-:.]+)", b)}
    return nodes, faces, names
def cut(P, tri, axis_n, c=0.0):
    """Segments where triangles (n,3 node idx) cross the plane x_axis_n = c."""
    d = P[tri][:, :, axis_n]-c; s = np.sign(d)
    cross = (s.max(1) > 0) & (s.min(1) < 0); T = tri[cross]; dd = d[cross]
    segs = []
    for (i, j) in ((0, 1), (1, 2), (2, 0)):
        pass
    pts = []
    for (i, j) in ((0, 1), (1, 2), (2, 0)):
        m = dd[:, i]*dd[:, j] < 0
        t = dd[:, i]/(dd[:, i]-dd[:, j]); X = P[T[:, i]]+(P[T[:, j]]-P[T[:, i]])*t[:, None]
        pts.append(np.where(m[:, None], X, np.nan))
    pts = np.stack(pts, 1)               # (n, 3 edges, 3)
    out = []
    for k in range(len(T)):
        q = pts[k][~np.isnan(pts[k][:, 0])]
        if len(q) == 2: out.append(q)
    return np.array(out)
