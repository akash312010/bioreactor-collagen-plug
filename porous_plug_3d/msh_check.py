"""Independent check of a Fluent ASCII 3D tetra .msh: orientation (RH normal toward c0), zone areas,
cell volumes and equi-volume skewness (Fluent definition)."""
import sys, re, numpy as np
txt = open(sys.argv[1]).read()
hdr = re.compile(r"\((10|12|13|45) \(([^)]*)\)(\(\n)?")
nodes = None; faces = []; cellz = []; names = {}
pos = 0
for m in hdr.finditer(txt):
    sid = m.group(1); h = m.group(2).split()
    if sid == '45':
        names[int(h[0])] = (h[1], h[2]); continue
    if m.group(3) is None or int(h[0], 16) == 0:
        if sid == '12' and int(h[0], 16) != 0: cellz.append([int(x, 16) for x in h])
        continue
    end = txt.index("))", m.end())
    block = txt[m.end():end]
    if sid == '10': nodes = np.array(block.split(), float).reshape(-1, 3)
    if sid == '13':
        z = int(h[0], 16); a = np.array([int(x, 16) for x in block.split()]).reshape(-1, 5); faces.append((z, int(h[3], 16), a))
allf = np.vstack([a for _, _, a in faces]); nc = int(allf[:, 3:].max())
# reconstruct tets from faces: collect node sets per cell
from collections import defaultdict
cn = defaultdict(set)
for col in (3, 4):
    for r in allf[allf[:, col] > 0]:
        cn[r[col]].update(r[:3])
T = np.array([sorted(cn[c]) for c in range(1, nc+1)]) - 1
assert T.shape[1] == 4
P = nodes; cc = P[T].mean(1)
fn = allf[:, :3]-1; c0 = allf[:, 3]-1
nrm = np.cross(P[fn[:, 1]]-P[fn[:, 0]], P[fn[:, 2]]-P[fn[:, 0]]); fcen = P[fn].mean(1)
ok = np.einsum('ij,ij->i', nrm, cc[c0]-fcen) > 0
print("orientation (RH normal toward c0):", ok.mean())
for z, bt, a in faces:
    A = 0.5*np.linalg.norm(np.cross(P[a[:, 1]-1]-P[a[:, 0]-1], P[a[:, 2]-1]-P[a[:, 0]-1]), axis=1).sum()
    print(f"  zone {z} {names.get(z)} bc={bt} faces={len(a)} area={A*1e6:.4f} mm^2")
p = P[T]; a = p[:, 1]-p[:, 0]; b = p[:, 2]-p[:, 0]; c = p[:, 3]-p[:, 0]
V = np.abs(np.einsum('ij,ij->i', np.cross(a, b), c))/6
# circumradius
A = np.stack([a, b, c], 1); rhs = 0.5*np.stack([np.einsum('ij,ij->i', a, a), np.einsum('ij,ij->i', b, b), np.einsum('ij,ij->i', c, c)], 1)
cen = np.linalg.solve(A, rhs[..., None])[..., 0]; Rc = np.linalg.norm(cen, axis=1)
Vid = 8*np.sqrt(3)/27*Rc**3; skew = (Vid-V)/Vid
print(f"cells {nc}, total volume {V.sum()*1e6:.4f} mL, min volume {V.min():.3e} m^3")
print(f"equi-volume skewness: max {skew.max():.3f}, 99.9 pct {np.percentile(skew,99.9):.3f}, 99 pct {np.percentile(skew,99):.3f}, mean {skew.mean():.3f}, cells > 0.9: {(skew>0.9).sum()}")
