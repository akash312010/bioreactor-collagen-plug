"""3D tetrahedral mesh of the bioreactor (inlet conduit, porous plug chamber, outlet conduit) written in
Fluent ASCII .msh format, like the 2D pore meshes. Units: metres in the file.
Face orientation follows Fluent's convention, verified on the Fluent-written 1M device case: the
right-hand-rule normal of each face's node order points toward cell c0 (boundary faces: c1 = 0).
Usage: python make_plug3d_mesh.py M1|M2"""
import sys, time, numpy as np, gmsh
level = sys.argv[1] if len(sys.argv) > 1 else 'M1'
f = {'M1': 1.0, 'M2': 0.5}[level]                       # M2 halves the sizes near the ports
ri, R, ro = 1.93, 16.664, 2.4995
y_in, y_roof, y_floor, y_out = 0.0, -31.68, -58.553, -60.33
h_far, h_port, h_cond = 0.8, 0.10*f, 0.30              # mm
t0 = time.time()
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("plug3d")
occ = gmsh.model.occ
cin  = occ.addCylinder(0, y_in,    0, 0, y_roof-y_in,    0, ri)
cch  = occ.addCylinder(0, y_roof,  0, 0, y_floor-y_roof, 0, R)
cout = occ.addCylinder(0, y_floor, 0, 0, y_out-y_floor,  0, ro)
occ.fragment([(3, cch)], [(3, cin), (3, cout)]); occ.synchronize()
vol = {}
for d, t in gmsh.model.getEntities(3):
    _, ymin, _, _, ymax, _ = gmsh.model.getBoundingBox(3, t)
    vol['conduit_in' if ymin > y_roof-1e-6 else 'conduit_out' if ymax < y_floor+1e-6 else 'plug'] = t
# port rim curves (circles at the two interfaces)
rims = []
for d, t in gmsh.model.getEntities(1):
    x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(1, t)
    rad = max(abs(x0), abs(x1), abs(z0), abs(z1))
    if (abs(y0-y_roof) < 1e-6 and abs(rad-ri) < 1e-3) or (abs(y0-y_floor) < 1e-6 and abs(rad-ro) < 1e-3):
        rims.append(t)
mf = gmsh.model.mesh.field
fd = mf.add("Distance"); mf.setNumbers(fd, "CurvesList", rims); mf.setNumber(fd, "Sampling", 400)
th = mf.add("Threshold"); mf.setNumber(th, "InField", fd)
mf.setNumber(th, "SizeMin", h_port); mf.setNumber(th, "SizeMax", h_far)
mf.setNumber(th, "DistMin", 0.2); mf.setNumber(th, "DistMax", 10.0)
cyl = mf.add("Box")   # conduits: cap the size
mf.setNumber(cyl, "VIn", h_cond); mf.setNumber(cyl, "VOut", h_far)
mf.setNumber(cyl, "XMin", -ro-0.01); mf.setNumber(cyl, "XMax", ro+0.01); mf.setNumber(cyl, "ZMin", -ro-0.01); mf.setNumber(cyl, "ZMax", ro+0.01)
mf.setNumber(cyl, "YMin", y_out-0.01); mf.setNumber(cyl, "YMax", y_in+0.01)
box2 = mf.add("Box")  # only outside the chamber (conduit regions in y)
mf.setNumber(box2, "VIn", h_far); mf.setNumber(box2, "VOut", h_far)
mn = mf.add("Min"); mf.setNumbers(mn, "FieldsList", [th]); 
# apply the conduit cap only in the conduit y-ranges
cap_in = mf.add("Box"); [mf.setNumber(cap_in, k, v) for k, v in dict(VIn=h_cond, VOut=1e9, XMin=-ri-0.01, XMax=ri+0.01, ZMin=-ri-0.01, ZMax=ri+0.01, YMin=y_roof, YMax=y_in+0.01).items()]
cap_out = mf.add("Box"); [mf.setNumber(cap_out, k, v) for k, v in dict(VIn=h_cond, VOut=1e9, XMin=-ro-0.01, XMax=ro+0.01, ZMin=-ro-0.01, ZMax=ro+0.01, YMin=y_out-0.01, YMax=y_floor).items()]
# port refinement: full grading (to 10 mm) inside the plug; in the conduits it relaxes to the conduit size within 3 mm
rp = mf.add("Restrict"); mf.setNumber(rp, "InField", th); mf.setNumbers(rp, "VolumesList", [vol['plug']])
thc = mf.add("Threshold"); mf.setNumber(thc, "InField", fd)
mf.setNumber(thc, "SizeMin", h_port); mf.setNumber(thc, "SizeMax", h_cond)
mf.setNumber(thc, "DistMin", 0.2); mf.setNumber(thc, "DistMax", 3.0)
rc = mf.add("Restrict"); mf.setNumber(rc, "InField", thc); mf.setNumbers(rc, "VolumesList", [vol['conduit_in'], vol['conduit_out']])
mn = mf.add("Min"); mf.setNumbers(mn, "FieldsList", [rp, rc, cap_in, cap_out]); mf.setAsBackgroundMesh(mn)
for k, v in {"Mesh.MeshSizeExtendFromBoundary": 0, "Mesh.MeshSizeFromPoints": 0, "Mesh.MeshSizeFromCurvature": 0,
             "Mesh.Algorithm": 6, "Mesh.Algorithm3D": 10, "Mesh.Optimize": 1, "Mesh.OptimizeNetgen": 1,
             "General.NumThreads": 4}.items():
    gmsh.option.setNumber(k, v)
gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.OptimizeThreshold", 0.45)
gmsh.model.mesh.optimize("Netgen", niter=3)
gmsh.model.mesh.optimize("", niter=3)
print(f"meshed in {time.time()-t0:.0f} s")
tags, xyz, _ = gmsh.model.mesh.getNodes(); xyz = xyz.reshape(-1, 3)
idx = np.zeros(int(tags.max())+1, np.int64); idx[tags.astype(np.int64)] = np.arange(len(tags))
P = xyz * 1e-3                                         # mm -> m
cells, czone = [], []
for zname in ('conduit_in', 'plug', 'conduit_out'):
    et, _, nt = gmsh.model.mesh.getElements(3, vol[zname])
    k = list(et).index(4); T = idx[nt[k].astype(np.int64)].reshape(-1, 4)
    cells.append(T); czone.append(np.full(len(T), ('conduit_in', 'plug', 'conduit_out').index(zname)))
gmsh.finalize()
T = np.vstack(cells); cz = np.concatenate(czone)
# positive orientation
v = np.einsum('ij,ij->i', np.cross(P[T[:, 1]]-P[T[:, 0]], P[T[:, 2]]-P[T[:, 0]]), P[T[:, 3]]-P[T[:, 0]])
neg = v < 0; T[neg, 1], T[neg, 2] = T[neg, 2].copy(), T[neg, 1].copy()
vol6 = np.abs(v)
# outward-oriented faces of each tet (for a positive tet n0..n3)
loc = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])
F = T[:, loc].reshape(-1, 3); owner = np.repeat(np.arange(len(T)), 4)
key = np.sort(F, axis=1)
_, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
inv = inv.ravel()
order = np.argsort(inv, kind='stable'); inv_s = inv[order]
first = np.r_[True, inv_s[1:] != inv_s[:-1]]
fa = order[first]                                      # first occurrence of each unique face
nU = len(fa)
other = np.full(nU, -1, np.int64)
dup = order[~first]; other[inv[dup]] = owner[dup]
A = owner[fa]; nodesA = F[fa]                          # nodes ordered with normal pointing OUT of A
interior = other >= 0
# Fluent: RH normal points toward c0. Interior: normal out of A into B -> c0 = B, c1 = A.
# Boundary: reverse order so the normal points into A -> c0 = A, c1 = 0.
fn = nodesA.copy(); fn[~interior] = fn[~interior][:, ::-1]
c0 = np.where(interior, other, A); c1 = np.where(interior, A, -1)
fcen = P[fn].mean(1)
# classify faces
zone = np.empty(nU, object)
same = interior & (cz[np.maximum(c0, 0)] == cz[np.maximum(c1, 0)])
zone[interior & same] = np.array(['interior-' + ('conduit_in', 'plug', 'conduit_out')[z] for z in cz[c0[interior & same]]], object)
port = interior & ~same
zone[port & (fcen[:, 1] > (y_roof+y_floor)/2*1e-3)] = 'port_in'
zone[port & (fcen[:, 1] < (y_roof+y_floor)/2*1e-3)] = 'port_out'
b = ~interior
tol = 1e-9
zone[b & (np.abs(fcen[:, 1]-y_in*1e-3) < tol)] = 'inlet'
zone[b & (np.abs(fcen[:, 1]-y_out*1e-3) < tol)] = 'outlet'
rest = b & (zone == None)
zone[rest] = np.array(['wall_' + ('conduit_in', 'plug', 'conduit_out')[z] for z in cz[c0[rest]]], object)
assert not (zone == None).any()
# checks
assert (vol6 > 0).all()
nrm = np.cross(P[fn[:, 1]]-P[fn[:, 0]], P[fn[:, 2]]-P[fn[:, 0]])
ccen = P[T].mean(1)
chk = np.einsum('ij,ij->i', nrm, ccen[c0]-fcen)
assert (chk > 0).all(), "orientation check failed"
pin = zone == 'port_in'; pout = zone == 'port_out'
area = 0.5*np.linalg.norm(nrm, axis=1)
print(f"{level}: nodes {len(P)}, cells {len(T)} (conduit_in {np.sum(cz==0)}, plug {np.sum(cz==1)}, conduit_out {np.sum(cz==2)}), faces {nU}")
print(f"port_in area {area[pin].sum()*1e6:.4f} mm^2 (exact {np.pi*ri**2:.4f}), port_out {area[pout].sum()*1e6:.4f} (exact {np.pi*ro**2:.4f}); "
      f"port_in faces {pin.sum()}, port_out faces {pout.sum()}")
print(f"plug volume {vol6[cz==1].sum()/6*1e6:.4f} mL (exact {np.pi*R**2*(y_roof-y_floor)/1000:.4f})")
# tet quality: radius ratio
def quality(T):
    p = P[T]; a = p[:, 1]-p[:, 0]; b_ = p[:, 2]-p[:, 0]; c = p[:, 3]-p[:, 0]
    V = np.abs(np.einsum('ij,ij->i', np.cross(a, b_), c))/6
    edges = [p[:, i]-p[:, j] for i, j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))]
    Lrms = np.sqrt(np.mean([np.einsum('ij,ij->i', e, e) for e in edges], axis=0))
    return 6*np.sqrt(2)*V/Lrms**3          # 1 for regular tet
q = quality(T); print(f"tet quality (V/Lrms^3 normalised): min {q.min():.3f}, 1st pct {np.percentile(q,1):.3f}, median {np.median(q):.3f}")
# ---- write Fluent ASCII ----
zone_order = ['interior-conduit_in', 'interior-plug', 'interior-conduit_out', 'port_in', 'port_out',
              'inlet', 'outlet', 'wall_conduit_in', 'wall_plug', 'wall_conduit_out']
bc = {'inlet': ('mass-flow-inlet', 0x14), 'outlet': ('pressure-outlet', 0x5)}
for z in zone_order:
    bc.setdefault(z, ('interior', 0x2) if z.startswith('interior') or z.startswith('port') else ('wall', 0x3))
# renumber cells by zone (already ordered by construction)
out = f"bioreactor_plug3d_{level}.msh"
with open(out, 'w') as fo:
    fo.write(f'(0 "Bioreactor 3D porous-plug mesh {level}: gmsh tetra, Fluent ASCII, units m")\n(2 3)\n')
    fo.write(f"(10 (0 1 {len(P):x} 0 3))\n(12 (0 1 {len(T):x} 0 0))\n(13 (0 1 {nU:x} 0 0))\n")
    fo.write(f"(10 (1 1 {len(P):x} 1 3)(\n")
    np.savetxt(fo, P, fmt='%.9e'); fo.write("))\n")
    start = 1; cid = {}
    for zi, zn in enumerate(('conduit_in', 'plug', 'conduit_out')):
        n = int(np.sum(cz == zi)); cid[zn] = 2+zi
        fo.write(f"(12 ({cid[zn]:x} {start:x} {start+n-1:x} 1 2))\n"); start += n
    start = 1; fid = {}
    for k, zn in enumerate(zone_order):
        sel = np.where(zone == zn)[0]; fid[zn] = 10+k
        if len(sel) == 0: continue
        fo.write(f"(13 ({fid[zn]:x} {start:x} {start+len(sel)-1:x} {bc[zn][1]:x} 3)(\n")
        rows = np.column_stack([fn[sel]+1, c0[sel]+1, np.where(c1[sel] >= 0, c1[sel]+1, 0)])
        np.savetxt(fo, rows, fmt='%x'); fo.write("))\n"); start += len(sel)
    for zn in ('conduit_in', 'plug', 'conduit_out'): fo.write(f"(45 ({cid[zn]} fluid {zn})())\n")
    for zn in zone_order: fo.write(f"(45 ({fid[zn]} {bc[zn][0]} {zn})())\n")
print("wrote", out, f"in {time.time()-t0:.0f} s total")
