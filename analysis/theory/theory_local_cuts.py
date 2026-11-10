"""Local test of lubrication theory: CFD pressure gradient along the centreline vs 12 mu q / w^3,
and its relation to wall slope, width gradient and centreline curvature."""
import sys, json, numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from fluentbin import read_cas, read_dat
from smooth import load, smooth_wall
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree
from shapely.geometry import LineString, Point
U = os.environ.get("CFD_RUNS", "CFD_runs_new/")  # Fluent case/data NOT in this repository
CAD = os.path.join(HERE, "..", "..", "data", "geometry") + os.sep
mu, rho = 9.64e-4, 1012.0
runs = {"p28": (U+"D1/new/New Run_D1/collagen_p28_fine", 3.0), "p31": (U+"D2/new/New_D2/collagen_p31_fine", 3.0),
        "p35": (U+"B1_v2_p35_fine/New60k_fine_mesh/New60k_fine_mesh/fine_mesh", 3.0),
        "p40": (U+"D3/new/New_D3/D3_new", 3.0), "p45": (U+"D4/new/New_D4/collagen_p45_fine", 1.0)}
SIG_C = 2.0      # centreline smoothing (um), as in the paper
SIG_L = 1.0      # smoothing of local gradients (um)

def tangent_angle(line, pt, h=0.05):
    d = line.project(pt); p1 = line.interpolate(max(d-h, 0)); p2 = line.interpolate(min(d+h, line.length))
    return np.arctan2(p2.y-p1.y, p2.x-p1.x)

def geometry(a, b):
    LA, LB = LineString(a), LineString(b)
    xs = np.linspace(0, min(a[:, 0].max(), b[:, 0].max()), 3001)
    ya = np.interp(xs, *a[np.argsort(a[:, 0])].T); yb = np.interp(xs, *b[np.argsort(b[:, 0])].T)
    dx = xs[1]-xs[0]; yc = gaussian_filter1d(0.5*(ya+yb), SIG_C/dx, mode='nearest')
    C = np.c_[xs, yc]; t = np.gradient(C, axis=0); t /= np.linalg.norm(t, axis=1)[:, None]; nrm = np.c_[-t[:, 1], t[:, 0]]
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))]
    th_c = np.unwrap(np.arctan2(t[:, 1], t[:, 0])); kappa = np.gradient(th_c, s)
    w = np.full(len(xs), np.nan); ta = np.full(len(xs), np.nan); tb = np.full(len(xs), np.nan)
    for i in range(len(xs)):
        seg = LineString([C[i]-30*nrm[i], C[i]+30*nrm[i]]); hit = []
        for Lw in (LA, LB):
            I = seg.intersection(Lw); pts = [I] if I.geom_type == 'Point' else list(getattr(I, 'geoms', []))
            if not pts: break
            p = min(pts, key=lambda q: np.hypot(q.x-C[i, 0], q.y-C[i, 1])); hit.append((p, Lw))
        if len(hit) == 2:
            w[i] = sum(np.hypot(p.x-C[i, 0], p.y-C[i, 1]) for p, _ in hit)
            # wall angle relative to the centreline direction
            ta[i] = np.tan(np.angle(np.exp(1j*(tangent_angle(hit[0][1], hit[0][0])-th_c[i]))))
            tb[i] = np.tan(np.angle(np.exp(1j*(tangent_angle(hit[1][1], hit[1][0])-th_c[i]))))
    ok = ~np.isnan(w)
    for arr in (w, ta, tb): arr[~ok] = np.interp(s[~ok], s[ok], arr[ok])
    ta = np.clip(ta, -3, 3); tb = np.clip(tb, -3, 3)       # guard against near-normal walls at pinches
    return dict(C=C, s=s, w=w, kappa=kappa, ta=ta, tb=tb)

res = {}; local = {}
for lab, (base, blend) in runs.items():
    A, B = load(CAD+f"collagen_channel_interior_{lab}.csv"); a = smooth_wall(A, 0.75, blend=blend); b = smooth_wall(B, 0.75, blend=blend)
    G = geometry(a, b); s, w = G["s"], G["w"]; ds = s[1]-s[0]
    P, F, names = read_cas(base+".cas"); D = read_dat(base+".dat"); zid = {n: z for z, (t, n) in names.items()}
    zf = [z for z, (t, n) in names.items() if t == 'fluid'][0]; p = D[(1, zf)][1]; nc = len(p)
    Ac = np.zeros(nc); Cx = np.zeros(nc); Cy = np.zeros(nc)
    for z, (bc, f0, fa) in F.items():
        n0, n1, c0, c1 = fa[:, 0]-1, fa[:, 1]-1, fa[:, 2]-1, fa[:, 3]-1; x0, y0, x1, y1 = P[n0, 0], P[n0, 1], P[n1, 0], P[n1, 1]; cr = x0*y1-x1*y0
        for c, sg in ((c0, 1.), (c1, -1.)):
            m = c >= 0; np.add.at(Ac, c[m], sg*cr[m]/2); np.add.at(Cx, c[m], sg*(x0[m]+x1[m])*cr[m]/6); np.add.at(Cy, c[m], sg*(y0[m]+y1[m])*cr[m]/6)
    cen = np.c_[Cx/Ac, Cy/Ac]*1e6
    q = -D[(18, zid['inlet'])][1].sum()/rho
    # section-averaged pressure along s: assign each cell to the nearest centreline point
    _, idx = cKDTree(G["C"]).query(cen); sc = s[idx]
    bins = np.arange(0, s[-1]+0.25, 0.25); k = np.digitize(sc, bins)-1
    ps = np.array([np.sum(p[k == j]*Ac[k == j])/max(np.sum(Ac[k == j]), 1e-300) if np.any(k == j) else np.nan for j in range(len(bins)-1)])
    sb = 0.5*(bins[1:]+bins[:-1]); okb = ~np.isnan(ps)
    pS = np.interp(s, sb[okb], ps[okb])
    gc = np.gradient(gaussian_filter1d(pS, SIG_L/ds, mode='nearest'), s)*1e6          # Pa/m
    gl = gaussian_filter1d(12*mu*q/(w*1e-6)**3, SIG_L/ds, mode='nearest')             # Pa/m
    wp = np.gradient(w, s); e_w = gaussian_filter1d(wp**2, SIG_L/ds, mode='nearest')
    e_wall = gaussian_filter1d(G["ta"]**2+G["tb"]**2, SIG_L/ds, mode='nearest')
    e_k = gaussian_filter1d((w*G["kappa"])**2, SIG_L/ds, mode='nearest')
    inner = (s > 6) & (s < s[-1]-12)                                                # away from both ends
    cuts = {}
    for lo, hi in ((6, 12), (4, 8), (8, 16), (10, 20)):
        mm = (s > lo) & (s < s[-1]-hi); cuts["%d_%d" % (lo, hi)] = float(np.sum(gc[mm]*ds)/np.sum(gl[mm]*ds))
    wt = (w**-3)/np.sum(w[inner]**-3)                                                  # lubrication weight
    C_int = np.sum(gc[inner]*ds)/np.sum(gl[inner]*ds)
    res[lab] = dict(C_cuts=cuts, C_inner=float(C_int), w_eff=float((s[-1]/np.sum(w**-3*ds))**(1/3)),
                    eps_w=float(np.sum(wt[inner]*e_w[inner])), eps_wall=float(np.sum(wt[inner]*e_wall[inner])), eps_k=float(np.sum(wt[inner]*e_k[inner])),
                    mean_wall_slope2=float(np.mean(e_wall[inner])), mean_wk2=float(np.mean(e_k[inner])))
    local[lab] = dict(s=s[inner].tolist(), r=(gc[inner]/gl[inner]).tolist(), e_w=e_w[inner].tolist(), e_wall=e_wall[inner].tolist(), e_k=e_k[inner].tolist(), w=w[inner].tolist())
    print(lab, {k2: round(v, 4) for k2, v in cuts.items()}, flush=True)
json.dump(dict(summary=res), open(os.path.join(HERE, "theory_local_cuts.json"), "w"))
