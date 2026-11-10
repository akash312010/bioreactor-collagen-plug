"""Test of the curvilinear second-order lubrication theory (curvi.py) against CFD, for several reference
centrelines (Gaussian smoothing SIG_C = 1, 2, 3, 4 um) and comparison intervals.
Local theory (relative to the leading term -3 mu q /(2 h^3)), walls at n = c -+ h about the reference curve:
  1 - c k + c'^2 + 7/5 h'^2 - h h''/10 + 2/15 h^2 k^2 - c c'' + 6 c c' h'/h - 6 c^2 h'^2/h^2 + 3/2 c^2 h''/h
Pressure drop over [s1, s2] integrated directly (no integration by parts, so no end term is dropped);
CFD pressure interpolated onto the reference curve n = 0, which is where the theory gives dp/ds."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from fluentbin import read_cas, read_dat
from smooth import load, smooth_wall
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree
from scipy.interpolate import LinearNDInterpolator
from shapely.geometry import LineString
U = os.environ.get("CFD_RUNS", "CFD_runs_new/")  # Fluent case/data NOT in this package
CAD = os.path.join(HERE, "..", "..", "data", "geometry") + os.sep
mu, rho = 9.64e-4, 1012.0
runs = {"G1": ("p28", U+"D1/new/New Run_D1/collagen_p28_fine", 3.0), "G2": ("p31", U+"D2/new/New_D2/collagen_p31_fine", 3.0),
        "G3": ("p35", U+"B1_v2_p35_fine/New60k_fine_mesh/New60k_fine_mesh/fine_mesh", 3.0),
        "G4": ("p40", U+"D3/new/New_D3/D3_new", 3.0), "G5": ("p45", U+"D4/new/New_D4/collagen_p45_fine", 1.0)}
SIG_D = 1.0                      # smoothing of the derivatives (um), as before

def frame(a, b, sigc):
    LA, LB = LineString(a), LineString(b)
    xs = np.linspace(0, min(a[:, 0].max(), b[:, 0].max()), 3001)
    ya = np.interp(xs, *a[np.argsort(a[:, 0])].T); yb = np.interp(xs, *b[np.argsort(b[:, 0])].T)
    dx = xs[1]-xs[0]; yc = gaussian_filter1d(0.5*(ya+yb), sigc/dx, mode='nearest')
    C = np.c_[xs, yc]; t = np.gradient(C, axis=0); t /= np.linalg.norm(t, axis=1)[:, None]; nrm = np.c_[-t[:, 1], t[:, 0]]
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))]
    th = np.unwrap(np.arctan2(t[:, 1], t[:, 0])); kap = np.gradient(th, s)
    na = np.full(len(xs), np.nan); nb = np.full(len(xs), np.nan)
    for i in range(len(xs)):
        seg = LineString([C[i]-30*nrm[i], C[i]+30*nrm[i]]); vals = []
        for Lw in (LA, LB):
            I = seg.intersection(Lw); pts = [I] if I.geom_type == 'Point' else list(getattr(I, 'geoms', []))
            if not pts: break
            p = min(pts, key=lambda q: np.hypot(q.x-C[i, 0], q.y-C[i, 1])); vals.append((p.x-C[i, 0])*nrm[i, 0]+(p.y-C[i, 1])*nrm[i, 1])
        if len(vals) == 2: na[i], nb[i] = max(vals), min(vals)
    ok = ~np.isnan(na); na = np.interp(s, s[ok], na[ok]); nb = np.interp(s, s[ok], nb[ok])
    return C, s, kap, na, nb

def run(lab):
    pk, base, blend = runs[lab]
    A, B = load(CAD+f"collagen_channel_interior_{pk}.csv"); a = smooth_wall(A, 0.75, blend=blend); b = smooth_wall(B, 0.75, blend=blend)
    P, F, names = read_cas(base+".cas"); D = read_dat(base+".dat"); zid = {n: z for z, (t, n) in names.items()}
    zf = [z for z, (t, n) in names.items() if t == 'fluid'][0]; p = D[(1, zf)][1]; nc = len(p)
    Ac = np.zeros(nc); Cx = np.zeros(nc); Cy = np.zeros(nc)
    for z, (bc, f0, fa) in F.items():
        n0, n1, c0, c1 = fa[:, 0]-1, fa[:, 1]-1, fa[:, 2]-1, fa[:, 3]-1; x0, y0, x1, y1 = P[n0, 0], P[n0, 1], P[n1, 0], P[n1, 1]; cr = x0*y1-x1*y0
        for cc, sg in ((c0, 1.), (c1, -1.)):
            m = cc >= 0; np.add.at(Ac, cc[m], sg*cr[m]/2); np.add.at(Cx, cc[m], sg*(x0[m]+x1[m])*cr[m]/6); np.add.at(Cy, cc[m], sg*(y0[m]+y1[m])*cr[m]/6)
    cen = np.c_[Cx/Ac, Cy/Ac]*1e6; q = -D[(18, zid['inlet'])][1].sum()/rho
    out = {}
    for sigc in (1.0, 2.0, 3.0, 4.0):
        C, s, kap, na, nb = frame(a, b, sigc); ds = s[1]-s[0]
        g = lambda f: gaussian_filter1d(f, SIG_D/ds, mode='nearest')
        h = 0.5*(na-nb); cc = 0.5*(na+nb); hs, cs, ks = g(h), g(cc), g(kap)
        h1 = np.gradient(hs, s); h2 = np.gradient(h1, s); c1 = np.gradient(cs, s); c2 = np.gradient(c1, s)
        corr = (-cs*ks + c1**2 + 1.4*h1**2 - 0.1*hs*h2 + (2/15)*hs**2*ks**2 - cs*c2 + 6*cs*c1*h1/hs - 6*cs**2*h1**2/hs**2 + 1.5*cs**2*h2/hs)
        corr_cart = 1.4*h1**2 - 0.1*hs*h2      # without curvature and offset terms
        lead = 1.5*mu*q/(h*1e-6)**3          # Pa/m, h in um
        # CFD pressure ON the reference curve (n = 0), the quantity the theory predicts
        pS = LinearNDInterpolator(cen, p)(C)
        okp = ~np.isnan(pS); pS = np.interp(s, s[okp], pS[okp])
        res = {}
        for lo, hi in ((4, 8), (6, 12), (8, 16), (10, 20)):
            m = (s > lo) & (s < s[-1]-hi); i0, i1 = np.where(m)[0][[0, -1]]
            dp_cfd = pS[i1]-pS[i0]
            T = lambda f: np.trapezoid(f[m], s[m])*1e-6          # non-uniform s spacing
            dp_lub = T(lead); dp_th = T(lead*(1+corr)); dp_th0 = T(lead*(1+corr_cart))
            res["%d_%d" % (lo, hi)] = dict(cfd_over_lub=dp_cfd/dp_lub, theory_over_lub=dp_th/dp_lub, theory_over_cfd=dp_th/dp_cfd,
                                           no_curv_over_cfd=dp_th0/dp_cfd,
                raw_signed_metric_min=float(np.min((1-cc*kap-np.abs(kap)*h)[m])),
                filtered_signed_metric_min=float(np.min((1-cs*ks-np.abs(ks)*hs)[m])),
                raw_metric_nonpositive_fraction=float(np.sum(np.gradient(s)[m]*((1-cc*kap-np.abs(kap)*h)[m]<=0))/np.sum(np.gradient(s)[m])))
        m = (s > 6) & (s < s[-1]-12); dsi = np.gradient(s); wt = (h**-3*dsi)[m]/np.sum((h**-3*dsi)[m])
        res["terms_6_12"] = dict(hk2=float(np.sum(wt*(2/15)*hs[m]**2*ks[m]**2)), ck=float(np.sum(wt*(-cs[m]*ks[m]))), c1sq=float(np.sum(wt*c1[m]**2)),
                                 h1sq=float(np.sum(wt*1.4*h1[m]**2)), hh2=float(np.sum(wt*(-0.1*hs[m]*h2[m]))), max_abs_c_over_h=float(np.max(np.abs(cs[m]/hs[m]))))
        out["sig%g" % sigc] = res
        r = res["6_12"]; print(lab, "sigC", sigc, "| CFD/lub %.3f theory/lub %.3f theory/CFD %.3f (no curvature/offset terms %.3f)" % (r["cfd_over_lub"], r["theory_over_lub"], r["theory_over_cfd"], r["no_curv_over_cfd"]), "| terms", {k: round(v, 4) for k, v in res["terms_6_12"].items()}, flush=True)
    return out
if __name__ == '__main__':
    allres = {lab: run(lab) for lab in runs}
    json.dump(allres, open("curvi_test.json", "w"), indent=1)
