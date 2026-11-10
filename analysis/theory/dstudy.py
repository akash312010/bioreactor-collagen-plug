import sys, json, numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from fluentbin import read_cas, read_dat
from smooth import load, smooth_wall
from scipy.ndimage import gaussian_filter1d
from shapely.geometry import LineString
U = os.environ.get("CFD_RUNS", "CFD_runs_new/")  # Fluent case/data NOT in this repository
CAD = os.path.join(HERE, "..", "..", "data", "geometry") + os.sep
mu, rho, phi = 9.64e-4, 1012.0, 0.27
runs = {"p28": (U+"D1/new/New Run_D1/collagen_p28_fine", 3.0), "p31": (U+"D2/new/New_D2/collagen_p31_fine", 3.0),
        "p35": (U+"B1_v2_p35_fine/New60k_fine_mesh/New60k_fine_mesh/fine_mesh", 3.0),
        "p40": (U+"D3/new/New_D3/D3_new", 3.0), "p45": (U+"D4/new/New_D4/collagen_p45_fine", 1.0)}

def width_profile(a, b, sig_um):
    LA, LB = LineString(a), LineString(b)
    xs = np.linspace(0, min(a[:, 0].max(), b[:, 0].max()), 4001)
    ya = np.interp(xs, *a[np.argsort(a[:, 0])].T); yb = np.interp(xs, *b[np.argsort(b[:, 0])].T)
    dx = xs[1]-xs[0]; yc = gaussian_filter1d(0.5*(ya+yb), sig_um/dx, mode='nearest')
    C = np.c_[xs, yc]; t = np.gradient(C, axis=0); t /= np.linalg.norm(t, axis=1)[:, None]; nrm = np.c_[-t[:, 1], t[:, 0]]
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))]; w = np.full(len(xs), np.nan)
    for i in range(len(xs)):
        seg = LineString([C[i]-30*nrm[i], C[i]+30*nrm[i]]); d = []
        for Lw in (LA, LB):
            I = seg.intersection(Lw); pts = [I] if I.geom_type == 'Point' else list(getattr(I, 'geoms', []))
            if not pts: break
            d.append(min(np.hypot(p.x-C[i, 0], p.y-C[i, 1]) for p in pts))
        if len(d) == 2: w[i] = d[0]+d[1]
    ok = ~np.isnan(w); w = np.interp(s, s[ok], w[ok])
    I = np.r_[0, np.cumsum(0.5*(w[1:]**-3+w[:-1]**-3)*np.diff(s))]
    return xs, s, w, I

res = {}
for lab, (base, blend) in runs.items():
    A, B = load(CAD+f"collagen_channel_interior_{lab}.csv"); a = smooth_wall(A, 0.75, blend=blend); b = smooth_wall(B, 0.75, blend=blend)
    P, F, names = read_cas(base+".cas"); D = read_dat(base+".dat"); zid = {n: z for z, (t, n) in names.items()}
    zf = [z for z, (t, n) in names.items() if t == 'fluid'][0]; p = D[(1, zf)][1]; nc = len(p)
    Ac = np.zeros(nc); Cx = np.zeros(nc)
    for z, (bc, f0, fa) in F.items():
        n0, n1, c0, c1 = fa[:, 0]-1, fa[:, 1]-1, fa[:, 2]-1, fa[:, 3]-1; x0, y0, x1, y1 = P[n0, 0], P[n0, 1], P[n1, 0], P[n1, 1]; cr = x0*y1-x1*y0
        for c, sgn in ((c0, 1.), (c1, -1.)):
            m = c >= 0; np.add.at(Ac, c[m], sgn*cr[m]/2); np.add.at(Cx, c[m], sgn*(x0[m]+x1[m])*cr[m]/6)
    X = Cx/Ac*1e6; L = P[:, 0].max()*1e6
    mdot = -D[(18, zid['inlet'])][1].sum(); q = mdot/rho
    inl = F[zid['inlet']][2]; lf = np.hypot(*(P[inl[:, 1]-1]-P[inl[:, 0]-1]).T); Win = lf.sum()*1e6
    pin = float(np.sum(D[(1, zid['inlet'])][1]*lf)/lf.sum())
    Atot = Ac.sum()*1e12
    st = np.linspace(0.15, 1.0, 18)*L; st = st[st <= L-10.0]
    prof = []
    for x in st:
        m = np.abs(X-x) < 0.3; prof.append(np.sum(p[m]*Ac[m])/Ac[m].sum())
    prof = np.array(prof)
    r = dict(cells=nc, L_um=L, area_um2=Atot, inlet_width_um=Win, q_m2s=q, mdot=mdot, p_inlet=pin, stations_um=st.tolist(), p_station=prof.tolist())
    for sig in (1.0, 2.0):
        xs, s, w, I = width_profile(a, b, sig); S = s[-1]; weff = (S/I[-1])**(1/3); tau = S/xs[-1]
        plub = 12*mu*q*I*1e12; ratio = prof/np.interp(st, xs, plub)
        r[f"sig{sig:g}"] = dict(S_um=S, tau=tau, w_eff_um=weff, w_min_um=float(w.min()), w_mean_um=float(np.trapezoid(w, s)/S),
                                ratio_stations=ratio.tolist(), ratio_median=float(np.median(ratio)), ratio_min=float(ratio.min()), ratio_max=float(ratio.max()),
                                ratio_inlet=pin/plub[-1], k_lub=phi*(weff*1e-6)**2/(12*tau**2))
    res[lab] = r
    s1, s2 = r["sig1"], r["sig2"]
    print(f"{lab}: N {nc}, L {L:.1f}, A {Atot:.1f}, Win {Win:.2f}, mdot {mdot:.4e}, p_in {pin:.1f} Pa | sig2: tau {s2['tau']:.3f} weff {s2['w_eff_um']:.3f} wmin {s2['w_min_um']:.2f} ratio med {s2['ratio_median']:.3f} [{s2['ratio_min']:.3f}-{s2['ratio_max']:.3f}] inlet {s2['ratio_inlet']:.3f} | sig1: tau {s1['tau']:.3f} weff {s1['w_eff_um']:.3f} ratio med {s1['ratio_median']:.3f} [{s1['ratio_min']:.3f}-{s1['ratio_max']:.3f}]", flush=True)
json.dump(res, open(os.path.join(HERE, "D_study.json"), "w"), indent=1)
