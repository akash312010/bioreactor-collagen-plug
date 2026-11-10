"""Geometry-only small parameters of the curvilinear expansion for G1-G5 (no CFD needed).
Reports, over the 6/12 um comparison interval and reference smoothing sigma_C: weighted (h^-3) means and maxima of
|h'|, |c'|, |kappa h|, metric minimum 1-|kappa|(|c|+h), and the h^-3-weighted correction <C>."""
import json, numpy as np
import os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'curvi_test.py')).read().split('def run(lab):')[0]; exec(src)
out = {}
for lab, (pk, base, blend) in runs.items():
    A, B = load(CAD+f"collagen_channel_interior_{pk}.csv"); a = smooth_wall(A, 0.75, blend=blend); b = smooth_wall(B, 0.75, blend=blend)
    out[lab] = {}
    for sigc in (1.0, 2.0, 3.0, 4.0):
        C, s, kap, na, nb = frame(a, b, sigc); ds = s[1]-s[0]
        g = lambda f: gaussian_filter1d(f, SIG_D/ds, mode='nearest')
        h = 0.5*(na-nb); cc = 0.5*(na+nb); hs, cs, ks = g(h), g(cc), g(kap)
        h1 = np.gradient(hs, s); h2 = np.gradient(h1, s); c1 = np.gradient(cs, s); c2 = np.gradient(c1, s)
        corr = (-cs*ks + c1**2 + 1.4*h1**2 - 0.1*hs*h2 + (2/15)*hs**2*ks**2 - cs*c2 + 6*cs*c1*h1/hs - 6*cs**2*h1**2/hs**2 + 1.5*cs**2*h2/hs)
        m = (s > 6) & (s < s[-1]-12); w = (hs**-3)[m]; w = w/w.sum()
        out[lab]["sig%g" % sigc] = dict(
            max_abs_hp=float(np.max(np.abs(h1[m]))), wmean_hp2=float(np.sum(w*h1[m]**2)),
            max_abs_cp=float(np.max(np.abs(c1[m]))), max_abs_kh=float(np.max(np.abs(ks[m]*hs[m]))),
            max_abs_k_reach=float(np.max(np.abs(ks[m])*(np.abs(cs[m])+hs[m]))),
            max_abs_c_over_h=float(np.max(np.abs(cs[m]/hs[m]))), wmean_C=float(np.sum(w*corr[m])), max_abs_C=float(np.max(np.abs(corr[m]))),
            frac_absC_gt_05=float(np.mean(np.abs(corr[m]) > 0.5)))
        r = out[lab]["sig%g" % sigc]
        print(lab, sigc, {k: round(v, 3) for k, v in r.items()}, flush=True)
json.dump(out, open('geom_params.json', 'w'), indent=1)
