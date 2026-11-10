"""Where does the curvilinear metric 1 - kappa*n reach zero inside the channel? Fraction of the 6/12 um comparison
interval (by arc length and by h^-3 resistance weight) with |kappa|(|c|+h) >= 1, per channel and reference smoothing."""
import json, numpy as np
import os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'curvi_test.py')).read().split('def run(lab):')[0]; exec(src)
out = {}
for lab, (pk, base, blend) in runs.items():
    A, B = load(CAD+f"collagen_channel_interior_{pk}.csv"); a = smooth_wall(A, 0.75, blend=blend); b = smooth_wall(B, 0.75, blend=blend)
    out[lab] = {}
    for sigc in (1.0, 2.0, 3.0, 4.0):
        C, s, kap, na, nb = frame(a, b, sigc); ds = s[1]-s[0]; g = lambda f: gaussian_filter1d(f, SIG_D/ds, mode='nearest')
        h = 0.5*(na-nb); cc = 0.5*(na+nb); hs, cs, ks = g(h), g(cc), g(kap)
        m = (s > 6) & (s < s[-1]-12); reach = np.abs(ks)*(np.abs(cs)+hs); bad = reach >= 1
        dsi = np.gradient(s); w = hs**-3*dsi
        out[lab]["sig%g" % sigc] = dict(len_frac=float(np.sum(dsi[m & bad])/np.sum(dsi[m])), res_frac=float(np.sum(w[m & bad])/np.sum(w[m])))
        print(lab, sigc, {k: round(v, 4) for k, v in out[lab]["sig%g" % sigc].items()}, flush=True)
json.dump(out, open('metric_check.json', 'w'), indent=1)
