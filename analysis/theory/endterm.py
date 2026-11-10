import sys, json, numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from smooth import load, smooth_wall
ns = {"__file__": os.path.join(HERE, "theory_local.py")}; exec(open(os.path.join(HERE, 'theory_local.py')).read().split("res = {}; local = {}")[0], ns); geometry = ns["geometry"]
CAD = os.path.join(HERE, "..", "..", "data", "geometry") + os.sep
blend = {"p28": 3.0, "p31": 3.0, "p35": 3.0, "p40": 3.0, "p45": 1.0}
out = {}
for lab in blend:
    A, B = load(CAD+f"collagen_channel_interior_{lab}.csv"); a = smooth_wall(A, 0.75, blend=blend[lab]); b = smooth_wall(B, 0.75, blend=blend[lab])
    G = geometry(a, b); s, w, ta, tb = G["s"], G["w"], G["ta"], G["tb"]; h = w/2; hp = np.gradient(h, s)
    rows = []
    for lo, hi in ((6, 12), (4, 8), (8, 16), (10, 20)):
        m = (s > lo) & (s < s[-1]-hi); ds = s[1]-s[0]
        cp = 0.5*(ta+tb); hq = 0.5*(ta-tb)
        I0 = np.sum(h[m]**-3)*ds; I2 = np.sum(h[m]**-3*(1.2*hq[m]**2+2*cp[m]**2))*ds
        i0, i1 = np.where(m)[0][[0, -1]]
        E = -0.1*(hp[i1]/h[i1]**2 - hp[i0]/h[i0]**2)
        rows.append(dict(cut=(lo, hi), factor_no_end=1+I2/I0, end_over_I0=E/I0, factor_with_end=1+(I2+E)/I0))
    out[lab] = rows
    print(lab, " | ".join("cut %s: F %.3f, end %.4f -> %.3f" % (r["cut"], r["factor_no_end"], r["end_over_I0"], r["factor_with_end"]) for r in rows), flush=True)
json.dump(out, open("endterm.json", "w"), indent=1, default=str)
