"""Parameter-free prediction of the CFD/lubrication factor from second-order lubrication theory.
Local correction derived here (sympy, elt2.py) for walls at y = c(x) -+ h(x):
   dp/dx = -(3 mu Q / 2 h^3) [1 + 7/5 h'^2 + 2 c'^2 - 1/10 h h'']
Integrated with weight h^-3 (parallel ends): factor = 1 + < 6/5 h'^2 + 2 c'^2 >_w  (+ curvature terms, not included)."""
import sys, json, numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from smooth import load, smooth_wall
import types
T = types.SimpleNamespace(); _src = open(os.path.join(HERE, 'theory_local.py')).read().split("res = {}; local = {}")[0]
_ns = {"__file__": os.path.join(HERE, "theory_local.py")}; exec(_src, _ns); T.geometry = _ns["geometry"]          # reuse geometry() without rerunning the CFD loop
CAD = os.path.join(HERE, "..", "..", "data", "geometry") + os.sep
meas = json.load(open(os.path.join(HERE, "theory_local.json")))["summary"]
blend = {"p28": 3.0, "p31": 3.0, "p35": 3.0, "p40": 3.0, "p45": 1.0}
out = {}
for lab in ("p28", "p31", "p35", "p40", "p45"):
    A, B = load(CAD+f"collagen_channel_interior_{lab}.csv"); a = smooth_wall(A, 0.75, blend=blend[lab]); b = smooth_wall(B, 0.75, blend=blend[lab])
    G = T.geometry(a, b); s, w, ta, tb, kap = G["s"], G["w"], G["ta"], G["tb"], G["kappa"]
    ds = s[1]-s[0]; inner = (s > 6) & (s < s[-1]-12)
    cp = 0.5*(ta+tb); hp = 0.5*(ta-tb)
    wt = w[inner]**-3; wt = wt/wt.sum()
    th2 = np.sum(wt*(1.2*hp[inner]**2 + 2*cp[inner]**2))                    # epsilon^2 theory
    # partial resummation: exact for a uniformly tilted channel (1+c'^2)^2, times the width-variation term
    res = np.sum(wt*((1+cp[inner]**2)**2*(1+1.2*hp[inner]**2)))-1
    curv = np.sum(wt*(w[inner]*kap[inner])**2)
    out[lab] = dict(C_cfd=meas[lab]["C_inner"], C_theory_eps2=1+th2, C_theory_resummed=1+res, mean_hp2=float(np.sum(wt*hp[inner]**2)),
                    mean_cp2=float(np.sum(wt*cp[inner]**2)), mean_wk2=float(curv))
    o = out[lab]; print(f"{lab}: CFD {o['C_cfd']:.3f} | theory eps^2 {o['C_theory_eps2']:.3f} | resummed {o['C_theory_resummed']:.3f} | <h'^2> {o['mean_hp2']:.3f} <c'^2> {o['mean_cp2']:.3f} <(w kappa)^2> {o['mean_wk2']:.3f}", flush=True)
json.dump(out, open(os.path.join(HERE, "theory_predict.json"), "w"), indent=1)
