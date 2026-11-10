import numpy as np, json
from spread import solve
R, H = 16.664, 26.873
out = {"a_over_R": [], "F": [], "F_holm": []}
for aR in (0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0):
    a = aR*R; F, _ = solve(R, H, a, a, 0.08 if aR > 0.06 else 0.04)
    holm = 1+(2/(4*a))/(H/(np.pi*R**2))
    out["a_over_R"].append(aR); out["F"].append(F); out["F_holm"].append(holm); print(aR, round(F, 3), round(holm, 3), flush=True)
dev = [solve(R, H, 1.93, 2.4995, h)[0] for h in (0.1, 0.05, 0.025)]
out["device"] = dict(a1=1.93, a2=2.4995, F=dev)
json.dump(out, open("fcurve.json", "w"), indent=1)
print("device", dev)
