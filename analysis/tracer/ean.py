import numpy as np, json
import os
U = os.environ.get("CFD_RUNS", "CFD_runs_new/")  # Fluent case/data NOT in this repository
D = 6.7e-10
geo = {"p28": dict(A=477.97e-12, weff=4.53e-6, wbar=6.24e-6, L=71.3e-6, Win=2.715e-6),
       "p35": dict(A=768.43e-12, weff=6.11e-6, wbar=7.66e-6, L=90.3e-6, Win=4.6088e-6),
       "p45": dict(A=972.68e-12, weff=7.81e-6, wbar=9.15e-6, L=96.7e-6, Win=12.195e-6)}
runs = {"E11": ("p35", 0.000653, "E11_run/E11/", ""), "E12": ("p35", 0.00392, "E12_run/E12/", ""), "E1": ("p35", 0.00769, "E1/E1_5ml/", ""), "E2": ("p35", 0.0140, "E2_run/E2/", ""), "E3": ("p35", 0.0209, "E3_run/New/New_E3/", ""), "E4": ("p35", 0.0292, "E4_run/E4/", ""),
        "E6": ("p35", 0.0469, "E6/E6_30Ml_fine_mesh_1..7E-6/", ""), "E7": ("p28", 0.00769, "E7_run/E7/", "_1_1"), "E8": ("p28", 0.0469, "E8_run/New2/New_E8/", "_1_1"), "E9": ("p45", 0.00769, "E9_run/E9/", ""), "E10": ("p45", 0.0469, "E10_run/E10/", "")}
def rd(f): d = np.loadtxt(f, skiprows=4); return d[:, 2], d[:, 1]
out = {}
for lab, (g, u, path, sfx) in runs.items():
    G = geo[g]; q = u*G["Win"]; tid = G["A"]/q
    t, cm = rd(U+path+f"report-def-0-rfile{sfx}.out"); _, ca = rd(U+path+f"area_weigted_avg-rfile{sfx}.out"); _, cv = rd(U+path+f"vol_avg-rfile{sfx}.out")
    i0, i1 = np.searchsorted(t, 0.05*tid), np.searchsorted(t, 0.4*tid)
    slope = (cv[i1]-cv[i0])/(t[i1]-t[i0])*tid
    tf = lambda c, f: float(t[np.argmax(c >= f)]/tid) if c.max() >= f else float("nan")
    Pe = q/D; Lam = G["A"]*D/(q*G["weff"]**2)
    r = dict(geom=g, u_in=u, q=q, t_ideal_ms=tid*1e3, vol_slope=slope, Pe_q=Pe, Pe_w=(q/G["wbar"])*G["weff"]/D, Lambda=Lam, end_t_over_tid=float(t[-1]/tid),
             first=float(t[np.argmax(cm > 1e-3)]/tid))
    for nm, c in (("mass", cm), ("area", ca), ("vol", cv)):
        for f in (0.5, 0.9, 0.99): r[f"{nm}_t{int(f*100)}"] = tf(c, f)
    out[lab] = r
    print(f"{lab} {g} u={u}: t_id {tid*1e3:.3f} ms, slope {slope:.4f}, Pe_w {r['Pe_w']:.1f}, Lambda {Lam:.3f}, end {r['end_t_over_tid']:.2f} | mass t50/90/99 {r['mass_t50']:.3f}/{r['mass_t90']:.3f}/{r['mass_t99']:.3f} | area t99 {r['area_t99']:.3f} | vol t90/99 {r['vol_t90']:.3f}/{r['vol_t99']:.3f} | first {r['first']:.3f} | t99 mass {r['mass_t99']*tid*1e3:.1f} ms")
json.dump(out, open("E_all.json", "w"), indent=1)
