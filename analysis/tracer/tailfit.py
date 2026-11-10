"""Late-time decay constant of the unreplaced volume 1 - c_V for every valid tracer run, and its ratio to t_ideal,
compared with the volume/outlet replacement ratio t99,V / t99,out. Uses existing Fluent report files only."""
import os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
U = os.environ.get("CFD_RUNS", "CFD_runs_new/")  # Fluent case/data NOT in this repository
E = json.load(open(os.path.join(HERE, "E_all.json")))
PATH = {"E11": ("E11_run/E11/", ""), "E12": ("E12_run/E12/", ""), "E1": ("E1/E1_5ml/", ""), "E2": ("E2_run/E2/", ""),
        "E3": ("E3_run/New/New_E3/", ""), "E4": ("E4_run/E4/", ""), "E6": ("E6/E6_30Ml_fine_mesh_1..7E-6/", ""),
        "E7": ("E7_run/E7/", "_1_1"), "E8": ("E8_run/New2/New_E8/", "_1_1"), "E14": ("../New run R1_E14/E14_tracer/", "_1_1"), "E9": ("E9_run/E9/", ""), "E10": ("E10_run/E10/", "")}
G = {"p28": "G1", "p35": "G3", "p45": "G5"}
out = {}
for r, (p, s) in PATH.items():
    d = np.loadtxt(os.path.join(U, p, "vol_avg-rfile%s.out" % s), skiprows=4); t, c = d[:, 2], d[:, 1]
    y = 1 - c; tid = E[r]["t_ideal_ms"]*1e-3
    res = {}
    for lo, hi in ((1e-3, 3e-2), (3e-3, 1e-1)):
        m = (y > lo) & (y < hi) & (t > t[np.argmax(y < 0.5)])
        if m.sum() < 8: res["%g_%g" % (lo, hi)] = None; continue
        k = np.polyfit(t[m], np.log(y[m]), 1)[0]; res["%g_%g" % (lo, hi)] = -1/k
    tau = res["0.001_0.03"] or res["0.003_0.1"]
    out[r] = dict(geom=G[E[r]["geom"]], Pe=E[r]["Pe_w"], t_ideal_ms=tid*1e3, tau_ms={k: (v*1e3 if v else None) for k, v in res.items()},
                  tau_over_tid=tau/tid, ratio_V_out=E[r]["vol_t99"]/E[r]["mass_t99"], vol_t99=E[r]["vol_t99"], out_t99=E[r]["mass_t99"],
                  ymin=float(y.min()))
json.dump(out, open(os.path.join(HERE, "tailfit.json"), "w"), indent=1)
for r, v in sorted(out.items(), key=lambda x: x[1]["tau_over_tid"]):
    print("%-4s %s Pe %6.1f tid %6.2f ms tau %s ms  tau/tid %5.2f  t99V/t99out %.2f  (min 1-c %.1e)" % (r, v["geom"], v["Pe"], v["t_ideal_ms"],
          {k: (round(x, 2) if x else None) for k, x in v["tau_ms"].items()}, v["tau_over_tid"], v["ratio_V_out"], v["ymin"]))
