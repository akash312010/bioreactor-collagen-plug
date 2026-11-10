"""Recompute every headline number of the manuscript from the files in this package (no Fluent case/data needed).
Run:  python reproduce_headlines.py      -> prints a table and writes headline_numbers.json
Inputs used (all in this package):
  analysis/theory/D_summary.json        channel-informed permeabilities k_cfd (G1-G5)
  analysis/porous/fcurve.json           axisymmetric Darcy spreading factor F (equal ports + device ports)
  analysis/theory/curvi_test.json       curvilinear second-order theory vs CFD (5 channels x 4 smoothings x 4 cuts)
  analysis/theory/metric_check.json     where the curvilinear metric 1 - kappa*n vanishes
  analysis/theory/geom_params.json      size of the second-order correction
  analysis/tracer/E_all.json            t_ideal, Pe, outlet/volume t50/t90/t99 per run
  data/tracer_reports/<run>/*.out       raw Fluent report series (outlet mass-weighted, volume average)
  data/device/Part1_data_reconciled.csv  archived device losses, checked against the original workbook
The open-device dimensionless calculation uses the explicitly identified 1M
reference series. The 2M field panel is a separate 15 mL/min field.
The exact signed-metric audit is in data/checks/signed_metric_audit.json; the older
metric_check.json records a sufficient conservative bound only."""
import os, json, csv, numpy as np
from scipy.stats import spearmanr
H = os.path.dirname(os.path.abspath(__file__)); J = lambda *p: json.load(open(os.path.join(H, *p)))
out = {}
# ---------------- constants (as meshed) ----------------
mu, rho, D = 9.64e-4, 1012.0, 6.7e-10
R, delta, phi = 16.664e-3, 26.873e-3, 0.27
A = np.pi*R**2; MLMIN = 1e-6/60
# ---------------- device (open chamber) ----------------
rows=list(csv.DictReader(open(os.path.join(H,'data','device','Part1_data_reconciled.csv'))))
Ro = 2.5e-3
loss={int(r['Q_ml_min']):float(r['dP_Pa']) for r in rows if r['mesh']=='1M' and int(r['Q_ml_min']) in (5,30)}
dev = {}
for Q, dp in loss.items():
    U = Q*MLMIN/(np.pi*Ro**2); Re = rho*U*2*Ro/mu; dev[Q] = dict(U_o=U, Re_o=Re, EuRe=dp*2*Ro/(mu*U))
out["device"] = dev
out['device_provenance']={'reference_grid':'1M','source':'data/device/Part1_data_reconciled.csv',
 'reference_losses_Pa':loss,'all_grid_loss_range_Pa':[min(float(r['dP_Pa']) for r in rows),max(float(r['dP_Pa']) for r in rows)],
 'outlet_diameter_m':2*Ro,'velocity_convention':'mass-conserving Q divided by nominal circular outlet area'}
# ---------------- permeability, capacity, spreading ----------------
names = {"p28": "G1", "p31": "G2", "p35": "G3", "p40": "G4", "p45": "G5"}
k = {names[r["threshold"]]: r["k_cfd"] for r in J("analysis", "theory", "D_summary.json")}
slope = {g: kk*A/(mu*delta)*1e3/MLMIN for g, kk in k.items()}       # mL/min per kPa, uniform faces
FC = J("analysis", "porous", "fcurve.json"); Fdev = FC["device"]["F"]; F = Fdev[-1]
a1, a2 = FC["device"]["a1"], FC["device"]["a2"]; Rm, Hm = R*1e3, delta*1e3
Fcorr = 1+np.pi*Rm**2/(4*Hm)*((1-a1/Rm)**1.5/a1+(1-a2/Rm)**1.5/a2)
out["capacity"] = dict(k_m2=k, slope_uniform_mLmin_per_kPa=slope, F_device_grids=Fdev, F_finest=F, F_corr=Fcorr,
    F_corr_rel_err=Fcorr/F-1, slope_portfed_G3=slope["G3"]/F,
    dp_5mLmin_uniform_G3_kPa=5/slope["G3"], dp_5mLmin_portfed_G3_kPa=5*F/slope["G3"],
    dp_5mLmin_portfed_range_kPa=[5*F/max(slope.values()), 5*F/min(slope.values())],
    F_full_face_limit=FC["F"][-1])
# ---------------- dye / recalled reversal band consistency (40-80 mL/min) ----------------
out["dye_check"] = {Q: dict(uniform_traverse_s=delta/(Q*MLMIN/(A*phi)), dp_uniform_kPa=Q/slope["G3"], dp_portfed_kPa=Q*F/slope["G3"])
                    for Q in (40, 80)}
# ---------------- curvilinear theory ----------------
C = J("analysis", "theory", "curvi_test.json"); cuts = ["4_8", "6_12", "8_16", "10_20"]
th = {}
for g in ["G1", "G2", "G3", "G4", "G5"]:
    r = lambda sig, key="theory_over_cfd": [C[g]["sig%d" % s][c][key] for s in sig for c in cuts]
    th[g] = dict(all16=[min(r([1, 2, 3, 4])), max(r([1, 2, 3, 4]))], sig2to4=[min(r([2, 3, 4])), max(r([2, 3, 4]))],
                 sig4=[min(r([4])), max(r([4]))], no_offset_curv_sig2to4=[min(r([2, 3, 4], "no_curv_over_cfd")), max(r([2, 3, 4], "no_curv_over_cfd"))])
M = J("analysis", "theory", "metric_check.json"); G = J("analysis", "theory", "geom_params.json")
out["theory"] = dict(ratio=th, conservative_metric_bound_failure=M,
    filtered_signed_metric_audit=J('data','checks','signed_metric_audit.json'),
    raw_geometry_metric_audit=J('data','checks','theory_reproduction_summary.json')['metric_checks'],
    corr_wmean_G3G5_sig2to4=[min(G[g]["sig%d" % s]["wmean_C"] for g in ["G3", "G4", "G5"] for s in (2, 3, 4)),
                             max(G[g]["sig%d" % s]["wmean_C"] for g in ["G3", "G4", "G5"] for s in (2, 3, 4))],
    frac_absC_gt05_G3G5_sig2to4=[min(G[g]["sig%d" % s]["frac_absC_gt_05"] for g in ["G3", "G4", "G5"] for s in (2, 3, 4)),
                                 max(G[g]["sig%d" % s]["frac_absC_gt_05"] for g in ["G3", "G4", "G5"] for s in (2, 3, 4))])
try:   # Tavakol et al. (2017) Eq. 19 one-flat-wall limit
    import sympy as sp
    h, hs, hss, c, cs, css, kap, Gg, Gs, Gss = sp.symbols('h h_s h_ss c c_s c_ss kappa G G_s G_ss')
    Cc = (-c*kap+cs**2+sp.Rational(7, 5)*hs**2-hss*h/10+sp.Rational(2, 15)*h**2*kap**2-c*css+6*c*cs*hs/h
          - 6*c**2*hs**2/h**2+sp.Rational(3, 2)*c**2*hss/h)
    red = sp.expand(sp.simplify(Cc.subs({kap: 0, c: h, cs: hs, css: hss}).subs({h: Gg/2, hs: Gs/2, hss: Gss/2})))
    tav = sp.expand((sp.Rational(36, 5)*Gs**2/Gg**3+sp.Rational(6, 5)*Gss/Gg**2)/(12/Gg**3))
    const = sp.expand(sp.simplify(Cc.subs({c: 0, cs: 0, css: 0, hs: 0, hss: 0}).subs({h: Gg/2})))
    out["theory"]["tavakol_check"] = dict(reduction=str(red), tavakol_eq19_normalised=str(tav), identical=bool(sp.simplify(red-tav) == 0),
                                          constant_width_curvature=str(const))
except ImportError:
    out["theory"]["tavakol_check"] = "sympy not installed"
# ---------------- tracer: raw reports -> t99, tail fit, ratio, Spearman ----------------
E = J("analysis", "tracer", "E_all.json"); tr = {}
for run in ["E11", "E12", "E1", "E2", "E3", "E4", "E6", "E7", "E14", "E8", "E9", "E10"]:
    rd = lambda f: np.loadtxt(os.path.join(H, "data", "tracer_reports", run, f), skiprows=4)
    o = rd("outlet_massweighted.out"); v = rd("volume_average.out"); tid = E[run]["t_ideal_ms"]*1e-3
    t99 = lambda d: float(d[np.argmax(d[:, 1] >= 0.99), 2]/tid) if d[:, 1].max() >= 0.99 else float("nan")
    t, y = v[:, 2], 1-v[:, 1]; m = (y > 1e-3) & (y < 3e-2) & (t > t[np.argmax(y < 0.5)])
    tau = -1/np.polyfit(t[m], np.log(y[m]), 1)[0]
    m_alt=(y>0.003)&(y<0.1)&(t>t[np.argmax(y<0.5)])
    tau_alt=-1/np.polyfit(t[m_alt],np.log(y[m_alt]),1)[0]
    tr[run] = dict(geom=names[E[run]["geom"]], Pe=E[run]["Pe_w"], t_ideal_ms=tid*1e3, t99_out=t99(o), t99_V=t99(v),
                   tau_ms=tau*1e3, Lambda_e=tau/tid, tau_sensitivity_ms=tau_alt*1e3,
                   Lambda_e_sensitivity=tau_alt/tid, Re_w=E[run]["Pe_w"]*rho*D/mu)
for r in tr.values(): r["ratio_V_out"] = r["t99_V"]/r["t99_out"]
lam = [r["Lambda_e"] for r in tr.values()]; rat = [r["ratio_V_out"] for r in tr.values()]
cons = [1-r["ratio_V_out"] for r in tr.values() if r["ratio_V_out"] < 1]
out["tracer"] = dict(runs=tr, spearman=float(spearmanr(lam, rat).statistic), n_conservative=len(cons),
                     spearman_sensitivity=float(spearmanr([r['Lambda_e_sensitivity'] for r in tr.values()],rat).statistic),
                     conservative_margin=[min(cons), max(cons)], Sc=mu/(rho*D),
                     Re_w_range=[min(r["Re_w"] for r in tr.values()), max(r["Re_w"] for r in tr.values())])
# ---------------- solute balance (Eq. 8, Table S5) ----------------
bal = {}
for run in tr:
    o = np.loadtxt(os.path.join(H, "data", "tracer_reports", run, "outlet_massweighted.out"), skiprows=3)
    v = np.loadtxt(os.path.join(H, "data", "tracer_reports", run, "volume_average.out"), skiprows=3)
    tid = tr[run]["t_ideal_ms"]*1e-3; t, co, cv = v[:, 2], o[:, 1], v[:, 1]; y = 1-cv
    m = (y > 1e-3) & (y < 3e-2) & (t > t[np.argmax(y < 0.5)])
    bal[run] = dict(max_deviation=float(np.abs(np.gradient(cv, t)*tid-(1-co))[t/tid > 0.5].max()),
                    residual_ratio=float(np.median(y[m]/(1-co[m]))))
out["solute_balance"] = bal
# ---------------- plug scale: axisymmetric advective transit times and 3D tracer R1T ----------------
out["plug_transport_grids"] = J("analysis", "porous", "plug_transport.json")
import importlib.util
spec = importlib.util.spec_from_file_location("r1t", os.path.join(H, "porous_plug_3d", "tracer_R1T", "analyze_r1t.py"))
r1t = importlib.util.module_from_spec(spec); spec.loader.exec_module(r1t)
out["R1T"] = r1t.r1t_metrics()[0]
# ---------------- image dark-area fraction (porosity sensitivity) ----------------
from PIL import Image
img = np.asarray(Image.open(os.path.join(H, "data", "geometry", "collagen_SHG_Niestrawska2022_Fig4.png")).convert("L")).astype(int)
out["image_dark_fraction"] = {names_pct: float((img <= int(round(e["threshold_value"]*255))).mean())
                              for e in J("data", "geometry", "cad_info.json") if e["variant"] == "interior"
                              for names_pct in [{28: "G1", 31: "G2", 35: "G3", 40: "G4", 45: "G5"}[e["threshold_pct"]]]}
json.dump(out, open(os.path.join(H, "headline_numbers.json"), "w"), indent=1, default=float)
# ---------------- print ----------------
c = out["capacity"]; t = out["tracer"]
print("Device: Re_o %.1f-%.1f, Eu*Re %.0f-%.0f (archived losses)" % (dev[5]["Re_o"], dev[30]["Re_o"], dev[5]["EuRe"], dev[30]["EuRe"]))
print("k(G3) %.3e m2; slopes uniform G1-G5 %.3f-%.3f, G3 %.4f mL/min/kPa" % (k["G3"], min(slope.values()), max(slope.values()), slope["G3"]))
print("F grids %s, finest %.3f, F_corr %.3f (%+.1f%%), full-face %.3f" % (np.round(Fdev, 3), F, Fcorr, 100*c["F_corr_rel_err"], c["F_full_face_limit"]))
print("5 mL/min: uniform %.2f kPa, port-fed %.1f kPa (G3), port-fed G1-G5 %.0f-%.0f kPa; port-fed slope %.3f" % (
    c["dp_5mLmin_uniform_G3_kPa"], c["dp_5mLmin_portfed_G3_kPa"], *c["dp_5mLmin_portfed_range_kPa"], c["slope_portfed_G3"]))
for Q, d in out["dye_check"].items(): print("Q=%d: uniform traverse %.2f s, dp uniform %.0f kPa, port-fed %.0f kPa" % (Q, d["uniform_traverse_s"], d["dp_uniform_kPa"], d["dp_portfed_kPa"]))
for g, d in th.items(): print("%s theory/CFD: all16 %.3f-%.3f | sig2-4 %.3f-%.3f | sig4 %.3f-%.3f | no offset/curv %.3f-%.3f" % (g, *d["all16"], *d["sig2to4"], *d["sig4"], *d["no_offset_curv_sig2to4"]))
exact=out['theory']['raw_geometry_metric_audit']
print('Raw geometry nonpositive metric fractions (sig2, 6/12 cuts):', {g:round(exact[g]['sig2']['6_12']['raw_metric_nonpositive_fraction'],5) for g in exact})
print("correction <C> G3-G5 sig2-4: %.2f-%.2f; |C|>0.5 on %.0f-%.0f%%" % (*out["theory"]["corr_wmean_G3G5_sig2to4"], *[100*x for x in out["theory"]["frac_absC_gt05_G3G5_sig2to4"]]))
print("Tavakol:", out["theory"]["tavakol_check"])
for r, d in tr.items(): print("%-4s %s Pe %6.1f tid %7.2f ms  t99o %.3f  t99V %.3f  ratio %.2f  tau %.2f ms  Lambda %.2f" % (r, d["geom"], d["Pe"], d["t_ideal_ms"], d["t99_out"], d["t99_V"], d["ratio_V_out"], d["tau_ms"], d["Lambda_e"]))
print("Spearman %.3f; conservative in %d runs by %.0f-%.0f%%; Sc %.0f; Re_w %.4f-%.2f" % (t["spearman"], t["n_conservative"], 100*t["conservative_margin"][0], 100*t["conservative_margin"][1], t["Sc"], *t["Re_w_range"]))
print("Solute balance: max deviation %.1e; residual ratio vs Lambda_e:" % max(b["max_deviation"] for b in bal.values()),
      {r: (round(float(b["residual_ratio"]), 2), round(float(tr[r]["Lambda_e"]), 2)) for r, b in bal.items()})
g = out["plug_transport_grids"][-1]
print("Plug scale (axisymmetric, h=%.3f mm): first arrival %.2f, outlet t99 %.2f, interior t99 %.2f mean residence times" % (g["h_mm"], g["first_arrival"], g["outlet"]["t99"], g["interior"]["t99"]))
q = out["R1T"]
print("R1T (3D): plug residence %.1f s vs %.1f s; outlet t50/t99 %.2f/%.2f, interior t50/t99 %.2f/%.2f; tail %.2f; residual ratio %.2f" % (
    q["plug_residence_from_integral_s"], q["tbar_s"], q["out_t50"], q["out_t99"], q["vol_t50"], q["vol_t99"], q["tail_time_constant_out"], q["tail_residual_ratio_median"]))
print("Image dark-area fraction:", {g_: round(v, 3) for g_, v in out["image_dark_fraction"].items()})
