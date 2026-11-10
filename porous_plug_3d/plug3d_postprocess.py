"""Post-process the 3D porous-plug runs R1-R3 from the Fluent surface-integral reports.

    python plug3d_postprocess.py

Reads reports/<run>/*.txt and prints the port-to-port pressure drop, the feeding factor
F3D = dp * k * A / (mu * delta * Q), and the convergence and uniformity checks.
Q is the flow actually delivered in each run (inlet mass flow / rho)."""
import math, os, re, json

HERE = os.path.dirname(os.path.abspath(__file__))
mu, rho = 9.64e-4, 1012.0
R, delta, k = 16.664e-3, 26.873e-3, 5.566342952e-13
A = math.pi * R**2
F_axisym = 7.123813486          # axisymmetric Darcy calculation (analysis/porous/spread.py)

RUNS = {"R1": "mesh M1 (production), velocity inlet, nominal 5 mL/min",
        "R2": "mesh M2 (port-refined), velocity inlet, nominal 5 mL/min",
        "R3": "mesh M1 (production), mass-flow inlet, 30 mL/min"}


def read_report(path):
    """Return {surface: value} from a Fluent 'Surface Integral Report' text file."""
    out = {}
    for line in open(path):
        m = re.match(r"\s*([A-Za-z_]+)\s+([-+0-9.eE]+)\s*$", line)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def F3D(dp, Q):
    return dp * k * A / (mu * delta * Q)


results = {}
for run, desc in RUNS.items():
    d = os.path.join(HERE, "reports", run)
    r = {n: read_report(os.path.join(d, n + ".txt")) for n in
         ("area", "mass_flow", "static_pressure_area_weighted", "static_pressure_mass_weighted",
          "total_pressure_area_weighted", "total_pressure_mass_weighted")}
    mdot_in, mdot_out = r["mass_flow"]["inlet"], r["mass_flow"]["outlet"]
    Q = mdot_in / rho
    dp = {key: r[key]["port_in"] - r[key]["port_out"] for key in r if "pressure" in key}
    res = dict(description=desc, Q_mL_min=Q * 6e7,
               mass_imbalance_rel=abs(mdot_in + mdot_out) / mdot_in,
               dp_static_area_Pa=dp["static_pressure_area_weighted"],
               dp_static_mass_Pa=dp["static_pressure_mass_weighted"],
               dp_total_area_Pa=dp["total_pressure_area_weighted"],
               F3D_area=F3D(dp["static_pressure_area_weighted"], Q),
               F3D_mass=F3D(dp["static_pressure_mass_weighted"], Q),
               conduits_Pa=(r["static_pressure_area_weighted"]["inlet"] - r["static_pressure_area_weighted"]["outlet"])
                           - dp["static_pressure_area_weighted"])
    res["F3D_over_axisym"] = res["F3D_area"] / F_axisym
    results[run] = res
    print(f"{run}: {desc}")
    print(f"   Q = {res['Q_mL_min']:.3f} mL/min, mass imbalance {res['mass_imbalance_rel']:.1e}")
    print(f"   dp(port_in - port_out): static area {res['dp_static_area_Pa']/1e3:.3f} kPa, "
          f"static mass {res['dp_static_mass_Pa']/1e3:.3f} kPa, total area {res['dp_total_area_Pa']/1e3:.3f} kPa")
    print(f"   F3D = {res['F3D_area']:.3f} (area), {res['F3D_mass']:.3f} (mass); "
          f"F3D / F_axisym = {res['F3D_over_axisym']:.4f}; conduits {res['conduits_Pa']:.2f} Pa")

r1, r2, r3 = (results[x] for x in ("R1", "R2", "R3"))
print(f"Mesh change M1 -> M2 in F3D: {100 * (r2['F3D_area'] / r1['F3D_area'] - 1):+.2f} %")
lin = (r3["dp_static_area_Pa"] / r3["Q_mL_min"]) / (r1["dp_static_area_Pa"] / r1["Q_mL_min"])
print(f"Linearity, (dp/Q at 30 mL/min) / (dp/Q at 5 mL/min) = {lin:.4f} (Darcy: 1)")
json.dump(results, open(os.path.join(HERE, "plug3d_results.json"), "w"), indent=1)
