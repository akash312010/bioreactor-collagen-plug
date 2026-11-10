# 3D porous-plug simulation (runs R1-R3)

Tests the feeding factor F = 7.12 of the axisymmetric Darcy calculation in the full 3D
device, with the actual inlet and outlet ports, at the G3 permeability. Results are in
Table IV of the paper.

## Geometry and meshes

| File | Content |
|---|---|
| `bioreactor_plug3d_3bodies.step` | Inlet conduit, chamber (plug) and outlet conduit as three solids with shared faces (mm) |
| `make_plug3d_geometry.py`, `make_plug3d_mesh.py` | Geometry and tetrahedral mesh with gmsh, written as Fluent ASCII `.msh` (m) |
| `msh_check.py` | Check of a written mesh: face orientation, zone areas, volumes, skewness |
| `bioreactor_plug3d_M1.msh.gz` | Production mesh M1: 630 657 cells, port-edge size 0.10 mm |
| `bioreactor_plug3d_M2.msh.gz` | Port-refined mesh M2: 1 385 111 cells, port-edge size 0.05 mm |

Dimensions (as in the open-device grids and the axisymmetric model):

| Part | Radius (mm) | Length (mm) |
|---|---|---|
| Inlet conduit | 1.930 | 31.680 |
| Chamber / plug | 16.664 | 26.873 |
| Outlet conduit | 2.4995 | 1.777 |

Zones: `inlet`, `outlet`; cell zones `conduit_in`, `plug`, `conduit_out`; interior sampling
discs `port_in` (top of the plug) and `port_out` (bottom of the plug); walls
`wall_conduit_in`, `wall_plug` (side wall and sealed roof and floor annuli), `wall_conduit_out`.

## Solver setup (ANSYS Fluent 2025 R2)

- Pressure-based, steady, laminar; gravity off.
- ρ = 1012 kg m⁻³, μ = 9.64 × 10⁻⁴ Pa s.
- `plug`: porous zone, isotropic viscous resistance 1/k = 1.79651 × 10¹² m⁻²
  (k = 5.566 × 10⁻¹³ m²), inertial resistance 0, porosity 0.27, superficial-velocity
  formulation. The conduits are plain fluid.
- Outlet: pressure outlet, 0 Pa gauge. Walls: no slip.

| Run | Mesh | Inlet | Delivered flow |
|---|---|---|---|
| R1 | M1 | velocity inlet, 7.12 × 10⁻³ m s⁻¹ | 4.986 mL/min |
| R2 | M2 | velocity inlet, 7.12 × 10⁻³ m s⁻¹ | 4.987 mL/min |
| R3 | M1 | mass-flow inlet, 5.06 × 10⁻⁴ kg s⁻¹ | 30.00 mL/min |

The meshed inlet disc is 0.28 % smaller than the exact disc, so the velocity inlet delivers
slightly less than 5 mL/min. F3D uses the delivered flow of each run.

## Reports and post-processing

`reports/R1`, `reports/R2`, `reports/R3` hold the Fluent surface-integral reports of the
converged solutions: surface areas, mass flow, and area- and mass-weighted static and total
pressure on `inlet`, `outlet`, `port_in` and `port_out`.

```bash
python plug3d_postprocess.py
```

computes, for each run, the delivered flow, the mass imbalance, the port-to-port pressure
drop Δp and

F3D = Δp k A / (μ δ Q),  A = πR², δ = 26.873 mm,

and writes `plug3d_results.json`. Output:

| Run | Q (mL/min) | Δp (kPa) | F3D | Relative mass imbalance |
|---|---|---|---|---|
| R1 | 4.986 | 31.56 | 7.120 | 5 × 10⁻⁸ |
| R2 | 4.987 | 31.47 | 7.097 | 4 × 10⁻⁸ |
| R3 | 30.00 | 190.0 | 7.123 | 2 × 10⁻⁷ |

Axisymmetric value: F = 7.124. Mesh change M1 to M2: −0.3 %. Δp/Q at 30 mL/min over Δp/Q at
5 mL/min: 1.0004. Area- and mass-weighted port pressures agree within 0.01 %.

## 3D tracer run R1T (`tracer_R1T/`)

A passive scalar on the frozen flow of R1, testing the plug-scale replacement curves of the
axisymmetric advective calculation (`analysis/porous/plug_transport.py`).

| Setting | Value |
|---|---|
| Flow | converged R1 solution, flow equations switched off |
| Scalar | one user-defined scalar, flux function mass flow rate |
| Diffusivity | ρD = 6.78 × 10⁻⁷ kg m⁻¹ s⁻¹ (D = 6.7 × 10⁻¹⁰ m² s⁻¹), all zones |
| Porous zone | porosity 0.27 in the storage term |
| Boundaries | inlet value 1 from t = 0; zero diffusive flux at the outlet and walls |
| Schemes | second-order upwind; bounded second-order implicit in time |
| Time step, duration | 0.05 s, 420 s (8400 steps) |

Files: `reports/outlet_massweighted.out` and `reports/plug_volume_average.out` (Fluent report series:
time step, value, flow time), `reports/area.txt`, `reports/velocity_magnitude_area_weighted.txt`, and
`axial_plane/axial_plane_t*.txt`, node exports of the scalar on the plane x = 0 every 25 s.

```bash
python porous_plug_3d/tracer_R1T/analyze_r1t.py     # writes R1T_metrics.json and R1T_curves.npz
```

Times are referred to the plug by removing the mean transit times of the inlet (4.4 s) and outlet
(0.4 s) conduits and are scaled by the mean residence time φπR²δ/Q = 76.2 s.

| Quantity | 3D (R1T) | Axisymmetric, advective |
|---|---|---|
| Plug residence time from the outlet integral | 75.6 s (−0.7 %) | 76.2 s |
| Outlet 50 / 90 / 99 % | 0.75 / 2.03 / 3.28 | 0.77 / 1.94 / 3.55 |
| Interior 50 / 90 / 99 % | 0.55 / 1.71 / 2.90 | 0.54 / 1.69 / 3.30 |
| Tail time constant | 0.40 | 0.70 |
| Interior / outlet residual in the tail | 0.40 | — |

The figure scripts in this folder (`plane_to_vtu.py`, `render_r1t_planes_paraview.py`,
`compose_figure7.py`, `compose_r1t_planes.py`) produce panels (c–e) of Fig. 7 and Fig. S9; see the
main README.
