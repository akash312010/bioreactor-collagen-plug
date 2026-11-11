# Pressure demand and solute replacement in a collagen-plugged bioreactor

Data and analysis code for the paper

> M. M. H. Akash, M. Yeasin, K. Van Horsen, S. Mallik and S. Basu,
> "Pressure demand and solute replacement in a collagen-plugged bioreactor",
> submitted to *Physical Review Fluids* (2026).

The repository contains the reduced simulation records, the pore geometry, the 3D
porous-plug meshes and reports, and the scripts that turn them into every number,
table and chart in the paper. One command reproduces the headline numbers without
ANSYS Fluent.

## Quick start

```bash
python -m pip install -r requirements.txt
python reproduce_headlines.py                       # prints the headline numbers, writes headline_numbers.json
python porous_plug_3d/plug3d_postprocess.py         # 3D porous-plug check (Table IV)
python analysis/tracer/solute_balance.py            # solute balance of the tracer runs (Eq. 8, Table S5)
python analysis/tracer/outlet_reconstruction.py     # interior replacement from the outlet history (Eq. 19, Table S6)
python analysis/porous/plug_transport.py            # plug-scale speeds and transit times, three grids
python porous_plug_3d/tracer_R1T/analyze_r1t.py     # 3D tracer run R1T compared with the plug-scale curves
```

Python 3 (tested with 3.13). Field figures additionally need ParaView 6 (`pvpython`).

## The five calculations and where they are

| Calculation (name in the paper) | Dimension | What it gives | Files |
|---|---|---|---|
| Open-device simulations (fluid only, no plug) | 3D, steady laminar | pressure used by the housing; outlet speed | `data/device/`, `data/fields/Device_cut.vtk` |
| Pore simulations, channels G1–G5 | 2D, steady creeping flow | pore resistance C_CFD; permeability k | `data/geometry/`, `data/fields/G1.vtk`, `G3.vtk`, `G5.vtk`, `analysis/theory/` |
| Tracer simulations, runs E1–E12 (paper labels, see below) | 2D, time-dependent scalar on the frozen pore velocity | outlet and volume replacement times t99,o, t99,V; τe; solute balance | `data/tracer_reports/`, `data/fields/E*.vtk`, `analysis/tracer/` |
| Axisymmetric Darcy calculation | axisymmetric | feeding factor F = 7.12; plug-scale speed and transit-time distributions | `analysis/porous/` |
| 3D porous-plug simulation, runs R1–R3 | 3D, steady Darcy–Stokes | F3D and plug pressure | `porous_plug_3d/` |
| 3D tracer simulation, run R1T | 3D, time-dependent scalar on the frozen R1 flow | plug-scale outlet and interior replacement; axial-plane fields | `porous_plug_3d/tracer_R1T/` |

The open-device, pore, tracer, porous-plug and 3D tracer simulations were run in ANSYS
Fluent 2025 R2. The axisymmetric Darcy calculation (`analysis/porous/spread.py`) and the
independent tracer solver used for the time-step and grid study
(`analysis/tracer/tsolver.py`) are written in Python and included here.

## Contents

```
reproduce_headlines.py      recomputes the headline numbers from the files below
headline_numbers.json       its output
requirements.txt
rebuild_figures.ps1         rebuilds the chart figures (Windows PowerShell)

data/
  device/                   open-device pressure and speed records (1M reference series, original workbook)
  geometry/                 traced pore walls of G1-G5 (CSV, mm), channel properties, source image (SOURCE.md)
  fields/                   VTK fields used in the field figures: device cut, G1/G3/G5 pore flow,
                            tracer snapshots of E8, E11, E12; provenance of each panel
  tracer_reports/E*/        Fluent report series per tracer run (folder names = Fluent runs): outlet
                            mass-weighted and volume-averaged concentration at every time step
  verification/             grid convergence of the G3 pore pressure (meshes B1-B4, GCI)
  checks/                   reproduction checks: signed curvilinear metric, theory reproduction,
                            advective-limit transit-time curves

analysis/
  common/                   Fluent binary case/data reader; wall smoothing
  native/                   Fluent-to-VTK conversion and field plotting
  porous/                   axisymmetric Darcy solver (spread.py) and F curve; plug-scale stream
                            function, speed distribution and transit times (plug_transport.py) and
                            the plug-scale figure (make_plug_transport_figure.py)
  theory/                   lubrication and curvilinear second-order theory vs the pore simulations;
                            k3_transit.py, advective (zero-diffusion) transit times;
                            image_dark_fraction.py, dark-area fraction of the source image (porosity sensitivity)
  tracer/                   replacement times, tail fits, independent tracer solver and its studies;
                            solute_balance.py (Eq. 8, Table S5); outlet_reconstruction.py (Eq. 19, Table S6);
                            e14_vs_sweep.py (paper run E9, Fluent run E14, against the diffusivity sweep)
  make_*_v3.py, make_singlecol_figures.py, compose_figures_v3.py, render_fields_paraview.py
                            figure scripts (output goes to figures/)

porous_plug_3d/             3D porous-plug simulation: geometry (STEP), meshes M1 and M2,
                            mesh scripts, Fluent reports of R1-R3, post-processing
  tracer_R1T/               3D tracer run R1T: Fluent report series, axial-plane exports (16 times),
                            analysis, ParaView rendering and the composition of Fig. 7 and Fig. S9
```

## Tracer run labels

The paper numbers the twelve tracer runs E1–E12 in the order of its Table III. The folders in
`data/tracer_reports/` and the keys in the JSON files keep the names of the original Fluent runs.
`analysis/run_labels.py` holds the same map, and the figure and table scripts print the paper labels.

| Paper label | Channel | Pe | Fluent run (folder) |
|---|---|---|---|
| E1 | G3 | 3.6 | E11 |
| E2 | G3 | 21.5 | E12 |
| E3 | G3 | 42.2 | E1 |
| E4 | G3 | 76.8 | E2 |
| E5 | G3 | 115 | E3 |
| E6 | G3 | 160 | E4 |
| E7 | G3 | 257 | E6 |
| E8 | G1 | 22.6 | E7 |
| E9 | G1 | 69.0 | E14 |
| E10 | G1 | 138 | E8 |
| E11 | G5 | 119 | E9 |
| E12 | G5 | 729 | E10 |

The field files `data/fields/E8_*.vtk`, `E11_*.vtk` and `E12_*.vtk` are the paper's runs E10, E1 and E2.

## Figures of the plug-scale results

```bash
python analysis/porous/make_plug_transport_figure.py           # Fig. 7(a,b), with the R1T symbols
python porous_plug_3d/tracer_R1T/plane_to_vtu.py               # axial-plane exports -> .vtu
pvbatch porous_plug_3d/tracer_R1T/render_r1t_planes_paraview.py  # ParaView panels, Turbo map 0-1
python porous_plug_3d/tracer_R1T/compose_figure7.py            # Fig. 7 with panels (c-e)
python porous_plug_3d/tracer_R1T/compose_r1t_planes.py         # Fig. S9
```

Run `pvbatch` under a display (or `xvfb-run -a pvbatch ...` on Linux). The panel PDFs are
vector files from ParaView's GL2PS export; some PDF viewers warn about their graphics-state
entries, which can be removed by rewriting each panel with `pdftocairo -pdf in.pdf out.pdf`.

## What needs the original Fluent case and data files

`reproduce_headlines.py`, `porous_plug_3d/plug3d_postprocess.py` and the chart
scripts use only files in this repository. Scripts that re-read the Fluent solutions
(for example `analysis/theory/theory_local.py`, `analysis/tracer/tsolver.py`,
`analysis/rebuild_device_cut.py`) look for them in the folder named by the
environment variable `CFD_RUNS` (or `BIOREACTOR_RAW`). These case and data files
are large and are not in the repository; they are available from the corresponding
authors on reasonable request.

## Units and conventions

SI units unless a file states otherwise; the pore-wall CSVs and the 3D STEP geometry are in
millimetres, `cad_info.json` in micrometres, the Fluent meshes in metres. Fluid: ρ = 1012 kg m⁻³, μ = 9.64 × 10⁻⁴ Pa s.
Reference permeability (G3): k = 5.566 × 10⁻¹³ m². Pressures in the 3D reports are
gauge pressures relative to the outlet.

## Licence

Code: MIT licence (`LICENSE`). Data: CC BY 4.0 (`LICENSE-DATA.md`). The source image
in `data/geometry/` is reproduced from Niestrawska *et al.* (2022) under CC BY 4.0;
see `data/geometry/SOURCE.md`.

## Citation

Please cite the paper (see `CITATION.cff`). A DOI for this repository will be added
on publication.

## Contact

Mohammad Mehedi Hasan Akash, makash@eng.famu.fsu.edu
Saikat Basu, Saikat.Basu@sdstate.edu
