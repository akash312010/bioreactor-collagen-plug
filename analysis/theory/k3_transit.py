"""K3(ii): advective (Pe -> infinity) transit-time distribution on the supplied G3.vtk.
Streamlines = iso-lines of the stream function psi (point data, normalised 0..1 between the walls, from solver
face mass fluxes). For a level psi, T(psi) = sum over contour segments of (segment length)/|u|_cell, using the
cell-centred speed of the cell that generated the segment (vtkContourFilter passes input cell data to its output).
Outlet flux-weighted distribution: psi is uniformly distributed in flux, so quantiles of T over psi are the outlet
flux-weighted t50/t90/t99. Volume replaced (advection only): R(t) = (1/A) int min(t, T) q dpsi.
Recirculating (never-reached) area: A_rec = A - q int T dpsi.   Run: pvpython k3_transit.py G3.vtk out.json"""
import sys, json, numpy as np
import vtkmodules.all as vtk
from vtkmodules.util.numpy_support import vtk_to_numpy as v2n
fn, outj = sys.argv[1], sys.argv[2]
r = vtk.vtkUnstructuredGridReader(); r.SetFileName(fn); r.ReadAllScalarsOn(); r.ReadAllVectorsOn(); r.ReadAllFieldsOn(); r.Update()
ug = r.GetOutput(); fd = ug.GetFieldData()
q = fd.GetArray("q_m2s").GetTuple1(0); A = fd.GetArray("A_um2").GetTuple1(0)*1e-12; tid = A/q
ug.GetPointData().SetActiveScalars("psi")
psi_rng = ug.GetPointData().GetArray("psi").GetRange()
# levels clustered toward both walls (Chebyshev-like), plus uniform interior
N = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
u = (np.arange(N)+0.5)/N
levels = 0.5 - 0.5*np.cos(np.pi*u)                      # dense near 0 and 1
w = np.pi/2*np.sin(np.pi*u)/N                            # d(psi) weights for this mapping (sum = 1)
T = np.full(N, np.nan); nseg = np.zeros(N, int); nlines = np.zeros(N, int)
cf = vtk.vtkContourFilter(); cf.SetInputData(ug); cf.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "psi")
cf.ComputeScalarsOff()
for i, lv in enumerate(levels):
    cf.SetValue(0, lv); cf.Update(); o = cf.GetOutput()
    if o.GetNumberOfCells() == 0: continue
    pts = v2n(o.GetPoints().GetData()); sp = v2n(o.GetCellData().GetArray("speed"))
    lines = o.GetLines(); lines.InitTraversal(); idl = vtk.vtkIdList(); tt = 0.0; k = 0
    while lines.GetNextCell(idl):
        ids = [idl.GetId(j) for j in range(idl.GetNumberOfIds())]
        for a, b in zip(ids[:-1], ids[1:]):
            L = np.linalg.norm(pts[a]-pts[b])*1e-6; tt += L/max(sp[k], 1e-30)
        k += 1
    T[i] = tt; nseg[i] = o.GetNumberOfCells()
ok = np.isfinite(T)
Tn = T/tid
order = np.argsort(Tn[ok]); cdf = np.cumsum(w[ok][order])/np.sum(w[ok])
quant = lambda p: float(Tn[ok][order][np.searchsorted(cdf, p)])
mean_flux_T = float(np.sum(w[ok]*Tn[ok]))           # = (A - A_rec)/A in units of t_ideal
A_rec_frac = 1 - mean_flux_T
def R(tau): return float(np.sum(w[ok]*np.minimum(tau, Tn[ok])))
taus = np.logspace(-1, 3, 2000); Rv = np.array([R(t) for t in taus])
def tvol(p):
    target = p*(1-A_rec_frac)                         # fraction of the through-flow (reachable) area
    i = np.argmax(Rv >= target); return float(taus[i]) if Rv.max() >= target else None
res = dict(input=fn, levels=N, psi_range=list(psi_rng), q_m2s=q, A_m2=A, t_ideal_s=tid,
           first_arrival=float(Tn[ok].min()), outlet_t50=quant(0.5), outlet_t90=quant(0.9), outlet_t99=quant(0.99),
           recirculating_area_frac=A_rec_frac, volume_t90_of_reachable=tvol(0.9), volume_t99_of_reachable=tvol(0.99),
           volume_t99_total=(float(taus[np.argmax(Rv >= 0.99)]) if Rv.max() >= 0.99 else None), volume_t90_total=(float(taus[np.argmax(Rv >= 0.90)]) if Rv.max() >= 0.90 else None),
           R_inf=float(Rv[-1]), levels_failed=int((~ok).sum()))
json.dump(res, open(outj, "w"), indent=1); print(json.dumps(res, indent=1))
