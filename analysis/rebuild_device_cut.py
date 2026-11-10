"""Robust native tetrahedron z=0 cut, including plane-coincident mesh nodes.
Set BIOREACTOR_RAW to the hydrated archive root containing Manuscript/.
This replaces the legacy strict-sign cut which omitted plane-coincident cells.
No resampling or flow solution is performed; each polygon retains its cell scalar.
"""
from pathlib import Path
import os,sys,json,shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/native'))
from dev3d import read3d
from fluentbin import read_dat
raw=Path(os.environ.get('BIOREACTOR_RAW','bioreactor_raw'))
base=raw/'Manuscript/figures/2026-10-21_arial12/device_native/device_2M_Q15'
p,faces,names=read3d(base.with_suffix('.cas'));data=read_dat(base.with_suffix('.dat'))
fluid=next(z for z,(t,n) in names.items() if t=='fluid');vy=data[(112,fluid)][1]
pieces={};positive=set();tol=1e-12
for arr in faces.values():
 xyz=p[arr[:,:3]-1];zz=xyz[:,:,2]
 mask=(zz.min(1)<=tol)&(zz.max(1)>=-tol);xyz=xyz[mask];zz=zz[mask];adj=arr[mask,3:]-1
 coords=[np.where((abs(zz[:,i])<=tol)[:,None],xyz[:,i,:2],np.nan) for i in range(3)]
 for i,j in [(0,1),(1,2),(2,0)]:
  cross=((zz[:,i]>tol)&(zz[:,j]<-tol))|((zz[:,j]>tol)&(zz[:,i]<-tol))
  frac=np.divide(zz[:,i],zz[:,i]-zz[:,j],out=np.zeros(len(zz)),where=(zz[:,i]!=zz[:,j]))
  q=xyz[:,i]+frac[:,None]*(xyz[:,j]-xyz[:,i]);coords.append(np.where(cross[:,None],q[:,:2],np.nan))
 coords=np.stack(coords,1)
 for pp,neighbours,z in zip(coords,adj,zz):
  pp=pp[~np.isnan(pp[:,0])]
  if not len(pp):continue
  for cell in neighbours:
   if cell<0:continue
   pieces.setdefault(int(cell),[]).extend(pp)
   if z.max()>tol:positive.add(int(cell))
polys=[];values=[];owners=[]
for cell,points in pieces.items():
 if cell not in positive:continue # deterministic positive-side choice on coincident faces
 v=np.unique(np.round(points,12),axis=0)
 if len(v)<3:continue
 center=v.mean(0);v=v[np.argsort(np.arctan2(v[:,1]-center[1],v[:,0]-center[0]))]*1000
 area=.5*abs(np.dot(v[:,0],np.roll(v[:,1],-1))-np.dot(v[:,1],np.roll(v[:,0],-1)))
 if area<1e-12:continue
 polys.append(v);values.append(abs(vy[cell]));owners.append(cell)
target=ROOT/'data/fields/Device_cut.vtk'
archived=ROOT/'data/checks/legacy_strict_sign_Device_cut.vtk'
if target.exists() and not archived.exists():shutil.copy2(target,archived)
with target.open('w') as f:
 f.write('# vtk DataFile Version 3.0\nRobust original Fluent tetrahedron z=0 cut; mm and m/s\nASCII\nDATASET POLYDATA\n')
 vertices=np.concatenate(polys);f.write(f'POINTS {len(vertices)} double\n');np.savetxt(f,np.c_[vertices,np.zeros(len(vertices))],fmt='%.12g')
 f.write(f'POLYGONS {len(polys)} {len(vertices)+len(polys)}\n');off=0
 for poly in polys:f.write(str(len(poly))+' '+' '.join(map(str,range(off,off+len(poly))))+'\n');off+=len(poly)
 f.write(f'CELL_DATA {len(polys)}\nSCALARS axial_speed double 1\nLOOKUP_TABLE default\n');np.savetxt(f,values,fmt='%.12g')
from shapely.geometry import Polygon
from shapely import union_all
shape=union_all([Polygon(v) for v in polys]);total=sum(Polygon(v).area for v in polys)
expected=3.860*31.68+33.328*26.873+4.999*(60.33-58.553)
audit={'source_case':str(base.with_suffix('.cas')),'source_data':str(base.with_suffix('.dat')),
 'cut_cells':len(polys),'total_polygon_area_mm2':total,'union_area_mm2':shape.area,'overlap_mm2':total-shape.area,
 'as_meshed_rectangular_section_area_mm2':expected,'area_relative_difference':shape.area/expected-1,
 'bounds_mm':list(shape.bounds),'method':'include plane-coincident nodes; choose positive side on coplanar internal faces',
 'legacy_issue':'strict z sign test omitted cells at coincident nodes; legacy slice preserved in data/checks/'}
(ROOT/'data/checks/device_cut_verification.json').write_text(json.dumps(audit,indent=2));print(audit,flush=True)
assert abs(total-shape.area)<1e-6
assert abs(shape.area/expected-1)<.002
