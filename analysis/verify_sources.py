"""Check the original device workbook, native snapshot integrals, and Darcy grids."""
from pathlib import Path
import sys,json,csv,importlib.util
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
from openpyxl import load_workbook
from vtk_ascii import read
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/checks'
rows=list(csv.DictReader(open(ROOT/'data/device/Part1_data_reconciled.csv')))
wb=load_workbook(ROOT/'data/device/Part1_dimensionless_corrected.xlsx',read_only=True,data_only=True)
ws=wb['Corrected_Dimensionless'];values=list(ws.values);matched=[]
for row in rows:
 mesh=row['mesh'];q=float(row['Q_ml_min']);dp=float(row['dP_Pa'])
 candidates=[r for r in values if len(r)>=3 and str(r[0])==mesh and isinstance(r[1],(float,int)) and float(r[1])==q]
 assert len(candidates)==1,(mesh,q,candidates)
 assert abs(float(candidates[0][2])-dp)<1e-10,(mesh,q)
 matched.append({'mesh':mesh,'Q_mLmin':q,'loss_Pa':dp})
wb.close();print('Workbook: all 18 pressure values agree with CSV.',flush=True)
pro=json.loads((ROOT/'data/fields/original_panel_provenance.json').read_text())['Fig_Tracer_Snapshots'];fields=[]
for v in pro:
 name=f"{v['run']}_t{round(v['time_over_ideal'])}";p,cells,sc=read(ROOT/'data/fields'/f'{name}.vtk')
 area=np.array([.5*abs(np.dot(p[c,0],np.roll(p[c,1],-1))-np.dot(p[c,1],np.roll(p[c,0],-1))) for c in cells])
 avg=float(np.average(sc['uds0'],weights=area));difference=avg-v['internal_average_from_report']
 assert abs(difference)<5e-9,(name,avg,difference)
 fields.append({'field':name,'cell_mean':avg,'report_mean':v['internal_average_from_report'],'absolute_difference':abs(difference),'time_s':v['time_s']})
 print(name,'mean',avg,'difference',difference,flush=True)
spec=importlib.util.spec_from_file_location('spread',ROOT/'analysis/porous/spread.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
saved=json.loads((ROOT/'analysis/porous/fcurve.json').read_text())['device']['F'];darcy=[]
for h,expected in zip([.1,.05,.025],saved):
 f,_=mod.solve(16.664,26.873,1.93,2.4995,h);assert abs(f-expected)<2e-8,(h,f,expected)
 darcy.append({'h_mm':h,'F_recomputed':f,'F_archived':expected});print('Darcy',h,f,flush=True)
f,_=mod.solve(16.664,26.873,16.664,16.664,.2);assert abs(f-1)<1e-8
result={'workbook_pressure_matches':matched,'snapshot_integral_checks':fields,'darcy_grids':darcy,'darcy_full_face_F':f,
 'workbook_caution':'Other dimensionless/geometry cells in the historical workbook are superseded; current calculations use as-meshed dimensions and mass-conserving outlet mean.'}
(OUT/'source_verification.json').write_text(json.dumps(result,indent=2));print('Source and Darcy checks passed.',flush=True)
