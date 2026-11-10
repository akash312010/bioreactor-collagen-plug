import sys, json, time, numpy as np
from tsolver import *
import os
U = os.environ.get("CFD_RUNS", "CFD_runs_new/")  # Fluent case/data NOT in this repository
B={"B1":U+"B1_v2_p35_fine/New60k_fine_mesh/New60k_fine_mesh/fine_mesh","B2":U+"B2_v2_p35_medium/medium_mesh/medium_mesh/collagen_p35_medium",
   "B3":U+"B3_v2_p35_coarse/Coarse_mesh/Coarse_mesh/collagen_p35_coarse","B4":U+"B4_v2_p35_fine_y1_25nm/fine_y1_25nm/fine_y1_25nm/collagen_p35_fine_y1_25nm",
   "D1":U+"D1/new/New Run_D1/collagen_p28_fine"}
S30=0.0469/0.00769
jobs=[]
# A: time-step study, p35 fine, E1 and E6 conditions, BDF2 (Fluent's scheme) and CN
for case,sc,T in (("E1",1.0,5.0),("E6",S30,6.5)):
  for N in (250,500,1000,2000,4000):
    jobs.append(("dt",case,"B1",sc,6.7e-10,N,T,"bdf2"))
  for N in (500,2000):
    jobs.append(("dt",case,"B1",sc,6.7e-10,N,T,"cn"))
# B: grid study of the tracer, N=2000 BDF2
for case,sc,T in (("E1",1.0,5.0),("E6",S30,6.5)):
  for g in ("B3","B2","B4"):
    jobs.append(("grid",case,g,sc,6.7e-10,2000,T,"bdf2"))
# C: diffusivity at 37 C
for case,sc,T in (("E1",1.0,5.0),("E6",S30,6.5)):
  jobs.append(("D37",case,"B1",sc,9.0e-10,2000,T,"bdf2"))
# D: p28 at 30 mL/min (E8, slow pocket), time step
for N in (1000,2000,4000):
  jobs.append(("dt",("E8"),"D1",S30,6.7e-10,N,10.0,"bdf2"))
part=int(sys.argv[1]); jobs=jobs[part::2]
res={}; cache={}
for j in jobs:
  kind,case,g,sc,D,N,T,sch=j; key=f"{kind}|{case}|{g}|D{D:.1e}|N{N}|{sch}"
  ck=(g,sc,D)
  if ck not in cache: cache.clear(); cache[ck]=build(B[g],D,sc)
  t0=time.time(); out,tid,cmin,cmax=run(cache[ck],N,T,sch)
  r=times(out); r.update(t_ideal=tid,cmin=cmin,cmax=cmax,vol_end=float(out[-1,3]),secs=time.time()-t0,cells=len(cache[ck]["A"]))
  res[key]=r; np.save("curve_"+key.replace("|","_")+".npy",out[::max(1,N//500)])
  print(key,{k:round(v,4) for k,v in r.items()},flush=True)
  json.dump(res,open(f"study_{part}.json","w"),indent=1)
