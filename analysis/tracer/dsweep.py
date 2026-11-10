import sys, json, numpy as np
from tsolver import *
import os
U = os.environ.get("CFD_RUNS", "CFD_runs_new/")  # Fluent case/data NOT in this repository
base=U+"D1/new/New Run_D1/collagen_p28_fine"; S30=0.0469/0.00769
f=float(sys.argv[1]); D=6.7e-10*f
S=build(base,D,S30); out,tid,cmin,cmax=run(S,1000,12.0,"bdf2")
t=out[:,0]*tid; c=out[:,3]; res={"f":f,"D":D,"t_ideal":tid}
for lo,hi in ((4,8),(6,10),(8,12)):
    k=(out[:,0]>lo)&(out[:,0]<hi)&(1-c>1e-7); a,b=np.polyfit(t[k],np.log(1-c[k]),1); res["tau_ms_%d_%d"%(lo,hi)]=-1e3/a
r=times(out); res.update(vol_t99=r["vol_t99"], mass_t99=r["mass_t99"])
print(json.dumps(res)); json.dump(res,open("dsweep_%g.json"%f,"w"))
