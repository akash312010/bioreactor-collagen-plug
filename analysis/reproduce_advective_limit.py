"""Repeat the advective-limit contour/transit calculation (k3_transit.py).

Run in official ParaView pvpython; input is the same final G3 cell field.
"""
from pathlib import Path
import sys,runpy,json,hashlib
import numpy as np
R=Path(__file__).resolve().parents[1]
source=R/'data/fields/G3.vtk'
out=R/'data/checks/advective_limit_local.json'
sys.argv=['k3_transit.py',str(source),str(out),'4000']
ns=runpy.run_path(str(R/'analysis/theory/k3_transit.py'),run_name='__main__')
np.savez_compressed(R/'data/checks/advective_limit_curves.npz',T=ns['Tn'],weights=ns['w'],time=ns['taus'],internal=ns['Rv'])
d=json.loads(out.read_text());d['sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
d['scope']='Zero-diffusion streamline limit; near-wall internal t99 is not a converged quantitative result.'
out.write_text(json.dumps(d,indent=2))
