"""Read the project's ASCII legacy VTK arrays for geometry and scalar checks.
Field rendering is performed exclusively by ParaView in render_fields_paraview.py.
"""
import re
import numpy as np
def read(path):
 text=path.read_text()
 m=re.search(r'POINTS (\d+) \w+\s*\n',text);n=int(m[1]);p=np.fromstring(text[m.end():],sep=' ',count=3*n).reshape(n,3)
 m=re.search(r'(?:CELLS|POLYGONS) (\d+) (\d+)\s*\n',text);nc=int(m[1]);flat=np.fromstring(text[m.end():],sep=' ',count=int(m[2]),dtype=int)
 cells=[];i=0
 while i<len(flat):
  k=flat[i];cells.append(flat[i+1:i+1+k]);i+=k+1
 start=text.index('CELL_DATA');end=text.find('POINT_DATA',start);chunk=text[start:end if end>0 else None]
 scalars={}
 for m in re.finditer(r'SCALARS (\w+) \w+(?: \d+)?\s*\nLOOKUP_TABLE \w+\s*\n',chunk):
  scalars[m[1]]=np.fromstring(chunk[m.end():],sep=' ',count=nc)
 return p,cells,scalars
