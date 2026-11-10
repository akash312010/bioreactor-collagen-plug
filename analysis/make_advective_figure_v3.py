# v3: letters only.
"""Numerical curves of the advective (zero-diffusion) limit."""
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from figure_style import STYLE,WIDTH,CM
plt.rcParams.update(STYLE)
R=Path(__file__).resolve().parents[1];O=R/'figures'
A=np.load(R/'data/checks/advective_limit_curves.npz');E=json.loads((R/'analysis/tracer/E_all.json').read_text());H=json.loads((R/'headline_numbers.json').read_text())['tracer']['runs']
T=A['T'];w=A['weights'];idx=np.argsort(T);cdf=np.cumsum(w[idx])/w.sum()
fig,ax=plt.subplots(3,1,figsize=(WIDTH,15*CM),layout='constrained')
ax[0].loglog(T[idx],np.maximum(1-cdf,1e-9),'--',color='.45',label='Advection only')
ax[1].loglog(A['time'],np.maximum(1-A['internal'],1e-9),'--',color='.45',label='Advection only')
for run,col in [('E1','#171717'),('E6','#bc5618')]:
 for a,name in zip(ax[:2],['outlet_massweighted.out','volume_average.out']):
  v=np.loadtxt(R/'data/tracer_reports'/run/name,skiprows=4)
  a.loglog(v[:,2]/(E[run]['t_ideal_ms']*.001),np.maximum(1-v[:,1],1e-9),color=col,label=f"G3, Pe = {H[run]['Pe']:.0f}")
for a,title,ylabel in zip(ax[:2],['(a)','(b)'],['Outlet residual','Internal residual']):
 a.set(xlabel=r'Time / $t_{\rm ideal}$',ylabel=ylabel,xlim=(.3,100),ylim=(.001,1));a.set_title(title,loc='left',pad=32)
 a.axhline(.01,color='.7',lw=.7,ls=':');a.legend(frameon=False,ncols=3,loc='lower left',bbox_to_anchor=(0,1.01),borderaxespad=0)
for geom,marker,col in [('G3','o','#171717'),('G5','^','#bc5618'),('G1','s','#bc5618')]:
 runs=[k for k in H if H[k]['geom']==geom]
 ax[2].semilogx([E[k]['Lambda'] for k in runs],[H[k]['ratio_V_out'] for k in runs],marker,linestyle='none',color=col,mfc='white' if geom=='G5' else col,label=geom)
ax[2].set(xlabel='Transverse-mixing number',ylabel=r'$t_{99,V}/t_{99,o}$',xlim=(.012,7),ylim=(.7,3.25))
ax[2].set_title('(c)',loc='left',pad=9)
ax[2].set_xticks([.02,.1,.5,1,5],['0.02','0.1','0.5','1','5']);ax[2].xaxis.set_minor_formatter(NullFormatter())
ax[2].legend(frameon=False,ncols=3,loc='upper right');ax[2].axhline(1,color='.6',ls='--',lw=.8)
for ext in ['pdf','svg','png']:fig.savefig(O/f'Fig_S6_Advective_Limit.{ext}',dpi=220)
plt.close(fig)
