# v3: letters only.
"""Assumption and fit sensitivity from existing data, not new measurements."""
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from figure_style import STYLE,WIDTH,CM
from run_labels import PAPER_LABEL
plt.rcParams.update(STYLE)
R=Path(__file__).resolve().parents[1];O=R/'figures'
H=json.loads((R/'headline_numbers.json').read_text());C=json.loads((R/'analysis/theory/curvi_test.json').read_text());T=json.loads((R/'analysis/tracer/tailfit.json').read_text())
fig,ax=plt.subplots(3,1,figsize=(WIDTH,18*CM),layout='constrained')
phi=np.linspace(.1,1,201)
dp=5/H['capacity']['slope_uniform_mLmin_per_kPa']['G3']
for f,c,ls,title in [(1,'#171717','-','Uniform faces'),(H['capacity']['F_finest'],'#bc5618','--','Port feeding')]:
 ax[0].plot(phi,dp*f*.27/phi,ls,color=c,label=title)
ax[0].axvline(.27,color='.5',ls=':',lw=.8)
ax[0].set(xlabel=r'Assumed porosity, $\phi$',ylabel='Plug pressure demand\nat 5 mL/min (kPa)',xlim=(.1,1),ylim=(0,90))
ax[0].legend(frameon=False,loc='upper right',ncols=2)
ax[0].set_title('(a)',loc='left',pad=10)
for sig,col,marker in [(3,'#171717','o'),(4,'#bc5618','D')]:
 x=[C[g][f'sig{sig}']['terms_6_12']['h1sq']/1.4 for g in ['G1','G2','G3','G4','G5']]
 y=[C[g][f'sig{sig}']['6_12']['theory_over_cfd'] for g in ['G1','G2','G3','G4','G5']]
 ax[1].plot(x,y,marker,linestyle='none',color=col,mfc='white' if sig==3 else col,label=rf'$\sigma_c={sig}$ µm')
 if sig==4:
  for i in [0,1]:ax[1].annotate(f'G{i+1}',(x[i],y[i]),xytext=(8,1),textcoords='offset points')
ax[1].axhline(1,color='.5',lw=.8,ls='--');ax[1].text(.02,.90,'G3–G5',fontsize=12)
ax[1].set(xlabel='Weighted mean squared wall slope',ylabel=r'$\Delta p_{\rm theory}/\Delta p_{\rm CFD}$',xlim=(0,.25),ylim=(.87,1.26))
ax[1].legend(frameon=False,loc='upper left',ncols=2)
ax[1].set_title('(b)',loc='left',pad=10)
for i,(run,v) in enumerate(T.items()):
 x=[v['tau_ms'][k]/v['t_ideal_ms'] for k in ['0.001_0.03','0.003_0.1']];y=v['ratio_V_out']
 ax[2].semilogx(x,[y,y],color='.6',lw=.9)
 ax[2].plot(x[0],y,'o',color='#171717',mfc='white');ax[2].plot(x[1],y,'D',color='#bc5618',ms=4)
 if run in ['E7','E8','E10','E14']:ax[2].annotate(PAPER_LABEL[run],(x[1],y),xytext=(5,4 if run!='E10' else -14),textcoords='offset points')
ax[2].plot([],[],'o',mfc='white',mec='#171717',label='Residual 0.001–0.03')
ax[2].plot([],[],'D',color='#bc5618',label='Residual 0.003–0.1')
ax[2].axhline(1,color='.5',lw=.8,ls='--');ax[2].axvline(1,color='.5',lw=.8,ls=':')
ax[2].set(xlabel=r'Exchange-to-flush ratio, $\Lambda_e$',ylabel=r'$t_{99,V}/t_{99,o}$',xlim=(.12,8),ylim=(.70,3.35))
ax[2].set_xticks([.2,.5,1,2,5], ['0.2','0.5','1','2','5'])
from matplotlib.ticker import NullFormatter
ax[2].xaxis.set_minor_formatter(NullFormatter())
ax[2].legend(frameon=False,loc='upper left')
ax[2].set_title('(c)',loc='left',pad=10)
for ext in ['pdf','svg','png']:fig.savefig(O/f'Fig_Model_Limitations.{ext}',dpi=220)
plt.close(fig)
