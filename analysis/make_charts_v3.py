# v3: panel titles removed; parameters stated in captions.
"""Rebuild numerical charts as editable PDF and SVG, with white backgrounds.

No CFD field is rendered by this script. Field panels are produced separately
with ParaView and assembled by compose_figures.py.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import json, re
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures'
OUT.mkdir(parents=True,exist_ok=True)
J=lambda p:json.loads((ROOT/p).read_text())
H=J('headline_numbers.json')
ORANGE='#bc5618'; BLACK='#171717'; GREY='#777777'
from figure_style import WIDTH, CM, STYLE
from run_labels import PAPER_LABEL
plt.rcParams.update(STYLE)

def save(fig,name):
    fig.savefig(OUT/(name+'.pdf'),bbox_inches=None)
    fig.savefig(OUT/(name+'.svg'),bbox_inches=None)
    fig.savefig(OUT/(name+'.png'),bbox_inches=None,dpi=200)
    plt.close(fig)
def label(ax,s): ax.text(-.13,1.05,s,transform=ax.transAxes,ha='left',va='bottom',fontsize=12)
def legend(ax,**kw): return ax.legend(frameon=False,handlelength=1.7,**kw)

G=[f'G{i}' for i in range(1,6)]; x=np.arange(5)
C=J('analysis/theory/curvi_test.json'); cuts=['4_8','6_12','8_16','10_20']
fig,axs=plt.subplots(1,2,figsize=(WIDTH,8*CM),layout='constrained')
for key,name,col,marker,shift in [
 ('leading','Leading',GREY,'o',-.18),
 ('no_curv_over_cfd','Terms omitted',BLACK,'s',0),
 ('theory_over_cfd','Second order',ORANGE,'D',.18)]:
    val=[1/C[g]['sig4']['6_12']['cfd_over_lub'] if key=='leading' else C[g]['sig4']['6_12'][key] for g in G]
    axs[0].plot(x+shift,val,linestyle='none',marker=marker,color=col,
                markerfacecolor='white' if key!='theory_over_cfd' else col,label=name)
axs[0].axhline(1,color=GREY,lw=.8,ls='--');axs[0].set_xticks(x,G)
axs[0].set(ylabel=r'$\Delta p_{\rm model}/\Delta p_{\rm CFD}$',ylim=(.65,1.3))
legend(axs[0],loc='upper right',fontsize=12);label(axs[0],'(a)')
for key,name,col,shift in [('theory_over_cfd','Second order',ORANGE,-.08),('no_curv_over_cfd','Terms omitted',BLACK,.08)]:
    vals=np.array([[C[g][f'sig{s}'][cut][key] for s in (3,4) for cut in cuts] for g in G])
    for i,v in enumerate(vals):
        axs[1].plot(np.full(len(v),i+shift)+np.linspace(-.025,.025,len(v)),v,'.',color=col,ms=3.3,alpha=.75)
        axs[1].vlines(i+shift,v.min(),v.max(),color=col,lw=1.3)
    axs[1].plot([],[],'.-',color=col,label=name)
axs[1].axhline(1,color=GREY,lw=.8,ls='--');axs[1].set_xticks(x,G)
axs[1].set(ylabel=r'$\Delta p_{\rm model}/\Delta p_{\rm CFD}$',xlabel='Channel',ylim=(.65,1.3))
legend(axs[1],loc='upper right',fontsize=12);label(axs[1],'(b)')
save(fig,'Fig_Resistance_Theory')

fig,axs=plt.subplots(1,2,figsize=(WIDTH,8*CM),layout='constrained')
P=np.geomspace(.5,100,200);s=np.array(list(H['capacity']['slope_uniform_mLmin_per_kPa'].values()))
F=H['capacity']['F_finest'];s0=H['capacity']['slope_uniform_mLmin_per_kPa']['G3']
for factor,col,style,title in [(1,BLACK,'-','Uniform faces'),(F,ORANGE,'--','Port feeding')]:
    axs[0].fill_between(P,P*s.min()/factor,P*s.max()/factor,color=col,alpha=.1,lw=0)
    axs[0].loglog(P,P*s0/factor,style,color=col,label=title)
axs[0].set(xlabel=r'$\Delta P_c$ (kPa)',ylabel=r'$Q_{\max}$ (mL/min)',xlim=(.5,100),ylim=(.03,200))
legend(axs[0],loc='upper left',fontsize=12);label(axs[0],'(a)')
fc=J('analysis/porous/fcurve.json');a=np.geomspace(.03,1,250);R,Ht=16.664,26.873
axs[1].loglog(a,1+np.pi*R/(2*Ht)*(1-a)**1.5/a,color=ORANGE,label='Correlation')
axs[1].loglog(a,1+np.pi*R/(2*Ht)/a,'--',color=GREY,label='Half space')
axs[1].loglog(fc['a_over_R'],fc['F'],'o',mfc='white',mec=BLACK,label='Darcy model')
axs[1].plot((fc['device']['a1']+fc['device']['a2'])/(2*R),F,'*',ms=10,color=ORANGE,mec=BLACK,mew=.5,label='Device ports')
axs[1].set(xlabel=r'$a/R$',ylabel=r'$F$',xlim=(.025,1.1),ylim=(.9,45))
legend(axs[1],loc='upper right',fontsize=12);label(axs[1],'(b)')
save(fig,'Fig_Operating_Map')

runs=list(H['tracer']['runs']);tr=H['tracer']['runs'];tail=J('analysis/tracer/tailfit.json')
fig,ax=plt.subplots(figsize=(WIDTH,9*CM),layout='constrained')
y=np.arange(len(runs))[::-1]
for yy,r in zip(y,runs):
    d=tr[r];xx=[d['t99_out'],d['t99_V']]
    ax.plot(xx,[yy,yy],color='#aaaaaa',lw=.8,zorder=1)
    ax.plot(xx[0],yy,'o',mfc='white',mec=BLACK,zorder=3)
    ax.plot(xx[1],yy,'D',mfc=ORANGE,mec=BLACK,mew=.55,zorder=3)
    if r=='E8':
        for value,dy in [(xx[0],.35),(xx[1],.35)]:ax.text(value,yy+dy,f'{value:.2f}',ha='center',fontsize=12)
ax.plot([],[],'o',mfc='white',mec=BLACK,label='Outlet')
ax.plot([],[],'D',mfc=ORANGE,mec=BLACK,mew=.5,label='Interior')
ax.set_yticks(y,[f"{tr[r]['geom']} | {tr[r]['Pe']:.1f}" for r in runs])
ax.set(xlabel=r'99% replacement time / $t_{\rm ideal}$',ylabel=r'Channel | Péclet number',xlim=(0,9.2),ylim=(-.65,10.7))
legend(ax,loc='upper right',ncols=2)
save(fig,'Fig_Replacement_Paired')

fig,axs=plt.subplots(2,2,figsize=(WIDTH,16*CM),sharex=True,sharey=True,layout='constrained')
for ax,sig,letter in zip(axs.ravel(),(1,2,3,4),'abcd'):
    for i,g in enumerate(G):
        vv=[C[g][f'sig{sig}'][cut]['theory_over_cfd'] for cut in cuts]
        ax.plot(i+np.linspace(-.12,.12,4),vv,'o',color=ORANGE,mec=BLACK,mew=.4,ms=3.6)
    ax.axhline(1,color=GREY,lw=.8,ls='--');ax.set_xticks(x,G)
    ax.set_ylim(.85,1.62)
    ax.set_ylabel(r'$\Delta p_{\rm theory}/\Delta p_{\rm CFD}$');label(ax,f'({letter})')
save(fig,'Fig_S3_Theory_Sensitivity')

fig,axs=plt.subplots(2,1,figsize=(WIDTH,16*CM),layout='constrained')
for run,col,style in [('E7',BLACK,'-'),('E8',ORANGE,'-')]:
    v=np.loadtxt(ROOT/'data'/'tracer_reports'/run/'volume_average.out',skiprows=4)
    axs[0].semilogy(v[:,2]*1000,np.maximum(1-v[:,1],1e-12),style,color=col,label=rf"G1, $\mathrm{{Pe}}={tr[run]['Pe']:.0f}$")
axs[0].axhline(.01,color=GREY,lw=.75,ls=':');axs[0].set(xlabel='Time (ms)',ylabel='Internal residual',ylim=(.001,1),xlim=(0,145));legend(axs[0]);label(axs[0],'(a)')
factors=[.5,1,2];tau=[J(f'analysis/tracer/dsweep_{f:g}.json')['tau_ms_6_10'] for f in factors]
axs[1].loglog(factors,tau,'o-',color=ORANGE,mec=BLACK,mew=.5)
axs[1].set(xlabel=r'Diffusivity / $D_0$',ylabel='Exchange time (ms)',xlim=(.4,2.5),ylim=(9,50));axs[1].set_xticks(factors,['0.5','1','2']);axs[1].set_yticks([10,20,40],['10','20','40']);label(axs[1],'(b)')
axs[1].xaxis.set_minor_formatter(NullFormatter());axs[1].yaxis.set_minor_formatter(NullFormatter())
save(fig,'Fig_S4_Exchange_Tails')

fig,ax=plt.subplots(figsize=(WIDTH,10*CM),layout='constrained')
for rr,col,marker,name in [(['E11','E12'],ORANGE,'D','Darcy interstitial mapping'),(['E1','E6'],BLACK,'o','Conduit peak mapping')]:
    ax.loglog([tr[r]['Pe'] for r in rr],[tr[r]['t_ideal_ms'] for r in rr],marker+'-',color=col,mec=BLACK,mew=.5,label=name)
    for r in rr:ax.annotate(PAPER_LABEL[r],(tr[r]['Pe'],tr[r]['t_ideal_ms']),xytext=(4,6),textcoords='offset points',fontsize=12)
ax.set(xlabel=r'G3 Péclet number',ylabel='Ideal flush time (ms)',xlim=(2.5,400),ylim=(2.5,400));legend(ax,loc='upper right',fontsize=12)
save(fig,'Fig_S5_Coupling')

fig,axs=plt.subplots(2,1,figsize=(WIDTH,16*CM),layout='constrained')
study=J('analysis/tracer/study_all.json')
for run,col,mark in [('E1',BLACK,'o'),('E6',ORANGE,'s'),('E8',GREY,'^')]:
    pts=[]
    for key,d in study.items():
        fields=key.split('|')
        if fields[0]=='dt' and fields[1]==run and fields[-1]=='bdf2':
            n=int(next(v[1:] for v in fields if re.match(r'^N\d+$',v)))
            pts.append((n,d['mass_t99']))
    pts=sorted(pts)
    axs[0].semilogx(*np.array(pts).T,mark+'-',color=col,mfc='white' if run=='E8' else col,label=PAPER_LABEL[run])
axs[0].set(xlabel='Steps per ideal flush time',ylabel=r'Outlet $t_{99,o}/t_{\rm ideal}$');legend(axs[0]);label(axs[0],'(a)')
axs[0].set_xticks([250,1000,4000],['250','1000','4000']);axs[0].xaxis.set_minor_formatter(NullFormatter())
co=J('analysis/tracer/courant.json');labs=list(co)
for key,name,col,marker in [('Co_mean_area','Area mean',BLACK,'o'),('Co_p99','99th percentile',ORANGE,'s'),('Co_max','Maximum',GREY,'^')]:
    axs[1].semilogy(np.arange(len(labs)),[co[l][key] for l in labs],marker+'-',color=col,ms=3.5,label=name)
axs[1].set_xticks(np.arange(len(labs)),labs);axs[1].set(xlabel='Frozen-flow mesh',ylabel='Cell Courant number');legend(axs[1],fontsize=12);label(axs[1],'(b)')
save(fig,'Fig_S2_Tracer_Verification')
pressure=J('data/verification/mesh_study_B.json');ss=pressure['stations']
fig,axs=plt.subplots(2,1,figsize=(WIDTH,16*CM),layout='constrained')
for row,col,mark in zip(ss,[BLACK,ORANGE,BLACK,ORANGE],['o','s','^','D']):
 axs[0].plot([21295,37045,52059],[row[k] for k in ('coarse','medium','fine')],mark+'-',color=col,mfc='white' if mark in ['^','D'] else col,ms=3.5,label=rf"$X={row['X_um']}\,\mu$m")
axs[0].set(xlabel='Cell count',ylabel='Interior pressure (Pa)');legend(axs[0],fontsize=12,ncols=2);label(axs[0],'(a)')
obs=[r['GCI_fine_pct'] for r in ss]
fixed=[1.25*abs(r['fine']-r['medium'])/abs(r['fine'])/(pressure['r21']**2-1)*100 for r in ss]
axs[1].plot(np.arange(4)-.07,obs,'o',color=ORANGE,label='Observed order');axs[1].plot(np.arange(4)+.07,fixed,'s',color=BLACK,mfc='white',label=r'Fixed $p=2$')
axs[1].set_xticks(range(4),[str(r['X_um']) for r in ss]);axs[1].set(xlabel=r'Station $X$ ($\mu$m)',ylabel='Fine-grid GCI (%)',ylim=(0,.4));legend(axs[1],fontsize=12);label(axs[1],'(b)')
save(fig,'Fig_S1_Pressure_Verification')

from vtk_ascii import read as read_vtk
p,cells,scalars=read_vtk(ROOT/'data/fields/Device_cut.vtk')
centers=np.array([p[c].mean(0) for c in cells]);m=(centers[:,1]>-60.2)&(centers[:,1]<-59.8)
fig,ax=plt.subplots(figsize=(8*CM,7.5*CM),layout='constrained');ax.plot(centers[m,0],scalars['axial_speed'][m],'.',color=BLACK,ms=1.7)
ax.set(xlabel='Outlet coordinate (mm)',ylabel='Axial speed magnitude\n(m/s)',xlim=(-2.6,2.6),ylim=(0,.046))
save(fig,'Device_Profile_Chart')
rows=list(csv.DictReader(open(ROOT/'data/device/Part1_data_reconciled.csv')))
fig,ax=plt.subplots(figsize=(WIDTH,6*CM),layout='constrained')
for mesh,col,marker in [('0.5M',GREY,'^'),('1M',BLACK,'o'),('2M',ORANGE,'s')]:
 rr=[r for r in rows if r['mesh']==mesh];ax.plot([float(r['Q_ml_min']) for r in rr],[float(r['dP_Pa']) for r in rr],marker+'-',color=col,mfc='white' if mesh=='0.5M' else col,ms=3.5,label=mesh)
ax.set(xlabel=r'Flow rate (mL min$^{-1}$)',ylabel='Device pressure\ndifference (Pa)',xlim=(3,32),ylim=(0,4.7));legend(ax,loc='upper left',fontsize=12)
save(fig,'Device_Loss_Chart')
print('Wrote eight numerical figures and device chart panels as PDF, SVG and PNG to',OUT)
