"""17 cm, Arial 12 pt geometry plate from measured meshes and source pixels.
v3 (1 Nov 2026): panel letters only, numeric scale bars, numeric channel length."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection,PolyCollection
from matplotlib.patches import Polygon,Rectangle,ConnectionPatch
from PIL import Image
from vtk_ascii import read
from figure_style import WIDTH, CM, STYLE
plt.rcParams.update(STYLE)
R=Path(__file__).resolve().parents[1];O=R/'figures'
p,cells,_=read(R/'data/fields/G3.vtk');p=p[:,:2]
edges=np.concatenate([np.column_stack([c,np.roll(c,-1)]) for c in cells])
unique,count=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
walls=p[unique[count==1]];polys=[p[c] for c in cells]
P,C,_=read(R/'data/fields/Device_cut.vtk');P=P[:,:2]
dev_edges=np.concatenate([np.stack([P[c],P[np.roll(c,-1)]],axis=1) for c in C])
fig,axs=plt.subplots(3,2,figsize=(WIDTH,18*CM))
fig.subplots_adjust(left=.04,right=.98,bottom=.035,top=.96,wspace=.12,hspace=.22)
a,b,c,d,e,f=axs.ravel()
ri,r,ro=1.93,16.664,2.4995;roof,floor,bottom=-31.68,-58.553,-60.33
outline=np.array([[-ri,0],[ri,0],[ri,roof],[r,roof],[r,floor],[ro,floor],[ro,bottom],[-ro,bottom],[-ro,floor],[-r,floor],[-r,roof],[-ri,roof]])
a.add_patch(Polygon(outline,fc='#f5f5f5',ec='black',lw=.8))
a.add_patch(Rectangle((-r,floor),2*r,roof-floor,fc='#e4c6ad',ec='none'))
a.add_patch(Polygon(outline,fill=False,ec='black',lw=.8))
a.annotate(r'$d_i$',xy=(ri,-3),xytext=(12,4),arrowprops={'arrowstyle':'-','lw':.7},va='center')
a.annotate('',xy=(-r,-43),xytext=(r,-43),arrowprops={'arrowstyle':'<->','lw':.7})
a.text(0,-39,r'$2R$',ha='center')
a.annotate('',xy=(-20,roof),xytext=(-20,floor),arrowprops={'arrowstyle':'|-|','lw':.7})
a.text(-24,-45,r'$\delta$',ha='center')
a.text(6,-16,r'$L_i$');a.annotate(r'$d_o, L_o$',xy=(ro,-59.5),xytext=(13,-63),arrowprops={'arrowstyle':'-','lw':.7},va='center')
b.add_collection(LineCollection(dev_edges,colors='.4',linewidths=.12))
b.add_patch(Polygon(outline,fill=False,ec='black',lw=.7))
def bar(ax,x,y,length,label,below=False):
 ax.plot([x,x+length],[y,y],color='black',lw=1.6,solid_capstyle='butt')
 ax.annotate(label,(x+length/2,y),xytext=(0,-3 if below else 3),textcoords='offset points',ha='center',va='top' if below else 'bottom',fontsize=12)
for ax in [a,b]:
 ax.set(xlim=(-29,29),ylim=(-74,8),aspect='equal');ax.axis('off')
 bar(ax,-25 if ax is a else -12,-64.5,10,'10 mm',below=True)
c.imshow(Image.open(R/'data/geometry/collagen_SHG_Niestrawska2022_Fig4.png'),interpolation='none');c.axis('off')
# Rotation is rigid: the narrow channel has the same height as the adjacent
# micrograph, without stretching its physical width.
theta=np.deg2rad(-60);rot=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
pp=p@rot.T;ww=walls@rot.T
d.add_collection(PolyCollection([pp[cc] for cc in cells],facecolor='#eddbc9',edgecolor='none'))
d.add_collection(LineCollection(ww,colors='black',linewidths=.7))
lo,hi=pp.min(0),pp.max(0)
d.set(xlim=(lo[0]-27,hi[0]+35),ylim=(lo[1]-15,hi[1]+12),aspect='equal');d.axis('off')
ex,ey=rot[:,0],rot[:,1]
origin=np.array([lo[0]-22,hi[1]-20])
for vector,name in [(ex,'X'),(ey,'Y')]:
 end=origin+14*vector
 d.annotate('',xy=end,xytext=origin,arrowprops={'arrowstyle':'->','lw':.8})
 d.text(*(end+4.5*vector),rf'${name}$',ha='center',va='center')
pin=p[np.isclose(p[:,0],p[:,0].min())].mean(0)@rot.T
pout=p[np.isclose(p[:,0],p[:,0].max())].mean(0)@rot.T
d.annotate('',xy=pin-9*ex,xytext=pin+5*ex,arrowprops={'arrowstyle':'->','lw':1,'color':'#bc5618'})
d.annotate('',xy=pout-5*ex,xytext=pout+9*ex,arrowprops={'arrowstyle':'->','lw':1,'color':'#bc5618'})
start=pin+16*ey;end=start+(p[:,0].max()-p[:,0].min())*ex
d.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'|-|','lw':.7})
d.text(*((start+end)/2+9*ey),f'{p[:,0].max()-p[:,0].min():.1f} µm',ha='center',va='center',rotation=np.rad2deg(np.arctan2(ex[1],ex[0])))
bar(d,lo[0]+2,lo[1]-6,10,'10 µm')
for ax,win in [(e,(41,53,-11.2,-3.8)),(f,(45.55,46.55,-9.66,-9.06))]:
 x0,x1,y0,y1=win
 mask=(p[unique][:,:,0].max(1)>=x0)&(p[unique][:,:,0].min(1)<=x1)&(p[unique][:,:,1].max(1)>=y0)&(p[unique][:,:,1].min(1)<=y1)
 ax.add_collection(LineCollection(p[unique[mask]],colors='.45',linewidths=.35))
 ax.add_collection(LineCollection(walls,colors='black',linewidths=.8))
 ax.set(xlim=(x0,x1),ylim=(y0,y1),aspect='equal');ax.set_xticks([]);ax.set_yticks([])
bar(e,41.6,-10.75,2,'2 µm');bar(f,45.6,-9.64,.2,'0.2 µm')
e.add_patch(Rectangle((45.55,-9.66),1,.6,fill=False,ec='#bc5618',lw=.9))
fig.add_artist(ConnectionPatch((46.55,-9.06),(45.55,-9.06),'data','data',axesA=e,axesB=f,color='#bc5618',lw=.7))
for ax,letter in zip(axs.ravel(),'abcdef'):
 ax.set_title(f'({letter})',loc='left',fontsize=12,pad=8)
for ext in ['pdf','svg','png']:fig.savefig(O/f'Fig_Geometry_Across_Scales.{ext}',dpi=220)
plt.close(fig)
