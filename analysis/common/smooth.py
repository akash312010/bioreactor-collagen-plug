import numpy as np, csv
from scipy.ndimage import gaussian_filter1d
def load(csvpath):
    a=np.array([[float(v) for v in r] for r in list(csv.reader(open(csvpath)))[1:]])*1e3   # um
    return np.c_[a[:,0],a[:,2]], np.c_[a[:,0],a[:,3]]        # wallA (Yl), wallB (Yr), each ordered X=0 -> X=L
def resample(P, ds):
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(P,axis=0),axis=1))]; n=max(int(np.ceil(s[-1]/ds)),2)
    t=np.linspace(0,s[-1],n+1); return np.c_[np.interp(t,s,P[:,0]),np.interp(t,s,P[:,1])]
def smooth_wall(P, sigma_um, ds=0.02, out_ds=0.1, blend=3.0):
    Q=resample(P, ds)
    if sigma_um>0:
        S=np.c_[gaussian_filter1d(Q[:,0],sigma_um/ds,mode="nearest"), gaussian_filter1d(Q[:,1],sigma_um/ds,mode="nearest")]
        # keep the exact end points, blend over 3 sigma so the end lines stay where they were
        n=len(S); k=int(blend*sigma_um/ds); wgt=np.ones(n)
        if k>0:
            r=np.linspace(0,1,k); wgt[:k]=r; wgt[-k:]=r[::-1]
        S=Q+(S-Q)*wgt[:,None]
    else: S=Q
    return resample(S, out_ds)
def min_radius(P):
    Q=resample(P,0.02); d1=np.gradient(Q,axis=0); d2=np.gradient(d1,axis=0)
    k=np.abs(d1[:,0]*d2[:,1]-d1[:,1]*d2[:,0])/np.maximum(np.linalg.norm(d1,axis=1)**3,1e-30)
    return 1/np.max(k[5:-5])
def poly_area(A,B):
    P=np.r_[A,B[::-1]]; x,y=P[:,0],P[:,1]; return 0.5*abs(np.dot(x,np.roll(y,-1))-np.dot(y,np.roll(x,-1)))
