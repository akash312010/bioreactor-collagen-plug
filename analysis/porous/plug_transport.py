"""Plug-scale speed and transit-time distributions of the port-fed Darcy plug.

Uses the axisymmetric finite-volume Darcy solution of spread.py (the calculation that gives F).
From the discrete face fluxes it builds the Stokes stream function psi (0 on the axis, Q on the
side wall), the superficial speed relative to Q/A, and the advective transit time of each stream
tube, T(psi) = phi dV/dpsi, where V(psi) is the plug volume enclosed by the stream surface psi.
Outlet and interior replacement curves follow from the advective-limit relations of the paper:
    c_o(t) = (1/Q) int_{T<=t} dpsi,      c_V(t) = (1/(phi V)) int min(t, T) dpsi.
Times are reported in units of the mean residence time tbar = phi V / Q, so they do not depend
on k, mu, Q or phi.    Run:  python plug_transport.py [h_mm ...]
"""
import sys, os, json, numpy as np
from spread import solve
HERE = os.path.dirname(os.path.abspath(__file__))

R, H, A1, A2 = 16.664, 26.873, 1.93, 2.4995          # mm: plug radius, thickness, inlet and outlet port radii
NSUB, NBIN = 8, 400

def fields(h):
    F, p = solve(R, H, A1, A2, h)
    nz, nr = p.shape; dr, dz = R/nr, H/nz
    rc = (np.arange(nr)+0.5)*dr; rf = np.arange(nr+1)*dr
    # downward flux through each horizontal face (k/mu = 1); row j=nz-1 is the top (inlet) row
    down = np.zeros((nz+1, nr))
    down[1:nz] = 2*np.pi*rc*dr*(p[1:]-p[:-1])/dz
    down[nz] = np.where(rc < A1, 2*np.pi*rc*dr*(1-p[-1])/(dz/2), 0.0)
    down[0] = np.where(rc < A2, 2*np.pi*rc*dr*(p[0]-0)/(dz/2), 0.0)
    Q = down[nz].sum()
    psi = np.zeros((nz+1, nr+1)); psi[:, 1:] = np.cumsum(down, axis=1)        # at cell corners
    # superficial velocity at cell centres
    uz = 0.5*(down[:-1]+down[1:])/(2*np.pi*rc*dr)
    out = 2*np.pi*rf[1:-1]*dz*(p[:, :-1]-p[:, 1:])/dr                          # radial face fluxes
    ur = np.zeros((nz, nr)); urf = out/(2*np.pi*rf[1:-1]*dz)
    ur[:, 1:-1] = 0.5*(urf[:, :-1]+urf[:, 1:]); ur[:, 0] = 0.5*urf[:, 0]; ur[:, -1] = 0.5*urf[:, -1]
    speed = np.hypot(uz, ur)/(Q/(np.pi*R**2))
    vol = (2*np.pi*rc*dr*dz)[None, :]*np.ones((nz, 1))
    return dict(F=F, p=p, psi=psi/Q, speed=speed, vol=vol, rc=rc, rf=rf, dr=dr, dz=dz, Q=Q, nz=nz, nr=nr)

def transit(f):
    """T(psi)/tbar on NBIN flux bins, from the volume enclosed by stream surfaces (sub-cell sampling)."""
    nz, nr, dr, dz, psi = f['nz'], f['nr'], f['dr'], f['dz'], f['psi']
    s = (np.arange(NSUB)+0.5)/NSUB
    P, Wt = [], []
    for a in s:            # radial sub-position
        for b in s:        # axial sub-position
            ps = ((1-a)*(1-b)*psi[:-1, :-1]+a*(1-b)*psi[:-1, 1:]+(1-a)*b*psi[1:, :-1]+a*b*psi[1:, 1:])
            r = f['rf'][:-1]+a*dr
            P.append(ps.ravel()); Wt.append(((2*np.pi*r*dr*dz/NSUB**2)[None, :]*np.ones((nz, 1))).ravel())
    P = np.concatenate(P); Wt = np.concatenate(Wt); V = np.pi*R**2*H
    # bins uniform in flux, refined toward the wall streamline psi -> 1 where T grows
    e = np.r_[np.linspace(0, 0.9, 3*NBIN//4, endpoint=False), 1-0.1*np.geomspace(1, 1e-4, NBIN//4), 1.0]
    vol_in, _ = np.histogram(P, bins=e, weights=Wt)
    dpsi = np.diff(e); T = (vol_in/V)/dpsi                                     # in units of tbar
    return e, dpsi, monotone(T, dpsi)

def monotone(T, w):
    """Flux-weighted pooling of adjacent bins so that T does not decrease from the axis to the wall.
    It removes the sampling noise of the narrow wall-side bins and conserves the plug volume (sum of T*w)."""
    Y, W, N = [], [], []
    for a, b in zip(map(float, T), map(float, w)):
        Y.append(a); W.append(b); N.append(1)
        while len(Y) > 1 and Y[-2] > Y[-1]:
            tw = W[-2]+W[-1]; Y[-2] = (Y[-2]*W[-2]+Y[-1]*W[-1])/tw; W[-2] = tw; N[-2] += N[-1]
            Y.pop(); W.pop(); N.pop()
    return np.repeat(Y, N)

def curves(dpsi, T, t):
    o = np.argsort(T); Ts, ds = T[o], dpsi[o]
    # flux with transit time below x, taken piecewise linear between bin values (bins are narrow in flux)
    cum = np.cumsum(ds)-0.5*ds
    Tu, first = np.unique(Ts, return_index=True)                                # pooled bins share one T
    last = np.r_[first[1:], len(Ts)]-1; cu = 0.5*(cum[first]+cum[last])
    co = np.interp(t, Tu, cu, left=0.0, right=1.0)
    cv = np.array([np.minimum(x, Ts).dot(ds) for x in t])
    return co, cv

def summary(h):
    f = fields(h); e, dpsi, T = transit(f)
    t = np.geomspace(0.05, 60, 4000); co, cv = curves(dpsi, T, t)
    q = lambda c, x: float(np.interp(x, c, t))
    w = f['vol']/f['vol'].sum(); sp = f['speed']
    frac = {str(x): float(w[sp > x].sum()) for x in (0.1, 0.5, 1, 1.2, 1.5, 2, 5, 9, 12, 26, 50)}
    res = dict(h_mm=h, F=float(f['F']), mean_T=float(T.dot(dpsi)), first_arrival=float(T[:5].mean()),
               outlet=dict(t50=q(co, .5), t90=q(co, .9), t99=q(co, .99)),
               interior=dict(t50=q(cv, .5), t90=q(cv, .9), t99=q(cv, .99)),
               ratio99=q(cv, .99)/q(co, .99), vol_frac_speed_above=frac,
               speed_axis_mid=float(sp[f['nz']//2, 0]), speed_wall_mid=float(sp[f['nz']//2, -1]))
    return res, f, (e, dpsi, T, t, co, cv)

if __name__ == "__main__":
    hs = [float(x) for x in sys.argv[1:]] or [0.1, 0.05, 0.025]
    allres = []
    for h in hs:
        res, f, c = summary(h); allres.append(res); print(json.dumps(res), flush=True)
    json.dump(allres, open(os.path.join(HERE, "plug_transport.json"), "w"), indent=1)
