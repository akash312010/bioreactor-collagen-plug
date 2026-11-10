"""Own finite-volume tracer solver on a frozen Fluent flow field (for the time-step / CFL / grid study).
Space: central (default) or first-order upwind advection, two-point diffusion. Time: BDF2 (Fluent's
second-order implicit) or Crank-Nicolson. Boundary conditions as in the Fluent runs: inlet c = 1, outlet
zero diffusive flux, walls zero flux. Records outlet mass-weighted, outlet area-weighted and volume averages every step."""
import sys, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common")); from fluentbin import read_cas, read_dat
RHO = 1012.0

def build(base, D, scale=1.0, adv="central"):
    P, F, names = read_cas(base+".cas"); Dd = read_dat(base+".dat"); zt = {z: (t, n) for z, (t, n) in names.items()}
    nc = len(Dd[(1, [z for z, (t, n) in names.items() if t == 'fluid'][0])][1])
    A = np.zeros(nc); Cx = np.zeros(nc); Cy = np.zeros(nc)
    for z, (bc, f0, fa) in F.items():
        n0, n1, c0, c1 = fa[:, 0]-1, fa[:, 1]-1, fa[:, 2]-1, fa[:, 3]-1; x0, y0, x1, y1 = P[n0, 0], P[n0, 1], P[n1, 0], P[n1, 1]; cr = x0*y1-x1*y0
        for c, s in ((c0, 1.), (c1, -1.)):
            m = c >= 0; np.add.at(A, c[m], s*cr[m]/2); np.add.at(Cx, c[m], s*(x0[m]+x1[m])*cr[m]/6); np.add.at(Cy, c[m], s*(y0[m]+y1[m])*cr[m]/6)
    C = np.c_[Cx/A, Cy/A]
    rows, cols, vals = [], [], []; diag = np.zeros(nc); b_in = np.zeros(nc); outlet = None
    outflux = np.zeros(nc)                      # sum of outgoing volumetric face fluxes per cell (for Courant number)
    pe_cell = []
    for z, (bc, f0, fa) in F.items():
        t, n = zt[z]; n0, n1, c0, c1 = fa[:, 0]-1, fa[:, 1]-1, fa[:, 2]-1, fa[:, 3]-1
        Q = Dd[(18, z)][1]/RHO*scale if (18, z) in Dd else np.zeros(len(fa))   # c0 -> c1, m2/s
        L = np.hypot(*(P[n1]-P[n0]).T); fm = 0.5*(P[n0]+P[n1])
        if t == 'interior':
            dist = np.hypot(*(C[c1]-C[c0]).T); g = D*L/dist
            np.add.at(outflux, c0, np.maximum(Q, 0)); np.add.at(outflux, c1, np.maximum(-Q, 0))
            pe_cell.append(np.abs(Q)/L*dist/D)
            if adv == "central":
                wl = wr = 0.5*np.ones_like(Q)
            else:                                  # first-order upwind
                wl = (Q > 0).astype(float); wr = 1-wl
            # flux out of c0 = Q (wl c0 + wr c1) ; into c1 the same
            np.add.at(diag, c0, Q*wl+g); rows += list(c0); cols += list(c1); vals += list(Q*wr-g)
            np.add.at(diag, c1, -Q*wr+g); rows += list(c1); cols += list(c0); vals += list(-Q*wl-g)
        elif n == 'inlet':
            dist = np.hypot(*(fm-C[c0]).T); g = D*L/dist
            np.add.at(diag, c0, g); np.add.at(b_in, c0, -Q+g)
        elif n == 'outlet':
            np.add.at(diag, c0, Q); np.add.at(outflux, c0, np.maximum(Q, 0)); outlet = (c0, Q, L)
    K = sp.csr_matrix((vals, (rows, cols)), shape=(nc, nc))+sp.diags(diag)
    return dict(K=K, b=b_in, A=A, outlet=outlet, outflux=outflux, pe_cell=np.concatenate(pe_cell), C=C)

def courant(S, n_per_tid):
    A = S["A"]; q = S["outlet"][1].sum(); tid = A.sum()/q; dt = tid/n_per_tid
    co = dt*S["outflux"]/A
    w = A/A.sum()
    return dict(t_ideal=tid, dt=dt, Co_max=float(co.max()), Co_p99=float(np.percentile(co, 99)), Co_p95=float(np.percentile(co, 95)),
                Co_mean_area=float(np.sum(w*co)), Co_median=float(np.median(co)),
                cellPe_max=float(S["pe_cell"].max()), cellPe_p99=float(np.percentile(S["pe_cell"], 99)))

def run(S, n_per_tid, nT, scheme="bdf2"):
    K, b, A = S["K"], S["b"], S["A"]; oc, oQ, oL = S["outlet"]
    q = oQ.sum(); Atot = A.sum(); tid = Atot/q; dt = tid/n_per_tid
    M = sp.diags(A/dt)
    if scheme == "cn":
        lu = spl.splu((M+0.5*K).tocsc()); R = (M-0.5*K).tocsr()
    else:
        lu1 = spl.splu((M+K).tocsc()); lu = spl.splu((1.5*M+K).tocsc())
    c_old = np.zeros(len(A)); c = np.zeros(len(A)); nsteps = int(round(nT*n_per_tid))
    out = np.zeros((nsteps+1, 4)); cmin, cmax = 0.0, 0.0
    for k in range(1, nsteps+1):
        if scheme == "cn":
            cn = lu.solve(R@c+b)
        elif k == 1:
            cn = lu1.solve(M@c+b)                  # first step backward Euler
        else:
            cn = lu.solve(M@(2*c-0.5*c_old)+b)
        c_old, c = c, cn
        out[k] = (k/n_per_tid, np.sum(oQ*c[oc])/q, np.sum(oL*c[oc])/oL.sum(), np.sum(A*c)/Atot)
        cmin = min(cmin, c.min()); cmax = max(cmax, c.max())
    return out, tid, cmin, cmax

def times(out, fr=(0.5, 0.9, 0.99)):
    r = {}
    for j, nm in ((1, "mass"), (2, "area"), (3, "vol")):
        y = out[:, j]
        for f in fr:
            i = np.argmax(y >= f)
            r[f"{nm}_t{int(f*100)}"] = float(np.interp(f, [y[i-1], y[i]], [out[i-1, 0], out[i, 0]])) if y.max() >= f and i > 0 else float("nan")
    return r
