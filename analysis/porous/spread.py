"""Darcy flow through a cylindrical plug (radius R, height H) fed through a central disk of radius a1 on top and
drained through a central disk of radius a2 at the bottom; all other plug faces impermeable.
Axisymmetric finite volumes, equipotential disks. Returns F = (actual resistance)/(1D resistance mu H/(k pi R^2))."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, sys
def solve(R, H, a1, a2, h):
    nr, nz = int(round(R/h)), int(round(H/h)); dr, dz = R/nr, H/nz
    rc = (np.arange(nr)+0.5)*dr; idx = lambda i, j: j*nr+i
    N = nr*nz; rows, cols, vals = [], [], []; diag = np.zeros(N); b = np.zeros(N)
    rf = np.arange(nr+1)*dr
    for j in range(nz):
        for i in range(nr):
            k = idx(i, j)
            # radial faces
            if i > 0: g = rf[i]*dz/dr; diag[k] += g; rows.append(k); cols.append(idx(i-1, j)); vals.append(-g)
            if i < nr-1: g = rf[i+1]*dz/dr; diag[k] += g; rows.append(k); cols.append(idx(i+1, j)); vals.append(-g)
            # axial faces
            g = rc[i]*dr/dz
            if j > 0: diag[k] += g; rows.append(k); cols.append(idx(i, j-1)); vals.append(-g)
            elif rc[i] < a2: diag[k] += 2*g          # bottom outlet disk, p = 0
            if j < nz-1: diag[k] += g; rows.append(k); cols.append(idx(i, j+1)); vals.append(-g)
            elif rc[i] < a1: diag[k] += 2*g; b[k] += 2*g*1.0   # top inlet disk, p = 1
    A = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))+sp.diags(diag)
    p = spl.spsolve(A.tocsc(), b).reshape(nz, nr)
    top = rc < a1; Q = np.sum(2*rc[top]*dr/dz*(1.0-p[-1, top]))*2*np.pi/2   # flux through inlet disk, times 2*pi (k/mu = 1)
    Q = np.sum(2*(1.0-p[-1, top])/dz*2*np.pi*rc[top]*dr)
    Rnum = 1.0/Q; R1D = H/(np.pi*R**2)
    return Rnum/R1D, p
if __name__ == "__main__":
    R, H, a1, a2 = 16.664, 26.873, 1.93, 2.4995
    for h in (0.2, 0.1, 0.05, 0.025):
        F, p = solve(R, H, a1, a2, h); print(f"h = {h} mm: F = {F:.4f}", flush=True)
    # Holm half-space estimate: R = mu/(4 k a) per equipotential disk
    print("half-space estimate F =", 1+(1/(4*a1)+1/(4*a2))/(H/(np.pi*R**2)))
