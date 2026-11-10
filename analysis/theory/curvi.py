"""Second-order lubrication theory in centreline (curvilinear) coordinates.
(s, n): arc length along a reference curve and normal distance (positive to the left); scale factor H = 1 - kappa n.
Walls at n = c(s) -+ h(s). Slow variable X = eps s, kappa = eps k(X). Stream function psi: u_s = psi_n, u_n = -psi_s/H.
Stokes: L(L(psi)) = 0, L f = (1/H)[d_s(H^-1 d_s f) + d_n(H d_n f)].  Momentum: (1/H) p_s = -mu omega_n, omega = -L psi.
Output: dp/ds at n = 0 relative to the leading lubrication value, to O(eps^2)."""
import sympy as sp
X, n, e, Q, mu = sp.symbols('X n epsilon Q mu', real=True)
h = sp.Function('h')(X); c = sp.Function('c')(X); k = sp.Function('k')(X)
H = 1 - e*k*n
def L(f):
    return (e*sp.diff(e*sp.diff(f, X)/H, X) + sp.diff(H*sp.diff(f, n), n))/H
eta = (n-c)/h
psi0 = Q/2*(3*eta - eta**3)/2
def solve_order(rhs):
    # d^4 psi / dn^4 = rhs, psi = psi_n = 0 at n = c -+ h
    Y = sp.symbols('Y'); r = sp.expand(sp.simplify(rhs.subs(n, Y+c)))
    P = r
    for _ in range(4): P = sp.integrate(P, Y)
    a0, a1, a2, a3 = sp.symbols('a0:4')
    P = P + a0 + a1*Y + a2*Y**2 + a3*Y**3
    eqs = [P.subs(Y, h), P.subs(Y, -h), sp.diff(P, Y).subs(Y, h), sp.diff(P, Y).subs(Y, -h)]
    s = sp.solve(eqs, [a0, a1, a2, a3], dict=True)[0]
    return sp.expand(sp.simplify(P.subs(s).subs(Y, n-c)))
psi1s, psi2s = sp.symbols('psi1s psi2s')
f1 = sp.Function('f1')(X, n); f2 = sp.Function('f2')(X, n)
full = L(L(psi0 + e*f1 + e**2*f2))
ser = sp.series(full, e, 0, 3).removeO()
o1 = sp.expand(ser.coeff(e, 1)); o2 = sp.expand(ser.coeff(e, 2))
# order 1: d4 f1 + (rest) = 0
rest1 = sp.simplify(o1 - sp.diff(f1, n, 4))
psi1 = solve_order(-rest1)
rest2 = o2 - sp.diff(f2, n, 4)
rest2 = sp.simplify(rest2.subs(f1, psi1).doit())
psi2 = solve_order(-rest2)
psi = psi0 + e*psi1 + e**2*psi2
omega = -L(psi)
ps = -mu*H*sp.diff(omega, n)            # dp/ds (since (1/H) p_s = -mu omega_n)
ps0 = sp.simplify(ps.subs(n, 0))
ser = sp.series(ps0, e, 0, 3).removeO()
lead = sp.simplify(ser.coeff(e, 0)); c1 = sp.simplify(ser.coeff(e, 1)/lead); c2 = sp.simplify(ser.coeff(e, 2)/lead)
print("leading dp/ds at n=0:", sp.factor(lead))
print("O(eps) relative:", sp.factor(c1))
print("O(eps^2) relative:", sp.expand(c2))
import pickle; pickle.dump(dict(lead=lead, c1=c1, c2=c2), open("curvi.pkl", "wb"))
