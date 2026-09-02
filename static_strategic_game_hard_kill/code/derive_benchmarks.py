"""Stage 4R-MC-G: B1/B2/B3 exact benchmark diagnostics.

Provenance:
- Python 3.13.5
- SymPy 1.14.0
"""
import sympy as sp

c, u, x1, x2, m, q, gamma, th = sp.symbols(
    "c u x1 x2 m q gamma theta", positive=True
)
d = c+u
a = 2*c+u

# B1: one public hub + metro.
th1 = (x1*d-q)/(x1*u)
A1 = sp.integrate(x1*(d-u*th), (th, 0, th1))
L1 = 1-th1
P1 = A1 + q*L1  # common -kappa omitted from derivatives
W1_b1_fixed = P1 + 2*A1 - gamma*x1**2/2
q_b1 = (m+x1*c)/2
W1_b1 = sp.factor(W1_b1_fixed.subs(q, q_b1))
F_b1 = sp.factor(sp.diff(W1_b1, x1))
SOC_b1 = sp.factor(sp.diff(W1_b1, x1, 2))

# B2: two public hubs, no metro, full participation.
th12 = (x1*d-x2*c)/(u*(x1+x2))
A1_b2 = sp.integrate(x1*(d-u*th), (th, 0, th12))
A2_b2 = sp.integrate(x2*(c+u*th), (th, th12, 1))
W1_b2 = sp.factor(3*A1_b2 + A2_b2 - gamma*x1**2/2)
C_b2 = sp.factor(sp.diff(W1_b2, x1, x2))
C_b2_target = x1*x2*a**2*(x1-5*x2)/(u*(x1+x2)**4)
assert sp.simplify(C_b2-C_b2_target) == 0

# symmetric B2 Nash and planner.
x = sp.symbols("x", positive=True)
F_b2_sym = sp.factor(sp.diff(W1_b2, x1).subs({x1:x, x2:x}))
xN_b2 = sp.solve(sp.Eq(F_b2_sym, 0), x)[0]
SW_b2 = 4*(A1_b2+A2_b2)-gamma*(x1**2+x2**2)/2
Fsp_b2_sym = sp.factor(sp.diff(SW_b2, x1).subs({x1:x, x2:x}))
xSP_b2 = sp.solve(sp.Eq(Fsp_b2_sym, 0), x)[0]

# B3 fixed metro offer, public-metro-public interior.
th1T = (x1*d-q)/(x1*u)
th2T = (q-x2*c)/(x2*u)
A1_b3 = sp.integrate(x1*(d-u*th), (th, 0, th1T))
A2_b3 = sp.integrate(x2*(c+u*th), (th, th2T, 1))
L = th2T-th1T
P = A1_b3+A2_b3+q*L
W1_b3 = sp.factor(P+2*A1_b3-gamma*x1**2/2)
F_b3 = sp.factor(sp.diff(W1_b3, x1))
C_b3 = sp.factor(sp.diff(W1_b3, x1, x2))
assert sp.simplify(C_b3) == 0

print("B1 FOC =", F_b1)
print("B1 SOC =", SOC_b1)
print("B2 cross =", C_b2)
print("B2 x_N =", sp.factor(xN_b2))
print("B2 x_SP =", sp.factor(xSP_b2))
print("B2 wedge =", sp.factor(xN_b2-xSP_b2))
print("B3 FOC =", F_b3)
print("B3 cross =", C_b3)
