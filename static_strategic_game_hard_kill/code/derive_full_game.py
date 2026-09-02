"""Stage 4R-MC-G: full strategic game G.

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
K = a*x1*x2/(x1+x2)
qstar = (m+K)/2

th1T = (x1*d-q)/(x1*u)
th2T = (q-x2*c)/(x2*u)
A1 = sp.integrate(x1*(d-u*th), (th, 0, th1T))
A2 = sp.integrate(x2*(c+u*th), (th, th2T, 1))
L = th2T-th1T
P = A1+A2+q*L
W1_fixed = sp.factor(P+2*A1-gamma*x1**2/2)

W1 = sp.factor(W1_fixed.subs(q, qstar))
F1 = sp.factor(sp.diff(W1, x1))
C12 = sp.factor(sp.diff(W1, x1, x2))
C_target = -a*(a*x1*x2*(x1+7*x2)+2*m*(x2**2-x1**2))/(4*u*(x1+x2)**4)
assert sp.simplify(C12-C_target) == 0

Wq = sp.factor(sp.diff(W1_fixed, q))
qx1 = sp.factor(sp.diff(qstar, x1))
decomp = sp.factor(sp.diff(W1_fixed, x1).subs(q, qstar)+Wq.subs(q,qstar)*qx1-F1)
assert sp.simplify(decomp) == 0

x = sp.symbols("x", positive=True)
B = 12*c**2+28*c*u+15*u**2
Fsym = sp.factor(F1.subs({x1:x,x2:x}))
Fsym_target = (3*B*x**2+4*a*m*x+4*m**2-32*gamma*u*x**3)/(32*u*x**2)
assert sp.simplify(Fsym-Fsym_target) == 0
Csym = sp.factor(C12.subs({x1:x,x2:x}))
assert sp.simplify(Csym + a**2/(8*u*x)) == 0

print("q* =", sp.factor(qstar))
print("F1^G =", F1)
print("W12^G =", C12)
print("F_sym =", Fsym)
print("W12_sym =", Csym)
