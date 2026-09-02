"""Stage 4R-MC-G: project-choice partitions.

Provenance:
- Python 3.13.5
- SymPy 1.14.0
- Frozen CSDS primitives only.
"""
import sympy as sp

th, c, u, x1, x2, q, kappa = sp.symbols(
    "theta c u x1 x2 q kappa", positive=True
)
d = c + u
a = 2*c + u

V1 = x1*(d-u*th)
V2 = x2*(c+u*th)

th12 = sp.solve(sp.Eq(V1, V2), th)[0]
th1T = sp.solve(sp.Eq(V1, q), th)[0]
th2T = sp.solve(sp.Eq(V2, q), th)[0]
th10 = sp.solve(sp.Eq(V1, kappa), th)[0]
th20 = sp.solve(sp.Eq(V2, kappa), th)[0]
K = sp.simplify(V1.subs(th, th12))

assert sp.simplify(th12-(x1*d-x2*c)/(u*(x1+x2))) == 0
assert sp.simplify(th1T-(x1*d-q)/(x1*u)) == 0
assert sp.simplify(th2T-(q-x2*c)/(x2*u)) == 0
assert sp.simplify(K-a*x1*x2/(x1+x2)) == 0

print("theta_12 =", sp.factor(th12))
print("theta_1T =", sp.factor(th1T))
print("theta_2T =", sp.factor(th2T))
print("theta_10 =", sp.factor(th10))
print("theta_20 =", sp.factor(th20))
print("K =", sp.factor(K))
