"""Stage 4R-MC-G: metro demand and interior price solution.

Provenance:
- Python 3.13.5
- SymPy 1.14.0
- Frozen CSDS primitives only.
"""
import sympy as sp

c, u, x1, x2, m, q = sp.symbols("c u x1 x2 m q", positive=True)
a = 2*c + u
d = c + u
K = a*x1*x2/(x1+x2)

th1T = (x1*d-q)/(x1*u)
th2T = (q-x2*c)/(x2*u)
ell = sp.factor(th2T-th1T)
D = sp.factor(2*ell)
profit = sp.factor((m-q)*D)
FOC = sp.factor(sp.diff(profit, q))
qstar = sp.solve(sp.Eq(FOC, 0), q)[0]
pstar = sp.factor(m-qstar)
SOC = sp.factor(sp.diff(profit, q, 2))

dpdx1 = sp.factor(sp.diff(pstar, x1))
dpdx2 = sp.factor(sp.diff(pstar, x2))

assert sp.simplify(ell-(x1+x2)*(q-K)/(u*x1*x2)) == 0
assert sp.simplify(qstar-(m+K)/2) == 0
assert sp.simplify(pstar-(m-K)/2) == 0
assert sp.simplify(dpdx1 + a*x2**2/(2*(x1+x2)**2)) == 0
assert sp.simplify(dpdx2 + a*x1**2/(2*(x1+x2)**2)) == 0

print("D_T =", D)
print("q* =", sp.factor(qstar))
print("p* =", pstar)
print("SOC =", SOC)
print("dp*/dx1 =", dpdx1)
print("dp*/dx2 =", dpdx2)
