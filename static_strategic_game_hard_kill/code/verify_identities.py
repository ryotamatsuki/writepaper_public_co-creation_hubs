"""Stage 4R-MC-G: exact identity and characteristic-space verification.

Provenance:
- Python 3.13.5
- SymPy 1.14.0
"""
import sympy as sp

c, u, x1, x2, theta, kappa, m, p = sp.symbols(
    "c u x1 x2 theta kappa m p", positive=True
)
d = c+u
Q1, Q2 = x1*d, x2*d
t1, t2 = x1*u, x2*u

U1_csds = x1*(c+u*(1-theta))-kappa
U2_csds = x2*(c+u*theta)-kappa
U1_char = Q1-t1*theta-kappa
U2_char = Q2-t2*(1-theta)-kappa
UT = m-p-kappa

assert sp.simplify(U1_csds-U1_char) == 0
assert sp.simplify(U2_csds-U2_char) == 0
assert sp.simplify(t1-(u/d)*Q1) == 0
assert sp.simplify(t2-(u/d)*Q2) == 0

# Same provider cutoffs follow mechanically.
th12_csds = sp.solve(sp.Eq(U1_csds,U2_csds),theta)[0]
th12_char = sp.solve(sp.Eq(U1_char,U2_char),theta)[0]
assert sp.simplify(th12_csds-th12_char) == 0

print("Exact characteristic-space mapping verified.")
print("U1 =", U1_char)
print("U2 =", U2_char)
print("UT =", UT)
print("theta12 =", sp.factor(th12_char))
