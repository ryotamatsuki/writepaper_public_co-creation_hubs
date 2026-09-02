"""Stage 4R-MC-G: strict frozen social-planner accounting.

Provenance:
- Python 3.13.5
- SymPy 1.14.0
- Metro m is project-side only, as frozen in Stage 3.
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

SW_fixed = 4*(A1+A2)+2*m*L-gamma*(x1**2+x2**2)/2
SW = sp.factor(SW_fixed.subs(q, qstar))
Fsp1 = sp.factor(sp.diff(SW, x1))

x = sp.symbols("x", positive=True)
B = 12*c**2+28*c*u+15*u**2
Fsp_sym = sp.factor(Fsp1.subs({x1:x,x2:x}))
Fsp_target = (B*x**2-8*gamma*u*x**3-4*m**2)/(8*u*x**2)
assert sp.simplify(Fsp_sym-Fsp_target) == 0

# Nash FOC substitution: solve gamma from G symmetric FOC.
gamma_N = (3*B*x**2+4*a*m*x+4*m**2)/(32*u*x**3)
Fsp_at_N = sp.factor(Fsp_sym.subs(gamma, gamma_N))
FspN_target = (B*x**2-4*a*m*x-20*m**2)/(32*u*x**2)
assert sp.simplify(Fsp_at_N-FspN_target) == 0

y = sp.symbols("y", positive=True)
psi = sp.factor((32*u*x**2*Fsp_at_N).subs(m,y*x)/x**2)
y_star = sp.solve(sp.Eq(psi,0), y)

print("F_SP sym =", Fsp_sym)
print("F_SP at Nash =", Fsp_at_N)
print("psi(y) =", psi)
print("positive roots =", y_star)
