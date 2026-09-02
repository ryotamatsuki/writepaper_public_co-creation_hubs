"""Stage 3M-4M symbolic verification for the productive-matching model."""

import sympy as sp

A, D, p, r, u, v, k = sp.symbols("A D p r u v k", positive=True, real=True)
mi, mj = sp.symbols("mi mj", integer=True)

# Match-probability primitives.
h = 2*p - p**2
k_ref = (1-p)**2 * r
t = h + k_ref

# Realized homogeneous-Cournot equilibrium.
qi = (A + 2*D*mi - D*mj) / 3
qj = (A + 2*D*mj - D*mi) / 3
Q = sp.expand(qi + qj)

# Bernoulli expectations, using E[m_i^2]=E[m_i]=u and independence E[m_i m_j]=uv.
E_qi2 = (
    A**2 + 4*A*D*u - 2*A*D*v + 4*D**2*u - 4*D**2*u*v + D**2*v
) / 9
E_Q2 = (
    4*A**2 + 4*A*D*u + 4*A*D*v + D**2*u + 2*D**2*u*v + D**2*v
) / 9
R = sp.factor(E_qi2 + E_Q2/4)

R_expected = (
    8*A**2 + 20*A*D*u - 4*A*D*v - 14*D**2*u*v + 17*D**2*u + 5*D**2*v
) / 36
assert sp.simplify(R - R_expected) == 0

D0 = sp.factor(R.subs({u:h, v:p}) - R.subs({u:0, v:0}))
D1 = sp.factor(R.subs({u:t, v:t}) - R.subs({u:p, v:h}))
Gamma = sp.factor(D1 - D0)

Gamma0 = sp.factor(Gamma.subs(r, 0))
Gamma0_target = -D*p*(8*A + D*(11 - 14*p**2 + 7*p**3))/18
assert sp.simplify(Gamma0-Gamma0_target) == 0

# Compact nested representation using k as the referral-induced increment.
Gamma_k = sp.factor(
    Gamma0_target + D*k*(8*A + D*(11 - 14*h - 7*k))/18
)
assert sp.simplify(Gamma - Gamma_k.subs(k, k_ref)) == 0

# Referral monotonicity derivative.
dGamma_dk = sp.factor(sp.diff(Gamma_k, k))

# Downstream-competition threshold wedge: nonstrategic matching threshold is k=p.
Gamma_at_p = sp.factor(Gamma_k.subs(k, p))
Gamma_at_p_target = -7*D**2*p**2*(p**2 - 4*p + 5)/18
assert sp.simplify(Gamma_at_p-Gamma_at_p_target) == 0

# Maximal referral.
Gamma_max = sp.factor(Gamma_k.subs(k, (1-p)**2))
Gamma_max_target = D*(
    8*A*(1-3*p+p**2) + D*(4-33*p+39*p**2-14*p**3)
)/18
assert sp.simplify(Gamma_max-Gamma_max_target) == 0

# Cross-region first-hub externality.
X0 = sp.factor(R.subs({u:p, v:h}) - R.subs({u:0, v:0}))
X0_target = D*p*(12*A+4*A*p+D*(27-33*p+14*p**2))/36
assert sp.simplify(X0-X0_target) == 0

# Exact sign witnesses.
cal = {A:sp.Rational(1), D:sp.Rational(1,2), p:sp.Rational(1,10)}
neg = sp.factor(Gamma.subs({**cal, r:sp.Rational(1,10)}))
pos = sp.factor(Gamma.subs({**cal, r:sp.Rational(1,5)}))
assert neg == -sp.Rational(761087, 72000000)
assert pos == sp.Rational(33521, 2250000)

print("R_i(u,v) =", R)
print("D0 =", D0)
print("D1 =", D1)
print("Gamma =", Gamma)
print("Gamma0 =", Gamma0)
print("dGamma/dk =", dGamma_dk)
print("Gamma(k=p) =", Gamma_at_p)
print("Gamma_max =", Gamma_max)
print("X0 =", X0)
print("negative witness =", neg)
print("positive witness =", pos)
print("SYMBOLIC CHECKS: PASS")
