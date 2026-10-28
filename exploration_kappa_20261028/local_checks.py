#!/usr/bin/env python3
"""Exact algebra and a copula counterexample. Standard library only.

The analytic proof is in theory_independent.md / REPORT.md. Rational grids
are regression checks, not substitutes for universally quantified proofs.
No calibrated primitive is changed by this file.
"""
from fractions import Fraction as F
import json
from pathlib import Path


# Sparse polynomials in (c, kappa), checked by coefficient equality.
def add(a, b):
    out = dict(a)
    for mon, val in b.items():
        out[mon] = out.get(mon, F(0)) + val
    return {mon: val for mon, val in out.items() if val}


def mul(a, b):
    out = {}
    for (i, j), v in a.items():
        for (k, l), w in b.items():
            mon = (i+k, j+l)
            out[mon] = out.get(mon, F(0)) + v*w
    return {mon: val for mon, val in out.items() if val}


def scale(a, x):
    return {mon: val*x for mon, val in a.items() if val*x}


mu = {(0, 0): F(1, 2), (0, 1): -F(1, 2), (1, 1): F(1)}
nu = {(0, 0): F(1, 3), (0, 1): -F(1, 3), (2, 1): F(1)}
fh = {(0, 0): F(1, 4), (0, 1): F(1, 4)}
mh = {(0, 0): F(1, 4), (1, 1): F(1), (0, 1): -F(3, 4)}
Q = {(0, 0): F(5, 24), (2, 1): F(1), (1, 1): -F(1, 4),
     (0, 1): -F(1, 3), (0, 2): F(1, 8), (1, 2): -F(1, 4)}
assert add(nu, scale(mul(mu, fh), -1)) == Q
variance = add(nu, scale(mul(mu, mu), -1))
assert add(variance, mul(mu, mh)) == Q  # d*q_h=Var+mu*d*m_h

checks = 2
for ci in range(1, 50):
    c = F(ci, 100)  # symbolic proof does not depend on this finite grid
    km = 1/(3-4*c)
    def qh(k):
        return F(5,24)+k*(c*c-c/4-F(1,3))+(F(1,2)-c)*k*k/4
    assert qh(0) > 0 and qh(1) == c*(c-F(1,2)) < 0
    mu_at = (1-km)/2+km*c
    var_at = (1-km)/3+km*c*c-mu_at*mu_at
    assert qh(km) == var_at > 0
    assert F(1,3) < km < 1
    checks += 4

# An alternative one-parameter copula with the SAME uniform margins,
# Pearson correlation kappa, independent endpoint 0 and diagonal endpoint 1:
# kappa*M+(1-kappa)*[1+beta(kappa)*(2theta-1)*g(z)] dtheta dz,
# beta=4*kappa*(1-kappa), g(z)=6z^2-6z+1.
# Integrals of g and z*g vanish; |g|<=1, beta in [0,1].
assert F(6,3)-F(6,2)+1 == 0
assert F(6,4)-F(6,3)+F(1,2) == 0
c, k = F(1,4), F(13,25)
beta = 4*k*(1-k)
g, gp = 6*c*c-6*c+1, 12*c-6
mix = F(1,4)+k*(c-F(3,4))
correction = (1-k)*beta*(g-(1-c)*gp)/12
alternative = mix+correction
assert mix == -F(1,100) < 0 < alternative
checks += 3

result = {
    'scope': 'Coefficient identities and exact copula sign counterexample; no full-game proof',
    'assertions': checks,
    'copula_counterexample': {
        'c': str(c), 'kappa': str(k), 'beta': str(beta),
        'd_times_mass_h_mixture': str(mix),
        'd_times_mass_h_alternative': str(alternative),
        'same_uniform_margins_and_corr': True,
        'same_independent_and_aligned_endpoints': True,
        'extra_free_parameters': 0,
    },
    'formal_lean_status': 'New results are NOT Lean formalized; original source untouched',
}
Path(__file__).with_name('local_checks.json').write_text(
    json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
print(json.dumps(result, ensure_ascii=False, indent=2))
