"""Exact-rational stress tests for the Stage11R-P1c uniform-tail reduction.

Checks necessary identities over representative cutoffs, not publication priority,
full-domain formal verification or a general-equilibrium uniqueness proof.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json


def phi(m):
    return m - m * m / 2


cuts = (F(0), F(1, 10), F(1, 4), F(1, 2), F(3, 4), F(9, 10), F(1))
checks = 0
for ell in cuts:
    m = 1 - ell
    q = (1 - ell * ell) / 2
    assert q == phi(m)
    assert q - m / 2 == m * (1 - m) / 2
    assert q >= m / 2
    checks += 3

# Absence: m=q=0. Public-only and public/private: m=1, q=1/2.
assert phi(F(0)) == 0
assert phi(F(1)) == F(1, 2)
checks += 2

for ell_a, ell_b in product(cuts, repeat=2):
    ma, mb = 1 - ell_a, 1 - ell_b
    qa, qb = (1 - ell_a**2) / 2, (1 - ell_b**2) / 2
    assert qa * qb == phi(ma) * phi(mb)
    d = qa * qb - ma * mb / 4
    assert d >= 0
    # Identity at full participation, or whenever one side abstains entirely.
    if ma == 0 or mb == 0 or (ma == 1 and mb == 1):
        assert d == 0
    checks += 2

result = {
    "stage": "11R-P1c",
    "frozen_specification": "uniform_0_1_upper_tail",
    "verdict": "PASS: single-index reduction and fixed-mass comparator identities",
    "rational_cases": len(cuts),
    "paired_cases": len(cuts)**2,
    "assertions": checks,
    "scope_limit": "Finite rational regression; symbolic identity follows algebraically from ell=1-m, not from enumeration.",
    "novelty": "No theorem-level priority or non-absorption conclusion",
}
Path(__file__).with_name("stage11r_p1c_reduction_checks.json").write_text(
    json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
