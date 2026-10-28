# Kappa-only exploration, stopped with a negative adoption decision

Base: PR #59 head `d2520c9afbf9b2b742a5f04fbba6159eb489db61`.
The frozen manuscript, Lean source, primitives, and existing checks are unchanged.

Read [REPORT.md](REPORT.md) first (Japanese final synthesis).
The decision is **NO-GO for adopting this extension as an upper-journal project; STOP exploration; preserve the JICT baseline and its existing HOLD**.
There are valid model-specific local propositions; theorem-level novelty and full extended-game certification are not established.

| File | Content |
|---|---|
| [theory_independent.md](theory_independent.md) | Moment formulas, global conditional fee, local three regimes, policy cross derivative, proofs |
| [literature.md](literature.md) | Original-text access ledger, proposition-level comparison, unresolved priority |
| [REFEREE.md](REFEREE.md) | Independent AI referee assessment and stopping decision |
| [NUMERICAL_REPORT.md](NUMERICAL_REPORT.md) | Exact vs numerical scopes, frozen regression and all-profile sweep |
| [moments.py](moments.py) | Exact half-plane moment integrals, corners, diagonal and independent cases |
| [local_checks.py](local_checks.py) | Exact polynomial identities and equal-correlation copula counterexample |
| [exact_threshold.py](exact_threshold.py) | Exact conditional BP/BB threshold with unchanged original costs |
| [numerics.py](numerics.py) | Joint P/P participation and posted-price best-response diagnostics |
| [results.json](results.json) | All 16 payoff configurations at 11 kappa values, three original scenarios |
| [frozen_regression.log](frozen_regression.log) | Original five-check regression output |
| [freeze_integrity.json](freeze_integrity.json) | Hashes and base commit for unchanged frozen assets |

The original pooled technology case is reproduced solely as the existing comparator; it is not an additional extension.
No external human referee, new empirical data, or new Lean certification is claimed.

From the repository root, standard-library exact checks:

```bash
python3 exploration_kappa_20261028/local_checks.py
python3 exploration_kappa_20261028/exact_threshold.py
python3 exploration_kappa_20261028/integrity_check.py
```

Numerical reproduction (Python 3.11+):

```bash
python -m pip install -r exploration_kappa_20261028/requirements.txt
python exploration_kappa_20261028/numerics.py
python exploration_kappa_20261028/numerics.py --seed-checks
```

The numerical checks do not prove global price equilibrium for kappa below one.
The participation fixed point is separately unique for all fees and kappa under an exact contraction bound.
The exact conditional public-entry threshold is certified for the original skill-half setting; it is not a uniform full-SPNE certificate.
See the dedicated GitHub Actions workflow for the unchanged five baseline checks plus these extension checks.

Dependencies: NumPy 2.3.5, SciPy 1.17.0 for numerical work; exact checks use only the Python standard library.
