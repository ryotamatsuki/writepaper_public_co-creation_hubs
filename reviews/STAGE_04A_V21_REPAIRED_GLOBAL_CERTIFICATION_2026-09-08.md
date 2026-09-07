# Stage 4A v2.1 — Independent Mathematical Adversarial Certification

Date: 2026-09-08 JST

Canonical workflow: `ryotamatsuki/research-paper-workflow` v2.1 @ `fe6fa5f2a94632be578d5e9fab39e7df86013113`

Canonical template: `templates/STAGE_04A_MATH_RED_TEAM.md`

Upstream repaired Stage-4 construction: `stage4_v21_global/`

## Executive verdict

**GO — MATHEMATICAL ADVERSARIAL CERTIFICATION PASS**

The repaired construction survives an implementation-independent all-regime attack. The clean-room evaluator does not import the Stage-4 repair solver, the old central-regime solver, or the old exact verifier. It reconstructs route choice, participation, private pricing, public welfare, finite deviations, and local derivatives from primitives.

The certified witness uses the existing model only, with

- `beta = 0.01`,
- `gamma = 0.825`,
- `tau = 0.35`,

and all other primitives unchanged.

The earlier canonical vector (`beta=.05, gamma=.9, tau=.05`) remains rejected because of profitable finite deviations. It is permanently non-authoritative for global-equilibrium claims.

## Independent reconstruction summary

The independent evaluator reconstructs the four project routes (`0`, `H1`, `H2`, `HT`) from the upper envelope of project utility, solves partner participation by fixed point from multiple starts, reoptimizes the private price after public deviations, and computes regional welfare directly from project surplus, partner surplus, and public investment cost.

Because `kL + tau = 0.62 > v + alpha = 0.60`, a nonresident rival public hub is analytically dominated by the outside option for every project type and every history. This removes the exact free-riding route that destroyed the old witness, without changing the model architecture.

## Certified on-path objects

Production Stage-4 construction:

- `x_G = 0.8371022382025995`
- `x_B3 = 0.8258903860237495`
- `p_G ≈ 0.0184036832`

Clean-room reconstruction:

- independent private price: `0.0184036796`
- independent global BR in G: `0.8371019830`
- G best-deviation gain over candidate: `7.47e-12`
- independent global BR in B3: `0.8258910554`
- B3 best-deviation gain over candidate: `-4.20e-14`

Full multistart recheck at the detected G best reply gives a welfare difference of `-6.84e-10`; B3 gives `-4.20e-14`. No profitable finite deviation is detected.

## Boundary / regime / continuation audit

The independent audit covers the full public strategy interval `[0,1]` by global sweep plus local refinement. It deliberately rechecks the old low-investment failure region, endpoints, high-investment corners, and both repaired candidate neighborhoods. Material private continuations are reoptimized rather than filtered when the preferred central branch fails.

Targeted full-multistart histories include `x = 0, .10, .18, .20, .35, .50, x_B3, x_G, .95, 1.0`. These histories are solved from seven participation starts in the intensive recheck. The result is `TARGETED_FULL_MULTISTART: PASS`.

No `UNRESOLVED` or material multiplicity state is treated as an unprofitable deviation.

## Independent strategic-sign audit

Finite-difference step-size attacks reproduce the same qualitative strategic relationship:

- G: `H11 < 0`, `H12 < 0`, therefore `BR_G' < 0`.
- B3: `H11 < 0`, `H12 > 0`, therefore `BR_B3' > 0`.

Representative independent values:

- G slope: about `-0.01986`.
- B3 slope: about `+0.003968`.

The sign survives `h = 6e-4, 4e-4, 3e-4`.

Thus the repaired witness satisfies

`BR_B3' > 0 > BR_G'`.

## Counterexample preservation

The previous Stage-4A counterexample remains a permanent regression artifact. It establishes that FOC/SOC and exact local stationary-root isolation are not sufficient for a global public Nash/SPNE claim. Future solvers must continue to fail closed on omitted route regimes and unresolved continuations.

## Scope discipline

PASS does **not** establish a primitive necessary-and-sufficient characterization, a global sign theorem over the full parameter space, or uniqueness for all parameter values. It establishes that the same model contains a repaired, all-regime computationally certified equilibrium witness with the headline public-public strategic sign reversal, together with the already established local analytic mechanism.

The manuscript must not continue to use the rejected old canonical vector as a global equilibrium witness.

## CI evidence

GitHub Actions run `34150133823` on head `03a045fe1e4ac6951946fdc95da4a5c887d3ea77`: **SUCCESS**.

Observed terminal markers include:

- `RIVAL_PUBLIC_NONRESIDENT_DOMINANCE: PASS`
- `TARGETED_FULL_MULTISTART: PASS`
- `STAGE4A_REPAIRED_GLOBAL_CERTIFICATION: PASS`

## Canonical routing

**GO TO STAGE 6 — NOVELTY RE-KILL.**

Stage 6 must evaluate the repaired surviving claim set, not the stale pre-Stage-4A manuscript wording or the rejected old numerical witness.
