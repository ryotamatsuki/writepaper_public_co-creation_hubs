# Theorem Certificate — Local Strategic Sign Reversal

Claim ID: `thm:reversal`

Stage-4A date: 2026-09-07

## Exact claim

Conditional on regular interior B3 and G branches existing at `beta=0`, nonsingular participation/equilibrium Jacobians, strict private/public SOCs, and smooth continuation, there exists `epsilon>0` such that for every sufficiently small positive `beta` on those continuing branches,

`BR_i^{B3 prime} > 0 > BR_i^{G prime}`.

The claim is local and sufficient. It does not assert a primitive-space iff characterization, global equilibrium existence, uniqueness over the full strategy set, or the absence of off-regime deviations.

## Quantifiers and domains

- existence of `epsilon>0`, conditional on the stated regular beta-zero branches;
- `beta in (0,epsilon)` only;
- central smooth allocation regime only;
- local stationary-response branches only unless global best-response status is separately certified;
- `x_i in [0,1]`, `p_T>=0` are the primitive strategy domains, but this theorem does not certify global optimality over those domains.

## Independent reconstruction

PASS.

A clean-room SymPy derivation in `stage4a_v2/code/independent_adversarial_audit.py` reconstructs from the primitive beta-zero system:

`M_G(0) = -3 T alpha^2 kappa_L^2 / [16 Delta (Delta+T)^3] < 0`,

and, writing `delta=alpha beta`,

`d M_B3 / d delta |0 = alpha^2 d^3/(Delta_i^3 Delta_j^2) > 0`,

hence

`d M_B3 / d beta |0 = alpha^3 d^3/(Delta_i^3 Delta_j^2) > 0`.

This derivation does not import the production solver or exact-certificate code.

## Assumptions actually used

- beta-zero central-regime formulas and positive quality gaps;
- `T>0`, `Delta_i>0`, `d>0`;
- smooth participation and decision equations;
- nonsingular local Jacobians;
- strict relevant SOCs;
- continuity of the continuing branches.

## Boundary/corner/globality audit

The local derivative identities survive. However, the canonical finite-beta stationary configurations fail global public-best-response tests in a rival-public-entry regime. That failure does not algebraically contradict this theorem because the theorem is conditional/local. It does prohibit interpreting `BR` here as a globally optimal best-response branch without additional Stage-4 certification.

## Counterexample search

No counterexample to the two beta-zero identities was found. A separate finite-deviation counterexample was found to global-equilibrium promotion; see `stage4a_v2/COUNTEREXAMPLE_CANONICAL_FREE_RIDE.md`.

## Maximum defensible manuscript wording

A local stationary-response sign reversal exists on regular continuing branches under the stated conditions.

Prohibited without a repaired Stage 4: calling the certified canonical pair a global public Nash equilibrium or SPNE solely because the local FOCs/SOCs and exact stationary-root certificate pass.

## Certificate state

`PASS — LOCAL/CONDITIONAL MATHEMATICAL CLAIM ONLY`.

This PASS does not rescue the failed global equilibrium/SPNE claim.