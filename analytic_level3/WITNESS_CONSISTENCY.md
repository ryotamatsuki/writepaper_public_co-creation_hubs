# Canonical Witness Consistency

The L3 track does not alter the frozen witness. It reconciles new identities against the production solver.

Canonical generated values:

- `G x = 0.6840281531716518`;
- `G p_T = 0.022746652886822923`;
- `G own second = -0.19938895349606398`;
- `G cross second = -0.004809167261951883`;
- `G BR slope = -0.024119527073232862`;
- `B3 x = 0.6560164129736359` in the Stage-4 strategic solver;
- `B3 own second = -0.11835102280915384`;
- `B3 cross second = 0.013543617462098548`;
- `B3 BR slope = 0.11443599844454423`;
- generated private price response `-0.016471065568401944` (frozen display `-0.016505`).

## Explicit participation branch at G

At the witness, the quadratic participation reduction has a strictly positive discriminant and two algebraic roots. The upper root reproduces the production `q_T` to numerical tolerance. The lower root implies `s=(kappa_T+p)/q_T>1` and is therefore excluded by the central-regime inequalities rather than by arbitrary root selection.

The symmetric Jacobian factorization `det J=uH` is strictly away from zero at the witness (approximately `0.11338`), confirming the same regular branch used by the production IFT argument.

## Verification status

`analytic_level3/code/verify_witness.py` independently loads the frozen Stage-4 solver and checks:

- strict B3/G SOC margins;
- opposite cross derivatives and BR signs;
- private price response against the frozen value;
- participation quadratic/root selection;
- nonsingular symmetric participation Jacobian.

This is a high-precision reconciliation of the new algebra with the production numerical certificate. It is **not** represented as an interval-arithmetic proof of a primitive Level-3 region.

## Relation to the small-beta theorem

The analytic small-beta theorem proves that reversal occurs on a nonempty neighborhood of any regular beta-zero B3/G branch satisfying its hypotheses. It does not prove that the canonical `beta=0.05` witness lies below an explicit primitive `epsilon`. The canonical witness therefore remains the production numerical existence certificate; no proof-status inflation is made.