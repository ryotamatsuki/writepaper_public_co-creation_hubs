# Stage 7.5A v2.1 — T3 Symmetry Quantifier Limited Reopen

Date: 2026-09-10 JST

Trigger: Stage 11 independent referee regression after Stage-10 merge `02bcd37ac1ab019dc6df2c0d29c8f312a1c992a4`.

Original Stage-7.5A authority: `reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md`.

Original Stage-8 freeze: `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md` @ merge `ad927ca783a6123ea4fc6f55f65598ebd6ab583b`.

## 1. Certification regression

**CERTIFICATION REGRESSION — T3 G-branch anchor was under-specified in the prose quantifier ledger.**

The proved beta-zero G sign used by T3 is T2:

`M_G(0) = -3 T alpha^2 kappa_L^2 / [16 Delta (Delta+T)^3] < 0`,

and T2 is established at a **symmetric regular beta-zero full-game stationary state**. The original T3 ledger said only that there exist continuing regular B3 and G stationary branches. Read literally, that statement permits a G branch through an asymmetric beta-zero state for which the displayed T2 sign formula was not proved.

The Stage-11 reviewer independently re-derived T1 and T2 from primitives and found no mathematical error in either formula. The defect is a missing anchor restriction in the T3 quantifier, not a failure of the mechanism or derivation.

## 2. Repaired exact T3 quantifier

The maximum certified T3 statement is:

> Suppose there exist (i) a regular interior B3 stationary branch through `beta=0`, and (ii) a regular interior G stationary branch through a **symmetric regular beta-zero full-game stationary state**. Assume the relevant participation and stationary-condition Jacobians are nonsingular, the private stationary price and both public stationary choices satisfy strict second-order conditions, and both branches continue smoothly for sufficiently small positive `beta`. Let B3 use the matched full-game price along the corresponding G branch. Then there exists `epsilon>0` such that, for every `beta in (0,epsilon)` on these continuing branches, `BR'_{B3}>0>BR'_G`.

Equivalently, choose local neighborhoods for B3 positivity, G negativity, strict SOCs, regular allocation ordering, nonsingular continuation, and matched-price continuity, and take their finite minimum radius. This yields a common positive `epsilon`.

## 3. What is unchanged

No change to:

- players, timing, strategy sets, utilities, objectives, route set, or welfare definition;
- T1 formula or scope;
- T2 formula or scope;
- the repaired primitive vector or numerical witness;
- Stage-4A global-deviation certification;
- B3 benchmark definition;
- W1/W2/W3;
- novelty envelope;
- robustness/generality exclusions.

No new theorem is added. The T3 domain is **narrowed to the state at which its negative-G anchor is actually proved**.

## 4. Earliest affected stage

Earliest affected stage: **Stage 7.5A — theorem quantifier certification**.

Stage 4/4A need not reopen because no equilibrium computation, continuation certificate, or global witness changes. Stage 7 need not reopen because no welfare object changes. Stage 6 need not reopen because novelty is not enlarged.

## 5. Downstream repair contract

The following downstream reader-facing and validation objects must be synchronized:

1. `sections/03_main_results.tex` theorem statement;
2. `paper/main.tex` abstract;
3. `sections/08_introduction.tex`;
4. `sections/05_robustness.tex` and `sections/09_discussion.tex` where T3 scope is summarized;
5. `sections/10_conclusion.tex`;
6. `scripts/stage75a_scope_audit.py` and `scripts/verify_freeze.py`;
7. the Stage-8 freeze amendment below.

All Stage-9/10 verification gates must be rerun after this bounded correction.

## 6. Corrected certificate

T1: PASS unchanged.

T2: PASS unchanged, symmetric beta-zero G state.

T3: **PASS AFTER QUANTIFIER NARROWING** to a G branch through the symmetric beta-zero state covered by T2.

This correction does not license any asymmetric-G theorem, global equilibrium theorem, uniqueness result, explicit `epsilon`, or primitive-space characterization.

## 7. Verdict

**LIMITED STAGE 7.5A REOPEN CLOSED — QUANTIFIER REPAIRED.**

Return to Stage 11 only after downstream manuscript and CI synchronization succeeds.
