# Preliminary theorem certificate — repaired global sign-reversal witness

Status: **Stage 4 construction PASS; Stage 4A certification pending**

Workflow authority: v2.1 @ `fe6fa5f2a94632be578d5e9fab39e7df86013113`.

## Object

Same public two-sided-platform game and same primitive definitions as the existing manuscript. Selected existing-primitive witness:

`beta=.01, gamma=.825, tau=.35`, with all other primitive values unchanged.

## Candidate G statement

At the symmetric candidate `x_1=x_2=x_G ~= .8371022382`, after the private intermediary observes `(x_1,x_2)` and globally chooses its access fee, the construction solver obtains a regular continuation with `p_G ~= .0184036832`. The local public best-response slope is negative, approximately `-.01986145`.

Dense all-regime unilateral-deviation search over `x_i in [0,1]`, re-solving the private continuation, detects no profitable public deviation above numerical tolerance; the best detected point is `x_i ~= .8370893604` with payoff difference `~9.4e-12` relative to the candidate.

Certificate maturity: **numerically supported global-equilibrium candidate, not yet independently certified**.

## Candidate B3 statement

Holding the private fee fixed at the matched all-regime G price, the symmetric candidate is `x_1=x_2=x_B ~= .8258903860`. The local public best-response slope is positive, approximately `+.00396813`.

Dense all-regime unilateral-deviation search over `x_i in [0,1]` detects no profitable public deviation above numerical tolerance; the best detected point is `x_i ~= .8258866749`.

Certificate maturity: **numerically supported global-Nash candidate, not yet independently certified**.

## Headline comparative object

At this witness the construction evidence supports

`BR_G' < 0 < BR_B3'`.

This is the repaired complements-to-substitutes comparison. It is not yet authorized as a certified theorem because Stage 4A has not independently verified globality, continuation completeness, numerical stability, and the local derivative signs.

## Domain and quantifiers

- public strategies: `x_i in [0,1]`;
- private price: nonnegative, globally optimized over the economically relevant compact interval implied by primitive demand bounds;
- project routes: full primitive choice set `0,H1,H2,HT`;
- partner participation: primitive `c~U[0,1]` with boundary clipping;
- equilibrium concept targeted: pure-strategy SPNE for G and public Nash equilibrium for B3;
- result type at Stage 4: **existence witness / construction-level**, not a primitive necessary-and-sufficient characterization and not a global parameter-space theorem.

## Required Stage 4A attacks

1. reconstruct allocation/payoffs without importing Stage 4 repair code;
2. independently optimize private price after material public deviations;
3. attack the full public strategy interval, boundaries, kinks, active-set changes and low-investment regions;
4. test multiple participation/private continuation candidates;
5. independently approximate/verify local own and cross derivatives around G/B3 candidates;
6. fail closed on unresolved histories;
7. retain any counterexample permanently.
