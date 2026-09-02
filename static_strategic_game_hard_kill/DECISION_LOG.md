# Decision Log

## Starting state

- canonical main: `3e525ecedfaa4fbac56652d7a40c3f6adf187d7e`
- PR #23: merged
- PR #24: merged
- open PRs at start: 0
- Stage 4 branch: `co-creation-theory/stage4r-mc-g-static-hard-kill`

## Gate decisions

### Gate 1 — Solvability
**PASS (key regimes).** The project partitions have an exact clipped-envelope representation; private price is piecewise quadratic and the central regime is closed form. Public objectives are analytically differentiable within regimes.

### Gate 2 — Matching specificity
**FAIL.** Exact transformation to a one-dimensional characteristic-space/Hotelling-type payoff system leaves all strategic objects unchanged.

Under the frozen decision tree, this is already dispositive:

**NO-GO — MATCHING STRUCTURE COLLAPSES.**

### Gates 3–6 — diagnostic only after fatal Gate 2

- private pricing is strategically relevant: PASS diagnostically;
- second public player is essential to the cross-government channel: PASS diagnostically;
- a triadic cross-effect exists: PASS algebraically;
- welfare thresholds exist: PASS locally/conditionally;
- none rescues the matching-specific contribution.

## Secondary structural failure

When metro is active, `q>kappa` implies full participation. Hence the static model cannot jointly represent active metro competition and extensive `0 -> H_i` additionality.

## Final routing

Stop the CSDS static route. Preserve dynamics only as a separately motivated future research question, not an automatic continuation.