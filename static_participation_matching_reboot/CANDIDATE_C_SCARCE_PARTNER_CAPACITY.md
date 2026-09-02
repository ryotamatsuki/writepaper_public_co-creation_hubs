# Candidate C — Scarce-Partner Competition / Capacity-Constrained Matching

## Skeleton

Projects `(k,z)` choose a Hub. Each Hub controls a finite vector of heterogeneous partner capacities. Compatibility is bilateral: project type `k` can be served by multiple partner types `s` with different surplus, and each partner has limited capacity. The Hub must solve an assignment/flow problem across participating projects and partners.

## Participation

PASS in principle. Low `z` projects choose `0`; different `k` regions can select H1/H2/HT.

## Matching state

A genuine candidate is `Xi_h=(remaining compatible capacities, assignment composition, shadow values lambda_{hs})`. One project's allocation changes another project's feasible assignment.

## Label stripping

This is the only family where a truly bilateral assignment constraint can survive provider-index reduction. However, two cases arise:

1. Dedicated partner classes / random rationing make the model tractable, but then `mu_h(k)` is merely class-specific capacity congestion and the model reduces to a multi-class service facility.
2. Cross-compatible scarce partners preserve matching-specific shadow values, but provider choice changes the composition of an endogenous assignment problem. With two strategic public investments plus a private follower price, the full game becomes high-dimensional and heavily piecewise; analytical best responses are no longer minimally credible.

Matching-with-contracts theory also shows how strongly tractability/existence depends on substitutability-type restrictions; allowing complementarities or composition interactions raises the complexity further.

**Verdict: KILL FOR BASELINE — K10. Strongest candidate, but the matching-specific version is not minimally tractable.**