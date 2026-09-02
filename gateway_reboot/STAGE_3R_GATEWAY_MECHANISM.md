# Stage 3R — Gateway Mechanism Specification

## Verdict

**PASS.**

## Minimal environment

Two peripheral regions, each mass 1 of firms with `c~U[0,1]`, access one shared metropolitan/national innovation ecosystem. Region `i` chooses binary gateway establishment `E_i in {0,1}` at fixed cost `F`.

Let `B=A_M-delta`; `A_M` and direct-access friction `delta` are otherwise redundant.

Shared-core value is `B + beta N` net of normalized direct-access friction, where `N` is the mass effectively connected to the common ecosystem.

Available access-route bonuses are:

- direct: `0`;
- own gateway: `rho`;
- other-region gateway: `rho-tau`.

Canonical overlap region: `0<tau<rho`.

A firm's realized utility is the best available route value minus `c`, or outside option 0. The outside option is essential for an endogenous participation margin; it is not an additional strategic feature.

Regional objective is aggregate realized utility of origin-region firms minus its own fixed gateway cost.

## Stage 3R kill tests

1. Both forces arise from the same allocation problem: **YES**.
2. Gateway entry changes connected mass `N`: **YES**.
3. Other-region gateway substitutes for own gateway when `rho>tau`: **YES**.
4. Shared-core network feedback is the only positive strategic channel capable of overturning overlap: **YES; beta=0 gives strict substitutes**.
5. Gateway overlap is the negative strategic channel: **YES**.
6. Delete beta: complementarity disappears. Delete cross-region usability by setting `tau=rho`: overlap disappears and the remaining network channel is positive.
7. Redundant parameter: `delta` with `A_M`; normalized into `B`.
8. Fixed point is not pathological: with `beta<1/2`, the aggregate participation map is a contraction in the interior.
9. Fixed point unique in the canonical domain: **YES**.
10. Outside option 0 required: **YES**; without it `N=2` mechanically and thickening disappears.

## Mechanism

Gateway entry simultaneously:

- brings additional firms into the shared ecosystem, increasing `N` and hence the common value `beta N`;
- provides an overlapping access route that reduces the incremental value of establishing another regional gateway.

The Stage 4R object is `Gamma^G=D1-D0`, the change in own regional entry benefit induced by neighbor entry.