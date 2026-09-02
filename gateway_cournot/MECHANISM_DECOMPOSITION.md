# Mechanism Decomposition and Hard Deletion Tests

The homogeneous-Cournot model has three economically distinct forces:

1. **Shared-core thickening.** Higher `x_j` enters firm `i`'s cost through `-beta x_j`.
2. **Gateway overlap.** A neighbor gateway lowers firm `i`'s engagement-cost coefficient from `kappa` to `kappa-rho+tau`, reducing the incremental value of an own gateway.
3. **Rival strengthening / product-market competition.** Higher `x_j` lowers rival `j`'s marginal cost and changes Cournot market shares.

Holding engagement choices fixed,

`dq_i/dx_j=(beta-1)/3`.

Thus the direct common-core effect and rival-strengthening effect cancel at `beta=1`. This is a firm-level diagnostic, not the local-government threshold.

## Exact network-deletion theorem for the regional-welfare baseline

Set `beta=0`. Write

`ell=kappa-rho`, `d=rho-tau>0`, `t=tau>0`.

The stability domain is `ell>4/3`. Define `y=ell-4/3>0`.

The exact regional-welfare cross effect simplifies to

`Gamma_W^C(0) = -8 A^2 R / [3(9y+8)(9d+9t+9y+8)(9ty+4t+9y^2+8y)^2]`,

where

`R = 189 d t^2 y +120 d t^2 +243 d t y^3 +432 d t y^2 +192 d t y +243 d y^4 +432 d y^3 +192 d y^2 +189 t^3 y +120 t^3 +189 t^2 y^2 +240 t^2 y +64 t^2`.

Every term in `R` is strictly positive for `d,t,y>0`; every denominator factor is positive. Therefore

**`Gamma_W^C(0)<0` throughout the canonical stable domain.**

Hence shared-core network feedback is necessary for gateway complementarity under the canonical regional-welfare objective.

## Exact positive and negative open-set witnesses

Use `A=1`, `kappa=5`, `rho=3/2`, `tau=6/5`.

At `beta=1/2`,

`Gamma_W^C = -27845669343 / 2525102237378 < 0`.

At `beta=6/5`,

`Gamma_W^C = 7916917186023346560 / 126835784788527792481 > 0`.

Both satisfy the strict stability domain. By continuity, both signs hold on nonempty open neighborhoods. Numerical local-perturbation tests confirmed the negative sign under +/-5% perturbations of every primitive and the positive sign under +/-2% perturbations of every primitive.

## Gateway-overlap limits

- `tau -> 0`: the neighbor gateway becomes maximally useful to the focal region, so duplication/overlap is strongest. In numerical admissible stress tests this pushes `Gamma_W^C` toward substitution; no global analytic sign theorem is claimed for all `beta` at this limit.
- `tau -> rho`: the neighbor gateway gives essentially no direct access-cost reduction to the focal firm. The overlap channel disappears, leaving the network-versus-rivalry interaction. Both signs remain possible depending on network strength and engagement feedback.

## Producer-only objective warning

For `W_i^P=pi_i-FE_i`, both signs also occur. But network feedback is not necessary for complementarity: at

`A=1, beta=0, kappa=3, rho=1/2, tau=1/10`,

`Gamma_P^C = 454818784/345792217681 > 0`.

This is strategic rent shifting/business stealing, not shared-core complementarity. Therefore the producer-only objective changes the economic interpretation qualitatively.