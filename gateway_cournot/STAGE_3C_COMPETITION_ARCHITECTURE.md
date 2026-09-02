# Stage 3C — Endogenous Product-Market Competition Architecture

## Verdict

**PASS.**

The model closes as a three-stage subgame-perfect game without unauthorized primitives.

## Timing

1. Regional governments simultaneously choose `E_i in {0,1}` and pay `F E_i`.
2. After observing gateways, firms simultaneously choose ecosystem engagement `x_i >= 0`.
3. Firms simultaneously choose homogeneous Cournot quantities `q_i >= 0`.

## Gateway technology

Firm `i` pays engagement cost

`(K_i(E_i,E_j)/2) x_i^2`,

with

- `K_i(0,0)=kappa`,
- `K_i(1,0)=kappa-rho`,
- `K_i(0,1)=kappa-rho+tau`,
- `K_i(1,1)=kappa-rho`,

and `0<tau<rho<kappa`.

Thus the other region's gateway is useful but weaker than an own-region gateway.

## Shared-core productivity primitive

`c_i = cbar - x_i - beta(x_i+x_j)`.

The private component is `x_i`. The common-core component is `beta(x_i+x_j)`. Rival engagement therefore simultaneously thickens the common core and strengthens a product-market rival.

## Product market

`P=a-q_1-q_2`.

Firm profit is

`pi_i = [P-c_i]q_i - (K_i/2)x_i^2`.

All firm profits, market shares, outputs and price are equilibrium objects.

## Regional-government baseline

The canonical regional-welfare objective is

`W_i = pi_i + CS/2 - F E_i`,

where `CS=Q^2/2` in the homogeneous linear-demand market and the two symmetric regions are assigned equal consumer shares.

The producer-only objective `pi_i-FE_i` is retained only as an instructed comparison.

## Kill-question answers

1. Firm profit endogenous: YES.
2. Gateway acts only through engagement cost: YES.
3. Engagement changes marginal cost: YES.
4. Marginal cost changes Cournot `q`, `P`, market share and profit: YES.
5. Rival engagement creates genuine competition: YES.
6. Shared core benefits both firms when `beta>0`: YES.
7. Common-core and rival-strengthening forces are separable by `beta=0` and by the sign of `beta-1`: YES.
8. Strategic profitability enters only through `A=a-cbar`; absolute price/cost level is otherwise redundant for the strategic system.
9. Engagement equilibrium is unique under the stated sufficient stability domain.
10. Unbounded engagement is excluded by strict concavity and joint-stability restrictions.
11. Separate direct participation is not needed in this Cournot architecture.
12. Consumer surplus is not needed for the existence of a sign reversal, but it changes mechanism identification materially.