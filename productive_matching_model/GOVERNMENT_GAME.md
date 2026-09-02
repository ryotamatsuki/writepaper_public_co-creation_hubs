# Decentralized Government Hub Game

## Regional objective

For profile `(E_i,E_j)`, regional government `i` maximizes

`W_i = R_i(mu_i,mu_j)-F E_i`,

where

`R_i(u,v)=[8A^2+20A Delta u-4A Delta v-14 Delta^2 u v+17 Delta^2 u+5 Delta^2 v]/36`.

The fixed cost affects the level of the entry decision but not strategic complementarity/substitutability because it cancels from `Gamma=D_1-D_0`.

## Gross entry thresholds

`D_0 = R_i(h,p)-R_i(0,0)`

is the gross gain from opening when the neighboring hub is absent.

`D_1 = R_i(t,t)-R_i(p,h)`

is the gross gain from opening when the neighboring hub is already present.

The binary hub game is:

- if `D_1<D_0` (strategic substitutes), intermediate fixed costs can support one-hub asymmetric equilibria;
- if `D_0<D_1` (strategic complements), intermediate fixed costs can support coordination multiplicity between no-hub and two-hub outcomes;
- exact equilibrium regions also depend on whether the gross thresholds are positive.

This topology is standard and is not a novelty claim.

## Public-decentralization wedge

Let `X_0` be the welfare effect on region `j` when region `i` opens from `00` to `10`. Exact calculation gives

`X_0 = Delta p[12A+4A p+Delta(27-33p+14p^2)]/36`.

The quadratic `27-33p+14p^2` is strictly positive because its discriminant is negative. Hence

`X_0>0`

throughout the admissible domain.

Thus the first hub creates a positive cross-region externality that the opening government does not internalize.

When the other hub already exists, define

`X_1 = R_j(t,t)-R_j(h,p)`.

`X_1` is not globally signed: referral can improve the neighboring firm's own match prospects, while the focal firm's strengthening changes downstream competition. Therefore the social marginal value of the second hub can lie above or below the decentralized regional marginal value.

## Planner comparison

Because `W_1+W_2` equals total producer surplus plus consumer surplus minus total hub costs, a single planner internalizes both regional changes.

Planner gross marginal benefits are

`S_0=D_0+X_0`,

`S_1=D_1+X_1`.

Decentralization is therefore theorem-relevant: the regional government uses `D_e`, not `S_e`.

## Classification

**Public decentralization has theorem-level content.** It creates uninternalized cross-region matching and product-market effects; it is not a relabeling of a privately priced intermediary.