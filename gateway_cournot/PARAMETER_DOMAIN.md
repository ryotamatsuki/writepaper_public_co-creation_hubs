# Canonical Parameter Domain

Define

`m=2+beta`, `n=beta-1`, `b=2mn`.

The lowest possible engagement-cost coefficient across all gateway profiles is

`L=kappa-rho`.

A clean sufficient domain for all four firm subgames is

- `A=a-cbar>0`,
- `beta>=0`,
- `0<tau<rho<kappa`,
- `L > [2m^2+|b|]/9`,
- `F>0`.

Equivalently,

`9(kappa-rho)-2(2+beta)^2 > |2(2+beta)(beta-1)|`.

## Why this is sufficient

For every gateway profile, `K_i,K_j >= L`. Thus

`a_i=9K_i-2m^2 > |b|`,

`a_j=9K_j-2m^2 > |b|`.

Consequently:

1. each firm's engagement objective is strictly concave in own engagement because `K_i>2m^2/9`;
2. `Delta=a_i a_j-b^2>0`, so the simultaneous engagement FOCs have a unique solution;
3. `a_i+b>0` and `a_j+b>0`, hence `x_i^*,x_j^*>0` for `A>0`;
4. the engagement FOC implies `q_i^*=3K_i x_i^*/(2m)>0`;
5. no explosive common-network feedback occurs inside the specified domain.

## Price and cost levels

The strategic system depends on `a` and `cbar` only through `A=a-cbar`. Absolute price and marginal-cost levels can therefore be shifted upward together without changing `(x,q,pi,Gamma)`. For every admissible strategic solution, choose a sufficiently large `cbar` and set `a=cbar+A`; then marginal costs and equilibrium price are nonnegative. This level normalization is not a strategic parameter and is not used to generate sign reversals.

## beta=0 deletion domain

When `beta=0`, the stability condition becomes

`kappa-rho > 4/3`.

This domain is used in the exact proof that the regional-welfare gateway interaction is strictly substitutable when the shared-core network effect is deleted.