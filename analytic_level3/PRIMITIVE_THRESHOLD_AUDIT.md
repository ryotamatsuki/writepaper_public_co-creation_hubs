# Primitive Threshold Audit

## Objective
The Level-3 target requires a primitive or transparent primitive-index characterization of the complete headline comparison, not merely a condition in endogenous Hessian blocks.

## Exact small-beta analytic fallback
A useful primitive-direction result can nevertheless be derived around `beta=0`.

Let `delta=alpha beta`, `Delta_i=A_i-A_T>0`, and keep the central fixed-price regime. At `delta=0`, B3 has `M_B3=0` identically. A first-order expansion of the exact participation system and regional objective gives

`partial M_B3 / partial delta |_{0} = alpha^2 d^3/(Delta_i^3 Delta_j^2) > 0`,

where `d=kappa_L-kappa_T-p>0`.

Since `delta=alpha beta`,

`partial M_B3 / partial beta |_{0} = alpha^3 d^3/(Delta_i^3 Delta_j^2)>0`.

Because `M_B3(beta=0,x,p)=0` identically, the derivatives of the B3 equilibrium path and of the matched price do not enter the first derivative at zero: their chain-rule multipliers are zero.

For the full game at `beta=0`, write `T=A_T=v+alpha rho_T` and, at a symmetric regular state, `Delta=A_i-A_T>0`. Substituting the exact private optimum yields

`M_G(0) = -3 T alpha^2 kappa_L^2/[16 Delta(Delta+T)^3] < 0`.

Hence, if the B3 and G equilibrium branches at `beta=0` are regular interior, satisfy the relevant private/public SOCs, and continue smoothly in `beta`, continuity implies the existence of an `epsilon>0` such that

`0<beta<epsilon  =>  BR_i^{B3′}>0>BR_i^{G′}`.

This is a genuine analytic **sufficient-condition theorem**. It is Level 2 because `epsilon` is implicit and the regular beta-zero equilibrium/SOC conditions have not been eliminated into a finite primitive threshold.

## Why the full Level-3 primitive threshold fails the hard kill

The headline comparison requires two distinct equilibrium roots:

1. `G`: solve the public equilibrium jointly with the endogenous private price policy;
2. `B3`: set `bar p=p^G` and solve a different public equilibrium.

Exploratory symbolic elimination on the frozen model produced the following complexity before reaching the required full cross derivative:

- symmetric private FOC after participation reconstruction: rational numerator total degree 8 in `(t,q_T,p)`, 36 monomials;
- symmetric B3 public FOC numerator: total degree 9, 58 monomials;
- general B3 cross derivative before aggressive elimination: roughly 26k SymPy operations;
- general private-policy first derivative after implicit reduction: roughly 115k operations;
- G public FOC with the exact first price-policy response substituted: roughly 124k operations, before differentiating again to obtain the G cross derivative and `p_{ij}`.

Aggressive symbolic `cancel/factor` of the second-derivative objects did not complete within the bounded exploratory simplification window. These counts are diagnostic instrumentation, not theorem statements, but they show that direct root elimination is not an economically transparent route.

## Threshold search result

No unique closed-form threshold such as `beta>beta^dagger(...)`, `gamma>gamma^dagger(...)`, or a private-elasticity threshold was found that is necessary and sufficient over all of `Theta^{RI,SOC}`.

The only exact one-dimensional statement obtained is local:

`0<beta<epsilon`

around a regular `beta=0` pair of equilibrium branches. It is sufficient, not necessary, and `epsilon` is not primitive-explicit.

## Classification

- Level-3A: FAIL.
- Level-3B: FAIL.
- Level-3C derivative-block equivalence: PASS.
- Strong Level-2 small-beta sufficient theorem: PASS.

This supports `NO-GO LEVEL 3 / GO LEVEL 2`.