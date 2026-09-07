# Permanent Counterexample — Rival-Public Free-Riding Deviation

Date: 2026-09-07

Primitive vector:

`v=.1, alpha=.5, beta=.05, rho=.15, rhoT=.05, kappa_L=.27, kappa_T=.02, tau=.05, gamma=.9`.

This is the paper's canonical finite-decimal vector.

## Why the old Stage-4 search could miss this

The old `public_two_sided_platform_hard_kill/code/numerical_hard_kill.py` solves only the smooth `0 -> H_T -> H_i` central regime. Its `profit` routine assigns `-1e9` whenever that central-regime solve is invalid. Thus a price/action state that moves into rival-public entry, public inactivity, clipping, or another valid regime is not re-solved from primitive route choice. Under workflow v2.0 this is fail-open/invalid continuation handling and cannot certify SPNE or global Nash behavior.

## Full game G

Canonical stationary configuration:

- `x_1=x_2 ~= 0.684028196361275`;
- independent all-regime private-price re-solve `p_T ~= 0.02274665143`;
- independent direct `W_1 ~= 0.22301946658`.

Finite deviation:

- hold `x_2 ~= 0.684028196361275`;
- set `x_1=0.1875`;
- independently re-solve downstream participation from the upper envelope of all four project routes and re-optimize the private price on its full positive-profit interval;
- obtain `p_T ~= 0.02518788380`;
- obtain project masses `(n_1^F,n_2^F,n_T^F) ~= (0,0.734535,0.622508)`;
- obtain `W_1 ~= 0.25444926939`.

Gain:

`Delta W_1 ~= +0.03142980281 > 0`.

Thus the canonical G stationary pair is not a Nash equilibrium of the stated public strategy game, and therefore cannot support the manuscript's global SPNE language.

## Matched-price B3

Canonical stationary configuration:

- `x_1=x_2 ~= 0.656020390747393`;
- `p_T=0.02274665145` fixed;
- independent direct `W_1 ~= 0.21091151931`.

The same deviation `x_1=0.1875`, holding `x_2` and `p_T` fixed, gives:

- project masses `(n_1^F,n_2^F,n_T^F) ~= (0,0.657340,0.747685)`;
- `W_1 ~= 0.23886596269`;
- gain `~= +0.02795444338 > 0`.

So the canonical B3 stationary pair is also not a global public Nash equilibrium.

## Economic reason

At `x_1=0.1875` with `n_1^F=0`,

`b_1=rho+x_1=0.3375`,

and therefore

`q_1=v+alpha b_1=.1+.5(.3375)=0.26875 < kappa_L=.27`.

Hence `H_1` gives negative utility even to `z=1`, so it is inactive. Region 1 free-rides on the stronger rival public hub for high project types. While `H_1` is inactive, changing `x_1` does not affect project routing or the private pricing problem; its direct `x_1` contribution is

`(rho+x_1)^2/4 - gamma x_1^2/2`.

Its derivative is

`(rho+x_1)/2 - gamma x_1`,

which is zero at

`x_1=rho/(2gamma-1)=.15/.8=.1875`,

with second derivative `1/2-gamma=-0.4<0`. The own-hub inactivity condition holds up to `x_1<=0.19`, so the free-riding optimum lies strictly inside that regime.

## Reproduction

Run:

`python stage4a_v2/code/independent_adversarial_audit.py`

The evaluator is independent of the production central-regime solver and preserves this case as a permanent Stage-4 regression test.