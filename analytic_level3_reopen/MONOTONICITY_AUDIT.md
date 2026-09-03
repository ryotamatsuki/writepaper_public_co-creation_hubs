# Monotonicity Audit

## `beta`

L3-0's exact small-`beta` result survives unchanged: B3 starts at zero cross
effect with a strictly positive first derivative, while G has a strictly
negative cross effect at the symmetric beta-zero state. This proves a local
sufficient reversal interval on continuing regular/SOC branches.

L3-1 did not find a global sign-stable derivative of either complete cross
effect with respect to `beta`. The equilibrium-manifold reduced polynomials have
mixed coefficients, and their implicit root-response terms are not sign-fixed
by the central regime alone.

## Parameter inversion

A stronger structural identity is

`delta = q_T(q_T-A_T)/(2(t q_T-a))`.

Thus the network parameter can be solved from the symmetric participation
manifold rather than solving the equilibrium variables from `delta`. This is a
promising parametric-boundary route, but monotonicity of the mapping over the
full regular branch was not certified.

## Other primitives

No rigorous global single-crossing theorem was found for `alpha`, `gamma`,
`rho_T`, or the cost wedge. Numerical monotonicity was not treated as proof.

## Verdict

- monotonicity threshold: **not obtained**;
- single-crossing threshold: **not obtained**;
- parameter inversion: **exactly obtained**;
- unique one-dimensional primitive threshold: **not established**.
