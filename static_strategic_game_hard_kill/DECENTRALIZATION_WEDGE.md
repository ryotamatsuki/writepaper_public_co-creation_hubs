# Decentralization Wedge — P5

## B2

In the full-public symmetric interior benchmark,

`x_N^{B2}-x_SP^{B2}=(8c^2+4cu-u^2)/(8 gamma u)`.

Hence B2 underinvests when `c/u < (sqrt(3)-1)/4` and overinvests above that threshold, subject to interior/KKT conditions.

## G: local planner marginal at the symmetric interior Nash

Let `a=2c+u`, `B=12c^2+28cu+15u^2`, and `y=m/x`.

Using the full-game Nash FOC to eliminate `gamma`, the planner marginal at the decentralized Nash point is

`F_SP^G(x_N)=[B x^2-4amx-20m^2]/(32u x^2)`.

Its sign is the sign of

`psi(y)=B-4ay-20y^2`.

The active metro regime implies `a/2<y<a/2+u`. `psi` is strictly decreasing. At the lower and upper boundaries,

`psi(a/2)=8(u^2-2c^2)`,

`psi(a/2+u)=-4(2c+3u)^2`.

Therefore, if `c<u/sqrt(2)`, there is a unique interior threshold

`y*=[-a+sqrt(a^2+5B)]/10`

such that the planner locally wants **more** public investment at the Nash point for `y<y*` and **less** for `y>y*`. If `c>=u/sqrt(2)`, the planner marginal is non-positive throughout the active regime.

This is a nontrivial welfare threshold, but it is not sufficient to rescue Stage 4: it is derived entirely from the characteristic-space representation, and B2 already contains a sign-changing decentralization wedge.

**P5 classification: CONDITIONAL/TRUE as a local marginal-wedge theorem; contribution KILL.**