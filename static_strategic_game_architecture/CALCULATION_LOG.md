# Calculation Log

Stage 3 computations are pre-screens only.

## Scoring
Fixed weights total 53 per score point, maximum weighted score 265. Candidate totals: A 216; B 173; C 230; D 252; E 168; F 195.

## Metro interior demand
With
`A=(1/u)(1/x_1+1/x_2)`, `B=1+2c/u`,

`D_T=2[A(m-p)-B]`.

Private profit:
`Pi_T=2p[A(m-p)-B]`.

Exact FOC/SOC pre-check:
`p*=m/2-(2c+u)x_1x_2/[2(x_1+x_2)]`,
`d2Pi/dp2=-4A<0`.

Exact comparative statics:
`dp*/dx_1=-(2c+u)x_2^2/[2(x_1+x_2)^2]`, mirror for `x_2`.

## Public-public threshold
`theta_PP=[x_1(c+u)-x_2 c]/[u(x_1+x_2)]`.

All expressions must be re-derived and audited with regime restrictions in Stage 4. Script: `code/stage3_precheck_sympy.py`.