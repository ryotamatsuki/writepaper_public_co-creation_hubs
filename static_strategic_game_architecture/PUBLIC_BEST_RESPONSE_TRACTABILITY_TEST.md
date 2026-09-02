# Public Best-Response Tractability Test

Without the metro wedge, H1/H2 indifference occurs at

`theta_PP=[x_1(c+u)-x_2 c]/[u(x_1+x_2)]`.

With an active metro interval, provider regions are separated by affine-in-`theta` boundaries. Regional welfare therefore reduces to integrals of affine functions over endogenous intervals plus quadratic investment cost.

The boundaries are rational functions of `x_1,x_2,p_T`; after substituting the closed-form interior `p_T^*`, Stage 4 can differentiate exact rational expressions using SymPy.

Risk: clipping and no-use corners create regimes. Stage 4 must solve an economically meaningful interior region first and then audit all boundary regimes; it may not report an FOC as a global equilibrium.

**Verdict: PASS WITH PIECEWISE-REGIME OBLIGATION.**