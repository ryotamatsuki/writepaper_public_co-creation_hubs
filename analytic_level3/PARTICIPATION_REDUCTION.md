# Participation Reduction

Inside the frozen central regime define

`delta=alpha beta`, `a=kappa_T+p`, `d=kappa_L-a`, `A_i=v+alpha(rho+x_i)`, `A_T=v+alpha rho_T`.

The original three-equation fixed point is

`F_i=t_i[A_i+delta(1-t_i)-q_T]-d=0`, `i=1,2`,

`F_T=q_T-A_T-delta[t_1+t_2-2a/q_T]=0`.

## Exact reduction

Multiplying the third equation by positive `q_T` yields

`q_T^2-[A_T+delta(t_1+t_2)]q_T+2 delta a=0`.

Thus, for `delta>0`,

`q_T^{±}={B ± sqrt(B^2-8 delta a)}/2`,

where `B=A_T+delta(t_1+t_2)`.

Given an admissible root and cutoffs, the public investments are recovered exactly from the first two equations:

`x_i=[q_T+d/t_i-delta(1-t_i)-v]/alpha-rho`.

Hence the central participation block is not intrinsically a three-dimensional numerical black box. It can be parameterized by `(t_1,t_2,p)` plus an explicitly selected quadratic branch.

## Root discipline

Both algebraic roots exist whenever the discriminant is positive. `Theta^RI` selects a root only by the economic inequalities `q_T>0`, `0<s<t_i<1`, positive masses, partner interiority, and the sorting-order restrictions. At the canonical witness the upper root reproduces the numerical `q_T`; the lower root violates central-regime project participation because `s=a/q_T>1`.

## Symmetric Jacobian factorization

At `t_1=t_2=t`, define

`u=q_i-q_T-delta t=d/t-delta t`,

`w=1-2 delta a/q_T^2`,

`H=u w-2 delta t`.

The participation Jacobian satisfies

`det J = u H`.

For an `x_1` perturbation on the symmetric branch,

`dt_1/dx_1=-(alpha t/2)(w/H+1/u)`,

`dt_2/dx_1=-(alpha t/2)(w/H-1/u)`,

`dq_T/dx_1=-alpha delta t/H`.

These expressions immediately show why cross-region participation effects vanish when `beta=0` (`delta=0`).

## Tractability assessment

**Gate 1: PASS.** The regular central participation system has an exact finite representation and a transparent root branch. This is substantially stronger than the production numerical solver, but it is valid only within the declared regime.