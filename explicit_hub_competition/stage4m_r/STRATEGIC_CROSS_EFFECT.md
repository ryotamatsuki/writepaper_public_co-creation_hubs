# Strategic Cross-Effect

Define

`Gamma=(W_i(1,1)-W_i(0,1))-(W_i(1,0)-W_i(0,0))`.

Fixed cost cancels. With referral,

`Gamma = K_B-2K_S+K_0 + sigma(mu_B-2mu_S)`.

Let

`dmu=mu_B-2mu_S = -p[(1-p)l+p(2-p)m] < 0`

and

`dnu=nu_B-2nu_S = -p[(1-p)l^2+p(2-p)m^2] < 0`.

Then

`Gamma=(4a/9+sigma)dmu+(11/18)dnu-(7/18)(mu_B^2-2mu_S^2)`.

## Exact global sign proof

Because

`dGamma/da=(4/9)dmu<0`
and
`dGamma/dsigma=dmu<0`,

the largest possible `Gamma` is approached at `a=h` and `sigma=0`.

Normalize `h=1`, let `Q=1-p`, `X=m/h`, `Y=l/h`, and impose `0<Y<X<1` by `X=Y+(1-Y)R`, `R in (0,1)`.

For referral define

`P=-36 Gamma/(p h^2)`

at the maximizing boundary. Before the `X=Y+(1-Y)R` substitution,

`P = -14Q^5X^2+14Q^4X^2-28Q^4XY+28Q^3XY-28Q^3X+14Q^3Y^2-22Q^2X^2+28Q^2XY+12Q^2X-14Q^2Y^2+7QX^2-28QXY+14QX+22QY^2+16QY-7Q+15X^2+2X+7`.

After mapping to the unit cube `(Q,R,Y)`, its exact tensor-product Bernstein representation has degree `(5,2,2)`, **no negative coefficient**, 3 zero coefficients, and minimum strictly positive coefficient `7/5`. Since every Bernstein basis function is strictly positive in the interior and at least one coefficient is positive, `P>0`.

Therefore

**`Gamma<0` globally for every admissible parameter vector.**

## No-referral benchmark

With `mu_B0=ph+p(1-p)l`,

`mu_B0-2mu_S=-p[(1-p)l+m]<0`.

The analogous exact Bernstein certificate has degree `(3,2,2)`, no negative coefficient, one zero coefficient, and minimum positive coefficient `7/3`.

Hence **no-referral `Gamma_0<0` globally** as well.

Referral cannot generate strategic complementarity.
