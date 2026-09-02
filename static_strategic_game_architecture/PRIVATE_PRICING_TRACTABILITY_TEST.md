# Private Pricing Tractability Test

In the unclipped interior metro-demand region define

`A=(1/u)(1/x_1+1/x_2)`, `B=1+2c/u`.

Then

`D_T=2[A(m-p_T)-B]`,
`Pi_T=2 p_T[A(m-p_T)-B]`.

FOC gives

`p_T^*=m/2-(2c+u)x_1x_2/[2(x_1+x_2)]`.

SOC is `-4A<0`.

Exact derivatives:

`dp_T*/dx_1=-(2c+u)x_2^2/[2(x_1+x_2)^2]<0`,
`dp_T*/dx_2=-(2c+u)x_1^2/[2(x_1+x_2)^2]<0`.

Required Stage 4 checks: `p_T^*>=0`, metro participation `m-p_T^*>=kappa`, `0<L_T<R_T<1`, clipping/corners, and whether price response changes public BR/welfare relative to B3.

**Verdict: PASS FOR STAGE-3 TRACTABILITY.**