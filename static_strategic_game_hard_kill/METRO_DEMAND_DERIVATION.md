# Metro Demand Derivation

There are two identical unit masses of projects (one per home region). In the central public–metro–public regime,

`theta_1T=(x_1 d-q)/(x_1 u)`,
`theta_2T=(q-x_2 c)/(x_2 u)`.

The metro share in either home region is

`ell_T = theta_2T-theta_1T`

`= (x_1+x_2)(q-K)/(u x_1 x_2)`.

Hence total metro demand is

`D_T=2 ell_T`.

The factor two from the two home-region populations affects profit levels but not the private argmax.

The exact global expression is

`D_T(q)=2 max{0, clip(theta_2T,0,1)-clip(theta_1T,0,1)}` for `q>kappa`,

and `D_T=0` for `q<kappa` (up to the stated tie rule at equality).

Thus demand is continuous, nondecreasing, and piecewise linear in the metro gross offer `q=m-p_T`.