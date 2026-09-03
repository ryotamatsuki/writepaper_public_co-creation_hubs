# Parameter Domain

Baseline primitives: `v>0, alpha>0, beta>=0, rho>=0, rhoT>=0, gamma>0, kappa_L>0, kappa_T>=0, tau>=0`.

Interior partner block requires each of `rho+x_i+beta n_i^F` and `rhoT+beta n_T^F` to lie strictly in `(0,1)`.

The main smooth regime uses, for each home region i,
`0 < s < t_i < 1`, `q_i>q_T`, and `kappa_T+p_T<kappa_L`, where
`s=(kappa_T+p_T)/q_T`,
`t_i=(kappa_L-kappa_T-p_T)/(q_i-q_T)`.

In this regime the ordering is `0 -> H_T -> H_i`; the rival public platform is locally dominated at a symmetric point when `tau>0` but remains feasible and becomes active off symmetry when its network-quality advantage exceeds the real access-cost wedge.