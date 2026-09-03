# Generic Distribution Test

Replace `z~U[0,1]` only for the generality audit by a `C^1` CDF `F` with density `f>0` at the relevant cutoffs. In the central regime,

`n_i^F=1-F(t_i)`,
`n_T^F=F(t_1)+F(t_2)-2F(s)`,

and cutoff equations remain

`s q_T = kappa_T+p_T`,
`t_i(q_i-q_T)=kappa_L-kappa_T-p_T`.

Together with partner fixed points these form a differentiable implicit system. If its Jacobian, the private-price SOC, and public SOCs are nonsingular, equilibrium derivatives vary continuously with `(F,f)`.

Because P2-R4 is a pair of strict inequalities at the uniform witness, it survives sufficiently small `C^1` perturbations of the uniform distribution. Verdict: **local structural generality established; arbitrary-distribution global generality not established.**
