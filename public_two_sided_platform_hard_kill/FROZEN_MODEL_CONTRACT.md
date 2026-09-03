# Frozen Model Contract

Players: `G_1,G_2,H_1,H_2,H_T`.

Each region has unit mass of projects `z~U[0,1]`; each project single-homes among `0,H1,H2,HT`.

`U_rh(z)=z[v+alpha n_h^P]-kappa_rh-1{h=T}p_T`, `U_0=0`.

`kappa_rr=kappa_L`, `kappa_rj=kappa_L+tau`. The metro primitive is taken symmetrically as `kappa_1T=kappa_2T=kappa_T`; it remains a real access/coordination cost and is not absorbed into the fee because private revenue is `p_T n_T^F`.

Partners have per-platform cost `c~U[0,1]`, may multihome, and satisfy interior equations `n_i^P=rho+x_i+beta n_i^F`, `n_T^P=rhoT+beta n_T^F`.

`x_i in [0,1]` is partner-side facilitation; cost `gamma x_i^2/2`. Private platform chooses one project-side fee `p_T>=0`.

Timing: `(x1,x2) -> p_T -> participation fixed point -> surplus`. Equilibrium: SPNE.