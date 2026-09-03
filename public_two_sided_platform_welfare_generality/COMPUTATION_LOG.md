# Computation Log

Baseline parameters: `v=.1, alpha=.5, beta=.05, rho=.15, rhoT=.05, kappa_L=.27, kappa_T=.02, tau=.05, gamma=.9`.

Frozen Stage-4 canonical certificate: G `x*=.684028`, `p*=.0227467`, `BR'_G=-.024118`; B3 `x*=.656020`, `BR'_B3=.114409`.

Stage-7 deterministic re-evaluation uses the same smooth `0->H_T->H_i` equations, SciPy root/minimize routines, centered finite differences, and the same fee benchmark. Small numerical differences in recalculated BR slopes are treated as solver/step-size effects; Stage 4 values remain canonical.

Stage-7 new values:
- G `W_i≈.223019`, `Pi≈.012603`, `W^N≈.458642`;
- B3 `W_i≈.210910`, `Pi≈.013891`, `W^N≈.435711`;
- at G `dW_j/dx_i≈.44471`, `dPi/dx_i≈-.02194`, `dW^N/dx_i≈.42277`;
- local `p^*_{x_i}≈-.0165`.

All scripts are deterministic; no random search is used for new Stage-7 claims.
