# Numerical Hard Kill

Deterministic witness primitives:
`v=.1, alpha=.5, beta=.05, rho=.15, rhoT=.05, kappa_L=.27, kappa_T=.02, tau=.05, gamma=.9`.

G: `x*=.684028`, `p*=.0227467`, `BR'_G=-.024118`.
B3 with `bar p=p*_G`: `x*=.656020`, `BR'_B3=+.114409`.

All active masses and SOCs are strict.

Local robustness experiment: seed `20260903`, 20 independent draws perturbing the positive primitives by ±0.5% around the witness. **20/20** retained an interior G equilibrium, an interior B3 equilibrium, `BR'_G<0`, and `BR'_B3>0`.

The numerical audit supports the analytic IFT/continuity open-set argument; it is not used as a universal proof over the parameter space.