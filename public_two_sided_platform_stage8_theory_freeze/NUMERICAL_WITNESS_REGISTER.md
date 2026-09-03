# Numerical Witness Register

Primitives:
`v=.1`, `alpha=.5`, `beta=.05`, `rho=.15`, `rho_T=.05`, `kappa_L=.27`, `kappa_T=.02`, `tau=.05`, `gamma=.9`.

## G
`x_1=x_2=.684028`, `p_T=.0227467`, own reduced second derivative `-.199391`, reduced cross derivative `-.00480885`, `BR'_G=-.024118`.

## B3
Fix `bar p=.0227467`. Stage-4 strategic solver reports `x_B3=.656020`, own second derivative `-.118358`, cross derivative `+.0135412`, `BR'_B3=.114409`.

Stage-7 welfare recomputation reports `x_B3=.656016`. This is recorded as a numerical precision/solver-display reconciliation, not a second theoretical equilibrium. For strategic proposition reporting use Stage-4 `.656020`; for reproducing the Stage-7 welfare table use `.656016`.

At G, `p^*_{x_i}=-.016505`.

No number in this register is a calibration to observed data; it is a deterministic existence witness.