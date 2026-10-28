# Numerical verification of the one-parameter extension

Source baseline: PR #59 head `d2520c9afbf9b2b742a5f04fbba6159eb489db61`.
All original numerical primitives are imported from the unchanged Stage11R exact-rational certificate.
Only the taste-skill copula is changed.
The pooled-half scenario remains an original robustness comparator, not an additional technology assumption for the extension.

## Reproduction

Use Python 3.11 or newer and `python -m pip install -r exploration_kappa_20261028/requirements.txt`.
The run used Python 3.12.14 with NumPy 2.3.5 and SciPy 1.17.0, identical to the repository dependency pins.
Then run:

```bash
python exploration_kappa_20261028/numerics.py
python exploration_kappa_20261028/exact_threshold.py
python exploration_kappa_20261028/numerics.py --seed-checks
```

The exact conditional-threshold certificate uses only the standard library and the frozen primitive definitions.
It also runs under a Python environment without NumPy or SciPy.
The numerical solver requires the two pinned packages.
`results.json`, `exact_threshold.json`, `numerics.log`, `price_seed_checks.json`, and `frozen_regression.log` record the results.

## Quantities and welfare accounting

`moments.py` computes the exact integrals of `(1, theta, z)` over the participation half-plane.
The independent component is integrated as a polynomial on each interval between the two square-edge crossings; the diagonal component uses its clipped cutoff.
Fraction inputs therefore return exact rational values.
The quantities are unnormalised population shares: `m=E[1_join]`, `q=E[theta 1_join]`, and `r=E[z 1_join]`.
Private-only direct service welfare uses `N*(u*m+Delta*r)`.
Co-provision has full total participation, `m=1,q=1/2`, but the premium upgrade moment is `E[z 1{z>=1/2}]=3/8`, independent of kappa.
Fee receipts are local transfers, not an extra welfare subtraction.

Exact rational checks cover empty/full participation, diagonal reduction `q=r=m-m^2/2`, both kappa endpoints, uniform marginal correlation, a square triangle with `m=1/2,q=r=1/3`, and residual ability selection at kappa zero.

## Price and participation calculations

For each P/P deviation, the other fee is held fixed and both regions' participation is recomputed simultaneously.
The program never substitutes a fixed-other-participation monopoly best response for this price-game best response.
Best responses search the entire feasible fee interval `[0,u+a+Delta+b/2]`, evaluate all grid cells and refine every detected local maximum.
In the original parameter sets the final candidate is sharpened using the exact central-strip derivatives and the implicit four-dimensional `(m_A,q_A,m_B,q_B)` demand Jacobian.
The final candidate is checked again with a denser 400-cell deviation grid and local refinements.
Demand seeds zero and one agree to floating-point precision; three price seeds `(0,0)`, `(.65,.65)`, `(1.3,1.3)` agree within `1.2e-16` on the nine endpoint/midpoint checks recorded separately.

The participation stage is separately covered by an all-fee contraction argument: each of `m,q` is separately `1/Delta`-Lipschitz in `U` and `h`, so the joint sup-norm participation map has Lipschitz bound `max_i(a_i+b_i)/Delta_i`. Exact original-scenario bounds, all strictly below one, are recorded in `participation_contraction.json`. Banach therefore gives a unique participation fixed point for every fee pair and every kappa in `[0,1]`. This theorem covers participation only. The price checks remain numerical diagnostics, not a general global-best-response, pure-price existence, or uniqueness theorem for kappa below one.
Undetected narrow optima or equilibria are not excluded mathematically by a finite grid.
Consequently the four-player entry results below are numerical equilibrium candidates whenever they rely on the off-path P/P continuation.

## Regression at kappa one

The four Stage11R checks and the original Stage33 reproduction all exit zero; see `frozen_regression.log`.
For every one of the 16 configurations in each of the three original scenarios, the extended solver's fees, participation masses, skill moments, local welfare, and private profits agree with the frozen exact-rational enclosure midpoints to maximum absolute error `8.9e-16`.
The discrepancy is floating-point rounding; exact equality of the model formulas at kappa one is also checked by Fraction identities.
The frozen exact P/P and entry certificates remain unchanged and retain their original scope.

## Sweep with original primitives

The grid is exactly `kappa=0,0.1,...,1`; there is no change in costs, technology, product markets, ownership shares or other numerical parameters.

| Original scenario | Four-player NE candidates | Restricted planner on the grid | Conditional B-government threshold at kappa 0 | Threshold at kappa 1 |
|---|---|---|---:|---:|
| Skill, lambda=0 | BP at every grid point | BB throughout | 0.0404250710071206 | 0.0404250710071206 |
| Skill, lambda=1/2 | BP at 0 through 0.9; BB at 1 | BB throughout | 0.0400434337024524 | 0.0405316701024489 |
| Pooled, lambda=1/2 | BP at every grid point | BB throughout | 0.0400432292331624 | 0.0400432292331624 |

The fixed government cost is `F_GB=0.0405`.
The skill-half case therefore provides a reversal under the original cost, without choosing a new cost to manufacture one.
Lower correlation reduces the value of adding public basic access in this example.
Its policy scale remains limited: the entire kappa-induced threshold range is about `0.00048824` in the original welfare units.
A numeric reversal in one calibrated scenario is not a general new economic theorem.

## Exact conditional threshold and its limits

`exact_threshold.py` strengthens the B-government comparison between BP and BB in the skill-half scenario.
With opposite A already full, write `U=u_B+a_B`, `h=b_B/2`, and `Delta=1`.
In the central price strip `h<f-U<1`, demand is affine in the fee with
`A=kappa/(1+h)+(1-kappa)` and `c=(1-kappa)*h/2`.
Thus the conditional fee is `(1+c+A*U)/(2*A)`.
For all kappa, prices below the strip have revenue at most `U+h`, and prices above the strip have revenue at most `(U+h+1)*h/(1+h)`.
The feasible trial fee `(1+U)/2` earns at least `(1+U)^2/4`, which strictly exceeds both bounds in exact rational arithmetic for the original parameters.
The central-strip candidate is interior throughout kappa in `[0,1]` and has strictly concave quadratic revenue there.
This establishes the conditional global optimum over all nonnegative fees, including corners and kinks.

Substitution yields a public-entry threshold that is a cubic-over-quadratic rational function of kappa. The numerator is constructed directly by exact polynomial operations, without interpolation: write the fee excess as `x=X/B` with affine `X,B`; each diagonal and independent moment has a numerator of degree at most two over `B^2`, and multiplying by `kappa` or `1-kappa` gives degree at most three. The mass numerator is affine, so the remaining mass terms multiplied by `B^2` are at most cubic as well.
Its derivative has a positive denominator and a numerator whose Bernstein coefficients on `[0,1]` are all strictly positive.
All coefficients and signs are checked with Fraction arithmetic and recorded in `exact_threshold.json`.
The threshold is therefore strictly increasing, with exactly one crossing of `F_GB=0.0405`:
`kappa*=0.9351610109823004` (an exact rational isolating interval is in the JSON).

This proves a conditional public-entry switch with the other region full.
It does not alone prove a full subgame-perfect equilibrium switch across the continuum, because all off-path P/P price continuations and the other entry deviations must also be certified uniformly.
The grid results provide corroboration, with the explicit numerical limitations above.
