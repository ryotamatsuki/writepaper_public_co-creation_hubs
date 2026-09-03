# Stage 4R-TP — Executive Verdict

## Canonical verdict

**GO**

Primary result: in the smooth interior regime `0 -> H_T -> H_i`, there is a nonempty open parameter region in which fixed-price B3 has **strategic complements** between regional public investments while full game G with endogenous private pricing has **strategic substitutes**.

A deterministic witness is
`v=.1, alpha=.5, beta=.05, rho=.15, rhoT=.05, kappa_L=.27, kappa_T=.02, tau=.05, gamma=.9`.

At the symmetric G equilibrium: `x*=0.684028`, `p_T*=0.0227467`, own second derivative `-0.199391`, cross derivative `-0.00480885`, hence `BR'_G=-0.024118`.

Fixing `p_T=0.0227467` in B3 gives the symmetric B3 equilibrium `x_B3=0.656020`, own second derivative `-0.118358`, cross derivative `+0.0135412`, hence `BR'_B3=+0.114409`.

Thus endogenous private pricing reverses the strategic relation from complements to substitutes. Twenty deterministic ±0.5% local perturbations (seed `20260903`) preserved the reversal.

P1-R is killed as a headline: in the sign-reversal region the national planner pushes investment toward a regime boundary, so no clean interior novel decentralization threshold was obtained.

Route: **GO -> Stage 6R-TP Proposition-Level Novelty Re-Kill**. Stage 5 is skipped.