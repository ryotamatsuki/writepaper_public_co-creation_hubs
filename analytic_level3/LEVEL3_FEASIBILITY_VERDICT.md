# Level-3 Feasibility Hard-Kill Verdict

## Baseline

- repository: `ryotamatsuki/writepaper_public_co-creation_hubs`
- immutable L3-0 base: `main@6132a4f063012aa34f473b3a9c0185b546486b38`
- Stage-12 PR at startup: `#38`, open/unmerged, head `c8ad54653ae19ae0fb11763ea0021779574d771b`, CI success
- production manuscript/theory changes in this track: **NONE**

## Gate results

1. **Regular interior representation — PASS.** The participation block reduces exactly to a quadratic in `q_T` plus explicit reconstruction of `x_i`; admissible-root selection is encoded by strict economic inequalities.
2. **B3/G analytic BR reduction — PASS IN BLOCK FORM.** IFT gives exact BR-slope formulas on regular SOC branches.
3. **SOC/denominator control — CONDITIONAL PASS.** Public SOCs must be included explicitly in `Theta^{RI,SOC}`; they are not primitive-global consequences.
4. **Economic full-game decomposition — PASS.** The G cross derivative equals a fixed-price/direct term plus four exact private-price-policy terms involving `p_i,p_j,p_i p_j,p_{ij}`.
5. **Necessary-and-sufficient reversal condition — PASS ONLY AT LEVEL-3C.** It is exact conditional on the two selected equilibrium branches and derivative blocks.
6. **Primitive reduction — FAIL HARD KILL.** The headline compares distinct B3 and G public equilibria; eliminating both roots, the matched G price, SOC restrictions, and first/second private-policy derivatives produces high-degree/high-operation rational systems without a transparent finite primitive factorization.
7. **Primitive nonempty reversal region — NOT PROVED.** The canonical numerical witness/open set remains valid; no primitive Level-3 region containing it was derived.
8. **Theory-value test — PASS FOR LEVEL 2.** A nontrivial analytic small-`beta` sufficient theorem is available.

## Strong Level-2 theorem obtained

Consider symmetric regular interior B3 and G equilibrium branches that exist at `beta=0`, satisfy strict private/public SOCs and nonsingularity, and continue smoothly for small `beta`.

At `beta=0`, B3 has zero public-public cross effect. The exact first derivative is

`d M_B3/d beta |0 = alpha^3 d^3/(Delta_i^3 Delta_j^2)>0`.

For the full game at a symmetric regular beta-zero state,

`M_G(0)=-3 T alpha^2 kappa_L^2/[16 Delta(Delta+T)^3]<0`.

Therefore by continuity there exists `epsilon>0` such that for every `beta in (0,epsilon)` on those continuing SOC branches,

`BR_i^{B3′}>0>BR_i^{G′}`.

This is analytic, economically interpretable and strictly stronger than a bare numerical mechanism decomposition. It is not Level 3 because `epsilon` and the beta-zero branch/SOC existence conditions have not been reduced to a necessary-and-sufficient primitive region.

## Why Level 3 is killed rather than merely unfinished

The hard-kill criterion is not whether a computer algebra system could eventually return a resultant. It is whether the frozen model admits an economically interpretable necessary-and-sufficient primitive characterization. Already before the final G cross expansion, the symmetric private and public stationary systems are degree-8/degree-9 rational polynomial systems, and the exact private-policy/G-FOC expressions expand into six-digit operation counts. Direct elimination would create a root-ordering/resultant condition rather than a usable economic threshold.

The model must not be simplified or altered solely to manufacture a Level-3A/B theorem.

# Final verdict

**NO-GO LEVEL 3 / GO LEVEL 2**

Level 3, as defined by this task, is not analytically tractable in an economically useful primitive necessary-and-sufficient form under the frozen model. A strong analytic Level-2 small-beta theorem and exact derivative decomposition are feasible and should be preserved for a later explicit upgrade decision.

## Next authorized research step

A separate Level-2 formal-proof stage may turn the small-beta result into publication-quality lemmas/theorem, audit all beta-zero regularity assumptions, and assess whether it materially improves the paper enough to reopen Stage 12 journal positioning. The production manuscript remains frozen until that later decision.