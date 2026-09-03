# Expression Complexity Register

Operation counts are from exact SymPy experiments and are method-dependent.

| Object | Raw / pre-reduction | Equilibrium-ideal reduced / final route | Degree / structure | Exact status |
|---|---:|---:|---|---|
| symmetric participation | small | quadratic | degree 2 in `q_T` | branch-dependent |
| generic symmetric private FOC numerator | ~182 ops | remainder ~486 ops | linear `q_T`, degree 5 `t`, degree 3 `p` | exact block |
| generic private resultant core | n/a | ~1,335 ops | degree 7 `t`, degree 4 `p` | branch-dependent |
| canonical G public FOC numerator | ~14,339 ops | q-remainder ~1,975 ops; private-ideal remainder ~420 ops in one ordering | cubic `p`, high `t` degree | equilibrium equation exact |
| canonical G `t` projection | degree 105 before factor selection | factors incl. degree 14 and 31 | witness on degree 31 | root exactly isolatable |
| canonical component Gröbner basis | large coefficients | 2 basis elements | degree-31 `t`; linear `p` | exact algebraic branch |
| B3 FOC numerator | hundreds | ~239 ops | linear `q_T`, degree 8 `t`, degree 4 `p` | exact block |
| B3 cross numerator | thousands | ~506 ops | linear `q_T`, degree 18 `t`, degree 8 `p` | exact block |
| B3 own/SOC numerator | thousands | ~838 ops | linear `q_T`, degree 18 `t`, degree 9 `p` | exact block |
| complete G cross expression | ~390k raw ops | direct common-denominator form >425k ops; degree extraction timed out | direct full expansion rejected | canonical sign closed exactly by reduced rational interval implicit differentiation |
| G-price/B3 matched coupling | high | direct elimination >60s; replaced by five-polynomial rational Krawczyk system | variables `(t_G,q_G,p_G,t_B,q_B)` with one common `p_G` | exact canonical coupling certified |

The central lesson is that raw operation count is a poor kill criterion. Several
large objects collapse dramatically after equilibrium-manifold reduction, and
direct expansion/resultant timeouts are method-specific failures rather than
unresolved proof obligations. L3-1 closes the canonical G cross-sign and
matched-price coupling through alternative exact certificates. The remaining
Level-3 blocker is the non-calibrated parametric semialgebraic projection and
uniform branch selection across primitive space.
