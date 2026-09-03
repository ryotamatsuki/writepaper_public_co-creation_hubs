# Method Matrix

| Frontier | Mathematical idea | Exact? | Primitive projection? | N&S potential | Result / blocker |
|---|---|---:|---:|---:|---|
| A | Reparameterize by `(t,q_T,p)` | yes | partial | high | **PASS**; removes `x` exactly and lowers degrees |
| B | Symmetry-first reduction | yes | partial | high | **PASS**; first serious branch |
| C | Resultants | yes | partial | high | **PASS/PARTIAL**; private core quartic in `p`, canonical G factorization; direct B3-G coupling resultant times out, but the canonical coupling obligation is closed by the rational Krawczyk route |
| D | Gröbner / elimination ideals | yes | partial | high | **PASS at canonical G component**; degree-31 component gives a linear `p(t)` relation |
| E | Discriminant / root geometry | yes | partial | medium | **PASS at witness**; exact rational Sturm intervals for G `t,p` |
| F | Quantifier elimination / CAD | in principle | yes | high | tool unavailable; full parameter CAD not computed |
| G | Sign at algebraic root | yes | partial | high | **PASS at canonical coupled root**; rational interval implicit differentiation closes G/B3 SOC and cross signs; uniform primitive sign projection remains open |
| H | Positivity certificates | potentially | partial | medium | no SOS/CAD stack; mixed-sign reduced polynomials give no immediate global primitive certificate |
| I | Monotonicity | exact target | no result | medium | small-`beta` local result survives; no global sign-stable derivative found |
| J | Single crossing | exact target | no result | medium | no certified single-crossing theorem |
| K | Homogeneity / dimensionless indices | exact | limited | medium | no legitimate normalization removing several primitives without changing economics |
| L | Higher-order `beta` expansion | exact | local | medium | first-order L3-0 theorem retained; not enough for N&S |
| M | Solve for network parameter | yes | yes | high | **PASS**: `delta=q(q-A_T)/(2(tq-a))` |
| N | Solve zero-crossing system | yes | partial | high | B3 boundary polynomial constructed; full primitive G boundary elimination remains large |
| O | Reduce on equilibrium ideal | yes | partial | high | **STRONG PASS**; very large G FOC drops to low `p` degree after sequential reduction |
| P | Envelope simplification | yes | partial | medium | augmented follower system is cleaner; four G policy terms do not vanish generically |
| Q | Schur complement / block Jacobian | yes | implicit | high | exact augmented-Jacobian representation available |
| R | Determinant / orientation | yes | implicit | high | participation responses governed by `u` and `H`; useful denominator control |
| S | Root-location transforms / Sturm | yes | partial | medium | Sturm isolation works directly; interval transform not needed for witness |
| T | Computer-assisted exact proof | yes | partial | high | **PASS at canonical calibration**; rational Krawczyk plus exact interval arithmetic close the coupled root/sign certificate; advanced CAD/QE tools unavailable for full primitive projection |

A method timeout is recorded as a method-specific failure only, not as evidence of
mathematical impossibility. The remaining Level-3 gap is parametric projection
and branch selection, not the canonical G/B3 certificate.
