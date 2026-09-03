# Proof Obligation Register

| ID | Statement | Required for L3? | Analytic status | Assumptions | Dependence | Method | Numerical confirmation | Unresolved issue | Blocking? |
|---|---|---:|---|---|---|---|---|---|---:|
| PO1 | central participation reduction | yes | PASS | central sorting | primitive + cutoffs | exact quadratic/reconstruction | yes | root admissibility by inequalities | no |
| PO2 | participation regularity | yes | CONDITIONAL | `det J != 0` | equilibrium | Jacobian/IFT | strict at witness | not primitive-global | yes for L3A/B |
| PO3 | private SOC | yes | CONDITIONAL | `R_p<0` | equilibrium | exact FOC/SOC | strict at witness | not primitive-global | yes for L3A/B |
| PO4 | private `p_i` derivative | yes | PASS FORMULA | private SOC | reduced derivatives | IFT | negative at witness | sign not primitive-global | yes for L3A/B |
| PO5 | private `p_ij` derivative | yes | PASS FORMULA | private SOC | reduced derivatives | second IFT | diagnostic reconciliation | sign not primitive-global | no for block theorem |
| PO6 | B3 public SOC | yes | CONDITIONAL | explicit SOC subregion | equilibrium | Hessian | strict at witness | nonconcave counterexamples exist | yes for primitive L3 |
| PO7 | G public SOC | yes | CONDITIONAL | explicit SOC subregion | equilibrium | reduced Hessian | strict at witness | nonconcave counterexamples exist | yes for primitive L3 |
| PO8 | B3 cross derivative | yes | PASS BLOCK | regular B3 | equilibrium derivative | implicit differentiation | positive at witness | primitive sign not global | yes for L3A/B |
| PO9 | B3 beta-zero derivative | no / L2 | PROVED | central `d,Delta_i>0` | primitives/state at beta 0 | first-order expansion | symbolic identity | epsilon not explicit | no for L2 |
| PO10 | G full cross decomposition | yes | PASS | smooth price policy | derivative blocks | exact chain rule | matches production sign | four price terms required | no |
| PO11 | G beta-zero cross sign | no / L2 | PROVED | symmetric regular beta-zero state | primitive state index | closed-form substitution | symbolic identity | branch existence separate | no for L2 |
| PO12 | headline N&S equivalence | yes | PASS LEVEL-3C | two regular SOC branches | two equilibria | BR/Hessian identity | canonical witness | not primitive | yes for full L3 |
| PO13 | primitive reduction | yes | FAIL HARD KILL | frozen model | primitives + two roots | elimination/factor search | n/a | expression/root explosion | fatal for L3 |
| PO14 | necessity in primitives | yes | NOT OBTAINED | — | — | — | — | only derivative-block necessity | fatal for L3 |
| PO15 | sufficient analytic theorem | fallback | PASS LEVEL 2 | regular beta-zero branches | local beta | exact beta-zero signs + continuity | consistent with mechanism | implicit epsilon | no |
| PO16 | nonempty primitive L3 region | yes | NOT PROVED | — | — | — | witness remains numerical | no explicit primitive region | fatal for L3 |
| PO17 | canonical witness inclusion in L3 region | yes | NOT APPLICABLE | no primitive L3 region | — | — | strict witness verified | no explicit threshold | fatal for L3 GO |

## Register conclusion
The model supports exact analytic reduction, exact response identities, and a nontrivial small-beta sufficient theorem. It does not support the requested primitive necessary-and-sufficient Level-3 characterization under the hard-kill tractability standard.