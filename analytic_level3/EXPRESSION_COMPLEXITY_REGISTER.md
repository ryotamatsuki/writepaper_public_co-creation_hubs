# Expression Complexity Register

Counts below come from bounded exploratory SymPy derivations on the frozen central-regime equations. They are diagnostic rather than mathematical invariants; simplifier choices can change operation counts. The purpose is to identify where primitive elimination ceases to be an intelligible proof strategy.

| Object | Representation | Diagnostic complexity | Degree / terms | Sign-normalized interpretation | Tractability |
|---|---|---:|---|---|---|
| participation `q_T` | explicit quadratic | low | degree 2 | select admissible root by regime inequalities | GOOD |
| reconstructed `x_i` | rational in `(t_i,q_T,p)` | low | simple rational | exact | GOOD |
| symmetric participation Jacobian | `det J=uH` | low | factored | regular iff `uH != 0` | GOOD |
| private FOC after symmetric participation reconstruction | rational | ~6k raw ops in exploratory form | numerator degree 8 in `(t,q_T,p)`, 36 monomials | private optimum condition | MODERATE-HIGH |
| B3 symmetric public FOC numerator | rational/polynomial | ~3.3k numerator ops | degree 9, 58 monomials | public stationary condition | HIGH |
| B3 general cross derivative | implicit rational block | ~26k raw ops | multi-parameter | fixed-price network force | HIGH |
| private `p_x` after full implicit reduction | rational derivative block | ~115k raw ops before full factorization | high-dimensional | repricing response | VERY HIGH |
| G public FOC with exact `p_x` substituted | rational derivative block | ~124k raw ops | high-dimensional | endogenous-price public FOC | VERY HIGH |
| G cross derivative | exact block formula retained | not expanded | requires `p_i,p_j,p_ij` | fixed term + four price-policy terms | EXPANSION REJECTED |

## Factorization result
The low-dimensional participation identities factor cleanly and are retained. Direct aggressive `cancel/factor` attempts on the second-derivative B3/full-G objects exceeded the bounded simplification window before producing an economically useful primitive factorization.

## Hard-kill implication
An enormous resultant or polynomial sign condition could in principle be called “analytic,” but it would fail the Stage L3-0 success standard because it would not characterize an economically interpretable primitive region. The proof route therefore stops before symbolic dumping.