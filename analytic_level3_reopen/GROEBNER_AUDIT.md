# Gröbner Audit

The highest-value Gröbner calculation was performed after factor selection, not
on the raw full system.

For the canonical G witness component, combine:

- the degree-31 univariate `t` factor,
- the symmetric private equilibrium core,
- the equilibrium-ideal-reduced public FOC.

Using lexicographic order with `p` before `t`, SymPy returns a two-polynomial
basis:

1. the degree-31 polynomial in `t`;
2. a polynomial **linear in `p`** and degree 30 in `t`.

Therefore, on this irreducible component,

`p = P_num(t)/P_den`

as an exact rational function in the algebraic root `t`.

This is a genuine algebraic-function characterization and demonstrates that the
canonical full-game equilibrium does not require radicals.

The raw multivariate basis was not used as a proof object because the reduced
component basis is dramatically smaller and economically cleaner.

Tool note: SageMath, Singular, Mathematica/Wolfram Engine, Maple, QEPCAD,
Redlog/Reduce, Z3 and dReal were unavailable in the execution environment.
