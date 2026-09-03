# Quantifier-Elimination Audit

## Formal semialgebraic status

After restricting to the central active regime and recording the sign of every
cleared denominator, the frozen model consists of finitely many polynomial
equalities and strict/non-strict polynomial inequalities in primitive parameters
and endogenous/algebraic variables. Adding regularity, SOC, matched-price
coupling and `M_B3>0>M_G` preserves semialgebraicity.

Hence the full reversal set in `(theta,z)` is semialgebraic, and its projection
onto primitive `theta` is semialgebraic by Tarski–Seidenberg.

Therefore L3-0's suggestion that a primitive necessary-and-sufficient
characterization may not exist is too strong. It exists in principle. The real
question is whether the projection can be computed and compressed into a
scientifically usable form.

## Actually computed

- generic symmetric private-block elimination;
- canonical full-game elimination to algebraic components;
- B3 equilibrium-manifold polynomialization;
- exact canonical G root isolation.

## Not computed

A full nine-primitive CAD / quantifier-free projection was not produced.

Environment: SymPy 1.14, SciPy and mpmath were available. SageMath, Singular,
Mathematica/Wolfram Engine, Maple, QEPCAD, Redlog/Reduce, Z3 and dReal were
unavailable. SymPy did not provide a practical full CAD route for this problem.

## Verdict

**Theoretical PASS; computational INCOMPLETE.** Tarski–Seidenberg settles
existence, not scientific usability.
