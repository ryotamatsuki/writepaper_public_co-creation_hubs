# Equilibrium-Manifold Reduction

Raw differentiation is not the relevant complexity metric. Every sign object was
reduced after imposing equilibrium identities.

## B3

At the canonical rational primitives, after eliminating `x` and reducing modulo
the symmetric participation quadratic:

- B3 public FOC numerator: linear in `q_T`, degree 8 in `t`, degree 4 in `p`;
- B3 public cross numerator: linear in `q_T`, degree 18 in `t`, degree 8 in `p`;
- B3 public own/SOC numerator: linear in `q_T`, degree 18 in `t`, degree 9 in `p`.

Eliminating `q_T` from the B3 FOC gives one non-boundary core polynomial of
degree 14 in `t` and degree 7 in `p`.

The own and cross derivatives share the same reduced denominator structure at
symmetry. This is useful because the SOC and cross-sign tests can be organized
as polynomial sign tests once denominator signs are included in the regular
branch definition.

## G

The raw exact G public FOC numerator at the canonical rational primitives has
roughly fourteen thousand operations before manifold reduction. After:

1. reconstructing `x`,
2. reducing modulo the participation quadratic,
3. eliminating the reduced private FOC,
4. reducing modulo the quartic private equilibrium core,

the relevant public-equilibrium relation is cubic in `p` (although high degree
in `t`). The subsequent resultant factors strongly.

This confirms Frontier O: large off-manifold formulas can have much smaller
normal forms on the equilibrium ideal.
