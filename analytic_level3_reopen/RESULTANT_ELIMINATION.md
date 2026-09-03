# Resultant Elimination

## Generic private block

The symmetric participation quadratic and the reduced private FOC were eliminated
exactly with SymPy resultants.

After reducing the private FOC numerator modulo the participation quadratic, the
result is linear in `q_T`, degree 5 in `t`, and degree 3 in `p`.

The generic resultant factors as

- `alpha`,
- `beta`,
- `(kappa_T+p)^2`,
- one non-boundary core polynomial of degree 7 in `t` and degree 4 in `p`.

The `(kappa_T+p)` factor is the `a=0` boundary and is extraneous to the target
regular-interior positive-access-price branch.

## Canonical full game

With the canonical primitive vector represented exactly as rationals, the
full-game public FOC was constructed through the augmented follower system
(participation equations plus the private FOC), then reduced sequentially:

1. eliminate `x` by exact reconstruction;
2. reduce modulo the participation quadratic;
3. eliminate `q_T` with the reduced private FOC;
4. reduce the public FOC modulo the private equilibrium core in `p`;
5. eliminate `p`.

The final univariate `t` resultant is degree 105 before factor selection but
splits heavily. Its exact factor profile contains low-degree boundary/extraneous
pieces and irreducible degree-14 and degree-31 components. The canonical G
witness lies on the degree-31 component.

This is an important correction to L3-0: raw operation count is not preserved
under equilibrium-manifold reduction.

## Matched B3 coupling

For fixed `p`, the B3 public FOC reduces modulo participation to an expression
linear in `q_T`, with degrees 8 in `t` and 4 in `p`. Eliminating `q_T` gives a
boundary factor `(50p+1)^2` under the canonical rationalization and one core
polynomial of degree 14 in `t` and 7 in `p`.

The next exact coupling step—eliminating the algebraic G price between the
degree-31 G price component and the B3 core—was attempted directly. SymPy did
not finish within the 60-second resource limit. This is recorded as a
method-specific computational obstruction, not mathematical nonexistence.
