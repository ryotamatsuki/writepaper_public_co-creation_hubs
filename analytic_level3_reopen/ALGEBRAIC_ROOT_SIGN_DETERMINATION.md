# Algebraic-Root Sign Determination

## What succeeded

The canonical G equilibrium root is exactly isolatable without radicals.
Sturm root counts identify unique rational intervals for both `t_G` and `p_G`.
A component Gröbner basis then makes `p_G` a rational function of the isolated
degree-31 algebraic `t_G`.

This is the correct architecture for exact sign determination:

1. isolate the target root by a Thom/Sturm encoding;
2. reduce the sign polynomial modulo the equilibrium ideal;
3. evaluate its sign at that encoded root through a Sturm/subresultant query.

## B3 cross and SOC blocks

For fixed matched price, exact B3 cross and own derivative numerators were
reduced to polynomials linear in `q_T`. Thus their zero boundaries can be
eliminated jointly with the participation and B3 FOC equations. This is a
tractable exact sign-query architecture in principle.

## Remaining G obstruction

The exact G reduced Hessian can be written from the augmented follower system,
but direct symbolic construction of the complete cross derivative produces an
expression of roughly 390k raw operations before common-denominator processing.
A `together` form exceeded 425k operations. Attempting full degree extraction /
expansion hit the imposed 60-second limit.

Therefore this stage did **not** complete the final Sturm/subresultant sign query
for `M_G` at the algebraic root. The obstacle is computational expression
management, not inability to define the sign exactly.

No claim of an exact G cross-sign certificate is made here.
