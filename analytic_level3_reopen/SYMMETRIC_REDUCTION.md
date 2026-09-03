# Symmetric Reduction

The first serious Level-3 target is the symmetric regular-interior branch, not a
fully asymmetric global theorem.

## Participation

With `x_1=x_2=x` and `t_1=t_2=t`, the exact participation system reduces to one
quadratic in `q_T` plus the reconstruction formula for `x`. No radical expression
is required; admissible roots can be selected by interval/Sturm conditions.

## Private price FOC

Using the symmetric Jacobian, the exact responses to `p` are rational:

`q_p=-2 delta(q_T+u)/(q_T H)`,

`t_p=(t q_p-1)/u`.

For private profit `Pi=2p(t-a/q_T)`, the FOC is

`R=2(t-a/q_T)+2p(t_p-s_p)=0`.

After reducing the numerator of `R` modulo the participation quadratic:

- degree in `q_T`: 1;
- degree in `t`: 5;
- degree in `p`: 3.

Eliminating `q_T` yields factors `alpha`, `beta`, `(kappa_T+p)^2`, and a single
non-boundary core polynomial with degrees:

- degree 7 in `t`;
- degree 4 in `p`.

This is an exact generic result, not a numerical calibration artifact.

## Implication

The symmetric private follower is algebraic of materially lower degree than the
raw L3-0 representation suggested. This reopens L3-R-style algebraic-root
characterization, but does not by itself solve the public FOC, matched B3
equilibrium, SOCs, or cross-sign projection.
