# Stage L3-1 — Multi-Frontier Analytic Solvability Reopen

This directory is an experimental stacked theory-upgrade track based on PR #39 HEAD
`c4d7936a98d685a6f1ddc0a24492ab50f129deae`. It does not modify the production
manuscript, Stage-8 frozen model, or Stage-12 positioning.

## Reopen result in one sentence

L3-0's statement that a primitive necessary-and-sufficient characterization may
"not exist" is too strong. After denominator clearing the central regime is
semialgebraic, so its primitive projection exists by Tarski–Seidenberg. The actual
question is computational and expository: whether that exact projection can be
computed and compressed into a scientifically usable finite condition set.

## Strongest new structural finding

On the symmetric regular-interior branch, use `(t,q_T,p)` instead of `(x,t,q_T,p)`.
The exact identities are

`q_T^2-(A_T+2*delta*t)q_T+2*delta*a=0`,

`x=(q_T+d/t-delta(1-t)-v)/alpha-rho`,

and

`delta=q_T(q_T-A_T)/(2(t q_T-a))`.

The symmetric private FOC, after reduction modulo the participation quadratic,
is linear in `q_T`; eliminating `q_T` yields boundary factors times one core
polynomial that is degree 7 in `t` and quartic in `p`. This is materially smaller
than the raw degree-8/9 systems reported in L3-0.

## Current status

The canonical full-game equilibrium can be projected to exact univariate algebraic
components; its economically relevant `t` and `p` roots are Sturm-isolatable.
However, the joint matched-price elimination linking that algebraic G price to
the distinct B3 public equilibrium remains the main bottleneck, and a full
primitive necessary-and-sufficient reversal region was not produced in this stage.

See `LEVEL3_REOPEN_VERDICT.md` for the final classification.
