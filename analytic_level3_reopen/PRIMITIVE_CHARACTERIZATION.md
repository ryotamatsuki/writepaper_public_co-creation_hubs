# Primitive Characterization

## Exact quantified formulation

Let `z` collect the symmetric G equilibrium variables, the matched B3 variables,
and any derivative-response auxiliaries. After denominator clearing, define
`Phi(theta,z)` as the conjunction of:

1. central regular-interior participation equations and inequalities;
2. private FOC and the declared private optimality conditions;
3. G public FOC and strict public SOC;
4. `p_B3=p_G`;
5. B3 public FOC and strict public SOC;
6. `M_B3>0`;
7. `M_G<0`;
8. nonsingularity / branch-isolation inequalities.

Then

`theta in Theta_REV  <=>  exists z : Phi(theta,z)`.

This is an exact quantified semialgebraic characterization. By
Tarski–Seidenberg it has an equivalent quantifier-free primitive form.

## Why this does not yet count as L3-S

The stage success standard requires the projected conditions themselves, or an
algebraic-root characterization with completed unique branch isolation and sign
queries. L3-1 has not computed the full projection and has not completed the G
cross/SOC sign query together with matched B3 coupling.

Therefore this file records mathematical existence, not a claimed Level-3
success.

## Status

- necessary conditions: identified;
- sufficient conditions: existing small-`beta` Level-2 theorem survives;
- primitive necessary-and-sufficient formula: **exists in principle, not
  explicitly computed**;
- primitive-only explicit condition set: **NO**;
- unique threshold: **NO**.
