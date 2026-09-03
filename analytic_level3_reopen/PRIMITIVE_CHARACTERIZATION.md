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
Tarski-Seidenberg it has an equivalent quantifier-free primitive form.

## Why this does not yet count as Level 3

At the canonical rational calibration, L3-1 does complete unique coupled branch
isolation and the required sign queries: a five-polynomial rational Krawczyk
certificate jointly isolates `(t_G,q_G,p_G,t_B,q_B)`, and rational interval
implicit differentiation certifies the G/B3 SOC and cross signs.

The remaining Level-3 requirement is different: the projected conditions must
hold as an explicitly computed necessary-and-sufficient characterization over a
declared non-calibrated primitive region. L3-1 has not computed that parametric
projection or uniform root-selection theorem.

Therefore this file records both mathematical existence and exact canonical
closure, without claiming a primitive-space Level-3 success.

## Status

- exact quantified semialgebraic characterization: **YES**;
- canonical coupled root/sign decision: **PASS**;
- necessary primitive conditions: identified as projection obligations;
- sufficient conditions: existing small-`beta` Level-2 theorem survives;
- primitive necessary-and-sufficient formula: **exists in principle, not explicitly computed**;
- primitive-only explicit condition set: **NO**;
- uniform primitive branch-selection theorem: **NO**;
- unique scalar threshold: **NO**.
