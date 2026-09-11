# Stage 4A v2.1 — Astra Global-Evidence Repair

Date: 2026-09-10 JST

Trigger: Astra Stage-11 hostile review of Stage-10 baseline `02bcd37ac1ab019dc6df2c0d29c8f312a1c992a4`.

Inherited canonical chain before this repair: Stage-11 merge `7fdb1020fe398e5bd6833f9faa5da56812ad345e`.

## Issue

The active Stage-4A evaluator clipped support participation mass to `[0,1]` and then used the interior shortcut `m_h^2/4` as each region's support surplus over the whole public-deviation domain. Under the actual primitive support utility `r_h-c`, with `c~U[0,1]`, that shortcut is valid only when `0<r_h<1`. In a saturation region, gross benefit `r_h` continues to rise while mass is capped at one, so the regional support surplus is one half of `r_h m_h-m_h^2/2`, not `m_h^2/4`.

This made some large-deviation public payoffs wrong even though the reported on-path candidate is support-interior.

A second issue is evidentiary. The Stage-4A procedure uses finite grids, local scalar refinement, and multiple fixed-point starts. It does not derive a rigorous regret bound covering all unsearched public/private points or prove enumeration of all continuation roots. Therefore its old label `computational global-equilibrium certification` exceeded what the implementation certified.

## Repair

`stage4a_v21_repaired/code/independent_repaired_audit.py` now:

- keeps gross support benefit `r_h` separate from clipped participation mass `m_h`;
- computes national support surplus as `r_h*m_h-m_h^2/2` and attributes one half to each region;
- regression-tests the `r<=0`, `0<r<1`, and `r>=1` branches and continuity at `0` and `1`;
- targets support-saturation neighborhoods in addition to boundaries and the old low-investment region;
- verifies participation residuals and numerical private/public FOC/SOC at the retained reported candidates;
- continues to reoptimize the private price after G deviations and holds the matched scalar fee fixed during B3 deviations;
- explicitly labels the result `SEARCH EVIDENCE` and records that no certified global regret bound is available.

## Claim effect

The repaired vector is not changed and no new vector is searched.

T1 and T2 are unchanged. They are local analytic identities in the support-interior central regime.

T3 remains the local stationary-branch theorem subject to the Stage-7.5A/Stage-8 symmetry amendment and the explicit support-interiority/routing/SOC conditions.

The old T4 label `global-equilibrium existence witness/certification` is withdrawn. The replacement statement is:

> At the existing repaired vector, a primitive-consistent all-regime numerical search detects no profitable public deviation under the documented grid, local-refinement, private-reoptimization/fixed-price, and multistart protocol and reproduces the local sign comparison. This is SEARCH EVIDENCE; no certified global regret bound, equilibrium-existence proof, or uniqueness proof is claimed.

Stage 12 remains blocked pending Astra limited recheck.
