# Stage 8 v2.1 — Astra Welfare / Numerical-Evidence Amendment

Date: 2026-09-10 JST

Historical Stage-8 freeze: `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md` at merge `ad927ca783a6123ea4fc6f55f65598ebd6ab583b`.

Subsequent T3-only amendment: `STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md`.

Trigger: Astra Stage-11 hostile review of Stage-10 baseline `02bcd37ac1ab019dc6df2c0d29c8f312a1c992a4`.

This amendment preserves the historical freeze and supersedes only the all-regime welfare-accounting authority and the evidence qualification of T4/W2/W3. It does not add a theorem, parameter vector, primitive, equilibrium-selection assumption, or welfare objective.

## A. Welfare-accounting amendment

The historical all-regime evaluator used the regional shortcut `m_h^2/4` after clipping support participation mass. The primitive model instead has gross support benefit `r_h`, mass `m_h=clip(r_h,0,1)`, and national support surplus

`S_h = r_h*m_h - m_h^2/2`.

Each region receives `S_h/2`. The shortcut `m_h^2/4` is controlling only when support participation is interior and `m_h=r_h`. Off-path saturation histories are henceforth governed by the general primitive integral documented in `reviews/STAGE_07_V21_SATURATED_PARTNER_SURPLUS_CORRECTION_2026-09-10.md`.

## B. T4 evidence-level amendment

Historical wording that calls the repaired vector an `all-regime computational global-equilibrium existence witness` or `global-equilibrium certification` is no longer controlling.

The controlling T4 numerical evidence level is:

**ALL-REGIME COMPUTATIONAL SEARCH EVIDENCE — NO CERTIFIED GLOBAL REGRET BOUND.**

At the existing vector `(beta,gamma,tau)=(.01,.825,.35)`, the corrected evaluator searches the public interval, reoptimizes private price after G deviations, holds the matched scalar fee fixed during B3 deviations, reconstructs all routing alternatives, uses multiple participation starts, and targets boundaries, low-investment histories, and support-saturation neighborhoods. If no profitable deviation is detected, the permitted statement is only that no profitable deviation was detected under the documented search.

Finite grids, local scalar refinement, solver success, and convergence from multiple starts are not a proof of a global regret upper bound or exhaustive enumeration of all continuation roots. Unless a separate rigorous bound is later constructed through governance, the vector must not be described as a certified global equilibrium, verified equilibrium existence, exact global best response, or uniqueness result.

## C. Welfare evidence amendment

W1 remains the exact fee-transfer cancellation identity.

W2 remains a local directional numerical statement at the reported G stationary candidate. The own reduced derivative must be generated numerically rather than assumed zero; the national derivative must also be calculated directly and checked against

`dW_N/dx_i = dW_i/dx_i + dW_j/dx_i + dPi_T/dx_i`.

W3 remains a numerical G-versus-B3 welfare ranking at the reported state pair only. Because T4 is now search evidence rather than globally certified equilibrium evidence, W2/W3 must not acquire an `equilibrium` qualification from the old T4 label.

## D. Analytic claims unaffected

T1 and T2 are unchanged at their local support-interior scope.

T3 is controlled jointly by the historical freeze, `STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md`, and the manuscript repair that explicitly maintains support-side interiority, strict central-routing inequalities, nonsingularity, strict SOCs, smooth branch continuation, and the matched-price path. The G branch remains anchored at the symmetric regular beta-zero central-interior state covered by T2.

T5 remains the sufficient nonresident rival-public route-dominance condition only.

## E. Downstream gate

Stage 12 is blocked until Astra limited recheck of this Stage-11R repair. The repository may record successful reproducibility and search checks, but no downstream journal-positioning authorization follows automatically from those checks.
