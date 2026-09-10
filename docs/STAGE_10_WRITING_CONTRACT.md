# Stage 10 Writing Contract

Authority hierarchy:

1. `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md` at final Stage-8 merge `ad927ca783a6123ea4fc6f55f65598ebd6ab583b`, together with `theory_freeze_v21/STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md` for T3 and `theory_freeze_v21/STAGE_08_V21_ASTRA_WELFARE_EVIDENCE_AMENDMENT_2026-09-10.md` for all-regime welfare/T4/W2/W3 qualification;
2. `reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md`, together with `reviews/STAGE_075A_V21_T3_SYMMETRY_LIMITED_REOPEN_2026-09-10.md` for T3 and `reviews/STAGE_07_V21_SATURATED_PARTNER_SURPLUS_CORRECTION_2026-09-10.md` for capped support-side welfare accounting;
3. Stage-11R verified/generated result layer and search-evidence record;
4. retained verified literature and institutional evidence.

Stage 10 exposition remains subordinate to the amended theory/evidence authority. It may edit organization, notation presentation, citations, and generated-object inclusion only. It may not change players, timing, utilities, route set, objectives, parameter restrictions, benchmark definitions, theorem quantifiers, welfare object, novelty envelope, or upgrade numerical evidence beyond the controlling amendments.

Mandatory claim discipline:

- T1/T2/T3 are analytic **local sufficient-condition results on regular stationary branches**. T3 requires the G branch to pass through the **symmetric regular beta-zero central-interior full-game stationary state covered by T2**, with support-side interiority, strict routing inequalities, nonsingular local systems, the matched-price path, and strict private/public SOCs continuing on a common neighborhood. Do not call these global-equilibrium, uniqueness, asymmetric-G, or primitive-space characterization theorems.
- T4 is now **ALL-REGIME COMPUTATIONAL SEARCH EVIDENCE — NO CERTIFIED GLOBAL REGRET BOUND** at the existing primitive vector. Permitted wording is that no profitable deviation was detected under the documented search. Do not describe the reported states as certified global equilibria, verified equilibrium-existence witnesses, exact global best responses, or uniqueness results.
- B3 is the **matched-price fixed-price identification benchmark**. Along the G stationary path, it uses the same scalar fee selected in G and holds that fee fixed during each B3 public deviation; B3 public choices and allocations otherwise readjust, so `x_B3` and `x_G` generally differ. Distinguish this cross-environment stationary comparison from the same-state chain-rule decomposition. Do not describe B3 as a planner problem or literal price-regulation policy.
- Support-side welfare uses gross participation benefit `r_h` and clipped mass `m_h=clip(r_h,0,1)`. National support surplus is `r_h*m_h-m_h^2/2`, half of which is attributed to each region. The shortcut `m_h^2/4` is valid only under support-side interiority.
- W1 is an exact transfer-cancellation accounting identity.
- W2 is a **local coordination-wedge numerical statement at the reported G stationary candidate**. The own reduced derivative must be computed, and the direct national derivative must agree with the component decomposition. Do not claim first best, global social optimum, optimal subsidy, or global underinvestment.
- W3 is a **one-state-pair numerical G-versus-B3 welfare ranking only**. Do not qualify it as a ranking of certified or unique equilibria and do not claim welfare dominance.
- The sufficient route-dominance condition `kappa_L+tau > v+alpha` is not necessary and the repaired numerical exercise is not direct public-hub user-poaching evidence.
- The old vector `(beta,gamma,tau)=(.05,.9,.05)` is rejected as global-equilibrium authority; the old 20-draw robustness is non-authoritative; the old exact Krawczyk certificate is only a local stationary-root diagnostic.
- No arbitrary-distribution, arbitrary nonlinear-network, heterogeneous-region, global primitive-space, uniqueness, equilibrium-existence, certified-regret, or broad genericity claim may be introduced without a formally governed earlier-stage reopening.
- The surviving novelty claim is the **matched-price public-public local strategic sign reversal when private repricing is activated**, not follower pricing, two-sidedness, public-private competition, or downstream reactions in general.
- The private-price response remains a model prediction, not observed causal evidence.

Every quantitative value included in prose or tables must be sourced from `generated/results/canonical_results.json` or generated LaTeX tables. Any new numerical object requires a reproducible generator and test before inclusion.

The Stage-11 T3 certification regression narrowed T3 to the beta-zero G state for which T2 is actually proved. The Stage-11R Astra repair additionally corrects off-path support surplus and lowers the T4 numerical evidence level. Both corrections are controlling downstream; neither authorizes new theory.

Before downstream closeout, run `make clean && make all`, `python scripts/stage75a_scope_audit.py`, and `python scripts/verify_manifest.py`. Stage 12 remains blocked until Astra limited recheck resolves the Stage-11R repair.
