# Public Innovation Hubs — Stage 11R Astra Repair Report

Date: 2026-09-10 JST

## 1. Provenance and branch policy

Astra audited baseline: `02bcd37ac1ab019dc6df2c0d29c8f312a1c992a4` (Stage-10 merge).

At repair start, repository default `main` remained on the older production line at `3ee1f0aadd18d4cd5d787ce6bb9bf68460a40d71`; it was not reset or overwritten. The current v2.1 authority line was `decision/v21-stage75-full-theory-freeze` at Stage-11 merge `7fdb1020fe398e5bd6833f9faa5da56812ad345e`. Stage 11 had already repaired the T3 symmetric beta-zero G anchor and added Kim (2024) to the closest-prior-art discussion. Those repairs were preserved and not duplicated.

Dedicated repair branch:

`repair/v21-stage11r-astra-global-welfare`

This Stage 11R is not authorized to merge itself or to declare Stage 12 open. It returns to Astra limited recheck.

## 2. R1 — saturated support-side surplus

### Issue

The Stage-4A evaluator used clipped participation mass `m_h=clip(r_h,0,1)` and then applied `m_h^2/4` as each region's support-side surplus contribution at all histories. That expression is valid only where `0<r_h<1` and `m_h=r_h`.

### Repair

`stage4a_v21_repaired/code/independent_repaired_audit.py` now separates:

- gross support benefit `r_h`;
- clipped support mass `m_h=clip(r_h,0,1)`;
- primitive national surplus `r_h m_h-m_h^2/2`.

Each region receives one half of that national surplus. Thus the scalar national support surplus is

- `0` for `r<=0`;
- `r^2/2` for `0<r<1`;
- `r-1/2` for `r>=1`.

No gross benefit is clipped and the unit support-population cap is retained.

`public_two_sided_platform_welfare_generality/PARTNER_SURPLUS_DERIVATION.md` and `sections/04_welfare.tex` now state the general formula and identify the squared formula as the support-interior specialization only.

### Regression

The evaluator and independent Stage-11R implementation test the interior branch, the saturated branch, and continuity at `r=0` and `r=1`.

At the reported G candidate the support side remains interior, so the candidate's local welfare and local derivatives are not mechanically changed by this repair. Large upward deviations do enter saturation. In particular, the corrected evaluator reproduces approximately:

| own x, rival fixed at x_G | pre-repair coded W_i | corrected W_i |
|---:|---:|---:|
| `x_G≈0.8371022382` | `0.2948786590` | `0.2948786590` |
| `0.9` | `0.2552099352` | `0.2827613692` |
| `1.0` | `0.1768349352` | `0.2543863692` |

These values are regression diagnostics, not a proof of global optimality.

**R1 status: PASS subject to final CI regeneration.**

## 3. R2 — existing vector and global evidence level

Existing vector retained, with no new parameter search:

`beta=.01, gamma=.825, tau=.35`.

Stored candidate states retained for direct revalidation rather than assumed optimality:

- `x_G≈0.8371022382025995`;
- `x_B3≈0.8258903860237495`;
- matched `p_G≈0.018403679612460814`.

The repaired Stage-4A evaluator now explicitly reports:

- on-path participation residual;
- private-price numerical FOC and SOC;
- public own numerical FOC and SOC in G and B3;
- local cross derivatives / BR slopes;
- large off-path histories including the support-saturation boundary and saturated region;
- multiple participation starts;
- endpoint and low-investment targets;
- best detected deviation gain from the documented grid + local-refinement search.

### Certification decision

The current numerical architecture does **not** provide a mathematically justified bound on the supremum of deviation gain between searched grid/refinement points, nor an interval/exhaustive certificate for every private-price and participation continuation. Increasing grid density is not treated as a certificate.

Therefore the controlling evidence level is:

**SEARCH EVIDENCE**

and not:

- certified approximate equilibrium bound;
- verified exact equilibrium existence;
- computational global-equilibrium certification.

The generated layer must set `certified_regret_upper_bound` to `null` and state `NO CERTIFIED GLOBAL REGRET BOUND`.

The Stage-8 historical T4 wording is retained as historical evidence and is superseded by `theory_freeze_v21/STAGE_08_V21_ASTRA_WELFARE_EVIDENCE_AMENDMENT_2026-09-10.md` for T4/W2/W3 qualification. `reviews/STAGE_04A_V21_ASTRA_GLOBAL_EVIDENCE_REPAIR_2026-09-10.md` records the Stage-4A claim-level repair.

**R2 status: PASS as scope repair if final search/regression passes; T4 old certification claim is withdrawn.**

## 4. R3 — T3 quantifier / symmetric anchor

The prior canonical Stage 11 had already identified the T3 certification regression and created:

- `reviews/STAGE_075A_V21_T3_SYMMETRY_LIMITED_REOPEN_2026-09-10.md`;
- `theory_freeze_v21/STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md`.

Stage 11R preserves that repair and tightens reader-facing wording. The theorem now requires the G branch to pass through the symmetric regular beta-zero **central-interior** stationary state covered by T2, together with support-side interiority, strict routing inequalities, nonsingularity, and strict private/public SOCs.

The matched-price path is stated as

`bar p(beta)=p_G(x_G(beta), beta)`

with the scalar value fixed within each B3 deviation problem. The B3/G continuity argument takes the common neighborhood over branch existence, support/routing regime, SOCs, nonsingularity, smooth matching, and sign preservation; equivalently epsilon may be chosen as the minimum of the corresponding positive neighborhood radii.

No asymmetric theorem or stronger parameter-region result is added.

**R3 status: PASS subject to scope CI.**

## 5. R4 — B3 interpretation

`sections/02_equilibrium.tex` and downstream exposition now distinguish two different comparisons:

1. G and B3 use the same scalar private fee at the matched comparison, while public investments and allocations reoptimize within their respective environments, so generally `x_B3 != x_G`.
2. The same-state chain-rule decomposition switches the local price-policy response on/off at a fixed state; it is not identical to comparing the two distinct stationary states.

B3 remains an identification benchmark. It is neither a planner problem nor a literal price-regulation counterfactual. No alternative benchmark was introduced.

**R4 status: PASS subject to manuscript CI.**

## 6. R5 — welfare derivative generation

`scripts/generate_results.py` no longer constructs the national coordination wedge as `rival derivative + private-profit derivative` while silently setting the own derivative to zero.

It now numerically computes, along the direction that changes own public investment while fixing rival public investment and resolving the downstream private-price/participation continuation:

- own reduced regional-welfare derivative;
- rival reduced regional-welfare derivative;
- private-profit derivative;
- direct national-welfare derivative;
- decomposition sum;
- decomposition error.

It also records the finite-difference step. `scripts/verify_numerical.py` requires the own derivative to be close to zero as a stationarity residual and requires the direct national derivative to agree with the three-part decomposition within tolerance.

Expected scale at the reported G state remains approximately:

- own derivative: near `0`;
- rival derivative: `0.500526`;
- private-profit derivative: `-0.010080`;
- national derivative: `0.490447`.

These are generated quantities, not hard-coded returned values.

Because R2 downgrades T4, W2 is described as a local directional welfare statement at the reported numerical G state, and W3 as a numerical G/B3 state comparison. Neither is qualified as a globally certified equilibrium-welfare result pending any stronger numerical certificate.

**R5 status: PASS subject to generated-output CI.**

## 7. R6 — reproducibility / stale manifest

The build dependency chain remains `results -> tables/figures -> manifest`. Stage 11R adds `scripts/verify_manifest.py`, which recomputes byte counts and SHA-256 hashes for every manifest target and fails on any mismatch. `make all` now includes `manifest-verify`.

Dedicated Stage-11R CI performs:

1. `make clean && make all`;
2. `python scripts/stage75a_scope_audit.py`;
3. independent Stage-11R regression;
4. manifest integrity verification;
5. clean deterministic regeneration and manifest comparison;
6. committed-generated-file diff gate.

The initial PR run is intentionally allowed to expose stale generated files as a failure; the generated diff is then committed and the gate rerun until clean. Build success and numerical guarantee level are reported separately.

**R6 status: PENDING final generated-file synchronization and CI.**

## 8. Claims after repair

### Maintained analytically

- T1: beta-zero B3 cross effect zero and positive first beta derivative, at the stated local support-interior regular branch scope.
- T2: negative beta-zero full-game cross effect at the stated symmetric regular central-interior state.
- T3: local stationary-branch sign reversal for sufficiently small positive beta under the repaired symmetry/interiority/routing/nonsingularity/SOC assumptions.
- W1: aggregate private-fee transfer cancellation identity.
- T5: `kappa_L+tau>v+alpha` as a sufficient remote-public route-dominance condition.

### Changed / withdrawn / narrowed

- T4: **withdrawn as a computational global-equilibrium certification**. Replaced by all-regime computational **search evidence** at the existing vector, with no certified global regret bound.
- W2: remains a local numerical directional welfare statement at the reported state; it is not promoted by T4 to a certified equilibrium welfare theorem.
- W3: remains a numerical ranking of the two reported states only.
- B3 explanation is narrowed to the exact matched-scalar-fee identification role.

No new player, timing, primitive utility, strategy set, parameter vector, benchmark structure, or analytic theorem has been introduced.

## 9. Limited prior-art closure

See `reviews/STAGE_11R_LIMITED_PRIOR_ART_CLOSURE_2026-09-10.md`.

- Mun & Nakagawa (2010): accessible author-hosted/primary material inspected; no direct matched-price public-public sign-reversal object identified.
- He et al. (2026): accessible publisher material inspected; generic pricing/investment/network/platform-competition ingredients are absorbed, but the exact public-public matched-price repricing comparison is not identified.
- Mun (2019): published abstract and public 2016 precursor inspected; final published full text was not available through the legal access route used here. Residual proposition-level uncertainty remains and is carried to Astra. Lack of access is not treated as non-absorption proof.

## 10. Files changed by Stage 11R

Core implementation / verification:

- `stage4a_v21_repaired/code/independent_repaired_audit.py`
- `stage11_v21_independent/code/independent_stage11_regression.py`
- `scripts/generate_results.py`
- `scripts/generate_tables.py`
- `scripts/verify_numerical.py`
- `scripts/verify_freeze.py`
- `scripts/stage75a_scope_audit.py`
- `scripts/validate_manuscript.py`
- `scripts/validate_stage10_architecture.py`
- `scripts/verify_manifest.py`
- `tests/test_pipeline.py`
- `Makefile`

Model / authority records:

- `public_two_sided_platform_welfare_generality/PARTNER_SURPLUS_DERIVATION.md`
- `reviews/STAGE_04A_V21_ASTRA_GLOBAL_EVIDENCE_REPAIR_2026-09-10.md`
- `reviews/STAGE_07_V21_SATURATED_PARTNER_SURPLUS_CORRECTION_2026-09-10.md`
- `reviews/STAGE_11R_LIMITED_PRIOR_ART_CLOSURE_2026-09-10.md`
- `theory_freeze_v21/STAGE_08_V21_ASTRA_WELFARE_EVIDENCE_AMENDMENT_2026-09-10.md`
- `docs/CLAIM_SOURCE_MAP.md`
- `docs/STAGE_10_WRITING_CONTRACT.md`
- `docs/STAGE_10_EXPOSITION_ARCHITECTURE.md`

Reader-facing manuscript:

- `paper/main.tex`
- `sections/02_equilibrium.tex`
- `sections/03_main_results.tex`
- `sections/04_welfare.tex`
- `sections/05_robustness.tex`
- `sections/08_introduction.tex`
- `sections/09_discussion.tex`
- `sections/10_conclusion.tex`
- `sections/appendices.tex`

CI:

- `.github/workflows/stage11r-astra-repair.yml`
- `.github/workflows/stage11-v21-referee.yml` (historical Stage-11 job excluded for this superseding repair branch).

Generated objects will be listed after synchronization.

## 11. Final verification commands

Required final commands/gates:

- `make clean && make all`
- `python scripts/stage75a_scope_audit.py`
- `python stage11_v21_independent/code/independent_stage11_regression.py`
- `python scripts/verify_manifest.py`
- clean deterministic regeneration + manifest comparison
- committed generated-file diff gate

The final exit codes and run IDs will be appended after the draft-PR CI completes.

## 12. Astra limited-recheck scope

Astra should recheck only:

1. R1 primitive support-surplus formula and saturated off-path welfare;
2. the existing vector after corrected welfare accounting, including search results and candidate stationarity diagnostics;
3. the explicit distinction between best detected gain and absent certified regret upper bound;
4. T4/W2/W3 evidence-level downgrade and authority supersession chain;
5. T3 common-neighborhood/symmetric-anchor wording;
6. B3 matched-price interpretation;
7. directly generated national-welfare derivative decomposition;
8. regenerated manifest/determinism evidence;
9. remaining Mun (2019) full-text prior-art uncertainty.

No Stage-12 decision is made here.

## Provisional repair verdict

**PARTIAL REPAIR — R6 GENERATED-OUTPUT / CI CLOSURE PENDING.**

This verdict must be updated only after final PR-head CI and generated-output synchronization.
