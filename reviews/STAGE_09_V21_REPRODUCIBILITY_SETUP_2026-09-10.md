# Stage 9 v2.1 — Repository / Reproducibility Setup

Date: 2026-09-10 JST  
Workflow authority: `ryotamatsuki/research-paper-workflow` v2.1 @ `b27568172d4145a2b98825a5b66a071dbfe25f36`  
Template: `templates/STAGE_09_REPRODUCIBILITY_SETUP.md` @ blob `b22be041db59283ae9feb70b122c9eb73969c898`  
Final Stage-7.5A certified input: `eeb48a3dd76ab6f43d5de175b12c3f374746db0f`  
Final Stage-8 merge / Stage-9 parent: `ad927ca783a6123ea4fc6f55f65598ebd6ab583b`

## 1. Starting remote/main SHA

Default `main` at Stage-9 start was `3ee1f0aadd18d4cd5d787ce6bb9bf68460a40d71`, an older production chain and not the v2.1 theory authority. PR #45 had been merged into `decision/v21-stage75-full-theory-freeze`, producing the authoritative Stage-8 merge `ad927ca783a6123ea4fc6f55f65598ebd6ab583b`. The stale `pipeline/v21-stage9-reproducibility` branch descended from pre-final Stage-7.5A input `91849e2893a5599fcd64e8bbbb38724dcf4e5c67`; it was preserved as history and not reset. Open legacy PRs were not overwritten.

Canonical Stage-9 branch: `pipeline/v21-stage9-reproducibility-final`.  
Stage-9 PR: #46.

## 2. Repository tree

Production sources are organized under `paper/`, `sections/`, `references/`, `scripts/`, `tests/`, `generated/results/`, `generated/tables/`, `generated/figures/`, `stage4a_v21_repaired/`, `analytic_level3/`, `reviews/`, `theory_freeze_v21/`, `docs/`, and `.github/workflows/`, with root `Makefile` and `requirements.txt`.

## 3. Build system

Canonical command:

```sh
make clean && make all
```

The dependency graph avoids repeated execution of the expensive all-regime audit and generated-result solve. A clean build regenerates the active result layer once, then verifies and consumes that generated object downstream.

## 4. Verification scripts / tests

- `scripts/verify_freeze.py` — final Stage-8 identity, repaired parameters/witness, claim ceiling.
- `scripts/verify_symbolic.py` — v2.1 small-beta theorem identities and fee-transfer identity.
- `stage4a_v21_repaired/code/independent_repaired_audit.py` — independent all-route, fail-closed global-equilibrium witness audit.
- `scripts/generate_results.py` / `scripts/verify_numerical.py` — repaired generated result layer and freeze-consistency checks.
- `scripts/stage75a_scope_audit.py` — theorem quantifier/globality/overclaim regression.
- `scripts/validate_bibliography.py` / `scripts/validate_manuscript.py` — source and manuscript integrity.
- `tests/test_pipeline.py` — provenance, repaired-witness, generated-table, scope, and rejected-evidence regressions.

The manuscript lint accepts semantically equivalent non-calibration wording rather than depending on one exact English phrase; this changes no manuscript claim.

## 5. Theorem-certificate / claim-scope locations

- Stage-8 freeze: `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md`.
- Stage-7.5A scope ledger: `reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md`.
- Repaired Stage-4A certification: `reviews/STAGE_04A_V21_REPAIRED_GLOBAL_CERTIFICATION_2026-09-08.md` and `stage4a_v21_repaired/`.
- Analytic theorem verification: `analytic_level3/code/derive_small_beta.py` and `analytic_level3/code/verify_symbolic_identities.py`.
- Claim/source crosswalk: `docs/CLAIM_SOURCE_MAP.md`.

## 6. Counterexample / regression-test locations

The rejected low-friction vector `(beta,gamma,tau)=(.05,.9,.05)` remains permanent negative evidence. The independent Stage-4A evaluator attacks the dangerous low-investment region and treats `UNRESOLVED` and `MULTIPLE_EQUILIBRIA` as fail-closed states. The old 20-draw perturbation exercise is excluded from active robustness authority. The old exact Krawczyk object remains only a local stationary-root diagnostic.

## 7. Environment / dependencies

`requirements.txt` pins NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0, and pytest 9.0.2. Canonical Stage-9 CI uses Python 3.12 with `latexmk` and standard TeX Live LaTeX packages. No external dataset, secret, or private service is required.

## 8. Figure / table pipeline

`scripts/generate_tables.py` creates four LaTeX tables from `generated/results/canonical_results.json`; numerical cells are not manually maintained. `scripts/generate_figures.py` records explicitly that no production figure is required by the frozen Stage-8 result set. `scripts/generate_manifest.py` records SHA-256 hashes for generated results, tables, and the figure registry.

## 9. CI / local-equivalent gate status

Final substantive PR-head audited: `a7d1c75ddb5a298fa68bc8603b8edd9f951cfab6`.

Stage-9 dedicated workflow:
- run `34411458515`;
- job `102666573365`;
- final Stage-8 ancestry gate: PASS;
- Python environment: PASS;
- LaTeX environment: PASS;
- `make clean && make all`: PASS;
- second-pass generated manifest byte comparison: PASS;
- post-regeneration Stage-7.5A scope audit: PASS;
- conclusion: **SUCCESS**.

Generic `reproducibility` workflow run `34411458537` also completed `make clean && make all` successfully. Stage-7.5A quantifier workflow run `34411458601` completed successfully. Legacy Stage-13 and Stage-14 workflows are now guarded to their corresponding later-stage branches and correctly skip this Stage-9 PR.

The clean gate independently reproduced the repaired all-regime witness, including `STAGE4A_REPAIRED_GLOBAL_CERTIFICATION: PASS`, opposite local slopes in G and B3, and the permanent rival-public dominance check. No theory or scope failure remains.

## 10. Provenance locations

- `docs/PROVENANCE.md`
- `docs/REPRODUCIBILITY.md`
- `docs/REPRODUCIBILITY_DECISION_LOG.md`
- `docs/CLAIM_SOURCE_MAP.md`
- `generated/results/verification_report.json`
- `generated/results/manifest.json`
- this Stage-9 closeout record

## 11. Remaining blockers

**NONE for Stage 9.** No theory, provenance, claim-scope, numerical, build, bibliography, deterministic-generation, or CI blocker remains.

## 12. Exact Stage-10 writing contract

Stage 10 may write and revise exposition/presentation only against the final Stage-8 freeze, the Stage-7.5A quantifier ledger, and the verified Stage-9 generated layer. It must:

1. preserve the distinction between the local sufficient-condition theorem on regular stationary branches and the one-vector repaired computational global-equilibrium existence witness;
2. source reported quantitative values from generated objects rather than hand-edited numbers;
3. keep B3 as a matched-price fixed-price identification benchmark, not a planner or regulation experiment;
4. retain the rejected vector only as negative/regression evidence;
5. not restore the old 20-draw exercise as active robustness evidence or the old Krawczyk object as global-equilibrium authority;
6. not introduce arbitrary-distribution, nonlinear-network, heterogeneous-region, uniqueness, global primitive-space, general welfare-dominance, first-best, or optimal-subsidy claims absent a formally reopened earlier stage;
7. preserve journal selection as deferred to Stage 12.

Canonical contract file: `docs/STAGE_10_WRITING_CONTRACT.md`.

## Theory delta audit

No player, timing, route, utility, objective, benchmark definition, theorem statement/quantifier, equilibrium concept, welfare object, novelty claim, repaired parameterization, or approved evidence envelope changed in Stage 9. Base-to-substantive-head changes are reproducibility engineering, generated-output synchronization, documentation, and CI routing only.

## Final verdict

**REPRODUCIBILITY BASELINE READY**

**STAGE 10 — PAPER BUILD AUTHORIZED UNDER THE CONTRACT ABOVE.**
