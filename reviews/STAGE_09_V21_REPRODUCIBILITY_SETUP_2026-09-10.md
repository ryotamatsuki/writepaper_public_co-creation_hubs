# Stage 9 v2.1 — Repository / Reproducibility Setup

Date: 2026-09-10 JST  
Workflow authority: `ryotamatsuki/research-paper-workflow` v2.1 @ `b27568172d4145a2b98825a5b66a071dbfe25f36`  
Template: `templates/STAGE_09_REPRODUCIBILITY_SETUP.md` @ blob `b22be041db59283ae9feb70b122c9eb73969c898`  
Final Stage-7.5A certified input: `eeb48a3dd76ab6f43d5de175b12c3f374746db0f`  
Final Stage-8 merge / Stage-9 parent: `ad927ca783a6123ea4fc6f55f65598ebd6ab583b`

## 1. Starting remote state

- default `main` at Stage-9 start: `3ee1f0aadd18d4cd5d787ce6bb9bf68460a40d71`; it is an older production chain and is not the v2.1 theory authority;
- PR #45 is merged into the Stage-7.5A decision branch, producing final Stage-8 merge `ad927ca783a6123ea4fc6f55f65598ebd6ab583b`;
- open PR #43 preserves the earlier negative Stage-4A record; open PR #42 preserves the older pre-v2.1 Stage-14 package; neither is overwritten;
- stale branch `pipeline/v21-stage9-reproducibility` is retained as history. It diverges from the final freeze and descends from pre-final Stage-7.5A input `91849e2893a5599fcd64e8bbbb38724dcf4e5c67`;
- canonical Stage-9 branch is recreated directly from the final Stage-8 merge.

## 2. Repository structure / production sources

Existing modular production structure is retained and re-authorized only where consistent with the final freeze:

- `paper/` — modular LaTeX entry point and preamble;
- `sections/` — manuscript sections and appendix;
- `references/` — bibliography source;
- `scripts/` — freeze, theorem, numerical, scope, bibliography, generation, build, report gates;
- `tests/` — deterministic pipeline regressions;
- `generated/results/` — machine-readable canonical result/report/manifest layer;
- `generated/tables/` — machine-generated LaTeX tables;
- `generated/figures/` — explicit figure registry;
- `stage4a_v21_repaired/` — independent all-regime global-equilibrium certification;
- `analytic_level3/` — analytic small-beta derivations and symbolic identities;
- `reviews/` — certification and claim-scope records;
- `theory_freeze_v21/` — final Stage-8 freeze;
- `docs/` — provenance, environment, reproducibility and Stage-10 contract;
- `.github/workflows/` — clean CI gate;
- `Makefile` / `requirements.txt` — one-command build and pinned Python environment.

## 3. Build system

Canonical command:

```sh
make clean && make all
```

The `all` target runs freeze consistency, analytic symbolic identities, the independent all-regime global audit, repaired numerical verification, Stage-7.5A claim-scope regression, bibliography/manuscript audits, pytest, deterministic tables/figure registry, LaTeX build, verification report, and generated-object SHA-256 manifest.

## 4. Verification scripts / tests

Canonical gates:

- `scripts/verify_freeze.py` — final Stage-8 SHA/parameter/claim-ceiling identity;
- `scripts/verify_symbolic.py` — local small-beta theorem identities + fee-transfer identity;
- `stage4a_v21_repaired/code/independent_repaired_audit.py` — all-route fail-closed global-equilibrium witness audit;
- `scripts/generate_results.py` — repaired machine-readable result layer;
- `scripts/verify_numerical.py` — repaired witness, slopes, welfare witness, local coordination wedge;
- `scripts/stage75a_scope_audit.py` — theorem quantifier/globality/overclaim regression;
- `scripts/validate_bibliography.py` and `scripts/validate_manuscript.py` — source/manuscript integrity;
- `tests/test_pipeline.py` — provenance, generated signs/tables, scope, and rejected-witness exclusion.

## 5. Theorem-certificate / claim-scope locations

- Stage-8 freeze: `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md`;
- Stage-7.5A scope ledger: `reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md`;
- repaired Stage-4A certification: `reviews/STAGE_04A_V21_REPAIRED_GLOBAL_CERTIFICATION_2026-09-08.md` and `stage4a_v21_repaired/`;
- analytic theorem suite: `analytic_level3/code/derive_small_beta.py` and `analytic_level3/code/verify_symbolic_identities.py`;
- claim-to-source crosswalk: `docs/CLAIM_SOURCE_MAP.md`.

## 6. Counterexample / regression preservation

The rejected low-friction vector `(beta,gamma,tau)=(.05,.9,.05)` remains permanent negative evidence. The independent Stage-4A evaluator explicitly attacks the old dangerous low-investment region and treats `UNRESOLVED` / `MULTIPLE_EQUILIBRIA` as fail-closed states. The old 20-draw perturbation exercise is explicitly excluded from active robustness authority. The old exact Krawczyk object remains local stationary-root evidence only.

## 7. Environment / dependencies

`requirements.txt` pins NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0, and pytest 9.0.2. Canonical CI uses Python 3.12 and installs `latexmk` plus standard TeX Live LaTeX packages. No external data, secrets, or private service is required.

## 8. Figure / table pipeline

`scripts/generate_tables.py` regenerates four LaTeX tables from `generated/results/canonical_results.json`; no numerical cell is maintained manually. `scripts/generate_figures.py` records that the frozen result set requires no production figure at Stage 9 rather than introducing a non-substantive graphic. `scripts/generate_manifest.py` hashes generated results, tables, and figure registry.

## 9. CI / local-equivalent gate

GitHub Actions workflow `stage9-v21-reproducibility-final` checks:

1. ancestry from `ad927ca783a6123ea4fc6f55f65598ebd6ab583b`;
2. pinned Python / LaTeX environment;
3. `make clean && make all`;
4. a second clean generated-layer regeneration with manifest byte comparison;
5. post-regeneration Stage-7.5A scope audit.

Status at initial Stage-9 commit: **PENDING PR-HEAD CI**.

## 10. Provenance locations

- `docs/PROVENANCE.md`
- `docs/REPRODUCIBILITY.md`
- `docs/REPRODUCIBILITY_DECISION_LOG.md`
- `docs/CLAIM_SOURCE_MAP.md`
- `generated/results/verification_report.json`
- `generated/results/manifest.json`

## 11. Remaining blockers

No theory blocker is open. The only initial closeout dependency is successful execution of the clean PR-head CI gate. A platform/package-install failure must be distinguished from a repository verification failure.

## 12. Exact Stage-10 writing contract

See `docs/STAGE_10_WRITING_CONTRACT.md`. Stage 10 is limited to exposition/presentation against the final Stage-8 freeze and Stage-7.5A quantifier ledger. It must preserve the distinction between the local analytic theorem and the one-vector repaired computational global-equilibrium witness, must not revive the rejected vector/20-draw robustness/exact stationary certificate as global authority, and must source all reported numerical values from generated objects.

## Provisional verdict

**CONDITIONAL GO — PR-HEAD CI REQUIRED FOR FINAL `REPRODUCIBILITY BASELINE READY`.**
