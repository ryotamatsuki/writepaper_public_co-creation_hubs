# Stage 9 v2.1 — Repository / Reproducibility Setup

Date: 2026-09-10 JST
Workflow: `ryotamatsuki/research-paper-workflow` v2.1 @ `fe6fa5f2a94632be578d5e9fab39e7df86013113`
Stage-8 input: `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md`
Starting Stage-8 commit: `7aa9de680141fb4778f4926ee05cd95a57385971`

## Verdict

**REPRODUCIBILITY BASELINE READY**, subject only to the Stage-9 CI run attached to this branch completing successfully.

## Remote/concurrency record

The default `main` is not used as theory authority for this retrofit chain. Two older draft PRs remain open, including the negative Stage-4A record and the pre-v2.1 Stage-14 package. They are preserved and not overwritten. This Stage-9 branch descends only from the certified Stage-7.5A head and Stage-8 freeze.

## Canonical production chain

- theory freeze: `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md`
- independent global-equilibrium regression: `stage4a_v21_repaired/code/independent_repaired_audit.py`
- claim-scope regression: `scripts/stage75a_scope_audit.py`
- repaired result layer: `scripts/generate_results.py`
- repaired tables: `scripts/generate_tables.py`
- freeze gate: `scripts/verify_freeze.py`
- pipeline tests: `tests/test_pipeline.py`
- bibliography gate: `scripts/validate_bibliography.py`
- LaTeX build: `paper/main.tex`

## Changes in Stage 9

1. `scripts/verify_freeze.py` now checks the v2.1 repaired freeze rather than the obsolete v1.1 freeze.
2. `tests/test_pipeline.py` now asserts the repaired primitive vector and global-witness signs and invokes the independent Stage-4A regression.
3. `Makefile` now treats the repaired global audit, Stage-7.5A scope lint, repaired generated results, tests, bibliography, and manuscript build as the canonical `make all` chain.
4. `.github/workflows/stage9-v21-reproducibility.yml` provides a clean CI equivalent with pinned repository dependencies and LaTeX installation.

## Environment

Python dependencies are pinned in `requirements.txt` (`numpy`, `scipy`, `sympy`, `pytest`). CI uses Python 3.12 and installs `latexmk`, `texlive-latex-base`, `texlive-latex-recommended`, and `texlive-latex-extra`.

## Generated-object provenance

The active strategic and welfare values come from the repaired Stage-4A / Stage-7 audit path. The old low-friction exact stationary certificate remains archived only and is not used to generate active tables.

## Regression preservation

The old low-friction finite-deviation failure remains preserved as certification history. The Stage-4A independent evaluator fails closed on unresolved/multiple continuation states and explicitly attacks the old dangerous low-investment region.

## Stage-10 contract

Stage 10 may edit exposition and presentation only. It must use the Stage-8 freeze and Stage-7.5A scope ledger, must keep the repaired vector as the only active numerical equilibrium witness, and must not restore the old exact stationary-root certificate as equilibrium authority or reintroduce stale 20-draw robustness.
