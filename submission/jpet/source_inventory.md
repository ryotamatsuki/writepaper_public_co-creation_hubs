# JPET Source Inventory — QA Candidate

Classification is for the **submission package**, not for repository deletion.

## REQUIRED — submission / compilation

| Path | Role |
|---|---|
| `paper/main.tex` | Main LaTeX document |
| `paper/preamble.tex` | Required LaTeX preamble/macros |
| `sections/*.tex` | Main text and Appendix source included by `main.tex` |
| `references/references.bib` | Bibliography database |
| `generated/tables/*.tex` | Generated tables included by manuscript |
| `Title_Page.docx` | Existing title-page / author-metadata authority |
| `submission/jpet/qa_candidate_manuscript.pdf` | Stage-14 peer-review PDF candidate, generated from clean build |
| `submission/jpet/cover_letter.md` | Cover-letter content; convert/map as live portal requires |
| `submission/jpet/declarations.md` | Submission declaration source |
| `submission/jpet/metadata.md` | Portal metadata source |

## REQUIRED / SUPPORTING — reproducibility subset

These files support claims made in the paper and should be preserved in the Stage-15 frozen research record. Whether each is uploaded as Supporting Information or retained as an external reproducibility record must follow the live portal and final data/code statement.

- `scripts/verify_freeze.py`
- `scripts/verify_symbolic.py`
- `scripts/verify_numerical.py`
- `scripts/generate_results.py`
- `scripts/generate_tables.py`
- `scripts/generate_figures.py`
- `scripts/validate_bibliography.py`
- `scripts/validate_manuscript.py`
- `scripts/validate_build_log.py`
- `tests/`
- `requirements.txt`
- `analytic_level3/code/derive_small_beta.py`
- `analytic_level3/code/verify_threshold.py`
- `analytic_level3/code/verify_symbolic_identities.py`
- `analytic_level3/code/verify_witness.py`
- `analytic_level3_reopen/code/verify_l31.py`
- canonical machine-readable generated result/manifest/verification outputs created by the clean Stage-14 build

## OPTIONAL — submission convenience only

- `submission/jpet/README.md`
- `submission/jpet/source_inventory.md`
- `submission/jpet/reviewer_candidates.md`
- `docs/STAGE_14_JPET_REQUIREMENTS_LEDGER.md`

These are QA/provenance aids and normally should **not** be uploaded to JPET as manuscript supporting information.

## INTERNAL — EXCLUDE FROM JOURNAL UPLOAD

- `reviews/`
- `docs/STAGE_12_*`
- `docs/STAGE_13_*`
- other Stage reports / decision logs / journal ladders
- hostile-referee or editor-attack reports
- killed-journal analyses
- `.github/`
- internal workflow configuration and CI metadata
- exploratory Level-3 prose reports under `analytic_level3/` and `analytic_level3_reopen/` that are not specifically identified above as verification code
- branch/PR metadata
- local build logs and temporary artifacts not explicitly frozen for provenance

## SENSITIVE — MUST EXCLUDE

No submission archive may contain:

- API keys, PATs, access tokens, credentials, cookies or secrets;
- local machine paths or environment secrets;
- private correspondence;
- unrelated personal information beyond the author/contact fields required by the journal title page;
- unrelated repository or research-project material.

Stage-14 automation inspects `Title_Page.docx` for required-field presence but deliberately does not print personal values into logs or this public QA inventory.

## Figure note

The current production pipeline does not include a standalone research figure in the manuscript; `generated/figures/README.md` is a generated provenance marker. Tables are generated and included from `generated/tables/*.tex`. If the final clean build regenerates a figure used by the manuscript, it must be reclassified REQUIRED before Stage 15.
