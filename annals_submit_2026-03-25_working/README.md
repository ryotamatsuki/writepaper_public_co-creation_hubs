# Annals Submission Working Folder

- Created on: 2026-03-25
- Target journal: The Annals of Regional Science
- Purpose: Clean working area for revising the former JRS submission toward an Annals submission

## Folder structure

- `manuscript_source/`
  Base manuscript source copied from `jrs_submit_2026-03-10/`.
  This is the main working directory for the Annals revision.
- `notes/`
  Revision notes copied from the JRS desk-reject record and the Annals revision TODO.
- `verification/`
  Python / SymPy checks for the model equations, propositions, appendix algebra, and the Section 08 parameter-region figure.
- `reference_docs/`
  Reference-only files from the previous JRS submission stage.
  These are not Annals-ready submission files and should be treated as source material only.

## Main files to edit first

- `manuscript_source/main.tex`
- `manuscript_source/sections/01_intro.tex`
- `manuscript_source/sections/03_model.tex`
- `manuscript_source/sections/05_short_run_dynamic_ineff.tex`
- `manuscript_source/sections/06_hub_alone_vs_refresh.tex`

## Notes

- The manuscript source was copied without LaTeX build byproducts such as `.aux` and `.log` files.
- `reference_docs/Title_Page.docx` and the old cover letter are retained only for reference from the JRS stage.
- The current revision strategy is documented in `notes/02_次投稿先向け改稿TODO_ID_5129307.md`.
- Python / SymPy verification is documented in `verification/README.md` and the latest successful run is recorded in `verification/latest_check_2026-03-25.txt`.
