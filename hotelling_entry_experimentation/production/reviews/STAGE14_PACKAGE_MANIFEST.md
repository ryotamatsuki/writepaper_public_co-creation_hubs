# Stage 14 — JE package manifest and clean-build plan

**Input SHA:** 7308052c0f6299597eba4605ad9b4906966eb8f8 (PR #73 head at commencement). **Status:** preview prepared; submission freeze is gated on author facts and verified clean build.

## Source-of-truth
- Main: `manuscript/main.tex` and the actual `manuscript/sections/*.tex` inputs.
- Supplement: `manuscript/supplement.tex`, `manuscript/appendix/appendix_a_derivations.tex`.
- Bibliography: `manuscript/references.bib`.
- Prior review PDFs: `revision_je_upgrade/je_revised_manuscript.pdf` (37 pages), `revision_je_upgrade/online_supplement.pdf` (8 pages).
- Audit/verification: `reviews/STAGE14_*.*`, previous `JE_UPGRADE_*`, GitHub Actions checks on original commit.
- Draft journal letter: `submission_je/cover_letter.md`.
- Confirmations list: `submission_je/AUTHOR_CONFIRMATION_REQUIRED.md`.
- Rebuild script: `submission_je/prepare_je_package.py`; preview archive will include only necessary submission materials, not internal audit logs/JSON/old draft PDFs.

## Current executable build workflow
Run from the repository root:
`python3 hotelling_entry_experimentation/production/submission_je/prepare_je_package.py`

This uses a temporary isolated directory, recursively flattens active `\\input` sources, runs `latexmk -pdf -bibtex` for main then supplement with `main.aux` available to `xr-hyper`, copies current bibliography, validates compilation logs and PDF pages, and produces `JE_STAGE14_PREVIEW_NEEDS_AUTHOR_CONFIRMATION.zip`. The supplement is packaged as `ESM_1.pdf`. This preview is explicitly NOT a submission-authorised final archive.

**Environment:** Python 3 standard library, `latexmk`, `pdflatex`, `bibtex`, TeX packages used by preambles, `pdfinfo`. GitHub CI builder is declared in `.github/workflows/stage14-je-package.yml` and does not require any new scientific computations.

## Stage 15 entry criteria
1. Author attests title-page email, ORCID, selected affiliation and all declarations/AI use; resolves institution/consent question.
2. Fresh clean CI build successfully generates *both* PDFs from the frozen current source, with no undefined citations/references, important box warnings, or missing files. Review the resulting job/artifact and note exact SHA and actual PDF hashes.
3. Inspect every page of both newly built PDFs; typography, symbols, tables, bibliography, supplement headers and page breaks. Verify no author identity is concealed/inserted contrary to actual portal instructions.
4. Confirm every cited bibliography item/published/accepted status and final JE formatting; verify portal anonymity, any charges and each requested file type.
5. Only after the above, enter Stage 15, record immutable SHA and generate the final submission archive. No PR merge or journal submission is part of this work.

- JE official instructions (accessed 2026-10-29 JST): https://link.springer.com/journal/712/submission-guidelines
- JE publication options (accessed 2026-10-29 JST): https://link.springer.com/journal/712/how-to-publish-with-us
- PR #73 head (read 2026-10-29 JST): 7308052c0f6299597eba4605ad9b4906966eb8f8
- Workflow v1.1 tag: 488e5ab06c207909296a7564eaf9066f7f94319c

