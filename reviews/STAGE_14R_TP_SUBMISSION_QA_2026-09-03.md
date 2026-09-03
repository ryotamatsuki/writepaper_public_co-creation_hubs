# Stage 14R-TP — JPET Submission QA

Date: 2026-09-03

Target: Journal of Public Economic Theory (JPET)

Canonical workflow: `ryotamatsuki/research-paper-workflow` v1.1 @ `488e5ab06c207909296a7564eaf9066f7f94319c`

Template: `templates/STAGE_14_SUBMISSION_QA.md`

Checklist: `checklists/SUBMISSION_CHECKLIST.md`

## Executive verdict

**CONDITIONAL PASS — non-substantive metadata completion only.**

No theoretical, mathematical, reproducibility, bibliographic, source-package, or PDF-readability defect requires reopening an earlier research stage. Two factual author-side metadata items remain before Stage 15 submission freeze:

1. enter/confirm the required ORCID in Wiley Authors / Research Exchange and the final title-page record;
2. supply a factual funding declaration (including an explicit no-funding statement if that is the fact).

The redacted title-page audit confirms the current manuscript title, a contact email, and conflict-of-interest language. It deliberately does not expose personal values in logs or public QA files.

## Startup authority

- Startup `main`: `3ee1f0aadd18d4cd5d787ce6bb9bf68460a40d71`
- Startup open PR count: 0
- Production manuscript: `paper/main.tex`
- Theory authority: Stage 8 frozen model
- Current proof authority: L3-0 + L3-1 + merged Stage 13 integration
- Journal positioning authority: Stage 12, with JPET as primary
- Stage 13 authority: PR #41 merged; `THEORY INTEGRATION — GO`

## Current JPET requirements

Current official JPET/Wiley guidance was re-opened on 2026-09-03. Operational conclusions are recorded in `docs/STAGE_14_JPET_REQUIREMENTS_LEDGER.md`.

- Submission platform: Wiley Research Exchange / Wiley Authors
- Article type: Original Article
- Review: single-anonymized
- Initial submission: Free Format
- LaTeX: accepted; upload `.tex` plus compiled review PDF and referenced support files
- Submission fee: none
- Mandatory APC under standard subscription route: none
- Optional OA APC: only if voluntarily electing OA, subject to Wiley terms/waivers
- ORCID: required
- Abstract: current Author Guidelines request a structured abstract
- Keywords: maximum seven
- Data Availability Statement: required for Original Article
- Generative-AI use: disclose transparently when use exceeds spelling/grammar-only editing

Legacy JPET pages that still describe ScholarOne and older formatting are treated as historical/secondary authority. The live Research Exchange form must be rechecked at Stage 15 for field-level requirements.

## Manuscript/package changes allowed and made

Stage 14 made only submission-QA changes:

- converted the existing abstract content into the current structured form without changing the economic claim;
- added seven keywords and candidate JEL codes;
- added submission-facing Data Availability and AIGC disclosure wording;
- verified/corrected bibliographic metadata for the closest cited literature;
- prepared JPET metadata, cover-letter, declarations, reviewer-status, and sanitized source inventory;
- refreshed stale title-page manuscript-title metadata while preserving author-side factual authority;
- added Stage-14 clean-build/package/PDF preflight CI.

No new theorem, proposition, primitive, robustness result, welfare result, novelty claim, strategic choice, player, timing, or model element was added.

## Reproducibility result

Substantive Stage-14 HEAD `3bd8fcb4f5f7f2bb3f4f88582cc3922cc9e5329a` passed:

- reproducibility run `33758486879`: **SUCCESS**;
- Stage-14 submission-QA run `33758486871`: **SUCCESS**.

The Stage-14 job passed all of the following steps:

- clean production reproducibility gate (`make clean && make all`);
- small-beta analytic identities;
- L3-0 symbolic and witness checks;
- exact canonical G/B3 certificate (`verify_l31.py`);
- title-page current-title refresh;
- package and redacted title-page audit;
- PDF structural preflight;
- sanitized source-candidate build;
- QA artifact upload.

The clean production gate regenerates numerical results, tables, figure provenance, bibliography validation, manuscript validation, tests, LaTeX PDF, and verification/manifest outputs.

## Mathematical integrity

**PASS.**

The Stage-8 model remains frozen. The local small-`beta` theorem, canonical numerical comparison, and L3-1 exact certificate regenerate. No stale pre-Stage-13 proof-status wording was found in the submission manuscript. The manuscript continues to distinguish:

- local analytic theorem vs global theorem;
- exact canonical certificate vs general primitive theorem;
- numerical robustness vs proof;
- strategic sign reversal vs welfare ranking.

## Bibliography

**PASS.**

All manuscript citation keys resolve. Closest-reference author/title/year/journal/working-paper metadata were normalized from verified publication/institutional records. No fabricated citation was identified by the audit. Initial reference-style churn is unnecessary under current Free Format rules.

## Submission package

`submission/jpet/` is explicitly a **QA CANDIDATE — NOT FROZEN / NOT SUBMITTED** package containing:

- `README.md`
- `metadata.md`
- `cover_letter.md`
- `declarations.md`
- `reviewer_candidates.md`
- `source_inventory.md`
- `qa_candidate_manuscript.pdf`

The CI also builds an ephemeral sanitized source candidate and excludes internal reviews, journal-strategy documents, hostile-referee notes, GitHub CI metadata, and exploratory Level-3 prose from the journal upload set.

## PDF visual QA

**PASS with one non-material layout warning.**

The 20-page QA candidate was rendered and inspected page by page. No clipped text, clipped equation, unreadable table, broken glyph, unresolved citation/reference, blank page, or overlap was observed. Page 20 contains only the final compact Appendix parameter table and therefore has substantial white space; this is visually suboptimal but not a readability or submission-integrity defect.

See `reviews/STAGE_14R_TP_PDF_VISUAL_QA_2026-09-03.md`.

## Metadata/declaration status

- Current title: PASS
- Contact email presence in title-page authority: PASS
- Conflict-of-interest language: PASS
- Data Availability Statement in manuscript/package: PASS
- Ethics/IRB applicability: N/A — theoretical/computational study
- AIGC disclosure: PASS, subject to final factual author review
- ORCID: **WARNING — missing from redacted title-page pattern audit; must be entered/confirmed before Stage 15**
- Funding: **WARNING — no funding language detected; factual statement required before Stage 15**
- Reviewer suggestions: not mandatory in current Author Guidelines; recheck live portal at Stage 15

## Hard-gate disposition

- A Current-base integrity: PASS
- B Workflow integrity: PASS
- C No theory drift: PASS
- D Mathematical reproducibility: PASS
- E Current journal rules: PASS
- F Official-source conflict resolution: PASS
- G Fee gate: PASS
- H Article-type fit: PASS
- I Claim consistency: PASS
- J Bibliographic integrity: PASS
- K Metadata completeness: **WARNING — ORCID pending factual completion**
- L Disclosure completeness: **WARNING — funding statement pending factual completion**
- M LaTeX package integrity: PASS
- N Package sanitization: PASS
- O PDF visual QA: PASS
- P Editorial readability: PASS
- Q No unresolved substantive issue: PASS
- R No actual submission: PASS

## Stage 15 contract

Do not freeze or submit until the two factual warnings above are resolved and the live Research Exchange form is re-opened to confirm current operational fields. These are administrative/non-substantive completions. They do not require theory reopening or another Stage-13 cycle unless a substantive manuscript change is made.

Final Stage-14 classification:

**CONDITIONAL PASS — non-substantive fixes only.**
