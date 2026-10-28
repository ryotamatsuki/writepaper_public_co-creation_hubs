# Stage 10 — English paper construction and closest-literature differentiation: closeout

**Date:** 2026-10-28 JST
**Theory authority:** `STAGE33-2026-10-28-THEORY-FREEZE-v1`
**Workflow:** `research-paper-workflow v2.9`
**Decision:** **STAGE 10 DRAFT GO → STAGE 11 INDEPENDENT REFEREE/ROBUSTNESS ATTACK**, not journal submission or journal-fit approval.

## Provenance / scope

The entire Stage 10 source is on the stacked draft PR **#52**:
https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/pull/52
Its base is Stage 9 draft branch `stage9/stage33-reproducibility-manuscript-setup-20261028`, not main. Stage 9 PR #51 separately remains a draft and the separate original public two-sided-platform theory is untouched.

## Completed artifacts

1. 14-page modular English manuscript: updated abstract, Introduction, Model, Backward Induction, Propositions P1–P3/C1 with proof bridges, local-versus-national welfare, P4 calibrated 16-profile example, adverse robustness limits, detailed related-literature comparison matrix, international policy context, conclusion, supplementary derivations and numerical/formal certificate distinction.
2. Source-level bibliography including official UK/EDIH/NIST references and directly relevant economics predecessors.
3. `production/STAGE10_PRIOR_ART_DIFFERENTIATION_AUDIT.md`: twelve comparator rows with already-known mechanisms, the exact additional joint object, evidence-access limits, overlap risks, forbidden novelty claims.
4. `production/EXPOSITION_ARCHITECTURE.md`: theory-first reader-arrival and main-text necessity gate; no invented exhibit quota.
5. `production/STAGE10_REVIEWER_VERIFIABILITY_MAP.md`: every P1–P4/L/C proof bridge, pricing kink, fixed-point limitation, and reviewer blockers.
6. `.github/workflows/stage33-stage10-manuscript.yml`: independently runs Python baseline and LaTeX/BibTeX manuscript build; uploads reproducible CSV and built PDF.

## Verification of final source code and PDF

- Stage 10 GitHub Actions run **40624011620** at source SHA `e05c23636429cf40450f6ac004e71cf72125c564`: completed **SUCCESS** for both `latex-and-citations` and `numerical-regression`.
- Build produced a **14-page A4 PDF**. PDF page-image inspection identified and corrected one right-edge overrun caused by an unbreakable code path; after correction, independent text-box preflight found **no off-page text lines**. The final PDF was re-rendered and visually checked at the formerly affected page.
- Stage 9/10 source numerical code enumerated **48 profile records** across three scenarios (16 each): skill lambda 0 -> equilibrium (B,P), skill lambda 0.5 -> (B,B), pooled lambda 0.5 -> (B,P); restricted planner (B,B) all three. This is precisely the frozen baseline, not new calibration.
- The entire 18-theorem Lean proof-critical certificate is inherited unchanged. No full-game SPNE proof or universal global P/P price certificate is asserted.

## The substantive novelty sentence (allowed)

> The paper integrates endogenous privately priced founder sorting by both participation mass and technical ability, cross-jurisdiction complementary venture commercialisation with differentiated-product incumbent/consumer effects, and four independent local public/private provision incentives, yielding conditional analytical thresholds and a restricted, reproducible numerical ownership-allocation illustration.

## Explicitly rejected priority claims

Not first for sorting, public-paid competition, public/private coexistence, complementary R&D, downstream business-stealing, local public-good underprovision, or four geographically distributed government/private roles. No universal ownership switch, physical-facility necessity, first-best planner, or unconditional source-to-theorem non-absorption. Existing Beta(2,2) and Bertrand negative portability remains.

## Next-stage referee gates

1. **Closest-paper absorption**: verify precise theorem-level text (especially Hua--Chen--Zhang 2016, Kim 2024, Peitz--Sato 2026, Berliant--Fujita, Cremer--Marchand--Pestieau, Long--Wong) rather than relying on abstract-level differences.
2. **Economic value**: demand a clear economically material addition beyond simply composing known sorting plus regional spillovers. Test the necessary ablations and robustness `q` vs `m` and the nontrivial business-stealing channel.
3. **Equilibrium audit**: independent verification of every off-path fee regime and the 16 first-stage incentive comparisons, explicitly preserving price selection/multiplicity.
4. **Numerical portability**: reproduce the frozen Beta/Bertrand adverse tests and check narrow welfare margins; do not label them robust.
5. **Literature/factual audit**: validate exact BibTeX records and full theorem statements, and separate 2026 UK/EU policy claims from empirical identification.
6. **Submission policy**: author declarations, target journal and journal requirements are Stage 12–15 items; nothing has been submitted.

Stage 8 is unchanged and frozen. Stage 10 manuscript is a legitimate full-paper development draft with bounded, explicitly differentiating claims. No theory-fix or proof expansion was smuggled into this prose stage.
