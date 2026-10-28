# Stage 9 — Production/reproducibility setup decision

Date: 2026-10-28 (JST)
Frozen theory identifier: STAGE33-2026-10-28-THEORY-FREEZE-v1
Decision: SETUP IMPLEMENTED ON REVIEW BRANCH; NOT A STAGE 10 MANUSCRIPT GO
Repo: ryotamatsuki/writepaper_public_co-creation_hubs
Dedicated PR base: research/hotelling-entry-experimentation, never main.

## Immutable anchors
- hotelling_entry_experimentation/STAGE8_CANONICAL_THEORY_FREEZE_V1.md (Git blob bafd537937a5abb40c57ac3af3b59df78ce5f702).
- hotelling_entry_experimentation/STAGE8_THEORY_FREEZE_MANIFEST.md.
- formal/Stage33Formal.lean; Lean 4 v4.32.1 + pinned Mathlib (SHA 520045ab14e26149ee970e2e617ca04b09bde5d6).
- theorem_certificates/STAGE33_LEAN_FORMAL_VERIFICATION_CERTIFICATE.md.
- hotelling_entry_experimentation/INTRO_REAL_WORLD_HANDOFF_STAGE9.md.
- hotelling_entry_experimentation/INTRO_REAL_WORLD_PRIMARY_SOURCE_LEDGER.md.
No frozen file has been modified on this branch.

## Reproduction protocol and result-to-claim map
- P1: interior differentiated Cournot incumbents and entrant; M_i, Pi, R0, S_i enter welfare as in frozen source. Current regression checks baseline levels, not derivative proof.
- P2: posted fees under fixed opposite-region participation; compare B-only private fee, participant mass and first skill moment at baseline. Numeric value and local cutoff are checked; analytic comparative statics remain governed by frozen P2.
- P3: compare (B,B) with (B,P) for region B, holding the other region fully covered. Check heterogeneous and pooled local government welfare differences; national wedge is computed from both regional W vectors, not silently asserted as a first-best.
- P4: enumerate sixteen provision configurations; best-reply binary flips for both G_i and P_i; restricted-planner maximum. Baseline equilibria: lambda=0 -> (B,P), lambda=0.5 skill -> (B,B), lambda=0.5 pooled -> (B,P). Planner (B,B) in all three. The local B public-addition root approximately 0.3511840318401. This is a conditional numerical witness with price-continuation selection, NOT a universal theorem.
- L1--L3/C1: frozen analytic results; no additional global numerical certificate created by Stage9.
- Negative portability: Beta(2,2) and differentiated Bertrand remain recorded failures for the same switch. Not reverified by this baseline-only script.

## Source and institutional scope
UK parliamentary report 13 March 2026 + government response published 16 June 2026:
actual disagreement over Catapult regional market overlap and hybrid finance; they do
NOT represent two independently entering providers G and P.
EU EDIH overview + operational FAQ D.12/E.1/E.2: regional/cross-region services,
free or discounted tariffs and a list-price aid-accounting instrument. They do
NOT establish real-world crowd-out effects, all-free hubs, or pure municipal owners.
NIST Montana MEP start-up case: optional concrete matching vignette, no causal
welfare inference. References are provided in manuscript/references.bib.

## Explicit Stage9 acceptance boundaries
PASS in the previous conversation's independent Python reconstruction: benchmark
equilibria, restricted planner, 2.407438368575148 full-coverage welfare and
B-private fee 0.6410074280605 were reproduced.
The committed standalone script and GitHub CI results must be inspected before
declaring the repo-level regression gate PASS.
LaTeX full-paper completion, journal compliance, full-proof insertion, robustness,
closest-prior-theorem article-by-article verification, final Introduction and
third-party reproducibility are Stage10--14 tasks.
Any CI, manuscript, bibliography or numerical disagreement must be documented as
an OPEN blocker rather than retroactively changing the frozen result.
