# Stage 9 production scaffold: public and private innovation intermediation

Date: 2026-10-28 JST. Workflow: research-paper-workflow v2.9.
Theory freeze: STAGE33-2026-10-28-THEORY-FREEZE-v1.
Canonical authority: ../STAGE8_CANONICAL_THEORY_FREEZE_V1.md.
Source branch: research/hotelling-entry-experimentation.
This scaffold is independent of the unrelated main-branch two-sided-platform theory.

## How to reproduce

From this directory with Python 3.12+ (standard library only):
1. make check
2. make data

The first command evaluates all 16 entry profiles in each of three scenarios:
- lambda=0, heterogeneous;
- lambda=0.5, heterogeneous;
- lambda=0.5, pooled success advantage (same fully covered mean).
It checks the Nash candidates against all four binary unilateral entry flips; evaluates
the restricted sixteen-profile planner maximum; checks private price, endogenous
participants, consumer-and-incumbent welfare, the local threshold at lambda=0.5,
and the local government's indifference root. The second command writes 48
deterministically ordered CSV rows to output/stage33_profiles.csv.

The P/P pure-price continuation is computed by iterating best replies, with
independently recomputed participation and finite-grid no-profitable-deviation
diagnostics. These are baseline numerical checks, NOT an unrestricted price
equilibrium existence, uniqueness or formal certificate. The complete Stage 4A
sufficient-condition analysis remains in the frozen external research records.
Do not infer a global comparative static from the checked calibration.

## Manuscript setup

- manuscript/main.tex: compileable modular article scaffold
- manuscript/sections/01_introduction_seed.tex: policy motivation using official UK/EU material
- manuscript/sections/02_model_summary.tex: expository, frozen-scope primitives
- manuscript/sections/03_results_map.tex: bounded theorem register and numerical evidence
- manuscript/references.bib: official institutional sources used for the introduction

Run make paper where latexmk, LaTeX, natbib, amsmath, lmodern, geometry and
related TeX dependencies are installed. These are drafting templates, not a
Stage 10 completed paper, refereed proof exposition or submission-ready manuscript.

## Key limitations and change control

- Public G means free, publicly financed basic service access, not necessarily
  a physically built or publicly operated building.
- The sixteen-profile social comparison is a restricted-instrument optimum,
  not first-best policy design.
- Lean's 18 proof-critical checks do not certify all 16 SPNE or unrestricted
  P/P tariffs; see ../STAGE8_THEORY_FREEZE_MANIFEST.md and
  ../../theorem_certificates/STAGE33_LEAN_FORMAL_VERIFICATION_CERTIFICATE.md.
- Beta(2,2) types and differentiated Bertrand negative portability tests remain
  part of the freeze and are not neutralised by this reproducibility example.
- No payoff, timing, benchmark, quantifier, theorem or participation domain is
  altered here. Any scientific modification requires the recorded change-control
  process and re-opening relevant theory gates.
- Sources are illustrations of institutional arrangements, not observations
  establishing sorting, counterfactual welfare, or the four-player equilibrium.

Stage 10 order: model -> continuation -> propositions -> welfare -> sensitivity
and adverse tests -> institutional bridge -> literature -> figures -> introduction
-> discussion -> conclusion. The present introduction is provisional.

See STAGE9_SETUP_DECISION.md for the result-to-source traceability and pending QA.
