# Stage 10 v2.1 — Section Map

Date: 2026-09-10 JST

Stage-9 parent / reproducibility authority: `63c5f8c2627e97be6874438ca6af16e1df1c338a`.

Stage-8 theory authority: `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md`.

Stage-7.5A quantifier authority: `reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md`.

The paper is journal-neutral at Stage 10. Journal selection remains deferred to Stage 12.

| Order | File | Function | Frozen inputs used | Stage-10 status |
|---|---|---|---|---|
| 1 | `sections/01_model.tex` | Players, routes, two-sided participation, timing, central smooth regime | model / timing / objective freeze | PASS — no theory change |
| 2 | `sections/02_equilibrium.tex` | Participation fixed point, private pricing, B3 and G, BR-slope identity | equilibrium definitions; matched-price B3 definition | PASS — no theory change |
| 3 | `sections/03_main_results.tex` | T1–T4: B3 first-order complementarity, repricing decomposition, local sign reversal, repaired all-regime witness | analytic theorem register; repaired Stage-4A certificate | PASS — no theory change |
| 4 | `sections/04_welfare.tex` | W1–W3: transfer cancellation, local coordination wedge, one-witness G/B3 welfare comparison | welfare register; generated numerical layer | PASS — no theory change |
| 5 | `sections/05_robustness.tex` | Separates proved local scope from unproved extension directions | Stage-7.5A scope ledger | PASS — rejected 20-draw exercise remains non-authoritative |
| 6 | `sections/06_institutional_empirical.tex` | Maps model primitives to public/private innovation-intermediation analogue class and states empirical predictions | Stage-7 institutional evidence | PASS — prediction/evidence boundary explicit |
| 7 | `sections/07_related_literature.tex` | Positions known components and surviving matched-price public-public reversal | Stage-6 novelty audit; verified bibliography | PASS after Stage-10 literature wording repair |
| 8 | `docs/STAGE_10_EXPOSITION_ARCHITECTURE.md` | Mandatory figure/table architecture gate before Introduction reauthorization | Stage-8 / Stage-9 evidence hierarchy | PASS when architecture validator passes |
| 9 | `sections/08_introduction.tex` | Motivation, model, analytic theorem, repaired witness, welfare, contribution | all preceding sections and architecture gate | PASS — claim ceiling aligned |
| 10 | `sections/09_discussion.tex` | Limitations, novelty boundary, interpretation, policy non-claims | Stage-7.5A scope ledger | PASS — no hidden extension |
| 11 | `sections/10_conclusion.tex` | Restates result at frozen proof ceiling | T1–T4 / W1–W3 | PASS — no overclaim |
| 12 | `sections/appendices.tex` | derivations, local theorem proof, archived rejected witness, repaired all-regime verification, generated parameter table | analytic scripts; Stage-4A audit; generated tables | PASS — authority hierarchy explicit |

## Section-level verification rule

Every section must satisfy all of the following before Stage 10 closes:

1. compile within `paper/main.tex`;
2. preserve the Stage-8 model, timing, objectives, restrictions, theorem quantifiers, welfare objects, and novelty envelope;
3. use generated quantitative objects rather than hand-maintained numbers when a generated object exists;
4. preserve the distinction between the local analytic theorem and the one-vector computational global-equilibrium witness;
5. preserve B3 as a matched-price identification benchmark, not a regulation policy or planner problem;
6. keep the rejected draft vector, its old 20-draw perturbation exercise, and the old exact Krawczyk object outside the active evidentiary chain;
7. resolve citations and cross-references in a clean build.

## Stage-10 manuscript delta

Substantive theory delta: **NONE**.

The only reader-facing Stage-10 repair identified in the section audit is in `sections/07_related_literature.tex`: López and Vives (2019) is described accurately as a simultaneous R&D model with a two-stage strategic-commitment extension, rather than as a direct private-follower-pricing predecessor. Private repricing is instead positioned against the public/private platform and sequential platform-pricing literature. This is a literature-description correction, not a novelty expansion.
