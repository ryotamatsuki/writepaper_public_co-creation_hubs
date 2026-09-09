# Stage 10 v2.1 — Section-by-Section Paper Construction

Date: 2026-09-10 JST

Workflow authority: `ryotamatsuki/research-paper-workflow` v2.1 @ `b27568172d4145a2b98825a5b66a071dbfe25f36`.

Canonical template: `templates/STAGE_10_PAPER_BUILD.md` @ blob `c0d4498d24c3706e73eb0b34ce1ee90496f735a3`.

Stage-9 reproducibility merge / exact Stage-10 parent: `63c5f8c2627e97be6874438ca6af16e1df1c338a`.

Stage-8 theory freeze: `ad927ca783a6123ea4fc6f55f65598ebd6ab583b` / `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md`.

Stage-7.5A quantifier authority: `reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md`.

Working title: **Strategic Interaction among Public Innovation Hubs under Private Repricing**.

Journal target: **NOT SELECTED — deferred to Stage 12**.

## 1. Stage objective

Reauthorize and complete the existing full manuscript against the final v2.1 theory chain, section by section, without reopening theory. The repository already contained a full journal-neutral manuscript inherited from the earlier production chain; Stage 10 therefore performed a clean reconstruction/audit against the final Stage-8 freeze and Stage-9 generated layer rather than treating the old Stage-10 verdict as authority.

## 2. Allowed / prohibited delta

Allowed: exposition, organization, citation description, generated-object inclusion, manuscript architecture, CI/validation.

Prohibited: changes to players, timing, routes, utilities, private objective, public objectives, parameter restrictions, benchmark definitions, theorem quantifiers, equilibrium/globality claims, welfare objects, or novelty envelope.

Theory delta in this Stage: **NONE**.

## 3. Section-by-section result

Detailed section map: `docs/STAGE_10_SECTION_MAP.md`.

### Model — PASS

File: `sections/01_model.tex`.

Frozen inputs: two regional governments, public hubs `H_i`, shared private intermediary `H_T`, project single-homing, partner multihoming, cross-side participation, convex public facilitation cost, private fee-profit objective, four-step timing.

Key verification: no change to route set, support-side equations, objective functions, timing, or central-regime inequalities.

### Equilibrium / benchmarks — PASS

File: `sections/02_equilibrium.tex`.

Frozen inputs: participation fixed point, private FOC/SOC, B3 matched-price fixed-price identification benchmark, full game G, local BR-slope identity.

Key verification: B3 remains `bar p_T = p_T^G` on path and suppresses only the price-policy response; it is not represented as regulation or a planner benchmark.

### Main results — PASS

File: `sections/03_main_results.tex`.

Frozen inputs: T1–T4.

Key claims:

- B3 cross effect is zero at beta zero with a strictly positive first derivative in beta;
- full-game beta-zero cross effect is strictly negative at the stated symmetric regular state;
- sufficiently small beta on continuing regular stationary branches gives `BR_i^{B3 prime} > 0 > BR_i^{G prime}` under strict local SOCs;
- the repaired all-regime result is a one-vector computational global-equilibrium existence witness, not a general theorem.

### Welfare — PASS

File: `sections/04_welfare.tex`.

Frozen inputs: W1–W3.

Key claims:

- private fee cancels as an aggregate transfer;
- the repaired G witness has a positive local coordinated `+x_i` national-welfare derivative;
- G exceeds B3 aggregate welfare at the repaired witness only.

No first-best, global social optimum, optimal subsidy, or welfare-dominance claim is made.

### Robustness / scope — PASS

File: `sections/05_robustness.tex`.

The section correctly separates the local analytic theorem from the repaired one-vector global witness. The rejected old vector, old 20-draw perturbation exercise, arbitrary distributions, arbitrary nonlinear network functions, and heterogeneous regions remain outside certified robustness.

### Institutional / empirical bridge — PASS

File: `sections/06_institutional_empirical.tex`.

Institutional examples are used only to establish an analogue class. The negative private-price response is explicitly a model prediction, not observed causal evidence. Primary-route single-homing and private profit maximization remain reduced-form abstractions.

### Related literature — PASS AFTER BOUNDED EXPOSITION REPAIR

File: `sections/07_related_literature.tex`.

Stage-10 source checking found one wording issue: López and Vives (2019) had been described too directly as a follower-response precedent. The description was corrected to the supported two-stage strategic-commitment interpretation. Public/private fee responses and sequential platform pricing are positioned against Liu et al. (2026), Bontems et al. (2025), and Sánchez-Cartas (2026) without expanding the novelty claim.

Detailed check: `docs/STAGE_10_LITERATURE_CHECK.md`.

### Introduction — PASS

File: `sections/08_introduction.tex`.

The Introduction is reauthorized only after the exposition architecture gate. It states the analytic theorem at local stationary-branch scope, separates the repaired global witness, keeps B3 as matched-price identification, and distinguishes strategic interaction from welfare.

### Discussion — PASS

File: `sections/09_discussion.tex`.

The section records the proof ceiling, old-witness failure, remote-public friction interpretation, single-homing and private-objective abstractions, empirical non-validation, and absence of policy-optimality claims.

### Conclusion — PASS

File: `sections/10_conclusion.tex`.

The conclusion restates the local theorem, the repaired computational existence witness, and local welfare result without upgrading them.

### Appendix — PASS

File: `sections/appendices.tex`.

The appendix carries the small-beta derivation/proof, archived exact stationary-root diagnostic with negative authority statement, repaired all-regime computational verification, and generated baseline-parameter table.

## 4. Mandatory Figure/Table Architecture Gate

Canonical map: `docs/STAGE_10_EXPOSITION_ARCHITECTURE.md`.

Decision: **NO REQUIRED FIGURE**.

Reason: no certified primitive formula for the small-beta endpoint or nontrivial comparative-static path exists. A threshold/path plot would risk implying globality or a regime map beyond the frozen theorem. A four-step timing diagram also fails the cognitive-load test because the timing is already short and linear.

Required reader-facing tables:

1. `tab:strategic` from `generated/tables/strategic_results.tex`;
2. `tab:welfare` from `generated/tables/welfare_comparison.tex`;
3. `tab:parameters` from `generated/tables/canonical_parameters.tex`.

Architecture is enforced by `scripts/validate_stage10_architecture.py` and is part of `make all`.

## 5. Citation / source check

- BibTeX parse/integrity gate remains in `scripts/validate_bibliography.py`.
- manuscript citation resolution remains part of the LaTeX build and build-log validation.
- Stage-10 external source spot-check is recorded in `docs/STAGE_10_LITERATURE_CHECK.md`.
- institutional claims remain sourced to the retained Stage-7 evidence class and keep evidence-level qualifiers.

## 6. Reproducibility / final CI

Stage-9 reproducibility gates remain intact. Stage 10 adds `.github/workflows/stage10-v21-paper-build.yml` requiring Stage-9 ancestry, the clean full production gate, explicit exposition-architecture validation, final Stage-7.5A claim-scope regression, and a no-theory-drift diff guard.

Substantive PR-head tested: `840712b2ae329c1102514c46a2a3e518c1880472`.

Dedicated Stage-10 workflow: run `34417804951`, job `102686459768` — **SUCCESS**.

Verified steps:

- Final Stage-9 ancestry gate — PASS;
- Python / LaTeX environment — PASS;
- `make clean && make all` — PASS;
- explicit Stage-10 architecture gate — PASS;
- final Stage-7.5A claim-scope regression — PASS;
- Stage-8 freeze / repaired Stage-4A / analytic theorem / Stage-7.5A ledger diff guard — PASS.

Independent Stage-7.5A workflow run `34417804940` also completed **SUCCESS** on the same substantive head.

Stage-13 and Stage-14 legacy workflows correctly skipped on the Stage-10 branch.

## 7. Remaining issues

Unresolved theory blocker: **NONE**.

Unresolved proof-scope blocker: **NONE**.

Unresolved exposition-architecture blocker: **NONE**.

Unresolved citation-description blocker identified by Stage 10: **NONE after bounded repair**.

The strongest remaining substantive risk is proposition-level / combination novelty, which belongs to the Stage-11 hostile referee gate rather than Stage-10 construction.

## Final verdict

**FULL DRAFT READY FOR REFEREE GATE**

**STAGE 11 — ROBUSTNESS / REFEREE ATTACK GATE AUTHORIZED.**
