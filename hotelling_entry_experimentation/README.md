# Public–Private Innovation Hubs: Hotelling, Entry and Experimentation

Independent exploratory track (created 2026-10-27).
Working branch: research/hotelling-entry-experimentation.

## Research question
Under what conditions should a local government establish a public co-creation hub, rather than rely on private provision? Conditional on establishment, what service orientation and matching scope should public and private hubs choose?

## Scope and relation to existing repository
This track extends a quadratic Hotelling spatial competition benchmark with public-welfare maximization, private profit maximization, relationship-specific experimentation, private entry and local-government establishment decisions.

It is **not** a revision of the canonical public two-sided-platform model. The established research artifacts in public_two_sided_platform_stage8_theory_freeze/ remain unchanged. This branch is to be treated as a separate research hypothesis, not an approved replacement manuscript or a novelty-certified proposition.

## Saved work
- STAGE_1_TO_3.md: chronology, formulas, numerical claims and validation status.
- MODEL_AUDIT.md: exact maintained assumptions and limitations.
- PRIOR_ART_AND_GATES.md: closest literatures and required kill tests.
- code/innovation_hubs_stage3.py: runnable fixed-location Stage 3 calculation supplied in the conversation (Python standard library).

## Current verdict
**EXPLORATORY / NOT SUBMISSION READY.** A fixed-location entry-policy benchmark has runnable calculations. Stage 1's fully endogenous hub-location numerical claims are not preserved in runnable code. The proposed endogenous-location, endogenous-experimentation, private-entry subgame has not been globally solved. Novelty and non-absorption by prior results are unresolved.

## Reproduce
Run from repository root:

    python hotelling_entry_experimentation/code/innovation_hubs_stage3.py

No pull request is opened for merging into main. Commit future work only to this research branch unless the owner explicitly approves a change in direction.
