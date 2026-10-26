# Public–Private Innovation Hubs: Hotelling, Entry and Experimentation

Independent exploratory track (created 2026-10-27).
Working branch: research/hotelling-entry-experimentation.

## Research question
Under what conditions should a local government establish a public co-creation hub, rather than rely on private provision? Conditional on establishment, what service orientation and matching scope should public and private hubs choose?

## Scope and relation to existing repository
This track extends a quadratic Hotelling spatial competition benchmark with public-welfare maximization, private profit maximization, relationship-specific experimentation, private entry and local-government establishment decisions.

It is **not** a revision of the canonical public two-sided-platform model. The established research artifacts in public_two_sided_platform_stage8_theory_freeze/ remain unchanged. This branch is to be treated as a separate research hypothesis, not an approved replacement manuscript or a novelty-certified proposition.

## Saved work
- STAGE_4_ENDOGENOUS_LOCATION.md: conditional endogenous-location solution, welfare and entry expressions, revised calibration, and global-equilibrium caveats.
- STAGE_5_PRICE_SUBGAME_AUDIT.md: formally established off-path price-equilibrium multiplicity, a profitable high-price deviation under one continuation, and limitations of prior SPNE/entry claims.
- STAGE_6_CONSTRUCTIVE_SPNE.md: global location best-reply proof and an SPNE construction under an explicitly selected low-price continuation, with v >= 4t and a restricted spillover domain.
- STAGE_7_BACKWARD_INDUCTION_AUDIT.md: market-structure-by-market-structure backward induction, profits, prices, welfare, entry thresholds, realized location corrections, exact rational verification.
- STAGE_8_FULL_MODEL_SETUP_DRAFT.md: full two-sided platform / feasible matching / bilateral Nash-in-Nash / sunk innovation investment / private entry / public establishment setup (PROPOSAL ONLY; not an equilibrium or theory freeze).
- STAGE_9_FULL_MODEL_BACKWARD_INDUCTION_OBSTRUCTION.md: closed-form bilateral bargaining, investment and contact-choice conditions; exact two-sided participant fixed point; mathematically demonstrated multiple equilibria; upstream price/location/entry equations remain set-valued.
- code/stage9_full_model_multiplicity.py: independently runnable pure-Python regression of two distinct equilibria in the same fixed-price matching subgame.
- code/stage4_interior_location.py: executable Stage 4 conditional entry-policy diagnostics.
- code/stage5_price_subgame_audit.py: executable Stage 5 price-subgame counterexample and low-price selection checks.
- code/stage6_constructive_spne.py: executable Stage 6 continuation and policy diagnostics.
- code/stage7_backward_induction_audit.py: exact rational regression checks of the complete restricted government/entry/location/price accounting, plus independent location-deviation grid tests.
- STAGE_1_TO_3.md: chronology, formulas, numerical claims and validation status.
- MODEL_AUDIT.md: exact maintained assumptions and limitations.
- PRIOR_ART_AND_GATES.md: closest literatures and required kill tests.
- code/innovation_hubs_stage3.py: runnable fixed-location Stage 3 calculation supplied in the conversation (Python standard library).

- STAGE_10_TRACTABLE_MODEL_OPTIONS.md: literature-grounded simplification options, bilinear-matching closed-form two-sided demand/price/entry/welfare candidate, and binary partner-scope extension (design proposal, not an SPNE certification).
- STAGE_11_BILINEAR_GLOBAL_SPNE.md: global price Nash proof inclusive of monopoly and duopoly corners, exact private-entry thresholds, and a three-regime public construction rule under explicit sufficient restrictions.
- code/stage11_bilinear_spne_audit.py: exact-rational regression and 25,921 alternative two-sided fee combinations (grid supplements full analytic proof).

## Current verdict
**STAGE 11 FIXED-LOCATION / ONE-LINK BILINEAR TWO-SIDED MODEL: CONDITIONAL SPNE CERTIFIED BY ANALYTIC GLOBAL BEST RESPONSES. ORIGINAL FULL MODEL: STILL OPEN.** Theorem under A1-A4 yields a pure-strategy SPNE of the restricted game with endogenous bilateral match volumes, two joining prices per provider, private entry and public build. It rules out **all** one-provider partial-coverage fee deviations and **all** two-provider corner fee deviations. It does NOT allow firms to relocate when both are active; does not endogenize contact number above one; assumes full registration of each market side when at least one hub exists; has no subsidy/digital counterpart. This is a tractable foundation, not a certified novel economic contribution. Stage 6 is a DIFFERENT single-sided model and should never be mixed algebraically with Stage 11.

## Reproduce
Run from repository root:

    python hotelling_entry_experimentation/code/stage11_bilinear_spne_audit.py

For earlier checkpoints, run their correspondingly numbered scripts.

No pull request is opened for merging into main. Commit future work only to this research branch unless the owner explicitly approves a change in direction.
