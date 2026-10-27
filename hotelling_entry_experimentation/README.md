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
- STAGE_12_ENDOGENOUS_LOCATION_REOPEN.md: analytical proof that Stage 11's fixed profile is NOT a location-game equilibrium; general interior two-sided location/fee equations; symmetric-side cubic candidate; global location SPNE unresolved.
- STAGE_13_BACKWARD_INDUCTION_LOCATION_AUDIT.md: independently audited backward-induction interior cubic, private profit, welfare and conditional entry thresholds; rigorous demonstration of widely differing off-path Nash fee continuations, hence **NO globally certified location SPNE yet**.
- code/stage13_endogenous_location_price_multiplicity.py: pure-Python exact inequality-based price Nash counterexample and stationary-root numerical reproduction.
- code/stage12_location_reopen.py: exact rational profitable relocation counterexample and analytic-cubic consistency diagnostic.
- code/stage11_bilinear_spne_audit.py: exact-rational regression and 25,921 alternative two-sided fee combinations (grid supplements full analytic proof).

- STAGE_14_PRICE_SUBGAME_REGIME_AUDIT.md: exact 9-case FULL-COVERAGE membership correspondence, interval-based finite-v opt-out fixed-point map, all-private price-NE continuum and rational PURE-PRICE NONEXISTENCE witness in a full-coverage off-path location (not a finite-v SPNE nonexistence theorem).
- code/stage14_price_subgame_regimes.py: exact rational branch enumeration, global fee-BR inequalities, full-coverage nonexistence certificate, and finite-v interval demand-map diagnostics.

- STAGE_15_PURE_SPNE_NONEXISTENCE_OPEN_SET.md: strict open-set theorem excluding pure fee Nash at feasible off-path locations and therefore pure SPNE in the MANDATORY-FULL-COVERAGE location game (NOT in original finite-v opt-out game); comprehensive close-prior-art audit.
- code/stage15_nonexistence_open_set_check.py: independently executed exact-rational strict-condition verification.

- STAGE_16_PUBLIC_BACKSTOP_PROPOSITION_SCREEN.md: mechanism-first SCREEN (not a finished Stage 12 SPNE): credible idle public fallback changes private investment and can crowd IN private matching-platform entry; four-zone exact public-installation cutoff, strict both-alone-unviable joint-viability example, and crucial direct 2005 prior-art absorption warning.

- STAGE_17_PRIVATE_INCUMBENTS_PUBLIC_UNDERSUPPLY_AUDIT.md: with existing private hubs, distinguishes number-of-hubs vs missing cross-hub matching FUNCTION, simple underprovision window, private procurement no-go, and closest 2026 OECD and platform competition priors.

- STAGE_19_SIMULTANEOUS_ENTRY_PREMATCH_INVESTMENT_ZONED_SPNE.md: fully constructed **ZONED** simultaneous public/private entry SPNE, global 2-sided fee Nash, analytic pre-match investment, endogenous locations with exact Bernstein global-BR and contraction certificates, payoffs in all 4 entry structures, over/under/efficient entry examples.
- code/stage19_zoned_spne_bernstein.py: symbolic exact rational Bernstein verification, numerical location root and stage-1 entry/second-best welfare examples.

## Current verdict
**STAGE 19: COMPLETE CONSTRUCTIVE SPNE IN A NEW, EXPLICITLY ZONED AND NORMALIZED SIMULTANEOUS-ENTRY GAME.** Public and private entry are simultaneous; neither preexisting hub is assumed; each selects its location endogenously within its own restricted area [0,.4] or [.6,1]; both choose endogenous two-sided fees; matching and innovator effort BEFORE realized matching are endogenous; all entry structures and payoffs are solved. A globally valid two-sided price continuation exists at ALL zoned locations; unique zoned location Nash is certified with exact rational Bernstein inequalities and best-response contraction; off-path consumer enrollment is selected using a compact fixed-point correspondence. **This is NOT the unrestricted spatial model, NOT a general comparative-static M/R theorem, and NOT a novelty-certified publication result.** Stage 15 unrestricted-location pure-SPNE impossibility pertains to a different mandatory-full-coverage game. Both under-entry and over-entry emerge by varying fixed costs in Stage19, but those are classical mechanisms until shown otherwise.

## Reproduce
Run from repository root:

    python hotelling_entry_experimentation/code/stage19_zoned_spne_bernstein.py

For earlier checkpoints, run their correspondingly numbered scripts.

No pull request is opened for merging into main. Commit future work only to this research branch unless the owner explicitly approves a change in direction.
