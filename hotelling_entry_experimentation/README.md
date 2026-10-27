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

- STAGE_20_NOVELTY_GATES_AND_LOCAL_THEOREM.md: exact low-matching-investment strength comparative static (public-private overentry wedge increases), proof that Stage19 underentry interval is entirely unrelated to matching k, and full kill tests against 2026 matching-investment efficiency / generalized Hosios results.

- STAGE_23_TWO_REGION_SIMULTANEOUS_OWNERSHIP_SPNE.md: asymmetric urban/rural two-sided pre-match matching investment with THREE simultaneous entry decision makers, 16 ownership/location histories; PP/PG/GP/GG pure entry NE and planner, closed-form PP→PG→GG welfare phase cuts, ownership-count inefficiency; fixed commission / fixed local user-split restriction explicit.

- STAGE_24_STRUCTURAL_AUDIT_ENDOGENIZATION_ROADMAP.md: audited Stage23 endogeneity and imposed assumptions, derived ownership-specific access investment marginal externality, highlighted unverified participation/commission/bargaining, proposed minimal M1–M4 nested extensions and prior-art publication update (Peitz–Sato RAND 2026).

- STAGE_25_CORRECTED_TWO_MUNICIPALITIES_NO_HUB_COURNOT_BASELINE_PRIOR_ART.md: CORRECTED four independent local public/private entrants, two asymmetric urban/rural regions, positive firm output/competition and R&D without hubs, and six precise closest economic-game models (AJ 1988; Cremer et al 1997; Long&Wong 2009; Haaland&Kind 2008; Dur&Staal 2008; Goyal&Joshi 2003). Stage23 one-prefecture model is NOT this game.

- STAGE_28_PARSIMONIOUS_FOUR_PLAYER_ENTRY_BACKWARD_INDUCTION_THEOREMS.md: SOLVED restricted four-local-actor game with independent nonzero Cournot baseline, public-free/private-endogenous one-sided tariffs, consumer sorting, cross-regional producer spillovers; exact four local NE threshold cases, 16-profile-equivalent SPNE conditions, closed-form social planner optimum, private underentry and inverted-U government co-location proofs, exact rational unique NE/planner divergence.

## Current verdict
**STAGE 28 FULLY SOLVED RESTRICTED FOUR-ACTOR BENCHMARK.** Four INDEPENDENT municipality/private operators decide simultaneous local hub establishment; producer Cournot sector exists with zero hubs; G free and P optimal one-sided membership fee (not full two-sided), heterogeneous project users choose G/P/out, region-to-region positive links and knowledge-spillover cost reduction. Analytic all-subgame BRs and all 16 entry profiles yield full necessary/sufficient pure NE conditions via local 2x2 inequalities, and planner global max via 4 macro cases. Exact rational sample: rho=.4 unique NE (urban BOTH, rural GOV), global planner (urban BOTH, rural BOTH), gap .075. Private under-coexistence threshold F in (NΔ/4,3NΔ/8), free-G addition inverted-U in foreign cross value u. Still **NOT** endogenous I/R matching, bilateral Nash bargaining, interconnection choice, endogenous own R&D, location or novel-publication guarantee. Stage23 numerical ownership results cannot substitute for user-corrected Stage28.

## Reproduce
Run from repository root:

    python hotelling_entry_experimentation/code/stage19_zoned_spne_bernstein.py

For earlier checkpoints, run their correspondingly numbered scripts.

No pull request is opened for merging into main. Commit future work only to this research branch unless the owner explicitly approves a change in direction.
