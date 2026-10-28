# Stage 11R-P0b: Contribution focus and near-prior update
Date: 2026-10-28 JST. Theory freeze: STAGE33-2026-10-28-THEORY-FREEZE-v1.
Verdict: P0B EDITORIAL REVISION EXECUTED; NOVELTY STILL CONDITIONAL; STAGE 11 OVERALL HOLD. No theory or equilibrium primitives changed.

## Single claim
The combination of fee-selected founder skill composition, cross-regional matched innovation, differentiated Cournot effects, and four independent provider entry decisions permits a specific change in public/private provision configuration. Three frozen numerical cases are: heterogeneous lambda=0 NE (B,P); heterogeneous lambda=0.5 NE (B,B); pooled/full-population-equal-mean lambda=0.5 NE (B,P). The restricted planner selects (B,B) in all three. This does not establish an all-parameter comparative-static theorem.

## Closest research and claim limits
- Hagiu (2004b): vertical developer heterogeneity and variety/quality welfare reversal already known; no novel generic screening claim.
- Hua, Chen and Zhang (2016): free public/paid private two-tier provision and fiscal policy already known; published full theorem-level verification outstanding.
- Long and Wong (2009): two regional governments with private local price competition already known.
- Peitz and Sato (2026): asymmetric platforms, price equilibrium, entry and different cross-side welfare effects already known.
- Berliant and Fujita (2008, 2011): knowledge heterogeneity/diversity and innovation already known.
- Lee (2026), Journal of Public Economic Theory 28(4), e70128, doi:10.1111/jpet.70128, open full article checked 2026-10-28: optimal tax/subsidy R&D policy can change sign through local network spillovers versus global Cournot rivalry. This absorbs any general novelty claim about sign-changing rivalry/externality effects; this paper instead treats local binary intermediation entry with endogenously selected founder capability moments, not optimal R&D tax/subsidy.
- Häckner and Herzing (2026), JICT 26:16, doi:10.1007/s10842-026-00475-y, publisher open article checked: strategic adoption of cost-reducing technology supplied competitively or by monopoly; comparable applied Cournot publication style, but not a priority certificate for the current paper.
- Cremer, Marchand and Pestieau (1997): familiar local-versus-social public-investment spillover wedge.

A more general predecessor could still nest the joint mechanism. Full theorem-level absorption is NOT certified. See PRIOR_THEOREM_ABSORPTION_REGISTER_2026-10-28.md for prior work and its evidence-maturity limitations.

## Actual manuscript changes
1. main.tex: unified 171-word Abstract; bounded configuration result and adverse robustness.
2. 01_introduction.tex: early question, causal economic channels, exact baseline finding, nearest research and limits.
3. 06_related_literature.tex: closest prior mechanisms grouped; explicit Lee 2026 comparison; table maps known elements and additional integrated result.
4. references.bib: Lee and Häckner-Herzing added with verified metadata.
5. 03_continuation.tex and 05_numerical.tex: P/P baseline global result replaces stale finite-grid-only statement.
6. appendix_a_derivations.tex: analytic fee best-response sufficient conditions, unique global nonnegative price Nash on the three baselines, rational-root enclosure and all-profile entry verification.
7. 07_institutions.tex: concise supported conclusion, no universal physical-hub prescription.
8. stage11r-p0-contribution-manuscript.yml: rerun P0a exact script and canonical calibration, compile LaTeX/BibTeX draft.

## Status matrix
- Baseline P/P global price Nash / four-agent entry / restricted planner: PASS as scoped in prior P0a report; separate referee attack still needed.
- Independent portable novel theorem: NOT CERTIFIED; P1-P3 are largely familiar; P4 is model-specific.
- Numerical configuration reversal: narrow public-cost domain (~1.21% of F_GB) and failure under Beta(2,2) remain.
- Independent differentiated Bertrand replication: NOT DONE in this revision.
- Empirical predictions and measured policy impact: NOT DONE.
- Physical facility, contracts, subsidy superiority: NOT MODELED.
- External human mathematical or journal peer review: NOT DONE.
- Stage 13, 14, 15: NOT CLOSED; no submission.

## Next bounded step
P1: derive which exact binary entry inequalities bind the baseline public-cost interval; verify strict open-region signs without changing primitives or picking more favorable calibration. Then add 2-3 falsifiable observational predictions and policy measurement requirements for JICT. Scope remains model-specific unless genuinely stronger mathematics is proved and all affected stages are re-certified.

All changes remain isolated on a draft PR. No merge, model change or publisher submission.
