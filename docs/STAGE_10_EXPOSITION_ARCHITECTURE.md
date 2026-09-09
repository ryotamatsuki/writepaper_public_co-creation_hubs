# Stage 10 v2.1 — Figure / Table / Exposition Architecture

Date: 2026-09-10 JST

Authority: Stage-8 freeze + Stage-7.5A quantifier ledger + Stage-9 reproducibility baseline.

This map is the mandatory Stage-10 exposition architecture gate required before the Introduction is reauthorized.

| Headline result | Economic object | Primary vehicle | Why this vehicle | Verified source / generator | Required in final paper? |
|---|---|---|---|---|---|
| T1 — first-order B3 complementarity | `M_B3(0)=0` and positive derivative in `beta` | Proposition + displayed equation | Exact sign and local quantifier are more precise than a visual | `analytic_level3/code/derive_small_beta.py`; `verify_symbolic_identities.py` | YES |
| T2 — private repricing channel | chain-rule contribution `P_price` | Displayed decomposition + proposition + concise mechanism prose | Reader needs the decomposition, not a plotted index | symbolic derivation / manuscript equation `repricingdecomposition` | YES |
| T3 — local strategic sign reversal | `BR_i^{B3 prime}>0>BR_i^{G prime}` for sufficiently small positive `beta` on continuing regular branches | Theorem + proof | The result is an existence/continuity theorem; no certified primitive threshold path is available to plot without implying more globality than proved | analytic theorem suite + Stage-7.5A scope audit | YES |
| T4 — repaired global-equilibrium existence witness | one-vector G/B3 public best responses, matched private price, local BR slopes | Numerical illustration + Table `tab:strategic` | Exact certified magnitudes and sign comparison are more informative than a single-point figure | `generated/results/canonical_results.json`; `generated/tables/strategic_results.tex`; independent Stage-4A audit | YES |
| T5 — remote-public route dominance at repaired witness | sufficient condition `kappa_L+tau > v+alpha` | Concise prose + displayed inequality | Binary analytical scope condition; a figure would be redundant | Stage-8 theorem register / Appendix derivation | YES |
| W1 — aggregate fee transfer cancellation | `-p_T n_T^F + Pi_T = 0` | Displayed identity + prose | Exact accounting identity | welfare symbolic verification | YES |
| W2 — local coordination wedge | national marginal derivative at repaired G witness | Displayed decomposition + numerical illustration | Single local derivative decomposition; chart would suggest unsupported comparative statics | `canonical_results.json`; `verify_numerical.py` | YES |
| W3 — repaired-witness G/B3 welfare ranking | aggregate welfare at the two certified witness environments | Table `tab:welfare` | Exact two-case comparison is the object of interest | `generated/tables/welfare_comparison.tex` | YES |
| Robustness / scope ceiling | local theorem versus one-vector global witness; excluded arbitrary distributions/nonlinearities/heterogeneity | Concise prose | The main information is logical scope, not a quantitative surface | Stage-7.5A scope ledger | YES |
| Rejected old witness | profitable finite deviations invalidate former global-equilibrium interpretation | Appendix prose | Negative/regression evidence should not compete visually with active results | Stage-4A historical audit / archived exact diagnostic | YES, appendix only |
| Institutional bridge | public facilitation + priced metropolitan alternative + project-level route substitution | Concise prose | No causal dataset exists; a quantitative figure would overstate empirical content | retained Stage-7 primary-source ledger | YES |
| Novelty boundary | known ingredients versus matched-price public-public sign comparison | Related-literature prose | Conceptual distinction is proposition-level, not numerical | Stage-6 novelty audit + verified bibliography | YES |

## Figure decision

**NO REQUIRED FIGURE.**

This is a deliberate Stage-10 architecture decision, not a missing-output exception. The central analytic result is a local existence theorem without a certified primitive formula for the endpoint `epsilon`; plotting an arbitrary continuation path or normalized sign index would risk communicating a global threshold/regime map that the freeze does not establish. The repaired global result is one certified point, for which a table is the lower-cost and more faithful vehicle. The welfare result is likewise a local derivative plus a two-environment numerical comparison.

A game-timing diagram was considered and rejected as nonessential: the four-step timing is short, linear, and already stated explicitly in the Model section. No diagram materially reduces cognitive load relative to the prose.

A future figure becomes admissible only if an earlier stage is formally reopened and verifies an actual economic object over a nontrivial domain (for example a certified comparative-static path or threshold/regime map). Stage 10 must not manufacture such a domain for exposition.

## Required tables

1. `tab:strategic` — generated from `generated/tables/strategic_results.tex`; required because it reports the repaired witness magnitudes and opposite local slopes.
2. `tab:welfare` — generated from `generated/tables/welfare_comparison.tex`; required because it reports the one-witness G/B3 welfare comparison.
3. `tab:parameters` — generated from `generated/tables/canonical_parameters.tex`; required in the appendix to make the repaired witness reproducible.

`generated/tables/proof_status.tex` remains a generated reproducibility object but is not required in the reader-facing paper because the theorem/witness proof hierarchy is stated directly in the text and appendix.

## Visual-scope rule

No figure or table may communicate globality, uniqueness, robustness, empirical calibration, optimal policy, or primitive-space coverage beyond the Stage-8 freeze. All numerical tables must remain generated from the Stage-9 result layer; no hand-entered quantitative cell is authoritative.
