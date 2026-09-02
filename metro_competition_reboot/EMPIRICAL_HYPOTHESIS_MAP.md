# Empirical Hypothesis Map

The theory is useful only if the surviving objects map to observable moderators and counterfactual outcomes.

## Core variables

| Theory | Observable counterpart |
|---|---|
| `x_i` matching-evaluation capacity | coordinator/specialist FTE, specialist coverage, active partner assessments per case, matching-program expenditure |
| project–ecosystem fit | industry/technology alignment, expertise similarity, capability tags |
| metro access | membership/usage fee, travel/time, prior metro network, existing national intermediary relationship |
| local overlap | Jaccard/member overlap, capability-vector similarity, industry-specialization similarity between hubs |
| creation | first new collaboration where no alternative intermediary relationship existed |
| displacement | switching/rerouting from an existing intermediary |
| match quality | follow-on project, contract, patent, sales/productivity, repeat collaboration |
| firm heterogeneity | size, age, innovation capability, prior network degree |

## H1 — Local-fit heterogeneity

The effect of a local hub on new productive collaborations is larger for projects whose needs align with capabilities present in the local partner ecosystem.

Empirical interaction: `HubExposure × LocalFit`.

## H2 — Metropolitan-access attenuation

Better pre-existing metropolitan access lowers the extensive-margin additionality of local public matching capacity, while raising the fraction of observed local-hub users who are switchers rather than newly matched projects.

Interaction: `HubExposure × MetroAccess`; outcomes must separately measure new collaboration and intermediary switching.

## H3 — Rival-overlap duplication

Conditional on own capability, higher overlap with a neighboring public hub reduces the incremental number of newly created collaborations from additional local matching capacity.

Interaction: `Delta x_i × RivalOverlap`.

## H4 — Triple condition for high additionality

The largest creation effect should occur where local fit is high, metropolitan access is weak/costly, and rival local overlap is low.

Conceptual regression:

`Outcome = beta0 + beta1 Hub + beta2 Hub×LocalFit + beta3 Hub×MetroAccess + beta4 Hub×RivalOverlap + beta5 Hub×LocalFit×MetroAccess×RivalOverlap + ...`.

## H5 — Prior-network heterogeneity

Projects with low prior network degree should experience a larger `0 -> productive match` response to local screening capacity than already well-connected firms; well-connected firms should exhibit more displacement/intensive-margin use.

## H6 — Private price response

A credible increase in nearby public matching capacity changes the price/fee or effective offer of metropolitan private intermediaries. The sign is not frozen; Stage 1 will characterize conditions. This is an indirect competitive channel, not a direct treatment effect.

## H7 — Productive outcome versus venue substitution

If observed use increases but follow-on contracts/patents/sales do not, the event is consistent with displacement rather than productive additionality. A valid evaluation should therefore combine utilization with downstream match-quality outcomes.

## H8 — Regional incidence is separate

A local hub may increase the share of collaborations with local partners without increasing total collaboration value. Partner-location effects should be estimated separately from productive outcomes to avoid mislabeling redistribution as additionality.

Recent Kohsetsushi evidence on heterogeneous effects by firm size, industry, spillover channel, agglomeration and knowledge base, and 2026 EDIH evidence on selective take-up, make these moderators empirically plausible rather than purely theoretical.