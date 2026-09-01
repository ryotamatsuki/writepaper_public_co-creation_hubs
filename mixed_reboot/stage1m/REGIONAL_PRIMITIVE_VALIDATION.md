# Regional Primitive Validation

## Purpose
`prefecture-year` remains the regional-state unit. It is not an architecture outcome. Stage 1R-M asks whether a decision unit can be linked to predetermined regional conditions without recreating the failed regional union.

## Standard controls
The following are **VALIDATED / HIGH FEASIBILITY** at prefecture-year scale from standard Japanese statistics or NISTEP regional indicators:
- establishment / employment density;
- industrial structure;
- university/research capacity;
- R&D and patent activity;
- broad accessibility.

## Local-network candidate
### Route A — science–industry local linkage density
Use the **absolute count/value of university–industry collaborations with firms in the same prefecture**, normalized by a predetermined local research base (e.g. university/researcher count). NISTEP's regional science and technology indicators aggregate MEXT university–industry collaboration microdata by prefecture and distinguish same-prefecture from other-prefecture firms.

Verdict: **VALIDATED WITH CHANNEL BOUNDARY**.

Interpretation authorized for Stage 2R: `local science–industry linkage density` (or its inverse as a *science–industry fragmentation proxy*).

Interpretation NOT authorized: general ecosystem/network fragmentation.

### Route B — local co-patent network density
Conceptually stronger for general network structure, but a reproducible public prefecture-year extraction was not completed here.

Verdict: SECONDARY / NOT YET VALIDATED.

### Route C — official project/repeat-tie networks
Potentially useful for selected regions; not harmonized nationally.

Verdict: SECONDARY.

## External-connectivity candidate
### Route A — external science–industry linkage intensity
Use the **absolute count/value of university–industry collaborations with firms outside the university's prefecture**, normalized by the same predetermined research base.

This corrects the failed Stage 1R design: local and external measures are not entered as shares that mechanically sum to one. They are separate activity levels and may both be high or both be low.

Verdict: **VALIDATED WITH CHANNEL BOUNDARY**.

### Route B — outside-prefecture co-patent/project links
Conceptually broader, but not yet a validated harmonized panel.

Verdict: SECONDARY.

## Independence test
The two validated primary proxies:
- arise from distinct same-prefecture vs other-prefecture collaboration counts;
- are normalized separately rather than converted to complementary shares;
- therefore are not mechanically `L + X = 1`;
- nevertheless share the same university–industry data-generating channel and must not be treated as a complete description of the regional innovation network.

## Temporal rule
For a decision episode at `t`, use the nearest defensible pre-decision NISTEP/MEXT edition or multi-year average ending before the choice when available. For older launch episodes, use the historically corresponding regional-indicator edition; if only post-choice data are available the variable is flagged and excluded from causal/stylized-fact interpretation.

## Stage 1R-M conclusion
The previous blocker is narrowed rather than hidden. There is now one validated **channel-specific local linkage measure** and one separately measured **channel-specific external linkage measure**. Stage 2R is authorized to test them as constrained proxies, not as generic local fragmentation/external-connectivity primitives. Broader co-patent/project-network measures remain a robustness objective.