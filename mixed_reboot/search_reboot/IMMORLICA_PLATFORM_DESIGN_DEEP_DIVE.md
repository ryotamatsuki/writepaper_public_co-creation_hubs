# Immorlica–Lucier–Manshadi–Wei Deep Dive

Reference: Nicole Immorlica, Brendan Lucier, Vahideh Manshadi & Alexander Wei, *Designing Approximately Optimal Search on Matching Platforms*, Management Science 69(8):4609-4626. Published online 2022; issue 2023. DOI 10.1287/mnsc.2022.4601.

## Why this paper was not adequately captured at Stage 0R-T
Stage 0R-T concentrated on search intermediaries, local/international search, directed-search segmentation and public intermediation. The decisive neighboring field is broader: **matching-platform search design**. This paper does not need a local/external label because it gives the designer a strictly richer control space.

## Players and dynamics
The market is two-sided and dynamic. Each side contains finitely many types. Agents arrive exogenously and may exit unmatched. For any pair of opposite-side types, match value contains a type-specific distribution and is learned when a meeting occurs. Agents strategically accept or reject meetings based on continuation value.

## Designer's control
The platform sets meeting rates between each pair of types. These rates are subject to feasibility/capacity/flow constraints. Hence the designer controls a bipartite meeting graph together with edge intensities.

Economically:
- edge absent / rate zero = type pair is not made operationally accessible;
- edge present / positive rate = pair is accessible;
- larger positive rate = greater meeting intensity/productivity for that pair.

## Objective
For any feasible meeting-rate system, the authors characterize a unique stationary equilibrium. The platform's optimal directed-search problem chooses meeting rates to maximize **equilibrium social welfare**. Congestion and cannibalization make the optimal design nontrivial, and optimal/approximate designs may deliberately restrict choice.

## Exact mapping of the Stage 0R-S candidate
Partition types by geography or policy relevance:
- `L_B, L_S`: local seeker/provider types;
- `E_B, E_S`: external seeker/provider types.

Then:
- `H`: increase meeting rates on selected local-local edges;
- `P`: allow selected local-external edges that were previously zero, or increase those rates;
- `HP`: do both;
- `0`: keep baseline matrix.

Thus the candidate four-way menu is a constrained subset of the meeting-rate matrix design problem.

## Why a separate 'access-cost' parameter does not preserve novelty
One could insist that P changes a discrete access cost whereas Immorlica changes a meeting rate. But if the economic outcome of P is to make a previously unavailable type pair meetable, zero-versus-positive meeting intensity is the operational access margin. Teh (2024), Banerjee et al. (2017), and Mekonnen (2026) independently establish accessibility/visibility/segmentation as equivalent or closely related platform controls. Re-parameterizing an edge activation as an access cost would not create a new strategic interaction absent an additional primitive.

## Classification
**EXACT PRIOR ART — GENERAL-MODEL SUPERSET** for the stripped designer problem.

This classification means the exact discrete institutional labels are absent, but the proposed endogenous margins and welfare design problem are nested in a more general formal model.

## Binding consequence
Do not proceed to a bespoke Search minimal model. Any future theory contribution must introduce an evidence-supported primitive that is not reducible to pairwise meeting-rate/visibility design.
