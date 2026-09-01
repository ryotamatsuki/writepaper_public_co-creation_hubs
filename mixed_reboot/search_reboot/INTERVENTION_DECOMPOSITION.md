# Intervention Decomposition

Stage 0R-S decomposes the surviving candidate into two formal margins rather than policy labels.

## T — Within-market matching productivity
Candidate H was intended to increase the rate/efficiency with which already-accessible local seeker-provider pairs obtain meetings.

This margin is already standard in several literatures:
- meeting/matching productivity parameters in matching functions;
- public-employment-agency matching effectiveness/contact-arrival rates;
- search targetability/quality;
- pair-specific meeting rates chosen by matching platforms.

Most important boundary: Immorlica et al. let the platform select `lambda_{ij}`-type meeting rates for each pair of types. Raising an already-positive local rate is an H-like intensive intervention.

## X — Market accessibility / segmentation
Candidate P was intended to make a previously segmented external pool reachable.

This margin is already formalized as:
- zero vs positive visibility/meeting links in visibility or meeting graphs;
- market segmentation into submarkets;
- meeting technology/accessibility choices;
- directed-discovery menus;
- foreign-market/trade search costs.

Most important boundaries:
- Teh-Wang-Watanabe: platform meeting technology maps to accessibility;
- Mekonnen: planner chooses segmentation and buyer allocation/tightness;
- Banerjee et al.: platform controls buyer-seller visibility graph.

## The joint T_L versus X_E problem
The Stage 0R-T hope was that the public actor would face an unresolved trade-off between:
1. improving matching inside incumbent market L; and
2. opening market E.

That distinction does not survive general matching-platform design. If L and E are represented as type subsets, a platform that controls all type-pair meeting rates controls both margins jointly:
- `T_L`: increase positive L-L rates;
- `X_E`: turn selected L-E/E-L rates from zero to positive and/or increase them.

Thus the discrete H/P control set is nested in a continuous/high-dimensional meeting-rate matrix.

## Why this is stronger than component overlap
The relevant designer already anticipates equilibrium responses and optimizes welfare. Congestion/cannibalization also enter the general problem. The proposed endogenous interaction `P -> allocation/tightness -> marginal value of H` therefore does not supply a new strategic feedback by itself.

## Conclusion
`T_L` and `X_E` are useful policy diagnostics but are not a distinct theoretical control pair under modern search-platform design.

**Verdict: SEARCH ABSORBED.**
