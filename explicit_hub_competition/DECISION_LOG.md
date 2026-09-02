# Stage 3M-R Decision Log

## Input provenance

- Repository: `ryotamatsuki/writepaper_public_co-creation_hubs`
- Input PR: #17
- Input state at Stage start: open, unmerged, mergeable
- Input head: `2af59774b3b9bd7fcc383989f7d08ddd342e6252`
- New branch: `co-creation-theory/stage3m-r-explicit-hub-competition`
- Workflow SHA: `07466bcb1a6d3bc654b52945f21b034b38e45281`

## Rollback

Decision: **limited rollback to Stage 3M only**.

Reason: the previous model had government strategic interaction but no endogenous choice between simultaneously available hubs. Explicit competition for intermediation demand was essential to the original research motivation.

## Selections

- Competitive object: **C4 — joint beneficiary-firm/project representation**.
- Hub differentiation: **D1 — partner-set differentiation**.
- Allocation: **A1 — expected productive-match value**.
- Project population: **Option L — region-specific beneficiary firms may use either hub**.
- Incidence: **I3 — beneficiary-firm and complementary-partner components accrue to their respective home regions**.

## Canonical competition mechanism

Overlapping project needs for which both partner pools can generate productive matches are endogenously allocated between hubs.

## Canonical cooperation mechanism

An unsuccessful first-hub search may be referred to non-overlapping capabilities in the rival hub's partner pool.

## Same-structure result

Competition and cooperation are both generated from the same underlying pair of partner pools `S_1,S_2`:

- overlapping feasible capability coverage -> competition;
- non-overlapping capabilities -> referral gains.

## Rejected architecture choices

- C1 alone: too coarse relative to project-level intermediation.
- C2 alone: loses direct mapping to downstream firm identity.
- C3: adds unnecessary second endogenous side.
- D2 as independent quality parameter: unnecessary if quality follows from pool composition.
- D4 Hotelling: ad hoc for baseline.
- A2 separate intermediation costs: unnecessary at architecture stage.
- A3 logit: rejected as smoothing device.
- A4 random splitting: tie breaker only.
- common unlocated project pool: rejected because it risks collapsing hub choice into location competition.

## Prior-art warning

The architecture remains under strong threat from competing-matchmaker theory (Damiano–Li; Caillaud–Jullien; Gal-Or), productive-collaboration matching (Banal-Estañol et al.), and regional public-input/location competition (Walz and related work).

The structural object that survives for Stage 4M-R hard kill is **split regional incidence of a collaboration link without beneficiary-firm relocation**.

## Stage verdict

**PASS — EXPLICIT HUB COMPETITION HAS A COHERENT MINIMAL MATCHING ARCHITECTURE.**

This is not a novelty verdict.

## Next route

Proceed to **Stage 4M-R — Explicit Hub Competition Minimal Model Hard Kill** under `STAGE_4M_R_CONTRACT.md`.