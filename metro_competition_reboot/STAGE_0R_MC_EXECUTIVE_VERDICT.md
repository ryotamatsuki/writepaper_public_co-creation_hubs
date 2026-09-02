# Stage 0R-MC Executive Verdict

Date: 2026-09-02

Starting main SHA: `d013afe9c14792abc7491ad0de2d8d416e846cb8`
Workflow SHA: `07466bcb1a6d3bc654b52945f21b034b38e45281`
PR #21 status: **NEGATIVE BENCHMARK — NOT NEW BASE**.

## Verdict

\[
\boxed{\textbf{CONDITIONAL GO — ONE ARCHITECTURAL ISSUE REMAINS}}
\]

The generic ingredients are heavily absorbed by prior art. Damiano & Li (2008) absorb heterogeneous endogenous choice among competing matchmakers; Takahashi (2004) absorbs strategic public-facility investment; Kim (2024) absorbs public/private two-sided platform competition with investment incentives; Liu, Reshidi & Rivadeneyra (2026) explicitly absorb the extensive-margin participation versus switching/displacement decomposition in public/private platform competition; Calcagnini et al. (2019) and related work absorb endogenous search intensity in innovation matching.

Accordingly, the paper cannot claim novelty from any of the following alone: public versus private objectives, endogenous intermediary choice, matching/search investment, or market-creation versus displacement.

The only architecture that survives Stage 0R-MC is **MC-A with a narrowed matching primitive**:

1. two regional governments choose costly local matching-evaluation capacity `x_1,x_2`;
2. the capacity operates over regionally anchored heterogeneous partner sets `S_1,S_2`;
3. overlap is derived from the intersection/common component of those partner opportunity sets, not from an arbitrary cross term;
4. a metropolitan private intermediary has an exogenous partner technology in the baseline and chooses access price `p_T`;
5. heterogeneous immobile projects choose among `H_1,H_2,H_T,0` by expected net productive match surplus;
6. the no-hub option generates a genuine uncovered-project margin;
7. the local government's marginal investment effect is decomposed into creation, match-quality improvement, displacement/duplication, and regional-incidence terms.

## The single unresolved issue

Stage 1R-MC must show that the overlap mechanism is not merely generic type-dependent facility quality in disguise. The required formal result is an **overlap-additionality lemma**: when local partner opportunity sets become more overlapping, a marginal increase in local screening/evaluation capacity must shift a nontrivial share of its effect from previously uncovered match creation toward duplication/displacement, under primitive restrictions on partner capability distributions rather than a reduced-form competition coefficient.

If this result collapses to a generic cross-elasticity or to `quality high enough -> use hub`, the branch is **NO-GO — ABSORBED**.

## Routing

Proceed only to `Stage 1R-MC — Minimal Formal Model Hard Kill`. Do not add referral, Cournot, taste shocks, an endogenous `x_T`, or an arbitrary overlap parameter to rescue the result.