# Stage 7.5 v2.1 — Full-Theory Freeze Decision

Date: 2026-09-10 JST

Canonical workflow: `ryotamatsuki/research-paper-workflow` v2.1 @ `fe6fa5f2a94632be578d5e9fab39e7df86013113`

Canonical template: `templates/STAGE_075_FREEZE_DECISION.md`

Upstream authorities:

- Stage 4 repaired construction: GO at the repaired witness.
- Stage 4A: `GO — MATHEMATICAL ADVERSARIAL CERTIFICATION PASS`.
- Stage 6: `GO — GO TO WELFARE / GENERALITY`, with broad novelty claims permanently killed.
- Stage 7: `GO TO STAGE 7.5`, with welfare/generality claims narrowly classified.

Working title: `Strategic Interaction among Public Innovation Hubs under Private Repricing`

Target journal retained for value assessment: Journal of Public Economic Theory (JPET), with a lower theory journal as fallback if referees view the result as too narrow.

## 1. Executive freeze-decision verdict

**GO TO STAGE 7.5A GENERALITY / QUANTIFIER RED-TEAM.**

The project still warrants full-paper investment, but only under a narrow proposition-level contribution. It is not a new theory of public platforms, two-sided markets, downstream feedback, or intergovernmental competition in general.

The full-paper case rests on five elements taken jointly:

1. a clean comparative object — the same two public investment choices under a matched private price level;
2. a qualitative strategic reversal — `BR_B3' > 0 > BR_G'`;
3. isolation of the endogenous repricing-response margin rather than a price-level difference;
4. a repaired all-regime global-equilibrium witness independently certified at Stage 4A;
5. a non-mechanical welfare implication — a positive cross-regional coordination wedge at the repaired full-game equilibrium after the private fee transfer is removed.

The result is narrow, but it is more than a parameter exercise because the mechanism is economically interpretable without notation, the nested benchmark identifies one strategic margin, and the welfare consequence comes from an externality rather than from the private transfer itself.

JPET remains a defensible primary target, but not a comfortable one. The contribution is now best described as a compact public-economic-theory paper with platform/IO ingredients. A skeptical referee may still regard it as an application of a generic multistage feedback principle. Journal of Economics or a comparable general theory outlet remains a credible fallback if JPET judges the novelty breadth insufficient.

## 2. Mechanism card

| Element | Frozen Stage-7.5 statement |
|---|---|
| Phenomenon | Two decentralized public innovation-hub investments can be strategic complements when a shared private alternative's price is held fixed, but strategic substitutes when that private alternative can reoptimize its price. |
| Friction / environment | Both jurisdictions interact with the same priced private two-sided alternative; public facilitation changes project routing and support-side participation. |
| Strategic response | Stronger public facilitation reduces demand for the private route; the profit-maximizing private intermediary cuts its access price at the repaired witness (`dp_T^*/dx_i < 0`). |
| Equilibrium effect | The private repricing feedback adds a negative cross-effect large enough to reverse the sign of the public-public best-response slope relative to the price-matched B3 benchmark. |
| Welfare effect | At the repaired G equilibrium, each government is locally at its own optimum but does not internalize the positive effect of its investment on the other region; the aggregate marginal welfare effect of increasing one public investment is positive. |
| Empirical implication | Jurisdictional public-investment responses should depend on the price flexibility of shared commercial intermediation options; stronger local facilitation should be associated with a lower optimal commercial access price in the model. |

Minimal causal chain, without notation:

> One region strengthens its public intermediation. That takes demand from a shared commercial intermediary. The commercial intermediary responds by cutting its price. The cheaper commercial option changes routing and network participation in both regions. Because the other government's marginal return to its own public investment is now different, the sign of government-to-government strategic interaction can reverse.

This explanation does not depend on institutional labels such as prefecture, incubator, or coworking hub.

## 3. Essential versus tractability assumptions

### Economically essential for the headline mechanism

1. **Two decentralized public decision makers.** The strategic object is the response of one public investment to the other's choice.
2. **A shared third-party private alternative.** Both jurisdictions must be linked through the same commercial option for one jurisdiction's policy to affect the private response relevant to the other.
3. **Endogenous private behavior in G and passive matched behavior in B3.** Without switching on the private response while matching the on-path price level, the identification comparison disappears.
4. **A channel from public investment to private demand and hence to the private optimum.** At the repaired witness this gives `dp_T^*/dx_i < 0`.
5. **A positive fixed-price public cross-effect somewhere in the admissible region.** Otherwise there is no complement-to-substitute reversal to explain.
6. **A sufficiently strong opposing repricing effect in G.** This is the sign-reversal condition in economic terms.

### Important for the repaired global witness but not the conceptual local mechanism

- `kappa_L + tau > v + alpha` at the repaired witness. This globally dominates nonresident use of the rival public hub and removes the free-riding deviation that killed the old vector. It is a transparent sufficient condition for the certified witness, not a general institutional fact and not the source of the repricing mechanism itself.
- The repaired numerical values `beta=.01`, `gamma=.825`, `tau=.35`. They construct one globally certified equilibrium witness; they are not empirical calibration.

### Normalization / tractability devices

- uniform project types;
- linear support-side participation feedback;
- quadratic public facilitation cost;
- zero real operating cost for the private intermediary;
- symmetric regional primitives at the baseline;
- equal attribution of national support-side surplus across regional welfare objectives;
- scalar access price and primary-route single-homing abstraction.

These devices cannot be called innocuous globally until Stage 7.5A completes its quantifier/generality attack.

## 4. Core propositions and certification status

| Claim | Status entering 7.5 | Maximum defensible statement |
|---|---|---|
| Repaired global-equilibrium witness exists in G and B3 | Stage 4A computationally certified via independent all-regime evaluator | existence at the repaired baseline vector; no global parameter-space uniqueness claim |
| `BR_B3' > 0 > BR_G'` at repaired witness | Stage 4A independently reproduced with multiple finite-difference step sizes | baseline-witness strategic sign reversal |
| small-positive-`beta` analytic reversal on regular branches | analytic local theorem from earlier proof work | local/branch-specific sufficient-condition theorem; not a global SPNE theorem and no explicit primitive `epsilon` characterization |
| `M_G = M_fixed + P_price` mechanism decomposition | algebraic/mechanism identity in the maintained smooth branch | mechanism exposition; repricing term itself is not novel |
| private price response is negative at repaired witness | Stage 7 numerical derivative | model prediction at the certified witness only |
| private access fee cancels in aggregate welfare | exact accounting identity under zero real private operating cost | proved under baseline accounting |
| positive local national-welfare derivative at repaired G equilibrium | Stage 7 numerical verification plus accounting identity | local under-provision in the coordinated `+x_i` direction at the repaired baseline |
| `W_G^N > W_B3^N` | Stage 7 one-witness numerical comparison | numerical witness only; no welfare-dominance theorem |

No headline mathematical claim used for the full-paper value case lacks a Stage-4A-compatible certification status. Broad generality claims remain outside the frozen contribution until Stage 7.5A.

## 5. Baseline versus robustness versus intended-general-theorem classification

### Baseline-specific certified results

- repaired all-regime G/B3 equilibrium witness;
- opposite public best-response slopes at that witness;
- negative private price response at that witness;
- positive local coordination wedge at that witness;
- G-versus-B3 aggregate-welfare ranking at that witness.

### Analytic sufficient-condition result

- sufficiently small positive network interaction yields the strategic sign reversal on regular interior continuations under the stated beta-zero regularity/SOC/nonsingularity conditions.

This is a **local sufficient-condition theorem**, not a global-equilibrium theorem over the full primitive parameter space.

### Robustness / generality not yet certified

- small `C^1` distribution changes;
- nonlinear network functions;
- regional asymmetry;
- a neighborhood of the repaired global-equilibrium witness;
- arbitrary positive `tau` or direct inter-hub project competition.

These are Stage-7.5A attack targets, not frozen contributions.

### Permanently stale evidence

The old `±0.5%` 20-draw exercise around the rejected pre-v2.0 canonical vector is non-authoritative and must not appear as current robustness evidence.

## 6. Closest-paper distinction

The broad ideas are prior art: fiscal-competition papers already show that downstream strategic policy responses can change or reverse public-investment conclusions; generic sequential-game work shows follower responses feeding back into earlier investment incentives; monetary-policy work shows a third-party feedback rule can turn strategic complements into substitutes; platform work already contains two-sided pricing, investment, and public/private competition. The surviving distinction is therefore not that feedback, repricing, strategic substitutes, or public/private platforms are new. The narrower result is that, for the **same two decentralized public investment choices**, one can hold the shared private intermediary's price at the **same on-path equilibrium level** and obtain strategic complementarity, then activate only the intermediary's profit-maximizing repricing policy and obtain strategic substitutability. The strongest close papers examined in Stage 6 do not supply this matched-price public-public sign comparison as a direct theorem, restriction, or relabeling.

That paragraph is the maximum defensible novelty framing unless Stage 7.5A requires further narrowing.

## 7. Welfare and generality case for a full paper

The welfare content is substantive enough to support, but not lead, the paper.

At the repaired G equilibrium:

- the government's own reduced marginal welfare is approximately zero;
- the rival-region marginal welfare effect is about `+0.50053`;
- the private-profit effect is about `-0.01008`;
- the aggregate effect is about `+0.49045`.

Thus the decentralized equilibrium exhibits a real positive cross-regional coordination externality after the private-price transfer is cancelled. This prevents the paper from being merely a statement about the slope of two reaction functions.

However:

- the equal regional attribution of support-side surplus materially influences the welfare wedge;
- no first best or global social optimum is solved;
- strategic substitutability is not itself socially desirable;
- G's higher aggregate welfare than B3 at the witness is not a theorem.

Institutional portability is credible across at least two innovation-intermediation settings already documented at Stage 7, but it is portability of interpretation, not a general theorem across sectors.

## 8. Fatal and major referee risks

### R1 — Generic-feedback novelty objection — **MAJOR, not fatal**

A referee can say: “This is another example of a known multistage feedback mechanism.” This is the strongest novelty threat. The response must rely on the matched-price public-public comparative object and not on the generic feedback idea.

### R2 — Strong `tau` condition weakens the literal hub-competition story — **MAJOR**

The repaired witness makes nonresident use of the rival public hub globally dominated. Therefore the certified global result does not represent direct public-hub project poaching. The intergovernmental strategic interaction is mediated through the shared private route and support-side network. Manuscript motivation must be aligned with this fact.

### R3 — Global existence evidence is constructive rather than a primitive characterization — **MAJOR but acceptable for target level**

Stage 4A certifies a global-equilibrium witness, not an iff parameter region. The earlier Level-3 projection remained economically impractical. For a compact JPET-type paper this is potentially acceptable if the local analytic theorem and computational global witness are clearly separated.

### R4 — Positive B3 slope is small — **MODERATE**

At the repaired witness `BR_B3' ≈ +.00397`, much smaller than the old rejected witness. It is separated from zero in the independent derivative checks, but Stage 7.5A should attack whether manuscript prose suggests a quantitatively strong effect or unjustified robustness.

### R5 — Welfare result depends on regional surplus attribution — **MODERATE**

The local coordination wedge is not invariant to arbitrary welfare-accounting rules. It should remain a baseline welfare implication rather than a universal policy theorem.

### R6 — Current manuscript is stale relative to repaired authority — **PROCESS BLOCKER, not research-value blocker**

The current paper still contains the rejected old canonical vector and stale robustness/welfare numbers. This prevents submission but does not affect the Stage-7.5 value decision. If 7.5A passes, Stage 8 must freeze the repaired theory before downstream manuscript integration.

No fatal attack currently requires a Stage-3/4 model pivot.

## 9. Full-paper value assessment

### Can the mechanism be explained without notation?

**Yes.** A local public investment weakens demand for a shared commercial intermediary; the intermediary cuts price; this common price response changes the other government's marginal return to public investment and can reverse the government-to-government strategic relation.

### Is the result more than model-specific algebra?

**Yes, narrowly.** The economic logic is a third-party market-response mediation of intergovernmental competition, and the matched-price benchmark isolates that mediation. The exact sign reversal is parametric/baseline-specific, but the causal architecture is economically recognizable.

### Does a credible alternative formulation already survive?

**Partially, not fully.** The mechanism is supported by two distinct verification routes — a local analytic small-network theorem and an independently reconstructed all-regime global witness — and is institutionally portable across more than one innovation-intermediation setting. A fully nonbaseline functional-form theorem is not yet certified. This limitation routes to Stage 7.5A rather than killing the full paper.

### Is welfare/organizational relevance substantive?

**Yes.** The coordination wedge is not merely a transfer artifact and is conceptually distinct from the reaction-function sign.

### Is the closest-paper difference cosmetic?

**No, but it is narrow.** The matched-price switch in a shared-private-intermediary environment defines a distinct comparative object. The contribution cannot be defended simply by saying the ingredients have not previously been combined.

### Would a skeptical referee see more than a parameter exercise?

**Probably, if the paper is rebuilt correctly.** The local theorem, mechanism decomposition, nested benchmark, independent global witness, and welfare wedge together create a coherent theory paper. If the manuscript instead leads with the repaired numerical vector, the paper will look like a parameter exercise and should be rejected.

## 10. Recommended journal level

### Primary

**Journal of Public Economic Theory — defensible but borderline.**

Why retain it:

- the strategic object is intergovernmental/public investment;
- there is an analytic mechanism theorem plus a globally certified constructive witness;
- the welfare externality is genuinely public-economic;
- the paper is compact and theory-forward.

Why it is not safe:

- novelty is proposition-level rather than building-block novelty;
- generic feedback-induced reversals are known;
- global-equilibrium existence is constructive at one repaired vector;
- the `tau` dominance condition limits the literal public-hub competition interpretation.

### Fallback

A general theory outlet around the **Journal of Economics** level remains a strong fallback if JPET views the contribution as too narrow. The project does **not** justify escalating to a more theory-heavy outlet before submission; additional algebra has poor expected journal-value payoff unless requested by referees.

### Research-note downgrade

Not recommended at this stage. A research note would become preferable only if Stage 7.5A shows that the local theorem/generality prose cannot be stated cleanly without making the result essentially “one numerical example,” or if the matched-price novelty is later absorbed by a direct predecessor.

## 11. Exact Stage-7.5A input package

Stage 7.5A must receive the following frozen inputs:

### Model

- same players, timing, utilities, participation system, welfare objects, and strategy domains as the repaired model;
- no new variable or strategic margin;
- repaired witness `beta=.01`, `gamma=.825`, `tau=.35`, with all other primitives as in the current baseline.

### Mathematical certificates

- `reviews/STAGE_04A_V21_REPAIRED_GLOBAL_CERTIFICATION_2026-09-08.md`;
- Stage-4A clean-room all-regime evaluator and successful CI run `34150133823`;
- previous Stage-4A counterexample to the old vector as a permanent regression target;
- local analytic small-`beta` theorem and mechanism decomposition from the existing proof artifacts.

### Surviving claims

1. baseline repaired global-equilibrium witness in G and B3;
2. `BR_B3' > 0 > BR_G'` at the repaired witness;
3. price-matched benchmark isolates the repricing-response margin;
4. local small-positive-`beta` sign-reversal theorem under explicitly stated regularity conditions;
5. transfer cancellation in aggregate welfare under zero real private cost;
6. local positive coordination wedge at repaired G witness;
7. institutional interpretation only as plausibility/portability, not causal validation.

### Benchmark register

- `G`: endogenous-pricing full game;
- `B3`: matched-price fixed-price identification benchmark;
- no first-best claim;
- no solved global social optimum;
- local coordination wedge only.

### Baseline / robustness / theorem classifications

Carry forward exactly the classifications in Sections 5 and 7 above and in Stage 7. In particular, `C^1` distribution robustness, nonlinear-network robustness, asymmetry, and repaired-neighborhood robustness are **not yet certified**.

### Mandatory Stage-7.5A attacks

1. local versus global quantifiers in the small-`beta` theorem;
2. equilibrium versus stationary-branch wording;
3. existence versus uniqueness wording;
4. whether any continuity argument improperly upgrades local strict signs into global SPNE claims;
5. distribution/network/asymmetry generality claims;
6. exact meaning of “matched price” and prohibition on interpreting B3 as a regulation experiment;
7. benchmark terminology (`first best`, `social optimum`, `underinvestment`);
8. `tau` scope and prohibition on claiming certified direct inter-hub user competition;
9. stale old-vector robustness/welfare statements;
10. quantitative wording around the small positive B3 slope.

## 12. Canonical verdict

**GO TO STAGE 7.5A GENERALITY / QUANTIFIER RED-TEAM.**

This is a full-paper value GO, not a theory freeze. Stage 8 remains blocked until Stage 7.5A independently certifies the exact quantifier, generality, benchmark, and prose scope of every headline claim.
