# Stage 4 — Prior-Art Isomorphism Matrix

Date: 2026-09-01

## 1. Label-stripped Stage 4 model

Remove `regional`, `public`, `innovation`, `co-creation`, `hub`, `local`, `outside`, and `knowledge`.

The model becomes:

> A principal chooses a continuous bias parameter for an intermediary. Bias gives the principal a direct favored-side benefit but lowers the expected access payoff of a second side. Agents on the second side enter endogenously under heterogeneous participation costs and congestion. The principal values the thickness of the second-side pool, so excessive bias can lower the principal's own payoff. A planner adds second-side participant surplus and chooses less bias.

This stripped description is the decisive novelty object.

## 2. Comparators

### de Cornière & Taylor (2019), RAND Journal of Economics

DOI: `10.1111/1756-2171.12298`

Their model studies consumers relying on a biased intermediary when choosing among strategic sellers. Bias changes seller incentives; depending on payoff congruence or conflict, the favored seller can improve or worsen its offer and consumer welfare can rise or fall. Policies toward bias are analyzed.

### Zennyo (2022), Journal of Industrial Economics

DOI: `10.1111/joie.12301`

The platform intermediates consumers and third-party sellers while also selling first-party products. It chooses fair versus biased search. Bias changes commission and equilibrium prices, which changes consumer participation; consumer participation then changes seller participation through indirect network effects.

### Rysman (2009)

DOI: `10.1257/jep.23.3.125`

Canonical two-sided-market benchmark: participation on one side changes the value of participation on the other side, and platform choices internalize some but not all cross-side effects.

### Current bias frontier

Reimers & Waldfogel, NBER Working Paper 31766, revised 2026, develop an equilibrium/welfare framework for platform ranking bias.

Chi, Choi, Hahn & Kim (2026), CESifo Working Paper 12749, jointly characterize self-preferencing and transaction fees and extend to consumer search.

Xu et al. (2025) model a continuous matchmaker-bias control factor and nonlinear effects on matching fairness/platform revenues.

## 3. Matrix

| Element | Stage 4 model | de Cornière–Taylor (2019) | Zennyo (2022) | Exact mapping? | Material difference? |
|---|---|---|---|---|---|
| Decision maker | sponsor/principal | biased intermediary environment; seller responses central | profit-maximizing integrated platform | No | Form differs, but a generic principal/platform can be substituted for the sponsor |
| Intermediary | implements `alpha` | advises consumers among sellers | search platform | Partial | Stage 4 intermediary has no independent objective |
| Favored side | local/favored activity | favored seller | first-party product | Partial | institutional identity only |
| Other side | potential outside partners | consumers/sellers depending mapping | third-party sellers/consumers | Partial | Stage 4 has only one endogenous entry side |
| Bias instrument | continuous `alpha` | biased advice/ranking | fair vs biased search | Partial | continuity is not economically distinctive |
| Source of bias | geographically bounded sponsor payoff | intermediary bias/payoff conflict | vertical integration / own content | No exact | This is the main interpretive difference, but it produces no new strategic term in the Stage 4 equations |
| Participation | cost-cutoff entry on second side | not the core abstract mechanism | consumer and seller participation | No exact to DCT; strong to Zennyo | Stage 4 uses congestion cutoff rather than price/commission channel |
| Cross-side feedback | bias lowers entry; pool thinning lowers principal value | strategic seller response to bias | bias -> price/consumer entry -> seller entry | Structural overlap | Generic feedback from bias to participation/value is established |
| Congestion | `1/(1+n)` under fixed attention | not central | indirect network effects rather than this congestion form | No | Different functional form, but this is a standard participation technology rather than a new institution-specific strategic margin |
| Principal payoff | favored-side return `A alpha` + pool value `B n` | intermediary/consumer/seller payoff structure | platform profit from first-party/third-party business | No exact | Stripped economic trade-off is still bias benefit versus participation value |
| Planner | adds entrant surplus | bias-policy welfare comparison | consumer/seller/platform welfare comparisons | Partial | welfare accounting is standard |
| Equilibrium | bias choice -> entry cutoff -> matching value | intermediary/seller strategic response | platform search/commission -> prices -> participation | No one-to-one | Stage 4 is simpler, not demonstrably more novel |
| Main regional result | interior self-limiting bias because excess bias destroys own pool value | bias can help/harm through seller responses | bias changes participation and platform-market outcomes | Structural overlap | Local-own-welfare reversal is a generic demand/participation discipline result |
| Planner wedge | planner chooses less bias | policy efficacy depends on congruence/conflict | welfare effect of search bias | Structural overlap | No unique regional term remains after stripping labels |
| Policy implication | neutrality/fair-access restraint may preserve participation | regulate intermediary bias conditionally | regulate platform self-preferencing/search | Structural overlap | application-specific interpretation only |

## 4. Exact-isomorphism verdict

The model is **not literally one-to-one isomorphic** to de Cornière & Taylor or Zennyo.

Important non-identities are:

- de Cornière & Taylor use strategic seller quality/offer responses rather than a participation-cost cutoff;
- Zennyo uses vertical integration, commissions, endogenous prices and two-sided participation rather than a single entry side with congestion;
- Stage 4 has a non-profit regional principal and no intermediary profit problem.

Therefore the strict `NO-GO — MECHANISM ISOMORPHIC TO PRIOR ART` label would overstate the evidence.

## 5. Structural-novelty verdict

Despite the absence of literal isomorphism, the economically active loop is generic after label stripping:

`bias benefit -> reduced participation/access -> lower value of the intermediary's participant pool -> self-limiting bias`.

The regional/public origin of bias enters only through the coefficient `A` and interpretation of the principal. It creates **no additional strategic state, best response, information problem, or equilibrium term**.

Replacing the regional sponsor with a private platform or any biased principal that values both favoritism and participation leaves the mathematical propositions unchanged.

This is fatal to the intended contribution.

## 6. Innovation-intermediary institutional literature

The institutional interpretation is credible but not novel enough to rescue the theory.

- Klerkx & Leeuwis (2008), `10.1016/j.technovation.2007.05.005`, explicitly discuss governance/revenue tensions and protecting intermediary credibility and impartiality.
- Kant & Kanda (2019), `10.1016/j.jclepro.2019.04.213`, identify neutrality as a core dimension affecting intermediary survival.
- Feser (2023), `10.1007/s11846-022-00593-x`, reviews neutrality/trust and participation/inclusion mechanisms in innovation intermediation.
- Eldebo & Hjelm (2026 volume), `10.1080/14479338.2025.2514462`, explicitly describe publicly funded open-innovation matchmakers as neutral intermediaries serving both sides.

Thus even the claim that neutrality is a productive property of innovation intermediation is established institutional background.

## 7. Final prior-art decision

- Literal one-to-one isomorphism: **NO**.
- Structural overlap with biased-intermediary/platform participation economics: **VERY HIGH**.
- Public/regional sponsorship mathematically essential: **NO**.
- Innovation-intermediary neutrality conceptually novel: **NO**.
- Defensible theorem-level contribution after stripping labels: **NO**.

Recommended Stage 4 classification:

**NO-GO — RESULT TOO WEAK.**

The reason is not algebraic failure. The minimal mechanism works, but its strongest theorem is a generic biased-principal participation trade-off rather than a distinct economics of public regional innovation intermediation.