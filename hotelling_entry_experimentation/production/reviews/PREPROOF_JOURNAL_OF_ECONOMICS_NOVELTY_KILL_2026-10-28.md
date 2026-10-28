# PRE-PROOF Stage 2/6 novelty and Journal of Economics value kill-test
Date: 2026-10-28 JST
Object: PR #69 frozen manuscript, head a463d3340d8e6c2307136a6d64db484a5fefc5fa.
User instruction: **DO NOT prove a chosen result before a literature-based independent novelty and journal-value gate.**
Status: RESEARCH-GO **NOT AUTHORIZED** for any candidate below. This report does not edit the paper, model, theorem, Lean, or numerical witnesses.

## 1. Verdict and immediate correction

Previous assistant proposed working on a "private price fixed: public access value falls, private price reoptimized: value rises" analytical result, and claimed it could push the manuscript toward *Journal of Economics*. This prematurely elevated a **transparent decomposition of a previous proven result** into the next research objective. The proper **pre-proof priority decision is REJECT AS HEADLINE**. It may become a short explanatory derivation if beneficial; no new main-proposition proof workflow should be opened solely for it.

There is no guarantee of acceptance even when an unpublished theorem appears not in the surveyed literature. But a credible ex-ante scientific novelty/priority audit can and should eliminate weak candidates *before* formalization and extensive computations.

## 2. Mechanical reducibility to a one-line chain-rule decomposition

On the exact PR69 conditional P->B comparison, partner full, interior private-only fee:
\`\`\`
d=Delta_i, C=gamma*N_i*N_j,
D=d+eta*lambda, eta=gamma*N_j*omega_i*Pi/2,
U=u_i+gamma*N_j*omega_i*R0,
ell=(f-U)/D when f fixed;
f*=(U+D)/2 when the private seller chooses an optimal interior price;
K=M_i+omega_i*Pi, S0=p_b*M_i+omega_i*R0.
T(ell,lambda)=N_i*(u_i*ell-d/8+d*ell^2/2)+C*(S0*ell+lambda*K*ell^2/4).
T_ell=N_i*(u_i+d*ell)+C*(S0+lambda*K*ell/2).
ell_lambda|f=-eta*ell/D,
df*/dlambda=eta/2,
dell*/dlambda=ell_lambda|f + (1/D)*df*/dlambda = eta*U/(2D^2).
Thus dT*/dlambda = dT/dlambda|fixed fee + T_ell*eta/(2D).
\`\`\`
The decomposition follows by the chain rule, *not* an independent equilibrium theorem. If K<0 and S0<0 and gamma small, T_ell>0; at the same local optimal fee, the fixed-fee slope is negative, and the strategically repriced slope can be positive in the already known profit-incidence band. This is a useful **within-model counterfactual identification**, but not a sufficient standalone originality claim.

Do not overstate that the fixed fee is a feasible equilibrium response to changed technology; it is a controlled mechanism diagnostic only. Do not equate the derivative of governmental welfare with an envelope theorem for the firm's private revenue.

## 3. Closest **primary** predecessors — mechanism-level tests, not keyword matching

| Paper; publication details | Actual established price/selection/welfare result | Risk to proposed A as independent JE contribution |
|---|---|---|
| Damiano & Li, **Competing Matchmaking**, JEEA 2008, https://doi.org/10.1162/JEEA.2008.6.4.789 (official publisher abstract) | Price-sorting heterogeneous matching participants under complementarity; high-quality matchmaking can be insufficiently exclusive and inefficient despite broad participation. | Core skill complementarity + fee-selected matching already old; current A is a short derivative inside a more restricted fee schedule. |
| Damiano & Li, **Price Discrimination and Efficient Matching**, Economic Theory 2007, https://doi.org/10.1007/s00199-005-0058-2 (bibliographic abstract only) | Conditions for whether a profit-maximizing matchmaker's incentive-compatible sorting is efficient. | Plain private-vs-social selection wedge is known. |
| Anderson & Bedre-Defolie, **Hybrid platform model**, RAND 2024, https://doi.org/10.1111/1756-2171.12478 (open original) | Better platform own-product quality increases its seller fee and lowers consumer welfare; platform decides fees and direct product rivalry. | Quality improvement -> *higher* strategic fee -> *worse* participation/welfare already in a rigorous close market-design theory. |
| Peitz & Sato, **Asymmetric Platform Oligopoly**, RAND 2026, https://doi.org/10.1111/1756-2171.70008 (open original) | Analyses endogenous platform fees, entry, welfare, partial coverage and **zero-priced outside platform** in appendix; network and participation comparative statics. | Public/free outside alternative + strategic re-pricing is not virgin territory. |
| Hua, Chen & Zhang, **Competition and Coordination in Two-Tier Public Service Systems**, POM 2016, https://doi.org/10.1111/poms.12555 (publisher abstract, final theorem full text NOT inspected) | Free government and priced private service coexistence, equilibrium strategies, conditions, budget/tax-subsidy coordination of misaligned objectives. | Policy-instrument comparison itself is known; no unsupported claim of theorem independence where inaccessible. |
| Mekonnen, **Luck Out or Outpay? Competing with a Public Option**, AEJ: Microeconomics forthcoming as of 2026-10-28, https://doi.org/10.1257/mic.20250244 (official forthcoming abstract) | Profit-max monopolist strategically congests a capacity-constrained free public option; expanding free capacity can raise private price and reduce consumer welfare, with conditions for gains. | Public option + perverse strategic price/welfare comparison is explicitly known; exact rationing mechanism differs. |
| Liu, Rivadeneyra & Reshidi, **Public vs. Private Payment Platforms**, Bank of Canada WP 2026-10, https://doi.org/10.34989/swp-2026-10 (official full working paper) | Free/welfare-oriented public network can induce private platform's higher fees and harm retained merchants; second-best policy under differentiation/switching. | Strategic public/private price-feedback generic priority claim is already weak. |
| Moraga-González, Hoàng & Motchenkova, **Pricing and Investment by a Two-Sided Monopoly Platform**, Tinbergen 26-015/VII (2026), https://papers.tinbergen.nl/26015.pdf (open working paper) | Prices, private investment, underparticipation, second-best optimal comparisons and commission caps; network benefit vs stand-alone investment distortions. | Joint endogenous price/investment and social optimal gap is crowded. *Working paper*, not yet a peer-reviewed theorem priority verdict. |
| Ikegami et al., **Caps vs. Prices: Regulating Matching Markets with Distributional Disparities**, UTMD-151 (2026), https://www.mdc.e.u-tokyo.ac.jp/en/news-en/9768/ (open working paper report) | Regional matching policy: quantity caps vs targeted taxes/subsidies; counterfactual welfare and feasible matching interventions. | Instrument-level interregional matching welfare comparisons are already directly researched; a new subsidy alternative needs a **new question**. |
| Benassi, Chirco & Colombo, **Income, Demand and Privatization**, *Journal of Economics* 145 (2025), https://doi.org/10.1007/s00712-025-00899-w (open original) | Adapts canonical mixed Cournot model to reservation-price distributions and studies stochastic-income-distribution changes; proves general directions of privatization threshold via demand scale/elasticity. | Actual JE benchmark: clear new microfoundation/threshold not a fixed-vs-optimized-price decomposition of a known equation. |
| Angelini, Castellani & Zirulia, **Platform investment and seller competition in two-sided markets**, *Journal of Economics* 144 (2025), https://doi.org/10.1007/s00712-024-00874-x (open original) | Price+investment decisions, strategic feedback, internal seller competition, private participant welfare, hybrid platform and additional forms of investment/competition. | Actual JE benchmark: genuinely new strategic choice/link and nontrivial comparative-static theorem; supplemental numerical analysis is acceptable when mechanism is independent. |

A source *not containing this exact formula* is not sufficient evidence that the economic mechanism is new. No general embedding/impossibility theorem is established here; avoid categorical claims "all literature absorbs" without a theorem-by-theorem map.

## 4. PREPROOF gate applied to the previously proposed A–F expansion plan

| Proposal | Scientific pre-proof verdict | Reason |
|---|---|---|
| **A fixed private fee vs repricing sign reversal** | **KILL AS JE HEADLINE** | Direct chain rule and old PR69 \(T'(\lambda)\); strong close strategic-fee predecessors; no standalone new economic decision. |
| **B parameter-domain + economic magnitude calculation** | **KILL AS INDEPENDENT HEADLINE**, retain ancillary robustness | Necessary and useful if a new theorem exists; does not itself make a new question. |
| **C connect sign to 4-entrant Nash entry** | **HOLD** | A general, open-domain *new* comparative-static/full equilibrium result could matter, but known mixed entry/public option prior and tiny cost margin make existing three-entry witness weak; no candidate new theorem has been independently established as plausible and unabsorbed. |
| **D weaken ability / premium valuation correlation** | **HOLD; NOT PRIORITY CLEARED** | Sorting with correlated types and complements already mature (Damiano-Li; Eeckhout-Kircher). Must identify nonstandard public/region welfare result NOT a generic sorting condition. |
| **E alternative distributions / Bertrand** | **SUPPORTING ONLY** | Robustness tests strengthen an accepted original contribution; not an independent JE theorem and existing negative original Beta/Bertrand evidence must remain disclosed. |
| **F compare free access, subsidies, vouchers, private procurement** | **HOLD; NOT PRIORITY CLEARED** | Hua et al 2016, Bergh 2026, Ikegami et al 2026, Liu et al 2026 already study instrument effects. An exact unanswered policy-equilibrium question must be isolated from their results before specifying/proving a new model. |

**NO positive Stage2/Stage6 originality-and-JE-fit gate has passed** on A–F. The action "begin formal proof A, then B/C, then D..." must be **halted**. Do not manufacture a desired sign by parameter engineering to find a theorem, and do not add a model solely to fit JE.

## 5. Research admission criteria before touching equations

For every prospective theorem, write **one page** containing:
1. **Three nearest predecessor exact claims** (players, objectives, choices, timeline, endogenous price/participation, theorem conclusions and assumptions), with accessible primary evidence. Unknown full text = explicit UNVERIFIED, not "the authors did not prove it".
2. **New economic choice/question**, not merely new parameter or existing welfare threshold's derivative. Identify a *different decision*, policy optimality ranking, inability/possibility theorem or robust reversal that was not resolved by the nearest predecessors.
3. **Formal reduction attempts**: map old demand, type distribution, fee rule, public objective and downstream payoff into predecessors. If candidate follows by substitution + standard chain rule: **KILL**.
4. **Meaningful primitive-domain existence and effect size**: a theorem requiring an artificially minuscule cost interval or a welfare-effect asymptotic to zero is at best a diagnostic, not automatically JE headline.
5. **Editorial JE comparator**: identify at least one actually published Journal of Economics paper with equally model-specific proof/scope and argue why our result has *comparable novelty*. If absent: **HOLD**.
6. **GO/STOP/HOLD before proof**: GO requires concrete pre-identified unresolved theorem question + no immediate absorption; this **cannot guarantee** post-proof publication. A failed proof -> rollback; lack of originality -> stop before coding.

## 6. Closing verdict

- **Fixed-price/repricing candidate A: KILL as main JE proposition.**
- **Other A–F programme: RESEARCH STOP/HOLD**, no math/formalization authorized solely on assumption of Journal of Economics value.
- Existing PR69 JICT-target applied-theory manuscript stays untouched; original Stage12 journal ranking unchanged.
- Next scientific step, if aiming at JE, is **a literature-first open-problem search and candidate novelty triage**, not proof-first exploration or six incremental features. A meaningful answer may be that none is available without changing the economic question; this is preferable to another polished but low-priority paper.

This is an ex-ante priority triage, NOT a proof the proposed theorem has been published verbatim, NOT a guarantee that JE will reject or accept any future result.
