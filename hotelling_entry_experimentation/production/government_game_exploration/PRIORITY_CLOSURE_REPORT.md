# Stage 11R P1h — high-priority audit closure and independent JICT novelty verdict

**2026-10-28 JST.** Model: frozen **Draft PR #61** `cb92ec3cc5bb6fb94b20420bc5085605064d00cc`, built on PR #62–#64 exploratory code. No changes to the baseline manuscript, actor set, technology, posterior sorting, timing, monetary instruments, Cournot conduct, or original baseline parameters. Original `p0=0.05` remains the manuscript benchmark. **User rule:** when a paper's complete original text is inaccessible, terminate retrieval, DO NOT infer absence of priority.

## Executive decision

**Reject the neighboring-private-price-rise/participation-decline result as the primary stand-alone JICT contribution at present; retain it as a rigorously certified optional supplementary proposition.** The mathematical evidence is now substantially stronger: not only the weak-`gamma` first-order effect but also **finite `gamma=0.1`**, global P/P price Nash, **two distinct unique first-stage four-entrant SPNE** and seven certified `p0` points. However:
- `p0=0.05`, the original frozen calibration, has a **price increase and participation INCREASE**; the decrease uses different `p0`, e.g. `0.001`.
- The actual neighboring participation effect is `~2.8e-5` of founders (0.0028 percentage points), economically minute; national match volume **increases** after A public entry. No national innovation decrease is established.
- The high-priority published original **Mekonnen (AEJ Micro forthcoming)** proves public-option changes can raise private prices, the private supply cutoff, and lower welfare under rationed capacity. Our **cross-border two-intermediary** mechanism is different, but common generic claims of price rises/crowding-out are **already known**.
- The paywalled Hua (2016), Zhou (2023), and SSRN-blocked Anderson–Bedre Defolie (2025) cannot be checked proposition-by-proposition. Under user instruction, acquisition efforts are **CLOSED, NOT VERIFIED**. This prevents certifying an independent priority claim.
- The four-entry cost phase diagram is algebraically complete for **pure** equilibrium with the already globally verified downstream pricing, but is the familiar consequence of binary participation and separable entry costs; no distinct general economic mechanism is created.

**Stage11 scholarly novelty: HOLD. Model and PR #61 unchanged.** This is a research-closure verdict, not a formal journal decision or a claim of actual acceptance probability.

## 1. Closure of the high-priority table

| Audit ID | User-facing status | Evidence | What is not claimed |
|---|---|---|---|
| L01 Hua et al. 2016 | **STOP: FULL TEXT UNAVAILABLE** | Restricted SAGE publisher original: abstract + existence-of-supplement metadata | No assertion that its uninspected theorems lack our result |
| L02 Zhou et al. 2023 | **STOP: FULL TEXT UNAVAILABLE** | INFORMS + university author abstract | No proof-to-proof comparison; **IMPORTANT CORRECTION**: their abstract compares one public/one private with **two profit-seeking private providers**, NOT simply with **one private provider only**. The broad crowdout phenomenon, not the exact comparator, is prior art |
| L03 Mekonnen 2026 forthcoming | **FULL-TEXT COMPARISON COMPLETED** | [v4 44-page original](https://arxiv.org/pdf/2507.12779), dated 7 Jul 2026, especially **Propositions 2, 3, 4**, and [AEA official accepted listing](https://topcat.aeaweb.org/articles?from=f&id=10.1257%2Fmic.20250244) | No mathematical nesting claimed: its public option is capacity-rationed, substitutable; ours is an uncapacitated basic service that may coexist with premium and change cross-region network quality |
| L04 Anderson–Bedre Defolie 2025 | **STOP: FULL TEXT BLOCKED** | SSRN full paper PDF returns 403; CEPR and SSRN official abstracts examined | No assertion that its two-sided heterogeneity and endogenous-fee propositions do not absorb ours |
| C05 other-region participation change | **FINITE-GAMMA + COMPLETE STRICT PURE SPNE CERTIFIED** | `finite_matching_spne.py` and seven `p0` points in `finite_matching_p0_scan.py`, exact rational intervals, 16 profiles and global P/P price theorem | Not all-`p0` or all-`gamma` comparative static |
| P11 neighboring premium fee increase | **FINITE-GAMMA + GLOBAL FEES CERTIFIED** | Same tests: B fee strictly increases at all seven distinct feasible `p0` points | Not an academically new generic public-option pricing result |
| E11 private fixed costs above low-cost dominance bound | **GENERAL COST CLASSIFICATION (FROZEN CONTINUATION)** | `entry_cost_phase_diagram.py`; if `F_{Pi}>N_i Delta_i/4`, coexistence state `B_i` cannot be pure NE, while `F_{Pi}<N_i Delta_i/4` implies private-entry strict dominance | No general closed-form for every mixed equilibrium |
| E12 all 16 states | **EXACT NECESSARY/SUFFICIENT PURE-SPNE COST BOXES** | All 16 profiles, 64 inequalities, 4 independent cost axes; examples include `GG` equilibrium above private dominance bound | Any 4-player mixed equilibrium other than the specially analysed low-`F_P` cases remains unclassified; **study stopped** rather than claimed solved |
| R05 original cost sensitivity | **SEVEN POINTWISE ROBUST FINITE-GAMMA SPNE CHECKS** | Same `F_{GB}=0.041`, `F_{GA}=0.1155` versus `0.1156` across 7 original-form `p0` settings | Not a verified continuous interval for the entire `p0` range |

Closure means **mathematically certified** or **explicitly terminated with a known limitation**. It does not mean all possible general theorems or unavailable prior articles have become accessible.

## 2. The strongest new model-internal result: two exact SPNE and opposite region effect

Hold the complete frozen structure fixed: two local governments, two fee-setting intermediaries, matching `gamma=0.1`, `lambda=0.5`, uniform common skill/valuation, original downstream differentiated Cournot, `F_{PA}=0.65,F_{PB}=0.24`, profit share `omega_A=0.8,omega_B=0.2`. Change only the **existing scalar primitive** `p0` from baseline `0.05` to `0.001`, and choose common `F_{GB}=0.041`.

At both public-cost vectors below, a prior original-model lemma gives **strictly dominant private entry** for A and B: `F_{Pi}<N_i Delta_i/4`. Also, at this `p0`, `F_{GB}=0.041` **strictly exceeds both B-government public-entry thresholds** conditional on either A-government action. So B-government nonentry is strict dominant in the *reduced two-government game*. Then A-government has a strict best response, excluding any other **mixed** first-stage equilibrium.

| Fixed `F_{GA}` | Strict **unique first-stage SPNE**, including mixed | Exact fee / founder participation contrast |
|---|---|---|
| `0.1156` | `(P,P)` | B private fee ~0.6407184; B private-only participant cutoff ~0.36018909 |
| `0.1155` | `(B,P)` | B private fee ~0.64084476; B private-only participant cutoff ~0.36021723 |

The strictly positive private-fee difference has outward rational enclosure `(0.0001263554,0.0001263575)`. The strictly positive B participation **cutoff** difference has outward rational enclosure `(0.0000281368,0.0000281389)`, hence the unnormalised mass of B participants falls by that amount.

The government A conditional public-entry threshold is strictly bounded in `(0.1155414128,0.1155414166)`. Because `F_{GA}` straddles it, its best response changes; B government stays out. This is a complete finite-`gamma` SPNE comparison, not merely an isolated smooth perturbation to off-equilibrium outcomes.

**Important:** A public entry raises national match volume from ~0.3959047 to ~0.6397828 (in the normalized mass-product units `m_A m_B`). Thus the opposite-region reduction does not imply lower aggregate matching/innovation.

**Magnitude:** B's participation change is ~0.00002814 of its unit-normalized population, about **0.002814 percentage points**, or only ~0.00440% of its private-only participation. It is a mathematically strict but small effect in this parameterization. `F_{GA}=0.1155` is significantly above the original benchmark `F_{GA}=0.09`.

Proof protocol: frozen global P/P nonnegative-fee **best-response contraction** sufficient conditions, four exact rational Poincare–Miranda stationary-face root checks, single-private provider global fee closed form, 16 profile welfare and private payoff intervals, 4 actor deviations per profile, 15 non-equilibrium deviation witnesses per cost vector, and B non-entry strict dominance. No new fee rule was introduced.

## 3. Same-government-cost 7-point sensitivity

Both A government costs `0.1155` and `0.1156`, the same B cost `0.041`, and **all other primitives** besides `p0` are maintained at every case.

| `p0` | `B private fee after A public entry` | B participation effect | Both strict SPNE |
|---:|---|---|---|
| `0.001` | increases by ~0.00012636 | **decreases** ~0.00002814 | BP vs PP |
| `0.005` | increases by ~0.00013142 | decreases ~0.00002310 | BP vs PP |
| `0.010` | increases by ~0.00013776 | decreases ~0.00001679 | BP vs PP |
| `0.020` | increases by ~0.00015042 | decreases ~0.00000419 | BP vs PP |
| `0.024` | increases by ~0.00015549 | **increases** ~0.00000084 | BP vs PP |
| `0.030` | increases by ~0.00016309 | increases ~0.00000840 | BP vs PP |
| `0.050` **original** | increases by ~0.00018844 | **increases** ~0.00003357 | BP vs PP **at the counterfactual government costs**, not the paper's baseline cost vector |

All signs, roots, provider fees and **all 16 first-stage profile checks** are exact rational interval certificates, not the float-only earlier exploratory result. Because only 7 points were checked, this is **not** a certificate for every `p0` in a continuous interval or for a unique `p0` crossing.

The observed fee increase itself is NOT a new general fact, while the regional participation sign switch is genuinely model-specific. The willingness-to-pay and technical-skill identity `z=theta` and low baseline success `p0` remain central constraints.

## 4. Full four-entrant entry-cost analysis (E11/E12)

Fix a downstream technology/market primitive vector for which **every** entry-status continuation is well-defined and the P/P fee Nash is unique and globally valid (this is certified for PR #61's baseline `lambda=1/2`). For each complete configuration `s in {0,G,P,B}^2`, for each binary player `a in {G_A,P_A,G_B,P_B}`, define the cost-free unilateral entry threshold `T_a(s_{-a})` as the payoff with player a active less the payoff with a inactive, **with the full follower fee and founder participation continuations re-optimised in both states**. With separable cost `F_a`, the profile is a pure first-stage Nash **if and only if**:

```
for each actor a:
   if actor a active at s: F_a <= T_a(s_-a)
   else:                   F_a >= T_a(s_-a)
```

There are `16 * 4 = 64` necessary-and-sufficient half-line conditions. Each one-profile feasible cost region is a Cartesian product of **four** half-lines restricted to nonnegative costs—an **axis-aligned box**. All sixteen boxes are recorded with **globally certified** follower-game prices in `all_16_cost_boxes.json`. The union is the set of cost vectors supporting at least one pure first-stage SPNE; overlaps give multiplicity; uncovered zones require mixed equilibrium analysis.

These are exact general structural facts under the specified game and well-defined continuation, not a deep novel game-theoretic theorem.

The provider fee behavior gives sharper regime exclusion:
- `F_{P_i} < N_i Delta_i/4`: private entry is strictly dominant (baseline regime), hence `G_i` alone cannot be a first-stage equilibrium.
- `F_{P_i} > N_i Delta_i/4`: whenever `G_i=1`, private premium entry is strictly loss making, so **coexistence `B_i` cannot occur in any pure equilibrium**; a `G_i`-only state may occur.
- With `G_i=0`, the private-only gross revenue always exceeds `N_i Delta_i/4` under the maintained `u_i>0`, so in the intermediate range own private entrant reacts in the opposite direction to its government (one-to-one service crowdout).
- Example using the original baseline technology: `F_{PA}=0.8,F_{PB}=0.3` (only existing private-cost primitives are changed) and original government costs `F_{GA}=0.09,F_{GB}=0.0405` generates **`(G,G)`** as the unique strict pure equilibrium, certified against all 15 other profiles.

**Boundary caution:** equalities `F_{P_i}=N_i Delta_i/4` can create indifference. Cost-box supports weak NE inequalities with equalities; the CI sample searches strict interiors. We do not claim all mixed four-player equilibria have been enumerated.

## 5. Theorem-by-theorem closest original literature (not mere bibliography)

### 5.1 Mekonnen, `Luck Out or Outpay? Competing with a Public Option`

Full available [2026 revised manuscript (44 pages)](https://arxiv.org/pdf/2507.12779), version 4 dated 2026-07-07, read including numbered results and proof text:
- **Proposition 2 (printed p.11):** private monopolist's buyer cutoff increases continuously and strictly with public-option capacity. Public-option expansion may crowd out private activity.
- **Proposition 3 (p.12):** buyers relying on the public option benefit from expansion, in terms of own allocation/payoff, despite the endogenous cutoff adjustment.
- **Proposition 4 (p.13):** gives the elasticity criterion under which public expansion **reduces** monopoly price; outside it, price increases and high types may be harmed. This is not an all-distributions unconditional price-increase claim.
- Uses *rationing/congestion of a free alternative* and private services substituting for public, versus our *always available basic public access plus optional premium* and cross-region pairing. It does not directly prove a decline in **another autonomous jurisdiction's** founder participation after government entry, or a four-service-entry mixed/pure SPNE comparison.
- Nevertheless **higher private posted fees, endogenous exclusion and welfare deterioration in public/private markets are prior phenomena**. Our exact sign cutoff is a novel expression **within this particular model**, not certified as new economic theory.

### 5.2 Zhou, Huang, Hsu & Guo (2023)

Publisher/university abstracts only. The article's principal comparison is public/private **mixed duopoly versus two profit-maximising private service providers** under heterogeneous quality taste, waiting-time congestion and balking. P1g's earlier wording that this was *single-private-only* versus *public-plus-private* should be **corrected**. It still establishes a general **public welfare objective can yield lower total welfare** phenomenon. No exact theorem nesting without main text, and per user instruction **stop further full-text search**.

### 5.3 Hua, Chen & Zhang (2016)

Official SAGE published page is restricted. The original publisher abstract verifies two-tier free government vs toll private provider, demand affected by quality and delays, mixed duopoly Nash and fiscal coordination; main proof not obtained. **STOP/NOT CERTIFIED**, not a permanent task to chase. Do not infer lack of coverage from access failure.

### 5.4 Anderson & Bedre Defolie (2025)

SSRN shows a 58-page paper but full PDF retrieval returns 403; CEPR/SSRN public abstracts establish endogenous seller-quality stock and consumer participation value, endogenous participation fees and commissions, and participation possibly falling when fee structures change. **STOP/NOT CERTIFIED** on actual theorem nesting. There is close generic quality-price-participation feedback; one may not assert independent priority over this work.

### 5.5 Additional predecessors already screened

Hagiu (2004b), Veiga–Weyl–White (2017), Peitz–Sato (2026), Monaco–Sabarwal (2016), Kalamov–Staal (2023), Zennyo (2020). Generic selection, bilateral network pricing, fee endogenous participation, region public-good strategic complement/substitute and simple 2×2 mixed-equilibrium classifications all predate this model. The exact regional venture-incidence mapping has not been proved a genuinely non-absorbed main contribution.

## 6. The JICT primary contribution test

Official [JICT Aims & Scope](https://link.springer.com/journal/10842/aims-and-scope): applied theory must make clear empirical and economic-policy implications, and industrial competitiveness/innovation/firm strategy are included. Applying to candidate **cross-border private-fee increase and regional participation decline**:

| Independent referee criterion | Evidence | Verdict |
|---|---|---|
| Mathematical validity of finite-gamma proposition | Global PP price Nash + rational first-stage certification | **PASS** at frozen/scanned numerical points |
| Primitive-derived condition | Weak matching formula `lambda Pi u_j ell_i > 2R0 Delta_j` | **PASS local**, not full global primitive theorem |
| Generic novelty of price rise/crowdout | Mekonnen, platform and service literatures already have mechanisms | **NO** |
| Structural distinctiveness of two-region pricing/innovation | Two autonomous regional matching stocks, founder sorting, lower neighboring participation | **DIFFERENT MODEL**, no theorem-level nonabsorption proven |
| Robustness | 7 distinct p0 points; 4 sign-decrease + 3 sign-increase, common SPNE cost pairs | **CONDITIONAL** |
| Economic magnitude | B participation ~0.0028 pp change at p0=.001; aggregate match volume increases | **WEAK** |
| Policy prediction | Cross-region fees and service users would need joint data; no calibration or causal identification | **WEAK** |
| Match to journal and concise manuscript | Topic fits JICT, but 27-page main paper already exceeds typical <20-page recommendation | **CAUTION** |
| Literature priority | Mekonnen fully read; three critical others inaccessible | **NO CERTIFICATION** |

**Final editorial-style decision: No, do not elevate the cross-border price/participation result to the independent PRIMARY JICT contribution on this evidence.** It is at most a **supplementary, model-specific existence/counterexample theorem**, unless a separate broader robust magnitude/structural priority proof or empirically grounded material effect emerges. It must not be described as a theorem of public options generally.

**Recommendation:** Keep PR #61 exactly frozen; archive this P1h stage as a transparent negative novelty gate. Do not continue chasing inaccessible full texts, and do not claim the 81-item ledger was globally exhaustive beyond its declared seven-scope domains. Future JICT submission should be based on PR #61's narrowly stated access-selection-incidence application and independently assessed on that contribution, not on the supplementary price paradox.

## 7. Reproduction / audit trail

- `finite_matching_spne.py`: exact rational price root, 16 state first-stage profiles, full entry deviations, two strict SPNE, private fee/participation deltas, match-volume direction.
- `finite_matching_p0_scan.py`: same certification at 7 ex-ante existing parameter points, both unchanged cost vectors.
- `entry_cost_phase_diagram.py`: 16×4 necessary/sufficient cost-axis conditions, including cost values beyond strict private dominance.
- `EQUILIBRIUM_COVERAGE_MATRIX.csv`: user-visible closure updates for the nine top-priority rows and limitations, 81 overall IDs unchanged.
- `.github/workflows/stage11r-p1h-priority-closure.yml`: regression checks original global P/P fees, finite-`gamma` proof, entry-cost boxes, additional theoretical identities and coverage checks.
- **No manuscript edits, PR merges, submission, change to baseline cost/technology in PR #61, or actual external human-peer referee evaluation**.
