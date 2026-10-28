# Stage 11R P1k — independent referee report: can the incidence-screening inversion be the MAIN economic proposition?

**Date:** 2026-10-28 (JST). Exact scope: frozen PR #61 original manuscript (`cb92ec3cc5bb6fb94b20420bc5085605064d00cc`), PR #67 mathematical addendum and newly verified public literature. No model/parameter/pricing/equilibrium or manuscript modifications here.

## Executive answer — three separate questions

1. **Relative internal theorem quality: YES, this is stronger than the current Propositions 1–2, Corollary 1 and generic parts of Proposition 3 for economic content.** It puts an independently interpretable and sharp **primitive incidence condition** on an unexpected opposite sign between the social value of successful ventures and the government's willingness to restore basic access. A *conditional* central proposition is editorially defensible.
2. **Logical independence: NO.** The exact `dT/dlambda` identity and qualitative `K_i<0, T'_i>0` observation are *already in frozen PR61 §4*. The new claim is a first-order explicit `omega` region, a sharp weak-matching sign condition, and a nondegenerate finite-`gamma` rigorously certified continuous parameter rectangle. It must be stated as a **strengthening/corollary of current Proposition 3**, not an unrelated discovered general theory.
3. **Novelty against economic literature/JICT acceptance: CONDITIONAL, NOT CLEARED.** Damiano–Li (2007/08) have sorting-by-fee with complementary bilateral matching; Hagiu has selection/quality; Takalo–Tanayama–Toivanen (2013) has public/private R&D financiers and cost/spillover-sensitive support; Atal–Cuesta–González–Otero (AER 2024) has localized public-option entry, segmented clientele and private price responses; Liu–Reshidi–Rivadeneyra (2026) has public/private network platform fees. None directly establishes the precise *regional downstream adverse-product-incidence versus improved public-basic-access threshold* under frozen two-region model in accessible texts, but structural non-nesting is not proved. **JICT submission main-novelty HOLD persists pending independent economist referee; a revised editorial main result is supportable, not a certified original general theorem.**

**Recommended publication structure:** front-load a concise **main Proposition: “Incidence–screening inversion under weak matching”** with **explicit quantifier and primitive band** plus a one-paragraph exact finite-`gamma` robust illustration (the user asked to assess candidacy; actual integration into frozen PR61 is NOT done). Relocate old P1/P2 to preliminary lemmas / current Proposition 3 to the threshold construction only if author decides to reorganize. Do NOT append an independent Proposition 4 after current §4 without replacing redundant prose.

## I. Precisely stated candidate theorem and proof

Maintain the frozen economic model: the opposite region is in G or B, so `(m_j,q_j)=(1,1/2)`, local private-only service has an interior cutoff, `theta=z~U[0,1]`, `0<u_i<d_i=Delta_i`, `N_i,N_j,Pi>0`, `M_i<0`, `p_b>0`, `omega_i>=0`, `omega_A+omega_B=1`. Hold the own free basic utility `u_i=v_i+rho*w_i` and `rho` constant while `gamma=rho*mu->0+`; maintain `p_b+lambda<1` or appropriate admissible success bound on a fixed compact `lambda` domain.

Let
```
C=gamma*N_i*N_j,   eta_i=gamma*N_j*omega_i*Pi/2,
U_i=u_i+gamma*N_j*omega_i*R0,   D_i=Delta_i+eta_i*lambda,
ell_i=(D_i-U_i)/(2D_i),    K_i=M_i+omega_i*Pi,
S_i^0=p_b*M_i+omega_i*R0.
```
The true P-to-B **conditional** willingness-to-pay ceiling is
```
T_i(lambda)=N_i[u_i*ell_i-Delta_i/8+Delta_i*ell_i²/2]
          +C[S_i^0*ell_i+lambda*K_i*ell_i²/4].
```
The cutoff's exact derivative is `ell_i'=U_i*eta_i/(2D_i²)>0` when `omega_i>0`. Hence
```
T_i'(lambda)=[N_i(u_i+Delta_i ell_i)+C S_i^0]*ell_i'
            + (C K_i/4)*(ell_i²+2lambda ell_i ell_i').
```
Expansion at `gamma=0+` gives, **uniformly over compact admissible lambda ranges**,
```
T'_i(lambda) =
gamma*N_i*N_j/(16*Delta_i²) * 
 [M_i*(Delta_i-u_i)² + omega_i*Pi*(Delta_i²+3u_i²)]
 + O(gamma²).
```

Define `r_i=u_i/Delta_i in (0,1)`, `omega_i^crit=(-M_i/Pi)*(1-r_i)²/(1+3r_i²)`. If
```
omega_i^crit < omega_i < min{1,-M_i/Pi}
```
with the top endpoint strict and `omega_i>=0`, then **for all sufficiently small positive gamma** and every lambda in the specified compact interval,
```
K_i<0,
S_i^0<0,
T'_i(lambda)>0.
```
The region is nonempty when `omega_i^crit<min{1,-M_i/Pi}`; `0<u_i<Delta_i` implies the economic interval `omega_i^crit<-M_i/Pi` but does NOT by itself guarantee `omega_i^crit<1`. Strict inequalities imply an **open model-primitive neighborhood**. This theorem is NOT a global 4-player SPNE; it compares two local service states with the neighboring provider(s) held in full coverage.

**Independent stronger implication (previous P1j report left implicit): `K_i<0` implies the base venture-welfare coefficient `S_i^0<0` automatically in the original optimal-effort model.** Because
```
p_b*Pi-R0=zeta²*Pi²/(2k) >= 0,
S_i^0=p_b*K_i - omega_i*(p_b*Pi-R0) < 0
```
for `p_b>0` and `omega_i>=0`. Thus **every additional successful-matching term in the P-to-B welfare contrast is strictly regionally harmful**, not only the skill-complement increment. This stronger conclusion depends on net R&D spending being correctly counted in `R0`.

Two counterfactuals that identify mechanism rather than word overlap:
- With `omega_i=0` and `M_i<0`, `eta_i=0`, `ell_i'=0`, so `T'_i=(C M_i/4)*ell_i²<0` at every feasible finite gamma: absent private venture-profit internalization, there is no endogenous price-selection inversion.
- At `u_i=0`, the leading-order sign bracket reduces to `Delta_i² K_i<0` whenever `K_i<0`; positive baseline basic service valuation is essential for **the weak-matching inversion**. It is NOT asserted impossible for every finite gamma outside the weak-coupling regime.
- The critical appropriation threshold strictly falls with `r_i`: `d[(1-r)²/(1+3r²)]/dr=-2(1-r)(1+3r)/(1+3r²)²<0`.

These three restrictions make the result substantially more **economically explanatory** than an existence numerical example.

## II. Stronger positive government valuation: finite gamma proof (not just signs)

On PR61 baseline technology and market primitives, for every point in the **entire closed rectangle** `omega_B in [0.13,0.27]`, `lambda in [0,0.5]`, with `omega_A=1-omega_B`, `gamma=0.1` and A full:
- `K_B<0` and `S_B^0<0`, so **B region's additional venture-related welfare from P->B is strictly NEGATIVE**;
- `H_B>0`, and the full *fixed-cost ceiling* `T_B>0`;
- `dT_B/dlambda>0`, with certified conservative lower bound `>3.34036e-05`;
- the strictly positive `T_B` lower bound is `>0.04036779`;
- conditional extra venture-related gain has certified upper bound `<-1.66554e-05`.

Evidence: independent exact-fraction 100 closed boxes for both derivative and total valuation, integrated with the unchanged Cournot primitives. This is a true uniform-continuum mathematical **certificate**, not a simulation; actual four-entrant pure SPNE persistence over entire omega/lambda rectangle **NOT** certified. This finite-gamma example answers the key logical objection that `T'_i>0` might hold only as matching vanishes, and that government value might be negative despite the derivative.

At fixed original `omega_B=0.2`: `T_B(lambda=0)≈0.04042507100712`, `T_B(.5)≈0.04053167010245` (difference ~0.0001065991). The benchmark cost `F_GB=.0405` lies strictly between them, yielding the previously certified separate strict SPNE outcomes `BP` versus `BB`. The **wide derivative-region proof is NOT a wide discrete ownership-switch theorem**; the threshold effect is small relative to F_G and the baseline ownership-switch band is ~0.2632% of F_G.

## III. Relative evaluation of old vs new propositions

| Existing item | Original actual theorem | Literature prior absorption | Proposed new primitive inversion vs old | Research contribution judgment |
|---|---|---|---|---|
| P1 regional incidence | Under exogenous **full coverage**, `dW_i/dlambda=C K_i/4`, cross-region and aggregate signs may differ | Standard distributive welfare/spillover sign | **New**: welfare can have a negative local technical sign while *the policy threshold derivative* is positive | Far stronger policy mechanism than P1 |
| P2 price screening | Conditional `d f*/d lambda>0`, `ell'>0`, average skill up but number down | Damiano–Li + Hagiu sorting + platform pricing | **New**: the screening channel offsets a negative local commercialization term, with a sharp ownership-share condition | Genuine economic synthesis, not a stand-alone innovation of screening |
| C1 selected relative to random | Upper-tail algebra | Pre-existing sorting integrals | Candidate gives actual government monetary decision threshold, not just statistical sorting | Candidate superior |
| P3 local/national fixed-cost thresholds | `F_G<=T_i` and national-local wedge `Delta W_j` | Standard local public-good & fiscal federalism | **New**: sharp `T'_i>0,K_i<0` sign condition identifies inverse policy response on the P3 threshold | Economically sharper, **but mathematically nested** inside P3 |
| Original diagnostic `dT/dlambda` (§4) | Has **exact same derivative**, already says `K_i<0` and possible `T'_i>0` | Broad sorting + local externality theory | **New only:** primitive necessary/sufficient FIRST-ORDER condition, open feasible profit-share band, uniform finite-gamma certification | A first-order theorem-strengthening, **not independent core math** |
| Three-scenario first-stage switch | `(BP)->(BB)`, versus pooled `BP` | Composite price/entry already studied; exact 3-pattern not matched | **New** robust policy-value monotonicity across a broad existing-primitive rectangle, rather than narrowly chosen cost switch | Better central explanatory pillar than single numeric counterexample; DO NOT claim broad SPNE switch |

**Mathematical novelty: moderate relative to current P1–3, low relative to the already stated exact derivative. Economic expository value: high within current paper. Journal-level stand-alone priority: not established.**

## IV. Mechanism-level prior comparison, not keyword overlap

1. **Damiano & Li (2007), Economic Theory**, [2007 abstract](https://ideas.repec.org/a/spr/joecth/v30y2007i2p243-263.html): single profit-max **matchmaker sets type-dependent entrance-fee SCHEDULES across both sides**; complementary matching and optimal assortative sorting; efficient and privately optimal fee sorting may diverge. **Direct overlap of fee-driven ability sorting + complementarity**. Our distinct choices: autonomous public local free-basic supply, two separate region-specific private monopolists, no assortative platform-managed matching, R&D success and Cournot business stealing into local W, and policy value comparison under P->B. Published theorem text restricted; do not infer original non-absorption.
2. **Damiano & Li (2008), JEEA**, [abstract](https://academic.oup.com/jeea/article-abstract/6/4/789/2295772): two competing **private** matchmakers, costly sorting and more mass versus sorting efficiency. Our two regional intermediaries are cross-side **complementary** rather than competing for the same pool; the government can supply basic access alongside premium. Again do not claim fee/skill matching new.
3. **Hagiu (2011)**: platform chooses quality exclusion and participation; match quantity/quality is a previously known tradeoff. Our private operator chooses only **one linear posted fee**, and ex-post innovation incidence is weighted by local `M_i` and `omega_i`, producing the sign band.
4. **Takalo, Tanayama & Toivanen (2013), IJIO**, [abstract](https://ideas.repec.org/a/eee/indorg/v31y2013i5p634-642.html): strategic public/private **R&D financing**, fixed R&D costs, financial frictions, spillovers and optimal subsidy. They explicitly derive that the effect of spillovers and finance costs on optimal funding/“additionality” is conditional. This is the strongest direct **industrial policy / R&D public support prior** and prevents claiming that “greater spillovers can make support more valuable in otherwise adverse situations” is generally new. But the **instrument** is financing/subsidy rather than a zero-priced basic two-region intermediation option; the marginal privately sorted **access margin** and `M_i` business stealing are not the same mechanism.
5. **Atal, Cuesta, González & Otero (2024), AER**, [official full text and abstract](https://www.aeaweb.org/articles?id=10.1257/aer.20211547): Chilean decentralized public pharmacy market; public option's lower-price/lower-quality supply induces **market segmentation** and **higher private prices** for consumers staying with the private pharmacy. A local government may find public competition beneficial while some private customers lose. No two-region R&D successes and incumbent-versus-venture Cournot market incidence. **Public option + market selection/price response is already empirical prior art.**
6. **Liu, Reshidi & Rivadeneyra (2026), Bank of Canada**, [official](https://www.bankofcanada.ca/2026/03/staff-working-paper-2026-10/): zero-priced/public-payment and private competing **two-sided payment networks** with network fragmentation, private merchant fee increase; public pricing partly disciplined by market-wide welfare. No technology-dependent screened ventures or region-specific downstream product-market damage. Generic price crowdout not novel.
7. **Peitz & Sato (2026), RAND J Econ**: network-fee equilibrium, partial coverage, entry in private platforms; fee and endogenous participation general architecture known.
8. **Cremer et al. (1997) local public investment**: local-national welfare gap and sign-flexible externalities generally known.

The combination is a **different structural game**, not a mere relabeling of a specific published model. But `q=phi(m)=m-m²/2` reduces quality-weighted ability to a single nonlinear participation index, so a sufficiently general nonlinear network-participation/fiscal externality framework might produce analogous conditions. No formal proof of non-embedding in that broad class. **Do not advertise the mere inverse sign as untouched economics.**

## V. JICT editorial gate and defensible presentation

[JICT official scope](https://link.springer.com/journal/10842/aims-and-scope) directly welcomes *applied theory* that highlights **empirical and economic policy implications**, with industrial competitiveness, innovation and firm strategy relevant. The candidate improves the real application to **three separately measurable channels**:
- the heterogeneous founders' **participation cutoff or take-up**, and the provider's **posted fees**;
- ex-ante venture success as a function of **both founders' capabilities**; and
- **regional incumbent displacement, venture profit repatriation/allocation** vs direct basic-service benefits.

Still, the user cannot infer “build a physical hub”, “free services dominate all cheaper instruments”, “local R&D damage must be subsidised”, or “more complementarities always favour government entry” from the result. Distribution `z=theta` and 100% full partner region are strong assumptions; initial baseline JICT paper's Beta(2,2) three-way equilibrium switch fails, and the new theorem has not been certified over that alternative distribution.

**Final decision:**
- **YES, elevate to the principal CONDITIONAL ECONOMIC PROPOSITION / thesis of the paper, ahead of P1/P2,** if the paper is deliberately reorganised around “why socially adverse local innovation can increase the value of free public access via price screening.” This is the strongest *internally derived and mathematically certified economic explanation* so far.
- **NO, do not assert an independent brand-new economic law or automatic JICT-ready novelty.** The complete effect is already implicit in old PR61 `diagnostic-derivative`, and close mechanisms have prior research. Treat it as an analytical strengthening of that previous diagnostic, with one precise primitive-level Proposition and existing P3 as a preceding lemma/threshold identity. Avoid duplicate theorem statements.
- **NO, do not claim the all-primitive parameter result from the finite box proof.** Analytical sign is sufficient/necessary *to leading order as gamma->0*, and finite gamma is rigorously checked only on the original two-parameter rectangle. State that every claim conditions on partner full, local private-only monopoly, uniform aligned skills/willingness-to-pay, and strict model feasibility.
- **Authoring risk:** absent stronger portable result / endogenous selection mechanism under another distribution / empirically meaningful magnitude, remains an applied-theory priority claim of **moderate** originality, not a clear high-novelty result over predecessors.

## Files
- `incidence_screening_sign_certificate.py`: complete rational 100-box proof of `dT/dlambda>0, K<0, S0<0`.
- `incidence_screening_valuation_check.py`: second independent rational 100-box check of `H>0, T>0, venture component<0`.
- New CI `.github/workflows/stage11r-p1k-headline-evaluation.yml`: checks both and original price+entry baseline. No manuscript changes.
