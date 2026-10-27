# Stage 26 — Where do private co-creation HUB operators earn profits? Evidence-based revenue model

Date 2026-10-27. User's UPDATED architecture = region A urban with independent municipality a and prospective private HUB operator 1; region B rural with independent municipality b and prospective private HUB operator 2. These four actors only establish hubs within their territories. DISTINCT incumbent downstream manufacturers in A and B produce/compete/research even with NO hub. **This distinction is the baseline recommendation.** Optional later vertically integrated alternative has operator itself also produce and internalize Δπ_product; must not quietly use manufacturing profits as broker fee revenue.

## Verified real-world evidence
1. Nesta, Incubation for Growth, §7.2.1: four revenues rent, business-support fees, equity/royalties, sponsorship from public/private; most mixed; due to liquidity-constrained startups revenue capture imperfect. https://www.nesta.org.uk/report/incubation-for-growth/
2. OECD (2026) Incubation in Entrepreneurial Ecosystems, financing, policy mapping: rents and service fees, equity, public support and tendered operating contracts. Of 34 listed financial support policies, grants and procurement regimes differ, many report performance; hence public procurement is a POLICY INSTRUMENT distinct from direct public ownership, should not be automatically credited as costless private operator income. https://www.oecd.org/en/publications/incubation-in-entrepreneurial-ecosystems_21be9d03-en/full-report/mapping-of-incubation-support-policies-in-oecd-countries_932b6122.html
3. SHIBUYA QWS OFFICIAL personal 16,500 JPY/mo (project), 27,500 JPY/mo (standard), business/public 110,000 JPY/mo; per-room ancillary charges. Good two-tier subscription evidence. https://shibuya-qws.com/entry
4. QUINTBRIDGE by NTT WEST OFFICIAL membership free; co-investment / PoC projects and business creation. A CORPORATE STRATEGIC return model, not evidence that free co-working is commercially self-sustaining from subscription. https://www.quintbridge.jp/join/ and https://www.quintbridge.jp/program/detail/202504211300.html
5. CIC Tokyo OFFICIAL office community plus programme operations with Tokyo Prefecture and others; identifies plausible private contracted hub assistance revenue, does NOT disclose profitability / exact contract amounts. https://jp.cic.com/innovation-program/
6. Japanese Diet 2026 Jun 19 answer about another planned state innovation campus: rent/usage to cover administrative overhead, competitive govt grants, sponsors and program revenue to cover research/venture creation. Context-specific proposed plan, NOT universal incumbent hub economics. https://www.shugiin.go.jp/Internet/itdb_kaigiroku.nsf/html/kaigiroku/000222120260619021.htm
7. Madaleno, Nathan, Overman, Waights (2022) Urban Studies: coworking rent; incubator rent or services + public grants; accelerator often investment/equity. https://journals.sagepub.com/doi/10.1177/00420980211004209
8. 2025 international business-model survey "Exploring the business model of innovation spaces": membership, donations, rental, products and services, institutional support, costs facilities/utility/equipment. DOI 10.1080/2157930X.2025.2466930

## Theory precedents: operator pricing not an ad hoc accounting identity
- Rochet–Tirole (2003) JEEA, Platform Competition in Two-Sided Markets DOI 10.1162/154247603322493212. Fee allocation across two sides matters to enrollment, profit and welfare.
- Armstrong (2006), RAND Competition in Two-Sided Markets, DOI 10.1111/j.1756-2171.2006.tb00037.x: monopoly, two rival platforms and competitive bottleneck; membership and transaction fee; single/multi-home.
- Caillaud–Jullien (2003), RAND Chicken & Egg: Competition among Intermediation Service Providers (2003): registration fee and successful-match transaction fee; multiple registrations can change dominant platform rent.
- Hagiu (2006), RAND Pricing and commitment by two-sided platforms: fixed prices vs variable royalty fees with commitment and early one-side entry.
- Nesta/OECD direct private equity and grant models are empirical institutional alternatives, not substitute for endogenous fee choice.

## Main mathematical recommendation: independent firm-owned local hub operator
i∈{A,B}; private h=P_i if entered; locally available innovator mass I_h≥0 and potential partner mass R_h≥0; operating matched-project intensity L_h≥0 (allocate cross-region matches to hosting operator with contractually enforced shares, and do not double-count). Public hub G_i (if entered) may compete for local participation, no mandated no-fee unless modelled.
Private posted two-sided participation prices f_h^I,f_h^R, optional company sponsor/program package B_h (NOT assumed costless), optional contractible match success fee τ_h. Venue, technology and service costs C_h(I_h,R_h,L_h,x_h); real fixed cost F_Pi.

Profit:
 Π_{P_i}= f_i^I I_{P_i} + f_i^R R_{P_i}
       +τ_i L_i^{contractible}+B_i
       − C_i(I_{P_i},R_{P_i},L_i,x_i)−F_{P_i}.
- In ORIGINAL minimal model, take τ=0,B=0 and cost C=v_I I+v_R R+κ_i x_i²/2, while I,R are ENDOGENOUS membership after fee posting and include no-hub outside option; commission and sponsors are added as robustness / institution-specific.
- More policy-relevant empirical lean version: basic startup access free (f^I=0), corporate/partner fee f^R>0, plus paid custom demo-project contracts B_i; demonstrate equilibrium feasibility and no-hub market baseline. Differential charges across two sides reflect cross-side network effect / liquidity constraints, not necessarily private opportunism.
- If one wants pure success-fee, specify who signs, exact fee base and no double-counting (cross-site fees split per contract). Returns from noncontractible knowledge spillovers are not operator revenue.

Useful microfoundations for sponsor/program revenue:
Corporate client j chooses pay B_j for staff programme / innovation challenge; utility V_j^{programme}(match quality, selected firms)−B_j, participation B_j≤V_j. Operator optimizes menu B_j and matchmaking investment; creates genuinely third group of fee-paying clients if necessary (complexity risk).
Public procurement T_i chosen by municipality as a policy instrument; local W_i subtracts financing and any marginal deadweight public-funds cost; T_i is a transfer in consolidated sum π_h+W_taxpayer, not new social output.

Firms WITHOUT hubs still earn production profits π_k^0 and baseline R&D link spillovers θ0. A firm only joins if incremental network/use value−fees−access/travel exceeds preexisting NO-HUB choice. A market without hubs cannot artificially be normalized to zero industry production / research.

## Optional "operator is also manufacturer" alternate
Corporate-owned hub operator 1 is an incumbent product-market manufacturer; its true objective:
 Π_{P_i}^{integrated}=Π_{P_i}^{platform}
    +[π_i^{product}(g,p,f,x)−π_i^{product,nohub}]
    − any cross-division double-counted transfers.
If zero membership (e.g. QUINTBRIDGE), it may still build based on Δπ_product, corporate client access, capability and new venture equity. BUT this changes downstream firm's participation incentives (own and rivals) and creates vertical foreclosure/selectivity; treat as separate, not silently merge with independent-broker model.

## Fee structure versus sunk pre-match experimental investment: mathematically critical
- Lump-sum memberships paid BEFORE investment, given attendance, do not directly change e_i* conditional on q, although they change attendance/matching probabilities and thus indirectly e.
- Fixed per-match success charges independent of realised output and e can shift participants' participation but do not alter effort FOC, when charged as an additive constant per match and match probability independent of e.
- Proportional contractual success commission τ on gross Ae BEFORE ex-post sharing implies expected investor return β(1−τ) q Ae. Ex ante best-response
   e*=β(1−τ)Aq/c
  instead of βAq/c; direct hold-up-style **investment wedge**.
  Hub expected proportional-fee revenue per innovator = τ q A e*=τβ(1−τ)(Aq)²/c.
  Social project surplus per innovator = [β(1−τ)−β²(1−τ)²/2](Aq)²/c (exclude pure fee transfers). Differentiating τ gives welfare investment distortion for β∈(0,1), τ>0.
- If success commission instead only comes out of the PARTNER'S ex-post Nash share (1−β)Ae with innovator share unchanged, current innovator effort FOC is unchanged; partner entry and bargaining endogeneity may still shift. Need precisely define contract order. Stage23 implicitly posits fixed σ partner fee independent of β without fully resolving off-equilibrium contracting and participation.

## Comparative statics and novel question
Urban greater participating mass can support private fee revenues, rural less. With strong cross-hub NONCONTRACTIBLE spillovers, cross-regional knowledge gains can make public investment socially beneficial while independent private broker does not privately internalize them. However **this underappropriation is established platform/externality theory**; claims of physical govt-owned hub need procurement/subsidy/incomplete-contract failure.
If each region has 1 active hub, gross revenue could rise from cross access, but one-hub local competition cannibalization yields total contribution redistribution. Effect of operator choice on both region's production output also must enter local government welfare.

## Phased recommended research baseline and next-stage gate
M0: no-hub DAJ Cournot+R&D firms exist; 4 independent local build actors; private platform Π=f_I I+f_R R−c_I I−c_R R−F; local participants opt out to baseline market, either region may build 0/1/2 hubs. G chooses local welfare; firms privately optimize platform payoffs.
M1: network links from within and other region; endogenous fee(s) with users making attendance decisions, baseline if none; compare PP/PG/GP/GG social optimum vs 16-profile four-player simultaneous Nash.
M2: platform-endogenous cross-link effort x with real cost, optional paid corporate program contract or success commission; test social/private access wedge and whether private gains can be appropriated via menu/two-sided pricing.
M3: procurement/subsidy and vertically integrated corporate-owned alternative; compare direct ownership not necessarily socially preferable.

**No Stage26 equilibrium asserted.** Sources and motivation only. Publication novelty requires nonabsorbed ownership/entry result with operational revenue/accounting separated from real activity and price/participation global equilibria.
