# Stage 25 — CORRECTED GAME: two asymmetric regions, TWO INDEPENDENT municipalities and TWO local private hub candidates; firm competition exists without hubs
Date 2026-10-27. STATUS: source-verified research-model recommendation, a constructed zero-hub Cournot benchmark; NOT a solved stage-25 SPNE.

## User correction (MANDATORY)
Region A urban with municipality a and private operator 1; region B rural with municipality b and private operator 2. FOUR distinct stage-1 simultaneous strategic actors each choose {build in OWN territory, out}; NO participant allowed to build its hub in the other region. Profile (g_a,p_1,g_b,p_2)∈{0,1}^4 yields 16 subgames, including PP, PG, GP, GG, double occupancy. a optimizes local W_A, b local W_B (not a single prefecture optimizing W_A+W_B). Planner values W_A+W_B (plus any unassigned spillover correct attribution). Local governments' behavior is DIFFERENT from Stage23 single-prefecture g=(g_A,g_B) choosing jointly.

Crucial: ZERO HUBS must NOT mean zero firms or zero economic activity. Producers sell output and choose own R&D WITHOUT hubs. Opening hubs alters incremental innovation opportunity/knowledge links, potentially competition, not existence of production market.

**Operator identity must be fixed in manuscript**: (1) private hub candidates 1,2 are also the two product-market producers, integrated ownership; or (2) they are standalone for-profit hub operators, while n_A and n_B incumbent downstream production firms exist regardless of hubs. If one firm only per region and firms cannot sell across regions, LOCAL COMPETITION CANNOT EXIST; model needs >1 producer per region or interregional cross-selling in the product market. Prohibit interregional HUB establishment but do not automatically ban interregional GOODS sales or cross-region R&D.

## Closest verified models (ranked by exact component)
1. d'Aspremont & Jacquemin (1988), AER 78(5) 1133–1137, **Cooperative and Noncooperative R&D in Duopoly with Spillovers**. Formally p=a-b(q1+q2), C_i=[c0-x_i-βx_j]q_i, R&D cost γ x_i². R&D first, Cournot output second. DOES provide NO-HUB nonzero product competition and endogenous independent R&D, and socially compare R&D arrangements. Hub is an additional institutional R&D collaboration tech that may alter probability, cost or transmission of links. Critical absorbs generic claim "cooperation changes output, investment and welfare".
   https://ideas.repec.org/a/aea/aecrev/v78y1988i5p1133-37.html
   Original author PDF https://www.di.ens.fr/~aspremon/Claude/PDFs/dAsp88a.pdf
2. Long & Wong (2009), Review of International Economics 17(2) 261–279, **A Tale of Two Ports: The Economic Geography of Inter-City Rivalry**. EXACT geographic provider-actor pattern: TWO city governments and TWO local private intermediary/distribution companies; each local government invests in its own public infrastructure; private service operators Bertrand-compete for producer clients; governments compete, sometimes preemptively overinvest. It differs because government provides complement to local private services, rather than independent competing G/P facility whose entry is simultaneous. Crucial template for regional private/public welfare functions, city competition, and private pricing.
   https://doi.org/10.1111/j.1467-9396.2009.00822.x
3. Cremer, Marchand & Pestieau (1997), Journal of Public Economics 65(1) 23–35, **Investment in Local Public Services: Nash Equilibrium and Social Optimum**. TWO municipalities simultaneously decide discrete/in\-divisible public investment with interregional spillovers, compare Nash and social optimum, conditional grants. Not private R&D or 4-owner simultaneous entry; strong absorption of "cross-border benefits cause public underprovision".
   https://doi.org/10.1016/S0047-2727(97)00004-2
4. Haaland & Kind (2008), Journal of International Economics 74(1) 170–187, **R&D policies, trade and process innovation**. TWO countries each has OWN company and public authority; firms compete across markets, conduct cost-reducing R&D; govs choose R&D subsidy. Closest downstream two-producer regional competition skeleton plus local-government welfare. Not four binary public/private hub build choices.
   https://doi.org/10.1016/j.jinteco.2007.04.001
5. Dur & Staal (2008), RSUE 38(2) 160–173, **Local Public Good Provision, Municipal Consolidation, and National Transfers**. CITY/VILLAGE asymmetric size local govt provision, interjurisdictional spillovers; over-/underprovision and public coordination already known. No endogenous downstream Cournot/PP/GP firms.
   https://doi.org/10.1016/j.regsciurbeco.2008.01.005
6. Goyal & Joshi (2003), Games and Economic Behavior 43(1) 57–85, **Networks of collaboration in oligopoly**. Firms form costly pairwise R&D links prior to output-market oligopoly; stable/efficient network structures. Link creation incremental to a perfectly alive non-hub product market; strategically stable asymmetric network despite ex-ante firm symmetry. Direct kill test for link-building as novelty.
   https://doi.org/10.1016/S0899-8256(02)00562-6
7. "Geography of Knowledge and R&D-led Growth" (2022), J Econ Geography 22(6) 1149..., DOI 10.1093/jeg/lbab019. URBAN/lagging regions, local vs cross-region knowledge acquisition, endogenous firm number/market power and innovation growth. Great relevance for asymmetric baseline but dynamic growth model too heavy for first publication candidate.
8. Optional regional cross-market competition note: "A higher-cost region excessively attracts firms" (2020), Journal of International Trade & Economic Development DOI 10.1080/09638199.2020.1799427. First choose geography THEN Cournot competing in each region; user FIXES location of hub providers, so importing strategic firm locations prematurely deviates from user goal.

## Proposed MINIMAL new formal downstream baseline (two producers sell in both markets, hub cross-siting prohibited)
For market k∈{A,B}, inverse regional demand P_k(Q_k)=a_k-b_k Q_k; a_k>0,b_k>0. Big city's market demand can be represented by a_A>a_B or a flatter demand b_A<b_B, keeping city/rural differences independent from competition intensity.
A-based producer 1 and B-based producer 2 exist from outset even if zero hubs. Their base costs c_1,c_2, and bilateral GOODS shipment transport τ>=0, i.e.
  C_1A=c_1, C_2A=c_2+τ; C_1B=c_1+τ, C_2B=c_2.
In both local goods markets k each firm chooses quantities q_{1k},q_{2k} simultaneously. At strict interior:
  q_{1k}^0=(a_k-2 C_1k+C_2k)/(3b_k)
  q_{2k}^0=(a_k-2 C_2k+C_1k)/(3b_k)
  Q_k^0=q_{1k}^0+q_{2k}^0
  CS_k^0=b_k (Q_k^0)²/2
  π_i^0=Σ_{k} b_k(q_{ik}^0)² - (own independent baseline R&D cost).
If quantities nonpositive, resolve complementary slackness/global nonnegative q — do NOT just use interior formulas.
These formulas PROVE positive ordinary economic activity with no hub, and DIFFERENT local production market conditions due to demand/cost/trade. **If one producer per region and no cross-selling**, local markets are monopolies, NOT intra-regional competition. To study competitive intensity n_A>n_B independent of market size, use n_A>=2 local producers and n_B>=1 or segmented oligopoly with cross-selling.

Add independent pre-hub R&D e_i>=0 and knowledge leakage baseline θ0∈[0,1], costs c_i(e)=c_i^0-e_i-θ0 e_j, real convex investment cost h_i e_i²/2. Then backward induction Stage output -> Stage R&D yields nonzero investment and positive CS and producer profits at zero hubs. Hub action vector h=(g_a,p_1,g_b,p_2) alters **incremental** R&D link formation and diffusion (not base market existence), e.g. own-region experimental site availability, search-cost reduction, and realized innovation link probability η(h). DO NOT simply add payoff term +λ when government owned, or social public welfare superiority becomes built in.

Suggested stages:
 0) Nature fixes region population/demand a_k,b_k, marginal costs c_i, investment efficiency h_i, trade cost τ, background knowledge transmission θ0, public and private site real fixed costs and social incidence.
 1) 4 ACTORS choose local hub investment SIMULTANEOUSLY, 16 choices (both local authorities independently).
 2) Firms choose link formation/participation (conditional on existence) and possibly local hub use, offsite R&D collaboration; if no hub can still use baseline conventional collaboration at θ0. Pricing fixed minimal first; then endogenous.
 3) Producers choose R&D e_i to cut costs and benefit from links; baseline R&D remains when h=0.
 4) Producers play Cournot in both regional product markets; calculate CS_i, producer profits and operator fee revenue.
 5) Realized consumption/incomes; planner computes W_A+W_B.
Solve Stage4→3→2, then four-way simultaneous Stage1 NE.
Local government a: W_A=CS_A+profits of A-resident firms (including profits from sales in B)+local taxes/transfers/knowledge spillovers−F_GA, whereas b maximizes W_B. Private operator P_i: choose build if PRIVATE total operating payoff after hub minus baseline private continuation payoff≥F_Pi (if integrated producer, include change in its production profit). If operators separate from participating manufacturing firms, their own payoffs are PLATFORM fee receipts net F, not product Cournot profits. Economic planner maximizes W_A+W_B over all 16, keeping downstream Nash continuation fixed: second-best institution vs global planner control of R&D/quantities must not be conflated.

## Publication-oriented thesis and falsifiable gate
"Do two region-specific public and private hub-investment decisions generate socially too many, too few, or wrong-owner sites when private product-market competition and R&D pre-exist, and cross-regional knowledge-sharing opportunities and asymmetric local objectives make local government and hub firms strategically substitutes/complements?"
This question is plausible but 100% NOT novelty certified:
- d'Aspremont/Jacquemin already R&D spillover and product competition.
- Cremer et al already two discrete local public facilities externalities / planner gap.
- Long/Wong already TWO govs + TWO private local firms strategic.
- Dur/Staal already urban/rural asymmetric public under/over.
- Goyal/Joshi already endogenous collaborative links before product competition.
Only a sharp NEW region-specific **ownership/entry** condition that cannot be derived by each parent model separately can clear absorption gate.

## Most IMPORTANT model-design decisions
(1) Are firms 1,2 themselves downstream competing PRODUCERS as well as private hub operators, or only private HUB intermediaries with a distinct set of incumbent producing firms? The question changes the game fundamentally.
(2) Can the producer goods be sold to the other region? "No other-region entry" here concerns hubs; cross-region goods sales may be allowed. If NO, create n_A,n_B≥2 or differentiated local competitors for actual local rivalry.
(3) Local government objective local W_i (could include firm's extraregional profit) vs universal W_A+W_B, and national planner W_A+W_B, must be separated.
(4) Are matched innovations product-market rivals (cost-reducing collaboration) or firms from complementary industries co-inventing NEW products? AJ/Goyal network assumes rival-product cost effects; joint-new-product models require different demand structures.
(5) How public hubs differ from private hub owners: pricing/appropriation, neutral link governance/contractibility, cost, cannot assume γ_G=1>γ_P without derivation.

Recommended **minimal start**: combine d'Aspremont/Jacquemin (1988) for zero-hub endogenous R&D/Cournot, Cremer et al (1997) for decentralized binary local governmental choices, and Long/Wong (2009) for two regional private intermediaries and governmental external interactions. Introduce link-investment mechanism next (Goyal/Joshi). Preserve all 16 pure entry profiles and government objectives local rather than prefectural. Never assert already solved Stage23 applies to this game.
