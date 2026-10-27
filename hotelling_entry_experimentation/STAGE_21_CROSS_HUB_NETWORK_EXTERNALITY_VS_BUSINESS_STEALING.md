# Stage 21 — Cross-hub collaboration vs business stealing: proposition-first audit

Date: 2026-10-27
Status: EXACT identities and CONDITIONAL balanced-assignment welfare threshold; NO new SPNE and NO novelty certification.

## Research idea
In the simultaneous-entry game with G and P, a second physical hub (i) competes for participants, (ii) geographically covers otherwise unserved distant innovators/partners, (iii) enables cross-hub collaborations whose complementarities and resulting knowledge spillovers may be privately nonappropriable. The research question is which effects prevail at equilibrium, not merely whether positive externalities exist.

## Critical no-go: interoperability alone does NOT create more pairs than one fully integrated hub
Unit mass I,R; P shares X,Y; G shares 1-X,1-Y, with full single-homing/coverage. Within-hub matches:
 L_within=mu*[XY+(1-X)(1-Y)].
Cross-hub matched links at per-pair efficiency rho∈[0,1]:
 L_cross=mu*rho*[X(1-Y)+(1-X)Y].
Total:
 L2=mu*[1-(1-rho){X(1-Y)+(1-X)Y}] <=mu=L1.
At balanced shares X=Y=.5, L2=mu*(1+rho)/2, compared with L1=mu.
Therefore cross-hub matches at rho<=1 restore losses from segmentation, but cannot by themselves prove that two hubs make MORE matches than a one-hub full-coverage technology with access to exactly the same potential population. An arbitrary positive "extra two-hub innovation bonus" by assumption is not an economic explanation.

## Pre-match investment microfoundation
For j=G,P let matching success potential from own side partners R_j and opposite site partners R_-j be
 q_j=mu*(R_j+rho*R_-j), with q_j<=1.
Gross project payoff Ae; innovator receives bargaining share beta, effort paid ex ante c e²/2, and therefore:
 e_j*=beta*A*q_j/c.
Real expected per-innovator experiment surplus:
 B_j=(beta-beta²/2)*(A²/c)*q_j².
For beta=2/3, k=2 A² mu²/(9c), this equals
 B_j=2k*(R_j+rho R_-j)².
Aggregate experiment surplus with P shares X,Y is
 S2(X,Y,rho)=2k*[X(Y+rho(1-Y))²+(1-X)((1-Y)+rho Y)²].
In a balanced full-coverage allocation X=Y=.5:
 S2=k*(1+rho)²/2.
Single fully integrated monopolist: S1=2k.
Difference S2-S1=-(k/2)(1-rho)(3+rho) <=0.
Under global resource feasibility, one monopolist enjoys at least the gross matching-investment payoff of fragmented fully covered hubs, unless the two-node technology also changes quality, participation or spillovers.

## Conditional location/transport comparison
For two sides uniform[0,1] with efficient distance assignment and quadratic disutility t, optimal one center=.5, two centers (.25,.75). Mismatch cost saving from one to two = t/8. Under the aforementioned X=Y=.5 **conditional** partition:
 ΔB_2-1(rho)=t/8-(k/2)(1-rho)(3+rho).
A two-hub addition with real extra fixed resource cost F is conditionally welfare beneficial iff
   t/8 - (k/2)(1-rho)(3+rho) > F.
If F>=t/8, no rho∈[0,1] justifies the second hub in this fully covered symmetrical benchmark. If F<=t/8-3k/2, all rho∈[0,1] justify it (assuming F≥0). Otherwise threshold
  rho> -1+sqrt(4 -2(t/8-F)/k).
This is NOT an endogenous location/price/entry equilibrium; exact only for the imposed socially efficient balanced allocation and equal-region types. Stage19's artificially zoned single-hub 0.4 vs dual .25/.75 makes transport saving 29t/200 instead of t/8. With k=1/450 the largest cross compatibility welfare gain at balanced full coverage is 3k/2=1/300=0.003333..., much too small by itself to displace Stage19's ~0.35 private-social wedge. This underscores need for real new participation/quality effects.

## What must be added to explain NEW distant cross-hub innovative projects?
Option A (preferred physical-hub rationale): two geographic clusters each have residents with finite benefit v and travel/experimental access costs. If only one hub opens, remote firms OPT OUT; if both hubs enter, two local populations enroll and cross-hub link rho makes NEW previously inaccessible cross-region matchings. One-hub reference no longer full population M=R=1. This is a substantive change from Stage19's v=4/full cover, and may resurrect price-subgame nonexistence so solvable subgames must be reworked.
Option B (complementary capabilities): cross-cluster partner project gross value A_cross > A_within, but importantly a single facility must be offered all feasible cross-region matches to avoid a tautological special benefit to "two buildings." Cross-hub connectivity effect should come from access and contracting, not arbitrary 2-hub contingent bonus.
Option C (nonappropriable spillover): cross projects create external local knowledge spillover chi>0. Distinguish true social nonappropriable chi from direct private payoffs captured in platform fees.
Endogenize interoperability investment/choice rho via public neutral governance vs private pricing/contracting and compare subsidy/procurement/digital coordination with direct build.

## Simultaneous entry outcome, rigorous strategic complementarity conditions
Let B_G, B_P be gross social surplus under public/private SOLO (potentially asymmetric due geographic participation), B_2(rho) both-active gross social value after valid continuation equilibrium; H0 private monopoly gross fee revenue; H1(rho) private duopoly gross revenue. With respective fixed costs FG and FP,
  W00=0, W10=B_G-FG, W01=B_P-FP, W11=B_2(rho)-FG-FP.
  Private entry payoff without G H0-FP, with G H1(rho)-FP.
Two structures both (0,0) AND (1,1) are STRICT simultaneous-entry pure Nash iff:
  B_G<FG<B_2(rho)-B_P
  H0<FP<H1(rho).
This pair of open cost intervals requires
  B_2(rho)>B_G+B_P,
  H1(rho)>H0.
It is STRONG: public and private entry are BOTH strategic complements; if B_2(rho)-FG-FP>0, no entry is an inefficient coordination equilibrium although both-entry is feasible.
Stage19 baseline **FAILS BOTH**: B2≈7.959 <2*7.818, H1≈.500 <H0≈7.284. Adding rho lightly to the Stage19 full-coverage game will NOT magically make the coordination result hold. The remote-absence/quality mechanism needs quantitative power.
If no strict both-side complementarity, there can still be over-/under provision relative to planner and the effect of rho on entry distortions is analytically interesting. Need solve full downstream price, location, entry and private service incentives, not assume rho shifts total benefit without changing H1.

## Economic mechanism distinction
"Business stealing" is primarily a redistribution of private operating revenue; it is NOT automatically a real welfare cost. Real effects are match breakup/dilution, insufficient pre-match investment, duplicate fixed costs, lost scale, transport and project-quality changes. If links boost competitors' users but hubs cannot charge those users, partial appropriation is an entry externality; show contractual or participation friction microfounded.
Potential publishable question: *When do public and private hubs' entries become strategic complements due to endogenous access/novel cross-region matches, whereas their fees remain strategic rivals?* Derive a non-monotone entry-distortion threshold in geographic remoteness, user complementarity, investment sunkness and market thickness that collapses if local catchment exists but no cross-projects (rho=0) OR cross-projects exist without pre-match effort (c→∞). High absorption risk versus network interconnectivity and interoperability literature; **do not claim originality yet**.

## Mandatory close literature
- Carballa-Smichowski, Duch-Brown, Martens, Gomez-Losada (2026), "Inter-platform ecosystems", Strategic Management Journal 47(7):1840–1877, DOI 10.1002/smj.70070: empirical 241 European platforms, cross-platform matching/co-opetition, four interconnection modes. Cross-hub spillover as such is NOT NEW.
- Zhu, Li, Valavi & Iansiti (2021), "Network Interconnectivity and Entry into Platform Markets", Information Systems Research DOI 10.1287/isre.2021.1010: inter-market buyer mobility, entry profitability and incumbent response, sometimes nonmonotone.
- Kamada & Kojima (2026), "Choice or competition: Does integration benefit everyone?", Theoretical Economics DOI 10.3982/TE6468: integrated matching broadens choice but strengthens competition, monotonic individual gains impossible for well behaved mechanisms.
- Roche, Oettl & Catalini (2022), "(Co-)Working in Close Proximity: Knowledge Spillovers and Social Interactions", NBER wp30120; physical co-location knowledge transfers and dissimilar partners.
- Cespa & Vives (2018), "Exchange Competition, Entry, and Welfare", CESifo: excessive or insufficient entry, fragmentation, shared/liquidity network benefits.
- Bourreau et al. (2023), "Interoperability between Ad-Financed Platforms with Endogenous Multi-Homing", CESifo 10332: private interoperability can deviate from social optimum.

Recommended next step: fully specify 2-region participation (real local nonparticipation), match quality for distinct capabilities, public and private jointly decided site entry, endogenous hub fees (or controlled reduced-form as precursor); prove off-path continuation existence, then derive exact thresholds where network gains overwhelm displacement/duplication and whether public directly owned hub is warranted vs subsidy/contract.
