# Stage 17 — Public co-creation hubs WITH existing private incumbents: distinguishing facility undersupply from function undersupply
Date 2026-10-27. Status: RESEARCH DESIGN / NECESSARY-SUFFICIENT THRESHOLDS FOR AN ILLUSTRATIVE BRIDGING PROJECT / NO NOVELTY CLAIM / NO FULL GAME SPNE.

## Research question
Why build a publicly operated co-creation hub **when private co-creation hubs already exist**? Can the private sector underprovide hub NUMBER, or only some socially valuable FUNCTIONS? Why public direct provision rather than subsidy, procurement, private coordination or purely digital matching?

## Critical verified closest sources
- Mankiw–Whinston 1986 RAND "Free Entry and Social Inefficiency": entry may be inefficient; business stealing and product diversity.
 https://mankiw.scholars.harvard.edu/sites/g/files/omnuum5931/files/mankiw/files/free_entry.pdf
- Basak & Petrakis 2021 JEMS "Social efficiency of entry: Implications of network externalities": free entry insufficient OR excessive; entry differentiation, network effects, compatibility threshold.
 https://doi.org/10.1111/jems.12431
- Cespa & Vives 2018 / revised 2021 CEPR "Exchange Competition, Entry, and Welfare": platform entry excessive OR insufficient; fragmentation and underinvestment.
 https://cepr.org/publications/dp13415
- OECD (June 25, 2026) "Incubation in Entrepreneurial Ecosystems", ch.3 and ch.11, explicitly say incubation is undersupplied owing to innovation spillovers, startup finance/info, service information, ecosystem connectivity/network failures, including outside major hubs, and that direct public provision is NOT necessarily superior to funding private providers.
 https://www.oecd.org/en/publications/incubation-in-entrepreneurial-ecosystems_21be9d03-en/full-report/why-does-policy-support-incubation_7d9c307d.html
 https://www.oecd.org/en/publications/incubation-in-entrepreneurial-ecosystems_21be9d03-en/full-report/key-policy-choices-in-supporting-business-incubation_023e1a23.html
- OECD (Oct 20, 2026) "Data Spaces between Co-ordination Challenges, Financial Viability and Public Policy", specifically identifies MULTI-PARTY neutral collective governance, interoperability and non-excludable external value as implementation problems; therefore simply framing neutral interoperability as new is absorbed.
 https://www.oecd.org/en/publications/data-spaces-between-co-ordination-challenges-financial-viability-and-public-policy_0066f11b-en.html
- Kamada & Kojima 2026, "Choice or competition: Does integration benefit everyone?", Theoretical Economics: merging fragmented matching markets can harm some participants due to competition; do not assume networking helps everyone.
 https://doi.org/10.3982/TE6468
- Bourreau et al. 2023 CESifo, "Interoperability between Ad-Financed Platforms with Endogenous Multi-Homing": private interoperability socially too low OR too high; mandated connection not always desirable.
 https://www.ifo.de/en/cesifo/publications/2023/working-paper/interoperability-between-ad-financed-platforms-endogenous-multi
- Murillo-Rojas et al. 2026, Strategic Entrepreneurship Journal, meta-analysis of 39 studies, >55,000 firms; bridging support stronger than space subsidies on average, public incumbents do not systematically dominate private.
 https://doi.org/10.1002/sej.70032
- OECD 1999 "Business Incubation": warns even if private services undersupplied, DIRECT government funding or supply need not be most efficient vs regulation/subsidies.
 https://www.oecd.org/content/dam/oecd/en/publications/reports/1999/09/business-incubation_g1ghg2c2/9789264173781-en.pdf

## First mathematical audit of OLD Stage 11 mechanism
For a previously integrated one-hub facility handling innovator/partner masses M,R, total bilinear links L0=mu*M*R. If G enters and steals fractions x,y of the two sides, all within-hub links sum to
 L1=mu*M*R*[x*y+(1-x)*(1-y)].
So L1-L0=-mu*M*R*[x*(1-y)+(1-x)*y] <=0; a symmetric half-half split loses one half of bilateral links. Public installation as a purely duplicative second co-creation facility may HARM a thin matching market.

For the restricted earlier Stage 11 policy game with PUBLIC INCUMBENT and private challenger, gross welfare gain from private challenger = T/16-Z/2, and entrant gross operating profits = T/4-Z/2. Hence entrant profit minus social gain=3T/16>0. This is an OVER-ENTRY force, not underentry in the original restricted benchmark. NOTE carefully: this is NOT an unrestricted private-vs-private free-entry theorem, since the incumbent is public and sets welfare-maximizing tariffs.

## Minimal new mechanism class: cross-private-platform public connector, NOT another segregating competitor
Suppose existing private hubs h=1,2; their enrolled innovator masses M_h and partner masses R_h are positive. Internal match mass L0=mu*(M_1 R_1+M_2 R_2).
Define potential cross-platform pair mass Omega=M_1 R_2 + M_2 R_1. A publicly operated cross-hub test/connector can make a fraction lambda of those pairs feasible at match efficiency mu; INCREMENTAL match mass Q=lambda*mu*Omega, capacity-constrained within these types if necessary.
Each incremental success yields real surplus S*>0 after optimal relationship-specific effort and ex-post bargaining, plus public knowledge/spillover benefit chi≥0 counted ONCE and not captured by private hub fees. Public cross-hub matching may cannibalize existing profitable matches; let D_real≥0 denote real resulting net displacement cost after opportunity costs, not foregone transfer revenue. Public fixed resource cost F>0.
Social net addition:
 deltaW_G=Q*(S*+chi)-D_real-F.
Private FOR-PROFIT BRIDGE with identical technology and cost F can privately capture no more than theta*Q*S* in operating fees, theta∈[0,1] — theta MUST later be microfounded from bargaining/fee competition; not an assumed intrinsic private inferiority.
Candidate functionality UNDERPROVISION (not necessarily number-of-hubs underprovision) occurs iff:
   theta*Q*S* < F < Q*(S*+chi)-D_real.
The region is nonempty if Q*((1-theta)*S*+chi)>D_real.
If D_real=0 this is an intermediate-thickness condition:
   F/(S*+chi) < Q < F/(theta*S*)
(for theta>0).
Meaning: very thin markets justify NEITHER platform, intermediate markets justify public bridge but not private bridge, sufficiently thick markets can privately support bridge. This is a CONDITIONAL simple algebraic example, NOT an SPNE of endogenous platform/price game and NOT novel versus general public-good/network underprovision.
Example M1=M2=R1=R2=1, mu=.1,lambda=1 -> Q=.2; S*=3/8, chi=1/8, theta=2/5,F=.06,D_real=0:
 private max fee revenue=.03<.06; social value=.10>.06, public net gain=.04.
This verifies a STRICT interval of functional underprovision even WITH TWO incumbent privately supplied hubs.

**Critical implementation-equivalence no-go lemma**:
If a municipality can procure exactly the same connector service, credible equitable access, matching quality, post-match effort and spillover from a private provider for total real resource cost F_private<=F_public (with transfers correctly accounted for and no administrative distortions), direct public establishment is never strictly more efficient on this dimension. Funding a private operator can realize the efficient connection even when it would not enter voluntarily. Government-specific control, credible neutral access, ex-post noncontractibility, public data/test assets or verifiable public-purpose effort must be ENDOGENOUS DIFFERENCES, not normative assumptions, to establish direct-ownership optimality.

## Research propositions and absorption kill tests
Target H1: Endogenously chosen private hubs can be TOO NUMEROUS yet invest TOO LITTLE in CROSS-HUB MATCHING; sign reversal of physical count distortion vs functional connecting distortions. Must microfound θ from revenue pricing and dynamic partner relationship-specific effort; test Aoyagi/Yoo 2022, Kanoria/Saban 2021, network interoperability papers.
Target H2: Public facilitator optimally COMPLEMENTS existing private hubs over an intermediate THICKNESS interval (not monotonically declining/increasing with scarcity); extensive participation and platform investment must be endogenous to claim originality. For Q→0 and F fixed >0, NO public build — refutes "thin regions always need a hub".
Target H3: Neutral public physical space beats subsidized private consortium only under explicit renegotiation/quality-verification or commitment constraint unavailable to private service contracts. Must compare digital/traditional brokerage as separate feasible instruments. Examining local government's privileged ability to offer administrative data/testbed access is promising but needs no-go and non-absorbed mechanism.
Target H4: Public bridge's improved matchability increases private hubs' own service profits/entry, rather than cannibalizing incumbents, AND generates higher welfare; sign threshold for crowd-in vs crowd-out as a separate proposition.

Ranking: best research goal is NOT generic numerical undersupply. Best is boundary between "another rival publicly owned building" (could fragment matching) vs "neutral public common testbed connecting incumbent private platforms" (could expand net feasible matches), followed by an OWNERSHIP/INSTRUMENT IMPLEMENTATION theorem.

ABSORPTION RISKS HIGH: The general mechanisms and some structural conditions are prior art, especially OECD2026 reports, network compatibility and integration models, public innovation intermediaries literature, and market failure public ownership. No novelty is certified at Stage17. Need full peer-reviewed-paper comparison and complete equilibrium derivation before journal-grade designation.
