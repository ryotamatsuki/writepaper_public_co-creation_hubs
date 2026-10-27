# Stage 16 — Identify an economic PUBLIC-ESTABLISHMENT proposition first (prior art / algebra / kill tests)
Date: 2026-10-27
Status: **NEW MODEL HYPOTHESIS / EXACT CONDITIONAL ALGEBRA PROVED / NOVELTY NOT CERTIFIED / NOT FULL STAGE12 SPNE**.

## Why change mechanism
Stage 15's off-path pure-SPNE nonexistence is a technical equilibrium-property result, NOT a portable rationale for municipal installation. Existing literature already studies mixed-duopoly location nonexistence, endogenous matching technology, optimal restrictions on partner search, directed matching and platform interoperability. We instead target **a credible PUBLIC FALLBACK as an unused bargaining outside option**, and ask whether it complements private entry even when the government facility has ZERO realized routine demand. Do not call this novel before closest contract-theory/public-option literature is audited.

## Explicit two-institution minimal candidate game
- There are L=μMR potentially successful private-facilitated experimental relationships if the incumbent private matching firm P enters. The scalar L is ENDOGENOUS to masses M,R and matching efficiency μ; this stripped-down first lemma does NOT endogenize initial platform sorting, the service location, or partner-scope n=1 vs 2.
- G chooses whether to establish a neutral alternative experimental site at real fixed cost FG. Its potential ability to support alternative experimentation upon bargaining breakdown provides an outside-option value z A e for the innovating firm, z∈[0,1]. z=η q_G, with q_G a feasibility probability that depends on available public reserve partner thickness and η investment redeployability. There is **no actual use of G on equilibrium path when P operates**: G is an option invoked only when the P match breaks down.
- Private P observes public establishment, chooses to enter and pay real FP. If P operates, L private projects occur. Innovator invests e into a project yielding gross Ae; noncontractible e costs c e²/2. Ex-post generalized Nash bargaining on the expected disagreement payoffs (innovator zAe, partner zero) gives innovator fraction phi(z)=beta+(1-beta)z and partner fraction (1-beta)(1-z). beta in (0,1).
- Fee/revenue microfoundation: P can charge the matched partner a preannounced ad valorem commission sigma∈(0,1) on her realized ex-post bargaining receipt. Commission is a TRANSFER and leaves the Nash division of gross project surplus unchanged (equivalently, it is paid by the partner out of the receipt, and after-fee participation is nonnegative). P does NOT choose sigma in this first lemma; allowing negotiated/dynamic commissions is mandatory next robustness.
- If G is the only facility and P stays out, G can directly facilitate L matches with the same A,c, matching efficiency, and beta but no independent second outside option, z=0. Total gross welfare is thus the same as private-only benchmark before fixed costs. A spare public test/partner capacity makes the alternative relationship credible, even if not taken up on path; the fixed public cost finances maintaining that capability.
- Outside no facility, no mediated relationship (baseline social welfare 0).
- No transfer of public tax revenue into private projects and no output subsidies. Funding has shadow resource cost 1.

## Exact post-investment bargains and efforts
Let C=L A²/c, phi=beta+(1-beta)z.
Ex-post firm payoff before effort is phi A e; investor optimal e(z)=phi A/c.
Social project value per match:
  Ae(z) - c e(z)²/2 = (A²/c)[phi - phi²/2].
Aggregate GROSS welfare if P operates:
  B(z)=C[phi-phi²/2].
Private P gross fee operating revenue:
  P(z)=sigma C(1-beta)(1-z)phi.
With no G, phi=beta, B0=C[beta-beta²/2], P0=sigma C (1-beta)beta.
With G, Bz=B(z), Pz=P(z).

**Proposition 16A (public outside option / private-entry COMPLEMENTARITY):**
  P(z)-P(0)
  =sigma C(1-beta)z[(1-2beta)-(1-beta)z].
Therefore public G raises PRIVATE gross entry revenue iff
  beta<1/2 AND 0<z<(1-2beta)/(1-beta).
The revenue curve is inverted U: argmax z*=(1-2beta)/(2(1-beta)) when beta<1/2.
Mechanism = investor's bigger effort makes the privately matched partner's bargaining rent larger at first, more than offsetting rent dilution from the investor's outside option; eventually rent dilution dominates. It is NOT direct subsidization, public price-cutting, or higher on-path match counts.
Private entry occurs with public iff FP<=Pz, without public iff FP<=P0. Thus public crowds in P entry iff P0<FP<=Pz. For constant z, this is an INTERMEDIATE matched-market-thickness window:
 FP/[sigma (A²/c)(1-beta)(1-z)phi] <= L
   < FP/[sigma (A²/c)(1-beta)beta].
This demonstrates public co-creation facility can crowd in a private intermediary even if neither facility is used jointly on path.

**Proposition 16B (exact public-installation cutoff conditional on a pre-specified private-entry continuation):**
  W(no G) = 1{FP<=P0}(B0-FP).
  W(G)    = 1{FP<=Pz}(Bz-FP-FG)
           +1{FP>Pz}(B0-FG).
Municipality builds iff W(G)>=W(no G). This generates FOUR regimes:
  I  FP<=min(P0,Pz):      build iff FG<=Bz-B0
       = C(1-beta)² (2z-z²)/2.
  II P0<FP<=Pz:         public CROWDS IN private entry; build iff FG<=Bz-FP.
  III Pz<FP<=P0:        public DETERS private entry; build iff FG<=FP.
  IV FP>max(P0,Pz):      neither has private entry; build iff FG<=B0.
At exact thresholds, entrant enters when indifferent.
Notice public installation can be optimal even if **G is NEVER USED ON THE EQUILIBRIUM PATH** when P operates. This is option value realized by higher sunk investment.
If P enters in both regimes, the threshold is a quadratic in z and MATCHING MASS L: FG< L A²/c (1-beta)²(2z-z²)/2.
In the crowd-in case (II), the threshold adds an extensive private entry effect and can exceed the value of G-only supply.

**Proposition 16C (strict joint viability — a diagnostic witness):**
beta=1/10, z=1/2, sigma=4/5, C=L A²/c=1,
 FP=13/100, FG=1/10.
Then phi=11/20,
 B0=19/200=.095,
 Bz=319/800=.39875,
 P0=9/125=.072,
 Pz=99/500=.198.
Hence without G the private firm does NOT enter (FP>P0), while with G it DOES enter (FP<Pz).
Public-only welfare B0-FG=-1/200<0, and private-only entry is unprofitable since P0<FP. Yet joint G+P welfare Bz-FP-FG=27/160=.16875>0.
=> Neither G-only nor P-only is privately/economically viable in its subgame, while the municipality BUILDING G and inducing P to enter is socially desirable. This is strict (holds in an open neighborhood).
Reproduction of rational witness performed independently in container at /mnt/data/stage16_outside_option_audit.py.

## Critical closest-neighbor discovery: IDLE public fallback is ALREADY KNOWN
Martinez-Giralt & Pita Barros (2005), "Bargaining and idle public sector capacity in health care", Economics Bulletin 9(5):1–8, **explicitly rationalize IDLE publicly owned capacity via the value of the outside option during negotiations with a privately efficient provider**.
https://ideas.repec.org/a/ebl/ecbull/eb-05i10002.html
Therefore "idle public facility has bargaining value without on-path use" is **NOT NEW, and must NOT be presented as a novel standalone economic mechanism**.
The candidate possible additional contribution is only the JOINT result:
- Noncontractible pre-match experimentation effort endogenously responds to public outside option;
- private intermediary derives endogenous partner-side commission profit that is inverted-U in fallback quality;
- public infrastructure **induces** private intermediary entry in an intermediate matched-market-thickness band, with neither public-only nor private-only viable but joint configuration viable;
- an implementability boundary separates real public provision from contractible private standby provision.
These elements must ALL survive direct comparison with prior art and a less contrived fee/entry model. Strong ABSORPTION RISK persists. Also Wilson (2025) already models municipal infrastructure, private investment crowd out/in and public entry, so simple private-entry crowd-in is NOT NEW in general.
The exact paper-specific contribution is NOT yet established.

## Critical no-go and absorption tests
1. **NOT NEW BY ASSERTION**: outside options and relation-specific investment are classic Grossman-Hart-Moore/hold-up logic; Goldlücke & Schmitz (2014 JET) directly study investments that affect outside options; Hwang & Li (2017 JET) study transparency and bargaining outside options; Bernhardt/Hwang/Krasa (2025 working paper) study joint and outside-option investments; Wilson (2025 AEJ Micro) shows municipal public competition and private investment incentives.
2. **Prove that public BUILDING, not public financing, is needed**: If a contracting authority can purchase the exact same fallback availability/quality from an existing private laboratory at resource cost <= FG, with enforceable open access and observed match quality, then own-construction is never uniquely desirable. This is an immediate implementation-equivalence / no-go theorem. Need formal unverifiability/commitment friction or a second-mover credible-neutrality advantage to justify ownership. Shleifer (1998 JEP) is especially relevant.
3. **Contractual credibility**: Post-effort public alternative zAe is a stipulated credible disagreement payoff. Need microfound from an equilibrium of public reserve-partner matching (and posterior information). Must check whether the Nash outside-option bargaining convention is appropriate if public option is realized probabilistically after breakdown.
4. **Private platform revenue**: P receives fixed ex-ante commission share sigma of partner receipts. Allow P to choose fees/commission and compare its revenue maximization with public fallback. If P re-optimizes ad valorem charges, crowd-in may disappear. Verify rather than assert.
5. **Double capacity and option**: public reserve match capacity can deliver fallback if P bargaining fails; when G alone it must be capable of matching the same L projects (or explicitly account for unequal capacities). Zero on-path use is consistent with credible off-path availability, but public-maintenance fixed cost FG must pay for reserve readiness.
6. **Platform choice/location**: This is a small stand-alone THEOREM TARGET, not yet nested in Stage 12's endogenous Hotelling locations and participant fixed points. Do not copy Stage 13's cubic. Require robustness with endogenous location or prove option mechanism after adding it.
7. **Alternative sources**: partner switching to a private backup, grants contingent on investment, private consortium, subsidy, tendering, public digital clearinghouse. The distinctive physical site must provide additional experimentally verified/noncontractible fallback quality.
8. **Public own outside options**: For public-only mode with no P assumed beta share, independent public alternative unavailable. If an alternative can be secured without public building, the public option may add zero incremental benefit.
9. **Prior art matching scope**: Matsushima & Shinohara (2014) partner numbers/quality; Aoyagi & Yoo (2022 GEB) optimal private matching mechanism; Kanoria & Saban (2021 Management Science) restricting contact; Teh, Wang & Watanabe (2024 JET) optimal restricted accessibility; Economides & Katsamakas (2003 IEP) two-sided platform compatibility; Kim (2025 Theoretical Economics) priority search entry intervention; Kamali Shahdadi & Mozafari-Vanani (2026 J Mathematical Econ) integrated unbalanced matching; McLaren (2003 JEEA) thickness and specific investments; Schmitz & Hoppe (2008/2019) public/private ownership under incomplete contracts.
10. **Novelty status**: Candidate exactly ties PUBLIC UNUSED FALLBACK with private operating entry complementarity and a four-branch installation threshold, but each constituent mechanism exists. No statement that neighbors have failed to derive this exact theorem without full-paper comparison; current status CONDITIONAL NEW-TO-THIS-TRACK, not new-to-literature proven.

## Preferred research question/title
"Public Innovation Hubs as Credible Outside Options: Matching, Relationship-Specific Investment, and Private Platform Entry."

**Pass to mathematical construction stage only if** absent two-hop public backup, a public hub can commit to maintain a neutral alternative relationship in a way that private intermediation/subsidy cannot; exact four-region build rule survives endogenous user sorting and chosen price/commission, and nearby full models do not already subsume the result.

Source links:
https://www.sciencedirect.com/science/article/pii/S0022053113002056
https://www.sciencedirect.com/science/article/pii/S0022053116301016
https://www.aeaweb.org/articles?id=10.1257/mic.20210110
https://pubs.aeaweb.org/doi/10.1257/jep.12.4.133
https://www.sciencedirect.com/science/article/pii/S0899825622001063
https://research.monash.edu/en/publications/strategic-limitation-of-market-accessibility-search-platform-desi/
