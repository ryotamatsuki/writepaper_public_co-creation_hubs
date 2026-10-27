# Stage 27 — Occam-parsimonious two-region FOUR-player simultaneous-entry model

Date: 2026-10-27.
Status: RESEARCH DESIGN AND EXACT CONDITIONAL ENTRY-INEQUALITY CHARACTERIZATION; no empirical or solved general SPNE claims. User's immutable premises: urban region A with independent government a and one independent private HUB operator 1, rural B with independent government b and private operator 2. Four players are TERRITORIALLY RESTRICTED and choose local hub entry simultaneously. BOTH public and private may operate in the SAME town. DOWNSTREAM producer firms exist/compete even if no hub.

## Nonnegotiable combinatorics and true simplification
2^4=16 first-stage complete entry profiles cannot be reduced WITHOUT losing action feasibility or imposing a restriction. But the model can be represented by TWO LOCAL 2×2 mixed-ownership games coupled by a single interregional network reach index, instead of a printed 16-row payoff matrix. In each region allowed {0,G,P,GP}. Conditional on other region's reach, local Nash conditions are merely 4 pairs of inequalities. Global Nash are LOCAL BEST RESPONSES CONSISTENT across the regions. Pure NE potentially multiple. Stage23 previous single prefecture, 3-player game and zero-without-hub demand is not the current model.

## Baseline production and product market: no hubs NOT zero activity
Existing n_A>n_B>=1 independently owned manufacturing firms located in respective regions, separate from the 2 private HUB intermediaries. Consumers in town i with inverse market price P_i=a_i-b_iQ_i, n_i symmetric Cournot manufacturers have baseline marginal cost c_i<a_i:
 q_i0=(a_i-c_i)/[b_i(n_i+1)],
 π_firm,i0=(a_i-c_i)^2/[b_i(n_i+1)^2],
 CS_i0=n_i²(a_i-c_i)²/[2b_i(n_i+1)²].
Business-as-usual R&D / private collaboration may also be added as θ0, but endogenous R&D effort omitted from CORE. Urban scale can be greater through a_i,b_i/n_i; urban scale and competition intensity n_i distinct.
If interregional PRODUCT trade desired, two-market Cournot also available; do not confuse territorial HUB ban with export ban.

## Minimally credible G-free vs P-priced co-presence
For local representative potential user of type θ~U[0,1], incremental utility OVER its no-hub production payoff:
 U0=0, UG=v_i>=0, UP=v_i+θΔ_i−f_i.
Private advantage θΔ_i is sector specialism, network quality, service bundle; must not just be unexplained automatic superiority in empirical interpretation.
If both public and private build and private has per-user variable cost c_i^hub∈[0,Δ_i):
 private participation N_P=N_i(1−f_i/Δ_i), public rest, interior f∈[0,Δ_i].
 Private gross contribution max_f (f−c_i^hub) N_P occurs at
 f*=(Δ_i+c_i^hub)/2,
 Π_P,operating=N_i(Δ_i−c_i^hub)²/(4Δ_i).
Hence coexistence IS FEASIBLE despite zero public fee iff F_Pi<N_i(Δ_i−c_i^hub)²/(4Δ_i), abstracting network effects and assuming G is already in. This is NOT by itself full simultaneous-entry SPNE: govt G decision and monopoly P-profit when no G must also be checked. If Δ_i=0 and c_i^hub>=0, profitable coexistence at positive private fee is impossible under identical free public service.
When P only, optimal fee differs and user participation may be partial; do not substitute coexistence formula into P-only history.

## Compact reduced-form local-entry threshold model
Let g_i,p_i∈{0,1}; define effective local interregional contact surface s_i=g_i+p_i−δ g_i p_i, 0<=δ<=1 (δ captures partial duplicated network contact when two hubs in SAME region), s_j similarly. Local gross welfare, net of real fixed costs, above Cournot no-hub benchmark W_i0:
  W_i = W_i0 + B_i^G g_i + B_i^P p_i −D_i g_i p_i
        +ρ Ω_i s_i s_j −FG_i g_i −FP_i p_i.
Here B^G,B^P are incremental real benefits from local hub services, innovation and producer competition net of endogenous consumer/producer transfers; D_i>=0 is local duplication/real competition inefficiency, distinct from private business stealing; Ω_i captures external benefits in region i arising from interregional links and may differ by urban/rural density. Private intermediary net entry payoff if builds (zero if not):
  π_Pi = p_i [R_i−χ_i g_i +ρ K_i(1−δ g_i)s_j−FP_i].
R_i is endogenous-price optimized operating receipt in P-only market, χ_i>=0 captures public free service's rent displacement, K_i represents private appropriable cross-network revenue. These reduced-form primitives MUST be derived from explicit demand/fee/product-R&D subgame for actual publishable proof; used here only as analytical shortcut and proposition finder. R_i and K_i not manufacturing profit.

Let z_i=s_j, a conditional observed scalar. Government's unilateral entry gain:
  ΔGi(p_i;z_i)=B_i^G−FG_i−D_i p_i+ρ Ω_i(1−δ p_i) z_i.
Private net gain:
  ΔPi(g_i;z_i)=R_i−FP_i−χ_i g_i+ρ K_i(1−δ g_i) z_i.
Local 2×2 conditions (weak inequalities permit ties):
 NONE: ΔGi(0;z)<=0 & ΔPi(0;z)<=0.
 G ONLY: ΔGi(0;z)>=0 & ΔPi(1;z)<=0.
 P ONLY: ΔGi(1;z)<=0 & ΔPi(0;z)>=0.
 BOTH: ΔGi(1;z)>=0 & ΔPi(1;z)>=0.
Globally require both local conditions simultaneously, with z_A=s_B and z_B=s_A. Does not drop any of the 16 original profiles. No sequence assumption.
Potential public/private coexistence in region i iff both last inequalities, requiring net G benefit after P and P residual rent after G. Public zero fee alone doesn't rule this out if specialized private value/ability to capture fee positive.
Social planner maximizes W_A+W_B across the 16 profiles AFTER downstream behavior; local government maximizes W_i only. Cross-region benefits cause marginal planner benefit of site's increase greater than local own-government benefit, e.g. rho Omega_j s_j times local Δs_i; standard Cremer et al. 1997 underprovision mechanism. No journal novelty from this threshold itself.

## Strict caveat around functional form
Use s_i as contact OPPORTUNITY capacity, not number of uniquely matched participants when both local hubs serve overlapping market, else double count. More honest structural derivation: seller/user participation mass by hub + compatible cross-hub link set. Do not call hub count automatically new unique contacts. D_i nonnegative only duplication real costs; neither D_i nor χ_i directly equals pure transfer.
Reduced-form formula may be used to SEARCH parameter candidates, not as new fundamental general structural theory.

## Publication novelty kill tests
Cremer–Marchand–Pestieau (1997) already 2 independent municipal indivisible public facilities plus cross-jurisdictional spillovers.
Long–Wong (2009) already 2 local governments/private firms and public infrastructure/price competition.
Kim (2024) mixed public/private two-sided markets price/quality, high overlap.
d'Aspremont–Jacquemin (1988) baseline innovation and Cournot; review 1990 erratum when using old propositions.
Zhu et al. (2021) cross-platform market connectivity and entry.
Thus merely proving over/under/multiple equilibrium or ownership PPP is insufficient.
The strongest possible central proposition: no-hub productive baseline + genuinely ENDOGENOUS client choice/fees + geographical asymmetry and interregional nonappropriable links causes ownership COUNT-correct but TYPE-wrong equilibrium, with sharp robust threshold NOT absorbed by nearest parent. Need public-broker subsidy/contract equivalence to claim direct construction necessity.

## Suggested nested minimal models
- Core M0: independent regional Cournot no-hub benchmark plus binary local authority/private building; ONE network interregion parameter rho; SINGLE paid private participation fee (the public free), TWO local free/paid differentiated choice; fixed local site, no match-specific hold-up, no dynamic R&D or endogenous access investment. Derive 16-profile payoff continuations but present via four local best-response inequalities.
- Extension M1: endogenous bilateral network-link investment, test market creation vs private cannibalization; possibly R&D investment before product-market Cournot.
- Extension M2: two-sided I/R participant populations/fees only IF M0/M1 paper theorem needs it.
- Extension M3: municipal subsidy, service contract, neutral data platforms vs physical public operation.
No move to sequential public leader or exogenous exclusive ownership unless explicitly acknowledged that the research question changed.

## Literature source pages
https://www.sciencedirect.com/science/article/abs/pii/S0047272797000042
https://onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9396.2009.00822.x
https://onlinelibrary.wiley.com/doi/10.1111/manc.12467
https://www.sciencedirect.com/science/article/pii/S0899825602005626
https://ideas.repec.org/a/aea/aecrev/v78y1988i5p1133-37.html
https://www.di.ens.fr/~aspremon/Claude/PDFs/dAsp90b.pdf
https://doi.org/10.1287/isre.2021.1010
