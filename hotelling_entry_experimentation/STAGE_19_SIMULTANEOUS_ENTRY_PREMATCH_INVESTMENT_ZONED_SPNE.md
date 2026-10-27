# Stage 19 — Constructive full SPNE of a **ZONED** simultaneous-public/private-entry, pre-match-investment, bilinear-matching game

Date: 2026-10-27
Verdict: **ANALYTICALLY VERIFIED GLOBAL PRICE NASH FOR EVERY FEASIBLE ZONED LOCATION AND EXACT-BERNSTEIN-CERTIFIED UNIQUE ZONED LOCATION NASH; FOUR-STRUCTURE PAYOFFS AND ALL PURE ENTRY NASH CONDITIONS.** NOT the unrestricted-location game of Stages 12–15; novel economic contribution NOT established.

## Model, important scope restriction
One mass of innovators x~Uniform[0,1] and partners y~Uniform[0,1]. No pre-existing hub (N0=0). Government G and private P SIMULTANEOUSLY choose whether to pay FG,FP >=0 at Stage 1; Stage 2 each surviving provider selects its service location SIMULTANEOUSLY: G a∈[0,2/5], P b∈[3/5,1]. **These site restrictions apply EVEN when only one provider enters**, so G-only optimally chooses a=2/5 and P-only chooses b=3/5 (NOT 1/2). Zoning/segmented service licenses are an EXPLICIT ADDITIONAL ASSUMPTION to avoid collocation-related pure-price-Nash nonexistence at off-path pairs. Both sites remain ENDOGENOUS within their entire intervals. No claim for unrestricted a,b∈[0,1].

Stage 3: each provider posts separate nonnegative enrollment fees p_j^I,p_j^R SIMULTANEOUSLY. Stage 4: each population chooses hub or outside (utility zero). Stage 5: actual matching links L_j=μ N_j^I N_j^R, 0<μ<=1, each participant matches at most once. Stage 6 (BEFORE match uncertainty resolves) innovators choose noncontractible effort e_j>=0 at per-innovator cost c e_j²/2. Stage 7: if match succeeds Nash bargain shares β to innovator, 1−β to partner of gross project return A e_j. Realized transfers cancel from W. Parameter β=2/3 symmetrizes the investor's NET expected benefit and the partner's expected receipts conditional on equal side masses; this is an analytically meaningful but restrictive choice. v=4,t=1,A=c=1,μ=1/10 for the exact illustrative SPNE. Fees transfer within W; no public-funds excess burden.

## Match/effort subgames
Innovator's success probability q_j=μ N_j^R, investor optimizes q_j β A e−c e²/2, so
  e_j*=β A μ N_j^R/c.
Define k=β² A² μ²/(2c)=2 A² μ²/(9c). Under β=2/3, partner expected remuneration per partner equals k N_j^I N_j^R, innovator's ex ante net surplus is k (N_j^R)^2. Hence aggregate REAL experiment surplus in j:
  S_j=2k N_j^I(N_j^R)^2.
With both hubs and P-side population fractions X,Y, aggregate experiment value
  S(X,Y)=2k[X Y²+(1−X)(1−Y)²].
For specified μ=.1,A=c=1, k=1/450.

## Demand and globally valid 2-sided price Nash for ALL zone location pairs
Fix a∈[0,.4], b∈[.6,1]; set d=b−a>=.2, h=t d, m=(a+b)/2, C=(1−a−b)/2.
For full coverage, given public two-side fees (gI,gR), private (pI,pR),
  2h(X−(1−m))=gI−pI+k(2Y−1),
  2h(Y−(1−m))=gR−pR+k(X+Y−1).
For arbitrary price vectors and opt-out zero, whenever any ONE rival's fees are held at prescribed equilibrium values, that rival's fallback utility is strictly positive for every type, so no user opts out under unilateral price deviations. The complete full-coverage shares solve the simultaneous clipped fixed-point equations corresponding to above. Their fixed-point map is a contraction, because d>=1/5 and k=1/450 imply max row-sum Jacobian <=k/h<=1/90<1; UNIQUE subscription equilibrium at all fees.

Both-active social surplus, before entry fixed costs:
  W0(a,b;X,Y)=2v − T(a,b;X,Y)+2k[XY²+(1−X)(1−Y)²],
with T the exact sum over the two sides' quadratic losses assigning the [0,1−X] innovators to a, [1−X,1] to b and similarly for partners Y.
Public can implement the socially optimal interior shares by its two side fees; global W0 strict concavity in (X,Y) follows from
  Hessian=[[-2h,4k],[4k,−2h+4k]],
negative definite for h>k(1+sqrt(5)).

The unique efficient subscriber fractions are explicit RATIONAL functions of a,b:
 D=h(h−2k)−4k² >0,
 X=1/2 + h² C/D,
 Y=1/2 + h(h+2k) C/D.
Private's inverse-fee revenue envelope Q over ANY shares X,Y∈[0,1]:
 Q=(gI+2h(1−m)−k)X+(gR+2h(1−m)−k)Y
       −2h(X²+Y²)+3kXY+kY².
The Q Hessian=[[-4h,3k],[3k,−4h+2k]] is strictly negative definite for h>=1/5,k=1/450. This is an upper bound when demand hits corners and exact when interior, hence controls ALL non-interior fee changes. The unique private fee BEST RESPONSE associated to the public efficient shares is:
  pI=2hX−kY,
  pR=−2kX+(2h−k)Y.
To implement W0-efficient shares the public charges
  gI=pI+k(2Y−1),
  gR=pR+3k(X+Y−1).
Then public has no fee deviation (global social optimum), and private no fee deviation (Q global strict maximum). Thus these rates constitute a fee Nash profile at EVERY feasible a,b, with unique demand. There may be other fee Nash profiles, not claimed unique.

All expected equilibrium enrollment fees remain strictly positive and <5/2, all X,Y strictly in (0,1), for every a,b in the zoned rectangle. These facts, plus D>0, were independently certified with EXACT RATIONAL Bernstein-polynomial bounds over a∈[0,2/5], b∈[3/5,1], for polynomial numerators D, X D, (1−X)D, YD,(1−Y)D,pI D,pR D,gI D,gR D,(5/2−gI)D,(5/2−gR)D. The bounded fees guarantee all outside options irrelevant under any unilateral fee deviations (baseline v=4 and max mismatch t=1).

Public's objective at location stage is the real gross W0 after the above price Nash, while private maximizes Π(a,b)=X pI+Y pR=2h(X²+Y²)−3k XY−k Y².

## Global location Nash: exactly verified uniqueness in zones
Define algebraic derivatives W_a=dW0(a,b)/da and Π_b=dΠ(a,b)/db after replacing the subscriber shares by the efficient rational X,Y; use social envelope for W_a:
  W_a=(1−X)²+(1−Y)²−2a(2−X−Y).
For P: differentiate Π's rational formula w.r.t. b (the derivative includes the fee/share re-optimization for each new b).

An exact rational **bivariate Bernstein certificate** for k=1/450 on whole rectangle verifies:
  W_aa<0 and Π_bb<0 everywhere,
  W_a(0,b)>0, W_a(2/5,b)<0 for all b∈[3/5,1],
  Π_b(a,3/5)>0, Π_b(a,1)<0 for all a∈[0,2/5].
Thus each provider has a UNIQUE **GLOBAL** location best reply interior to its zoned feasible set.

Further exact rational Bernstein signs show the best-reply correspondences lie first in
  a_BR ∈[.19,.35], b_BR∈[.65,.85],
and after reapplication within
  a_BR∈[.21,.29], b_BR∈[.72,.80].
On the latter rectangle the cross partials W_ab and Π_ba are strictly positive and the Bernstein sign inequalities imply
  0< W_ab/(−W_aa)<1/2,
  0< Π_ba/(−Π_bb)<1/2.
Hence successive scalar best replies compose into a global contraction with modulus less than 1/4 on the invariant rectangle. All simultaneous location Nash equilibria must lie there, so the equilibrium is UNIQUE. This argument certifies GLOBAL location best responses over WHOLE permitted location sets; not merely the FOCs or finite grids.

Unique solution of W_a=Π_b=0:
  a*=0.2488406355..., b*=0.7465841115...
  X*=0.50230842287..., Y*=0.50232903521...
  G fees (I,R)=(0.49893554535,0.49674414791)
  P fees (I,R)=(0.49892519408,0.49671323152)
  G operating fee revenue≈0.49553115792 (TRANSFER, not G's objective)
  P operating fee revenue Π_D≈0.50012780574
  P matches μXY≈0.02523241054; G matches μ(1−X)(1−Y)≈0.02476866474
  P innovator effort=(2/3) μY≈0.03348860235; G effort=(2/3) μ(1−Y)≈0.03317806432
  GROSS joint real W_D≈7.95943671247.

## Single-provider and no-provider histories — separately solved
No hub: W00=0.
G-only G chooses nearest-to-central feasible a=2/5, fee zero (full coverage). The exact gross social use welfare:
  N=2v−2[t/12+t(.1)²]+2k=1759/225≈7.8177777778.
P-only P chooses nearest-to-central feasible b=3/5, sets each side's full-coverage monopoly fee
  p_m=v−t(3/5)²+k=4−9/25+1/450=1639/450≈3.6422222222.
Private P-only GROSS operating fee revenue:
  H=2p_m=1639/225≈7.2844444444.
Private monopoly full-coverage is GLOBALLY optimal, not just local: for any privately feasible location b>=.6 and enrollment share z, the radius R_b(z) of the z nearest customers obeys R_b(z)∈[0,1], weak derivative R_b'(z)≤1 where defined; the marginal revenue of increasing market coverage on either side is ≥v−(R_b(z))²−2zR_b(z)R_b'(z)≥v−3>0; network match gains are nonnegative. b=.6 first-order stochastically dominates other b>=.6 in distances to all consumers, so the center-most P site is optimal. G-only social optimum is a=.4.

For entry (g,p), real societal payoffs:
 W00=0
 W10=N−FG
 W01=N−FP
 W11=W_D−FG−FP.
Private gross revenue if P enters: H if g=0, Π_D if g=1; private NET entry payoffs H−FP and Π_D−FP. G's pricing receipts are transfers already counted in W.
Define C_D=W_D−N≈0.14165893469.
Numerical hierarchy: 0<C_D<Π_D<H<N.

## Simultaneous stage1 binary-game equilibria (allow weak inequalities and multiplicity)
 (0,0): FG>=N and FP>=H.
 (1,0): FG<=N and FP>=Π_D.
 (0,1): FG>=C_D and FP<=H.
 (1,1): FG<=C_D and FP<=Π_D.
At boundaries multiple pure or mixed Nash; these conditions collectively exhaust PURE entry NE. Every selected pure entry NE plus specified lower-stage continuation constitutes a FULL PURE SPNE of the **ZONED** game, including all price/location histories within feasible strategy sets. Some entry parameter regions admit both public-only and private-only NE; equilibrium selection matters.

Compare equilibrium facility count to SECOND-BEST socially optimal count obtained by maximizing W00,W10,W01,W11 with exactly the same downstream technology and zoned location game. This is NOT the unconstrained first-best planner's joint location-and-entry optimum; both entrants' locations in W11 are Nash. Unrestricted planner can obtain 2v−1/24+k/2≈7.95944444444 at a=.25,b=.75,X=Y=.5, only ≈7.7e-06 higher than zoned decentralized W_D.

Examples:
- FG=.05, FP=.05 -> UNIQUE (1,1) NE, W11=7.8594367125; second-best chooses (1,1): efficient TWO.
- FG=.1, FP=.3 -> UNIQUE (1,1) NE, W11=7.5594367125; second-best chooses G-only W10=7.71777778: OVERENTRY of P, two instead of one.
- FG=8, FP=7.5 -> UNIQUE (0,0) NE; planner would choose P-only W01=7.81777778−7.5=.31777778>0: UNDERENTRY, zero instead of one.
- FG=.4,FP=1 -> BOTH (1,0) and (0,1) entry NE, count one either way; planner prefers G-only because FG<FP. Equal facility COUNT need not equal same SOCIAL WELFARE.
All numeric values are hypothetical and illustrative; no empirical calibration.

## Mathematical/comparative-scientific caveats
1. ZONING restriction is STRONG, NEW, EXPLICIT, nontrivial and required for proof. It is NOT an SPNE for unrestricted [0,1] location choice; full unrestricted location/fee game (Stages 12–15) remains OPEN/nonexistence-prone.
2. Normalize masses M=R=1, A=c=t=1, μ=.1, beta=2/3. This is one strictly admissible calibration. The rational fee system generalizes to a DOMAIN of positive k and large spacing d, but the Bernstein global location certification given is at this k, not arbitrary thick/thin markets. No general M/R comparative-static theorem.
3. Outside options are allowed in formal game but baseline v=4 ensures no opting out in any equilibrium OR unilateral fee deviation against the constructed price continuation. For arbitrary OFF-PATH prices, some users may opt out; demand equilibrium needs to be selected. With prices high on both sides, no-participation may have multiple fixed points. Existing proof of unique demand via clipped fixed point only applies to FULL COVERAGE, not universal high-high fee profiles. For a TRUE SPNE, **every off-path price vector** must induce a user-selection Nash equilibrium; use a fixed well-defined selection such as smallest masses from Tarski (cross-side positive feedback) when applicable. This requires careful construction not just contraction. Location/price Nash checks apply to unilateral fee changes only when fallback from the other hub guarantees full coverage. Do not mistake universal demand uniqueness for full opt-out game.
4. Unlike Stage 11 post-match investment (constant e), ex ante e_j depends on membership and social experiment value is cubic. It creates a genuine amplification force but has NOT been compared sufficiently with prior pre-match investment and network entry papers, so novelty NOT established.
5. Welfare calculated gross 2v is independent of price transfers. Fixed-cost resource weights assumed unit. No municipal direct provision vs private subsidies/procurement assessment. Mixed NE of entry game not characterized (straightforward 2x2 formula).
6. The client asked first for location/profit/welfare; next task should parameterize M,R,μ and derive EXISTING literature differential under/over thresholds in general, not announce journal contribution.

IMPORTANT TYPO AUDIT: Social planner ideal unrestricted a=.25,b=.75 gross W=2v−1/24 + k/2=7.9594444444 (since k=1/450); single N=1759/225, H=1639/225.

## Closest benchmarks
Kim (2024), Mixed duopoly in two-sided markets, Manchester School, DOI 10.1111/manc.12467; the full-coverage mixed public-private platform pricing/quality skeleton.
Basak & Petrakis (2021), JEMS, DOI 10.1111/jems.12431; simultaneous private entry and network externality over/under-entry.
Ghosh & Morita (2007), IJIO, DOI 10.1016/j.ijindorg.2007.02.002; bargaining and socially insufficient entry.
The formulas above are derived for a NEW RESTRICTED benchmark, not falsely attributed to those papers.
