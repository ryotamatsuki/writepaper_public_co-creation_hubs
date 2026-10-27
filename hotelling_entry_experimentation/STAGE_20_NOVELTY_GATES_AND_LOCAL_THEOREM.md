# Stage 20 — What theorem would constitute a publishable contribution? An independent local derivative audit

Date: 2026-10-27.
Status: **CONCRETE LOCAL THEOREM IN THE ZONED MODEL, NOT NOVELTY-CERTIFIED; STRONG ABSORPTION SCREEN FROM MATCHING-INVESTMENT LITERATURE.**

## 1. Baseline is Stage 19 (not unrestricted Hotelling game)
Symmetric unit-mass two-sided subscription, full market cover on the equilibrium path, β=2/3, t=A=c=1, matching success μ N_R, innovation investment BEFORE match realization. Define
 k=β² A² μ²/(2c)=2A²μ²/(9c).
In region where the zoned location continuation is unique and interior, define a(k),b(k), x(k),y(k), public/private 2-sided fees gI,gR,pI,pR, P gross duopoly fee revenue P_D(k) and two-hub gross social welfare V_D(k). Locations G∈[0,.4],P∈[.6,1] (including monopolist subgames); base valuation v=4. N(k)=public-only gross social value, H(k)=private-only gross operating revenue.

## 2. NEW exact low-k Taylor result — proof verified with SymPy and numerical finite-differences
At k=0, the interior Nash location satisfies a0=1/4,b0=3/4,x0=y0=1/2.
The location FOCs are public social envelope W_a=(1−x)^2+(1−y)^2−2a(2−x−y)=0 and P re-optimized revenue derivative ∂Π_D/∂b=0.
Differentiating these simultaneous FOCs wrt k at k=0 yields Jacobian
  J=[[-3/2,1/2],[1,-3]]
and direct shock partial ∂(W_a,Π_b)/∂k=[0,-4]'.
Hence
  a'(0)=−1/2, b'(0)=−3/2,
  x'(0)=y'(0)=1,
  pP_I'(0)=−1/2, pP_R'(0)=−3/2.
The private operating revenue derivative is zero:
  Π_D(0)=1/2, Π_D'(0)=0,
  Π_D(k)=1/2+O(k²).
The two-hub social welfare:
  W_D(0)=2v−1/24,
  W_D'(0)=+1/2,
  W_D(k)=2v−1/24+k/2+O(k²).
Single-hub public gross welfare:
  N(k)=2v−14/75+2k.
Marginal SOCIAL incremental gross value of adding P to public-alone market:
  Δ(k)=W_D(k)−N(k)
      =29/200−(3/2)k+O(k²).
PRIVATE-v-SOCIAL entry wedge:
  Π_D(k)−Δ(k)=71/200+(3/2)k+O(k²).

**Precisely interpretable proposition**: In the neighborhood of k=0 with zoned mixed-duopoly price/location continuation, increasing pre-match matching-investment complementarity INCREASES the private-entry-versus-social-entry wedge even though total social surplus increases in both market structures. The link amplification raises ONE-hub experimental surplus at rate 2 compared with TWO-hub at rate 1/2, while private entrant operating profit has zero first-order change. This is a rigorous **local comparative static**, not proof of a global parameter-region crossing or journal novelty.

## 3. The existing under-entry numerical example has NOTHING to do with matching investment
The monopoly gross private revenue:
 H(k)=2v−18/25+2k.
Single-hub social welfare:
 N(k)=2v−14/75+2k.
Therefore,
  N(k)−H(k)=8/15  (IDENTICALLY, all k in benchmark).
The possible inefficient no-hub equilibrium arises with FG≥N and H<FP<N. Its under-entry interval width 8/15 is COMPLETELY independent of μ, β, A,c and k.
This reflects consumers' inframarginal surplus / private monopoly fee extraction, not a new matching-investment mechanism. A paper claiming that Stage19 demonstrates matching-caused both over-/underentry would be scientifically unjustified.

## 4. Direct absorption literature found and VERIFIED (2026-10-27)
- Kim (2024) Manchester School mixed public/private two-sided platforms with endogenous quality investment.
- Basak & Petrakis (2021) JEMS both entry insufficiency and excess under network externalities and compatibility.
- Tan & Zhou (2021) ReStud multi-sided platform entry over/under and sign reversals of competition.
- Ghosh & Morita (2007) IJIO bargaining and insufficient entry.
- Mangin & Julien (2021) JET "**Efficiency in search and matching models: A generalized Hosios condition**": efficient market-side entry depends on sum of matching elasticity AND endogenous match-surplus elasticity. This is EXTREMELY close to Stage19's main pre-match investment force: higher market thickness increases both match PROBABILITY and joint match VALUE, so a generic statement on investment-amplified over-/underentry is ABSORBED unless the public-private two-platform equilibrium entails a qualitatively new ownership/number/selection effect. https://www.sciencedirect.com/science/article/pii/S0022053121000259
- Bidner, Roger & Moses (2016) AEJ Micro "**Investing in Skill and Searching for Coworkers: Endogenous Participation in a Matching Market**": pre-search skill investment, coordination failures, AND separate matching marketplaces. Preinvest/search/participation by themselves not novel. https://www.aeaweb.org/articles?id=10.1257/mic.20140110
- Atakan, Richter & Tsur (2026) Theoretical Economics "**Efficient investment, search, and sorting in matching markets**": some matching+investment+Nash bargaining structures are constrained efficient. Thus generic statements "ex-ante investment + matching implies inefficiency" are FALSE without qualifications. https://doi.org/10.3982/TE6786
- Aoyagi & Yoo (2022) GEB matching strategic agents on two-sided platform: optimal matching mechanism can deliberately destroy allocative matches to extract rents.

## 5. Rank / KILL criteria for publishable primary propositions
P1 (candidate HIGHER): **A genuinely new entry-distortion sign threshold arising from two-sided pre-match investment, with M/R imbalance and municipal player's simultaneous strategic participation; show an interval where an *actual* socially underprovided second (or public) hub would not appear in the post-match effort or compatible-network comparator.** Needs new general M,R and an implementable pure/mixed equilibrium, analytic sharp condition, welfare; no such theorem in Stage19.
P2 (candidate HIGHER): **Public entry has a true strategic complementarity with profitable private entry from verifiable institutional neutral access/physical testbed, despite two-sided business stealing; find exact parameter space where government and private both build but neither alone can enter, and direct public provision dominates cost-equivalent procurement/subsidy under precise incomplete contracting.** NOT possible in Stage19 current payoff ranking N≫Δ, H≫P_D. Requires genuine new link/spillover/commitment primitive; do not call it a mere "extension".
P3 (candidate MEDIUM): **Entry sign reversal across market sizes** due to imbalance / investment complementarity while incumbent private hubs are profit maximizers and public hub internalizes extra-project knowledge spillovers. Strong absorption risk from Basak/Petrakis and Mangin/Julien; must produce nonmonotone unusual thresholds or new policy optimum unavailable in parents.
P4 (potential supplementary): **Government social planner cannot decentralize its first-best number of facilities using single-sided per-subsidy at any fixed schedule when partner-side investments are unverifiable; a public, neutral hub implements it**, with precise implementation failure. But subsidy taxes and complex instruments must be modeled, otherwise not credible.
P5 (NOT publishable alone): Pure SPNE existence/nonexistence, entry examples under fixed costs, cubic location roots, \`a=.24884\`, "preinvestment generates cubic social surplus", or simply "too few and too many occur".

## 6. Immediate research direction
Before generalizing all of Stage19, derive a sharp comparative-statics theorem in a SIMPLE price equilibrium where M,R and μ vary and both population participation and investor effort are endogenous; explicitly compare:
(i) before-match vs after-match investment,
(ii) municipal entry vs private-only entry,
(iii) single-homing fragmented vs connected/partly interoperable matches,
(iv) optional public funding instruments.
A general matching-investment surplus elasticity is ALREADY studied (Mangin-Julien), so a new threshold needs public/private ownership or simultaneous entry to matter.
The zone restriction is a substantive implementation detail and cannot be the sole driver of the theorem.

*Novelty status remains HOLD.* No likelihood-of-acceptance guarantee for IJIO/JRS or another field journal.
