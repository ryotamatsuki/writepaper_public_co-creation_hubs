# Stage 23 — Simultaneous urban–rural public/private ownership: solved restricted entry game

Date: 2026-10-27. **Proved complete PURE first-stage Nash/SPNE conditions and socially optimal entry for a DISCRETE-site, pre-match-investment, FIXED-COMMISSION benchmark; NOT solved unrestricted strategic two-sided price competition and endogenous sorting. No novelty certified.**

## 1. Actors/timing and all 16 profiles
One prefectural social-welfare maximizer chooses g=(g_U,g_R)∈{0,1}², two private firms independently choose p_U,p_R∈{0,1}; all choose SIMULTANEOUSLY. Therefore 4×2×2=16 histories, including overlapping public and private sites in one region, 3-/4-site alternatives and one-/zero-site configurations. The one public actor may deviate in BOTH regions at once and must compare four g-vectors.

Urban potential investor/partner masses 1,1; rural masses r,r with 0<r≤1. Without a hub in their OWN region, agents do not enter the regional matching market (remote attendance prohibitively costly). One hub locally serves all residents; when G and P coexist in a region, two predetermined micro catchments split investors and partners 50:50. Each active site has constant use value v and minimized quadratic mismatch cost t/(12 n_i²) per local user, n_i=g_i+p_i. These restricted service assignments are NOT free price/sorting equilibria. Neither platform chooses a fee; private intermediaries earn the stipulated commission σ∈(0,1) on the partner's ex-post Nash bargaining receipt.

Public hub's cross-institutional access coefficient is γ_G=1; private's is γ_P=α∈(0,1]. A link between sites in DIFFERENT regions has efficiency ρ γ_h γ_h', ρ∈[0,1]. Thus PP can interconnect but less effectively if α<1; do NOT postulate that only public sites can interconnect. **The α difference is REDUCED-FORM contractual/technical friction; public ownership not yet endogenously justified.**

## 2. Backward induction
For local hub h in region i, enrolled innovator mass I_ih, partner mass R_ih, weighted remote partner pool Rhat_j=Σ_{h' in j} γ_{h'} R_jh':
 q_ih=μ[R_ih+ρ γ_h Rhat_j]≤1.
Innovator ex ante (before match outcome) chooses effort to maximize q_ih β A e−c e²/2, hence uniquely:
 e_ih*=β A μ[R_ih+ρ γ_h Rhat_j]/c.
Each successful project ex post shares Ae with innovator β, partner 1−β. Define κ=(β−β²/2)A²μ²/c.
 Gross real experiment value = Σ_{i,h} κ I_ih[R_ih+ρ γ_h Rhat_j]².
 Gross social B(g,p) = Σ_{i:n_i>0}(M_i+R_i)[v−t/(12 n_i²)] + experiment value.
 Net welfare = W(g,p)=B(g,p)−Σ_i(F_Gi g_i+F_Pi p_i). Price commissions are transfers.
 A private firm's GROSS intermediary revenue is
 Π_i^op = σ(1−β) A μ R_iP [I_iP e_iP* + ρ α Σ_{h' in j} γ_h' I_jh' e_jh'*].
 Its NET entry payoff = Π_i^op−F_Pi, or 0 if absent.

## 3. EXACT exhaustive pure Nash characterization (not merely four pairwise cases)
(g,p) is a simultaneous-entry PURE NE iff
  (A) W(g,p)≥W(g',p) for ALL FOUR g'∈{0,1}²;
  (B) if p_i=1, Π_i^op(g,p)≥F_Pi;
  (C) if p_i=0, Π_i^op(g,p_i=1,p_-i)≤F_Pi,
for i=U,R. Every deviation's continuation recomputes user split, all matching probabilities, pre-match investment, commissions and total welfare. This is backward induction; entry is solved SIMULTANEOUSLY, not sequentially. The planner maximizes W over ALL sixteen profiles, not four named forms.

## 4. CLOSED-FORM four central ownership forms with exactly one site in each region
With z=ρ γ_U γ_R, the investors' experimental REAL surplus is exactly
 S(r,z)=κ[1+r³+r(1+r)(2z+z²)].
Let H(r)=2(v−t/12)(1+r).
 W_PP=H+S(r,ρ α²)−F_PU−F_PR;
 W_PG=H+S(r,ρ α)−F_PU−F_GR;
 W_GP=H+S(r,ρ α)−F_GU−F_PR;
 W_GG=H+S(r,ρ)−F_GU−F_GR.
PG is URBAN PRIVATE, RURAL PUBLIC; GP is URBAN PUBLIC, RURAL PRIVATE.

Threshold conditions (pairwise, not in isolation sufficient for 16-profile planner optimum):
 PG>PP iff
 κr(1+r)[2ρα(1−α)+ρ² α²(1−α²)] > F_GR−F_PR.
 GG>PG iff
 κr(1+r)[2ρ(1−α)+ρ²(1−α²)] > F_GU−F_PU.
 PG>GP iff F_GU−F_PU > F_GR−F_PR.
When α=1 all two-site configurations have IDENTICAL real gross benefits, so if F_Gi>F_Pi at both locations, the PP one-site-per-region structure strictly dominates the other three, for ANY r and ρ. **No direct-public-benefit result follows from simple interconnection absent specific public/private technological or contract differences.**

Let φ=σβ(1−β)μ² A²/c. For one P per relevant region:
 Π_U^op=φ[1+z r(1+r)+z²r];
 Π_R^op=φ[r³+z r(1+r)+z²r²].
Off-path simultaneous-entry NE still requires all 16 profiles, including deviations that form collocated public/private sites.

## 5. EXACT RATIONAL calibration and phase transitions
Assumptions, not empirical data:
 v=1/20, t=1/50, β=2/3, A=c=1, μ=2/5, σ=7/10, α=1/2.
Costs: F_PU=.018,F_PR=.012,F_GU=.110,F_GR=.037.
Thus κ=16/225, φ=28/1125. q_ih≤0.8 across r,ρ≤1.

At r=.7, the social planner's global maximum over ALL SIXTEEN profiles has the following unique ownership form except at boundaries:
 - PP: 0≤ρ<0.4978979545925306;
 - PG: 0.4978979545925306<ρ<0.7095688787169052;
 - GG: 0.7095688787169052<ρ≤1.
The cutoffs follow from exact-rational quadratic differences:
 W_PG−W_PP=(714ρ²+1904ρ−1125)/45000,
 W_GG−W_PG=(714ρ²+952ρ−1035)/11250.
Quadratic upper-envelope checking of all 16 profiles, not mere ρ sampling, certifies these transitions; adjacent plans tie at each cutoff.

At r=.7,ρ=.6, ALL FOUR entry configurations PP,PG,GP,GG are PURE Nash equilibria, but ONLY PG maximizes welfare. Every equilibrium has 2 sites; ownership, not count, is the inefficiency:
 PP: W=.2571262222, gross private revenue U=.0297235556,R=.0132539556;
 PG: W=.2632248889, gross private U=.0353422222;
 GP: W=.1962248889, gross private R=.0185198222;
 GG: W=.2448462222, no private.
At r=.2,ρ=.5, NE include urban-only P- and GG, while the planner chooses PP; the P- equilibrium has UNDER-provision in total hub number.

A complementary comparative-geography experiment changes only public fixed costs to F_GU=.035,F_GR=.115 and fixes r=.7,ρ=.6: GP becomes planner's UNIQUE optimum, while PP,PG,GP,GG are all pure Nash. Thus PG does not follow from model labels alone.

For owner-neutral connectivity α=1 with publicly more expensive costs, the planner never chooses government ownership in one-site-per-region comparisons — a critical implementation/novelty null.

## 6. Necessary qualifications
- This is a deliberately RESTRICTED fixed-fee/split, region-anchored game: no free within-region Hotelling location choice, no platform fee BR, no endogenous user provider sorting, no ex ante interoperability choice. The user's fuller Stage19 game is not solved by this.
- These are ALL PURE entry equilibria. Mixed Nash entry profiles and equilibrium selection are not analyzed. Never select PG unilaterally as the predicted market structure where PP/GP/GG also satisfy Nash.
- Higher public interface openness is assumed. Absent a reason private contracts, outsourcing, subsidies or mandated open access cannot implement γ=1 at lower social resource cost, this is NOT a proof of the necessity of public OWNERSHIP.
- Strategic public access choice/contract feasibility and endogenous local sorting are top priorities; price-subgame nonexistence of unrestricted stage14 may reappear if restored carelessly.
- All parameters normalized illustrative; journal originality unproved. Existing model absorption risk: Kawasaki (2026), Zhu et al (2021), Jullien-Pavan-Rysman (2021).

Local reproducible outputs available: \`stage23_two_region_ownership.py\` (Fraction EXACT 16-profile checks) and \`stage23_two_region_equilibria.csv\` (420 parameter-pair exhaustive NE/planner profiles).

## 7. Additional rural-thickness thresholds (holding rho=.6)
Using exactly the same cost vector and alpha=.5, the one-site-per-region pairwise differences yield two explicit positive-root scales:
 - PP vs PG: r*(1+r) = (.037-.012) / [κ*(2*.6*.5*.5+.6²*.5²*(1-.5²))],
   giving r_1=0.59846832137355.
 - PG vs GG: r*(1+r) = (.110-.018) / [κ*(2*.6*(1-.5)+.6²*(1-.5²))],
   giving r_2=0.81797912180628.
Both thresholds were checked against the **global 16-profile planner** at multiple neighboring and interior points. In this calibration at rho=.6, r=0.2,.4,.5 planner=PP; r=.6,.7,.8 planner=PG; r=.85,.9,1 planner=GG. Thus the MODEL DOES NOT imply \"public is more necessary the thinner the rural market\"; when the population is too thin cross-region gains cannot finance the higher government fixed cost. Potential NO-HUB or single-hub optimality must be rechecked if other fixed costs change.
At r=.7,rho=.6, replacing alpha=.5 by alpha=1 while retaining public costs above private makes planner PP, W_PP=.36184622. This is a core OWNER-NEUTRAL IMPLEMENTATION equivalence check: the public ownership advantage is exclusively supplied by the presumed difference in interoperable access.
