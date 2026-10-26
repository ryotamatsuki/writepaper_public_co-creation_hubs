# Stage 11 — Global SPNE for the restricted bilinear two-sided matching benchmark
Date: 2026-10-27
Status: **CONSTRUCTIVE SPNE EXISTENCE PROVED IN RESTRICTED FIXED-LOCATION / ONE-LINK MODEL**.
Do not promote to proof of endogenous-location or endogenous partner-scope full model, uniqueness, or publication novelty.

## 0. Players / sequential stages / technological restrictions
Local municipality chooses g∈{0,1} and bears F_G≥0 if g=1. Exactly one potential private provider then chooses entry and pays F_P≥0 if it enters. A public-only or private-only hub uses position ℓ=1/2; when BOTH are active, public G uses a=1/4 and private P uses b=3/4, by exogenous service-technology specification. **Locations are NOT strategic in this benchmark**. Each active hub simultaneously sets nonnegative joining fees (p_j^I,p_j^R). Each mass-M innovator (type x~U[0,1]) and mass-R partner (type y~U[0,1]) chooses one hub or outside. Actual within-hub matched mass is L_j=μN_j^I N_j^R; μR<=1, μM<=1 ensure match-capacity feasibility. Any successful match is one-to-one and initiates sunk effort e before ex-post bilateral Nash bargaining. There are no other operating costs and no direct public informational advantage; no subsidy/cost of public funds wedge. Prices and Nash transfers are internal transfers in aggregate W. All optimization in stages is conditional on actual downstream market structure, not hypothetical duopoly.

Postmatched bargaining and effort imply e*=βA/c, net innovator per-match profit u=β²A²/(2c), partner payoff τ=β(1-β)A²/c, total social match surplus S=u+τ=β(1-β/2)A²/c. Innovator utility at hub j is v_I-p_j^I-t_I(x-ℓ_j)^2+μuN_j^R; partner is v_R-p_j^R-t_R(y-ℓ_j)^2+μτ N_j^I. Users absent hub have outside zero. Hubs are single-homing. When an indifferent endpoint joins a hub it is a measure-zero tie; single hub complete-coverage endpoint is counted as participating. No direct participant search cost and scope fixed at one.

## 1. Sufficient parameter domain (nonempty)
Set a_I=t_I/2>0, a_R=t_R/2>0, B=Mv_I+Rv_R, T=Mt_I+Rt_R, Z=μMR S. Require:
(A1) a_I a_R > μ²MR S²  (STRICT strong-sorting / welfare-concavity).
(A2) p_I*=a_I-μR τ ≥ 0, p_R*=a_R-μM u ≥ 0  (feasible candidate joining fees).
(A3) v_s >= max(3 t_s/4, p_s*+9 t_s/16), s∈{I,R} (sufficiently high base valuation).
(A4) μ M <=1 and μ R<=1; F_G,F_P>=0; β∈(0,1), A,c,t_s,M,R,μ>0.

A1 implies product-of-two-side sorting feedback μ²MR uτ/(a_I a_R) <1, since S²>uτ. This ensures UNIQUE full-coverage allocation for any unilateral price deviations with the other platform fixed at equilibrium fees. A3 guarantees that at the other platform's fees (p_s*), all users strictly prefer at least that platform to outside regardless of its positive match benefit, because max distance to a fixed endpoint (.25 or .75) is .75. Therefore unilateral fee deviations CANNOT reduce total platform attendance. For other off-path price profiles, compactness/continuity of continuum demand ensures at least one user participation equilibrium (pick any well-defined selection, e.g. lexicographically smallest vector of masses). The proposed SPNE is an existence profile, not necessarily unique.

## 2. Conditional two-hub market: closed-form fees and global price Nash
With both hubs, let x=N_P^I/M, y=N_P^R/R. Full coverage implies N_G^I=M(1-x), N_G^R=R(1-y).
The private-side indifference cutoff, in the interior, yields:
p_P^I-p_G^I=a_I(1-2x)-μR u(1-2y)
p_P^R-p_G^R=a_R(1-2y)-μM τ(1-2x).
At either market-side corner x or y ∈{0,1}, the same equation supplies an upper bound on the private fee revenue required to retain that mass; a side with zero customers contributes zero revenue.

Candidate public/private fees are:
p_G^{I*}=p_P^{I*}=a_I-μR τ,
p_G^{R*}=p_P^{R*}=a_R-μM u.
At these fees the unique full-coverage mass allocation is x=y=1/2, due to A1.

**GLOBAL PRIVATE FEE BEST RESPONSE**:
For any unilateral nonnegative private fee vector, full coverage follows from public fees remaining fixed at p_G*. Let (x,y)∈[0,1]^2 be ANY resulting demand equilibrium. Since a private side with positive membership cannot charge above the inverse-demand retention bound, its revenue is bounded by:
Q_P(x,y)=2 M a_I x(1-x)+2R a_R y(1-y)-Z(x+y)+2Zxy.
The private candidate profit is Π^D=Q_P(1/2,1/2)=(M a_I+R a_R-Z)/2=T/4-Z/2.
Exactly:
Π^D-Q_P(x,y) =
  2M a_I(x-1/2)^2+2R a_R(y-1/2)^2-2Z(x-1/2)(y-1/2) >=0.
The RHS is positive definite when a_I a_R>Z²/(4MR), which is implied by A1. It includes EVERY corner x=0,1 and y=0,1, unilateral fee cuts, and abandonment of a market side. Therefore NO profitable private fee deviation exists even when the price change induces tipping/corners. This is an analytic global proof, not a finite grid inference.

**GLOBAL PUBLIC FEE BEST RESPONSE**:
For any unilateral public fee vector, all users can still join the private platform at the candidate fixed P fees by A3. The global maximized surplus for given allocations is achieved by assigning the lowest types to G and higher types to P on EACH side; actual price-induced allocation has exactly that cutoff structure. Holding the active facilities fixed, the welfare loss relative to (x,y)=(1/2,1/2) is:
W_D(1/2,1/2)-W_D(x,y)
 =M a_I(x-1/2)^2+R a_R(y-1/2)^2
   -2Z(x-1/2)(y-1/2) >=0.
The RHS is positive definite by A1. Thus G also has no profitable price deviation, including cases in which one platform receives all users on one/both sides.

This establishes a genuine Nash equilibrium of the two-sided fee subgame with all noninterior price deviations excluded, under A1-A4.

## 3. One-hub and no-hub continuations
For ONLY one public G hub, place its service at 1/2 (a fixed service capability), choose p_G^I=p_G^R=0, and all innovators and partners join because v_s>=3t_s/4. This yields W_G=B-T/12+Z-F_G.

For ONLY one private P hub, position 1/2 and complete-coverage candidate fees:
p_P^{I,m}=v_I-t_I/4+μR u,
p_P^{R,m}=v_R-t_R/4+μM τ.
Gross private monopoly revenue is Π^M=B-T/4+Z. Both sides fully join with endpoint tie-breaking.

**GLOBAL PRIVATE MONOPOLY FEE PROOF**:
Any mass fractions x,y ∈[0,1] served by a sole platform at center must be sourced from the center-nearest segments, or else the maximal fee maintaining each membership mass is even lower. Thus firm revenue is at most
Q_M(x,y)= M x[v_I-t_I x²/4+μR u y]
        + R y[v_R-t_R y²/4+μM τ x].
For v_s>=3t_s/4, each partial derivative satisfies
∂Q_M/∂x >= M(v_I-3t_I/4) >=0,
∂Q_M/∂y >= R(v_R-3t_R/4)>=0
on the WHOLE unit square, because matching cross-terms are nonnegative.
Therefore Q_M maximizes at x=y=1, where the stated monopoly prices attain it. This rules out global partial-coverage fee deviations, unilateral abandonment of either side, and corner monopoly strategies. Center 1/2 weakly dominates any other solitary characteristic if the solitary provider is permitted to choose its location, by symmetric-uniform rearrangement.

No hub gives W_empty=0 and no matching. At zero public prices public-only full participation maximizes gross welfare, since per-user intrinsic social value v_s-t_s(x-1/2)² is positive for every user and S>0.

For arbitrary off-equilibrium price profiles, user allocation fixed points exist by continuity of the continuum enrollment map (each type's tie sets are of measure zero with t_s>0 and distinct duopoly positions). Select a fixed continuation allocation from the equilibrium correspondence. The payoff bounds at equilibrium rivals' fees are independent of which off-path user allocation is selected.

## 4. Private entry: two closed-form thresholds
Let
  D := Π^D =T/4-Z/2,
  H := Π^M =B-T/4+Z,
  C := T/16-Z/2,
  N := B-T/12+Z.
A3 implies H>D>0 (indeed B>=3T/4 implies H-D >=T/4+3Z/2>0).

Private entry AFTER government build occurs iff F_P<=D, with indifference resolved by entry.
Private entry when NO government build occurs iff F_P<=H, with indifference resolved by entry.

Realized SPNE social welfare for the four market structures:
- no hubs: 0
- G only: N-F_G
- P only: N-F_P
- G and P: B-T/48+Z/2-F_G-F_P
This last term is N+C-F_G-F_P.

## 5. Municipality setup choice: FULL piecewise closed-form policy rule
At stage 0 municipality chooses g=1 if W(g=1)>=W(g=0), with tie in favor of building:
CASE (i) 0<=F_P<=D: P enters irrespective of G.
  W(g=1)=N+C-F_G-F_P; W(g=0)=N-F_P.
  **Build iff F_G<=C**.
CASE (ii) D<F_P<=H: private entrant is deterred by public establishment.
  W(g=1)=N-F_G; W(g=0)=N-F_P.
  **Build iff F_G<=F_P**.
CASE (iii) F_P>H: no private facility in either contingency.
  W(g=1)=N-F_G; W(g=0)=0.
  **Build iff F_G<=N**.

All entry decisions above reflect complete backward induction and all market structures, not a conditional duopoly payoff reused in the monopoly continuation.

The incremental gross social welfare from P entering conditional on G present is C =T/16-Z/2. The private gross fee revenue is D=T/4-Z/2. Hence private over-entry wedge D-C=3T/16, independent of μ (though both entry levels separately depend on μ).
This wedge reproduces generic business-stealing and DOES NOT establish originality. The participation-extensive margin, endogenous location/scope and physical-v-digital policy comparison remain outside this restricted theorem.

## 6. Rational parameter witness and sample SPNE
M=R=t_I=t_R=1; μ=β=1/2; A=c=1; v_I=v_R=4; F_P=1/10; F_G=1/50.
All A1-A4 strict where required.
u=1/8; τ=1/4; S=3/8; Z=3/16; T=2; B=8.
Two-hub fees: p_I*=3/8, p_R*=7/16; each platform hosts half of EACH side, each makes μMR/4=1/8 matches, and private gross revenue D=13/32.
Private-only fees: p_I^m=61/16, p_R^m=31/8; all join, private gross revenue H=123/16.
Social use welfare for a single platform before fixed cost N=385/48; for two N+C=773/96.
C=1/32>F_G=1/50.
P enters with and without G. Municipality BUILDS; equilibrium social welfare with G+P is 773/96-1/10-1/50=19037/2400 ≈7.9320833; without G is 385/48-1/10=1901/240 ≈7.9208333.
Public-only (but NOT the realized outcome under private entry) would yield 385/48-1/50≈8.000833, illustrating socially excessive private entry if the public setup has already occurred.
Independent regression script: code/stage11_bilinear_spne_audit.py (exact rational payoff identities, 51×51 all-demand-share corners and 161×161 two-price deviations).

## 7. Science / publication limits
This Stage 11 CONDITIONAL SPNE is for the fixed positions and ONE possible match per innovator; both entire potential groups join an active hub because of deliberately high v. Matching realization is endogenous because it depends on endogenous choices of platforms by both sides, but it does NOT endogenize the overall participation masses, two-match scope, relative hub position, or the ownership-versus-private-implementation instrument.
A1-A4 are sufficient, not necessary. No uniqueness claim for other equilibrium branches or outside this region. Proof does not establish new scholarly mechanism beyond Armstrong (2006) cross-group feedback and classical free-entry distortions.
Next: introduce the smallest nontrivial partner-count / hold-up decision while retaining the stage-0/1 entry decisions and maintaining verified equilibrium existence.
