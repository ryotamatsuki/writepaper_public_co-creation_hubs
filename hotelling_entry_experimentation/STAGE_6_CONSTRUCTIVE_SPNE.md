# Stage 6 — Constructive SPNE under an explicit low-price continuation
Date: 2026-10-27
Status: **RESTRICTED EXISTENCE THEOREM (with a specified equilibrium selection); UNIQUENESS and NEW ECONOMIC CONTRIBUTION NOT ESTABLISHED.**

## 0. Changes from previous stages
- The Stage 4/5 numerical exercise had baseline user value v=2t, which was not sufficient for an immediate uniform bound on off-path full participation. For this Stage 6 existence proof impose v >= 4t; illustrative parameters use v=4, t=1.
- Stage 5 found a continuum of price Nash equilibria in the off-path all-private regime; we keep that correct observation.
- We now specify the **lowest-public-price Nash equilibrium** in that regime as a CONTINUATION SELECTION, not as a general behavioral refinement. A lexicographic taste for low public prices conditional on a fixed private price is not enough to select it across Nash price equilibria.
- We PROVE unilateral location optimality for the chosen continuation in the domain below, instead of relying only on finite price/location grids.
- No claim that every continuation supports the same location, entry or establishment decisions.

## 1. Assumptions and game
Types x uniform on [0,1], mass M>0; mismatch t(x-l)^2, t>0. Public hub G objective total surplus and private hub P objective its profit. Both prices are nonnegative and marginal cost is zero. Utilities: U_G=v-p_G-t(x-a)^2, U_P=v+delta-p_P-t(x-b)^2, outside option zero.

Five-stage game:
0. Municipality selects g in {0,1}; if g=1 selects matching scope n in {1,2} and pays F_G.
1. An associated pilot firm chooses relationship-specific effort e>=0, maximizing r_n e -e^2/2, r_n=alpha*n/[1+sigma(n-1)]. Optimal e_n=r_n; spillover delta_n=gamma*n*e_n; program direct social value H_n=e_n^2/2+B*n*e_n-k*n (separate from downstream user spillover).
2. One potential private operator chooses enter/no enter at fixed cost F_P.
3. Active hubs choose positions simultaneously.
4. Active hubs choose prices simultaneously.
5. Users choose provider or nonparticipation. Private gets customers at exact ties when a=b, delta>0 and p_P=p_G+delta. All other ties are zero-mass except at exactly equal locations.

When g=0, no pilot occurs and no delta exists; private monopoly optimizes its own location/price. When g=1 and private stays out, public monopoly optimizes social allocation.

## 2. Selected price continuation for EVERY a,b in [0,1]
Let d=delta/t in (0,1/4) and D=|b-a|.

(a) b>a:
 z=(a+b)/2-delta/(2tD).
 If z>0 (interior efficient allocation), select p_G=p_P=2tD(1-z); P share=1-z.
 If z<=0 (all-private efficient allocation), define Delta=delta-tD(a+b)>=0 and select
 q_low=max{0,2tD-Delta}; p_G=q_low; p_P=q_low+Delta; P share=1.

(b) b<a:
 z=(a+b)/2+delta/(2tD).
 If z<1, select p_G=p_P=2tD*z; P share=z.
 If z>=1, define Delta=delta-tD(2-a-b)>=0 and select
 q_low=max{0,2tD-Delta}; p_G=q_low; p_P=q_low+Delta; P share=1.

(c) b=a:
 Select p_G=0,p_P=delta; users choose P on indifference.

In each efficient interior regime, public best responds by implementing efficient consumer allocation, which requires p_G=p_P; the private price first-order condition yields the stated p. In each all-private region the selected q and q+Delta are Nash prices: P optimally sets the highest all-market price if q+Delta>=2tD, and G obtains the efficient allocation. In exact collocation, P has quality advantage delta and p_G=0,p_P=delta is Nash with P-favoring tie-break.

All selected q are <=2t and all interior selected prices are <=t+delta<1.25t. If v>=4t, the public service alone is accessible at equilibrium q to every customer even at the worst taste position; hence no one opts out at these selected prices, and any private unilateral price change is disciplined by a viable public alternative. This proves the chosen price continuation is a valid Nash equilibrium in every post-location price subgame. Other Nash equilibria exist, especially for all-private regions.

## 3. Main restricted existence theorem (conditional on the selected continuation)
Let d=delta/t and
 d_min=(sqrt(5)-2)/4 = approximately 0.059016994.
If v>=4t and d_min<=d<1/4, then under the above selected price equilibria the following pair is a Nash equilibrium in the endogenous location stage:

  a*=1/4-d,   b*=3/4-d,
  z*=1/2-2d,  p*=t/2+2delta.

Both hubs remain interior, b*-a*=1/2, and their shares are z* and 1-z*.
Private operating profit / M = (t+4delta)^2/(4t).

Proof of G's GLOBAL best response:
Given b=b*=3/4-d, efficient social welfare equals its value when all users use P plus the gain from serving efficiently with G. For a<b define D=b-a, z=b-D/2-d/(2D). If z>0, welfare improvement /t is D*z^2; otherwise zero.
Derivative of D*z^2 has the sign of z[b-3D/2+d/(2D)]. The second factor has exactly one positive root D=1/2 when b=3/4-d, and changes from + to -. Hence G's global optimum in a<b is a=b-1/2=1/4-d; gain /t is (1/2)(1/2-2d)^2.
For a>b define D=a-b<=1-b=c=1/4+d. When G can attract a positive right-hand tail, its gain /t is D*[c-D/2-d/(2D)]_+^2. Since D/2+d/(2D)>=sqrt(d), the gain is at most c(c-sqrt(d))^2. This is STRICTLY less than (1/2)(1/2-2d)^2 for all 0<=d<1/4. When a=b the gain is zero. Thus no G deviation in [0,1] is profitable.

Proof of P's GLOBAL best response:
Fix a=a*=1/4-d. For b>a write s=b-a. The all-private regime is 0<s<=2d. The lowest-price continuation gives P profit/(Mt) = max{Delta/t, 2s} <= max(d,4d)=4d. But at b*=a+1/2 the profit/(Mt)=(1/2+2d)^2 >=4d, strictly for d<1/4.
For s>2d (interior), P share is h(s)=3/4+d-s/2+d/(2s); profit/(Mt)=2s h(s)^2. Its derivative has the sign of -(s-1/2)(s-2d/3), and since s>2d, the unique maximum is s=1/2.
For b<a, let s=a-b<=a=1/4-d. P profit/(Mt)<=max(d,2s)<=max(d,1/2-2d). The candidate profit (1/2+2d)^2 is >=both terms provided d>=d_min=(sqrt(5)-2)/4. At b=a, P profit/(Mt)=d, also lower. Thus no P deviation in [0,1] is profitable.

This proves a global LOCATIONAL equilibrium under a COMPLETE chosen continuation of Nash price equilibria. It does NOT imply uniqueness of the equilibrium of the original extensive-form game.

## 4. Lift to a full SPNE of the restricted entry and establishment game
For any n in {1,2}, assume d_n=delta_n/t satisfies d_min<=d_n<1/4.

- Under g=1 and P entry, use the above continuation and location equilibrium. Expected private operating profit:
  Pi_n(M)= M*(t+4delta_n)^2/(4t).
- Under g=1 and no P entry, public locates at a=1/2 and picks p_G=0; all users served and gross welfare:
  W_public_alone=M*(v-t/12)+H_n.
- Under g=0 and P entry, private locates at b=1/2 and charges p=v-t/4 (full coverage, global monopoly optimum for v>=4t). Operating profit M*(v-t/4); total gross surplus M*(v-t/12).
- Under g=0 and no P entry, zero service and zero gross welfare.

Entry (tie breaks toward entering):
with public pilot n: enter iff F_P<=Pi_n(M);
without public: enter iff F_P<=M*(v-t/4).

The incremental gross social surplus of private entry in the selected mixed equilibrium, relative to one public hub at 1/2, is:
  DeltaS_n(M)=M*(t+4delta_n)^2/(16t).
The policy welfare after g=1,n is:
  W_n=M*(v-t/12)+H_n-F_G+1{FP<=Pi_n(M)}*[DeltaS_n(M)-FP].
No-public welfare is W_0=1{FP<=M(v-t/4)}*[M(v-t/12)-FP].
Stage 0 selects argmax{W_0,W_1,W_2}. At equality select the indicated fixed tie-break.

Every continuation (pilot investment, entry, position, price and adoption) is optimal on its own subgame under the explicitly specified off-path price selection. Hence this CONSTRUCTS AT LEAST ONE SPNE of the restricted complete game (not the unqualified full empirical model).

## 5. Recalibrated illustration
t=1, v=4, alpha=.09, sigma=1.4, gamma=1, B=.5, k=.001, F_P=.55, F_G=.4, M variable.
n=1: e=.09, delta=.09, H=.04805; (a,b)=(.16,.66); gross private profit/M=.4624; DeltaS/M=.1156.
n=2: e=.075, delta=.15, H=.0758125; (a,b)=(.10,.60); gross private profit/M=.64; DeltaS/M=.16.

Policy intervals away from indifference:
- M<~.08277128: no public hub
- ~.08277128<M<.859375: build, n=2
- .859375<M<~1.18944637: build, n=1
- ~1.18944637<M<2.026171875: no public hub
- M>2.026171875: build, n=2
At any threshold, tie-breaking and right-/left-continuity must be specified.

Other useful threshold: the no-government private monopoly enters at M=.1466666667. This is correctly included in W_0. All values are hypothetical, not measured policy estimates.

## 6. Interpretation and STRICT limits
1. This is an **existence** result for a restricted SPNE with a fully specified equilibrium-selection convention. It does NOT prove SPNE uniqueness or robustness to other off-path Nash-price continuations. The low-price selection is an EXOGENOUS equilibrium selection; lexicographic public price minimization at each best reply alone does not pin it down.
2. The market size comparative statics rely on exactly one prospective private entrant and H_n not scaling with M.
3. The "scope" n is the potential number of experiment partners, not an endogenous negotiated matching network; the investment recovery coefficient r_n remains reduced form.
4. Matching and spillovers do NOT depend on geographic/service-characteristic distance and cannot yet determine the advantage of physical public facilities over subsidies, procurement or digital platforms.
5. No new economics is claimed from pricing multiplicity, free-entry inefficiency or conditional nonmonotone policy; those mechanisms are close to Cremer–Marchand–Thisse (1991), Zhang–Li (2013), Kang–Kim (2022), Mankiw–Whinston (1986) and Matsushima–Shinohara (2014).
6. The next critical stage is policy-instrument dominance (direct public ownership versus a subsidy/contract yielding the same pilot spillover) and a genuine bargaining/matching microfoundation. If equivalent contracted provision can deliver the pilot's H_n and delta_n at no higher total social cost, the claim that a publicly owned physical hub is NECESSARY fails.

## Reproducibility
- Executable script: code/stage6_constructive_spne.py
- Tests include all private/public unilateral location deviations on grids, 21x21 arbitrary location profiles with sampled price BR checks, and reported stage-0 policy outcomes.
- Analytic proof above, not the grids, is the equilibrium certificate.
- **Scientific verdict**: mathematical conditional GO; independent referee and prior-art novelty gates HOLD.
