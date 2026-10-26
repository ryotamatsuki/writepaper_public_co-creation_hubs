# Stage 13 — Endogenous location, backward induction, and the off-path price-continuation obstruction
Date: 2026-10-27.
Status: **POSTMATCH BARGAINING / MATCHING / CONDITIONAL INTERIOR FEES / LOCATION FOCs SOLVED. FULL LOCATION/ENTRY/ESTABLISHMENT SPNE NOT YET CERTIFIED.**
Do NOT misrepresent the location cubic as a global SPNE. This stage actively tests the Stage 12 candidate with rigorous OFF-PATH price equilibria.

## 0. Game (precisely specified)
Stage0: government build G at cost F_G or no G.
Stage1: one potential private P enters at cost F_P or stays out, observing G.
Stage2: existing hubs choose service location(s) simultaneously; G maximizes social W, P maximizes operating profit, each location in [0,1].
Stage3: both hubs (if active) choose nonnegative fees on the innovator side and partner side simultaneously.
Stage4: mass-1 each of innovators and partners, taste uniformly distributed on [0,1], independently choose the highest utility active hub or outside.
Stage5: within-hub realized matching link mass L_j=mu*N_j^I*N_j^R (one link per participant maximum), with mu<=1.
Stage6: upon matched link innovator makes e>=0 investment, e costs c*e²/2.
Stage7: Nash bargaining with innovator share beta of Ae gross project value.
Fee and Nash transfers cancel from social welfare. Match net surplus enters once. Fixed costs are resources.
The equilibrium analysis below selects t_I=t_R=t, M=R=1, beta=2/3, A=c=1, mu=1/10, baseline v_I=v_R=4. Then u=tau=w=2/9, total match surplus S=4/9, network coefficient r=mu*w/t=1/45. Both sides are symmetrical. The participant outside-option constraint is slack in all reported examples.

## 1. Backward induction to subscription and prices
e*=beta*A/c; u=beta²A²/(2c), tau=beta*(1-beta)*A²/c, S=u+tau.
For two active hubs a<b, d=b-a, m=(a+b)/2, let x denote the PRIVATE share of EACH of the two sides. Conditional on interior efficient subscription and full coverage:
 (d-2r)(x-1/2)=d*(1/2-m)
 x=1/2+d*(1-a-b)/(2(d-2r)).
Two-sided Nash fees:
 p_P^I=p_P^R=2t*(d-r)*x,
 p_G^I=p_G^R=2t*d*x-t*r.
Private gross fee revenue:
 Pi_P^op(a,b)=4t*(d-r)*x².
These equations apply ONLY when efficient fractions 0<x<1, d>2r, all candidate fees nonnegative and all-price deviations are unprofitable. For this strictly interior domain, the Stage 11 inverse-fee-revenue quadratic bound generalizes: firm P profit shortfall for any shares X,Y is 2t d[(X-X*)²+(Y-Y*)²]-2mu*S*(X-X*)(Y-Y*) when M=R=t=1, which is positive definite if d>2r. Public W shortfall is t d[(X-X*)²+(Y-Y*)²]-2mu*S*(X-X*)(Y-Y*), positive definite when d>2r. With high v, the fixed rival fee covers the entire population even under unilateral fee changes. Hence the formulas define a GLOBAL price Nash equilibrium for the given interior a,b, not merely a local FOC solution.

## 2. Endogenous-location interior candidate
At fixed private b, public location FOC in the efficient interior:
 a=(1-x)/2.
Private location FOC given fixed public a:
 (d-r)*(1-2d)=r*x, once joint public FOC imposed.
Eliminating x:
 x=[d(1-d)-2r]/(d-4r),
 -2d³+(1+11r)d²-(6r+8r²)d+6r²=0.
At r=1/45 this polynomial has three real roots approximately
 d=.0287094683948, .1058066254652, .48770612836218.
Only the largest root respects d>2r, a,b∈[0,1], x∈(0,1), both prices>=0.
Candidate:
 d=.4877061283621769
 x=.5150339452401376
 a=.2424830273799312
 b=.7301891557421081
 p_P=.4794800252501016, p_G=.4801482005941077, per side
 Pi_P^op=.4938969781368012 (pre-fixed cost)
 Mixed gross W_use=7.9803069527047180 (pre fixed costs).
The PUBLIC's globally efficient location given b fixed, allowing the entire [0,1] deviation domain and planner-optimal allocation, is numerically a≈.2424830, agreeing with the FOC candidate. This grid/numerical bound is NOT a symbolic global location best-reply proof and cannot itself assert a SPNE.

## 3. Crucial non-interior price Nash equilibria (RIGOROUS obstruction)
Fix public a≈.24248302738. Let private DEVIATE in stage2 to b_dev=a+d0 with d0=1/100=.01<r and b_dev<1/2, so private's new location is slightly more central.
Define m0=(a+b_dev)/2 and Delta=r-2*d0*m0=r-d0*(a+b_dev)>0.
For any public fee q in {0,1}, consider on BOTH sides:
  p_G^I=p_G^R=q;
  p_P^I=p_P^R=q+Delta.
Then private hub serves EVERY user in both sides; there is no nonparticipation because v=4.
This is an exact Nash price equilibrium, not just a candidate:
- SOCIAL optimality of all-P for public: For any possible allocations X,Y of the two sides, write s=(X+Y)/2,h=(X-Y)/2. Holding a,b fixed, social W is a quadratic with s² coefficient +4*(2r-d0)*t [positive since d0<r<2r], and h² coefficient negative. Hence global planner-optimal allocation is either all-P or all-G; private b_dev closer to 1/2 than a, so all-P strictly maximizes W. G's q is a best response because all-P obtains globally maximum W.
- PRIVATE global price BR: Given G fees q on both sides, any unilateral private fee vector inducing private shares (X,Y)∈[0,1]^2 has total revenue bounded by
  Q(X,Y)=[q+2d0(1-m0)-r]*(X+Y)-2d0*(X²+Y²)+4r*X*Y.
  Let s=(X+Y)/2,h=(X-Y)/2:
  Q=2[q+2d0(1-m0)-r]*s+4(r-d0)*s²-4(r+d0)*h².
  With d0<r, this is convex in s and maximized in h at h=0. On s∈[0,1], global maximum occurs at s=0 or s=1; Q(1,1)=2(q+Delta)>Q(0,0)=0. So NO private fee change on either or both sides can increase profit, including loss of one side, tipping, and all boundary shares. Thus (q,q+Delta) is a Nash price profile.
- Baseline value v=4 ensures the public alternative at q≤1 covers all types, so the bound remains valid if P deviates even at high fees.
For q=0: Delta=.0172725616746236, P operating revenue=.0345451233492472.
For q=1: P operating revenue=2.0345451233492472.
The old interior candidate P operating profit is .4938969781368012. Thus q=1 yields a PROFITABLE location deviation, and q=0 does NOT. This is DIRECT evidence that the original cubic/entry threshold are NOT automatically robust to off-path price Nash equilibrium selection.
This does NOT prove no pure SPNE exists: SPNE uses one price Nash continuation at each off-path location, not all price equilibria, and a suitable low-price continuation might deter this particular deviation. All possible off-path continuation regimes must be addressed.

## 4. Monopoly histories closed form and CONDITIONAL higher stages
With only one hub, unconstrained optimum hub location is 1/2. For both sides v=4,t=1:
 each side monopoly full-coverage joining fee v-t/4+mu*w=3.7722222222;
 gross private monopoly fee profit
   H=2*(v-t/4)+mu*S = 8-1/2+2/45=7.5444444444.
One-hub gross social welfare
   N=2*v-2*t/12+mu*S=8-1/6+2/45=7.8777777778.
With two hubs at interior FOC candidate the private gross operating profit is Dpi=.4938969781368012, gross mixed social welfare 7.9803069527047180, and marginal mixed SOCIAL benefit vs a single public hub:
   C=7.9803069527047180-N=.1025291749269403.
These values can be used in a CONDITIONAL proposed backward-induction structure:
  FP<=Dpi: private would enter with and without public; government builds if FG<=C.
  Dpi<FP<=H: private would enter ONLY without public; government builds if FG<=FP.
  FP>H: private never enters; government builds if FG<=N.
**DO NOT promote these conditions as proved SPNE thresholds** until the location equilibrium has passed ALL off-path price and location best responses. The high-price Nash continuation constructed above invalidates a naive continuation-independent assertion.

## 5. What is solved and what remains
SOLVED: ex post Nash sharing; sunk effort; bilinear matching; inner participation shares; inner global price BR; location FOCs/cubic and feasible local candidate; symmetric monopoly histories; a mathematically certified multiple-equilibrium off-path price counterexample.
UNSOLVED: complete Nash PRICE equilibrium correspondence for all arbitrary locations (especially close colocation and negative candidate fees); public/private LOCATION global best responses under a fully specified off-path Nash PRICE continuation; full SPNE; unconditional entry and public establishment policy thresholds.
Status: **STAGE 13 = CONDITIONAL LOCAL CANDIDATE + RIGOROUS OFF-PATH OBSTRUCTION; NOT SPNE CERTIFIED.**
Next stage: select or microfound public fee rule, or characterize the whole price equilibrium correspondence and construct/prove at least one full location SPNE. Do not change g/private entry/matching or silently freeze a local solution.

Reproduction: code/stage13_endogenous_location_audit.py
