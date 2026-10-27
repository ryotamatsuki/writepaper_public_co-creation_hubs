# Stage 14 — Arbitrary-location fee-subgame correspondence and global SPNE obstruction
Date: 2026-10-27.
Status: **PIECEWISE DEMAND FULLY CLASSIFIED UNDER FULL COVERAGE; EXPLICIT PURE PRICE-NE FAMILY; PURE PRICE-NONEXISTENCE CERTIFICATE FOR AN OFF-PATH FULL-COVERAGE SUBGAME; ORIGINAL FINITE-v SPNE STILL OPEN.**
Do not write "global SPNE proved" or claim that the finite-v game has no SPNE without a separate opt-out branch proof.

## 0. Maintained model, clearly distinguishing different claims
Stage 12/13 retains: (i) public establishment; (ii) one prospective private entrant; (iii) simultaneous endogenous positions a,b∈[0,1] if both active; (iv) simultaneous nonnegative fees for TWO sides; (v) endogenous two-sided membership and matched links μ N^I_j N^R_j; (vi) ex-post bilateral bargaining and relationship-specific investment. For analytic subgame certificates specialize M=R=1, t_I=t_R=t, beta=2/3, A=c=1, so innovator and partner receive u=tau=w=2/9 per match and total match net surplus S=4/9. Let r=μw/t.

*Full coverage* below is the standard no-outside-option subgame (or a valid pricing region with high v and an active alternative). The original game has an outside option zero, v=4 in its illustrative calibration. Strictly distinguish the former's nonexistence certificate from a theorem about the latter.

## 1. Full-coverage allocation for arbitrary positions and fees
For a != b define d=|b-a|, midpoint m=(a+b)/2, and orientation c:
  c=1-m if b>a (P on right); c=m if b<a (P on left).
Let private shares x∈[0,1] of innovators and y∈[0,1] of partners, and public fees g_I,g_R≥0, private fees p_I,p_R≥0.
Let k=r/d, A=c+(g_I-p_I)/(2td)-k/2, B=c+(g_R-p_R)/(2td)-k/2.
Subscription Nash allocations are EXACT fixed points
   x=clip(A+k*y,0,1); y=clip(B+k*x,0,1).
The exhaustive combinations (some combinations coalesce at boundaries):
  II:  x=(A+kB)/(1-k²), y=(B+kA)/(1-k²), requiring 0<x,y<1 and k≠1.
  I0:  x=A∈(0,1), y=0, B+k*A<=0.
  I1:  x=A+k∈(0,1), y=1, B+k*(A+k)>=1.
  0I:  x=0, y=B∈(0,1), A+k*B<=0.
  1I:  x=1, y=B+k∈(0,1), A+k*(B+k)>=1.
  00:  x=y=0, A<=0, B<=0.
  01:  x=0,y=1, A+k<=0,B>=1.
  10:  x=1,y=0, A>=1,B+k<=0.
  11:  x=y=1, A+k>=1,B+k>=1.
At k=1 interior fixed points can form a continuum if A+B=0, subject to box. At d=0, taste-location differences cancel and full-side indifference/coordinated tipping must be separately defined; price Nash examples given below.
For d>r, k<1 yields a unique full-coverage subscription Nash equilibrium via contraction. For d<r, multiple subscription Nash equilibria can occur **AT THE SAME FEE PROFILE**; the fixed-point correspondence, not a single demand curve, is the primitive.

For finite v and arbitrary possibly unaffordable fees, these nine FULL-COVERAGE branches are not exhaustive: add a no-participation set to each population, determined exactly as the measurable argmax of {0,U_G,U_P}. An exact finite-v full-game price correspondence remains unresolved.

## 2. A whole family of pure two-sided fee Nash equilibria (both sides tip to P)
Write c above for fraction of types toward P's end, and Δ=t[r-2d*(1-c)]. This is the maximum P-vs-G posted fee premium compatible with all-P in each side; for P on the right Δ=t[r-d(a+b)].

Suppose P is closer to 1/2 than G (c>1/2) and 0<d<2r. Then the socially optimal FULL-COVERAGE assignment of both sides is entirely to P (because network matching efficiency is supermodular and W(s,h) is convex in the symmetric share s for d<2r, with P all outperforming G all).

For any q satisfying
   q>=0,
   q+Δ>0,
   if r<d<2r, additionally q+Δ>=2t(d-r),
and such that public at price q still covers the whole population when a private fee deviation is made,
the symmetric fee vector
    g_I=g_R=q,
    p_I=p_R=q+Δ
is a pure-strategy fee Nash equilibrium with x=y=1 and P gross profit 2(q+Δ). In particular, for 0<d<r and c>1/2, Δ>0, so all sufficiently low nonnegative q work.

Proof of P global best reply:
For any attainable P shares x,y (with public fallback full coverage) its gross fee revenue is bounded by
  Q= A_q*(x+y)-2td*(x²+y²)+4tr*x*y,
  A_q=q+2td*c-tr.
Set s=(x+y)/2, h=(x-y)/2. Then
  Q=2A_q*s+4t(r-d)*s²-4t(r+d)*h².
If d<=r, Q is convex (or linear) in s, with max at s∈{0,1}; since Q(1,1)=2(q+Δ)>0, all P is a global max.
If r<d<2r, Q is strictly concave in s, and its maximum over [0,1] is s=1 iff q+Δ>=2t(d-r). The h² term is always nonpositive. The chosen prices attain Q(1,1).
G already attains maximal social welfare by assigning both groups fully to the more central P, so cannot improve via any own-price deviation.
=> The fee equilibrium is analytically certified against ALL joint side-fee deviations, not just local FOCs.

If a=b, another family exists under the tie rule that sends all indifferent agents to P: g_I=g_R=q, p_I=p_R=q+tr, within the public-full-cover range, giving P profit 2(q+tr). At same locations any single matched hub is planner-efficient. Subscriber tie-selection must be explicitly stated; this demonstrates multiplicity rather than a uniquely defined continuation.

## 3. Exact rational witness: pure fee Nash NONEXISTENCE in a full-coverage off-path regime
This result is for the NO-OUTSIDE-OPTION full-coverage game, not an unconditional assertion for the original finite-v game.

Take
 t=M=R=1, r=1/45, a=99/200=.495, b=21/40=.525.
Then d=3/100=.03, m=51/100=.51:
 r<d<sqrt(3)*r<2r,
 G is closer to 1/2 than P.
Define
 A0=2d*(1-m)-r=323/45000 >0.
 s0=A0/[4(d-r)]=323/1400≈.230714.
 s_bar=d*(2m-1)/(2r-d)=27/650≈.0415385.
 s_zero_price=A0/(4d)=323/5400≈.0598148 >s_bar.
Here s=(x+y)/2 and h=(x-y)/2.

(A) Any best reply by P to any public nonnegative fee vector has s>s_bar:
 - If x=1 or y=1, then s>=1/2>s_bar.
 - If P has any zero fee on one side, fixed-point participation implies the private share on that side >=A0/(2d); hence s>=s_zero_price>s_bar.
 - If both private fees are positive and both shares interior, the global P FOCs for its strictly concave inverse-demand envelope Q imply
   x+y=[2A0+g_I+g_R]/[4(d-r)],
   so s>=s0>s_bar.
 - Any zero-P-demand side cannot be a P best reply: by setting a small strictly positive fee on that side it attracts a positive fringe because A0>0 and the other private side only reinforces it, without reducing the revenue earned on the other side.
Thus s>s_bar in every private best reply, including price corners.

(B) For ANY full-coverage allocation, relative to all-P welfare:
 W(x,y)-W(1,1)
 =2t(1-s)[d*(2m-1)-(2r-d)*s]
   -2t(d+2r)*h².
For s_bar<s<1 this is STRICTLY negative. If s=1 then h=0 and W equals W(1,1).

(C) In the no-outside full-coverage model, the public player can always raise both nonnegative own fees sufficiently to induce all consumers to choose P at fixed P fees, achieving W(1,1). Therefore no fee pair yielding s<1, with P playing a global best reply, can be a public best reply.

(D) If x=y=1 and P is globally optimal, its side prices must reach their all-P cutoffs p_s=g_s+Δ and, by global P optimality, satisfy
  p_s>=2t(d-r)=7/450≈.0155556.
At g_I=g_R=0, those fixed P prices strictly exceed A0=323/45000≈.0071778, causing ALL users to choose G in the UNIQUE subscription fixed point (since d>r). Public obtains W(0,0)>W(1,1) because G is closer to the center. This is a profitable G price deviation.
=> No pair of pure fee strategies forms a Nash equilibrium for this full-coverage location profile. This is a rigorous obstruction in the full-coverage specialization, not merely a failure of candidate first-order conditions.

Note that Stage 13's original finite-v=4 game retains an outside option. Raising public fees may produce OPT-OUT if P's fees are exorbitant. An additional global upper-bound/corner audit would be required before extending (C) to the exact finite-v game. Do NOT make that extension without proof.

## 4. What this implies for the requested global SPNE construction
- An SPNE in pure strategies must prescribe a Nash continuation in EVERY price subgame, including off-path locations. In the **full-coverage no-outside specialization** above, the explicit location with no pure fee Nash rules out pure-strategy SPNE of the ENTIRE location-then-fee game, regardless of on-path location cubic roots.
- In the CURRENT finite-v=4 game, we have an exact nine-branch characterization whenever all agents are served, and a provably nonunique all-P fee-NE family. The full finite-v correspondence (additional opt-out branches), existence/nonexistence across ALL off-path a,b, and a pure or mixed full-game SPNE are NOT yet settled.
- Hence an assertion that Stage 13's interior cubic and entry/setup thresholds are globally certified is NOT justified.
- A mixed-strategy continuation is the appropriate mathematical direction if pure pricing fails at an off-path node, but mixed existence and measurable off-path selection must be proved, not assumed.
- Another legitimate MODEL CHANGE is a publicly committed access fee or modified price timing/insulating tariffs; this changes the game and cannot be hidden inside "equilibrium selection".
- The economic interpretation is important: spatial gains compete with two-sided matching concentration; public welfare maximization can want one-sided concentration while private price best replies insist on partial fragmentation. Coordination failure can destroy pure pricing equilibrium.

## 5. Academic context
Armstrong (2006 RAND, doi:10.1111/j.1756-2171.2006.tb00037.x), Ambrus & Argenziano (2009 AEJ:Micro, doi:10.1257/mic.1.1.17), Fernández-Ruiz (2020 Letters in Spatial and Resource Sciences, doi:10.1007/s12076-020-00249-y), and Jullien/Pavan/Rysman (2021 HIO overview) cover competing two-sided market structures, allocation multiplicity and equilibrium difficulties. There is no certified independent novelty from the existence/nonexistence phenomenon.
Reproduction: code/stage14_price_subgame_regimes.py
