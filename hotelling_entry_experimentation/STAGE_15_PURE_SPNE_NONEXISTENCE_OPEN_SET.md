# Stage 15 — Pure-Strategy SPNE Nonexistence Theorem: Open Set, Not Just a Numerical Witness

Date: 2026-10-27
STATUS: **PROVED for the MANDATORY-FULL-COVERAGE specialization with simultaneous unrestricted two-sided nonnegative fees and endogenous locations. NOT PROVED for original finite-v with opt-out, nor for mixed-strategy SPNE. Novelty NOT ESTABLISHED.**

## Main proposition
Take the symmetric bilateral-matching/mixed-duopoly model of Stage 12/14:
- two unit masses of single-homing agents, independently uniform on [0,1], one side innovators, other side partners;
- active public G is a social-welfare maximizer, active private P a total posted-fee revenue maximizer;
- posted nonnegative fees on each side chosen simultaneously AFTER both hubs select arbitrary a,b in [0,1];
- all users on each side MUST subscribe to one of the active hubs (full coverage NO OUTSIDE OPTION);
- quadratic taste mismatch t>0, positive per-link Nash shares u=tau=w>0 (e.g. beta=2/3), feasible bilinear matching links μ N^I_j N^R_j, r=μw/t>0;
- no private variable operational cost; private/public establishment fixed costs are already sunk when the price subgame is reached.

For any a<b with d=b-a, m=(a+b)/2 satisfying the following STRICT OPEN inequalities
 (i) r<d<2r,
 (ii) 1/2<m<1,
 (iii) 0<2m-1 < [(2r-d)(d-r)]/[d(3d+2r)],
and 0<a<b<1,
the simultaneous nonnegative two-sided FEES subgame has NO pure-strategy Nash equilibrium.

Consequently, in an extensive game with public establishment, subsequent potential private entry, UNRESTRICTED simultaneous location choices and then that fee subgame, **NO PURE-STRATEGY SPNE exists at all**: the offending fee subgame is a proper subgame following feasible (possibly off-path) histories in which public builds, P enters and locations (a,b) are chosen. This is a pure-strategy statement ONLY. It does not rule out mixed-strategy SPNE and does not claim socially relevant on-path nonexistence.

The parameter domain is NOT a knife-edge: all inequalities are strict, hence there is an open set of (a,b,r). It holds even with arbitrarily weak positive network effects: for any sufficiently small r>0, choose d=3r/2 and m=101/200=0.505, because 2m-1=1/100 <1/39 (the RHS at d=3r/2). Positions a=m-3r/4 and b=m+3r/4 remain feasible for r<0.66.

## Proof
Normalize t=1 (the result is invariant to common multiplication of all service payoffs by t). P's fractions of each side are x,y∈[0,1]. Define s=(x+y)/2 and h=(x-y)/2. Let public and private posted fee vectors g=(gI,gR), p=(pI,pR), each in R²_+.

Since k=r/d<1, for each fee vector the full-coverage demand fixed point is unique by contraction.

If x,y are interior, the private inverse demand is
 pI=gI+A0-2d*x+2r*y,
 pR=gR+A0-2d*y+2r*x,
 where A0=2d(1-m)-r.
The same two expressions are UPPER BOUNDS on private feasible posted fees at all demand corners, with equality when that side's private share is interior.
The resulting PRIVATE REVENUE ENVELOPE is
 Q(x,y)=(gI+A0)*x+(gR+A0)*y -2d(x²+y²)+4r*x*y.
Hessian diag -4d and offdiag 4r is strictly negative definite because d>r.

(A) Show every private global best reply, to ANY g≥0, implies
 s>s_bar := d(2m-1)/(2r-d).
By assumption A0>0 (follows from (iii)). A private-side share 0 cannot be optimal: reducing that side's own price to small positive epsilon yields strictly positive clientele even if other private-side membership is zero, while the positive network response of the other side cannot reduce its revenue.
If either side has share 1, s≥1/2>s_bar.
If either private fee is zero, the inverse-demand bound gives the corresponding side share ≥ A0/(2d), hence s≥A0/(4d)>s_bar.
Otherwise 0<x,y<1 and both private fees >0. Interior FOCs for maximizing the strictly concave envelope give
 x+y=[2*A0+gI+gR]/[4(d-r)],
 hence s≥ A0/[4(d-r)]>A0/(4d)>s_bar.
The strict inequalities s_bar<1/2 and s_bar<A0/(4d) both follow from (iii).

(B) PUBLIC WELFARE relative to sending both unit masses to P is exactly
 W(x,y)-W(1,1)
 =2(1-s)[d(2m-1)-(2r-d)s]-2(d+2r)h².
For s_bar<s<1 this is strictly negative. The public player can always make its two own fees arbitrarily high, holding P fees fixed; since subscribers are FORCED to use one hub, this implements x=y=1 and improves social welfare. Thus at a pure fee Nash, P must serve BOTH whole sides: x=y=1.

(C) To serve all P at fee Nash, each private price equals its all-P cutoff:
 pI=gI+Delta, pR=gR+Delta,
 Delta=r-d(a+b)=r-2dm.
If P prices less than these cutoffs on a side, a small price increase raises profit with unchanged share.
A necessary PRIVATE BEST RESPONSE condition is
 pI≥2(d-r), pR≥2(d-r).
For pI,pR>0, this follows by allowing a small reduction in one private-side market share while keeping the other share exactly at one via a compensating (feasible) small discount on the OTHER side. The derivative of the revenue envelope at (1,1) with respect to that share is p_s-2(d-r), which must be nonnegative for a maximizer at the upper boundary.
If one of pI,pR were zero, the other positive, the same compensated deviation is feasible; if both were zero, a simultaneous small price increase would be profitable. Hence the necessary componentwise lower bound is valid.

Since m>1/2, d>r, we have
 2(d-r)>A0=2d(1-m)-r,
because 2dm-r>0.
Therefore if PUBLIC unilaterally changes its TWO prices to (0,0), the UNIQUE full-coverage demand equilibrium has x=y=0 (all-G): at P fixed fees p_s>A0 each side's best-response zero when the other P-side has mass zero; contraction yields uniqueness.
But W(0,0)-W(1,1)=2d(2m-1)>0: G prefers all-G to all-P. This is a profitable public deviation, contradicting price Nash. QED.

## Explicit exact-rational witness
t=M=R=1; β=2/3, A=c=1, thus w=u=tau=2/9. μ=1/10 => r=1/45.
a=99/200=.495, b=21/40=.525; d=3/100, m=51/100=.51.
A0=323/45000=.007177777...
s_bar=27/650=.04153846...
A0/(4d)=323/5400=.0598148... >s_bar;
A0/[4(d-r)]=323/1400=.230714...
2(d-r)=7/450=.0155555... >A0;
W(all-G)-W(all-P)=3/2500=.0012.
Hence an entire OPEN neighborhood of parameters/locations generates pure fee Nash nonexistence.

## Stronger interpretation / caution
In this MANDATORY-FULL-COVERAGE specialization, **no matter how small but strictly positive r is**, unrestricted endogenous locations admit a sufficiently close off-path (a,b) for which the price subgame lacks a pure NE. Thus zero-network-effect benchmark may have pure SPNE, whereas each positive r produces an off-path impossibility for PURE SPNE. This is a discontinuous feature of the pure-SPNE requirement across the entire strategy space, NOT a claim that network instability is necessarily observed on the equilibrium path.

The original Stage 12 game allows NONPARTICIPATION with finite v (e.g. v=4). That is a materially different game. Stage 14 derived its exact participation fixed-point operator but HAS NOT shown no pure fee NE at the example when outside options are allowed. **Do not assert the general full-model SPNE nonexistence until the finite-v price game is proved to lack pure equilibria at a feasible off-path location.**

The mechanism and existence status must be compared to:
- Sanjo (2007), *Hotelling's Location Model with Quality Choice in Mixed Duopoly*, Economics Bulletin: no SPNE in mixed-duopoly location-quality-price sequence.
- Lu (2006), *Hotelling's Location Model in Mixed Duopoly*, Economics Bulletin: changed price-equilibrium location conditions but no SPNE (linear travel).
- Matsumura & Matsushima (2009), *Cost Differentials and Mixed Strategy Equilibria in a Hotelling Model*, Annals of Regional Science: pure-strategy nonexistence, mixed strategies.
- Heywood, Li & Ye (2022), *Mixed duopoly under hotelling with convex production costs*, Annals of Regional Science: nonexistence under simultaneous locations, existence under leadership.
- Roger (2017), *Two-sided competition with vertical differentiation*, Journal of Economics: nonquasiconcave two-sided payoffs, pure equilibrium disappearance, mixed equilibrium existence.
- Kim (2024), *Mixed duopoly in two-sided markets*, Manchester School: public/private two-sided platform competition in a fully covered linear city.

**NOVELTY AS A STANDALONE NONEXISTENCE RESULT: HIGH ABSORPTION RISK.** A worthwhile paper would show new threshold architecture linked to matching effort/hold-up, an on-path policy consequence, or a credible institutional resolution (public tariff commitment or sequential fee move) compared with these precedents, and rigorously handle finite-v opt-out/mixed-equilibrium welfare.

Stage 14 proof was independently cross-checked with exact fractions and symbolic algebra on this date. Regime enumeration and finite-v formulas remain in Stage 14.
