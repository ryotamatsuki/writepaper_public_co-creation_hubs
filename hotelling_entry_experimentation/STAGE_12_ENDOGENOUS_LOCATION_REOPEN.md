# Stage 12 — Endogenizing location in the bilinear two-sided mixed-duopoly model
Date: 2026-10-27
Status: **PRELIMINARY INTERIOR ANALYSIS, NOT GLOBAL LOCATION SPNE**.

## Why this stage is necessary
Stage 11 rigorously proves the SPNE of a restricted GAME WITH EXOGENOUS positions a=1/4,b=3/4. It does NOT certify this pair for a location-choice game. Indeed, the Stage 11 parameterization has a profitable unilateral private location deviation. Thus, before restoring two-matching-partner scopes, endogenous location must be audited and solved.

## 1. Retain complete Stage 11 primitives
M innovators and R partners choose between public G at a and private P at b, quadratic taste losses t_I,t_R, link mass mu*N_j^I*N_j^R, match-surplus partition u, tau, S=u+tau. Two-sided joining fees and private entry/public setup stay endogenous. Here analyze a<b and full-coverage interior demand ONLY.
Let d=b-a>0, m=(a+b)/2, X=N_P^I/M, Y=N_P^R/R, Z=mu*M*R*S, h_I=t_I*d,h_R=t_R*d.

Sorting indifference (full coverage, interior) yields:
  h_I[X-(1-m)]-mu*R*u*(Y-1/2)=-(p_P^I-p_G^I)/2,
  h_R[Y-(1-m)]-mu*M*tau*(X-1/2)=-(p_P^R-p_G^R)/2.

## 2. Public-efficient demand shares, given a,b
Total real W (pre fixed costs and baseline intrinsic values) = -quadratic transport losses on both sides
  + Z*[(1-X)*(1-Y)+X*Y].
Public's fee best response (provided it can implement positive-price efficient shares) induces maximization over X,Y. First-order conditions:
  h_I[X-(1-m)]-mu*R*S*(Y-1/2)=0,
  h_R[Y-(1-m)]-mu*M*S*(X-1/2)=0.
Under q=h_I*h_R - mu²MR*S²>0, unique efficient interior fractions:
  X = 1/2 + (1/2-m)*h_R*(h_I+mu*R*S)/q,
  Y = 1/2 + (1/2-m)*h_I*(h_R+mu*M*S)/q.
The efficient shares coincide at X=Y=1/2 if a+b=1 but generally do NOT equal one-half.

Given these shares, public best-response fees must satisfy
  p_G^I=p_P^I +2*mu*R*tau*(Y-1/2),
  p_G^R=p_P^R +2*mu*M*u*(X-1/2).
For general locations these are NOT equal to private fees, contrary to a naive reuse of Stage 11's symmetric fees.

## 3. Private fee FOCs
On the interior, with rival posted fees held fixed, own-price demand derivatives follow from the inversion of the linear two-sided sorting system.
Let k=h_I*h_R -mu²MR*u*tau >0.
Private first-order conditions imply:
  p_P^I = 2*(h_I*X-mu*R*tau*Y),
  p_P^R = 2*(h_R*Y-mu*M*u*X).
Consequently, the conditional maximized private operating profit is
  Pi_P(a,b)=2*M*h_I*X²+2*R*h_R*Y²-2*Z*X*Y.
These equations are LOCAL to the interior, no-outsider, nonnegative-fee regime. Full noninterior price deviations must be independently certified for a location choice.

## 4. Public location condition (envelope)
When public fee choices achieve globally efficient sorting for given a,b, differentiating W conditional on b yields:
  a* = [M*t_I*(1-X)²+R*t_R*(1-Y)²] /
       [2*(M*t_I*(1-X)+R*t_R*(1-Y))].
This is an endogenous-location FOC, not a complete global best-reply proof.

Private location condition:
  ∂Pi_P(a,b)/∂b=0
with x,y being ENDOGENOUS equilibrium price/subscription responses at every b.
Noninterior b<=a, corner market shares, zero or forbidden negative prices, monopoly histories all require checks.

## 5. Exact counterexample against fixed Stage 11 locations
Parameters M=R=t_I=t_R=1, mu=beta=1/2, A=c=1, v_I=v_R=4:
 u=1/8, tau=1/4, S=3/8, Z=3/16.
At (a,b)=(1/4,3/4):
 X=Y=1/2; p_P^I=3/8, p_P^R=7/16,
 Pi_P=13/32=.40625.
Differentiating the endogenous two-sided fee-equilibrium profit along b, with a held at 1/4, yields ∂Pi_P/∂b=-3/10.
At (a,b)=(1/4, 37/50=.74):
 X=Y=0.508099173553719, p_P^I=.370912396694215,
 p_P^R=.434424793388430,
 p_G^I=.372937190082645, p_G^R=.435437190082645;
 Pi_P=.409191160713066 > 13/32.
The deviating location supports positive fees, strict interior demand and full coverage; the associated fee subgame is globally incentive-compatible under the same quadratic positive-definiteness argument used in Stage 11 with shifted, noncentral efficient optimum (to be detailed in a subsequent globally certified Stage 12 theorem). Thus the Stage 11 pair is NOT a location equilibrium when choices of positions are allowed.

## 6. Symmetric two-sided specialization producing an analytic CUBIC
For M=R=1, t_I=t_R=t, and equal Nash rewards u=tau=w (possible with beta=2/3), set r=mu*w/t>0, d=b-a.
By symmetry public-efficient participation fractions are X=Y=x.
Their interior allocation condition is:
   (d-2r)*(x-1/2)=d*(1/2-(a+b)/2).
Public location FOC: a=(1-x)/2.
Private fee profit: Pi_P=4*t*(d-r)*x².
Private location FOC, with a fixed in unilateral b derivative:
   (d-r)*(1-2d)=r*x.
Combining with public FOC implies:
   x=[d*(1-d)-2r]/(d-4r),
   (d-r)*(1-2d)*(d-4r)-r*[d*(1-d)-2r]=0.
Expanded cubic:
   -2d³+(1+11r)d²-(6r+8r²)d+6r²=0.
For beta=2/3, t=1, mu=1/10, A=c=1, u=tau=2/9, r=1/45, one feasible local root is
   d≈.48770612836,
   x≈.51503394524,
   a≈.24248302738,
   b≈.73018915574.
Other algebraic roots violate interior fractions or feasible positions. This is a LOCAL location equilibrium CANDIDATE. It is NOT yet a global location Nash equilibrium or entry/public setup SPNE.

## 7. Mandatory next gates
1. Prove fee subgame Nash for ALL relevant location deviations, including center crossing, colocation, tipping, monopoly and nonparticipation.
2. Prove both location best responses globally or identify parameter domain or equilibrium selection supporting them; derive all corner equilibria.
3. Recalculate entry with public under endogenous-location Π_P^D; no-public private monopoly remains center under the technology.
4. Recalculate public construction welfare for actual market structure; previously Stage 11 fixed-location thresholds are NOT transferable as location-game results.
5. Restore n∈{1,2} only after endogenous-location benchmark is certified; verify whether the central causal mechanism survives.

No claim of newly discovered research contribution follows from this cubic: public/private spatial differentiation and platform network effects have direct predecessors.
