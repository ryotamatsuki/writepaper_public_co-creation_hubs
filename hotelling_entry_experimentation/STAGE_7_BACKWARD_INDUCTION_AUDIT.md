# Stage 7 — Backward-Induction Equilibrium Ledger and Independent Audit
Date: 2026-10-27
STATUS: **RESTRICTED SPNE (ONE EXPLICIT PRICE CONTINUATION), CONSISTENTLY RECOMPUTED**.
This is not an equilibrium uniqueness result or a claim of novel economic mechanism.

## Objective
Resolve user concern: evaluate **all** necessary payoffs, price/location equilibria, private entry, social welfare, government setup and policy scope by backward induction, instead of plugging mixed-duopoly formulas into monopoly/no-entry histories.

## 0. Fully specified economic environment
Potential users mass M>0, type x uniform on [0,1], utility when served by:
 G: v-p_G-t(x-a)^2,
 P: v+delta_n-p_P-t(x-b)^2 if a pilot was established, otherwise v-p_P-t(x-b)^2.
Users may stay out for utility 0. Zero variable operating costs. Public G maximizes total surplus; private P maximizes profits. Both prices nonnegative. Private entry fixed cost F_P; public installation fixed cost F_G.
One pilot partner is **separate from the continuum of end users**; its net private experimentation payoff from effort e is r_n e-e^2/2; administrative spillover B*n*e; program matching operating cost k*n. Pilot scope n in {1,2}; r_n=alpha*n/[1+sigma(n-1)]. This treatment avoids counting the same firm twice. Spillover to PRIVATE hub users delta_n=gamma*n*e_n; it is a productivity/utility increment, not a transfer.
Timing: Government choose g and conditional n; pilot firm chooses e (when g=1); private firm decides entry; active facilities choose a,b; active facilities set prices; end users select facility/opt-out.
Selection: upon private indifference enter; government indifference choose build; on off-path all-P allocation pick the LOWEST feasible nonnegative G price Nash equilibrium. Other price continuations exist and may produce different SPNE.
Parameter domain: v>=4t and for EACH chosen n, d_n=delta_n/t in [(sqrt(5)-2)/4, 1/4).

## 1. Backward induction
### Step 5: user demand, price game
For a<b, private cutoff:
 z=(a+b)/2+(p_P-p_G-delta_n)/[2t(b-a)].
Public chooses efficient assignment, so p_G=p_P when both have positive demand. Private FOC yields
 p_G=p_P=2t(b-a)*(1-z_eff),
 z_eff=(a+b)/2-delta_n/[2t(b-a)].
When private efficiency-dominates all users, there are MULTIPLE price Nash equilibria: p_G=q>=max(0,2tD-Delta), p_P=q+Delta, Delta=delta_n-tD(a+b) for b>a; symmetric opposite orientation. Pick the minimum feasible q, with precise piecewise definitions in Stage 6.
At a=b select p_G=0,p_P=delta_n with P-favoring ties. v>=4t ensures full coverage at ALL selected off-path price equilibria; P price deviations are disciplined by G coverage.

### Step 4: service-characteristic equilibrium
For both facilities active under selected continuation:
 a*=1/4-d_n; b*=3/4-d_n;
 z*=1/2-2d_n, private share=1/2+2d_n;
 p_G=p_P=t/2+2delta_n;
 private GROSS operating profit = M*(t+4delta_n)^2/(4t);
 public fee revenue = M*(t^2-16delta_n^2)/(4t).
Stage 6 provides an analytic proof that each facility is a GLOBAL best response over [0,1] under the low-price continuation in the stated domain. It does not claim uniqueness.
For PUBLIC ONLY, optimize a=1/2 and select p_G=0 (many other W-equivalent prices exist).
For PRIVATE ONLY without pilot, optimize b=1/2 and p_P=v-t/4 when v>=4t (full coverage). Global proof: for each p, the centermost b maximizes the mass of users within distance sqrt((v-p)/t); with b=1/2, revenue (v-t*r^2)*min(1,2r) is maximized at r=1/2 for v>=4t.

### Step 3: exact firm payoff accounting
Pilot firm's PRIVATE payoff at optimal investment: r_n^2/2 (not to be conflated with platform's profit).
Private hub gross operating profit under public and P entry: Pi_n=M*c_n, c_n=(t+4delta_n)^2/(4t).
Private hub NET profit if enters: M*c_n-F_P. If stays out: 0.
Private-only gross profit under g=0: M*(v-t/4); net after entering M*(v-t/4)-F_P.
No private hub profit arises when the candidate does not enter.
Public hub revenue is a transfer in welfare (public price is not a separate resource benefit).

### Step 2: entry
Given build and n, P enters iff M>=M_n=F_P/c_n=4tF_P/(t+4delta_n)^2.
Without build, P enters iff M>=M_0=F_P/(v-t/4).
At equality P enters by convention.

### Step 1: pilot investment
e_n = r_n =alpha*n/(1+sigma(n-1)).
delta_n=gamma*n*r_n.
H_n=r_n^2/2 + B*n*r_n - k*n.
H_n includes the pilot firm's PRIVATE net payoff and the administrative benefit; private user's extra spillover delta_n is included separately in usage welfare.

### Step 0: government welfare comparisons across ACTUAL entry paths
Define A(M)=M*(v-t/12).
W_0(M) = A(M)-F_P IF private-only enters; otherwise 0.
W_n(M) = A(M)+H_n-F_G IF the private operator stays OUT;
W_n(M) = A(M)+H_n-F_G + M*c_n/4-F_P IF the private operator enters.
These are FOUR DIFFERENT market structure branches: no hub, private-only, public-only, public-private. Do not use (a*,b*,p*) from MIXED duopoly in public-only or private-only branches.
Municipality chooses argmax{W_0,W_1,W_2} with tie-break specified.
All-users consumer surplus + public price receipts + private GROSS operating receipts + direct H minus fixed fees equals the welfare formula.

Public-private gross use surplus:
W_GP,usage=M*(v-t/48+delta_n/2+delta_n^2/t).
Gross welfare increment against optimal centrally located public monopoly:
DeltaS_n=M*(t+4delta_n)^2/(16t)=M*c_n/4.
SOCIAL desirability of private entry (conditional on fixed n) requires M*c_n/4>=F_P, FOUR TIMES the market mass required for PRIVATE profitability in the specified equilibrium.
This 4:1 ratio is a model-specific overentry result, not inherently new.

## 2. INDEPENDENT symbolic and exact rational verification
Executed independently in Python using fractions.Fraction, sympy, scipy (2026-10-27). Sympy integral of aggregate mismatch utility yielded the W_GP expression EXACTLY; private gross profit =4*DeltaS_n EXACTLY. For each calibration n=1,2, grid of 20,001 deviations in [0,1] for EACH player plus scipy bounded optimizer recovered candidate a*,b* and no profitable unilateral location deviation under selected price continuation. Stage 6 also contains its own price-game BR grid; analytic arguments are required, numerical grids alone are not certification.

Calibration:
 t=1, v=4, alpha=.09, sigma=1.4, gamma=1, B=.5, k=.001, F_P=.55, F_G=.4.
n=1: e=.09, delta=.09, H=.04805; conditional mixed a=.16,b=.66,z=.32,p=.68, public revenue/M=.2176, private gross/M=.4624, entry threshold 1375/1156=1.18944636678; socially efficient entry threshold ~4.75778546713.
n=2: e=.075, delta=.15, H=.0758125; conditional mixed a=.10,b=.60,z=.20,p=.80, public revenue/M=.16, private gross/M=.64, entry threshold 55/64=.859375; socially efficient entry threshold=3.4375.
Private-only threshold (g=0): 11/75=.1466666667.
Government policy threshold INITIAL: 15561/188000=.08277127659574468.
Government policy reentry as M grows: 5187/2560=2.026171875.

Representative equilibria and welfare (proper realized market structures):
M=.10: no-build => no firm, W0=0; best public scope 2, P stays out, G location=.5 and pG=0, W2=.06747917.
M=.20: no-build => P enters alone, b=.5 pP=3.75, W0=.23333333; best public scope 2, P stays out, G=.5 pG=0, W2=.45914583.
M=1.00: no-build => P enters alone, W0=3.36666667; build n=1 => public ONLY (private net hypothetical -.0876), a=.5,pG=0,W1=3.56471667; build n=2 => both enter, a=.1,b=.6,p=.8,P net profit=.09,W2=3.20247917. Therefore choose build n=1.
M=1.20: no-build P only W0=4.15; build n=1 both enter W1=3.93677; build n=2 both enter W2=4.0178125. Therefore DO NOT BUILD.
M=2.10: no-build P only W0=7.675; build n=2 both enter W2=7.6868125. Therefore build n=2.

Interiors of optimal policy regions, with specified tie-breaking:
(0,15561/188000): no public;
(15561/188000,55/64): public scope 2;
(55/64,1375/1156): public scope 1;
(1375/1156,5187/2560): no public;
(5187/2560, infinity): public scope 2.
At each threshold, private indifference can cause jumps; clarify realized tie outcome. The no-g private entry threshold 11/75 lies within the broad public-scope-2 policy region, without changing the optimum.

## 3. Strict scientific caveats
- Existence is supported ONLY under a specified low-price continuation, NOT for arbitrary price Nash equilibria, and uniqueness remains unresolved.
- Scope n is selected by government, **not independently by public and private operators**; the pilot firm is one actor, not n realized separate investing firms. This is NOT yet the full original Stage 1 matching and firm investment specification.
- Pilot net welfare H_n does not scale with M, giving a positive direct pilot value even at negligible user mass. This is a substantive normalization that must be microfounded/replaced to interpret regional thin markets.
- In the current model the public hub directly creates H_n and delta_n by assumption; there is no explicit subsidy/procurement/private platform alternative. Thus the policy regions are NOT proof of a necessary physical public hub.
- Overentry and regime-switch nonmonotonicity likely reflect known mechanisms. Exact absorption/novelty has not been certified.
- If the pilot firm is drawn from the same user mass M, the separate-agent surplus formulation would need adjustment to prevent double counting.
- All welfare assumes quasi-linearity and no marginal public-funds cost; public fee revenue is a transfer. Multiple potential private entrants and endogenous network investment are not modeled.
- Claimed policy thresholds are for the SPECIFIC selection-dependent SPNE, not robust general predictions.
- This Stage is an *audit and consolidation*, not a new theorem.
