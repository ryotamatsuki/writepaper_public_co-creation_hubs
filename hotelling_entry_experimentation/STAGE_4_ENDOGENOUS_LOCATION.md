# Stage 4: Endogenous locations, entry and government establishment (2026-10-27)

## Verdict
**CONDITIONAL ANALYTICAL RESULTS + REPRODUCIBLE NUMERICAL DIAGNOSTICS; GLOBAL SUBGAME-PERFECT EQUILIBRIUM NOT CERTIFIED.**

This builds on the Stage 3 fixed-location benchmark. Crucial difference: when a single facility serves all potential participants, it chooses the center x=1/2. Old fixed-endpoint numerical policy comparisons must NOT be transferred unchanged to endogenous-location outcomes.

## Game (restricted specification)
1. Government chooses build g in {0,1}; conditional on build, pilot matching scope n in {1,2}, at cost F_G.
2. Pilot partner invests e_n before the private entry choice; investment return is r_n e - e^2/2; effort e_n=r_n; private-hub user spillover delta_n=gamma*n*e_n, where r_n=alpha*n/[1+sigma(n-1)].
3. A single prospective private operator decides whether to pay F_P to enter, observing g,n and pilot results.
4. Active public and private hubs choose service characteristics a,b, respectively. Public maximizes local total welfare; private maximizes operating profits.
5. Active hubs simultaneously choose prices.
6. Users with x uniform on [0,1] choose their preferred available service or the outside option. When both hubs exist, focus first on the interior, full-coverage pricing branch.

The pilot's direct social surplus is H_n=r_n^2/2+B*n*r_n-k*n. This is additive in the present diagnostic, not proportional to market mass M. Importantly, there is no exogenous upstream pilot when g=0.

Utility when both active:
U_G=v-p_G-t(x-a)^2
U_P=v+delta_n-p_P-t(x-b)^2.

## Interior price-game derivation for a<b
The consumer cutoff given prices:
z=(a+b)/2+(p_P-p_G-delta)/(2*t*(b-a)).
Socially efficient marginal reallocation yields
z_eff=(a+b)/2-delta/(2*t*(b-a)).
Thus public best response in the interior requires p_G=p_P.
Private pricing FOC gives p_P=2*t*(b-a)*(1-z_eff).
This is conditional on both hubs serving positive demand, no outside-option demand, and the appropriate second-order conditions.

## Endogenous-location first-order conditions
Under the same restricted price-game solution:
dW/da=0 implies a=z/2.
Private dPi/db=0 gives
1 + a/2 -3*b/2 -delta/[2*t*(b-a)]=0.
Solving yields the candidate
a=1/4-delta/t,
b=3/4-delta/t,
z=1/2-2*delta/t,
p_G=p_P=t/2+2*delta.

Necessary strict interior domain: 0<=delta<t/4.
The local second derivatives at the candidate are:
W_aa = - (t-4*delta)*(3*t+4*delta)/(4*t), excluding multiplier M;
Pi_bb = - (3*t-4*delta)*(t+4*delta)/(2*t), excluding M.
Both are strictly negative on the stated domain.

**Certification boundary:** These are local FOC/SOC statements. One still must solve the full price game for all off-path a,b including near-colocation, all-private allocation, partial coverage, and public price indifference. The present note does not assert that the displayed candidate is a global location Nash equilibrium or that its entry outcome is a globally valid SPNE.

## Operating profit and welfare of marginal private entry
At the stated interior candidate:
Pi_P^gross = M*(t+4*delta)^2/(4*t).
The public-only hub, conditional on no private entry, locates at x=1/2 and delivers W_G^gross=M*(v-t/12)+H_n.
When both hubs are active, pre-fixed-cost surplus is
W_GP^gross=M*(v-t/48+delta/2+delta^2/t)+H_n.
Thus incremental gross welfare of the private hub, when compared with a single centrally located public hub:
DeltaS_P=M*(t+4*delta)^2/(16*t).
In this restricted candidate, Pi_P^gross = 4*DeltaS_P.
Private entry criterion F_P<=Pi_P^gross. Socially desirable entry criterion F_P<=DeltaS_P.
This is a particular manifestation of familiar socially excessive entry, not a novel theorem.

## No-government private counterfactual
With V=2, t=1, single privately owned monopoly with free position chooses b=1/2 and maximizes revenue at price p=V-t/4=1.75; all users participate. Private entry when 1.75*M>=F_P; no-build social welfare then M*(V-t/12)-F_P, otherwise zero. This differs sharply from the fixed-endpoint monopoly comparator used in Stage 3.

## Validated numerical illustration (not a global equilibrium)
Parameters:
t=1, v=2, alpha=.09, sigma=1.4, gamma=1, B=.5, k=.001, F_P=.55, F_G=.4.

n=1: e=.090, delta=.090, H=.04805; (a,b,z)=(.16,.66,.32), price=.68, private profit/M=.4624, DeltaS/M=.1156, entry iff M>=1.18944637.
n=2: e=.075, delta=.150, H=.0758125; (a,b,z)=(.10,.60,.20), price=.80, private profit/M=.64, DeltaS/M=.16, entry iff M>=.859375.

For market masses M and welfare evaluated in the restricted branch:
| M | Best policy under the restricted branch |
| ---: | --- |
| 0.10 | No public hub |
| 0.20 | Public hub, n=2, private stays out |
| 0.50 | Public hub, n=2, private stays out |
| 0.90 | Public hub, n=1, private stays out |
| 1.00 | Public hub, n=1, private stays out |
| 1.20 | No public hub (private-only) |
| 1.50 | No public hub (private-only) |
| 2.10 | Public hub, n=2, private enters |

Illustrative threshold sequence, conditional on the branch:
- public build starts at approximately M=0.16914130;
- private-only entry threshold with no government is M=0.31428571;
- n=2 induces private entry at M=0.859375;
- n=1 induces private entry at M=1.18944637;
- if both n yield private entry, public establishment under n=2 becomes favorable again at M=2.026171875.
Ties at exact thresholds require a tie-breaking convention.

A grid search over ALL a in [0,1] for the PUBLIC actor, assuming welfare-efficient buyer assignment and the other firm's b held at its candidate, found the displayed a maximizes welfare for n=1,2. This is a useful unilateral allocation check but NOT a complete pricing-subgame/global equilibrium proof. It does not certify private deviations that lead to boundary price equilibria.

## Crucial limitations
- Stage 3's alpha=.3, sigma=1.4, gamma=1 gave delta in {.3,.5}, OUTSIDE delta<t/4, so the endogenous-location interior branch CANNOT replicate Stage 3 using those parameters.
- A single private entrant and perfect information are imposed.
- Pilot partner count n is NOT a realized endogenous matching network or Nash-bargaining outcome.
- Value H_n is independent of M; this is a strong scaling choice.
- The spillover delta is generated by pilot investment, but does NOT depend on service distance |a-b|; Zhang & Li (2013) already study location-dependent R&D spillovers.
- Monopoly with a public facility and monopoly without one must have their own location and price solutions; they cannot be computed as fixed-endpoint monopolies.
- Pricing outside-option and near-colocation regimes, price constraints, private matching scope, and global location SPNE remain unverified.
- Public provision versus equivalent subsidies, contracting, private matching infrastructure, and digital alternatives remains unresolved.

## Prior art and priority
Exact close precedents: Cremer, Marchand & Thisse (1991); Zhang & Li (2013); Matsushima & Shinohara (2014); Mankiw & Whinston (1986); Kang & Kim (2022); Kim (2024).
Sources: 
- https://doi.org/10.1016/0167-7187(91)90004-5
- https://doi.org/10.1007/s00168-013-0556-2
- https://doi.org/10.1016/j.jebo.2014.07.014
- https://doi.org/10.2307/2555627
- https://doi.org/10.1007/s00712-022-00790-y
- https://doi.org/10.1111/manc.12467

**Next gate (Stage 5):** Exhaustive piecewise price best responses with outside options, all location orientations, the full public/private deviation map, pure-strategy SPNE existence or precise non-existence, and a counterfactual with government subsidies. Only after that should this be used to claim local-government establishment criteria.
