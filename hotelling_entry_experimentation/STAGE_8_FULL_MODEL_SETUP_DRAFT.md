# Stage 8 — Proposed FULL MODEL setup (NOT a solved equilibrium)
Date 2026-10-27 | Status: SETUP ONLY / NOT THEORY FREEZE

## Research question
Given private hub entry, endogenous platform profiles, two-sided user/partner participation, feasible matching, relationship-specific sunk effort, post-effort bargaining, knowledge reuse, and endogenous matching scope: when should a local authority directly establish a co-creation facility? What specialization and coverage result?

**This supersedes the reduced-form Stage 3-7 game as a modeling PROPOSAL, not as an established derivation.** Stage 6/7 closed-form calculations are benchmarks with different simplifying assumptions. They must not be quoted as equilibria of the game below.

## Actors and populations
- Public service hub G, if local authority establishes it at fixed resource cost F_G.
- One potential privately owned hub P that may enter at fixed resource cost F_P.
- Innovator firms i of total mass M, independently typed x_i ~ U[0,1], outside option 0.
- Potential experiment partners r of total mass R, typed y_r ~ U[0,1], outside option 0. Each partner can serve at most one accepted project.
- Platform service characteristic l_j in [0,1] has mismatch costs t_I (x_i-l_j)^2 for innovators and t_R (y_r-l_j)^2 for partners. There is no built-in public technical advantage. Marginal access prices are (p_j^I,p_j^R)>=0, two-sided.
- Financial transfers (hub charges, ex-post Nash transfers) cancel from social welfare; all physical, effort, proposal, and hub costs count.

## Strategy timing (perfect information about previous public choices; rational expectations of mean-field matching)
0. Municipality chooses g in {0,1}. If g=1 pay F_G and activate hub G.
1. Observing g, potential private hub P chooses enter/no enter; entering costs F_P.
2. Each active hub j simultaneously chooses service characteristic l_j in [0,1] and proposal-slot limit n_j in {1,2}. Public maximizes anticipated total social welfare; private maximizes anticipated operating profits net its fixed cost. A hub also pays scope cost f_j(n_j), with f_j(2)>=f_j(1).
3. Active hubs simultaneously choose two-sided nonnegative access prices (p_j^I,p_j^R). G maximizes welfare, P maximizes profit. No presupposed low-price selection: report the price-equilibrium correspondence when multiple equilibria.
4. Each innovator simultaneously chooses a hub j or outside and requested contacts k_i in {0,...,n_j}; each partner chooses a hub j or outside, accepts at most one project. Joining and k choices are Nash best responses to anticipated aggregate matching probabilities and spillovers. Innovators pay a real proposal/search resource cost kappa*k_i.
5. The matching protocol processes proposals into feasible links. Each innovator i has a realized integer K_i in {0,...,k_i}; each partner takes at most one link.
6. After the complete matched coalition size K_i is observed, the innovator chooses noncontractible relationship-specific effort e_i>=0, at cost c*e_i^2/2.
7. The innovator and its K_i matched partners bargain ex post over the realized project gross value, and output plus distinct re-usable knowledge benefits accrue. All actors receive payoffs.

Stage 4 is a two-sided participation equilibrium with endogenous congestion, so it may be nonunique. Stage 2 is a simultaneous location/scope game and Stage 3 simultaneous tariff game. A full SPNE requires equilibrium strategies following ALL possible histories.

## Matching: explicit proposal/capacity-feasible mean-field protocol
Define membership masses N_j^I, N_j^R and proposed link mass S_j = integral_{i joining j} k_i dmu_I(i).
Let R_j=N_j^R.
Each proposal identifies a feasible contact with probability
 a(R_j)=1-exp(-lambda*R_j), lambda>0;
if successful, it targets one uniformly chosen participating partner.
When a partner receives multiple proposals, it randomly accepts exactly one.
In the Poisson/continuum large-population approximation, completed match mass is
 L_j = R_j*[1-exp(-a(R_j)*S_j/R_j)] for S_j,R_j>0,
and L_j=0 if either is zero.
Proposal success rate q_j = L_j/S_j in (0,1].
Under mean-field sampling, innovator K_i is approximated by Binomial(k_i,q_j) (individual-dependence corrections need verification for exact finite matching). Every completed link consumes one partner capacity: L_j<=min(S_j,R_j).
This matching technology is BOTH supply constrained and sensitive to absolute partner market thickness; holding participant shares fixed while M,R ->0 makes matching rare.
Participation and matching must be jointly solved; q_j is NEVER an exogenous fitted "success probability".

## Ex post bilateral Nash-in-Nash bargaining and innovation effort (PRIMARY specification)
For a project with realized K=m in {1,2} matched partners and sunk innovator effort e, transferable total project value is V(m,e)=A*m^rho*e, A>0, 0<rho<=1, with V(0,e)=0. The investor incurs a resource effort cost C_m(e)=c*d_m*e^2/2, c>0, d_1=1, d_2 in [1,2]; d_2 is an effort-cost scope factor, preserving the possibility of shared or partner-specific quality investment.

Once effort is sunk, the innovator separately and simultaneously negotiates with EACH matched partner using bilateral Nash bargaining, with innovator bargaining weight beta in (0,1) and partner weight (1-beta). If one relation breaks, the remaining m-1 connections deliver value V(m-1,e) and existing transfers to other partners are maintained (specified no-renegotiation disagreement convention). Therefore each partner's marginal bilateral surplus is
   DeltaV_m(e)=V(m,e)-V(m-1,e)=A*[m^rho-(m-1)^rho]*e.
The Nash-in-Nash payment to EACH matched partner is
   T_m(e)=(1-beta)*DeltaV_m(e).
The innovator keeps the residual before paying investment cost:
   R_m(e)=V(m,e)-m*T_m(e)=r_m*e,
   r_m=A*[m^rho-m*(1-beta)*(m^rho-(m-1)^rho)].
The total partner transfers plus innovator pre-cost rents add EXACTLY to gross project value.

The innovator anticipates these separate ex post bargains when it makes its sunk investment. Interior effort:
   e*_m=r_m/(c*d_m), e*_0=0,
assuming r_m>0, checked over the admissible parameters.
Innovator private NET income:
   pi_i_proj(m)=r_m^2/(2*c*d_m).
Each successful partner's Nash transfer:
   T_m(e*_m)=(1-beta)*A*[m^rho-(m-1)^rho]*e*_m.
Real social project value net of investment (bargaining transfers cancel):
   V_social(m)=A*m^rho*e*_m -c*d_m*(e*_m)^2/2
              =[A*m^rho*r_m-r_m^2/2]/(c*d_m).

The key comparison is endogenous, not imposed:
  e*_2 < e*_1 IFF r_2/d_2 < r_1,
  e*_2 > e*_1 IFF r_2/d_2 > r_1.
With A=c=1, beta=.9, rho=.3, d_2=1.8:
 r_1=.9, r_2≈1.18491553, e*_1=.9, e*_2≈.65828641, so broader scope lowers effort.
With beta=.5, rho=.5, d_2=1.8: r_1=.5, r_2=1, e*_1=.5, e*_2≈.55555556, so broader scope raises effort.
Thus the sign is NOT hard-coded. More partners may raise or reduce investment, which is exactly a question to investigate after equilibrium membership and scope are solved. These are comparative statics at FIXED REALIZED K, not endogenously selected scopes.

The SOURCE 2014 Matsushima–Shinohara paper uses simultaneous separate Nash bargaining, whereas the earlier trial multilateral grand-coalition Nash split gives theta_m=beta/[beta+m*(1-beta)]. This alternative MAY be used in sensitivity analyses, but it is NOT the primary specification. The source's precise variable-production-cost technology and outside-option rules should still be audited against this stylized Nash-in-Nash formula.

## Endogenous platform usage
For each i at hub j, expected project net return:
 Uproj_j(k)=sum_{K=1}^k Binom(k,K)*q_j^K*(1-q_j)^(k-K)*r_K^2/(2*c*d_K).
Set Uproj_j(0)=0.

Let total reusable knowledge production
 Z_h=integral_{i joins h} E[K_i^rho*e*(K_i)] dmu_I(i).
Hub-specific value of externally reusable technology:
 delta_j = sum_{h active} gamma_{h,j} * Z_h / M, gamma_{h,j}>=0.
This benefit is distinct from direct bargaining project revenues. It enters innovators' product value (and W) exactly once. Cross-hub term gamma_G,P Z_G/M is public-pilot knowledge spillover to private customers, but is now endogenous to participation, matching and effort. If no G, its spillover contribution is identically zero.

Innovator expected utility:
 U^I_j(x,k) = v_I - p_j^I - t_I*(x-l_j)^2 - kappa*k
               + Uproj_j(k) + delta_j.
An innovator selects (j,k) or no participation 0.

Let aggregate partner bargaining receipts in hub j be:
 B_j=integral_{i joins j} E[K_i*T_{K_i}(e*_{K_i})] dmu_I(i).
If partner assignments are symmetric, average partner revenue is B_j/N_j^R (take 0 at zero partners). Partner expected utility:
 U^R_j(y)=v_R-p_j^R-t_R*(y-l_j)^2+B_j/N_j^R.
Each partner chooses (j or no participation 0).
All N_j^I,N_j^R,S_j,q_j,Z_j,delta_j,B_j and individual k choices form one fixed-point participation/matching equilibrium.

## Platform profit and social welfare
Operational facility resource cost: C_j=f_j(n_j)+c_m*L_j (and no additional costs by default).
Private hub net payoff:
 Pi_P = p_P^I*N_P^I+p_P^R*N_P^R -C_P -F_P if it entered; otherwise 0.
Public hub total social welfare objective:
 W =
   integral_{active innovators}[v_I-t_I*(x_i-l_j)^2-kappa*k_i + delta_j] dmu_I
 + integral_{active partners}[v_R-t_R*(y_r-l_j)^2] dmu_R
 + integral_{active innovators} E[V_social(K_i)] dmu_I
 + xi*sum_j Z_j
 -sum_{active j}[f_j(n_j)+c_m*L_j]
 -g*F_G -1{P enters}*F_P.
Here xi*Z values a DISTINCT nonappropriable policy benefit (e.g. public task learning) not captured in reusable private δ; xi can be set to zero if reuse already measures all additional value.
No fees or Nash bargaining shares are added to social welfare, because they are internal transfers (under quasi-linear utility and marginal cost of public funds 1).
Alternative subsidies, procurement and purely digital networks are POLICY COUNTERFACTUALS, not silently counted as already included choices.

## Required solution/verification sequence (reverse induction)
7. Solve transfers in K=1,2 coalitions.
6. Solve e*(K) and each side's monetary payoff, check resource accounting.
5+4. Given active hubs, l,n,p, solve the coupled membership/scope/matching/spillover fixed point, including outside options and feasible partner capacities; characterize multiplicity.
3. Solve two-sided fee equilibrium for EVERY possible upstream location/scope/entry history.
2. Solve BOTH hubs' global location AND scope deviations, including all noninterior participation/price branches.
1. Solve private entry based on its EXPECTED optimized net profit and no-public counterfactual.
0. Compute all actually realized no-hub/private-only/public-only/GP social welfare values and optimize g (including F_G).
At EVERY step verify optimizer feasibility, corners, global profitable deviations, and conservation of transfers / actual resource costs.

## Key scientific distinctions
The 2014 Matsushima–Shinohara paper chooses 1 or 2 buyers and then uses separate bilateral simultaneous Nash bargaining; the present multilateral alternative is tractable but must not be conflated with the original.
Kim (2024) already studies public/private competition in two-sided Hotelling platforms with investment, and Zhang–Li (2013) already study endogenous R&D/location in a mixed oligopoly. Entry and private-v-social disparity are also known.
It is NECESSARY to run stripped-label ablations and a direct policy alternative test before claiming a new mechanism for why public physical ownership is necessary.

## Stage status
PROPOSED SETUP COMPLETE, EQUILIBRIUM UNSOLVED. This full two-sided model is not nested mechanically by the older single-sided numerical results: the equilibrium objects and matching/bargaining timing have been changed and must be solved afresh.

Primary sources:
Matsushima & Shinohara (2014): https://doi.org/10.1016/j.jebo.2014.07.014
Kim (2024): https://doi.org/10.1111/manc.12467
Zhang & Li (2013): https://doi.org/10.1007/s00168-013-0556-2
