# Stage 9 — Full Model Backward-Induction: Analytical Solution Boundary and Verified Nonuniqueness
2026-10-27. Status: **BARGAINING / INVESTMENT / CONTACT CHOICE SOLVED CONDITIONALLY; TWO-SIDED PARTICIPATION MULTIPLE; FULL SPNE NOT SOLVED**.

## Scope
Do not reuse Stage 6/7 unilateral quadratic-Hotelling formulas for this full two-sided, capacity-constrained matching, relationship-specific investment, spillover and public/private entry game. This is not a reduced model.

## 1. Bargaining conditional on realized m (m=1,2)
V(m,e)=A m^rho e, C_m(e)=c d_m e^2/2, d_1=1. In simultaneous independent bilateral Nash bargaining with no renegotiation of other agreements:
 partner payment T_m(e)=(1-beta)A[m^rho-(m-1)^rho]e.
 innovator's coefficient r_m=A[m^rho-m(1-beta)(m^rho-(m-1)^rho)].
 Therefore r_1=beta A, r_2=A[2(1-beta)+(2beta-1)2^rho].
 Unique sunk investment e_m=r_m/(c d_m), private innovator net payoff u_m=r_m^2/(2c d_m);
 per-partner Nash receipt tau_m=(1-beta)A[m^rho-(m-1)^rho]e_m;
 project social net surplus S_m=[A m^rho r_m-r_m^2/2]/(c d_m).
 **Resource accounting: S_m=u_m+m*tau_m exactly**.

These results are conditional on fixed realized matched partner count m, given the contractual disagreement convention. Do NOT confuse realized partner count m with chosen requested proposals k or hub matching limit n.

## 2. Innovator endogenous proposal count (conditional on q)
If per-proposal completion probability q, q is not exogenous:
 U0=0
 U1=q*u1-kappa
 U2=2*q*(1-q)*u1+q^2*u2-2*kappa.
For n=2, k=2 iff U2>=U1 and U2>=0, that is
 q*u1+q^2*(u2-2*u1)>=kappa,
 2*q*(1-q)*u1+q^2*u2>=2*kappa.
k=1 iff U1>=U0 and U1>=U2, k=0 iff both alternatives <=0. For n=1 restrict to k=0,1.
Within a hub all innovator types have the SAME payoff differences across k, since x affects only the hub mismatch term; absent indifference all innovators at a given hub choose the same k. This follows from the original full-model assumptions; no homogeneous-contact assumption is added.

## 3. Exact partner receipts / knowledge for given k and q
per innovator expected partner transfer b(k,q):
 b(0,q)=0;
 b(1,q)=q*tau1;
 b(2,q)=2*q*(1-q)*tau1+2*q^2*tau2.
per innovator reusable knowledge z(k,q):
 z(0,q)=0;
 z(1,q)=q*e1;
 z(2,q)=2*q*(1-q)*e1+q^2*2^rho*e2.
Then for hub j:
 S_j=N_j^I*k_j;
 L_j=R_j[1-exp(-{1-exp(-lambda R_j)}*S_j/R_j)] if R_j,S_j>0, else zero;
 q_j=L_j/S_j if S_j>0, else 0;
 B_j=N_j^I*b(k_j,q_j);
 Z_j=N_j^I*z(k_j,q_j);
 delta_j=sum_h gamma_{hj} Z_h/M;
 innovator quality Phi_j=max_{k=0..n_j}U_k(q_j)+delta_j;
 partner quality Psi_j=B_j/N_j^R if N_j^R>0, else defined zero.
Social project surplus per innovator with k=2:
 2q(1-q)S_1+q^2 S_2, and with k=1: q S_1.

## 4. Coupled two-sided market fixed point
For both hubs active, positions a<b, full coverage and interior shares, define z_I=N_G^I/M, z_R=N_G^R/R_total.
 z_I=(a+b)/2+[p_P^I-p_G^I+Phi_G-Phi_P]/[2 t_I(b-a)].
 z_R=(a+b)/2+[p_P^R-p_G^R+Psi_G-Psi_P]/[2 t_R(b-a)].
 N_G^I=M*z_I, N_P^I=M*(1-z_I);
 N_G^R=R_total*z_R, N_P^R=R_total*(1-z_R).
Together with S,L,q,B,Z,delta and contact optimality, these constitute a piecewise transcendental FIXED-POINT SYSTEM. Do not apply cutoff formulas when some users opt out, service positions cross, a hub is absent, or a side has zero membership; the general exact equation is membership as the integral of the indicator that its utility dominates all alternatives and zero outside.

For fixed upstream tariffs/locations there need NOT be a UNIQUE solution of this system. Therefore downstream Nash payoffs and upstream fee/location/entry thresholds need equilibrium selection or an equilibrium correspondence; there is no single generic closed-form vector for the full SPNE under the current specification.

## 5. Constructive counterexample to uniqueness in a single-hub fixed-price subgame
Consider a permissible history with one active hub at l=0, scope n=1, M=R_total=1, lambda=4, A=c=1, beta=1/2, rho=1/2, kappa=.02, d_1=1, v_I=5, p^I=0, t_I=.1, v_R=.02, p^R=.10, t_R=.1, gamma_{hj}=0. Innovators of any x join the hub because of high base utility.

Ex post m=1: r_1=.5, e_1=.5, u1=.125, tau1=.25, S1=.375.

(A) ZERO-activity membership Nash equilibrium:
  R_h=0, all innovators choose k=0, S_h=0, q=0, L=0.
  Innovators will not unilaterally propose because no partner exists and kappa>0; partners will not unilaterally join because their base net value is -.08-.1*y^2<0 and there are no proposals. Thus every agent best responds.

(B) POSITIVE-activity membership Nash equilibrium:
  All innovators choose k=1. Let r denote partner mass (those with types 0<=y<=r).
  h(r)=1-exp(-4r);
  q(r)=r[1-exp(-h(r)/r)], because S=1.
  Partner average transfer=0.25 q(r)/r.
  Their membership cutoff must satisfy:
   r=sqrt((0.25*[1-exp(-(1-exp(-4r))/r)]-.08)/.1).
  The right minus left side is positive at r=.90 and negative at r=.95, so by continuity at least one root exists in (.90,.95). Numerically r=.915277703237985, q=.599592299833587, L=.599592299833587.
  Innovators: q*u1-kappa=.0549490374791984>0, hence k=1 strictly preferred.
  Partner type y=r is indifferent, participants y<r strictly prefer participation and others prefer nonparticipation. Therefore this second Nash equilibrium is valid.

The same fixed tariffs, service position and all parameter values support two distinct equilibria; uniqueness is FALSE already in a strict lower subgame. This does NOT imply full SPNE nonexistence, and it does NOT mathematically prove non-expression in every special function, but refutes a unique universal closed-form equilibrium claim absent selection/refinement.

## 6. Upstream backward-induction correspondence (not a solved scalar expression)
Let E(g,z) denote the SET of participation/matching equilibria for each history and choices z=(locations, scopes, tariffs).
- Private operating profit in an equilibrium e in E:
  pi_P^op(e)=p_P^I*N_P^I+p_P^R*N_P^R-f_P(n_P)-c_m*L_P.
  Entry net profit = pi_P^op - F_P.
- Public gross welfare:
  W = integrated baseline benefits of participating innovators/partners
      - taste mismatches - proposal resource costs
      + integrated expected project net social surpluses S_{K_i}
      + integrated reusable external benefits delta_j + xi sum Z_j
      - all hubs' real scope and matching costs
      - fixed private entry cost F_P (if entered)
      - fixed public installation cost F_G (if installed).
  Nash transfers and fees are NOT social resources.
- Stage 3 private price conditions on a smooth unique selected e branch:
  N_P^s + sum_{r=I,R} p_P^r * (dN_P^r/dp_P^s)
        - c_m*dL_P/dp_P^s=0, s=I,R.
  Public price choices instead maximize W; FOCs may not exist at corners or participation jumps.
- Stage 2 positions and scope are GLOBAL best replies of pi_P^op and W after solving all downstream price/participation subgames; discontinuities and multiplicity matter.
- Stage 1 private entry is decided from equilibrium-conditioned expected operating profits relative to F_P.
- Stage 0 public build requires equilibrium-conditioned comparison W(g=1) >= W(g=0), respecting realized private entry and facilities. If multiple continuations, thresholds are a correspondence and depend on equilibrium selection.

## 7. Scientific status and next decision
**NOT SOLVED IN GENERAL / NOT THEOREM-FROZEN.** Bargaining, investment and contact decision inequalities are exact. Multiple lower-subgame equilibria and coupled exponential participation equations block unique closed-form backward induction for general parameters.
Next options for actual rigorous full-game results:
(A) maintain full two-sided matching and derive an equilibrium CORRESPONDENCE (existence, extremal equilibria, parameter-dependent closed-form bounds, and verified numerical SPNE);
(B) retain all strategic margins but alter the matching rule / add a microfounded coordination mechanism, explicitly declaring this a modified—not identical—model.
Before choosing either, the original current model is now properly diagnosed; no novel public-facility theorem has been proven.
Bibliographic context: Matsushima & Shinohara (2014), JEBO 106:428–441, doi:10.1016/j.jebo.2014.07.014; Kim (2024), Manchester School 92:281–295, doi:10.1111/manc.12467.
