# Stage 28 — Solved parsimonious four-player simultaneous local hub entry with preexisting Cournot competition
Date: 2026-10-27. STATUS: MODEL PROPOSITIONS AND PROOFS VERIFIED IN RESTRICTED SETTING; PUBLICATION NOVELTY NOT ESTABLISHED.

## 1. Game and primitive assumptions
Urban A and rural B. Four independent simultaneous binary entrants: local government G_i (maximize OWN region's consumer surplus+resident producer profit+resident hub consumer surplus+resident hub operator profit, less real local hub fixed cost), independent for-profit local intermediary P_i (maximize OWN tariff profit net cost), for i=A,B. Operators CANNOT establish in another territory, but cross-region cooperation is allowed. Each region can have neither, G only, P only, or BOTH; 16 global entry configurations.

For each region n_i>=2 pre-existing identical Cournot product manufacturers, and separate mass N_i>0 of heterogeneous potential hub-project users. Product demand price a_i-b_i Q_i, baseline cost c_i^0<a_i. Urban has n_A>n_B,N_A>N_B; potential further asymmetry in demand. Local existence indicator h_i=max{g_i,p_i}∈{0,1}. Hub enables PROCESS TECHNOLOGY SPILLOVER to existing manufacturers:
 c_i=c_i^0-h_i(d_i+rho xi_i h_j)
with d_i,xi_i>=0, c_i^0>=d_i+xi_i. Thus Cournot markets LIVE even at h_A=h_B=0.

Potential incremental hub users type theta~Uniform[0,1], local mass N_i. With local h_i=1, baseline hub common project return u_i=v_i+rho w_i h_j, 0<v_i, w_i>=0 and Delta_i>v_i+w_i for all rho∈[0,1]. A public G site gives incremental utility u_i FREE (fee=0). Private P site gives u_i+theta Delta_i−f_i with endogenous scalar f_i>=0; no hub user outside option incremental zero, despite ongoing product-market activity. Hub-project users are distinct from incumbent manufacturers, so project benefit and product-market CS/profit are not counted twice. Local public & private can coexist. If no local hub, no incremental project service, not zero product activity. Fixed resource costs FG_i,FP_i>0; fee revenue is LOCAL transfer without lump-sum public tax distortion.
Timing (1) FOUR entry actors move SIMULTANEOUSLY; (2) private P operators (if active) choose posted local fee; (3) uniform hub users choose public/private/none; (4) incumbent n_i Cournot firms choose quantities; (5) payoffs. Solve backwards.

## 2. Solved price, participation, product subgames
For fixed h_i,h_j and producer surplus margin A_i=a_i−c_i^0>0:
 q_i^*=[A_i+h_i(d_i+rho xi_i h_j)]/[b_i(n_i+1)].
 Γ_i=n_i(n_i+2)/(2 b_i(n_i+1)^2)=[1−1/(n_i+1)^2]/(2b_i).
 W_i^{production}=Γ_i[A_i+h_i(d_i+rho xi_i h_j)]².
 No-hub W_i^0=Γ_i A_i²>0. This counts CS and resident manufacturing profits.
 Let z=h_j and T_i(z)=Γ_i{[A_i+d_i+rho xi_i z]²−A_i²}, u_i(z)=v_i+rho w_i z.
 When P alone and G absent, tariff f_M=(u_i(z)+Delta_i)/2; private subscriber share (u+Delta)/(2Delta), gross fee receipts N(u+Delta)²/(4Delta); gross project utility N*3(u+Delta)²/(8Delta). When P & G coexist, P tariff f_D=Delta_i/2, subscriber fraction 1/2, gross fee receipts N Delta/4; public takes rest; project utility N(u+3Delta/8). When G alone all N attend G at zero charge and enjoy u. Price best responses are GLOBAL given Delta>u, including 100% demand, partial demand, zero-demand corner cases.

Net LOCAL welfare increments beyond W_i^0:
 X_i(0,z)=0;
 X_i(G,z)=T_i(z)+N_i u_i(z)−FG_i;
 X_i(P,z)=T_i(z)+3N_i(u_i(z)+Delta_i)²/(8Delta_i)−FP_i;
 X_i(B,z)=T_i(z)+N_i[u_i(z)+3Delta_i/8]−FG_i−FP_i,
 B denotes BOTH local hubs.
Private NET intermediary profit if alone pi_i(P,z)=N_i(u_i(z)+Delta_i)²/(4Delta_i)−FP_i; if coexisting pi_i(B,z)=N_iDelta_i/4−FP_i; otherwise 0.

## 3. Theorem 1 (NE iff conditions, all sixteen histories)
Define unilateral deviations at fixed foreign active status z:
 DGi0(z)=T_i(z)+N_i u_i(z)−FG_i,
 DGi1(z)=N_i[u_i(z)/4−3u_i(z)²/(8Delta_i)]−FG_i,
 DPi0(z)=N_i(u_i(z)+Delta_i)²/(4Delta_i)−FP_i,
 DPi1(z)=N_iDelta_i/4−FP_i.
Then local (g_i,p_i) is best response:
 0 iff DGi0<=0 and DPi0<=0;
 G iff DGi0>=0 and DPi1<=0;
 P iff DGi1<=0 and DPi0>=0;
 B iff DGi1>=0 and DPi1>=0.
A FOUR-player stage1 PURE Nash exists at profile (s_A,s_B) exactly when the local condition holds in BOTH regions with z_A=h_B,z_B=h_A. The selected downstream price, user, Cournot subgames are all unique, so each first-stage pure NE lifts to a pure SPNE. This 2×local-2x2 characterization is equivalent to testing all 16 cases. Finite stage1 game admits mixed NE, not generally analyzed.
Proof: expressions are exact payoffs under fully solved continuations. Each binary actor has a single unilateral alternative; sign inequalities exhaust 4 deviations; conversely all signs imply no unilateral gain. Subgame-perfect continuation optimality already established.

## 4. Theorem 2 (global constrained planner closed form)
Define M_i(z)=max{X_i(G,z),X_i(P,z),X_i(B,z)}. An integrated planner who selects ONLY hub establishment vector but respects decentralized continuation equilibrium achieves
 W*=W_A0+W_B0 + max{0, M_A(0), M_B(0), M_A(1)+M_B(1)}.
Proof: partition 16 states by two binary indicators (h_A,h_B): (0,0),(1,0),(0,1),(1,1). All states within last class have z_A=z_B=1, so welfare maximization over local ownership is additively separable; same holds first three classes. Max the four category scores. This is a second-best planner, not output/price first best.

## 5. Proposition 3 (sharp private underentry relative to socially desirable P co-location)
Conditional on a G_i already being present, P_i earns incremental NET profit N_iDelta_i/4−FP_i, but its presence increases CONSOLIDATED social welfare by 3N_iDelta_i/8−FP_i. The product spillover h_i and cross spillover to other region remain unchanged, so comparison is EXACT. If
 N_iDelta_i/4<FP_i<3N_iDelta_i/8,
private P will NOT build beside G even though the national planner would ADD P to G. Equivalent thickness interval 8FP/(3Delta)<N_i<4FP/Delta. The fixed-cost interval width is N_iDelta_i/8. This is a known user-surplus appropriability wedge; NOT automatically a novel field-journal theorem. Full entry equilibrium depends on all other players, but rational example in §8 proves feasibility of UNIQUE pure NE/planner divergence.

## 6. Proposition 4 (nonmonotone free public co-location incentive)
If P_i and a FOREIGN hub already exist, the own municipal incentive to add FREE G_i is
 DGi1(rho)=N_i[u_i/4−3u_i²/(8Delta_i)]−FG_i, u_i=v_i+rho w_i.
 dDGi1/drho=N_i w_i[1−3u_i/Delta_i]/4.
 For w_i>0, second derivative=−3N_i w_i²/(4Delta_i)<0; peak at rho=(Delta_i/3−v_i)/w_i, if interior.
 The highest pre-fixed-cost gain is N_i Delta_i/24. G enters over incumbent P exactly for u in [u_minus,u_plus], when FG<N Delta/24, where
 u_±=(Delta±sqrt(Delta²−24Delta FG/N))/3,
intersected with feasible u range; otherwise never enters (except a knife-edge).
Proof: welfare P-only vs B difference; quadratic formula/differentiation. Nonmonotonicity arises because improving private standalone service coverage reduces benefits of a free public alternative and crowding replaces differentiated value. Must separately check private P's own simultaneous entry deviation (could exit). Conditional finding does not imply nonmonotonic optimum of PUBLIC-ONLY first build.
For rural B parameters §8: u_B=.1+.45rho, FG_B=.04,N_B=Delta_B=1: u interval [4/15,2/5], rho∈[10/27,2/3].

## 7. Proposition 5 (cross-border underprovision conditional on active other site)
Suppose region i empty, other j active with owner state s_j∈{G,P,B}. Adding G_i yields local municipal benefit DGi0(1). Its NATIONAL marginal welfare gain is DGi0(1)+E_j(s_j), where
 E_j=[T_j(1)−T_j(0)]+N_j[V_{s_j}(v_j+rho w_j)−V_{s_j}(v_j)],
 with V_G(u)=u, V_P(u)=3(u+Delta_j)²/(8Delta_j), V_B(u)=u+3Delta_j/8.
 E_j>=0, >0 when rho>0 and w_j>0 or xi_j>0. If -E_j<DGi0(1)<0, first public hub is not locally desirable but IS nationally desirable holding foreign site state fixed. This is classic Cremer–Marchand–Pestieau 1997 fiscal-spillover logic, not novelty.
Proof via additive regional welfare and positive derivatives; local govt does not internalize E_j.

## 8. Exact rational 16-history example with unique NE and socially suboptimal ownership-count
Urban A:
 n=5,N=3,a=2,b=1,c0=1,d=.05,xi=.1,v=.15,w=.22,Delta=1,FG=.09,FP=.65.
Rural B:
 n=2,N=1,a=1.7,b=1,c0=1,d=.05,xi=.18,v=.1,w=.45,Delta=1,FG=.04,FP=.30.
At rho=.4, checking ALL FOUR independent entry deviations for ALL SIXTEEN states using exact fractions:
  unique PURE NE=(s_A,s_B)=(B,G), i.e. urban co-located public+private, rural free public only;
  unique global second-best social optimum=(B,B), BOTH public+private at both sites;
  no-hub Cournot production gross SW total=35/72+49/225=0.7038888889;
  FULL welfare NE=2.216852611111..., planner=2.291852611111..., exact welfare gap 3/40=.075.
The rural P entry gain under G is 1/4−.3=−.05, while socially it adds 3/8−.3=.075. B municipality's G addition to P under h_A=1 at rho=.4 returns +.0006>0; thus rural P-only is NOT another NE. Other 3-player deviations excluded by exact rational enumeration (script).
At rho=.6 same unique NE(B,G) and SO(B,B), gap .075. At rho=.0 pure NE include (P,G),(P,P) and planner (P,P); rho=.2 NE (B,G),(B,P), planner (B,P); rho=.8 and 1 NE (B,G),(B,P) and planner (B,P).
Proof for computed illustration by exact Fraction exhaustive 16-profile payoff vectors and four flip deviations; avoid general claim of NE uniqueness for all parameters.

## 9. Competition-thickness lemma
Γ_i=[1−(n_i+1)^(-2)]/(2b_i), with n treated continuously ∂Γ/∂n=1/[b_i(n+1)^3]>0. At unchanged product-market scale A_i, direct+cross efficiency δ>0, Cournot regional W gain T_i=Γ_i(2A_iδ+δ²) rises with incumbent competitor count. This is an aggregate surplus statement only (individual producer profit may differ), not a universal claim more competitor density warrants more public spending.

## 10. Scientific limits & prior art
This is deliberately ONE-SIDED uniform-taste differentiated hub user market, not two-sided I/R endogenous matching. No direct producer's decision to use a hub, no endogenous R&D, no Nash bargaining, and no match formation probabilities. The product efficiency improvement upon any hub establishment is exogenous reduced-form and separate from hub users' additive project values. This is a BENCHMARK not a complete microfoundation of co-creation innovation. Cross externality rho affects both project value and manufacturer productivity but not privately chosen hub interoperability. Both public and private can supply innovation equally; private has user-type service specialty. Venue costs fixed only; local governments count local platform profit and firms; no social cost of public funds. Do not treat existing public-free assumption as optimal-fee proof.
Prior models: d'Aspremont & Jacquemin (1988) product Cournot+R&D; Cremer et al. (1997) two local public goods; Long & Wong (2009) two government/private port service tiers; Kim (2024) mixed two-sided market; Zhu, Li, Valavi, Iansiti (2021) platform interconnectivity+entry. Neither generic underentry nor positive cross externality is novelty. Need newer nonabsorbed theorem with endogenous links/participation and public vs subsidy/procurement instrument comparison before publication claim.
