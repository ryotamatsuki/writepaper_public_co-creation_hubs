# Stage 30: Differentiated new product from cross-regional joint R&D — exact four-entrant result
2026-10-27. **NEW-TO-BASELINE RESULT, NOT a publication novelty certification.** Full local reproducible report and code exported to conversation files: STAGE_30_DIFFERENTIATED_NEW_PRODUCT_EQUILIBRIUM.md; stage30_differentiated_product_hubs.py; stage30_16_profiles.csv.

## Model: corrected 4 independent actors
A urban, B rural, municipality G_i and independent private hub provider P_i in each. All 4 decide build 0/1 at once, own-region hubs ONLY; GP collocation allowed => 16 sites states (s_A,s_B)∈{0,G,P,BOTH}². Governments maximize OWN local welfare, private intermediaries maximize own membership revenues less own fixed costs. Positive Cournot incumbent goods output/profit without hubs. Local hub members are DISTINCT from the final-goods incumbents and JV investors; representative venture commercialization rewards do not directly enter the members' utility (major explicit simplifying assumption).

Timing backward: 4-player simultaneous siting -> P chooses scalar posted fee, G free -> people select nonparticipant/G/P -> joint A/B R&D venture optimizes effort -> success draws differentiated new good, afterward n_A and n_B incumbent Cournot competitors and new venture producer compete in EACH product market.

## Ex-post new differentiated product
In market i, incumbent old good Q_i=n_i q_i, new good y_i; quasilinear representative utility a_i Q_i+a_vi y_i−(b_i Q_i²+2χ_i Q_i y_i+dvi y_i²)/2, χ_i²<b_i dvi. Let A_i=a_i−c_i>0, D_i=a_vi−c_vi>0.
Cournot after commercialization with n_i identical resident incumbents and ONE joint venture new-good supplier:
H_i=2b_i dvi(n_i+1)−n_i χ_i².
q_i^1=(2dvi A_i−χ_i D_i)/H_i, y_i=(b_i(n_i+1)D_i−n_iχ_i A_i)/H_i.
No innovation: q_i^0=A_i/[b_i(n_i+1)]>0.
CS_i^1=[b_i(n_i q_i^1)²+2χ_i n_i q_i^1 y_i+dvi y_i²]/2.
incumbent total π_i^1=n_i b_i(q_i^1)².
Joint venture operating profit Π_E=Σ_i dvi y_i²; residential ownership shares ω_A=ω,ω_B=1−ω.
Per-success LOCAL social value relative pre-success local Cournot:
V_i=ΔCS_i+Δπ_inc_i+ω_i Π_E; V=Σ V_i.
R&D financing cost shares ω_i.
Prob successful p(e)=p0+η e, η=η0+ρ λ m_A m_B; conventional η0>0,p0>0 even without hubs, effort cost k e²/2. Joint venture optimizes p(e)Π_E−k e²/2, giving e*=η Π_E/k and p*=p0+η²Π_E/k, with interior global success probability cap.
Local R&D change coefficient K_i=(Π_E/k)(V_i−ω_iΠ_E/2), global K=Σ K_i=(Π_E/k)(V−Π_E/2); K_i can be NEGATIVE despite gross-success V_i>0 due to competition/ownership/research-cost incidence.

## Endogenous local paid/private fee and attendance (Stage28 unchanged)
Uniform users θ∈[0,1], mass N_i, no hub incremental U0=0, public free UG=u_i, private UP=u_i+θΔ_i−f_i; u_i=v_i+ρw_i 1{foreign active}, 0<u_i<Δ_i. If P alone optimal fM=(u_i+Δ_i)/2, participation m_i^P=(u_i+Δ_i)/(2Δ_i), fee gross N_i(u_i+Δ_i)²/(4Δ_i), real user service value 3N_i(u_i+Δ_i)²/(8Δ_i). If G+P, fD=Δ_i/2, P/G half each, TOTAL attendance m_i=1, gross fee N_i Δ_i/4, gross user service N_i(u_i+3Δ_i/8). G only: all users attend m_i=1, benefit N_i u_i. No hub m_i=0. Let Y_i(state,z) be local net real user value minus all local fixed hub costs.

Exact regional continuation welfare W_i=W_i^Cournot,0 +Y_i(state,z)+p0 V_i +K_i η(state_A,state_B)². Private fee operating profit unchanged by venture commercialization in arms-length operator setup.

## THEOREM 1: complete pure first-stage SPNE
At other-region state s_j, define government unilateral gains
D_Gi^0=N_i u_i−FG_i + K_i[η(G,s_j)²−η(0,s_j)²],
D_Gi^1=N_i[u_i/4−3u_i²/(8Δ_i)]−FG_i +K_i[η(BOTH,s_j)²−η(P,s_j)²].
Private unconditional gross fee gains (given local G=0 or G=1, respectively)
D_Pi^0=N_i(u_i+Δ_i)²/(4Δ_i)−FP_i, D_Pi^1=N_iΔ_i/4−FP_i.
Local state 0 NE BR iff D_Gi^0≤0,D_Pi^0≤0; G iff D_Gi^0≥0,D_Pi^1≤0; P iff D_Gi^1≤0,D_Pi^0≥0; BOTH iff D_Gi^1≥0,D_Pi^1≥0. **For the 4-player stage-1 pure NE these inequalities must hold for BOTH regions simultaneously**, with subgame fee, attendance, JV R&D and Cournot continuations. Proof: exactly exhaust all four independent players' one binary deviation.

## THEOREM 2: planner and composition mismatch
Second-best planner controls entry site flags but respects decentralized fee/membership/R&D/Cournot continuation. Let W_base=Σ_i W_i0+p0 V+Kη0². With M_i(0)=max_{s_i∈{G,P,BOTH}}Y_i(s_i,0):
W^SO=W_base+max{0, M_A(0),M_B(0),max_{s_A,s_B ∈{G,P,BOTH}}[Y_A(s_A,1)+Y_B(s_B,1)+K(η(s_A,s_B)²−η0²)]}.
Cannot reduce both-active branch to M_A(1)+M_B(1) because interaction m_A m_B depends on ownership through total attendance. Nine both-active states reduce to 2×2 coverage types since G and BOTH give attendance 1 and P alone gives t_i<1.

## THEOREM 3: new-to-Stage28 public free co-entry condition
Fix P_i existing, the other region has positive members m_j, hence h_i remains=1 when G_i adds alongside it. Define a=ρλ m_j>0, t_i=(u_i+Δ_i)/(2Δ_i)∈(1/2,1),
E_i=(η0+a)²−(η0+a t_i)²=a(1−t_i)[2η0+a(1+t_i)]>0.
Base Stage28 local GOV marginal co-entry H_i=N_i[u_i/4−3u_i²/(8Δ_i)]−FG_i.
New extended local GOV marginal: H_i+K_i E_i.
Global planner marginal: H_i+K E_i.
If K_i<0<K and −K E_i< H_i <−K_i E_i, OWN municipality declines free G co-entry yet aggregate second-best planner desires it. Unlike Stage28, new venture formation does NOT cancel when adding G to preexisting P: extra free site expands total local *participation*, thus bilateral match feasibility. Private operator P's own tariff entry thresholds D_P stay as baseline.

## Exact 16-history numerical counterexample
Urban A n=5,N=3,a=2,c=1,b=dv=1,D=1,χ=.9,Δ=1,v=.15,w=.22,FG=.09,FP=.65.
Rural B n=2,N=1,a=2,c=1,b=dv=1,D=1,χ=.9,Δ=1,v=.1,w=.45,FG=.0405,FP=.24.
rho=.4, λ=.05, η0=.005, k=.0004, p0=.05, ω_A=ω_B=.5. Exact rational solver computes all 16 profiles, each actor's four unilateral BR checks, ALL second-best planner alternatives. Probability min .0569,max .22291 (feasible). Post-success JV profit Π_E=.11066084465; per-success regional V_A=.039532985516, V_B=.024555417635, both POSITIVE; K_A=+3.2832448351,K_B=−.8603309463,K=+2.4229138888.
Stage28 benchmark with PRODUCT innovation disabled, and PROCESS innovation d=xi=0: unique pure NE=unique planner=(BOTH,BOTH), W=2.40405555556.
Stage30 product innovation ON: unique pure NE=(BOTH,P), planner=(BOTH,BOTH); W_NE=2.40792765175, W_SO=2.40877429689, gap .000846645144.
For P_B existing while A both, u_B=.28,m_B^P=.64, η_BP=.0178,η_BB=.025, E_B=.00030816:
old B-gov co-entry H_B=+.0001,
new B-gov co-entry H_B+K_B E_B=−.0001651195844,
planner co-entry H_B+K E_B=+.000846645144.
Every pure first-stage deviation at (BOTH,P) is strictly loss-making; no other pure NE in all 16 histories; planner (BOTH,BOTH) uniquely best.

## Limits + related literature
This is a minimal one-sided fee hub SERVICE model plus ONE representative venture's stochastic effort, NOT a two-sided endogenous bilateral participant matching or actual founders receiving IP royalties. Neither χ product differentiation nor venture/production market assignment is endogenously chosen; public-free a primitive, government fiscal funding lump-sum, JV funds exactly owned shares; possibility of paid private subsidies/procurement not compared. These can all materially alter ownership equilibria and should be the next robustness test before publication claim.
The mechanism of regional externalities, business stealing and innovation joint ventures already exists in Cremer et al. 1997, Capuano & Grassi 2019, Lambertini et al 2002, Choné & Linnemer 2020, Kim 2024, and Zhu et al 2021. This proves a theorem for OUR extended game and demonstrates Stage28 conditions genuinely change; not yet proof the theorem is UNABSORBED by general existing models.
References: https://doi.org/10.1016/j.ijindorg.2020.102663 ; https://doi.org/10.1080/10438599.2018.1461333 ; https://doi.org/10.1016/S0167-7187(00)00082-5 ; https://doi.org/10.1111/manc.12467 .
