# Stage 32 — MEMBER-OWNED VENTURE COMMERCIALIZATION, FOUR-ENTRANT PURE SPNE
Date 2026-10-27. Actual detailed proofs and exact Fraction Python outputs provided in same ChatGPT conversation as artifacts STAGE_32_MEMBER_OWNED_VENTURES_BACKWARD_INDUCTION.md, stage32_member_owned_venture_entry.py and stage32_16_ownership_profiles.csv.

## User-corrected market actors
A urban / B rural, each with INDEPENDENT local government G_i and independent LOCAL private hub entrant P_i. Simultaneous binary establishment by four actors, no cross-region hub construction, G and P co-locate, total 16 configurations. Product incumbents n_A=5,n_B=2 continue positive Cournot output with no venues. Member founders N_A=3,N_B=1 are existing *venture-development companies*, owners of the profits of ventures they co-create, but distinct from incumbent homogeneous good sellers (important simplification). One matched regional A/B pair innovates a DIFFERENTIATED new product competing with n_i incumbents in each market. Different pairs serve separate independent niche markets so project welfare sums and private commercialization income not counted multiple times.

## Stage 5 differentiated-goods Cournot
Final-goods consumer utility in town i: a_i Q_i+a_vi y_i−(b_i Q_i²+2χ_i Q_i y_i+d_i y_i²)/2 +money. Let A_i=a_i−c_i, D_i=a_vi−c_vi, H_i=2b_i d_i(n_i+1)−n_i χ_i²; a patent success yields old-good q_i=[2d_i A_i−χ_i D_i]/H_i and new-good y_i=[b_i(n_i+1)D_i−n_iχ_i A_i]/H_i. Without success old incumbent q_i0=A_i/[b_i(n_i+1)]>0. JV revenue Π=Σ_i d_i y_i²; old-goods CS+resident incumbent profit change M_i (normally negative). ω_i>0 private JV share sum to 1.

## Stage 4 matched founders' R&D
Expected private net after matched project p(e)Π−k e²/2, p(e)=p0+ζe and p<1. e*=ζΠ/k. Result p*=p0+ζ²Π/k, R=p*Π−k(e*)²/2=p0Π+ζ²Π²/(2k)>0. **OWN joined company per-match value** ω_i R. Local social impact from one match S_i=p*M_i+ω_iR, global S=Σ_i S_i. S_i can be NEGATIVE while ω_iR>0 and global S>0 when product cannibalization M_i is large. If M_i<0, S_i<0 iff ω_i<−p*M_i/R.

## Stage 3 joining/revenue consistency
Local firm mass N_i, θ uniform[0,1], G publicly fee zero, P charges f_i and offers incremental specialty θΔ_i. Extra base participation payoff u_i=v_i+ρw_i h_j. An own local firm, independent of G or P ownership, obtains matched JV private return a_i m_j with a_i=ρμ N_j ω_iR. No hub outside option is incremental zero (pre-existing production survives). Its utilities (0,G,P)=(0,u_i+a_i m_j,u_i+a_i m_j+θΔ_i−f_i).
m_i=0 with no local hub; m_i=1 with G or BOTH; m_i=clip(1+(u_i+a_i m_j−f_i)/Δ_i) if P only.
Incremental new matches L=ρμ N_A N_B m_A m_B, assume ρμ N_j≤1. Aggregate business founder profits = L(ω_A+ω_B)R=LR; no duplicate reward.

## Stage 2 price Nash (correct bilateral feedback)
If P only and other region 0 or G or BOTH (foreign member mass fixed), optimal posted fee f_i=(Δ_i+u_i+a_i m_j)/(2), and member share m_i=(Δ_i+u_i+a_i m_j)/(2Δ_i).
With local BOTH, total membership1 and P specializes top-half θ: fee f_i=Δ_i/2, customer share1/2.
With both regions P only, define κ_i=a_i/Δ_i, t=κ_Aκ_B, z_i=1+u_i/Δ_i. Under small-network sufficient restriction (Δ_i−a_i)(1−t)>u_i+a_i, unique equilibrium in BOTH jurisdictions:
  f_i=Δ_i[(2−t)z_i+κ_i z_j]/(4−t),
  m_i=[(2−t)z_i+κ_i z_j]/[(1−t)(4−t)].
PROOF: simultaneous membership fixed point m_i=z_i+κ_i m_j−f_i/Δ_i, so dm_i/df_i=−1/[Δ_i(1−t)], own FOC f_i=Δ_i(1−t)m_i. Solve two linear equations (2−t)m_i−κ_i m_j=z_i. Check every own-price deviation globally via partner fee interior and piecewise-linear self demand. Outside weak-network region, uniqueness isn't asserted.

## Stage 1 and second-best planner
Local welfare W_i=W_i0 + Z_i(s_i)−FG_i g_i−FP_i p_i+ L S_i; Z_i=0 if 0, N_i u_i if G, N_i[u_i m_i+Δ_i(m_i−m_i²/2)] if P only, N_i[u_i+3Δ_i/8] if BOTH. Hub fees are transfers. Local private profit =N_i f_i m_i−FP_i if P-only, N_iΔ_i/4−FP_i if BOTH, otherwise0.
Define fully recomputed *all downstream equilibrium* sign changes:
 DG_i0(s_j)=W_i(G,s_j)−W_i(0,s_j); DG_i1(s_j)=W_i(B,s_j)−W_i(P,s_j).
 DP_i0(s_j)=π_Pi(P,s_j)−π_Pi(0,s_j); DP_i1(s_j)=π_Pi(B,s_j)−π_Pi(G,s_j).
Local pure 2x2 NE states:
  0: DG0<=0, DP0<=0;
  G: DG0>=0, DP1<=0;
  P: DG1<=0, DP0>=0;
  B: DG1>=0, DP1>=0.
Both towns must simultaneously satisfy for 4-entrant PURE subgame-perfect equilibrium (all lower subgames unique under restrictions). The planner maximizes W_A+W_B over all sixteen states with the same downstream behavior. Unlike Stage28, if the OTHER town is P-only, changing one area's G can also change foreign P's OPTIMAL FEE and membership: NEVER hold foreign mass constant when computing DG.

## Key PUBLIC co-entry theorem
Fix foreign town j=G or B with m_j=1 and local town i=P only. Let m_i=(u_i+a_i+Δ_i)/(2Δ_i), 0<m_i<1. Local G add causes m_i→1 (other town mass unchanged). Define ΔL=ρμN_i N_j(1−m_i)>0 and
 H_i=N_i[u_i(1−m_i)+3Δ_i/8−Δ_i(m_i−m_i²/2)]−FG_i.
Local municipal marginal value of G add = H_i+ΔL S_i; national marginal social welfare=H_i+ΔL(S_A+S_B).
Therefore if S_i<0<S_A+S_B and −ΔL(S_A+S_B)<H_i<−ΔL S_i, G refuses to add but planner prefers co-located G+P. Private return to local participant ω_iR>0 meanwhile. Note conditional local discrepancy is an application of regional spillovers; general novelty vs Cremer 1997, Greenlee/Cassiman 1999, Goyal/Joshi 2003, Rochet/Tirole 2003 not certified.

## Exact Fraction all 16 ownership profiles counterexample
Use Stage30 market primitives: A n=5,N=3,u base v=.15,w=.22,FG=.09,FP=.65; B n=2,N=1,v=.10,w=.45,FG=.0405,FP=.24; a=2,c=1,b=d=1,χ=.9,Δ=1, D=1 for both. rho=.4,μ=.25,p0=.05,ζ=.005,k=.0004,ω_A=.8,ω_B=.2.
Cournot venture Π≈.110660844653, e≈1.383260558,p≈.05691630279, private pair profit R≈.005915724187. M_A≈−.01579743681,M_B≈−.03077500469; S_A≈+.00383344765,S_B≈−.00056845465, S_total≈+.00326499300. All probability and cap restrictions satisfied.
If venture social returns exist but venue users are *passive*, private P tariff ignores the commercial return: unique pure 4-player NE=(A BOTH,B BOTH), planner=(BOTH,BOTH).
When venture founders THEMSELVES own commercialization net profits and JOIN EXPECTING THEM: unique pure NE=(A BOTH,B P ONLY), planner=(BOTH,BOTH).
W_NE=2.4046961742, W_SO=2.4050350535, gap .00033887925505. For local B public G add next to existing P, with A both:
 private B tariff fB≈.64017747, participation mB≈.64017747; ΔL≈.10794675848;
 H_B≈−.00001356616; local municipal ΔWB≈−.00007492899; global planner ΔW≈+.00033887926.
With commercial returns NOT in member utility (passive project investors), local municipal ΔWB≈+.0000386069>0 so the sign FLIPS upon giving local members ownership of the project profits. For BOTH P-only pricing, closed-form fees fA≈.61915141655,fB≈.64010982813.
All 16 states' 4 unilateral flips and all planner scores enumerated with Python stdlib fractions.

## Explicit limits / future agenda
Membership good-type firms themselves hold JV stakes but not the n_i old-good incumbents; a founder who also sells old goods must subtract its OWN post-success incumbent loss from participation benefit. Each project in different niche; no overlapping innovation/product rivalry among newly established ventures. ρ and μ exogenous, quality θ affects member venue value but NOT the quality of startup match or product; public budget financed lump sum; fees transfer local. Research matching/effort is per pair and the fraction ω fixed (joint contracting assumed). Baseline firms continue Cournot without any hub. Multiple platform/member price equilibria outside small-network parameter restrictions cannot be ruled out. No direct-public-vs-private-subsidy comparison, and results are NOT literature-verified publication novelty. 
