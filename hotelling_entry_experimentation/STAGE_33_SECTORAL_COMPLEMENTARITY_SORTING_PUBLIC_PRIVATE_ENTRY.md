# Stage 33 — Urban R&D × rural manufacturing: endogenous skill-complementary JV innovation, sorting, 4-player ownership entry

2026-10-27. STATUS: **FULL BACKWARD MODEL SETUP, CLOSED-FORM PRICES EXCEPT TWO PRIVATE HUBS (FOUR EQUATIONS), EXACT CONDITIONAL PROPOSITIONS, AND ALL-16-OUTCOME NUMERICAL CHECK.** Publication novelty not verified. Full 20KB Japanese manuscript plus Python script and 32-row CSV generated in chat as Stage33 files.

## Empirical restriction
Do NOT say "rural areas are mostly manufacturing": OECD Enhancing Rural Innovation in Japan says trade/services are 38% of rural employment, manufacturing 19%. Regional manufacturing concentrations are documented by 2025 Japanese Cabinet Office report. Interpret the assumption as *comparative advantage in the venture-participating firms*: A's skill θ_A is R&D, AI, service design; B's skill θ_B is manufacturing/implementation. Berliant & Fujita (2008 IER, 2011 SEJ) already analyze heterogeneous complementary knowledge and knowledge creation; this itself is not novel.

## Four entrant architecture
Two independent municipalities G_A,G_B and two independent LOCAL private hub intermediaries P_A,P_B all simultaneously choose binary entry in OWN jurisdiction; public free, private endogenous fee. G+P co-location in same region permitted, total 16 profiles. Existing old-goods n_A=5,n_B=2 Cournot industries operate and generate profits and CS with ZERO hubs. Venture founders are separate active companies; each A/B matched founder pair can create differentiated new goods in separate independent niches, earning its OWN project profits (no multiple counting). Municipal governments maximize LOCAL welfare; independent P operators maximize own fee receipts−fixed cost.

Backward order 5 differentiated-product Cournot ->4 matched JV R&D ->3 founders' participation ->2 private posted fee game ->1 4 simultaneous entry flags.

## Fifth/fourth stages and new quality kernel
Use Stage32 exact 2-market differentiated Cournot for old vs new good:
consumer U_i=a_i Q_i+a_i^v y_i−[b_iQ_i²+2χ_iQ_i y_i+d_i y_i²]/2.
If new good succeeds, incumbent quantity x_i=[2d_i A_i−χ_iD_i]/[2b_id_i(n_i+1)−n_iχ_i²], new good y_i=[b_i(n_i+1)D_i−n_iχ_i A_i]/same denominator. Old-product baseline x_i0=A_i/[b_i(n_i+1)]>0. New goods JV operating profit Π=Σd_i y_i² and CS + incumbent profit change M_i.

Founder types θ_A,θ_B uniform[0,1]. Complementarity in new-product chance p(e,θ_A,θ_B)=p0+ζe+λ θ_A θ_B where type A is research skill and B is implementation skill, λ>=0 with p(e*,1,1)<1. JV effort cost ke²/2, optimally e*=ζΠ/k (true backward induction; effort is constant across types under this separable success function). Define base p_b=p0+ζ²Π/k, private per-pair net R0=p_bΠ−k(e*)²/2, R1=λΠ. A/B founder private payoff shares ω_i(R0+R1θ_Aθ_B). Per-match LOCAL social effect S_i(θ_A,θ_B)=S_i0+S_i1θ_Aθ_B, S_i0=p_b M_i+ω_i R0, S_i1=λ(M_i+ω_iΠ), possibly S_B0<0,S_B1<0 even though each founder private profit positive.

## Endogenous participant ability distributions
MASS m_i = ∫join dθ_i, quality first moment q_i=∫join θ_i dθ_i. Facility direct net value u_i=v_i+ρw_i h_j. P special value Δ_iθ_i, P charges fee f_i, G is free.
 a_i=ρμ N_j ω_i R0, b_i=ρμN_jω_i R1.
 Founder utilities U_0=0, U_G=u_i+a_i m_j+b_i θ_i q_j, U_P=U_G+Δ_iθ_i−f_i.
Zero facility m_i=q_i=0. G-only or Both m_i=1,q_i=1/2. P-only founder cutoff ℓ_i satisfies:
 m_i=1−ℓ_i; q_i=(1−ℓ_i²)/2; f_i=u_i+a_i m_j+(Δ_i+b_i q_j)ℓ_i.
NB use b_i^{match} to distinguish from demand slope b_i.

## Stage2 prices
Other market G or Both implies m_j=1,q_j=1/2 fixed; other absent implies m_j=q_j=0. P only global monopoly price f_i*=[(u_i+a_i m_j)+(Δ_i+b_i q_j)]/2; cutoff ℓ_i=[(Δ_i+b_i q_j)−(u_i+a_i m_j)]/[2(Δ_i+b_i q_j)] in assumed 0<U<T region.
With G+P co-located, private extra advantage Δ_iθ_i means f_i=Δ_i/2, half local users to P, total full skill mix m_i=1,q_i=1/2.
When BOTH regions P only, all-membership and fee interdependence:
 T_i=Δ_i+b_i(1−ℓ_j²)/2, C_i=a_i+b_i ℓ_iℓ_j,
 f_i=u_i+a_i(1−ℓ_j)+T_iℓ_i
 f_i=(1−ℓ_i)(T_i−C_i C_j/T_j).
Four equations for ℓ_A,ℓ_B,f_A,f_B. This accounts for PRICE Nash while other posted price is held fixed, including feedback in membership. Unique PARTICIPATION fixed point for fixed posted fees if max_i (a_i+b_i)/Δ_i<1. GENERAL P/P PRICE equilibrium global uniqueness NOT formally proved, but is verified via multiple roots and direct price deviations in reproducible numerical instance.

## Social welfare and first-stage entry
L0=ρ μ N_A N_B m_A m_B, L1=ρ μ N_A N_B q_A q_B. Local real hub direct service value Z_i=0 if 0; =N_i u_i if G; =N_i(u_i m_i+Δ_i q_i) if P; =N_i(u_i+3Δ_i/8) if Both. Local W_i=W_i0+Z_i−FG_i g_i−FP_i p_i+L0 S_i0+L1 S_i1. Private operator π_i=N_i f_i×(m_i if P-only;1/2 if Both)−FP_i.
For each i, define D_G0=W_i(G,s_j)−W_i(0,s_j), D_G1=W_i(Both,s_j)−W_i(P,s_j); D_P0=π_i(P,s_j)−π_i(0,s_j), D_P1=π_i(Both,s_j)−π_i(G,s_j). All 16 pure entry SPNE completely characterized by BOTH local 2x2 sign systems: state0 iff DG0<=0,DP0<=0; G iff DG0>=0,DP1<=0; P iff DG1<=0,DP0>=0; Both iff DG1>=0,DP1>=0. SOCIALLY OPTIMAL second-best entry = highest W_A+W_B over all 16 downstream equilibrium configurations.

## PROVED PROPOSITIONS

### P1 Sorting with complementarity (different sign from uniform average productivity)
If other locality G or Both, m_j=1,q_j=.5, P-only destination i:
 T_i(λ)=Δ_i+[ρμN_jω_iΠ/2]λ; U_i=u_i+ρμN_jω_iR0.
 f_i*(λ)=(T_i+U_i)/2, ℓ_i=.5−U_i/(2T_i), m_i=1−ℓ_i.
Therefore ∂f_i/∂λ>0, ∂ℓ_i/∂λ>0, ∂m_i/∂λ<0. Yet when success advantage is *homogenized* to identical λ/4 for every pair (same expected success conditional on all firms attending), T_i=Δ_i, U_i=u_i+ρμN_jω_i(R0+λΠ/4), so ∂m_i/∂λ>0. Same population N, R&D technology at full participation and competitive market, but sorting creates different paid-platform price AND participation responses.

### P2 Additional free public venue recruits lower ability founders
Given existing local P-only cutoff ℓ_i, adding G leads full participation. New enrolled mass ℓ_i, quality moment ℓ_i²/2; average newcomer ability ℓ_i/2 compared with existing premium-platform average (1+ℓ_i)/2. Thus quality-weighted new matches rise with ℓ_i² but aggregate count rises ℓ_i. This is NOT equivalent to renaming volume matching parameter ρ, which would preserve a constant quality/volume ratio. Proof uniform truncated moments.

### P3 Closed conditional PUBLIC entry thresholds
Other region j has G or Both => m_j=1,q_j=.5. Adding local G to existing P-only changes local real amenity welfare by H_i^u=N_i[u_i ℓ_i−Δ_i/8+Δ_iℓ_i²/2].
Define C=ρμN_iN_j.
Municipal welfare difference D_Gi^1=H_i^u−FG_i +C(S_i0ℓ_i+S_i1ℓ_i²/4).
Constrained planner's welfare difference D_SO=H_i^u−FG_i +C[(S_i0+S_j0)ℓ_i+(S_i1+S_j1)ℓ_i²/4].
Govt opens local free G iff fixed cost below threshold F_loc*=H_i^u+C(S_i0ℓ_i+S_i1ℓ_i²/4); second-best planner prefers addition iff FG_i<=F_SO* with added foreign terms. **These are exact conditional comparisons**; simultaneous first-stage equilibrium separately requires 4 independent BRs for all profiles.

### P4 Counterexample in full 4-entrant game: owner allocation depends on skill sorting despite identical full-participation innovation mean
Urban: n_A=5,N_A=3,v=.15,w=.22,FG=.09,FP=.65. Rural n_B=2,N_B=1,v=.1,w=.45,FG=.0405,FP=.24. Both product markets a=2,c=1,b=d=1,new-product net demand intercept=1, substitution χ=.9. ρ=.4,μ=.25,ω_A=.8,ω_B=.2,p0=.05,ζ=.005,k=.0004,λ=.5. Maximum matched new-good success probability p_b+λ=.55691630279<1. All costs/firms constant in comparisons.

Under **skill complementary project quality** p=p_b+λθ_Aθ_B, fully-solved B private-only tariff with A Both: f_B=.6410074280605, m_B=.639945174957, q_B=.435180261481. Marginal B GOV coentry welfare +.00003167010245, national +.00079932699126. Across 16 candidates numeric unique pure NE=(Both,Both); global 16-profile planner=(Both,Both).
Under **homogenized productivity** p=p_b+λ/4, which has the EXACT SAME expected success probability and aggregate welfare when both towns fully participate, B private fee=.6405924498931, m_B=.6405924498931, q_B=.435413106463. Marginal B GOV coentry welfare −.00045677076684, national +.00093681610672. Numeric unique pure NE=(Both,P-only), global planner=(Both,Both).
Both models' fully subscribed (Both,Both) welfare SAME, 2.407438368575148. Local B conditional government entry critical fixed cost: skill .04053167010245 vs pooled .04004322923316, whereas actual FG_B=.0405 lies between. At λ=0 Stage32 count-only comparison unique NE=(Both,P-only), second-best planner=(Both,Both); under complementary λ=.5 owner outcome flips to public-private coexistence in rural site. Root critical λ≈.35118403184.
This is a **new-to-Stage32 theorem/counterexample**, NOT proof unabsorbed by existing papers. Numeric PP pricing subgame has global deviations searched/checked; no global all-parameter general proof of PP price equilibrium.

## What remains for journal
Crucial need to distinguish type θ for technical production ability from premium-venue taste; equality/positive monotonic correlation is a structural assumption. Mechanism may disappear if uncorrelated. All nonhub firms' background production remains, but product launch effort is identical across types due separability; pair-specific commercialization output quality or endogenous profit splitting/investor location would deepen. Need policy alternative private subsidy/fee cap/contracting vs taxpayer G before saying *government direct ownership* is necessary. Berliant/Fujita knowledge complementarity and conventional selection in platforms/public service must be absorbed/given sharp non-absorption comparison.

Full source/manuscript and exact computational artifacts in chat:
STAGE_33_SECTORAL_COMPLEMENTARITY_BACKWARD_INDUCTION.md
stage33_skill_complementarity.py
stage33_16_complementarity_profiles.csv

References: https://www5.cao.go.jp/j-j/cr/cr25-1/chr25-1_01-02.html ; https://www.oecd.org/en/publications/enhancing-rural-innovation-in-japan_8d3b682a-en/full-report/understanding-innovation-in-rural-japan_f09da41a.html ; https://doi.org/10.1111/j.1468-2354.2008.00512.x ; https://www.rieti.go.jp/en/publications/summary/10040016.html.
