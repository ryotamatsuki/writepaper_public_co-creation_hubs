# Stage 22 — Urban/rural asymmetric TWO REGION endogenous ownership compositions PP/GP/GG
Date 2026-10-27. STATUS: research architecture / equilibrium DEFINITIONS exact; no computed PP/GP/GG equilibria and no new theorem.

## User goal
Rather than presuppose one G and one P (Stage19) or private incumbents (Stage17), let the two-region co-creation hub network's number, urban/rural location and ownership composition emerge in a SIMULTANEOUS entry game. Compare no hub, monopoly, PP, GP or PG, GG, and possible regionally overlapping or 3-/4-hub configurations. Urban U has larger potential innovators/partners than rural R: M_U>M_R>0; N_U>N_R>0, with asymmetric valuations, costs, travel disutility and cross-regional contact. When public funding controlled by one prefecture, both regional public facilities are selected by ONE public decision-maker, NOT two independent public players. A GG structure means same government operates two sites; its sites maximize JOINT social welfare rather than competing as independent public entities. If two independent municipal governments, their objectives are different and add a separate institutional extension.

## Critical simultaneous entry game — DO NOT confuse market structure 3x3 with strategic action space
Potential entrants: private operator P_U in urban and P_R in rural, each chooses p_i∈{0,1}; ONE government chooses g=(g_U,g_R)∈{(0,0),(1,0),(0,1),(1,1)}. These 3 independent players act simultaneously. Stage-1 has 4×2×2=16 possible entry profiles, and there are 4 facilities maximum. PP, PG, GP, GG are merely some 2-provider market structures. Do NOT simply postulate "each region independently selects ∅/P/G" as one 3-action player — that would replace independent simultaneous municipal/firm entry with an exogenous proprietor-selection rule. To reduce to at most one entrant per region, must explicitly justify capacity/rationing/infeasible joint operation; otherwise 3-/4-provider states and overlapping within-region occupancy are indispensable OFF-PATH entry deviations.

Stage 2 for initial tractability: location discretely tied to U or R (the two accessible catchments); no extra continuous Hotelling optimization. Stage 3: two-sided tariff or tractable access/quantity competition (can start with posted fees fixed for benchmark, then endogenize fees). Stage 4: endogenous participation (including outside option), local vs cross-region search; stage5: innovators make sunk effort BEFORE realized matching; stage6: ex-post Nash bargaining. Then W(g,p), Π_U(g,p), Π_R(g,p) computed AFTER complete lower equilibrium and net real fixed costs. The welfare-maximizing public actor chooses a complete vector g and thus must check all FOUR g deviations, rather than four local playerwise unilateral deviations.

The simultaneous first-stage Nash iff:
  g*∈argmax_{g∈{0,1}²} W(g,p*).
  p_i*=1 => Π_i(g*,1,p_-i*)>=F_Pi;
  p_i*=0 => Π_i(g*,1,p_-i*)<=F_Pi.
Define Π_i here as private gross equilibrium continuation profit WHEN its entrant joins (all other entries fixed). W includes FG and FP as REAL resource losses. If downstream multiple equilibria, W and Π need a declared continuation selection, or a robust equilibrium-set treatment; do not write as single-valued without specifying that.

Examples:
PP at (g=00,p=(1,1)): both net private best responses enter; and W(00;11)>=max W(10;11),W(01;11),W(11;11).
GP at urban-private/rural-public (g=01,p=(1,0)): P_U profitable, P_R not, AND W(01;10)>=max_{g'} W(g';10).
GG at (g=11,p=(0,0)): neither private firm wants to enter anywhere, AND W(11;00)>=max_{g'} W(g';00).
These ALL involve off-path 3-/4-site comparisons.

## Minimal matching-investment primitives, without ad hoc two-hub bonus
Urban/rural potential innovators M_i and partners R_i endogenous active masses I_i ≤M_i, N_i≤R_i.
Within-region probability coefficients μ_i; cross-region connection efficacy ρ∈[0,1] and encounter parameter μ_c; within project productivity A_ii, cross project productivity A_ij. For each innovator of region i, expected gross project yield per unit effort:
  T_i=μ_i A_ii N_i + ρ μ_c A_ij N_j, i≠j.
Assume total match probabilities μ_i N_i+ρ μ_c N_j ≤1 or explicitly capacity-normalize and competing contacts. With ex-post Nash investor share β and ex-ante cost c_i e_i²/2:
  e_i* = β T_i / c_i.
  aggregate REAL experimental value in region i = I_i (β−β²/2)T_i²/c_i.
Thus new links operate at BOTH extensive and intensive margins; real project heterogeneity appears through A_ij rather than granting surplus to "G presence" by fiat. Cross-benefit exists even for PP if private-private linkage is technically feasible; to make ownership matter, specify endogenously chosen interoperability/neutral data access ρ_{s_U,s_R} and real provision costs, not a gratuitous PG quality bonus.

Baseline market asymmetry:
  M_U>M_R>0; N_U>N_R>0 (may differ in proportional size); local participation decline with transport cost to remote site, proper outside option.
Prefer a two-region DISCRETE location model with region-specific (not necessarily monopolistic) user fee/quality markets. Asymmetric location is not itself a new theorem.

Social welfare compare:
W_{PP}=S_{PP}(M_U,M_R,ρ,...)−FP_U−FP_R;
W_{GP}=S_{GP}−FG_R−FP_U;
W_{PG}=S_{PG}−FG_U−FP_R;
W_{GG}=S_{GG}−FG_U−FG_R.
S_o includes consumer benefits, real travel costs, investment costs, knowledge spillovers, not payment transfers.
This comparison is a constrained PLANNER CHOICE, NOT the simultaneous entry equilibrium; it is meaningful only after all downstream equilibria are computed for each ownership pattern.

## Candidate journal-grade questions (NOT proved)
1. Public-rural/private-urban ownership specialization is endogenous under intermediate rural market size & cross-region technologically complementary matches, but PP or GG reemerges outside sharp regions. Need strict equilibrium AND social-optimal comparisons; avoid hard coding cost efficiencies.
2. Cross-hub connection can cause PRIVATE operators to benefit from an unprofitable-for-private rural public facility even if privately operated urban hub poaches rural participants — quantify threshold for entrant profit difference and cross-spillover social value. Sign changes may require opt-out.
3. Number-of-hubs correct (2) yet TYPE composition wrong: e.g. PP equilibrium but GP social optimum, or GG equilibrium but GP optimum. Novel if endogenized pre-match effort, public commitment to neutral connection and urban/rural market thickness, and vanishes in reasonable benchmark models.
4. Geography (M_U/M_R) and interoperability (ρ) interact NONMONOTONICALLY in decentralization incentives; neither more connections nor bigger urban market need unambiguously improve rural efficient provision. Check against Zhu et al 2021 and Peitz/Sato 2024.
5. Direct public building vs cost-effective contracted private link: MUST run implementation-equivalence no-go.

## Closest priors (as of 2026-10-27)
- Inoue, Kamijo & Tomaru (2009), "Interregional mixed duopoly", RSUE 39(2):233–242: equilibrium endogenous location can co-locate or separate, local public objective; DOI 10.1016/j.regsciurbeco.2008.10.001.
- "Differentiated mixed duopoly and decentralization in a two-city model" (2018), Papers in Regional Science 97(4):1425–1441; DOI 10.1111/pirs.12307, municipality vs regional government objective and privatization.
- Akio Kawasaki (2026), "Private Firms' Location Choice Under Multi-Market Mixed Duopoly With Cost Differentials", Bulletin of Economic Research 78(4):950–967 DOI 10.1111/boer.70066. Strongest RECENT competitor: high-/low-cost region private entry, public ownership intensity, distorted socially suboptimal siting; baseline public already in both markets, quantities not two-side matching, no fixed entry costs. Has a policy-subsidy comparison, so do not claim this is unaddressed generally.
- Martin Peitz & Susumu Sato (2024), "Asymmetric Platform Oligopoly", CEPR DP19584: logit-based tractable two-sided pricing with within and cross-side network effects, entry and compatibility plus optional nonparticipation; plausible tractable stage3 pricing alternative to unstable Stage14 location-price.
- Zhu, Li, Valavi & Iansiti (2021), "Network Interconnectivity and Entry into Platform Markets", Information Systems Research DOI 10.1287/isre.2021.1010; interconnectivity may HELP OR HURT entrant profitability, selective entry/market targeting.
- Ambrus & Argenziano (2009), "Asymmetric Networks in Two-Sided Markets", AEJ Micro 1(1):17–52 DOI 10.1257/mic.1.1.17; endogenous asymmetric networks, so simply producing asymmetric ownership/geographical membership not new.
- Dong, Bárcena-Ruiz & Garzón (2024), "Cooperation between governments to set up public firms", Economic Systems 48(2) 101189: asymmetric cross-government ownership and welfare.

## Recommended phased scope
Part I: one-prefecture three-player simultaneous decision (4 governmental vectors ×2×2 private), urban vs rural locations FIXED to regions, realistic M_U>M_R, R_U>R_R, **ex ante investment endogenous**, cross region connections ρ perhaps exogenously fixed; solve all 16 entry-continuation profiles including collocated providers, tariffs kept simple enough to certify global best responses. Compare PP vs GP vs GG equilibrium ownership distributions and planner.
Part II: endogenize interoperability ρ and public neutral-access commitment; then test novel type-vs-count distorted equilibrium and welfare.
Part III: allow private candidate geographic relocation (choice U/R), dynamic pricing, municipal-vs-prefectural jurisdiction and policy tools as robustness. Avoid re-adding both continuous Hotelling positions and 4-player strategic site games simultaneously.

NEVER claim PP/GP/GG equilibria already computed. Stage 19 is one G one P with exogenous permitted zones and cannot itself realize PP or GG. Stage22 is new architecture, not a solved theorem.
