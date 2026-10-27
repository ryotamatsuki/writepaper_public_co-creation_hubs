# Stage 24 — Structural audit and minimal endogenous extensions of Stage 23

Date 2026-10-27. PURPOSE: explain every primitive, verify which decisions are actually endogenous, identify next solvable microfoundation, and define non-absorption/referee tests. **NO new Stage 24 equilibrium claimed.** All Stage 23 phase diagrams CONDITIONAL upon fixed local split, fee rate and owner-specific access coefficients.

## 1. What Stage 23 truly endogenizes
- Three simultaneous players: ONE prefecture choosing (g_U,g_R)∈{0,1}² jointly, TWO region-specific independent profit maximizers choosing p_U,p_R∈{0,1}; 16 stage-1 profiles, up to four sites including collocation and simultaneous 3-/4-site off-path deviations.
- Conditional on profile, bilateral match success chance for innovator in site h, region i:
  q_ih=μ [R_ih+ρ γ_h Σ_{h'∈other region} γ_h' R_jh'], with q_ih<=1.
- Every innovator chooses PRE-MATCH noncontractible effort e_ih=β A q_ih/c; ex-post Nash split β/(1−β) of A e; private h earns fixed share σ of partner's settlement. W net of entry resource fixed costs; transfers cancel.
- W, profits and ALL pure simultaneous entry Nash profiles, constrained planner over all 16; multiple pure Nash including PP/PG/GP/GG. Mixed entry Nash not resolved.

## 2. What Stage 23 DOES NOT endogenize
- Region sizes M_U=R_U=1, M_R=R_R=r, their 1:1 innovator/partner mix (no side-specific regional imbalance).
- Supply/location: each private only builds in its preassigned region; government regional flags; within-region continuous service locations absent. Transport cost t/(12n_i²) corresponds to EXOGENOUS optimally spaced n_i sites, not a location Nash.
- User enrollment: all residents enroll if any local site; NONE if their region has zero local hubs; if G/P both operate locally, each gets exactly half of each user group. No consumer outside-choice, prices, mobility, asymmetric brand/taste sorting. This 50:50 split is imposed, NOT an endogenous Nash outcome.
- Prices: no posted admission fees; private broker receives exogenously fixed ad-valorem charge σ on partner's Nash share, without bargaining/commission feedback into participation and β. No private fee/commission best response.
- Cross connectivity ρ exogenous common physical/communication opportunity; government γ_G=1, private γ_P=α<=1 EXOGENOUS owner advantage. No governance or access investment costs, access-rights bargaining or multi-homing.
- Productivity A uniform within/cross; cost c uniform; valuation v uniform across towns; government spending has real cost equal its face value; one prefecture internalizes both regions, not local-government rivalry.
- Link capacity and probability: q<=1 restricted calibration; bilateral μ I R matching treats multiple matches as EXPECTED mass and needs capacity checks if regions grow.
- Subsidy/procurement/digital broker alternatives; partnership types/quality heterogeneity; match formation after investment success single-shot (no endogenous rejection or match persistence).

## 3. Derivations
Let (I_ih,R_ih) be fixed site user masses and γ_h ownership openness, Rhat_j=Σγ_h'R_jh'. q_ih=μ[R_ih+ρ γ_h Rhat_j]. Innovator solves
  max_e β A q e - c e²/2 -> e*=β A q/c.
Gross social real experimentation per innovator:
 q A e* -c(e*)²/2 = (β−β²/2) A² q²/c.
Set κ=(β−β²/2) A²μ²/c.
Gross welfare
  B(g,p)= Σ_i 1{n_i>0}(M_i+R_i)[v−t/(12n_i²)]
           +κ Σ_i,h I_ih [R_ih+ρ γ_h Rhat_j]².
W=B−Σ(F_Gi g_i+F_Pi p_i).
Private P_i (gross) operational commission formula as in Stage23:
  σ(1−β) A μ R_iP [I_iP e_iP* + ρ α Σ_{h'∈j}γ_h' I_jh' e_jh'*].
This fee is a TRANSFER; if bargainers take fee into account, the fixed β formula may need reformulation.

For exactly ONE hub in each region and symmetric within-region cohorts 1,r, effective connectivity z=ρ γ_U γ_R:
 e_U=β A μ(1+zr)/c,
 e_R=β A μ(r+z)/c,
 S(r,z)=κ[1+r³+r(1+r)(2z+z²)].
Within + cross effects expand as:
 baseline κ(1+r³);
 NEW pair access first-order 2κ r(1+r)z;
 larger PRE-MATCH effort contribution κ r(1+r)z².
Cross term not necessarily all "externality": some private platform participants appropriate benefits. Separate private spillovers and ownership fees in welfare analysis.

With stage23 α=.5,ρ=.6,r=.7, β=2/3,A=c=1,μ=.4:
 z_PP=.15, z_PG=.3, z_GG=.6;
 e_U PP=.2946667 PG=.3226667 GG=.3786667;
 e_R PP=.2266667 PG=.2666667 GG=.3466667.
The public/private composition differences follow the EXOGENOUS γ multipliers. At α=1 and public fixed costs above private, public advantage disappears in same-site-configuration comparisons.

## 4. Nontrivial strategic lesson in Stage23
At r=.7,ρ=.6 ALL four PP/PG/GP/GG are simultaneous entry Nash, while PG is unique constrained planner optimum. The planner can REPLACE the rural private with public (PP→PG), but a unilateral government deviation under PP can only ADD a public rural site while private rural P's entry stays fixed; G CANNOT force P_R to exit. Likewise GP/GG can be locked-in by jointly required changes, creating an inefficient ownership-composition equilibrium even though every observed profile has exactly two sites. This is a coordination/entry externality, not itself new without deeper primitives.

## 5. What to endogenize FIRST (highest value, reasonable analytical effort)
### A: Endogenous private/public openness and interoperability effort x_h
Replace exogenous γ_G=1,γ_P=α with x_h∈[0,1] chosen simultaneously AFTER entry at real convex cost C_h(x)=k_h x²/2. Cross success z=ρ x_U x_R. Public decides x for owned sites to maximize W; private decides x for profit; for GG single public actor JOINTLY optimizes both x_U,x_R. For local collocated sites, define all pairwise cross links x_h x_h'.
For exactly one site per region:
 S(r,z)=κ[1+r³+r(1+r)(2z+z²)].
 Given private P_U private revenue
  Π_U^op(z)=φ[1+z r(1+r)+z²r],
  φ=σβ(1−β)A²μ²/c.
Hence
  ∂S/∂x_U=2κ r(1+r)(1+z)ρ x_R,
  ∂Π_U/∂x_U=φ r[(1+r)+2z]ρ x_R.
For the parameter domain β∈(0,1),σ∈(0,1],r>0, the social marginal return > private intermediary marginal return for identical marginal resource cost; hence private may underinvest in interoperability. Both sides have complementarities in x and equilibrium multiplicity possible. For strict global concavity in OWN x require k_G > 2κ r(1+r)ρ² (public) and k_P>2φ rρ² (urban P), etc. However high k may imply low x; prove existence and do not promise owner superiority.
Remove γ privilege before optimizing. Only if optimal access levels differ ENDOGENOUSLY and the public direct-production is efficient vs subsidized private provision do we obtain a rationale for government owner.
### B: Endogenous participation and sorting incl remote market
Each side's utility when enrolling hub h:
  u_s(z,h)=v_{is}(z)-\ell_h^s-t_s distance(z,site)^2 + expected private Nash payoffs conditioned on hub memberships; outside utility 0.
Users can enroll own/remote/none, firms pay fees endogenous later. Choose logit or uniform reservation values to obtain manageable share equilibrium; with high network strength fixed points may multiply, so explicitly characterize consumer participation and off-path NE selection. Big payoff: remote market may be completely absent without a rural hub; cross-local access may newly include users; now "two sites create MORE opportunities than one" is possible.
### C: Endogenous broker fees
Let private operator choose \ell_P^I,\ell_P^R (or σ itself) jointly with other platforms and fee constraints; test rival-platform price BR including corners. Cannot take uniform local split 50:50 if fees are chosen. For initial tractability, separate Step C from endogenous location: discrete region anchor. Peitz/Sato "Asymmetric Platform Oligopoly" is now in RAND Journal of Economics vol.57(1), pp.60-77 (issue 2026; first published Aug 5 2025), doi:10.1111/1756-2171.70008. The 2024 CEPR DP was the PREPRINT; cite journal version in new manuscripts.
### D: Bilateral match quality / cluster specialization
Let A_{ij} vary by technological complements and travel costs; do not hard-code “public match inherently higher output.” To establish increased TOTAL matches from dual hubs, remote region agents must opt out without local hub or sites offer nonredundant capabilities. Cross complementarity A_UR>A_UU must be tied to skill distributions/assortativity, not named alone.
### E: Location choices
After A/B/C have tractable equilibrium, allow P_U/P_R location choice U or R (or continuous intra-region d), expand firm actions to {out,U,R}, stage1 profiles 4×3×3=36 with two potential private entrants; government initially vector still 4; beware multi-site overlapping 3/4 sites and stage14 nonexistence if unrestricted continuous Hotelling.
### F: Ownership/contract implementation test
If exactly same efficient x levels, quality, impartiality and access can be contracted from private provider for resource cost not exceeding govt direct supply, publicly owned facilities aren't strictly optimal. Need endogenous enforceability/incomplete contracting or verifiable public testbed/data required; benchmark Hart/Shleifer/Vishny (1997 QJE). Do not infer physical hub requirement merely from market failure in connectivity.

## 6. Research priority ordering / suggested nesting
M0 current exact 16-profile restricted game (benchmark).
M1 = M0 + endogenous x (cross-access effort), keep masses, fixed commissions, fixed within-region micro-catchments. Goal: derive PP/PG/GG access-investment FOCs and global incentives, private underprovision and planner thresholds, no exogenous owner advantage.
M2 = M1 + endogenous rural/urban user ENROLLMENT with outside option (could first keep fee fixed), test two-site new participation and match growth; subgames & multiplicity.
M3 = M2 + optimal two-sided platform fees + matching mass reoptimization, choose a solvable discrete-choice demand technology (e.g. Peitz/Sato), not full Hotelling/free location concurrently.
M4 = M3 + private siting U/R or heterogeneous sector-specific capabilities and explicit welfare policy alternatives (access subsidy, public procurement, public data grants, digital-only hub).
If M1 doesn't have a new nonabsorbed ownership threshold then do NOT inflate M2/M3 as a remedy; prior theory (e.g. Zhu et al 2021, Mangin & Julien 2021) already contains network-entry and endogenous surplus logic.

## 7. New publication-target hypotheses, NOT results
H1 Same access technology and costs can yield endogenously PUBLIC higher cross-access investment in a GP network than private in PP, but a sufficiently efficient private consortium can dominate PG, with an analytically sharp boundary. Need see whether subsidizing x_P replicates before interpreting public ownership.
H2 In asymmetric regional sizes, urban-private/rural-public PG is optimal only on an INTERMEDIATE market-thickness and connectivity range, even where decentralized PP/GG persist as ownership-inefficient NE; derive actual equilibrium under endogenous x and user enrolment. Stage23's PP→PG→GG transition is NOT novel unless robust to replacing exogenous γ hierarchy.
H3 Interoperability and participation interact to cause BOTH market-creation and cannibalization, making private entry externality switch sign at a tractable cutoff. Compare prices, costs, matching/innovation effects and fees.
H4 Ownership inefficiency can occur with socially correct NUMBER of sites, and be Pareto-ambiguous; public access commitments or access subsidies reduce it. Must compare Mankiw-Whinston, Peitz-Sato, Zhu et al, mixed-duopoly ownership papers.

## 8. Specific methodological / accounting checks
- Govt is single prefecture optimizing JOINT g_U,g_R; if municipalities independently own sites, this is another game.
- Fees are TRANSFERS and commission after bargain can change effective bargained division unless contract timeline specified.
- Potential overmatching when masses/μ increase; enforce q<=1 and partner capacity or saturating technology/matching CRS.
- Stage23 transportation savings from n_i increase are ASSUMED ideal spacing and do not certify strategic spatial location.
- At most 2 hubs per area currently, no public/private composition effects on local user mass except 1/2 split.
- If \rho=0, public premium connection yields no project benefits; government can still enter solely if lower fixed cost or other primitives, not from the connection.
- At r→0 fixed positive incremental rural build costs exceed any cross-site match gains; “more rural scarcity implies more public building” is false in the baseline.
- Stage 23 input numbers are calibrated ILLUSTRATIONS, not observations from Ehime or any region.

Verified external literature (as of Oct 27 2026):
Kim (2024) DOI 10.1111/manc.12467.
Peitz & Sato 2026 RAND DOI 10.1111/1756-2171.70008 (first online Aug 2025).
Zhu et al (2021) ISR DOI 10.1287/isre.2021.1010.
Mangin & Julien (2021) JET DOI 10.1016/j.jet.2021.105208.
Atakan, Richter, Tsur (2026) Theoretical Economics DOI 10.3982/TE6786.
Hart, Shleifer, Vishny (1997) QJE DOI 10.1162/003355300555448.
OECD (2026) Incubation in Entrepreneurial Ecosystems DOI 10.1787/21be9d03-en.
