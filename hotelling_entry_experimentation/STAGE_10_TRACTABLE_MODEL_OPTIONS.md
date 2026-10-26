# Stage 10 — Tractable theoretical specifications after direct prior-art review
Date 2026-10-27. Status: **MODEL DESIGN OPTIONS, NOT FREEZE OR FULL SPNE**.

## What the primary literature actually simplifies
1. Matsushima & Shinohara (2014, JEBO) analyzes 1 supplier, 2 possible buyers, 1-vs-2 relationships, and quadratic sunk quality cost f*q^2 or d*f*q^2 with d in [1,2]; bilateral Nash bargaining. https://doi.org/10.1016/j.jebo.2014.07.014
2. Armstrong (2006, RAND) uses two fixed-endpoint platforms, two uniform mass-one sides, full coverage, single homing and linear cross-side externalities; under weak network feedback obtains linear demand systems and closed-form platform prices. See Section 4, equations (7)-(12), and discussion of the drawbacks of full coverage. https://doi.org/10.1111/j.1756-2171.2006.tb00037.x
3. Kim (2024, Manchester School) uses fixed-demand full-coverage mixed duopoly two-sided platform Hotelling competition. https://doi.org/10.1111/manc.12467
4. Cremer, Marchand, Thisse (1991, IJIO) exploits quadratic Hotelling preferences for spatial mixed oligopoly SPNE. https://doi.org/10.1016/0167-7187(91)90004-5
5. Caillaud & Jullien (2003, RAND, Chicken & Egg: Competing Matchmakers) shows both-sided matching intermediaries and potential coordination/multiple equilibria. Need verify full version provenance. https://cepr.org/publications/dp2885
6. Petrongolo & Pissarides (2001, JEL) surveys tractable matching-function restrictions and empirical issues. https://doi.org/10.1257/jel.39.2.390

## Essential variables to keep endogenous
- g in {0,1}: local government direct establishment.
- private entry choice, fixed cost FP.
- realized matching mass and probabilities, determined by how BOTH sides sort among existing operators. For the minimal closed-form model, full participation of each side is imposed but destination choice is endogenous.
- location and investment need not BOTH be fully endogenous in the first theoretical proposition. The key issue is the physical-facility causal margin, not maximum choice-set cardinality.
- Ex-post bargaining and relationship-specific investment can remain with exactly one possible link in a base equilibrium, while the scope-1-vs-2 margin belongs to a separate theorem, not to an unexamined claim that base scope one suffices for the final paper.

## Proposed analytical baseline ("A"): match production bilinear, full coverage
M innovators and R partners, independently uniform types x,y in [0,1]; each belongs to one available hub G or P (single-homing) while at least one exists. All are available ex ante but actual link realization remains endogenous.
Per-hub MATCHED LINK MASS:
  L_j=mu*N_j^I*N_j^R, mu>0, mu*M<=1, mu*R<=1.
These restrictions imply feasibility L_j <= min(N_j^I,N_j^R).
An innovator link probability is q_I,j=mu*N_j^R and partner's is q_R,j=mu*N_j^I.
The restriction is a stated technological reduced form, NOT the same Poisson matching model as Stage 8.

One match (m=1) ex post bilateral Nash:
  effort e*=beta*A/c;
  innovator net payoff u=beta^2*A^2/(2*c);
  partner transfer tau=beta*(1-beta)*A^2/c;
  social match surplus S=u+tau=beta*(1-beta/2)*A^2/c.
Set proposal cost kappa=0 in the BASE model (at q>0 all proposers attempt one link); knowledge spillover beyond partner and innovator project income xi=gamma=0 at this stage.
Innovator joining hub j receives
  U_I,j(x)=v_I-p_I,j-t_I*(x-l_j)^2 +u*mu*N_j^R.
Partner joining hub j receives
  U_R,j(y)=v_R-p_R,j-t_R*(y-l_j)^2+tau*mu*N_j^I.
The entire mass M and R is on one side or other (full coverage); no-participation is retained only if both hubs absent. This loses the participation-extensive margin, which matters in thin markets: restore in the paper's next lemma or robustness.
Two active hubs have fixed symmetric positions a=1/4,b=3/4 (distance D=1/2). A lone hub may center its own position 1/2. The roles and positions are a deliberately restricted physical service-design technology. Two-sided posted joining fees p_j^I,p_j^R>=0. Zero variable matching administration cost is a parametric specialization; each hub instead incurs fixed cost FG or FP.

## Explicit linear two-sided membership system
Let dI=N_G^I/M-1/2, dR=N_G^R/R-1/2; DI=p_P^I-p_G^I, DR=p_P^R-p_G^R.
Define a_I=t_I*D, a_R=t_R*D, w_I=u*mu*R, w_R=tau*mu*M.
Interior both-side demand satisfies
  a_I*dI-w_I*dR=DI/2,
  a_R*dR-w_R*dI=DR/2.
Set determinant Q=a_I*a_R-w_I*w_R.
If Q>0, exactly:
  dI=(a_R*DI+w_I*DR)/(2Q),
  dR=(w_R*DI+a_I*DR)/(2Q).
This is a CLOSED-FORM endogenous two-sided matching and membership demand; after those shares, matched volume L_G+L_P=mu*M*R*(1/2+2*dI*dR).
Zero-participation/tipping corners are outside this interior branch.

## Conditional mixed-duopoly fee equilibrium and profit
Assume t_I t_R D^2>S^2*mu^2*M*R (strong weak-feedback condition), both computed fees nonnegative, and large enough reservation values to sustain full coverage. Then public welfare is strictly concave in shares and uniquely maximized at dI=dR=0 (because positions symmetric); public can implement this by price matching p_G^I=p_P^I and p_G^R=p_P^R.
Private FOCs for unilateral posted fee deviations against FIXED G prices imply:
  p_P^I=p_G^I=t_I*D-mu*R*tau,
  p_P^R=p_G^R=t_R*D-mu*M*u.
This solves the smooth interior price game. Global deviations/tipping and boundary price regimes must still be certified before an unconditional SPNE existence claim.
Market shares N_j^I=M/2,N_j^R=R/2, so L_j=mu*M*R/4 and combined matching=mu*M*R/2.
Private pre-entry revenue/gross profit:
  Pi_P^op=(M*(t_I*D-mu*R*tau)+R*(t_R*D-mu*M*u))/2.
With D=1/2 and T=M*t_I+R*t_R:
  Pi_P^op=T/4-mu*M*R*S/2.
Under the no-marginal-cost assumptions private enters after public setup iff FP<=Pi_P^op, conditional on the selected branch.

## Complete realized-market-structure welfare for conditional branch
Let B=M*v_I+R*v_R, T=M*t_I+R*t_R.
No hub: W_empty=0.
One hub at 1/2, public only: W_G=B-T/12+S*mu*M*R-FG.
One hub at 1/2, private only: W_P=B-T/12+S*mu*M*R-FP.
Two hubs at 1/4,3/4, equal shares: W_GP=B-T/48+S*mu*M*R/2-FG-FP.
Incremental total welfare from a PRIVATE entry conditional on existing public:
  W_GP-W_G = T/16 - S*mu*M*R/2 - FP.
The private's PRE-FIXED-COST entry revenue equals T/4-S*mu*M*R/2. Hence
  Pi_P^op - [W_GP-W_G+FP] = 3*T/16 > 0.
This is an explicit PRIVATE vs SOCIAL entry threshold wedge. It appears to reproduce a known business-stealing mechanism; NOT a novelty claim.
Compared with a PRIVATE single-hub outcome (when P enters even without G), adding G while P remains yields
  W_GP-W_P=T/16-S*mu*M*R/2-FG.
The tradeoff is between diversity reduction in mismatch cost (T/16) and fragmented matching loss (S*mu*M*R/2), plus public setup cost FG.
No-public private entry gross operating revenue in a full-coverage centered-monopoly candidate:
  Pi_P^mono=M*(v_I-t_I/4+u*mu*R)+R*(v_R-t_R/4+tau*mu*M).
The one-hub monopoly must be checked against profitable partial coverage across BOTH sides, and prices/participation off-path must be verified. With sufficiently high v_I,v_R the formula can be supported, but not claimed globally without proof.
To finish Stage 0: no-G P enters iff FP<=Pi_P^mono, with-G P enters iff FP<=Pi_P^op; municipality compares the ACTUALLY induced outcomes, not hypothetical mixed-duopoly welfare if P stays out.

## Sanity-check parameter witness (interior candidate only)
Take M=R=t_I=t_R=1, mu=1/2, beta=1/2, A=c=1. Then
 u=1/8, tau=1/4, S=3/8.
 D=1/2 implies price pair p^I=3/8,p^R=7/16.
 Private gross operating profit Pi=13/32.
 Gross social gain from private entry relative to public-only T/16-S*mu*M*R/2=1/32.
 Private market entry socially excessive for 1/32<FP<=13/32, given the candidate mixed continuation.
 All formulas can be independently checked as exact rational arithmetic.

## Scope-1 vs scope-2 extension preserving bargaining economics
Replace power gross project technology with A*[1+eta*(m-1)]*e for m=1,2, 0<eta<1, partner-specific disagreement loss as in Stage 8, and investment cost c*d_m e^2/2. Then
 r_1=beta*A,
 r_2=A*[1+(2beta-1)*eta],
 e_1=beta*A/c,
 e_2=A*[1+(2beta-1)*eta]/(c*d_2).
At beta=3/4, eta=1/2, d_2=2, A=c=1:
 r_1=3/4,r_2=5/4;
 e_1=3/4, e_2=5/8 (lower);
 u_1=9/32, u_2=25/64;
 total social project S_1=15/32, S_2=35/64 (higher).
This exemplifies the intensive-effort vs extensive-partner tradeoff WITHOUT presupposing a sign.
If two requests succeed independently with probability q and each has a distinct partner, k=2 vs k=1 expected innovator gain is
 F(2,q)-F(1,q)=q*(18-11*q)/64-kappa.
For q in [0,1] this is quadratic and gives a radical threshold in q, but endogenous q=mu*N_j^R then feeds back into participation. Full two-sided fee equilibrium with both scope choices does NOT mechanically inherit the linear baseline solution; it must be solved piecewise, and general closed forms might not exist.

## Design / novelty / identification gates
1. Verify all price-branch global deviations, no-G private monopoly and complete stage-0 SPNE before stating a final theorem.
2. Never confuse endogenous matching QUANTITY with endogenous extensive-margin PARTICIPATION; total masses fixed under full coverage.
3. Purely bilinear matching is an assumed tractable technology, subject to supply capacities; do not mislabel as empirically universal.
4. The extra public facility in the current conditional model can have value due to diversification, despite lost matches, but whether a PRIVATE additional facility or a subsidy can achieve identical results is unaddressed. Therefore the paper does NOT yet establish why public ownership or a physical hub is necessary.
5. Benchmark matches Armstrong's linear two-sided-demand structure; novelty will require a proposition not reducible to platform demand plus generic free-entry business stealing.
6. Restore either n=1 vs 2 as endogenous matching/hold-up source, or a verifiable public-only feature not contractible to a private alternative, before promotion to midtier-field-paper main model.

## Research recommendation
Choose Model A as the CLOSED-FORM SOLVABLE CORE (public build / private entry / endogenous realized bilateral matching and endogenous sorting preserved). Treat partner-count choice n=1 vs 2 as the NEXT LEMMA and test whether it actually delivers a non-absorbed proposition. If partner-count-dependent hold-up is essential to the paper's novelty, replace scope-one baseline with the minimal binary-scope main model, not a final one-link-only paper. Do not silently freeze.
