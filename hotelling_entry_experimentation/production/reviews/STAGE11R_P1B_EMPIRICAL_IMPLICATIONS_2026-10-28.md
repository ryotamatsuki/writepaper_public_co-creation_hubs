# Stage11R P1b — observable predictions, empirical identification, and policy interpretation

**Date:** 2026-10-28 JST. **Scope:** frozen STAGE33-2026-10-28-THEORY-FREEZE-v1, JICT-focused applied-theory manuscript. **Parent:** P1 exact entry-cost bands and 3D local open neighborhood in Draft PR #58. **Decision:** P1b EXPOSITION IMPLEMENTED, empirical/theoretical identities only. Stage11 overall remains HOLD for novelty, actual external referee review, and strong generality claims.

## The theory-to-observation bridge

### I. Fee selection (conditional prediction, strong testable sign)

Frozen service-only regime P, other region's participant moments fixed, service quality fixed, founders uniformly distributed and willingness to pay perfectly aligned with skill. Put U = u+a^m m_other, T = Delta+b^m q_other, 0 < ell = (f-U)/T < 1.

- Participation mass m=1-ell; dm/df=-1/T<0.
- Unnormalised skill moment q=(1-ell^2)/2; dq/df=-ell/T<0.
- Participating mean skill q/m=(1+ell)/2; d(q/m)/df=1/(2T)>0.
- ***NO claim that project completion count rises***. Jointly endogenous fees, opposite-region responses and non-aligned skill valuation invalidate direct unconditional causal extrapolation.

Observables: full eligible population denominator, ex-ante technical skill distribution, posted/offered net fees, each firm's take-up and project participation, provider identity, counterpart region. Identification: independent tariff reform or external fee shock, nonprice quality stability, propensity/sorting control. A fee vs participant quality correlation alone is not causal.

### II. Free basic service versus specialist fee revenues (conditional prediction)

When B co-provision is selected in a jurisdiction and uniform increment Delta*theta applies: premium fee f_B=Delta/2, paid premium take-up = 1/2 of eligible population, provider gross receipts N*Delta/4. Holding the other region's matching moments and technology fixed, the private-only optimum receipts exceed co-provision receipts under U>0 and T>=Delta. Positive co-provision profits when F_P<N*Delta/4 are compatible with ***lower revenue***, not proof of crowding-in.

Observables: cohort-specific free/basic access counts, paid-upgrade customers, posted and net fees, specialist revenue, actual fixed cost and exit. Interpretation: avoid an unqualified 'one-half upgrade share' prediction outside uniform tastes, quality-invariance, no capacity/rationing, and single premium-provider assumptions. Avoid comparing non-equivalent public and private providers.

### III. Geographic incidence of innovation (conditional accounting object, no universal sign prediction)

There is a prior-stage microeconomic check on the same channel: the success technology has partial cross-derivative d²p/(d theta_A d theta_B)=lambda>=0 at fixed research effort and interior probability. An exact cross-difference-in-differences in ex-ante skills is lambda*(theta_A^H-theta_A^L)*(theta_B^H-theta_B^L). Under positive lambda, higher technical skills are complements in conditional project success. Recovering this from observational data requires failed as well as successful matches, ex-ante skills, effort, project composition and exogenous or conditionally ignorable matching; a regression interaction is not an automatically causal parameter. The checker now also verifies this difference using exact frozen probabilities.

For j fully covered, the nonlocal increment from i's public addition is exactly
Delta W_j = C(S_j^0 ell_i + lambda K_j ell_i^2/4)
          = F_i^national - F_i^local,
with K_j=M_j+omega_j Pi. M_j is the downstream consumer-plus-incumbent external effect of successful products in market j, and Pi is combined venture return. Its sign depends on product-market displacement and incidence; the local/national difference is ***not*** automatically positive. The same result can be checked against the frozen W_j(B,B)-W_j(B,P) baseline by independent exact fractions.

Observables: geographic pairing of eligible applicants, ex ante skills, joint R&D investment and real outcomes, launch/new product market, incumbents' competitive prices/quantities/profits, consumer demand/welfare estimate, founder profit distribution, and government fixed costs. Identification requires a credible market-without-venture counterfactual. A regional firm-revenue or public-user-count regression alone cannot recover W_j or the planners' counterfactual.

### Numerical cost-band result (NOT a population-level prediction)

The certified three-scenario band T_skill0<F_GB<T_skill_half is 0.263% of the baseline cost. It proves nonisolated existence of a regional service-configuration contrast in the frozen uniform-skill Cournot model, not that empirical local governments will generally enter as lambda rises. Beta(2,2) and differentiated-Bertrand caution remain negative.

## Current nearest empirical references and scope

1. **Carpentier, D'Adda and Nepelski (2026)**, *Which firms receive digital innovation support? Evidence from the European Digital Innovation Hubs network*, **International Small Business Journal**, doi:10.1177/02662426261461491, published online August 29, 2026. Open journal abstract reports about 15,600 firms and systematic composition differences favoring larger/faster-growing/more digitally advanced EDIH customers; this is ***observational selection***, not causal evidence that posted prices screen technical ability or that public provision is insufficient. Supports measuring the eligible denominator, not simple registered-user shares.
2. **Eldebo and Hjelm (2025)**, *Acceleration of cleantech startups through open innovation matchmaking with established companies: The role of intermediaries*, **Innovation**, doi:10.1080/14479338.2025.2514462, published online June 10, 2025. Public abstract describes a dataset of 3,129 matchmaking meetings (2018–2022) involving 545 Swedish startups; illustrates possible process-level matching records, not an identification of new-good Cournot regional welfare. Full publisher text not fully accessible during the present check.
3. **Official JICT aims and scope**, https://link.springer.com/journal/10842/aims-and-scope (verified October 28, 2026): applied-theory submissions must clearly highlight implications for empirical analysis and economic policy. This is ***why*** the new manuscript section now explicitly states testable signs, record linkage and identification qualifications.

## Actual manuscript / reproducibility edits

- New section `production/manuscript/sections/06_empirical_implications.tex`, inserted between numerical and related literature in `manuscript/main.tex`, contains three formula-grounded predictions, measurement variables, identification issues and explicit policy boundaries.
- `manuscript/references.bib` adds two publication-supported empirical context references. Introduction roadmap updated, and link to Section `sec:numerical` verified.
- Independent `production/reviews/stage11r_p1b_empirical_implications_check.py` checks screening derivative algebra through exact rational examples, private-only vs co-provision revenue and survival, and equality between direct regional welfare-accounting differences and the cross-border wedge for all three frozen technology scenarios. It tests a **necessary manuscript consistency** condition, not actual data or identified causal effects.
- New CI workflow `.github/workflows/stage11r-p1b-empirical-implications.yml` reruns prior P/P and P1 analytic/rational certificates, the conditional identity checks and a complete LaTeX+BibTeX build, retaining the working paper and JSON audit trail.

## Gate and deliberate restrictions

| Item | Status |
|---|---|
| Three conditional empirical implications explicitly stated | DRAFT IMPLEMENTED |
| Sign derivation and model-to-observable mapping | ALGEBRAIC / exact rational audit, CI to be checked |
| Empirical data and causal inference | NOT PERFORMED |
| Generic effects of public hubs / universal cost-effective physical building | NOT MODELED |
| Whole-paper journal novelty | STILL CONDITIONAL, predecessor absorption not cleared |
| Economic robustness over skill distributions, alternative competition | NEGATIVE or NOT RE-CERTIFIED |
| External human mathematical/empirical review | NONE |
| Stages 13, 14, 15 (format, anonymous package, submission) | NOT CLOSED |

**Next bounded work:** Finish proposition-level absorption against nearest full-text predecessors and prepare a no-oversell Stage11 independent referee review. Separately, Stage13 manuscript integration and Stage14 JICT technical compliance require explicit later gates. No economic model alteration and no journal submission.
