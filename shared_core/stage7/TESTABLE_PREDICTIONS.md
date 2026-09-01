# Stage 7 — Testable Predictions and Observability

These are theory implications only. No empirical design is authorized at Stage 7.

## Comparative-statics audit

The strategic threshold is

\[
A^*=A_L-\tau-4\beta(1-\beta)-\Omega,
\]

where

\[
\Omega=\sqrt{2\beta\left[2\beta(1-2\beta)^2+(1-\beta)\tau^2\right]}.
\]

Exact derivatives:

\[
\boxed{\frac{\partial A^*}{\partial A_L}=1.}
\]

\[
\boxed{\frac{\partial A^*}{\partial\tau}=-1-\frac{2\beta(1-\beta)\tau}{\Omega}<0.}
\]

\[
\frac{\partial A^*}{\partial\beta}
=-(1-2\beta)\left[4+\frac{4\beta-16\beta^2+\tau^2}{\Omega}\right].
\]

The \(\beta\) derivative is globally negative under the transparent sufficient condition \(0<\beta\le1/4\). Stage 7 does not claim a global sign for the wider \(0<\beta<1/2\) domain.

## P1 — Neighbor-entry interaction changes with core attractiveness

Prediction:

\[
\Gamma(A_M)=D_1-D_0
\]

moves from negative to positive through the canonical crossing region. Thus the marginal establishment incentive is lower when the neighboring entrant exists below \(A^*\), but higher when it exists above \(A^*\).

### Observability

- Dependent variable: establishment/opening probability or investment commitment of peripheral option \(i\).
- Key explanatory variables: neighbor establishment \(E_j\), common-core attractiveness \(A_M\), and their interaction.
- Ideal unit: peripheral market/region × decision episode.
- Candidate environments: local marketplaces, retail centers, regional airport/gateway investments; co-creation hubs only if sufficient historical decisions exist.
- Identification feasibility: difficult but conceivable with exogenous core-quality/connectivity shocks.
- Confounders: correlated regional demand, common subsidies, macro shocks, strategic policy coordination.

## P2 — Anti-coordination becomes coordination

Prediction:

Below \(A^*\), establishment outcomes should disproportionately take the one-entry form when fixed costs are intermediate. Above \(A^*\), intermediate-cost environments should be more likely to display clustered symmetric outcomes—either both peripheral options or neither.

### Observability

- Dependent variable: joint outcome \((E_1,E_2)\).
- Explanatory variables: core attractiveness and proxies for fixed establishment cost.
- Ideal unit: paired peripheral markets sharing the same core × decision episode.
- Candidate data: mall openings around a downtown core, airport/gateway infrastructure, local platform launches.
- Identification feasibility: moderate for long panels; weak for contemporary public co-creation hubs because of few observations.
- Confounders: joint planning, common funding programmes, regional demand correlation.

## P3 — Core-side feedback, not entrant-side feedback, is decisive

Prediction:

The strategic reversal should be strongest in environments where loss of core users materially reduces core value. It should disappear or be substantially weaker when the core's value is insensitive to its installed base, even if peripheral options themselves exhibit local network effects.

### Observability

- Outcome: estimated interaction of neighbor entry with core attractiveness.
- Explanatory heterogeneity: measures of core network sensitivity, e.g. connectivity/frequency, user-side matching thickness, footfall/variety dependence.
- Ideal unit: market × core architecture.
- Identification feasibility: conceptually strong but data intensive.
- Confounders: endogenous compatibility, pricing, core investment response.

## P4 — Cross-regional friction lowers the required core-attractiveness threshold

Because

\[
\frac{\partial A^*}{\partial\tau}<0,
\]

higher friction between peripheral options makes the substitute-to-complement transition occur at a lower level of core attractiveness.

Economic interpretation: when cross-use of the neighboring peripheral option is more difficult, the direct business-stealing force between peripheral entrants is weaker, so less core strength is required for residual-core feedback to dominate.

### Observability

- Outcome: sign/strength of neighbor-entry interaction.
- Explanatory variable: geographic distance, travel time, eligibility differences, or other cross-peripheral friction.
- Ideal unit: pair of peripheral markets sharing one core.
- Identification feasibility: plausible using geographic variation, but distance may correlate with market segmentation.

## P5 — Social marginal entry interaction changes at the same threshold

Stage 7 derives

\[
\boxed{G_1-G_0=2\Gamma.}
\]

Therefore the same \(A^*\) that changes private strategic interaction also changes social marginal establishment benefits from decreasing to increasing.

### Observability

This is harder to estimate directly because it requires welfare measurement. Candidate proxies include user surplus, access/connectivity, realized matching value, and resource costs before and after the first and second entry.

- Ideal unit: paired-entry episodes with detailed user-level outcome data.
- Identification feasibility: low to moderate depending on environment.
- Major confounders: transfers, pricing, omitted consumer groups, induced investment.

## P6 — Coordination multiplicity has an asymmetric welfare ranking

In the canonical coordination region, whenever \((0,0)\) and \((1,1)\) are both Nash equilibria,

\[
SW_{00}>SW_{11}.
\]

Thus observed clustered full entry need not indicate socially beneficial coordination; it can be excessive duplication.

### Observability

- Outcome: aggregate user welfare/resource cost under zero versus full peripheral entry.
- Explanatory state: core attractiveness above the strategic threshold and intermediate fixed costs.
- Identification feasibility: difficult; best viewed as an institutional/policy prediction rather than an immediate reduced-form regression.

## Empirical observability summary

| Prediction | Observable in principle? | Causal identification | Best candidate environment |
|---|---|---|---|
| P1 interaction sign | Yes | Moderate/difficult | retail, airports, platform entry |
| P2 joint regime clustering | Yes | Moderate | paired infrastructure/market entry |
| P3 core-feedback dependence | Yes with architecture data | Difficult | digital platforms, airports, retail agglomerations |
| P4 friction shift | Yes | Moderate | spatial retail/airport markets |
| P5 social threshold coincidence | In principle | Difficult | detailed structural case studies |
| P6 excessive full-entry equilibrium | In principle | Difficult | public infrastructure / retail development |

At least P1–P4 are externally meaningful without adding model primitives.
