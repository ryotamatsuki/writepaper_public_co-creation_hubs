# Stage 3 — Candidate Mechanism Table

Date: 2026-09-01

Base: `daf7b9e196f5557e08378ea2e5f347acb4b6d589`

Workflow: `07466bcb1a6d3bc654b52945f21b034b38e45281`

## Scoring rule

Weights were fixed before ranking: novel strategic mechanism 25; closest-paper survival 20; mechanism clarity 15; minimal-model tractability 15; welfare content 10; institutional relevance 10; empirical bridge 5. Each dimension was scored 0–5 and converted to 100. Scores are structured research judgments, not statistical estimates. Fatal prior art or an assumption-driven result overrides the numerical score.

| ID | Mechanism | Strategic feedback loop | Endogenous margins | Minimum new primitive | Closest prior-art threat | Expected nontrivial result | Welfare wedge | Tractability | Referee risk | Score | Verdict |
|---|---|---|---|---|---|---|---|---|---|---:|---|
| M1 | Local-funder service-composition externality | local funding → service mix → outside benefits → local funding | local vs external effort | geographically bounded objective | Cremer et al. (1997) | local over-weighting / external under-weighting | omitted external surplus | HIGH | collapses to standard spillover underprovision | 57.5 | **KILL** |
| M2 | Shared intermediary / strategic regional funding | regional contributions → common intermediary task mix → regional payoffs → contributions | contributions and task allocation | common agent serving several principals | Siqueira & Sandler (2004); common-agency literature | distorted geographic/task allocation | free riding + competing principals | MEDIUM | already a common-agency/multiple-public-goods problem | 51.7 | **KILL** |
| M3 | Service-mix-induced participant sorting | service mix → participant composition → match productivity → service mix | intermediary mix and participation | endogenous participation by outside/local actors | Rysman (2009); two-sided-platform theory; Carayol & Sterzi (2021) | interior composition distortion / participation tipping | platform fails to internalize cross-side participation | HIGH | generic two-sided-market feedback unless innovation-specific primitive matters | 73.1 | **SECONDARY / TOP 3** |
| M4 | Public KPI / multitask intermediation | measured output weights → task effort → measured performance → contract weights | contract and task effort | differential measurability | standard multitask agency; Russo et al. (2019) | effort distortion | principal-agent measurement wedge | HIGH | standard multitask moral hazard with new labels | 58.9 | **KILL** |
| M5 | Public-private specialization and crowd-out | public entry → private specialization/exit → residual demand → public task mix | entry, specialization, service mix | private intermediary alternative | Carayol & Sterzi (2021); generic public/private crowd-out | non-monotone composition / entry regimes | appropriability vs externality | MEDIUM | generic crowd-out/channel-substitution model | 59.1 | **KILL** |
| M6 | Direct vs intermediated search | intermediary offer → private direct search → matching return → intermediary offer | direct/intermediated mode, search | outside option of direct search | Carayol & Sterzi (2021); Calcagnini et al. (2019) | selection across transfer modes | search/selection externality | HIGH | formal TTO literature already contains the core choice | 51.0 | **KILL** |
| M7 | Scarce boundary-spanning capacity with cross-side response | allocation of scarce brokers → outside participation → returns to allocation → broker allocation | capacity allocation and participation | fixed rival intermediary capacity | platform/resource-allocation theory; Rysman (2009) | interior service mix / corner switch | intermediary ignores participation spillovers | HIGH | capacity constraint alone is mechanical | 69.1 | **SECONDARY** |
| M8 | Disclosure-mediated external search | disclosure policy → outside participation/search → matching quality → disclosure incentives | disclosure intensity, participation | credible anonymization/disclosure service | Kawakami (2024); Boudreau & Lakhani (2015); Howells & Thomas (2022) | disclosure threshold / participation regime | information provision changes matching market thickness | MEDIUM | disclosure + matching + public service is already formalized | 73.2 | **SECONDARY / TOP 3** |
| M9 | Endogenous overembeddedness / network state | brokerage → network topology → future learning/competition → brokerage incentives | links, openness, brokerage | endogenous network state | König (2016); Dasaratha (2023) | sparse/dense network regimes / non-monotonicity | network externalities | LOW | close strategic-network prior art and too complex for minimal Stage 4 | 55.0 | **KILL** |
| M10 | **Funder-induced neutrality loss / regional home-bias feedback** | regional local-priority incentive → intermediary home bias → outside participation/trust ↓ → external match pool ↓ → relative return to local matching ↑ → stronger local tilt | public sponsor's bias/priority choice; outside participation; realized match composition | intermediary neutrality affects outsiders' expected matching payoff | de Cornière & Taylor (2019); Zennyo (2022); Reimers & Waldfogel (2023/2026) | composition wedge; participation threshold/tipping; possible non-monotone local welfare | geographically bounded sponsor does not value outsider surplus and fails to internalize loss of market thickness/neutrality | HIGH–MEDIUM | may collapse to standard biased intermediary/two-sided platform | **83.8** | **KEEP / PREFERRED** |
| M11 | Asymmetric regions / intermediary location | location → participation/returns → funding → location | location and regional funding | regional asymmetry | standard location/fiscal competition + Cremer-style regional games | asymmetric equilibrium / relocation threshold | jurisdictional competition | MEDIUM | asymmetry alone is not a mechanism | 51.8 | **KILL** |
| M12 | Service-specific capability formation | current service mix → future search capability → future returns → current service mix | service mix and capability | service-specific capability accumulation | Kamien & Zang (2000); Leahy & Neary (2007); Bae (2016) | dynamic composition distortion | agents under/overinvest in future capability | LOW–MEDIUM | generic endogenous absorptive-capacity model with dynamics | 61.7 | **SECONDARY** |

## First-pass kill result

Seven candidates are killed: **M1, M2, M4, M5, M6, M9, M11**.

Five survive at least the first pass: **M3, M7, M8, M10, M12**. M7 and M12 are retained only as secondary mechanisms because the incremental novelty is weak relative to standard platform/capability theory.

## TOP 3

1. **M10 — Funder-induced neutrality loss / regional home-bias feedback** — preferred.
2. **M8 — Disclosure-mediated external search with endogenous participation** — strong institutional fit but severely threatened by Kawakami (2024) and disclosure theory.
3. **M3 — Service-mix-induced participant sorting / market-thickness feedback** — tractable but severely threatened by canonical two-sided-market theory.

## Interpretation of the ranking

M10 ranks first because it gives a single compact feedback loop linking a geographically bounded public principal, intermediary neutrality, endogenous participation, and realized service composition. The preferred mechanism is **not certified novel**. Its high score only justifies one Stage 4 minimal-model kill test. If stripping the regional/innovation labels reduces the mechanism to de Cornière–Taylor or Zennyo-style biased intermediation/platform participation, Stage 4 must return `NO-GO`.