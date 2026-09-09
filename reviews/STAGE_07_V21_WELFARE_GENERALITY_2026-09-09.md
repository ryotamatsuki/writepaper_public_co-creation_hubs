# Stage 7 v2.1 — Welfare, Generality & Institutional Validation

Date: 2026-09-09 JST

Canonical workflow: `ryotamatsuki/research-paper-workflow` v2.1 @ `fe6fa5f2a94632be578d5e9fab39e7df86013113`

Canonical template: `templates/STAGE_07_WELFARE_GENERALITY.md`

Upstream authorities:

- Stage 4A repaired-global certification: `GO — MATHEMATICAL ADVERSARIAL CERTIFICATION PASS`
- Stage 6 v2.1 novelty re-kill: `GO — GO TO WELFARE / GENERALITY`
- Surviving contribution only: matched-price public-public strategic sign reversal caused by activating the shared private intermediary's optimizing repricing response.

## 1. Executive welfare/generality verdict

**GO TO STAGE 7.5.**

The repaired mechanism has a coherent and nontrivial welfare interpretation, but the welfare and generality claims must remain narrower than the strategic contribution.

At the repaired full-game witness (`beta=.01`, `gamma=.825`, `tau=.35`):

- `W_1 = W_2 ≈ 0.2948786590`
- `Pi_T ≈ 0.0070310097`
- `W^N ≈ 0.5967883277`

At the matched-price B3 equilibrium using the same scalar private price `p_G ≈ .0184036796`:

- `W_1 = W_2 ≈ 0.2893342840`
- `Pi_T ≈ 0.0072598538`
- `W^N ≈ 0.5859284219`

Hence `W^N_G-W^N_B3 ≈ +0.0108599058` at this witness. This is **NUMERICAL WITNESS EVIDENCE ONLY**. It is not a theorem that endogenous repricing raises welfare.

More importantly, at the decentralized G equilibrium the local national-welfare derivative in the `+x_i` direction is

`dW^N/dx_i ≈ +0.4904468182 > 0`.

The government's own reduced derivative is numerically zero (`≈7.1e-8`), while the rival-region effect is `≈+0.5005263571` and the private-profit effect is `≈-0.0100796096`. Thus the repaired baseline exhibits **local under-provision relative to a coordinated marginal increase in one public investment**, after private fee transfers are accounted for correctly.

This local coordination result is separate from the sign of the public best-response slope. Strategic substitutability is not itself a welfare ranking.

## 2. Exact project surplus and welfare derivation

For region `r`, let `A_{rh}` denote the set of project types assigned to route `h` by the upper envelope of

`U_{rh}(z)=z q_h-kappa_{rh}-1{h=T}p_T`,

with non-participation utility zero. Exact project surplus is therefore

`CS_r = sum_h integral_{A_rh} [z q_h-kappa_rh-1{h=T}p_T] dz`.

This is the exact utility integral; no triangular-CS shortcut is used. In the smooth central ordering `0 -> H_T -> H_r`, this reduces to the manuscript expression

`integral_s^{t_r}(z q_T-kappa_T-p_T) dz + integral_{t_r}^1(z q_r-kappa_L) dz`.

The welfare section already relies on uniform support-side participation costs. Making that implicit microfoundation explicit, let the participation index for route `h` be `m_h` and individual support-side cost be `c~U[0,1]`. In an interior regime the participating mass is `n_h^P=m_h` and exact support-side surplus is

`integral_0^{n_h^P}(n_h^P-c) dc = (n_h^P)^2/2`.

The model's reduced-form participation equations identify the corresponding indices as

- `m_i = rho+x_i+beta n_i^F`,
- `m_T = rho_T+beta n_T^F`.

Under the baseline symmetric regional attribution, each government counts one half of national support-side surplus, so

`W_i = CS_i + (1/4) sum_h (n_h^P)^2 - (gamma/2)x_i^2`.

National welfare is

`W^N = W_1+W_2+Pi_T`.

Because project surplus includes `-p_T n_T^F` and private profit is `Pi_T=p_T n_T^F`, the fee cancels exactly at the national level. Run `34334195347` verifies

`p_G n_T^F = Pi_T = 0.00703100969181901`.

Thus the private access fee is a transfer in aggregate welfare under the baseline zero-real-private-cost normalization. Price changes matter for welfare through allocation and participation, not mechanically through the transfer itself.

### Exposition gap to repair later

The current Model section states the support-side mass equations but does not state the uniform support-side cost microfoundation as explicitly as the Welfare section uses it. This is not a new strategic assumption; it is the already-used welfare interpretation and should be made explicit during the next paper-build integration.

## 3. Planner-objective / choice-set register

No first-best or global planner optimum is claimed at Stage 7.

| Object | Objective | Choice set | Downstream pricing | Correct label |
|---|---|---|---|---|
| Decentralized G | each government maximizes its own `W_i` | own `x_i in [0,1]` | private intermediary reoptimizes `p_T` | decentralized full-game equilibrium |
| B3 | each government maximizes its own `W_i` | own `x_i in [0,1]` | `p_T` fixed at the scalar equilibrium G price | matched-price fixed-price identification benchmark |
| National welfare accounting | `W^N=W_1+W_2+Pi_T` | no optimization problem is asserted | follows evaluated regime | aggregate welfare accounting identity |
| Local coordination wedge | derivative of `W^N` in one public-investment direction at G | infinitesimal change in `x_i`, with downstream G continuation | private intermediary reoptimizes | local coordinated marginal welfare effect |

`B3` is not price regulation, not a planner problem, and not first best. The local coordination derivative is not a globally solved constrained optimum.

## 4. Benchmark-definition audit

**PASS.**

The defensible benchmark language is:

- `B3`: matched-price fixed-price benchmark / identification benchmark;
- `G`: endogenous-pricing full game;
- `W^N`: aggregate welfare measure;
- positive `dW^N/dx_i` at G: local coordination wedge / local under-provision in that direction.

Prohibited language unless a new optimization problem is separately solved and certified:

- `first best`;
- `social optimum`;
- `optimal coordinated investment`;
- `welfare dominance of G over B3` as a theorem.

## 5. Private versus social decision map

At the repaired G equilibrium:

- own reduced marginal welfare: `dW_i/dx_i ≈ 7.06e-8`, consistent with the decentralized FOC;
- rival-region welfare effect: `dW_j/dx_i ≈ +0.5005263571`;
- private-profit effect: `dPi_T/dx_i ≈ -0.0100796096`;
- aggregate marginal effect: `dW^N/dx_i ≈ +0.4904468182`.

Therefore the government does not internalize the full positive external effect of its investment under the baseline regional-welfare attribution. The result is **local under-provision relative to the coordinated `+x_i` direction**, not a global statement about the distance to an unrestricted social optimum.

The magnitude of this wedge is strongly influenced by the baseline assumption that each region receives one half of national support-side surplus. That attribution is a modeling normalization, not an empirically estimated share. The sign and size of the welfare wedge must not be advertised as general to alternative regional attribution rules.

## 6. Welfare propositions / thresholds

### W7-1 — Exact transfer cancellation

**PROVED by accounting identity.** Under zero real operating cost for `H_T`, `-p_T n_T^F + Pi_T = 0` in national welfare.

### W7-2 — Local coordination wedge at repaired G witness

**NUMERICALLY VERIFIED AT THE CERTIFIED BASELINE.** `dW^N/dx_i>0` at the repaired global-equilibrium witness. Hence a marginal coordinated increase in one public investment raises national welfare locally.

### W7-3 — B3 versus G welfare ranking at repaired witness

**NUMERICAL WITNESS ONLY.** `W^N_G>W^N_B3` by about `0.01086`. No global or parameter-uniform welfare ranking is established.

### W7-4 — Strategic sign and welfare sign are distinct

**CONCEPTUAL/ACCOUNTING RESULT.** `BR_i'` describes strategic response; it does not determine `dW^N/dx_i` or a welfare ranking between regimes.

## 7. Institutional evidence table

Primary-source recheck date: 2026-09-09.

| Primitive / institutional link | Current primary-source evidence | Classification | Stage-7 interpretation |
|---|---|---|---|
| Regional public startup/innovation intermediation | Ehime Prefecture, EGF Startup Community: public support, exchange, mentoring, support institutions, network building; FY2026 operation procurement confirms continuing prefectural program | ESTABLISHED | supports a regional public-intermediation interpretation |
| Paid metropolitan innovation/community option | CIC Tokyo official site: coworking and startup community linking innovators, investors and firms | ESTABLISHED | supports a priced metropolitan outside option |
| Commercial operator can alter charges | SHIBUYA QWS membership page and membership agreement: paid plans; operator may alter plans/charges for operational or economic reasons | ESTABLISHED for discretion; UNVERIFIED for model-specific objective | supports price flexibility, not `max p n_T^F` literally |
| Another paid innovation-community setting | Level39 official membership: paid London base with events, mentors and investors | ESTABLISHED | shows the commercial-community primitive is not Japan-specific |
| Municipal/university public incubation and partner connection | City of Helsinki Campus Incubators: jointly funded incubation, city connects incubators to relevant ecosystem partners | ESTABLISHED | second public-institutional setting consistent with facilitation/network role |
| Municipal startup incubation/network support | Barcelona Activa Startup Lab: public incubator designed to foster cooperation networks | ESTABLISHED | further supports public incubation/network facilitation |
| Exact private profit maximization `p_T n_T^F` | no institutional source establishes this exact one-period objective | UNVERIFIED / REDUCED FORM | retain explicitly as a modeling abstraction |
| Private price response `p^*_{T,x_i}<0` | Stage-7 repaired-model computation gives `≈ -0.01203693` | MODEL PREDICTION, not observed fact | empirical prediction only |
| Project primary-route single homing | plausible abstraction, not directly established by cited institutional pages | SUGGESTIVE | do not describe as literal exclusivity of firms |
| Partner multihoming | networks contain investors/experts/support institutions, but exact simultaneous multihoming rate is not documented | SUGGESTIVE | structural abstraction only |
| Two public jurisdictions sharing the same private alternative | plausible in metropolitan-facing innovation ecosystems but not directly identified causally by sources | SUGGESTIVE | motivating architecture, not established empirical fact |

Primary sources used:

- https://www.pref.ehime.jp/page/112623.html
- https://www.pref.ehime.jp/site/nyusatsu/135305.html
- https://jp.cic.com/cic-tokyo/
- https://shibuya-qws.com/membership
- https://shibuya-qws.com/membership/agreement/
- https://level39.co/workspaces/community-membership/
- https://www.hel.fi/en/business-and-work/campus-incubators-programme
- https://www.barcelonactiva.cat/en/-/incubadora-glories

## 8. Generality / robustness evidence-classification table

| Claim | Baseline form/class | Evidence type | Assumptions used | Maximum defensible scope now | Stage-7.5A attack target |
|---|---|---|---|---|---|
| `BR_B3'>0>BR_G'` at repaired witness | uniform projects; linear support-side feedback; quadratic public cost; repaired parameter vector | Stage 4/4A all-regime computational certification + local derivative checks | full baseline model; `beta=.01`, `gamma=.825`, `tau=.35` | BASELINE FUNCTIONAL FORM; existence at certified witness | attempt boundary/regime/quantifier overstatement and alternate wording |
| small-positive-beta reversal | same baseline functional form on regular interior branches | analytic local theorem | regular beta-zero branches, strict SOCs, nonsingular Jacobians, smooth continuation | SUFFICIENT-CONDITION THEOREM, local/branch-specific | verify theorem does not imply global equilibrium or a known explicit epsilon |
| remote rival-public route irrelevance at repaired witness | baseline project utility | analytic dominance inequality | `kappa_L+tau>v+alpha`; here `.62>.60` | SUFFICIENT CONDITION FOR ROUTE DOMINANCE ONLY | ensure manuscript does not generalize this to arbitrary tau |
| G aggregate welfare exceeds B3 | repaired baseline witness | numerical Stage-7 audit | baseline welfare attribution and same matched scalar price | NUMERICAL ROBUSTNESS ONLY at one witness | search nearby parameters for ranking reversal; prohibit theorem wording |
| local under-provision in G | repaired baseline witness | numerical derivative + exact welfare identity | equal regional attribution of national partner surplus; endogenous price continuation | NUMERICAL BASELINE RESULT ONLY | attack alternative attribution/parameter values; prohibit global optimum language |
| small `C^1` changes in project distribution preserve local sign | nonbaseline distribution | continuity argument presently stated in manuscript, not independently red-teamed | positive density, same ordering, nonsingularity, strict signs | CONJECTURED/CONDITIONAL LOCAL GENERALITY pending 7.5A | construct admissible perturbations and separate local sign from globality |
| nonlinear network functions preserve result | nonbaseline network functions | continuity intuition only at present | same derivative signs, smoothness, same regular branch | CONJECTURED GENERALITY pending 7.5A | construct nonlinear counterexamples / quantify sufficient restrictions |
| small regional asymmetry preserves result | asymmetric primitives/actions | continuity intuition only | same regular branch and strict signs | CONJECTURED LOCAL GENERALITY pending 7.5A | distinguish perturbing a state from existence of an asymmetric equilibrium |
| old `±0.5%`, 20-draw robustness around rejected vector | old canonical vector | stale computation | rejected global-equilibrium witness | REJECTED / NON-AUTHORITATIVE | must not survive freeze unless rerun around repaired witness |

## 9. Generality across two institutional settings

The mechanism can be interpreted without changing its strategic architecture in at least two distinct public-innovation settings, but this is institutional portability, not a general theorem.

### Setting A — prefectural startup-support ecosystem

A regional government funds/operates entrepreneurship facilitation and network formation (Ehime EGF-type activity), while projects can also use paid metropolitan innovation communities such as CIC Tokyo or QWS. The private option has price discretion and partner-network functions. This mapping supports the model primitives. It does not establish that actual Japanese prefectures satisfy the repaired parameter values or the predicted strategic sign reversal.

### Setting B — municipal/university incubation ecosystem

A city coordinates and funds campus incubation and connects programs with ecosystem partners (Helsinki-type activity; Barcelona Activa provides another municipal incubation example), while firms can also use paid private technology communities such as Level39/CIC. The same public facilitation / partner network / paid private alternative architecture is coherent without changing the model mechanism.

The two settings differ in governance structure and delivery organization, but both remain innovation-intermediation applications. Stage 7 does **not** claim generality to arbitrary two-sided markets, transport, payments, or other industries merely by relabeling.

## 10. Important repaired-witness interpretation: direct rival-hub use is inactive

The repaired global witness satisfies

`kappa_L+tau = .62 > v+alpha = .60`.

Therefore a nonresident rival public hub is dominated by non-participation for every project type and every history. This condition removes the finite-deviation/free-riding route that killed the old witness.

This has an important interpretation consequence: the certified global witness does **not** rely on direct project poaching between the two regional public hubs. The public-public strategic interaction operates through the shared private route and cross-side participation. Later manuscript language must not describe the repaired witness as proving direct inter-hub user competition.

The magnitude `tau=.35` is constructive, not empirically calibrated. Institutional sources make geographic/organizational frictions plausible but do not establish this dominance condition quantitatively.

## 11. Empirical predictions

The following are model predictions, not established empirical facts:

1. At the repaired configuration, stronger regional public facilitation reduces the shared private intermediary's optimal access price: `dp_T^*/dx_i ≈ -0.01204`.
2. Price rigidity versus price flexibility should mediate the sign of cross-jurisdiction public investment responses: the matched-price fixed-price environment has a positive local public BR slope, while endogenous repricing produces a negative one at the certified witness.
3. A credible empirical design should separately measure public facilitation, project routing, support-side participation, and private access fees; observing only public spending cannot identify the repricing channel.
4. Because the B3 construction matches the on-path private price level, empirical work should distinguish the level of commercial access prices from the responsiveness of those prices to public investment shocks.
5. The model predicts a positive national coordination wedge at the repaired baseline even though the public actions are strategic substitutes in G; strategic-substitute behavior is therefore not evidence of socially excessive investment.

## 12. Result-to-exposition triage

| Headline result | Economic object | Candidate vehicle | Why this vehicle | Verified source | Stage-10 action |
|---|---|---|---|---|---|
| matched-price public-public sign reversal | `BR_B3'>0>BR_G'` | theorem/proposition + compact comparison table | signs and benchmark distinction are the claim; a BR-curve plot could misleadingly imply global shape | Stage 4/4A + Stage 6 | replace rejected canonical witness with repaired witness; keep scope narrow |
| repricing decomposition | `M_G=M_fixed+P_price` | displayed equation + concise prose | mechanism is analytic and does not require a graphic | analytic derivation | retain as mechanism, not separate novelty |
| repaired global-equilibrium witness | full-strategy all-regime BR check | compact certificate/parameter table | communicates globality evidence without implying a global theorem | Stage 4A | add repaired parameters, global-BR residual/gain, slope signs |
| transfer cancellation | fee as transfer in `W^N` | displayed identity | exact accounting point | Stage 7 CI | retain prominently in welfare section |
| local coordination wedge | `dW^N/dx_i>0` | proposition-style numerical result + decomposition table | separates rival welfare and private-profit channels | Stage 7 CI | replace stale old-witness numbers |
| B3-G welfare levels | regime-level welfare | one compact table | useful benchmark comparison, explicitly numerical | Stage 7 CI | regenerate `welfare_comparison.tex` from repaired witness |
| institutional validation | public facilitation / priced private alternative | concise prose or appendix table | avoids overstating motivating cases as evidence | current primary sources | update access dates and evidence labels |
| generality limits | local theorem vs computational witness vs conjecture | concise scope table / prose | prevents theorem inflation | Stage 7 classification | remove stale 20-draw and unqualified `C^1` robustness claims unless re-certified |

No figure is required at this stage.

## 13. Policy scope and limits

Defensible policy interpretation:

- decentralized public investment may leave a positive cross-regional coordination externality even when the public actions are strategic substitutes;
- public-investment analysis should account for endogenous commercial responses, not only other governments' direct choices;
- the relevant empirical object is partly the responsiveness of private access pricing.

Not defensible from the present model:

- a general recommendation to subsidize public hubs;
- a quantitative optimal subsidy;
- regulation of private platform prices;
- a claim that price flexibility is socially beneficial in general;
- a claim that strategic substitutes imply overinvestment or harmful competition;
- direct claims about actual Ehime/Kagawa/Tokyo causal responses without empirical identification.

## 14. Candidate counterexample targets for Stage 7.5A

1. Quantifier attack on the small-beta theorem: verify local/branch-specific wording never becomes global/SPNE wording.
2. Globality attack on any continuity statement: local strict-sign persistence must not be used to claim persistence of global equilibrium.
3. Welfare-ranking attack: search nearby admissible parameters for `W^N_G-W^N_B3` sign reversal.
4. Regional-attribution attack: identify exactly where the equal split of national partner surplus enters the positive coordination wedge; prohibit general welfare claims beyond that baseline.
5. Remote-friction attack: verify all repaired-witness prose acknowledges `kappa_L+tau>v+alpha` and inactive nonresident rival-public use.
6. Asymmetry attack: distinguish continuity of derivatives at perturbed states from existence/globality of asymmetric equilibria.
7. Functional-form attack: challenge the manuscript's current `C^1` distribution/network wording with admissible perturbations or downgrade it.
8. Benchmark-language attack: ensure B3 is never called first best, price regulation, or a planner solution.
9. Evidence-maturity attack: remove or quarantine the old `.05/.9/.05` canonical vector, old exact-equilibrium wording, and old 20-draw robustness from any claim of global equilibrium.

## 15. Remaining fatal / major concerns

### Fatal concerns

None at Stage 7.

### Major concerns to carry forward

1. **Current manuscript is stale.** It still reports the rejected old canonical vector, old welfare table, and old exact-certificate equilibrium language. Stage 7 GO applies to the repaired theory record, not to the current submission manuscript.
2. **Welfare attribution is baseline-specific.** The equal regional attribution of national support-side surplus materially drives the coordination wedge and is not empirically calibrated.
3. **Strong remote friction limits the global witness interpretation.** Direct rival-public project use is inactive at the repaired witness.
4. **Welfare ranking is numerical only.** `G>B3` in aggregate welfare is not a theorem.
5. **Generality statements in the current robustness section are too mature relative to current certification.** The old 20-draw exercise is stale; `C^1` distribution/network and asymmetry language must remain conditional until Stage 7.5A.

## 16. CI evidence

GitHub Actions run `34334195347` on head `529594048ee792b8ca619314d3e7109efa63a781`: **SUCCESS**.

Terminal markers:

- `TRANSFER_CANCELLATION: PASS`
- `LOCAL_COORDINATION_CLASSIFICATION: UNDERPROVISION IN THE COORDINATED +x_i DIRECTION`
- `PLANNER_BENCHMARK_REGISTER: NONE CLAIMED; NO FIRST-BEST LABEL USED`
- `STAGE7_WELFARE_NUMERICAL_AUDIT: PASS`

## 17. Canonical verdict and Stage 7.5 contract

**GO TO STAGE 7.5.**

Stage 7.5 must make a full-paper value/freeze decision without adding extensions. It must evaluate the paper using only:

- the narrow Stage-6 novelty claim;
- the repaired Stage-4/4A global-equilibrium witness;
- the exact transfer-cancellation identity;
- the baseline-only local positive coordination wedge;
- the institutional evidence classifications above;
- the strict generality limits above.

Stage 7.5 must not restore the rejected old vector, broad follower-response novelty, generic complements-to-substitutes novelty, direct rival-hub competition at the repaired witness, or unverified welfare/generality claims.
