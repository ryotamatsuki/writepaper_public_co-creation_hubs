# Stage 7.5A v2.1 — Generality / Quantifier Red-Team Gate

Date: 2026-09-10 JST

Canonical workflow: `ryotamatsuki/research-paper-workflow` v2.1 @ `fe6fa5f2a94632be578d5e9fab39e7df86013113`

Canonical template: `templates/STAGE_075A_GENERALITY_QUANTIFIER_RED_TEAM.md`

Checklist applied: `checklists/THEOREM_CERTIFICATION_CHECKLIST.md`

Upstream authorities:

- Stage 4A repaired-global certification: `GO — MATHEMATICAL ADVERSARIAL CERTIFICATION PASS`
- Stage 6 novelty re-kill: `GO`
- Stage 7 welfare/generality: `GO TO STAGE 7.5`
- Stage 7.5 full-paper value decision: `GO TO STAGE 7.5A GENERALITY / QUANTIFIER RED-TEAM`

## 1. Executive scope-certification verdict

**GO — GENERALITY / QUANTIFIER CERTIFICATION PASS.**

The paper can proceed to Stage 8 only under the narrowed claim set recorded here.

The red-team found one material scope defect in the pre-audit manuscript: the old draft vector (`beta=.05`, `gamma=.9`, `tau=.05`) still appeared in the abstract, introduction, main-results section, robustness section, appendix, and generated tables as though its exact rational stationary-root certificate established an equilibrium pair. Stage 4A has already shown that this vector admits profitable finite public deviations once omitted routing regimes are restored. Therefore the exact Krawczyk result cannot serve as global-equilibrium or SPNE authority.

This defect is **scope/evidence misclassification, not a new mathematical failure of the surviving mechanism**. It has been corrected inside Stage 7.5A without changing the model or deriving a new theorem:

1. the analytic theorem is explicitly restricted to **regular interior stationary branches** and local best-response geometry;
2. the old exact rational certificate is reclassified as an **archived local stationary-root diagnostic only**;
3. the repaired vector (`beta=.01`, `gamma=.825`, `tau=.35`) is the sole active numerical equilibrium/welfare witness;
4. the repaired result is described as an **all-regime computational global-equilibrium witness under the baseline functional form and documented tolerance**, not a general theorem;
5. stale `±0.5% / 20-draw` robustness around the rejected vector is removed from the evidentiary chain;
6. distribution, nonlinear-network, and asymmetric-region robustness are downgraded to **unproved extension directions**;
7. the local welfare wedge is labelled a **baseline-specific local under-provision result**, not first best, social optimum, or global welfare theorem;
8. the active generated-result pipeline is rerouted from the rejected central-regime solver to the repaired Stage-4A evaluator.

No substantive rollback to Stage 4 or Stage 7 is required.

## 2. Formal quantifier table

| ID | Claim | Formal scope / quantifiers | Domain / class | Current evidence | Status |
|---|---|---|---|---|---|
| Q1 | First-order B3 complementarity | For any regular interior B3 **stationary branch** through `beta=0` satisfying the stated gap, smoothness and strict-local-SOC conditions, `M_B3(0)=0` and `dM_B3/dbeta|0>0`; therefore there exists a local interval of positive `beta` on that continuing branch with positive local BR slope | Baseline functional form; central interior regime; local branch | Analytic derivation | PASS |
| Q2 | Full-game beta-zero negative cross effect | At a symmetric regular beta-zero full-game **stationary state** with `T>0`, `Delta>0` and stated regularity, `M_G(0)<0` by the displayed closed form | Baseline functional form; local stationary state | Analytic derivation | PASS |
| Q3 | Local strategic sign reversal | **There exists** `epsilon>0` such that for every `beta in (0,epsilon)` **on the assumed continuing regular stationary branches**, `BR_B3'>0>BR_G'` | Baseline functional form; local branch; strict SOCs and nonsingular Jacobians | Analytic sufficient-condition theorem | PASS after wording correction |
| Q4 | Repaired baseline global witness | **There exists one documented baseline vector** (`beta=.01`, `gamma=.825`, `tau=.35`, remaining primitives fixed) at which independent all-regime computation finds global public best responses in G/B3 within the certification tolerance and opposite local slopes | Baseline functional form; public strategy set `[0,1]`; all project routes; documented numerical protocol | Stage 4A independent computational certification | PASS |
| Q5 | Remote-public route dominance at repaired witness | If `kappa_L+tau>v+alpha`, then for all project types and all histories a nonresident rival public route is dominated by non-participation | Baseline project utility; sufficient inequality only | Analytic utility bound | PASS |
| Q6 | Welfare transfer cancellation | For all evaluated allocations under zero real private operating cost, `-p_T n_T^F + Pi_T = 0` in national welfare | Baseline accounting identity | Exact identity + Stage 7 numerical check | PASS |
| Q7 | Local national coordination wedge | At the repaired G witness, `dW^N/dx_i ≈ .49045>0`; hence a marginal coordinated `+x_i` change raises national welfare locally | One repaired baseline witness; equal regional attribution of partner surplus | Stage 7 numerical derivative audit | PASS as baseline numerical result |
| Q8 | G versus B3 welfare ranking | At the repaired witness only, `W_N^G≈.59679 > W_N^B3≈.58593` | One repaired baseline witness | Numerical comparison | PASS only as one-witness result |
| Q9 | Distribution / nonlinear-network / heterogeneous-region robustness | No affirmative theorem claimed | Outside baseline functional form | Not independently certified | PASS because removed/downgraded |
| Q10 | Old exact rational certificate | There is a unique coupled **stationary root in its isolating box** with exact local derivative signs at the old vector | Rejected old vector; local algebraic box only | Exact rational Krawczyk / interval derivative diagnostic | PASS only as archived diagnostic; REJECTED as equilibrium authority |

## 3. Assumption-dependence table

| Result | Economic assumptions actually used | Shape / parameter restrictions | Normalization / tractability devices | What the assumptions do **not** establish |
|---|---|---|---|---|
| Q1/Q3 local B3 sign | public investment affects support-side attractiveness; cross-side feedback `beta`; shared private alternative held at matched price | regular central ordering; positive quality gaps; `d>0`; strict public local SOC; smooth/nonsingular stationary system | uniform project types; linear partner feedback; quadratic public cost | global public best reply; uniqueness; arbitrary distributions/functions |
| Q2/Q3 G sign | profit-max private follower reprices after public investment | symmetric regular beta-zero stationary state; `T>0`, `Delta>0`; private/public strict local SOCs; nonsingularity | zero private real operating cost is irrelevant for strategic sign but used in welfare | global sign over parameter space; necessary condition; genericity |
| Q4 repaired witness | same baseline game plus all-route deviations and endogenous private continuation | repaired vector; fail-closed numerical protocol; strategy set `[0,1]` | selected finite-decimal vector is constructive | positive-measure parameter region; exact global proof for all primitives |
| Q5 route dominance | remote public route carries additional access friction `tau` | `kappa_L+tau>v+alpha` | upper bound uses max partner mass `<=1` | necessity of the inequality; empirical calibration of `tau` |
| Q7 welfare wedge | government internalizes own regional welfare but not all rival-region effect; private profit included nationally | repaired witness; equal half attribution of national partner surplus | zero private real operating cost makes fee a transfer | first-best distance; globally optimal subsidy; robustness to other attribution rules |
| Q8 welfare ranking | same welfare accounting in both environments | repaired witness; matched scalar private price in B3 | baseline functions | general welfare dominance |

## 4. Function-class counterexample / attack audit

### 4.1 Uniform project distribution

The pre-audit manuscript asserted that the uniform distribution was “not locally essential” and that sufficiently small `C^1` perturbations preserve the result. This is too strong for the certified evidence chain if “the result” includes global equilibrium status. Smooth local derivative continuation may be plausible under density/ordering/nonsingularity restrictions, but Stage 4A did not certify a global best-response separation theorem over a `C^1` function class.

**Disposition:** affirmative robustness claim removed. Maximum defensible wording: this is an open extension direction.

### 4.2 Nonlinear partner/network response

“Smooth nonlinear perturbations preserve the result” is likewise not certified. Merely preserving the signs of first derivatives does not by itself establish all required cross derivatives, curvature, downstream uniqueness, route ordering, and global public optimality.

**Disposition:** no robustness theorem claimed. A future theorem would need explicit derivative/shape restrictions and an equilibrium-separation argument.

### 4.3 Regional asymmetry

Perturbing the action state is not equivalent to proving existence of an asymmetric equilibrium under perturbed regional primitives. The earlier prose blurred this distinction.

**Disposition:** removed as an affirmative robustness claim. No arbitrary or small-heterogeneity equilibrium theorem is claimed.

### 4.4 Public-cost curvature / parameter perturbations

The old 20-draw exercise is invalid as evidence for the repaired paper because it is centered on a vector that fails global optimality. A new cloud around the repaired vector is not required for the current paper because the analytic theorem and the repaired global witness have separate roles.

**Disposition:** old exercise permanently non-authoritative; no replacement robustness claim.

### 4.5 Low-friction direct rival-hub access

Stage 4A supplies an actual counterexample to careless generality: at the old low-`tau` vector a government profitably reduces own investment and free-rides on the rival public hub. This is preserved as a regression artifact. It demonstrates that local FOC/SOC and stationary-root isolation do not imply global equilibrium.

**Disposition:** repaired witness uses the transparent sufficient dominance condition `kappa_L+tau>v+alpha`. The paper does not call that condition necessary for sign reversal.

## 5. Baseline / robustness / general-theorem classification

| Object | Classification | Maximum current scope |
|---|---|---|
| beta-zero closed forms and local reversal | `SUFFICIENT-CONDITION THEOREM` | local regular stationary branches under baseline functional form |
| repaired `BR_B3'>0>BR_G'` witness | `BASELINE FUNCTIONAL FORM` + `COMPUTATIONAL GLOBAL-EQUILIBRIUM EXISTENCE WITNESS` | one repaired primitive vector under documented all-regime protocol |
| remote-public dominance inequality | `SUFFICIENT-CONDITION THEOREM` | route dominance only |
| transfer cancellation | `GENERAL ACCOUNTING IDENTITY WITHIN BASELINE WELFARE DEFINITION` | all evaluated allocations under zero real private operating cost |
| local under-provision | `NUMERICAL BASELINE RESULT ONLY` | repaired G witness and current regional attribution |
| `W_N^G>W_N^B3` | `NUMERICAL BASELINE RESULT ONLY` | repaired witness only |
| arbitrary distribution/network robustness | `NOT CLAIMED / OPEN EXTENSION` | none |
| heterogeneous regions | `NOT CLAIMED / OPEN EXTENSION` | none |
| old exact Krawczyk object | `LOCAL EXACT STATIONARY-ROOT DIAGNOSTIC` | old isolating box only; not equilibrium authority |

No result in the current paper is classified as a global primitive-space theorem or a generic function-class theorem.

## 6. Benchmark-definition audit

**PASS.**

### G

Decentralized full game. Each government chooses only its own `x_i in [0,1]`; the private intermediary observes public investments and reoptimizes `p_T`; participation then resolves. The repaired baseline provides computational global best-response evidence.

### B3

Matched-price fixed-price identification benchmark. Governments reoptimize their own public investment while the private scalar price is held at the G on-path price. B3 suppresses the repricing response but holds the on-path price level fixed.

**B3 is not:**

- a price-regulation policy experiment;
- a planner problem;
- first best;
- a claim that private price is exogenous in reality.

### National welfare

`W^N=W_1+W_2+Pi_T` is an accounting object. No unrestricted planner problem is solved.

### Local coordination wedge

`dW^N/dx_i>0` at the repaired G witness is a local marginal welfare statement. Correct wording is “local under-provision in the coordinated `+x_i` direction.”

**Prohibited without a new certified planner problem:** `first best`, `social optimum`, `optimal coordinated investment`, `optimal subsidy`, or global underinvestment.

## 7. Claim-scope ledger

| Claim ID | Manuscript location | Verified source | Allowed wording | Prohibited stronger wording | Certificate |
|---|---|---|---|---|---|
| C1 | Abstract; Introduction; Main Results; Conclusion | analytic small-beta proof + Stage 4A identity check | “on regular continuing stationary branches, sufficiently small positive beta yields opposite local BR slopes” | “all nearby equilibria”; “global sign theorem”; “generic reversal” | PASS |
| C2 | Main Results; Discussion; Appendix | Stage 4A repaired audit | “the baseline model contains an all-regime computationally certified global-equilibrium witness with opposite local slopes” | “exact global theorem”; “unique equilibrium”; “global region” | PASS |
| C3 | Main Results; Discussion | route-dominance inequality | “`kappa_L+tau>v+alpha` is sufficient to eliminate nonresident rival-public use” | “necessary for reversal”; “empirically calibrated” | PASS |
| C4 | Welfare | accounting identity | “private fee cancels as an aggregate transfer under zero real private operating cost” | “private price is welfare irrelevant” | PASS |
| C5 | Welfare; Introduction; Discussion | Stage 7 numerical audit | “local national coordination wedge is positive at repaired G witness” | “first-best underinvestment”; “global social optimum requires more x” | PASS |
| C6 | Welfare | Stage 7 numerical audit | “G welfare exceeds B3 at this witness” | “endogenous repricing raises welfare”; “G welfare-dominates B3” | PASS |
| C7 | Robustness | none beyond baseline/local theorem | “not established; extension direction” | arbitrary `C^1`, nonlinear-network or asymmetry robustness | PASS by downgrade |
| C8 | Appendix only | exact L3-1 Krawczyk verifier + Stage 4A counterexample | “exact local stationary-root diagnostic for rejected vector” | “canonical equilibrium certificate”; “SPNE certificate” | PASS after reclassification |
| C9 | Literature / contribution | Stage 6 novelty re-kill | “matched-price public-public sign comparison under shared private repricing” | generic first result that feedback turns complements into substitutes | PASS |

## 8. Theorem Certification Checklist application

### Headline analytic theorem C1

- exact claim identity: PASS;
- quantifiers: local / exists epsilon / every beta on continuing regular stationary branches — PASS;
- parameter/function domain: baseline functional form and stated regularity — PASS;
- local/global distinction: corrected to local stationary geometry — PASS;
- sufficient versus necessary: explicitly sufficient only — PASS;
- uniqueness: not claimed — PASS;
- independent support: Stage 4A re-derived the key beta-zero identities independently — PASS;
- manuscript mapping: abstract/introduction/results/conclusion corrected — PASS.

### Repaired global witness C2

- exact claim identity: computational existence witness — PASS;
- complete public choice domain `[0,1]`: checked — PASS;
- corners/regimes/low-investment attacks: checked — PASS;
- downstream private continuation reoptimized after deviations: checked — PASS;
- unresolved states treated fail-closed: checked — PASS;
- independent evaluator distinct from Stage-4 construction: checked — PASS;
- numerical evidence not upgraded to general analytic proof: corrected — PASS.

### Welfare C5/C6

- national welfare accounting: PASS;
- fee transfer cancellation: PASS;
- planner objective: no planner problem claimed — NOT APPLICABLE with explicit reason;
- first-best label: prohibited and absent from affirmative claim — PASS;
- local versus global welfare: explicit — PASS;
- one-witness G/B3 ranking not generalized — PASS.

### Broad robustness

No headline broad function-class robustness claim remains. Therefore a counterexample-function theorem certificate is **NOT APPLICABLE** to the final claim set; the absence of such a theorem is stated explicitly rather than hidden.

## 9. Required wording downgrades completed in this stage

The following pre-audit formulations were removed or narrowed:

1. `equilibrium branches` in the analytic theorem → `regular interior stationary branches`;
2. old exact rational pair described as `equilibria` → archived `stationary-root diagnostic`;
3. old exact certificate as headline numerical authority → repaired all-regime Stage-4A witness;
4. “uniform distribution is not locally essential” → no certified function-class robustness claim;
5. nonlinear-network persistence → open extension;
6. small-asymmetry persistence → open extension;
7. old 20-draw perturbation evidence → non-authoritative;
8. old welfare numbers/wedge → repaired Stage-7 values;
9. first-best/social-optimum implications → explicitly prohibited;
10. direct public-hub user competition at repaired witness → explicitly disclaimed because the rival route is dominated under the repaired `tau` condition.

## 10. Active result-pipeline correction

The pre-audit `scripts/generate_results.py` still imported the rejected central-regime Stage-4 solver and old Stage-7 welfare script. This created a reproducibility hazard: a future `make tables` could silently restore rejected numbers into the manuscript.

Stage 7.5A reroutes the active result generator to `stage4a_v21_repaired/code/independent_repaired_audit.py` and updates the generated result/table layer to the repaired witness. A dedicated scope lint rejects old headline values and old exact-equilibrium wording in active manuscript sections.

This is a reproducibility/scope correction only; it introduces no new model or theorem.

## 11. Earliest-stage rollback requirement

**NONE.**

The Stage-4A counterexample and repaired witness were already incorporated upstream. Stage 7.5A found no new false theorem after the local/global distinction is stated correctly. All required repairs are wording, evidence-role classification, and active-result-pipeline routing allowed inside Stage 7.5A.

A future attempt to restore any of the following would reopen the indicated stage:

- claim the old exact root is global equilibrium → Stage 4A;
- claim arbitrary distribution/network/asymmetry robustness → Stage 4/7.5A depending on theorem content;
- claim first-best/global welfare optimum → Stage 7;
- broaden novelty to generic feedback-induced sign reversal → Stage 6.

## 12. Stage 8 input package

Stage 8 must freeze exactly the following hierarchy:

### Analytic headline

Local sufficient-condition theorem on regular interior stationary branches:

`exists epsilon>0 such that for beta in (0,epsilon) on the stated continuing branches, BR_B3'>0>BR_G'`.

### Baseline global witness

- `beta=.01`
- `gamma=.825`
- `tau=.35`
- `x_G≈.83710224`
- `x_B3≈.82589039`
- `p_G≈.01840368`
- `BR_G'≈-.01986`
- `BR_B3'≈+.003968`

Classification: all-regime computational global-equilibrium existence witness under baseline functional form and documented tolerance.

### Welfare

- transfer cancellation: exact baseline accounting identity;
- repaired local national coordination derivative: approximately `+.49045`;
- repaired `W_N^G-W_N^B3≈+.01086`: numerical witness only.

### Permanent exclusions

- old vector as global equilibrium authority;
- old 20-draw robustness;
- broad distribution/network/asymmetry theorem;
- first-best/global-optimum language;
- direct public-hub poaching interpretation of the repaired witness;
- generic novelty claim that third-party feedback can turn complements into substitutes.

## 13. CI evidence

Dedicated workflow: `.github/workflows/stage75a-v21-quantifier.yml`.

The workflow regenerates the repaired result/table layer, applies `scripts/stage75a_scope_audit.py`, and requires the checked-in generated files to reproduce exactly. Final run ID and conclusion are to be recorded after the report commit.

## 14. Canonical verdict and routing

**GO — GENERALITY / QUANTIFIER CERTIFICATION PASS.**

Route to **Stage 8 — Canonical Theory Freeze**.

Stage 8 may freeze only the claim hierarchy and repaired evidence package in Section 12. No broader theorem, globality, robustness, planner, welfare, or direct-hub-competition interpretation may be introduced without reopening the relevant earlier gate.
