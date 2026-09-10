# Stage 11 v2.1 — Independent Robustness / Referee Attack Gate

Date: 2026-09-10 JST

Workflow authority: `ryotamatsuki/research-paper-workflow` v2.1 @ `b27568172d4145a2b98825a5b66a071dbfe25f36`.

Canonical template: `templates/STAGE_11_REFEREE_GATE.md` @ blob `232a7efea21886b7bad6f061d6989a87501849db`.

Checklist: `checklists/REFEREE_ATTACK_CHECKLIST.md` @ blob `44b7e9918976de4fc157d932532bc45cb4d8123b`.

Stage-10 baseline / parent: `02bcd37ac1ab019dc6df2c0d29c8f312a1c992a4`.

Stage-8 freeze: `ad927ca783a6123ea4fc6f55f65598ebd6ab583b` plus the Stage-11-triggered T3 amendment `theory_freeze_v21/STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md`.

Stage-7.5A authority: `reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md` plus the limited-reopen correction `reviews/STAGE_075A_V21_T3_SYMMETRY_LIMITED_REOPEN_2026-09-10.md`.

Working title: **Strategic Interaction among Public Innovation Hubs under Private Repricing**.

Journal target: **NOT SELECTED — deferred to Stage 12**.

This audit is independent of the separately commissioned Astra hostile-referee review. No Astra conclusion is used below.

## 1. Executive referee-gate verdict

**PROVISIONAL: CONDITIONAL GO — PR-HEAD CI REQUIRED.**

No fatal mathematical, equilibrium, welfare, or novelty failure has been found. Two material bounded issues were found and repaired:

1. **CERTIFICATION REGRESSION — T3 symmetry anchor under-specified.** T2 proves the negative beta-zero G cross effect at a symmetric regular full-game stationary state, while the old T3 prose referred only to an arbitrary regular G stationary branch. Stage 7.5A was therefore reopened in a bounded way and T3 was narrowed to a G branch through the symmetric beta-zero state covered by T2.
2. **Closest-prior-art omission — Kim (2024).** `Mixed Duopoly in Two-Sided Markets` combines a welfare-maximizing public platform, a profit-maximizing private platform, two-sided demand, endogenous prices, and first-stage R&D/investment. It does not contain two decentralized public investors sharing a third private alternative or the matched-price public-public best-response sign comparison, but it is close enough that omission would invite a referee objection. It is now cited and distinguished directly.

No model primitive, T1/T2 formula, repaired global witness, B3 definition, welfare object, or novelty envelope is expanded.

Final success criterion after CI: **GO TO JOURNAL POSITIONING**.

## 2. Referee A — novelty / mechanism

### Attack A1 — classic-result / combination-novelty attack

**Attack.** Strip labels and compare the paper to two-sided platform competition, mixed public-private platform competition, sequential pricing, and strategic public investment.

**Severity.** MAJOR RISK, currently defensible.

**Evidence.** The literature already contains: two-sided strategic interaction; public-private mixed platform competition; endogenous fee responses; investment followed by downstream price competition; and qualitative reversals induced by follower responses. The paper therefore cannot claim novelty for any component.

**Surviving object.** The narrow remaining object is: two decentralized public investors share one private fee-setting outside option; B3 fixes that private fee at the G on-path value while governments reoptimize; restoring the shared private price policy reverses the sign of the same public-public local best-response slope.

**One-page-specialization kill test.** No inspected closest paper supplies this exact three-actor matched-price public-public sign comparison by direct relabeling or a one-step parameter restriction. Kim (2024) is especially close institutionally but has one public and one private platform choosing investment and then prices; it does not generate the present two-public strategic object.

**Result.** PASS, with **incremental/combination novelty** as the strongest surviving editorial risk.

### Attack A2 — no-new-mechanism attack

The negative private-price response and positive cross-side feedback are known ingredients. The paper survives only because the comparison changes the qualitative strategic relationship between two other actors rather than presenting either channel as new. Current Introduction and Related Literature respect this ceiling.

**Result.** PASS WITH NARROW CONTRIBUTION.

## 3. Referee B — assumptions / mathematics / globality

### T1 — First-order B3 complementarity

Independent symbolic reconstruction from primitive beta-zero participation and welfare confirms:

- `M_B3(0)=0` by separability at beta zero;
- `dM_B3/d(delta)|0 = alpha^2 d^3/(Delta_i^3 Delta_j^2)>0`;
- with `delta=alpha beta`, `dM_B3/dbeta|0 = alpha^3 d^3/(Delta_i^3 Delta_j^2)>0`.

The fact that `M_B3(0,x,p)` is identically zero on the beta-zero central regime makes the chain-rule multipliers on the stationary/matched-price path derivatives vanish at first order.

**Verdict: PASS.**

### T2 — Beta-zero G cross effect

An independent derivation solves the private beta-zero profit problem, substitutes the private optimum into reduced regional project surplus, and directly differentiates the reduced objective rather than calling the production chain-rule derivation. At symmetry it reproduces

`M_G(0) = -3 T alpha^2 kappa_L^2 / [16 Delta (Delta+T)^3] < 0`.

**Verdict: PASS.**

### T3 — Local strategic sign reversal

**Attack discovered a real quantifier defect.** The negative G anchor is proved only at a symmetric beta-zero state. The historical T3 prose did not explicitly require the G branch to pass through that state.

**Repair.** Stage 7.5A limited reopen and Stage-8 amendment now require a regular G branch through the symmetric beta-zero state covered by T2. With smooth continuation, strict SOCs, nonsingular systems, preserved regime inequalities, and matched-price continuity, take the minimum of the finite local radii supporting B3 positivity and G negativity. This gives a common `epsilon>0`.

**Verdict: PASS AFTER CERTIFICATION-REGRESSION REPAIR.**

No asymmetric-G theorem is certified.

### T4 — repaired computational global-equilibrium existence witness

The Stage-4A evaluator already performs an all-domain public sweep and targeted full multistart continuation. Stage 11 additionally reconstructs the model without importing that solver.

Independent reviewer-side results reproduce:

- private price near `0.01840368`;
- G public BR near `0.83710`, with detected gain over the candidate below the certification tolerance;
- B3 public BR near `0.82589`, with no detected profitable gain;
- opposite local Hessian/BR-slope signs.

**Verdict: PASS AS A COMPUTATIONAL EXISTENCE WITNESS.**

Surviving limitation: this is not an analytic proof that no machine-scale deviation exists, not a uniqueness certificate, and not a parameter-region theorem. The manuscript states these limitations.

## 4. Referee C — welfare / institution / benchmark

### W1 — fee transfer cancellation

With private real operating cost normalized to zero, project payments `-p_T n_T^F` and private profit `+p_T n_T^F` cancel in national welfare.

**Verdict: PASS.**

### W2 — local coordination wedge

Independent numerical reconstruction reproduces the decomposition at the repaired G witness:

- own reduced welfare derivative approximately zero;
- rival welfare derivative approximately `+0.50053`;
- private profit derivative approximately `-0.01008`;
- national derivative approximately `+0.49045`.

The manuscript calls this only local under-provision in the coordinated `+x_i` direction and explicitly rejects first-best/global-optimum language.

**Verdict: PASS.**

### W3 — G versus B3 welfare ranking

Independent reconstruction reproduces approximately

- `W_N(G)=0.5967883`;
- `W_N(B3)=0.5859284`;
- difference `0.0108599`.

The text limits this to one repaired witness.

**Verdict: PASS.**

### B3 benchmark attack

B3 is artificial by design but economically interpretable as a matched-price identification benchmark. It suppresses the private fee policy response while avoiding a mechanical on-path fee-level difference. It is not presented as regulation or first best.

**Verdict: PASS, with editorial caveat that the benchmark must remain diagnostic rather than policy language.**

### Remote-public friction attack

At the repaired vector, partner mass is bounded by one, so a remote public route has quality at most `v+alpha=.60`, while its access cost is `kappa_L+tau=.62`. For every `z<=1`, remote-public utility is at most `-.02<0`; outside option utility is zero.

**Verdict: PASS as a sufficient dominance condition.**

Editorial limitation: the repaired witness is intentionally a high-friction case and is not evidence of direct public-hub user poaching. Current manuscript discloses this.

### Institutional attack

Named hubs support plausibility of paid private access, network/intermediation functions, and public analogues. They do not validate the causal price response or exact objective/timing. The manuscript labels the price response as a model prediction and treats single-homing/profit maximization as reduced form.

**Verdict: PASS.**

## 5. Referee D — journal / exposition / claim scope

The paper is technically coherent after the T3 repair. The main editorial weakness is theorem strength: the analytic theorem is local and the global-equilibrium support is one computational witness. That is acceptable for journal positioning, but likely insufficient for the very top IO/public-economics outlets unless the mechanism is viewed as unusually sharp.

The Stage-10 `NO REQUIRED FIGURE` decision survives. A beta-path or threshold figure would require a certified path/domain not proved by T3. The strategic and welfare tables are more faithful to the actual evidence.

**Referee recommendation before journal fit is selected: SEND OUT / BORDERLINE depending on outlet.**

## 6. Candidate-deviation re-audit

Independent reviewer implementation: `stage11_v21_independent/code/independent_stage11_regression.py`.

It does not import `stage4a_v21_repaired`, `scripts/generate_results.py`, or the analytic derivation scripts.

Attacks:

- reconstruct all four project routes from utility upper envelopes;
- solve partner participation directly;
- optimize the private continuation after public deviations;
- scan the full public interval and locally refine candidate maxima;
- compare repaired candidate welfare to independently detected best replies.

Current result from reviewer-side development: no profitable G or B3 finite deviation detected beyond documented tolerance.

Final state: **PENDING CI RE-RUN OF THE COMMITTED INDEPENDENT SCRIPT**.

## 7. Alternative-equilibrium / multiplicity re-audit

The manuscript does **not** claim uniqueness. Stage 11 nevertheless attacks alternative equilibria separately from candidate deviations.

Reviewer-side sampled best-response maps use rival public choices spanning `0,.25,.50,.75`, the repaired candidate neighborhood, and `1`. No sampled fixed point outside the repaired neighborhood is detected.

This is **not a uniqueness proof** and is not recorded as one. Its only role is to satisfy the hostile search for an obvious alternative branch/corner equilibrium and to check that the paper is not silently relying on uniqueness.

**Verdict: PASS FOR EXISTENCE-SCOPE PAPER; uniqueness remains explicitly unclaimed.**

## 8. Indifference / zero-payoff re-audit

Private prices above the maximum willingness-to-pay bound can generate zero private demand and zero profit, but at the material repaired and adversarial public histories the private intermediary has a strictly positive-profit interior continuation. Hence those zero-profit prices are not payoff-equivalent best replies at the histories supporting the paper.

Project route ties occur at cutoff types of measure zero under the continuous uniform type distribution and do not change aggregate shares/welfare. No tie-breaking refinement is used to remove an inconvenient positive-mass equilibrium.

**Verdict: PASS at the claimed witness/continuation scope.**

## 9. Selection / refinement provenance and symmetry audit

No weak-dominance refinement, no-loss restriction, or asymmetric equilibrium-selection device is introduced to select the repaired witness. Remote-public route dominance follows directly from primitives at the repaired vector. Participation multiplicity is treated fail-closed by multistart comparison.

Symmetry is now used transparently in the only theorem place where it is mathematically required: the beta-zero G anchor for T2/T3.

**Verdict: PASS AFTER T3 REPAIR.**

## 10. Independent equilibrium / continuation re-audit

Adversarial public histories checked independently include

`x_i = 0, .05, .10, .18, .20, .35, .50, .75, .95, 1.0`

against the repaired rival choice, with private repricing and ten dispersed participation starts. No material participation multiplicity or nonconvergence was observed in the reviewer-side reconstruction.

The production Stage-4A audit independently checks the old dangerous low-investment region, boundaries, detected maxima, repaired candidates, route changes, and private repricing, and fails closed on unresolved/multiple continuations.

**Continuation verdict: PASS at the repaired witness certification scope.**

## 11. Welfare-selection regression audit

No relevant equilibrium multiplicity was found by the Stage-4A or Stage-11 witness audits. The welfare statements are one-witness statements and do not claim selection-free rankings over a broader equilibrium correspondence.

**Verdict: PASS.**

## 12. Independent quantifier / function-class re-audit

The manuscript makes no arbitrary-distribution, nonlinear-network, heterogeneous-region, or general function-class theorem. Such variations are explicitly extension directions. Therefore no unsupported function-class quantifier remains to kill.

The one quantifier defect actually found was the missing symmetric G anchor in T3; it has been repaired through the required Stage-7.5A limited reopen.

**Verdict: PASS AFTER REPAIR.**

## 13. Solver-failure / unresolved-continuation ledger

Production Stage-4A semantics: unresolved participation continuation and multiple private/participation equilibria are fail-closed rather than treated as unprofitable deviations.

Reviewer-side committed regression: exceptions/nonconvergence and multistart spread beyond tolerance fail the run.

Current observed material unresolved continuations in the independent development attack: **0**.

Final CI count/status: **PENDING**.

## 14. Evidence ledger for material PASS states

| Claim | Attack actually performed | Evidence / artifact | Result | Surviving limitation |
|---|---|---|---|---|
| T1 | beta-zero separability + independent first-order fixed-point expansion | `stage11_v21_independent/code/independent_stage11_regression.py` | PASS | baseline local branch only |
| T2 | solve private beta-zero optimum, substitute, differentiate reduced welfare directly | same | PASS | symmetric beta-zero G state |
| T3 | quantifier attack against T2 anchor | limited-reopen + freeze amendment + manuscript repair | PASS AFTER REPAIR | local; symmetric G anchor; no explicit epsilon |
| T4 | independent route/continuation/private/public reconstruction | reviewer regression + Stage-4A audit | PASS | computational existence witness, not uniqueness |
| remote route | primitive utility upper bound | reviewer regression + analytic inequality | PASS | sufficient; high-friction witness |
| W1 | accounting reconstruction | symbolic welfare identity / direct accounting | PASS | zero real private operating cost |
| W2 | central-difference component reconstruction | reviewer regression | PASS | one witness/local direction |
| W3 | independent national-welfare computation | reviewer regression | PASS | one witness only |
| novelty | proposition-level comparison including Kim (2024) | related-literature repair / prior-art audit | PASS, narrow | combination novelty risk |

## 15. Certification-regression ledger

### CR-1 — T3 symmetric G anchor

- **Failure:** Stage-7.5A/Stage-8 T3 prose allowed a reading broader than the actual T2 proof anchor.
- **Severity:** MAJOR BUT FIXABLE.
- **Earliest affected stage:** Stage 7.5A.
- **Would have been prevented by:** theorem-certificate check requiring every local sign used in a continuity theorem to identify the exact branch anchor state.
- **Repair:** `reviews/STAGE_075A_V21_T3_SYMMETRY_LIMITED_REOPEN_2026-09-10.md` and `theory_freeze_v21/STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md`, plus manuscript/validator synchronization.
- **Theory change:** NO new result; scope narrowing only.
- **Resolved:** YES subject to CI regression.

No Stage-4A certification regression was found in the repaired witness.

## 16. Consolidated severity table

| Attack | Severity | Status |
|---|---|---|
| T3 omitted symmetric G anchor | MAJOR BUT FIXABLE / CERTIFICATION REGRESSION | REPAIRED; CI pending |
| Kim (2024) closest-paper omission | MAJOR BUT FIXABLE positioning | REPAIRED |
| combination novelty | MAJOR editorial risk | DEFENSIBLE, unresolved as journal-fit risk |
| remote-public friction artificiality | MINOR-to-MAJOR editorial limitation | DISCLOSED / no correctness failure |
| T1/T2 mathematics | — | PASS |
| repaired global deviation audit | — | PASS at computational-witness scope |
| alternative equilibrium search | — | no obvious alternative found; uniqueness not claimed |
| welfare accounting/scope | — | PASS |
| institutional mapping | — | PASS as analogue class |
| figure/table architecture | — | PASS; no figure required |

Fatal issues: **NONE FOUND**.

## 17. Required fixes and earliest affected stage

Completed:

1. narrow T3 to a G branch through the symmetric beta-zero state — earliest Stage 7.5A;
2. propagate the corrected scope through Abstract, Main Results, Robustness, Introduction, Discussion, Conclusion, and validators;
3. add and directly distinguish Kim (2024) in Related Literature;
4. add an independent Stage-11 regression implementation.

Remaining: only CI execution and closeout evidence.

## 18. Theory-change implications

No new player, timing, utility, objective, route, parameter restriction, equilibrium concept, benchmark, welfare object, or novelty claim is introduced.

The T3 correction is a **theorem-domain narrowing mandated by the proof already in force**. It is handled through formal Stage-7.5A limited reopening because the base freeze's change-control explicitly routes theorem-quantifier corrections there.

No Stage 4/4A, Stage 6, or Stage 7 theoretical reopening is required.

## 19. Resolved versus unresolved attacks

Resolved subject to CI:

- T1/T2 independent derivation;
- T3 quantifier regression;
- repaired candidate-deviation attack;
- off-path continuation attack;
- remote-public route dominance;
- welfare reconstruction;
- literature omission;
- benchmark/institutional wording.

Unresolved substantive correctness attacks: **NONE currently identified**.

Unresolved editorial risk: whether the narrow combination novelty is sufficient for a given journal. That is the proper object of Stage 12.

## 20. Verdict and Stage-12 contract

Current verdict: **CONDITIONAL GO — PR-HEAD CI REQUIRED**.

If the dedicated Stage-11 workflow passes the clean full build, independent Stage-11 regression, amended Stage-7.5A scope audit, and freeze-amendment consistency checks, the final verdict becomes:

**GO TO JOURNAL POSITIONING**

Stage 12 must select an outlet for the actually surviving contribution. It may not remove the symmetric-G T3 anchor, promote the repaired witness to uniqueness/general theorem status, convert B3 into policy regulation, or suppress the high-friction/one-witness limitations to fit a preferred journal.
