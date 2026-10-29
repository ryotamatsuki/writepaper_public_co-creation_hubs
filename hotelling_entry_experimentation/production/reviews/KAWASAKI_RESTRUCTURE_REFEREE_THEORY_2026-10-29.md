# Independent economic-theory referee report

Date: 2026-10-29 JST. Role: economic theory referee, separate from manuscript construction. Reviewed the revised active `manuscript/main.tex` inputs, Online Supplement, frozen 37-item ledger, and research-workflow Stage 6/11 rules. This is a manuscript-facing proof/scope attack, not a replacement for the independent baseline execution report or a complete fresh proof of all price-game certificates. No frozen mathematics, code, primitives or numerical results were edited.

## Verdict

**NO-GO for submission readiness on the package as read.** Editorial reorganisation can proceed and is substantially clearer. The central derivative and primitive sign calculation reconstruct correctly, but a pre-existing strictness/hypothesis gap must be recorded and resolved at the appropriate mathematical gate before claiming a fully closed submission package. Theorem-level novelty also remains unclosed. Neither concern is repaired merely by rearranging sections.

This verdict does not invalidate the fixed positive-baseline numerical certificates. It separates (i) preservation of the original result, (ii) correctness on the proof's explicit positive-baseline domain, and (iii) the scope of the printed proposition.

## Independent reconstruction of the central result

1. With full partner participation, private willingness to pay has intercept `U=u+gamma N_j omega R0` and slope `D=d+eta lambda`, with `eta=gamma N_j omega Pi/2`. The interior monopoly cutoff is `(D-U)/(2D)`; hence its lambda derivative is `U eta/(2D^2)`. The direct-service value derivative must retain this private-price response: a social envelope argument would be incorrect.
2. Differentiating the printed threshold gives the two terms in `eq:diagnostic-derivative`. The direct-service/baseline term multiplies the cutoff derivative; the technological term contains both direct lambda variation and the squared-cutoff response.
3. At zero matching, the cutoff tends to `(d-u)/(2d)`. Keeping the first-order direct-service response and direct technological response gives the coefficient
   `M(d-u)^2 + omega Pi(d^2+3u^2)`
   multiplied by `gamma N_i N_j/(16d^2)`. Expanding the intermediate expression in the supplement verifies that the omega coefficient reduces to `d^2+3u^2`, not `(d-u)^2` alone. The `CS0 ell'` and lambda-times-ell-prime corrections are second order.
4. Strict positivity of that coefficient is equivalent to the lower profit-share bound. The upper bound makes `K<0`. These are two different economic restrictions, rather than a welfare sign assumed into the conclusion.
5. Compact feasible lambda range, fixed rho and fixed basic value, and denominators uniformly away from zero justify the uniform remainder. For fixed strict share margins, the positive first-order coefficient dominates the remainder at sufficiently small positive gamma. This is not uniform in profit shares approaching the lower boundary, arbitrary matching intensities, or expanding lambda intervals.
6. The limit of the public ceiling is the stated direct-service expression. The extra restriction `0<u<2d/3` makes it positive; `0<u<d` alone does not establish a positive ceiling.

The main bridge from primitives to conditional policy value is visible and reconstructable without opening production code. The broad research question is answered only in this conditional sense plus the separately certified examples.

## Major mathematical boundary issue: strict baseline incidence

**MAJOR BUT FIXABLE; pre-existing frozen issue, not introduced by this revision.**

The printed Proposition 1 lists nonnegative `R0` and feasible probabilities, and concludes strictly `S_i^0<0`. Its supplement proves this step explicitly **for `p_b>0`**, using

`S_i^0 = p_b K_i - omega_i zeta^2 Pi^2/(2k)`.

The model text as read does not explicitly impose `p0>0` or `zeta>0`. The boundary `p0=zeta=0` is compatible with the displayed success-probability feasibility restriction and nonnegative `R0`; it gives `p_b=R0=S_i^0=0`. With a feasible nonempty profit-share band, the private screening slope can still be positive and `K<0`. Thus strict negativity of `S_i^0` does not follow on that boundary. At lambda zero the corresponding added venture component is also zero rather than strictly negative.

This is a logical boundary test using existing primitives, not a proposed alternative calibration, simulation or extension. It leaves the main first-order derivative argument intact and does not affect the original baseline, whose `p0` is strictly positive. The theorem freeze prevents this referee from silently adding a hypothesis or weakening a conclusion. Record the mismatch, establish whether a positive-baseline restriction already exists in the authoritative frozen domain, and route any required statement/domain correction to Stage 4/4A with fidelity checks. Prose should not conceal the issue.

## Other numbered results and scope

| Result | Independent assessment |
|---|---|
| Lemma 1 | Full-coverage incidence follows from fixed moments and profit shares summing to one. It is not an endogenous participation comparative static. |
| Lemma 2 | Conditional monopoly derivatives follow from the interior optimum. Positive partner skill moment and appropriation are essential; no generic P/P equilibrium derivative follows. |
| Corollary 1 | Uniform-tail moment identity and fixed-cutoff technological difference are correct. The identity also demonstrates why mass and skill moment are not independent state dimensions here. |
| Lemma 3 | Local service/matching increments and the national-minus-local wedge reconstruct directly; real private fixed costs cancel between the same two configurations. This does not check all entry flips. |
| Proposition 1 | Derivative, weak-matching sign band and positive-ceiling restriction reconstruct. Strict baseline incidence needs the positive-baseline condition discussed above. |

A source comparison against commit `a463d3340d8e6c2307136a6d64db484a5fefc5fa` confirmed the four preparatory numbered blocks verbatim in their new destinations. The Proposition 1 mathematical statement and proof content remain intact; a subsequent prose correction from “open” to “closed” rectangle is an accurate description of the already printed inclusive interval. Stable original result identifiers cause Lemma 2 to appear before Lemma 1; this is harmless mathematically but should be disclosed as deliberate identifier preservation.

## Certificates, corners and negative evidence

- The private-only formula treats the full-participation kink, and the heterogeneous/pooled comparison reprices each technology and includes full-coverage corners. No interior fee formula should be extrapolated beyond its stated branch.
- The finite-gamma closed rectangle certifies a conditional slope and value, fixing partner coverage. The three baseline entry examples require distinct global private/private price and all-profile deviation certificates. The revised manuscript keeps these layers separate.
- Vanishing local/national wedge is only equality for the same conditional comparison. It is not sufficient for equality of the full decentralised and planner configuration choices. The revised welfare introduction expressly says so.
- The planner compares the existing sixteen configurations under decentralised continuations. It is not a first best or an optimisation over entry subsidies, fees, research allocations or buildings.
- The narrow cost band and tiny joint persistence box support strict-example persistence. They do not establish economically broad equilibrium comparative statics.
- Beta and Bertrand records are correctly presented as adverse portability evidence; the latter is not claimed to be independently rerun. These failures constrain configuration portability, and should not be misreported as a universal distributional refutation of the conditional slope theorem.
- Lean's limited theorem map must remain distinct from full-game verification. This referee did not rebuild Lean or re-execute the interval programs; those actions belong to the separately documented verification work.

## Novelty and economic importance

The application-neutral reduction is a fee-selected participation game with welfare terms in participation products and a nonlinear quality index. Eliminating uniform cutoffs gives `q=phi(m)`, and national-minus-local welfare is an ordinary externality increment. This absorbs claims of independent two-state network innovation or a new decentralisation-wedge principle.

The remaining plausible contribution is the primitive region combining adverse local commercialisation with a positive public-access-value response to complementarity after private repricing. Its algebraic derivation is more specific than the generic observation that private and social returns differ. Whether it is an independent theorem or an application-specific corollary of general platform/public-option results remains a Stage-6 source and absorption question. This review does not certify literature priority by failing to find an identical formula.

The broad market-structure question remains only partly answered: the Cournot primitives generate `M` and `Pi`, but no general signed incumbent-count, population, demand or cost comparative statics are supplied. The revised text accurately recognises this limit. The small numerical effect also weakens a claim that the discrete configuration switch alone is economically central. Diagnostic value is defensible; implementation or optimal-instrument value has not been established.

## Editorial disposition

| Dimension | Judgment |
|---|---|
| Structure | Improved: private allocation, local boundary, national comparison and bounded technology evidence have distinct roles. |
| Visibility of existing contribution | Improved without adding a theorem. |
| Actual novelty | Unresolved; no automatic gain from reordering. |
| Theoretical rigour | Central bridge reconstructs, with the positive-baseline strictness issue requiring explicit closure. |
| JE fit | Applied-theory topic is potentially relevant; this review cannot certify sufficient novelty, importance or acceptance. |
| Unresolved before submission | Close the printed-domain/proof-domain mismatch, complete evidence-scoped theorem absorption, and reconcile the diagnosis-only contribution with an editor's expected economic payoff. |

Minor presentation issues observed in the interim version: consecutive subsection headings immediately before Proposition 1; terminology should say closed certified rectangle; stable nonsequential lemma numbering should be intentional. These can be repaired editorially. The mathematical boundary issue must remain visible in the audit even if the frozen manuscript cannot be altered within this task.
