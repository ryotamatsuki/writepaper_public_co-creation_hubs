# Independent mathematical-scope closure for the minimal Proposition 1 correction

Date: 2026-10-29 JST. This follow-up responds to the author's explicit authorization to make the minimum necessary correction to Proposition 1 and verify it. Ownership: this report only; manuscript/code implementation is performed separately. This does not retrospectively change the earlier referee's reading or certify literature novelty.

## Recommended minimum correction

Add **strictly positive baseline success probability `p_b>0`** to Proposition 1's maintained hypotheses. `p_b` is the existing object in `eq:effort`, not a new parameter, function, agent, instrument or model extension. The correction narrows the theorem's printed domain to the positive-baseline domain already used explicitly by its proof. It leaves every displayed formula, weak-matching coefficient, cutoff, certified parameter value and numerical conclusion unchanged.

Exact proposed replacement in `sections/04_propositions.tex`:

> Maintain the model's uniform and perfectly aligned technical-skill/premium-value distribution, a positive interior Cournot continuation, an interior local private-only fee, partner full participation, `$M_i<0$`, `$0<u_i<d$`, `$\Pi>0$`, nonnegative `$R_0$`, and strictly positive baseline success probability `$p_b>0$` as defined in equation~`\eqref{eq:effort}`.

The subsequent explanation “`S_i^0=...<0` when `K_i<0`” should read “...when `K_i<0` and `p_b>0`.” This is hypothesis fidelity, not a new mathematical result. The supplement already states the required positivity at the relevant proof step; add a short explanatory sentence there identifying it as a now-explicit maintained hypothesis, without changing the existing equality or proof.

Suggested supplementary sentence:

> The maintained condition `$p_b>0$` is explicit in Proposition~1; without it, the boundary `$p_0=\zeta=0$` gives `$S_i^0=0$`, so strict negativity cannot be inferred from `$K_i<0$` alone.

The manuscript's abstract and broad research question need not be altered: the added restriction is a hypothesis of the specific analytical result, and the existing finite-intensity rectangle already lies within it. The revision record must disclose the authorized exception to the earlier complete-theory freeze rather than claiming literal unchanged hypotheses.

## Why a weak inequality alone is not the minimum coherent repair

Changing only `S_i^0<0` to `S_i^0<=0` would leave other strict claims unjustified on the zero-baseline boundary. At `p0=zeta=0`, both `p_b` and `R0` vanish. The venture component is then zero at lambda zero, even if `K<0`, while the existing conclusion and prose describe a strictly negative increment throughout a feasible lambda interval that can include zero. Retaining that larger domain would require coordinated changes to the strict-incidence conclusion, its interpretation and related statements. It is a legitimate different theorem version, but is a larger change than making the proof's existing `p_b>0` domain explicit.

## Complete proof chain under the corrected hypothesis

1. The profit-share lower bound makes the leading weak-matching derivative coefficient strictly positive. The upper bound gives `K_i=M_i+omega_i Pi<0`.
2. Existing primitives imply exactly `p_b Pi-R0=zeta^2 Pi^2/(2k)>=0`, with `k>0` and `Pi>0`.
3. Hence `S_i^0=p_b K_i-omega_i zeta^2 Pi^2/(2k)<0`: its first term is strictly negative because `p_b>0` and `K_i<0`; its second term is nonpositive because founder profit shares are nonnegative.
4. The private cutoff remains positive on the maintained interior branch. Thus the baseline venture increment is strictly negative. The lambda-dependent increment is nonpositive because lambda is nonnegative and `K_i<0`. Their sum remains strictly negative even at lambda zero.
5. The derivative proof does not otherwise change: `D=d+eta lambda`, `U=u+gamma N_j omega R0`, and `ell'=U eta/(2D^2)` produce the existing exact derivative and first-order coefficient. Fixed strict profit-share inequalities, a fixed compact feasible lambda interval, and denominators bounded away from zero yield the uniform positive slope for sufficiently small positive gamma.
6. The separate `0<u<2d/3` restriction yields a strictly positive limiting service ceiling and therefore a positive ceiling for sufficiently small gamma. The added positive-baseline condition does not alter that limit.

No extra sign condition on `S0`, no change to the matching technology and no adjustment of fixed costs is needed to close this specific gap.

## Remaining domain discipline

The correction does not remove the already maintained requirements: positive interior Cournot quantities; nonnegative feasible success probabilities; the adopted optimal-effort continuation; positive service and premium values; nonnegative profit shares summing to one; interior private-only fee; partner full coverage; a fixed compact feasible lambda interval; and weak matching obtained by varying mu while holding rho/basic value fixed.

For the narrow sign proof, no further new hypothesis is required after `p_b>0`, provided those existing continuation restrictions are enforced. In particular, `0<u<d` makes the gamma-zero fee cutoff strictly interior, which persists uniformly for small gamma on a fixed compact lambda interval. `D` stays bounded away from zero by positive `d` for nonnegative gamma, lambda and appropriation. No uniform claim is established as profit shares approach the threshold, lambda domains expand, or probability/effort branches change.

The effort formula is the maintained admissible continuation; it must not be extrapolated to a negative unconstrained effort or a probability-cap regime where that formula ceases to solve the constrained problem. The frozen supplement already excludes the success-probability corner when deriving effort. This audit finds no need to enlarge that restriction or change the effort model to repair the incidence claim. It does not certify all unrestricted primitive values outside the maintained continuation.

## Formal-statement fidelity

Whole-repository search located the formal project at:

- `formal/Stage33Formal.lean`
- `formal/lakefile.toml`
- `formal/lean-toolchain`
- `theorem_certificates/STAGE33_LEAN_FORMAL_VERIFICATION_CERTIFICATE.md`

The Lean file contains eighteen theorem declarations. They cover uniform moment identities, conditional fee/cutoff ordering, the local/planner accounting wedge, public-entry polynomial identities, binary-flip logic and exact fixed-baseline polynomial signs. It **does not encode** Proposition 1's general weak-matching expansion, primitive profit-share interval, `p_b`/`R0` definitions or strict `S0` consequence. The pre-existing certificate explicitly excludes end-to-end market-primitives mapping and the full game.

Accordingly, adding `p_b>0` does not alter or invalidate any encoded Lean theorem statement. No new Lean theorem is necessary merely to preserve the fidelity of those existing targets. The revised scope ledger must state: **Proposition 1 is analytically proved with a separate exact-Python finite-baseline certificate; it is not a Lean-verified general proposition.** An actual Lean rebuild is a separate execution claim and may only be reported if completed. Unchanged source hashes and the prior successful build establish preservation, not a fresh build.

## Existing Python coverage and exact audit

`government_game_exploration/incidence_screening_sign_certificate.py` fixes positive `p0`, reconstructs `Pi` from the original Cournot primitives, and computes `pb=p0+zeta*zeta*Pi/k`; its original certificate therefore already satisfies the corrected hypothesis. Its interval checks establish negative `K`, negative `S0` and positive conditional derivative across the fixed closed rectangle.

`incidence_screening_valuation_check.py` uses the same fixed positive baseline and certifies positive public value/direct service and negative venture increments on the same rectangle. `revision_p1d/validate_revision.py` verifies protected sources and runs older baseline/diagnostic checks; by itself it does not check the new printed hypothesis or encode a general Proposition 1 theorem.

This reviewer executed a non-writing Python exact-`Fraction` check importing the existing Cournot reconstruction. It verified positive baseline `pb`, the identity `pb*Pi-R0=zeta^2*Pi^2/(2k)`, and the incidence identity/negative signs at the fixed rectangle endpoints and original appropriation share. It also verified the zero-baseline boundary equals zero. Output:

> PASS exact fixed-baseline positive pb and incidence identity; p0=zeta=0 strictness boundary equals zero. No files generated.

That check validates the missing logical link at the unchanged baseline; it is not a new simulation, a new certified parameter rectangle or a general computer proof. The original full-rectangle scripts should be rerun unchanged by the verification owner, and their outputs compared with the frozen outputs.

## Provisional disposition pending implementation inspection

The minimum correction is mathematically sufficient for the identified strict-incidence gap, with no baseline-number changes. After checking the actual manuscript edit and proof-domain annotation, this specific Stage-4/4A scope objection can be closed. Submission readiness still depends independently on novelty absorption, editorial significance and the other reviewers' concerns. A narrow mathematical repair does not manufacture a new substantive economic contribution.

## Implemented-edit inspection and final scope verdict

The root's subsequent manuscript edits were inspected directly. Proposition 1 now explicitly assumes a positive baseline success intercept `p_b>0` and points to its existing definition. The supplement's proof introduction explicitly retains `p_b=p0+zeta^2 Pi/k>0` as a requirement of that proposition. A new explanatory paragraph excludes the degenerate zero-baseline case and states why the strict sign follows. The existing proof's strict-incidence step already lists `p_b>0`; its statement and use are now aligned.

An independent Python extraction compared every equation/align/gather/display block in the edited proposition file and supplement against the pre-repair `HEAD`. All displayed mathematics was identical. The formal source was byte-identical, with SHA-256 `abbe081e4e6c96f49271a7ab41e722b583bfdf59ffdb477e5b50d2c72704d32e`. This reviewer also rechecked the unchanged exact baseline positivity/identity as documented above. Full baseline regression, LaTeX builds and any actual Lean rebuild remain execution claims for the verification owner.

**GO — identified mathematical-scope repair closed.** The explicit hypothesis is sufficient for the previously unsupported strict-incidence assertion on the maintained theorem domain. This is a domain narrowing using an existing derived object, not a new model or a new economic result. The earlier referee objection is retained as historical diagnosis and superseded on this specific point by this closure. No unresolved mathematical defect from this particular boundary attack remains after the edit.

This is not a blanket Stage-4/4A recertification of every possible model branch, a Lean certification of Proposition 1, a novelty certificate, or a journal-acceptance recommendation. Separate submission and editorial judgments must consume this narrower closure without treating it as a change in their own evidence.
