# Stage 11 Referee B — Assumptions / Mathematics

## Overall recommendation
**MAJOR REVISION IN CERTIFICATE EXPOSITION; mathematical validity survives.** I do not find a fatal equilibrium or accounting error, but a computational existence result should disclose its verification boundary more explicitly.

## Serious attacks

### B1 — FOC versus equilibrium / SOC / interiority
Attack: The reported stationary point may fail the intended equilibrium conditions.
Severity: MINOR.
Evidence: The frozen verifier solves the participation fixed point and public/private problems, rejects invalid central-regime states, and checks negative public own second derivatives. Stage-11 independent rerun of the frozen code reconfirmed near-zero fixed-point/FOC residuals, interior cutoffs and partner masses, negative public and private SOCs, and nonsingular fixed-point/equilibrium Jacobians.
Can the paper answer now? YES.
Required fix: State the certificate checklist in the appendix rather than merely saying the code verifies the result.
Does the fix reopen theory? NO (Class V).
Resolved? YES mathematically; exposition pending pre-fix.

### B2 — Participation fixed-point / local uniqueness
Attack: Local uniqueness is assumed rather than justified.
Severity: MINOR.
Evidence: Stage-8 restrictions require nonsingular fixed-point/equilibrium Jacobians; the frozen computation satisfies the regularity gate. The claim is local, not global uniqueness.
Can the paper answer now? YES.
Required fix: Make the local nature explicit.
Does the fix reopen theory? NO.
Resolved? YES.

### B3 — Numerical-not-proof
Attack: Calling the central statement a proposition disguises a computational example as a theorem.
Severity: MAJOR BUT FIXABLE.
Evidence: Stage-8 status is explicitly `NUMERICALLY SUPPORTED ONLY`; analytic pieces are the derivative decomposition/essentiality/accounting, while existence is certified numerically and extended locally by continuity/IFT conditions.
Can the paper answer now? YES.
Required fix: Label the result as a computationally certified local existence result and spell out analytic versus computational components immediately around the statement and in the appendix.
Does the fix reopen theory? NO.
Resolved? NO pre-fix.

### B4 — Parameter exercise
Attack: One searched parameter vector is not a theory result.
Severity: MAJOR BUT FIXABLE.
Evidence: The signs are strict, the point is regular, 20/20 deterministic ±0.5% perturbations preserve the reversal, and conditional local C1 distribution/network and action-state continuity results are frozen. This is more than one point but less than a global characterization.
Can the paper answer now? YES with precise scope.
Required fix: Present the witness as a constructive certificate, place the analytic dominance condition before the numbers, and frame the perturbation exercise as computational confirmation of local persistence rather than a separate theorem.
Does the fix reopen theory? NO.
Resolved? NO pre-fix.

### B5 — Alternative demand / functional form
Attack: Uniform-linear structure may drive the result.
Severity: MINOR.
Evidence: Only local C1 robustness is established. That is sufficient to reject a literal knife-edge claim but not to support arbitrary demand/network language.
Can the paper answer now? YES if language remains local.
Required fix: None beyond local wording.
Does the fix reopen theory? NO.
Resolved? YES.

### B6 — Symmetry
Attack: The result may rely on symmetric regions.
Severity: MINOR.
Evidence: Small asymmetric action states survive by continuity; no asymmetric-primitives theorem exists.
Can the paper answer now? YES within stated scope.
Required fix: Keep the limitation explicit.
Does the fix reopen theory? NO.
Resolved? YES.

### B7 — Alternative contract/information set
Attack: Richer tariffs or simultaneous instruments could undo the result.
Severity: MINOR.
Evidence: The frozen research object has one public facilitation instrument and one private access price. The paper does not claim contract-set robustness.
Can the paper answer now? YES as a scope limitation.
Required fix: Make clear that richer contracting is outside the claim.
Does the fix reopen theory? NO unless pursued.
Resolved? YES as limitation.

### B8 — Timing
Attack: The effect is just Stackelberg sequencing.
Severity: MINOR.
Evidence: Sequential private repricing is deliberately the mechanism under study and is conceded as known in generic form. The contribution is the induced public-public sign reversal relative to the nested fixed-price benchmark.
Can the paper answer now? YES.
Required fix: Interpret the sequence as the maintained research object, not a universal timing theorem.
Does the fix reopen theory? NO.
Resolved? YES.

### B9 — Precision / tolerance
Attack: Finite differences could flip signs.
Severity: MINOR.
Evidence: Frozen BR slopes are materially separated from zero, derivative signs persist in the deterministic perturbation audit, and production verification is reproducible. Stage-11 independent rerun matches the signs and regularity gates.
Can the paper answer now? YES.
Required fix: Keep computational-protocol transparency in appendix.
Does the fix reopen theory? NO.
Resolved? YES.

## Referee B bottom line
No fatal mathematical defect is located. The publication risk is the **proof ceiling**, not correctness: a theory-oriented referee may demand a global analytic characterization that the frozen paper does not provide.