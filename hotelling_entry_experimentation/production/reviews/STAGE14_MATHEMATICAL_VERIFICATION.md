# Stage 14 — Mathematical and reproducibility verification

**Frozen source:** 7308052c0f6299597eba4605ad9b4906966eb8f8. No substantive change to any proposition, assumption, proof, numerical primitive or script.

## Formal statements and evidence boundaries
- Lemma 1: local/national regional-incidence accounting under full founder participation, fixed positive interior Cournot quantities and uniform independent skills; not the fee/entry response.
- Lemma 3: conditional public basic-access cost ceiling given partner full participation and local interior private-only fee; the national difference adds the foreign increment. Does not solve all four entry decisions.
- Proposition 1: positive success baseline `p_b>0`, positive interior Cournot continuation, aligned uniform skills/premium value, feasible compact `lambda` interval, interior local private-only fee, fully participating other region, `M_i<0, Pi>0, 0<u_i<d`, nonnegative `R_0`, nonempty open profit-share band and sufficiently weak positive matching. The strict conclusion `K_i<0, S_i^0<0, T_i'(lambda)>0` is conditional. The strictly positive threshold level additionally needs `u_i<2d/3`.
- Four-party configurations: three previously certified finite baseline examples and 16 configuration comparisons, not a general 16-configuration comparative-static theorem; cost band narrow, Beta(2,2) and differentiated-Bertrand limitations retained.
- Supplement proof mapping: `manuscript/appendix/appendix_a_derivations.tex` is included by `manuscript/supplement.tex`. No new proof was added.

## Independent spot-check with SymPy (Stage 14 execution)
With `U=u+gamma N_j omega R_0`, `eta=gamma N_j omega Pi/2`, `D=d+eta lambda`, `ell=(D-U)/(2D)`, `C=gamma N_iN_j`, `K=M+omega Pi` and the stated local threshold:
1. `simplify(diff(T,lambda)-reported_exact_derivative) == 0`.
2. `simplify(diff(diff(T,lambda),gamma).subs(gamma,0)-N_iN_j*(M*(d-u)^2+omega Pi*(d^2+3u^2))/(16 d^2)) == 0`.
3. `simplify(T.subs(gamma,0)) == N_i*u*(2*d-3*u)/(8*d)`.
These algebraic identities are a new *verification* of existing formulas, not new propositions or tests of every feasibility domain. The `p_b>0` strict sign was also previously checked by committed SymPy JSON (SymPy 1.14.0). No unrestricted SPNE inference.

## Recorded reproducibility evidence / CI
- Existing `revision_je_upgrade/verification_checks.json`: ten old scripts on baseline and revised source, zero exit codes; eight output SHA-256 hashes unchanged. The previous run's environment was described as Python 3 with SymPy 1.14.0 for the new scope regression. This Stage 14 did not rerun those ten scripts locally.
- Baseline commit checks completed: `verify` at https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/40780728568/job/123141457214 ; `hypothesis-and-baseline` at https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/40780728653/job/123141456809 ; `preserved-formal-core` at https://github.com/ryotamatsuki/writepaper_public_co-creation_hubs/actions/runs/40780728653/job/123141457418 .
- Existing Lean source: `formal/Stage33Formal.lean` contains eighteen theorem declarations, including selection, fee/cutoff, accounting wedge, polynomial and binary-flip results. The GitHub workflow rebuilds the pinned existing formal core with `lake build --wfail`, but **does not formalise Proposition 1**, its profit-share band or the full four-party equilibrium.
- Existing previous builds: `build_checks.json` reports 37/8 pages and no undefined/overfull/underfull warnings. A fresh independent compile and visual review are still required.

**Mathematical verdict:** no rollback-triggering defect found in the checked identities or current specified proof domains. This does not assert an exhaustive new theorem proof.

## Sources
- JE official instructions (accessed 2026-10-29 JST): https://link.springer.com/journal/712/submission-guidelines
- JE publication options (accessed 2026-10-29 JST): https://link.springer.com/journal/712/how-to-publish-with-us
- PR #73 head (read 2026-10-29 JST): 7308052c0f6299597eba4605ad9b4906966eb8f8
- Workflow v1.1 tag: 488e5ab06c207909296a7564eaf9066f7f94319c

