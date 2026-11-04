# JE contribution-closure revision: independent verification

Date: 2026-10-29 JST. Baseline: PR #73 commit `3d2eaf2` (the completed frozen-results reconstruction). Current scope: the user-authorised explicit positive-baseline hypothesis correction plus exposition and source positioning. No new economic model, simulation specification or numerical parameter is authorised by these checks.

## Verification verdict and actual change

**Corrected mathematical scope regression: PASS. Existing Python reproduction: PASS. Display/table/label preservation: PASS. Final LaTeX/reference/rendering QA: PASS. Local Lean rebuild: NOT RUN (toolchain unavailable).** This verdict concerns the stated repair and reproducibility; literature priority and journal acceptance remain independent editorial questions.

The diagnosed defect is the strict claim `S_i^0<0` at a degenerate zero baseline. The previous proposition admitted `p_0=zeta=0`, although the existing proof explicitly used `p_b>0`. The revision adds `p_b>0` to Proposition 1 and aligns the supplement's opening statement with that existing proof condition. This **narrows the theorem's stated domain**; it is a scientific scope correction, not merely an editorial reordering. The proof's displayed algebra and all reported numerical cases remain unchanged.

The new exact regression is `JE_UPGRADE_SCOPE_CHECK_2026-10-29.py` (Python 3, SymPy 1.14.0); its machine-readable results are in the adjacent JSON file. SymPy was installed in the runtime to execute this test; no frozen dependency or verification file was modified.

## Exact symbolic regression of the diagnosed boundary

The new test checks the primitive definitions directly, rather than relying on numerical substitution:

1. Substituting `p_b=p_0+zeta^2 Pi/k`, `R_0=p_0 Pi+zeta^2 Pi^2/(2k)` and `K=M+omega Pi` into `S^0=p_b M+omega R_0` gives exactly `S^0=p_b K-omega zeta^2 Pi^2/(2k)`; SymPy's simplified residual is exactly zero.
2. Substituting `p_0=zeta=0` gives exactly `S^0=0`, reproducing the failure of the former strict assertion. This is the regression case that the old statement failed.
3. Writing `K=-loss`, with `loss>0`, `p_b>0`, `omega>=0`, `Pi>0`, `k>0` and `zeta>=0`, SymPy confirms the resulting expression is strictly negative. It also verifies both primitive routes to a positive baseline: `p_0>0` or `zeta>0` under the other maintained nonnegative restrictions.

These are checks of the previously used proof bridge and corrected hypothesis. They add no independent economic proposition, new technology, or stronger four-party equilibrium claim.

## Mathematical and archive invariance

The checker recursively expands current main/supplement inputs and compares with PR #73. It verifies all **48 main** and **26 supplement** displayed mathematics instances as identical multisets after whitespace normalisation, all labels as identical multisets, and all tabular bodies as identical multisets. The formal Proposition 1 statement differs **exactly** by the authorised `p_b>0` hypothesis insertion: all other statement sentences and displayed conclusions are identical. The supplement adds the same condition inline at its proof opening; existing proof algebra is retained.

All **28** baseline-tracked scientific `.py`, `.lean`, `.csv` and `.json` files in production/formal were compared byte-for-byte and are unchanged. All **18** PR #73 archive/report files covered by the `revision_kawasaki_results/` and `reviews/KAWASAKI_` prefixes were likewise compared and remain unchanged. New JE reports and new PDF outputs are stored under distinct names.

## Existing checks rerun against both PR #73 and the revision

A detached PR #73 worktree and an isolated copy of the revised sources were used. This prevents regenerated verification outputs from replacing frozen output/archive files in the working branch. For each tree, `python3` executed each unchanged script below. All twenty executions (ten per tree) exited successfully.

| Existing script, relative to production | Reproduced evidence |
|---|---|
| `reviews/stage11r_pp_global_certificate.py` | Three baseline global P/P price certificates, sixteen profiles and all four unilateral flips per profile; restricted planner |
| `reviews/stage11r_p1_open_region_certificate.py` | Existing open cost band and small joint parameter neighbourhood, including price continuation checks |
| `reviews/stage11r_p1b_empirical_implications_check.py` | Existing three conditional model-to-observable mappings |
| `reviews/stage11r_p1c_reduction_check.py` | Existing single-index reduction and selection/composition identities |
| `verification/stage33_reproduce.py` | Original numerical regression, all profiles, original fees/welfare/thresholds |
| `government_game_exploration/incidence_screening_sign_certificate.py` | Existing exact-rational continuous two-parameter sign certificate |
| `government_game_exploration/incidence_screening_valuation_check.py` | Existing positive conditional public value and negative venture increment certificate |
| `revision_p1d/diagnostic_check.py` | Existing exact threshold/decomposition/cost-width checks |
| `verification/beta22_adverse.py` | Existing negative ability-distribution portability result |
| `verification/referee_sensitivity.py` | Existing 396-row sensitivity panel and unchanged baseline restoration |

Eight regenerated output artifacts are byte-identical between baseline and revised-source runs: P/P certificate JSON, open-region JSON, both incidence JSONs, diagnostic JSON, Beta JSON, sensitivity CSV and sensitivity SVG. The new package's `verification_checks.json` preserves script exit codes and output SHA-256 hashes.

The baseline three outcomes remain `(B,P)/(B,B)`, `(B,B)/(B,B)`, `(B,P)/(B,B)` (decentralised equilibrium/restricted planner). Beta outcomes remain `(P,P)/(P,P)` at both reported complementarity values. These reruns preserve the original evidence classifications: interval certificates are not simulated samples; the conditional rectangle is not a four-party equilibrium rectangle; sensitivity/Beta numerical continuations are not newly certified global price theorems. The differentiated-Bertrand record was not independently reconstructed.

## Lean scope

The unchanged `formal/Stage33Formal.lean` does not encode Proposition 1's screening-incidence sign theorem. No Lean theorem statement therefore needed to be changed by this hypothesis correction. Neither `lean` nor `lake` is available locally, so **no local Lean rebuild is claimed**. Byte-identical formal artifacts preserve the previous certification scope; the existing limited formal statements do not become a formal proof of the entire game.

## PDF compilation and rendering

Final build and visual QA results are recorded in the new `revision_je_upgrade/` package and the final verification update below. The previous PR #73 PDFs are preserved as an archive. Current sources must build the main before the supplement so cross-document references read the new `main.aux`.


### Final compilation and visual QA

The final source revision was flattened into a scratch build directory. `latexmk -pdf -bibtex -interaction=nonstopmode -halt-on-error main.tex` succeeded, followed by the same command for `supplement.tex` against the new main auxiliary file. The new PDFs are `revision_je_upgrade/je_revised_manuscript.pdf` (**37 pages**) and `revision_je_upgrade/online_supplement.pdf` (**8 pages**). The abstract contains **197** whitespace-delimited words. Neither document has an undefined/multiply-defined-reference warning, overfull box or underfull box. `build_checks.json` records these results.

All **45** pages were rendered and inspected in contact sheets; the changed Proposition 1 page, the amended supplement proof-opening page and the original numerical table were also inspected individually. No clipping, overlap, broken equation rendering or unreadable table layout was observed. End-of-document tables follow the inherited `endfloat` configuration. The scope checker was rerun after final prose changes and passed. `git diff --check` passed. The previous PR #73 PDF/report archive was not overwritten.

The source checker reports Lean as “NOT EXECUTED BY THIS CHECK” so it remains truthful on other hosts; local toolchain unavailability is an observation of this verification session, rather than a property asserted by the portable script. A separately triggered CI Lean job, if added, is additional evidence and is not assumed successful in this report.
