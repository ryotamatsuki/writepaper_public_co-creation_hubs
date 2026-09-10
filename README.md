# Strategic Interaction among Public Innovation Hubs under Private Repricing

Current production status: **Stage 11 — INDEPENDENT ROBUSTNESS / REFEREE ATTACK GATE IN PROGRESS**.

Stage-10 manuscript authority: merge `02bcd37ac1ab019dc6df2c0d29c8f312a1c992a4`, verdict **FULL DRAFT READY FOR REFEREE GATE**.

Stage-9 reproducibility authority: merge `63c5f8c2627e97be6874438ca6af16e1df1c338a`, verdict **REPRODUCIBILITY BASELINE READY**.

Historical Stage-8 freeze authority is `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md`, merged at `ad927ca783a6123ea4fc6f55f65598ebd6ab583b`. Stage 11 found one theorem-quantifier certification regression: T3 relied on the symmetric beta-zero G state certified by T2 but the historical T3 prose did not state that branch anchor explicitly. The bounded Stage-7.5A repair is recorded in `reviews/STAGE_075A_V21_T3_SYMMETRY_LIMITED_REOPEN_2026-09-10.md`; the controlling T3 freeze amendment is `theory_freeze_v21/STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md`. All other frozen objects remain unchanged.

Stage 11 independent audit record: `reviews/STAGE_11_V21_REFEREE_GATE_2026-09-10.md`.

Reviewer-side independent regression: `stage11_v21_independent/code/independent_stage11_regression.py`. It does not import the Stage-4A production audit, result generator, or analytic derivation scripts.

Install `requirements.txt` and a LaTeX environment with `latexmk`, then run:

```sh
make clean && make all
python stage11_v21_independent/code/independent_stage11_regression.py
python scripts/stage75a_scope_audit.py
```

The Stage-11 gate attacks the analytic T1/T2 derivations, repaired T3 quantifier, candidate deviations, off-path continuations, sampled alternative equilibria, welfare accounting, B3 interpretation, remote-public friction, novelty/prior art, institutional mapping, and reader-facing scope. It does not treat a sampled alternative-equilibrium search as a uniqueness proof.

Generated quantitative objects live in `generated/` and must not be hand-edited. The rejected low-friction witness `(beta,gamma,tau)=(.05,.9,.05)`, its old 20-draw perturbation exercise, and the old exact stationary-root certificate are not active equilibrium authority.

Stage-10 exposition conclusion remains: **no production figure is required**. The analytic result is carried by theorem/proof and the repaired one-vector strategic and welfare results by generated tables.

Current Stage-11 provisional verdict: **CONDITIONAL GO — PR-HEAD CI REQUIRED**. Fatal issues found so far: **NONE**.

Journal target: **NOT SELECTED — deferred to Stage 12**.
