# Strategic Interaction among Public Innovation Hubs under Private Repricing

Current production status: **Stage 11R — ASTRA FINDINGS REPAIR IN PROGRESS / STAGE 12 BLOCKED**.

Astra audited Stage-10 baseline `02bcd37ac1ab019dc6df2c0d29c8f312a1c992a4` and required an earlier-stage bounded repair. The canonical non-Astra Stage-11 merge `7fdb1020fe398e5bd6833f9faa5da56812ad345e` is retained as historical pre-Astra audit evidence, not current final authority.

Historical Stage-8 freeze: `theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md`, merged at `ad927ca783a6123ea4fc6f55f65598ebd6ab583b`.

Controlling amendments:

- T3 symmetry/quantifier: `theory_freeze_v21/STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md` plus `reviews/STAGE_075A_V21_T3_SYMMETRY_LIMITED_REOPEN_2026-09-10.md`.
- Astra welfare/evidence repair: `theory_freeze_v21/STAGE_08_V21_ASTRA_WELFARE_EVIDENCE_AMENDMENT_2026-09-10.md`.

Stage 11R repair record: `reviews/STAGE_11R_ASTRA_REPAIR_REPORT.md`.

The repair corrects support-side welfare when participation saturates. With gross benefit `r_h` and mass `m_h=clip(r_h,0,1)`, national support surplus is `r_h*m_h-m_h^2/2`; the familiar `m_h^2/2` expression is valid only in the interior where `m_h=r_h`.

The existing vector `(beta,gamma,tau)=(.01,.825,.35)` is retained and rechecked; no new vector is searched. The all-regime numerical procedure covers the full public interval with grids/local refinement, private repricing after G deviations, multiple participation starts, boundaries, the low-investment region, and support-saturation histories. Its active evidence level is deliberately **SEARCH EVIDENCE**. It does **not** provide a certified global regret upper bound or verified global-equilibrium existence certificate.

T1 and T2 remain analytic local results. T3 remains a local stationary-branch sign-reversal theorem with the G branch anchored at the symmetric regular beta-zero central-interior state and with the stated support-interiority, routing, nonsingularity, and SOC conditions. B3 remains a matched-scalar-fee identification benchmark, not a planner or regulation counterfactual.

Install `requirements.txt` and a LaTeX environment with `latexmk`, then run:

```sh
make clean && make all
python scripts/stage75a_scope_audit.py
python stage11_v21_independent/code/independent_stage11_regression.py
python scripts/verify_manifest.py
```

Generated quantitative objects live in `generated/` and must not be hand-edited. Stage 11R verifies their hashes/sizes and clean-regeneration determinism. The rejected low-friction vector `(beta,gamma,tau)=(.05,.9,.05)`, its old 20-draw perturbation exercise, and the old exact stationary-root certificate remain historical negative/local evidence only.

Stage-10 exposition conclusion remains: **NO REQUIRED FIGURE**.

Journal target: **NOT SELECTED — Stage 12 remains blocked pending Astra limited recheck**.
