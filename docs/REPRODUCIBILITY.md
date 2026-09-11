# Reproducibility

The current v2.1 repair line is governed by the historical Stage-8 freeze together with its explicit Stage-11 and Stage-11R amendments. From a clone of the repair branch or its eventual merged authority, install the pinned Python requirements and a LaTeX environment with `latexmk`, then run:

```sh
python -m venv .venv
# activate the environment for your shell
python -m pip install -r requirements.txt
make clean && make all
python scripts/stage75a_scope_audit.py
python stage11_v21_independent/code/independent_stage11_regression.py
python scripts/verify_manifest.py
```

`make all` performs the following chain without adding a new theorem or parameter vector:

1. verifies the historical Stage-8 freeze plus the controlling T3 and Astra welfare/evidence amendments;
2. verifies the analytic small-beta identities and exact fee-transfer identity;
3. runs the repaired all-regime Stage-4A **search-evidence** audit using primitive-correct support surplus under saturation;
4. regenerates the canonical numerical result layer;
5. checks the G/B3 local slope signs, reported-state welfare comparison, and local coordination derivative decomposition;
6. runs the Stage-7.5A quantifier/scope regression against active manuscript text;
7. validates bibliography/citations and Stage-10 exposition architecture;
8. runs deterministic pipeline tests including support-surplus boundary checks;
9. regenerates four LaTeX tables and the figure registry;
10. builds `paper/main.pdf` and validates the build log;
11. writes a verification report and SHA-256 manifest for generated objects;
12. recomputes and verifies every manifest byte count and hash.

Stage 11R additionally regenerates the result layer a second time from clean inputs and compares manifests for determinism. It also requires committed generated files to match the regenerated files exactly.

The numerical all-regime routine searches the full public interval with grids, local refinement, multiple participation starts, endpoint/low-investment/saturation checks, and private repricing after G deviations. These checks are **SEARCH EVIDENCE**. They do not provide a certified supremum bound on unsearched deviation gains. Accordingly the active result layer must state `NO CERTIFIED GLOBAL REGRET BOUND`; build success must never be reported as proof of global equilibrium existence.

No journal template, secret, external dataset, or hidden manual numerical step is required. Stage 12 remains blocked pending Astra limited recheck.
