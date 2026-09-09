# Reproducibility

From a clone of the canonical Stage-9 branch:

```sh
python -m venv .venv
# activate the environment for your shell
python -m pip install -r requirements.txt
make clean && make all
```

`make all` performs the following chain without changing theory:

1. verifies the final Stage-8 freeze and Stage-7.5A scope authority;
2. verifies the analytic small-beta identities and exact fee-transfer identity;
3. runs the independent repaired all-regime Stage-4A global-equilibrium audit;
4. regenerates the repaired canonical numerical result layer;
5. checks the repaired G/B3 slope signs, welfare witness, and local coordination wedge;
6. runs the Stage-7.5A quantifier/scope regression against active manuscript text;
7. validates bibliography/citations;
8. runs deterministic pipeline tests;
9. regenerates four LaTeX tables and the figure registry;
10. builds `paper/main.pdf` and validates the build log;
11. writes a verification report and SHA-256 manifest for generated objects.

The GitHub Actions Stage-9 workflow additionally checks that the branch descends from final Stage-8 merge `ad927ca783a6123ea4fc6f55f65598ebd6ab583b` and compares generated manifests across a clean second regeneration.

No journal template, secret, external dataset, or hidden manual numerical step is required.
