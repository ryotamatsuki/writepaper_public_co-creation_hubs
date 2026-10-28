# Anonymous review source and replication files

The two LaTeX files and bibliography are flat, with no subfolder dependencies.
No author identity or identifying repository link is included.
Build `main.tex` first, then `supplement.tex` in the same directory so that the
supplement can import the main equation labels using xr-hyper.

```bash
latexmk -pdf -bibtex -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -bibtex -interaction=nonstopmode -halt-on-error supplement.tex
```

Python 3.11+ standard library suffices for the mathematical checks:

```bash
python stage11r_pp_global_certificate.py
python stage11r_p1_open_region_certificate.py
python stage11r_p1b_empirical_implications_check.py
python stage11r_p1c_reduction_check.py
python stage33_reproduce.py
python diagnostic_check.py
```

The first two scripts certify the stated original benchmarks and small cost/
parameter regions. The original numerical reproduction is not a general
price-equilibrium proof. The diagnostic script checks the unchanged benchmark
thresholds and their service/matching/market decomposition. It does not certify
historical priority, all-lambda SPNE, or any new Lean theorem.

Funding, competing interests, author affiliation and corresponding-author
details belong on a separate title page and must be confirmed by the author
before submission. This anonymous package is a review draft, not a record of
journal submission or acceptance.
