# Environment

Python dependencies are pinned in `requirements.txt`:

- NumPy 2.3.5
- SciPy 1.17.0
- SymPy 1.14.0
- pytest 9.0.2

Canonical CI uses Python 3.12 on `ubuntu-latest` and installs `latexmk`, `texlive-latex-base`, `texlive-latex-recommended`, and `texlive-latex-extra` from the runner package repository.

The earlier local-equivalent environment used Python 3.13.5, GNU Make 4.4.1, latexmk 4.86, and pdfTeX/TeX Live 2025/dev. Stage 9 does not rely on unpinned Python packages, private services, secrets, external datasets, or machine-specific absolute paths.
