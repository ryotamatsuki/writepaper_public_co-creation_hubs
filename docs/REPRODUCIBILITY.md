# Reproducibility
From a clone of the Stage-9 branch: `python -m venv .venv`; activate it; `python -m pip install -r requirements.txt`; run `make clean && make all`. `make all` regenerates canonical JSON/tables, runs freeze/symbolic/numerical/bibliography/tests gates, writes the figure registry, and compiles `paper/main.pdf`. No journal template or secret is required.
