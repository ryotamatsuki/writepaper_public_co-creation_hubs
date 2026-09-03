PYTHON ?= python
.PHONY: help freeze symbolic numerical results tables figures bibliography test manuscript verify report manifest all clean
help:
	@echo 'make verify|tables|figures|test|manuscript|all|clean'
freeze:
	$(PYTHON) scripts/verify_freeze.py
symbolic:
	$(PYTHON) scripts/verify_symbolic.py
numerical:
	$(PYTHON) scripts/verify_numerical.py
results:
	$(PYTHON) scripts/generate_results.py
tables: results
	$(PYTHON) scripts/generate_tables.py
figures:
	$(PYTHON) scripts/generate_figures.py
bibliography:
	$(PYTHON) scripts/validate_bibliography.py
test: results tables
	$(PYTHON) -m pytest -q tests
verify: freeze symbolic numerical bibliography
manuscript: tables figures
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
manifest: tables figures
	$(PYTHON) scripts/generate_manifest.py
report: verify test manuscript
	$(PYTHON) scripts/generate_verification_report.py
all: report manifest
clean:
	-rm -f generated/results/canonical_results.json generated/results/manifest.json generated/results/verification_report.json generated/tables/*.tex generated/figures/README.md
	-cd paper && latexmk -C main.tex >/dev/null 2>&1 || true
