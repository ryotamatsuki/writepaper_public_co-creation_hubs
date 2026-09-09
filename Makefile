PYTHON ?= python
.PHONY: help freeze global scope results tables bibliography test manuscript report all clean
help:
	@echo 'make freeze|global|scope|tables|test|manuscript|all|clean'
freeze:
	$(PYTHON) scripts/verify_freeze.py
global:
	$(PYTHON) stage4a_v21_repaired/code/independent_repaired_audit.py
scope:
	$(PYTHON) scripts/stage75a_scope_audit.py
results:
	$(PYTHON) scripts/generate_results.py
tables: results
	$(PYTHON) scripts/generate_tables.py
bibliography:
	$(PYTHON) scripts/validate_bibliography.py
test: freeze results tables scope
	$(PYTHON) -m pytest -q tests
manuscript: tables scope bibliography
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
	$(PYTHON) scripts/validate_build_log.py
report: freeze global scope test manuscript
	@echo 'v2.1 reproducibility gate complete'
all: report
clean:
	-rm -f generated/results/canonical_results.json generated/results/manifest.json generated/results/verification_report.json generated/tables/*.tex generated/figures/README.md
	-cd paper && latexmk -C main.tex >/dev/null 2>&1 || true
