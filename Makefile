PYTHON ?= python
.PHONY: help freeze symbolic global results numerical scope tables figures bibliography manuscript-audit test verify manuscript manifest report all clean

help:
	@echo 'make freeze|symbolic|global|numerical|scope|tables|figures|test|manuscript|all|clean'

freeze:
	$(PYTHON) scripts/verify_freeze.py

symbolic:
	$(PYTHON) scripts/verify_symbolic.py

global:
	$(PYTHON) stage4a_v21_repaired/code/independent_repaired_audit.py

results:
	$(PYTHON) scripts/generate_results.py

numerical: results
	$(PYTHON) scripts/verify_numerical.py

scope: results
	$(PYTHON) scripts/stage75a_scope_audit.py

tables: results
	$(PYTHON) scripts/generate_tables.py

figures:
	$(PYTHON) scripts/generate_figures.py

bibliography:
	$(PYTHON) scripts/validate_bibliography.py

manuscript-audit: results scope
	$(PYTHON) scripts/validate_manuscript.py

test: freeze symbolic numerical scope tables
	$(PYTHON) -m pytest -q tests

verify: freeze symbolic global numerical scope bibliography manuscript-audit

manuscript: tables figures scope bibliography manuscript-audit
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
	$(PYTHON) scripts/validate_build_log.py

manifest: results tables figures
	$(PYTHON) scripts/generate_manifest.py

report: verify test manuscript
	$(PYTHON) scripts/generate_verification_report.py

all: report manifest

clean:
	-rm -f generated/results/canonical_results.json generated/results/manifest.json generated/results/verification_report.json generated/tables/*.tex generated/figures/README.md
	-cd paper && latexmk -C main.tex >/dev/null 2>&1 || true
