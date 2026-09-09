PYTHON ?= python
RESULT_JSON := generated/results/canonical_results.json
TABLE_STAMP := generated/tables/.stage9.stamp
FIGURE_REGISTRY := generated/figures/README.md
MANIFEST := generated/results/manifest.json

.PHONY: help freeze symbolic global results numerical scope tables figures bibliography manuscript-audit stage10-architecture test verify manuscript manifest report all clean

help:
	@echo 'make freeze|symbolic|global|numerical|scope|tables|figures|bibliography|manuscript-audit|stage10-architecture|test|manuscript|all|clean'

freeze:
	$(PYTHON) scripts/verify_freeze.py

symbolic:
	$(PYTHON) scripts/verify_symbolic.py

global:
	$(PYTHON) stage4a_v21_repaired/code/independent_repaired_audit.py

results: $(RESULT_JSON)

$(RESULT_JSON): scripts/generate_results.py stage4a_v21_repaired/code/independent_repaired_audit.py
	@mkdir -p generated/results
	$(PYTHON) scripts/generate_results.py

numerical: $(RESULT_JSON)
	$(PYTHON) scripts/verify_numerical.py

scope: $(RESULT_JSON)
	$(PYTHON) scripts/stage75a_scope_audit.py

tables: $(TABLE_STAMP)

$(TABLE_STAMP): $(RESULT_JSON) scripts/generate_tables.py
	@mkdir -p generated/tables
	$(PYTHON) scripts/generate_tables.py
	@touch $(TABLE_STAMP)

figures: $(FIGURE_REGISTRY)

$(FIGURE_REGISTRY): scripts/generate_figures.py
	@mkdir -p generated/figures
	$(PYTHON) scripts/generate_figures.py

bibliography:
	$(PYTHON) scripts/validate_bibliography.py

manuscript-audit: scope
	$(PYTHON) scripts/validate_manuscript.py

stage10-architecture:
	$(PYTHON) scripts/validate_stage10_architecture.py

test: freeze symbolic numerical scope tables stage10-architecture
	$(PYTHON) -m pytest -q tests

verify: freeze symbolic global numerical scope bibliography manuscript-audit stage10-architecture

manuscript: tables figures scope bibliography manuscript-audit stage10-architecture
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
	$(PYTHON) scripts/validate_build_log.py

manifest: $(MANIFEST)

$(MANIFEST): $(RESULT_JSON) $(TABLE_STAMP) $(FIGURE_REGISTRY) scripts/generate_manifest.py
	$(PYTHON) scripts/generate_manifest.py

report: verify test manuscript
	$(PYTHON) scripts/generate_verification_report.py

all: report manifest

clean:
	-rm -f generated/results/canonical_results.json generated/results/manifest.json generated/results/verification_report.json generated/tables/*.tex generated/tables/.stage9.stamp generated/figures/README.md
	-cd paper && latexmk -C main.tex >/dev/null 2>&1 || true
