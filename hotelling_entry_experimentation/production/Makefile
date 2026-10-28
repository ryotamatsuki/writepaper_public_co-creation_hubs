PYTHON ?= python3
LATEXMK ?= latexmk
.PHONY: check data paper clean
check:
	$(PYTHON) verification/stage33_reproduce.py
data:
	$(PYTHON) verification/stage33_reproduce.py --csv output/stage33_profiles.csv
paper:
	mkdir -p build
	$(LATEXMK) -pdf -bibtex -interaction=nonstopmode -halt-on-error -outdir=build manuscript/main.tex
clean:
	rm -rf build output
