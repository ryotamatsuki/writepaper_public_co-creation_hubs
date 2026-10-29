# JE contribution-closure revision

The authoritative revised sources remain `../manuscript/main.tex`, its active section inputs, `../manuscript/supplement.tex`, the appendix with the explicit baseline-success hypothesis annotation, and `../manuscript/references.bib`.

- `je_revised_manuscript.pdf`: complete revised English article (37 pages).
- `online_supplement.pdf`: rebuilt existing proof derivations with the matching positive-baseline domain with current main-text cross-references (8 pages).
- `build_checks.json`: final LaTeX/reference/layout warning checks.
- `verification_checks.json`: documented hypothesis-delta and scientific-content comparison, ten existing script runs and eight byte-identical baseline/revision outputs.
- Detailed independent verification: `../reviews/JE_UPGRADE_VERIFICATION_2026-10-29.md`.

These artifacts are for review. Proposition 1 now explicitly requires the positive baseline success already used in its proof. The mathematical scope correction and unchanged numerical evidence are checked separately from literature novelty and journal fit. Local Lean rebuild was not available. The earlier PR #73 archive remains unchanged.

## Rebuild PDFs

Run from `hotelling_entry_experimentation/production` with Python 3, `latexmk`, a TeX installation including the packages in the manuscript preambles, and BibTeX. The following creates flat sources in a separate build directory; it does not overwrite the legacy submission package.

```bash
python3 - <<'PY'
from pathlib import Path
import re, shutil
root=Path.cwd()
out=root/'revision_je_upgrade'/'build'
out.mkdir(exist_ok=True)
def flatten(path):
    return re.sub(r'\\input\{([^}]+)\}',
                  lambda m: flatten((root/m[1]).with_suffix('.tex')),
                  path.read_text())
for name in ('main','supplement'):
    text=flatten(root/'manuscript'/f'{name}.tex')
    (out/f'{name}.tex').write_text(text.replace(
        r'\bibliography{manuscript/references}',r'\bibliography{references}'))
shutil.copy2(root/'manuscript/references.bib',out/'references.bib')
PY
cd revision_je_upgrade/build
latexmk -pdf -bibtex -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -bibtex -interaction=nonstopmode -halt-on-error supplement.tex
```

The main must be built before the supplement so `xr-hyper` reads the current `main.aux`. Existing verification script paths and their evidence limits are listed in the detailed report. Generated verification outputs should be produced in an isolated copy to preserve frozen output files.
