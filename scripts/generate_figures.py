from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; p=ROOT/'generated/figures/README.md'; p.parent.mkdir(parents=True,exist_ok=True)
p.write_text('# Figures\n\nNo frozen Stage-8 result requires a production figure at the reproducibility-baseline stage.\n',encoding='utf-8'); print(p)
