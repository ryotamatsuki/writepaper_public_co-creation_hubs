from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
r = json.loads((ROOT / 'generated/results/canonical_results.json').read_text())
out = ROOT / 'generated/tables'
out.mkdir(parents=True, exist_ok=True)
ROW_END = " \\\\"


def write(name, body):
    (out / name).write_text('% DO NOT EDIT — GENERATED FILE\n' + body, encoding='utf-8')


p = r['parameters']
c = r['computed']
w = c['welfare']
rows = '\n'.join([f"{k} & {v:.6g}{ROW_END}" for k, v in p.items()])
write(
    'canonical_parameters.tex',
    "\\begin{tabular}{lr}\\toprule Parameter & Value \\\\ \\midrule\n"
    + rows
    + "\n\\bottomrule\\end{tabular}\n",
)
write(
    'strategic_results.tex',
    f"\\begin{{tabular}}{{lrr}}\\toprule Environment & reported $x$ & local BR slope \\\\ \\midrule\n"
    f"G & {c['G']['x']:.6f} & {c['G']['BR_slope']:.6f}{ROW_END}\n"
    f"B3 & {c['B3']['x']:.6f} & {c['B3']['BR_slope']:.6f}{ROW_END}\n"
    "\\bottomrule\\end{tabular}\n",
)
write(
    'welfare_comparison.tex',
    f"\\begin{{tabular}}{{lrrr}}\\toprule Environment & $W_i$ & $\\Pi_T$ & $W^N$ \\\\ \\midrule\n"
    f"G & {w['G']['W_i']:.6f} & {w['G']['Pi_T']:.6f} & {w['G']['W_N']:.6f}{ROW_END}\n"
    f"B3 & {w['B3']['W_i']:.6f} & {w['B3']['Pi_T']:.6f} & {w['B3']['W_N']:.6f}{ROW_END}\n"
    "\\bottomrule\\end{tabular}\n",
)
write(
    'proof_status.tex',
    "\\begin{tabular}{ll}\\toprule Claim & Status \\\\ \\midrule\n"
    f"Analytic reversal & LOCAL SUFFICIENT-CONDITION THEOREM{ROW_END}\n"
    f"Repaired numerical vector & ALL-REGIME SEARCH EVIDENCE{ROW_END}\n"
    f"Certified global regret bound & NOT AVAILABLE{ROW_END}\n"
    f"Broad robustness & NOT CLAIMED{ROW_END}\n"
    "\\bottomrule\\end{tabular}\n",
)
print('generated', len(list(out.glob('*.tex'))), 'tables')
