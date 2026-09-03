from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]; r=json.loads((ROOT/'generated/results/canonical_results.json').read_text()); out=ROOT/'generated/tables'; out.mkdir(parents=True,exist_ok=True)
def write(name,body): (out/name).write_text('% DO NOT EDIT — GENERATED FILE\n'+body,encoding='utf-8')
p=r['parameters']; f=r['frozen_display']; c=r['computed']; w=c['welfare']; rows='\n'.join([f"{k} & {v:.6g} \\\\" for k,v in p.items()])
write('canonical_parameters.tex',"\\begin{tabular}{lr}\\toprule Parameter & Value \\\\ \\midrule\n"+rows+"\n\\bottomrule\\end{tabular}\n")
write('strategic_results.tex',f"\\begin{{tabular}}{{lrr}}\\toprule Regime & $x$ & BR slope \\\\ \\midrule\nG & {c['G']['x']:.6f} & {c['G']['BR_slope']:.6f} \\\\ \nB3 & {f['B3_strategic_x']:.6f} & {c['B3']['BR_slope']:.6f} \\\\ \n\\bottomrule\\end{{tabular}}\n")
write('welfare_comparison.tex',f"\\begin{{tabular}}{{lrrr}}\\toprule Regime & $W_i$ & $\\Pi_T$ & $W^N$ \\\\ \\midrule\nG & {w['G']['W_i']:.6f} & {w['G']['Pi_T']:.6f} & {w['G']['W_N']:.6f} \\\\ \nB3 & {w['B3']['W_i']:.6f} & {w['B3']['Pi_T']:.6f} & {w['B3']['W_N']:.6f} \\\\ \n\\bottomrule\\end{{tabular}}\n")
write('proof_status.tex',"\\begin{tabular}{ll}\\toprule Claim & Status \\\\ \\midrule\nHeadline reversal & NUMERICALLY SUPPORTED ONLY \\\\ \nRobustness & LOCAL/CONDITIONAL \\\\ \n\\bottomrule\\end{tabular}\n"); print('generated',len(list(out.glob('*.tex'))),'tables')
