import json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def test_generated_results_and_signs():
    subprocess.run([sys.executable,str(ROOT/'scripts/generate_results.py')],check=True); r=json.loads((ROOT/'generated/results/canonical_results.json').read_text()); assert r['computed']['G']['BR_slope']<0<r['computed']['B3']['BR_slope']; assert r['computed']['local_perturbation_successes']==20; assert r['journal_target']=='NOT SELECTED'
def test_generated_tables_match_json():
    subprocess.run([sys.executable,str(ROOT/'scripts/generate_tables.py')],check=True); r=json.loads((ROOT/'generated/results/canonical_results.json').read_text()); t=(ROOT/'generated/tables/strategic_results.tex').read_text(); assert f"{r['computed']['G']['x']:.6f}" in t; assert f"{r['computed']['G']['BR_slope']:.6f}" in t; assert 'DO NOT EDIT' in t
def test_freeze_gate(): subprocess.run([sys.executable,str(ROOT/'scripts/verify_freeze.py')],check=True)
