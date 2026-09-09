import json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def run(path):
    subprocess.run([sys.executable, str(ROOT/path)], check=True)

def test_repaired_global_certificate():
    run('stage4a_v21_repaired/code/independent_repaired_audit.py')

def test_generated_results_and_signs():
    run('scripts/generate_results.py')
    r=json.loads((ROOT/'generated/results/canonical_results.json').read_text())
    assert r['parameters']['beta']==0.01
    assert r['parameters']['gamma']==0.825
    assert r['parameters']['tau']==0.35
    assert r['computed']['G']['BR_slope']<0<r['computed']['B3']['BR_slope']
    assert abs(r['computed']['G']['x']-0.8371022382025995)<1e-9
    assert abs(r['computed']['B3']['x']-0.8258903860237495)<1e-9
    assert r['proof_status']['repaired_witness'].startswith('ALL-REGIME')

def test_generated_tables_match_json():
    run('scripts/generate_tables.py')
    r=json.loads((ROOT/'generated/results/canonical_results.json').read_text())
    t=(ROOT/'generated/tables/strategic_results.tex').read_text()
    assert f"{r['computed']['G']['x']:.6f}" in t
    assert f"{r['computed']['G']['BR_slope']:.6f}" in t
    assert 'DO NOT EDIT' in t

def test_scope_gate():
    run('scripts/stage75a_scope_audit.py')

def test_freeze_gate():
    run('scripts/verify_freeze.py')
