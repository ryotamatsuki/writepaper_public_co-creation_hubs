import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(path: str) -> None:
    subprocess.run([sys.executable, str(ROOT / path)], check=True)


def results() -> dict:
    return json.loads((ROOT / "generated/results/canonical_results.json").read_text(encoding="utf-8"))


def test_final_freeze_and_scope_authorities_exist():
    assert (ROOT / "theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md").exists()
    assert (ROOT / "reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md").exists()
    run("scripts/verify_freeze.py")


def test_repaired_global_certificate_is_canonical_target():
    p = ROOT / "stage4a_v21_repaired/code/independent_repaired_audit.py"
    assert p.exists()
    text = p.read_text(encoding="utf-8")
    assert "STAGE4A_REPAIRED_GLOBAL_CERTIFICATION: PASS" in text
    assert "UNRESOLVED" in text and "MULTIPLE_EQUILIBRIA" in text
    assert "beta=.01, gamma=.825, tau=.35" in text


def test_generated_results_and_signs():
    r = results()
    assert r["provenance"]["stage8_merge_sha"] == "ad927ca783a6123ea4fc6f55f65598ebd6ab583b"
    assert r["parameters"]["beta"] == 0.01
    assert r["parameters"]["gamma"] == 0.825
    assert r["parameters"]["tau"] == 0.35
    assert r["computed"]["G"]["BR_slope"] < 0 < r["computed"]["B3"]["BR_slope"]
    assert abs(r["computed"]["G"]["x"] - 0.8371022382025995) < 1e-10
    assert abs(r["computed"]["B3"]["x"] - 0.8258903860237495) < 1e-10
    assert r["proof_status"]["repaired_witness"].startswith("ALL-REGIME")
    assert r["journal_target"] == "NOT SELECTED — DEFERRED TO STAGE 12"


def test_generated_tables_match_json():
    r = results()
    table = (ROOT / "generated/tables/strategic_results.tex").read_text(encoding="utf-8")
    assert f"{r['computed']['G']['x']:.6f}" in table
    assert f"{r['computed']['B3']['x']:.6f}" in table
    assert f"{r['computed']['G']['BR_slope']:.6f}" in table
    assert f"{r['computed']['B3']['BR_slope']:.6f}" in table
    assert "DO NOT EDIT" in table


def test_scope_gate():
    run("scripts/stage75a_scope_audit.py")


def test_old_rejected_witness_is_not_active_generated_authority():
    r = results()
    payload = json.dumps(r, sort_keys=True)
    assert '"beta": 0.05' not in payload
    assert '"gamma": 0.9' not in payload
    assert '"tau": 0.05' not in payload
    assert "20/20" not in payload
