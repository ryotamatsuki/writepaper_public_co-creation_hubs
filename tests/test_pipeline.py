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
    assert (ROOT / "theory_freeze_v21/STAGE_08_V21_ASTRA_WELFARE_EVIDENCE_AMENDMENT_2026-09-10.md").exists()
    assert (ROOT / "reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md").exists()
    run("scripts/verify_freeze.py")


def test_repaired_global_search_is_canonical_target():
    p = ROOT / "stage4a_v21_repaired/code/independent_repaired_audit.py"
    assert p.exists()
    text = p.read_text(encoding="utf-8")
    assert "SEARCH EVIDENCE" in text
    assert "PARTNER_SURPLUS_INTERIOR_SATURATION_BOUNDARIES: PASS" in text
    assert "certified regret bound" in text
    assert "beta=.01, gamma=.825, tau=.35" in text


def test_partner_surplus_primitive_formula():
    from stage4a_v21_repaired.code.independent_repaired_audit import partner_surplus_scalar

    assert partner_surplus_scalar(-0.2) == 0.0
    assert abs(partner_surplus_scalar(0.4) - 0.08) < 1e-14
    assert abs(partner_surplus_scalar(1.0) - 0.5) < 1e-14
    assert abs(partner_surplus_scalar(1.2) - 0.7) < 1e-14


def test_generated_results_and_signs():
    r = results()
    assert r["parameters"]["beta"] == 0.01
    assert r["parameters"]["gamma"] == 0.825
    assert r["parameters"]["tau"] == 0.35
    assert r["computed"]["G"]["BR_slope"] < 0 < r["computed"]["B3"]["BR_slope"]
    assert abs(r["computed"]["G"]["x"] - 0.8371022382025995) < 1e-10
    assert abs(r["computed"]["B3"]["x"] - 0.8258903860237495) < 1e-10
    assert r["proof_status"]["repaired_witness"] == "ALL-REGIME COMPUTATIONAL SEARCH EVIDENCE — NO CERTIFIED GLOBAL REGRET BOUND"
    assert r["computed"]["global_search"]["certified_regret_upper_bound"] is None
    assert r["journal_target"] == "NOT SELECTED — STAGE 12 BLOCKED PENDING ASTRA RECHECK"


def test_welfare_derivative_generation_is_direct():
    r = results()["computed"]["coordination"]
    assert "own_welfare" in r
    assert "rival_welfare" in r
    assert "private_profit" in r
    assert "national_direct" in r
    assert "decomposition_sum" in r
    assert "decomposition_error" in r
    assert abs(r["decomposition_error"]) < 1e-5
    assert abs(r["own_welfare"]) < 1e-3
    assert r["rival_welfare"] > 0
    assert r["private_profit"] < 0
    assert r["national_direct"] > 0


def test_generated_tables_match_json():
    r = results()
    table = (ROOT / "generated/tables/strategic_results.tex").read_text(encoding="utf-8")
    assert f"{r['computed']['G']['x']:.6f}" in table
    assert f"{r['computed']['B3']['x']:.6f}" in table
    assert f"{r['computed']['G']['BR_slope']:.6f}" in table
    assert f"{r['computed']['B3']['BR_slope']:.6f}" in table
    assert "DO NOT EDIT" in table
    proof = (ROOT / "generated/tables/proof_status.tex").read_text(encoding="utf-8")
    assert "SEARCH EVIDENCE" in proof
    assert "COMPUTATIONAL CERTIFICATION" not in proof


def test_scope_gate():
    run("scripts/stage75a_scope_audit.py")


def test_manifest_integrity():
    run("scripts/verify_manifest.py")


def test_old_rejected_witness_is_not_active_generated_authority():
    r = results()
    payload = json.dumps(r, sort_keys=True)
    assert '"beta": 0.05' not in payload
    assert '"gamma": 0.9' not in payload
    assert '"tau": 0.05' not in payload
    assert "20/20" not in payload
