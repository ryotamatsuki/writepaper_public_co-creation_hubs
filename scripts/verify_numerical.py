from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "generated/results/canonical_results.json"


def close(got: float, expected: float, tol: float, label: str) -> None:
    if abs(float(got) - float(expected)) > tol:
        raise AssertionError(f"{label}: {got} != {expected} within {tol}")


assert RESULTS.exists(), "canonical result layer missing; run make results first"
r = json.loads(RESULTS.read_text(encoding="utf-8"))
p = r["parameters"]
c = r["computed"]

expected_parameters = {
    "v": 0.1,
    "alpha": 0.5,
    "beta": 0.01,
    "rho": 0.15,
    "rhoT": 0.05,
    "kL": 0.27,
    "kT": 0.02,
    "tau": 0.35,
    "gamma": 0.825,
}
for key, expected in expected_parameters.items():
    close(p[key], expected, 1e-12, f"parameter {key}")

assert r["provenance"]["stage8_merge_sha"] == "ad927ca783a6123ea4fc6f55f65598ebd6ab583b"
assert r["provenance"]["stage75a_certified_input_sha"] == "eeb48a3dd76ab6f43d5de175b12c3f374746db0f"
assert r["provenance"]["stage11r_repair"] == "reviews/STAGE_11R_ASTRA_REPAIR_REPORT.md"

# The reported stationary candidates are retained because the R1 correction is inactive
# on path (support participation is interior there).  We nevertheless verify stationarity
# numerically rather than treating the stored values as proof of optimality.
close(c["G"]["x"], 0.8371022382025995, 1e-10, "G x")
close(c["B3"]["x"], 0.8258903860237495, 1e-10, "B3 x")
close(c["G"]["p_T"], 0.018403679612460814, 2e-5, "G p_T")
assert c["G"]["BR_slope"] < 0 < c["B3"]["BR_slope"], "headline repaired slope signs failed"
close(c["G"]["BR_slope"], -0.019870, 5e-3, "G BR slope")
close(c["B3"]["BR_slope"], 0.003968, 5e-3, "B3 BR slope")

s = c["stationarity"]
assert s["participation_residual_inf_G"] < 1e-8
assert s["participation_residual_inf_B3"] < 1e-8
assert abs(s["private_price_foc"]) < 2e-5
assert s["private_price_soc"] < 0
assert abs(s["G_public_own_foc"]) < 2e-5
assert abs(s["B3_public_own_foc"]) < 2e-5

search = c["global_search"]
assert search["evidence_level"] == "SEARCH EVIDENCE"
assert search["public_interval"] == [0.0, 1.0]
assert search["G_deviation_private_reoptimization"] is True
assert search["B3_fixed_matched_price"] is True
assert search["certified_regret_upper_bound"] is None
assert search["G_best_detected_gain"] < 3e-5
assert search["B3_best_detected_gain"] < 3e-5

w = c["welfare"]
close(w["G"]["W_N"], 0.596788, 7e-4, "G aggregate welfare")
close(w["B3"]["W_N"], 0.585928, 7e-4, "B3 aggregate welfare")
assert w["G"]["W_N"] > w["B3"]["W_N"], "one-state-pair welfare ranking failed"
# On-path support masses are interior, so the corrected primitive integral agrees with
# the historical interior shortcut at the reported states.
assert all(0.0 < float(v) < 1.0 for v in w["G"]["support_gross"])
assert all(0.0 < float(v) < 1.0 for v in w["B3"]["support_gross"])

coord = c["coordination"]
close(coord["national_direct"], 0.49045, 7e-3, "direct local coordination derivative")
close(coord["sum_components"], 0.49045, 7e-3, "decomposed local coordination derivative")
assert abs(coord["own_welfare"]) < 2e-5
assert coord["rival_welfare"] > 0
assert coord["private_profit"] < 0
assert abs(coord["decomposition_error"]) < 1e-8

assert "SEARCH EVIDENCE" in r["proof_status"]["repaired_witness"]
assert "NO CERTIFIED GLOBAL REGRET BOUND" in r["proof_status"]["repaired_witness"]
assert "NOT GLOBAL EQUILIBRIUM AUTHORITY" in r["proof_status"]["old_exact_certificate"]
assert r["journal_target"] == "NOT SELECTED — STAGE 12 BLOCKED PENDING ASTRA RECHECK"
print("PASS: corrected welfare, stationarity, search-evidence scope, and direct welfare derivative checks")
