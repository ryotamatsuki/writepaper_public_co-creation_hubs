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

close(c["G"]["x"], 0.8371022382025995, 1e-10, "G x")
close(c["B3"]["x"], 0.8258903860237495, 1e-10, "B3 x")
close(c["G"]["p_T"], 0.018403679612460814, 2e-5, "G p_T")
assert c["G"]["BR_slope"] < 0 < c["B3"]["BR_slope"], "headline repaired slope signs failed"
close(c["G"]["BR_slope"], -0.019870, 5e-3, "G BR slope")
close(c["B3"]["BR_slope"], 0.003968, 5e-3, "B3 BR slope")

w = c["welfare"]
close(w["G"]["W_N"], 0.596788, 7e-4, "G aggregate welfare")
close(w["B3"]["W_N"], 0.585928, 7e-4, "B3 aggregate welfare")
assert w["G"]["W_N"] > w["B3"]["W_N"], "one-witness welfare ranking failed"
close(c["coordination"]["national_wedge"], 0.49045, 7e-3, "local coordination wedge")

assert r["proof_status"]["repaired_witness"].startswith("ALL-REGIME")
assert "NOT GLOBAL EQUILIBRIUM AUTHORITY" in r["proof_status"]["old_exact_certificate"]
assert r["journal_target"] == "NOT SELECTED — DEFERRED TO STAGE 12"
print("PASS: generated repaired witness matches final freeze, slope signs, welfare witness, and coordination wedge")
