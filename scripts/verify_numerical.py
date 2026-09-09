from __future__ import annotations

from generate_results import build_results


def close(got: float, expected: float, tol: float, label: str) -> None:
    if abs(float(got) - float(expected)) > tol:
        raise AssertionError(f"{label}: {got} != {expected} within {tol}")


r = build_results()
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
print("PASS: repaired all-regime witness, slope signs, welfare witness, and coordination wedge")
