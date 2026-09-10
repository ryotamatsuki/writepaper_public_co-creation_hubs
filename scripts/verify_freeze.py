from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FREEZE = ROOT / "theory_freeze_v21" / "CANONICAL_THEORY_FREEZE_2026-09-10.md"
AMEND = ROOT / "theory_freeze_v21" / "STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md"
SCOPE = ROOT / "reviews" / "STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md"
REOPEN = ROOT / "reviews" / "STAGE_075A_V21_T3_SYMMETRY_LIMITED_REOPEN_2026-09-10.md"

for path, label in [
    (FREEZE, "canonical Stage-8 freeze"),
    (AMEND, "Stage-8 T3 amendment"),
    (SCOPE, "Stage-7.5A claim-scope ledger"),
    (REOPEN, "Stage-7.5A T3 limited-reopen record"),
]:
    assert path.exists(), f"missing {label}"

text = FREEZE.read_text(encoding="utf-8")
amend = AMEND.read_text(encoding="utf-8")
reopen = REOPEN.read_text(encoding="utf-8")

required = [
    "Workflow authority: `ryotamatsuki/research-paper-workflow` v2.1 @ `b27568172d4145a2b98825a5b66a071dbfe25f36`",
    "Certified input SHA: `eeb48a3dd76ab6f43d5de175b12c3f374746db0f`",
    "THEORY FROZEN — GO TO REPRODUCIBILITY SETUP",
    "beta=.01",
    "gamma=.825",
    "tau=.35",
    "0.8371022382025995",
    "0.8258903860237495",
    "0.018403679612460814",
    "local sufficient-condition theorem",
    "all-regime computational global-equilibrium existence witness",
    "old vector `(beta=.05, gamma=.9, tau=.05)` is rejected",
    "old `±0.5%`, 20-draw robustness around that rejected vector is non-authoritative",
    "No unresolved continuation is authorized for the headline repaired witness",
    "reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md",
]
for item in required:
    assert item in text, f"freeze drift/missing authority: {item}"

# Stage-11 certification regression: the historical freeze is preserved, while this
# amendment narrows T3 to the symmetric beta-zero G anchor already required by T2.
assert "supersedes **only the T3 G-branch anchor wording**" in amend
assert "symmetric regular beta-zero full-game stationary state" in amend
assert "CERTIFICATION REGRESSION" in reopen
assert "LIMITED STAGE 7.5A REOPEN CLOSED — QUANTIFIER REPAIRED" in reopen

for forbidden in [
    "arbitrary distribution robustness;\n- arbitrary nonlinear-network robustness;\n- heterogeneous-region equilibrium theorem;\n- global primitive-space reversal theorem;\n- uniqueness of the repaired equilibrium;",
]:
    # The listed items must occur only inside the explicit NOT-CLAIMED register.
    assert forbidden in text

print("PASS: v2.1 Stage-8 freeze + T3 symmetry amendment, repaired witness, and claim ceiling")
