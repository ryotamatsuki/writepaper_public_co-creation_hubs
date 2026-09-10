from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FREEZE = ROOT / "theory_freeze_v21" / "CANONICAL_THEORY_FREEZE_2026-09-10.md"
T3_AMEND = ROOT / "theory_freeze_v21" / "STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md"
ASTRA_AMEND = ROOT / "theory_freeze_v21" / "STAGE_08_V21_ASTRA_WELFARE_EVIDENCE_AMENDMENT_2026-09-10.md"
SCOPE = ROOT / "reviews" / "STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md"
T3_REOPEN = ROOT / "reviews" / "STAGE_075A_V21_T3_SYMMETRY_LIMITED_REOPEN_2026-09-10.md"
WELFARE_CORR = ROOT / "reviews" / "STAGE_07_V21_SATURATED_PARTNER_SURPLUS_CORRECTION_2026-09-10.md"
GLOBAL_REPAIR = ROOT / "reviews" / "STAGE_04A_V21_ASTRA_GLOBAL_EVIDENCE_REPAIR_2026-09-10.md"

for path, label in [
    (FREEZE, "historical canonical Stage-8 freeze"),
    (T3_AMEND, "Stage-8 T3 amendment"),
    (ASTRA_AMEND, "Stage-8 Astra welfare/evidence amendment"),
    (SCOPE, "historical Stage-7.5A claim-scope ledger"),
    (T3_REOPEN, "Stage-7.5A T3 limited-reopen record"),
    (WELFARE_CORR, "Stage-7 support-surplus correction"),
    (GLOBAL_REPAIR, "Stage-4A Astra evidence repair"),
]:
    assert path.exists(), f"missing {label}"

text = FREEZE.read_text(encoding="utf-8")
t3_amend = T3_AMEND.read_text(encoding="utf-8")
astra = ASTRA_AMEND.read_text(encoding="utf-8")
t3_reopen = T3_REOPEN.read_text(encoding="utf-8")
welfare = WELFARE_CORR.read_text(encoding="utf-8")
global_repair = GLOBAL_REPAIR.read_text(encoding="utf-8")

# Historical freeze is immutable evidence of what Stage 8 concluded at the time.
required_historical = [
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
]
for item in required_historical:
    assert item in text, f"historical freeze drift/missing authority: {item}"

# Later amendments are the controlling authority; they must not silently rewrite history.
assert "supersedes **only the T3 G-branch anchor wording**" in t3_amend
assert "symmetric regular beta-zero full-game stationary state" in t3_amend
assert "CERTIFICATION REGRESSION" in t3_reopen
assert "LIMITED STAGE 7.5A REOPEN CLOSED — QUANTIFIER REPAIRED" in t3_reopen

assert "supersedes only the all-regime welfare-accounting authority and the evidence qualification of T4/W2/W3" in astra
assert "ALL-REGIME COMPUTATIONAL SEARCH EVIDENCE — NO CERTIFIED GLOBAL REGRET BOUND" in astra
assert "must not be described as a certified global equilibrium" in astra
assert "Stage 12 is blocked until Astra limited recheck" in astra

assert "r_h*m_h - m_h^2/2" in welfare
assert "SEARCH EVIDENCE" in global_repair
assert "global-equilibrium existence witness/certification" in global_repair

print("PASS: historical Stage-8 freeze preserved; T3 and Astra welfare/evidence amendments are controlling")
