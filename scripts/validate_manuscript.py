from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
TEX_FILES = [ROOT / "paper" / "main.tex", *sorted((ROOT / "sections").glob("*.tex"))]
text = "\n".join(p.read_text(encoding="utf-8") for p in TEX_FILES)

# Keep internal workflow / decision jargon out of reader-facing manuscript prose.
forbidden_internal = [r"\bhard kill\b", r"\bre-kill\b", r"\bNO-GO\b", r"\bcanonical workflow\b"]
for pat in forbidden_internal:
    assert not re.search(pat, text, re.I), f"internal workflow language in manuscript: {pat}"

# Citation keys used in LaTeX must exist in bibliography.
bib = (ROOT / "references" / "references.bib").read_text(encoding="utf-8")
bib_keys = set(re.findall(r"@\w+\{([^,]+),", bib))
cite_blocks = re.findall(r"\\cite[tp]?\{([^}]+)\}", text)
used = {k.strip() for block in cite_blocks for k in block.split(",") if k.strip()}
missing = used - bib_keys
assert not missing, f"missing bibliography keys: {sorted(missing)}"

# Final v2.1 theorem / witness hierarchy must be visible.
required = [
    r"\begin{proposition}[First-order benchmark complementarity]",
    r"\begin{theorem}[Local strategic sign reversal]",
    r"\BR_i^{B3\prime}>0>\BR_i^{G\prime}",
    "regular interior stationary branches",
    "Repaired all-regime computational witness",
    "computationally certified all-regime witness",
    "not an exact interval proof of global optimality for all primitives",
    "not a solved first-best or global social-optimum comparison",
]
for needle in required:
    assert needle in text, f"required final-v2.1 manuscript scope missing: {needle}"

# The repaired vector must remain explicitly constructive/non-calibrated, but do not
# couple the gate to one exact English phrasing.
calibration_disclaimers = [
    "not empirically calibrated",
    "rather than empirically calibrated",
    "not an empirical calibration",
]
assert any(phrase in text for phrase in calibration_disclaimers), (
    "required repaired-witness calibration disclaimer missing"
)

# Permanent rejected-evidence discipline.
assert "Archived Exact Stationary-Root Diagnostic" in text
assert "not evidence that the old root is a global public Nash equilibrium or subgame-perfect equilibrium" in text
assert "earlier 20-draw perturbation exercise around the rejected vector is non-authoritative" in text

# Prevent unambiguously positive regressions to claims killed by the freeze.
for phrase in [
    "we characterize the complete primitive parameter region",
    "we prove global strategic sign reversal",
    "we establish uniqueness of the repaired equilibrium",
]:
    assert phrase.lower() not in text.lower(), f"killed/overclaim language found: {phrase}"

print(f"PASS: final-v2.1 manuscript audit ({len(used)} citation keys, {len(TEX_FILES)} tex files)")
