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

# Final Stage-11R theorem / numerical-evidence hierarchy.
required = [
    r"\begin{proposition}[First-order benchmark complementarity]",
    r"\begin{theorem}[Local strategic sign reversal]",
    r"\BR_i^{B3\prime}>0>\BR_i^{G\prime}",
    "symmetric regular beta-zero central-interior full-game stationary state",
    "Repaired all-regime computational search",
    "search evidence",
    "Support-side surplus with participation caps",
    "one-state-pair numerical comparison only",
]
for needle in required:
    assert needle in text, f"required Stage-11R manuscript scope missing: {needle}"

# Accept equivalent reader-facing negations of a certified regret bound.
regret_bound_disclaimers = [
    "no rigorous upper bound on regret",
    "neither a rigorous upper bound on regret",
]
assert any(phrase in text for phrase in regret_bound_disclaimers), (
    "required Stage-11R regret-bound disclaimer missing"
)

# The active manuscript may discuss historical certificates negatively, but it must not
# positively describe the current repaired vector as a certified global equilibrium.
positive_overclaims = [
    "computationally certified all-regime witness",
    "computational global-equilibrium existence witness",
    "all-regime computational certification at one baseline parameter vector",
    "This repaired computation is an existence witness",
    "Global best-response status is addressed separately",
]
for phrase in positive_overclaims:
    assert phrase not in text, f"stale positive global-certification wording found: {phrase}"

# Correct support-side primitive accounting must be reader-visible.
assert "r_hm_h-\\frac{m_h^2}{2}" in text
assert "m_h=\\min\\{1,\\max\\{0,r_h\\}\\}" in text

# The repaired vector remains constructive/non-calibrated.
calibration_disclaimers = [
    "not empirically calibrated",
    "rather than empirically calibrated",
    "not an empirical calibration",
]
assert any(phrase in text for phrase in calibration_disclaimers), (
    "required repaired-vector calibration disclaimer missing"
)

# Permanent rejected-evidence discipline.
assert "Archived Exact Stationary-Root Diagnostic" in text
assert "not evidence that the old root is a global public Nash equilibrium or subgame-perfect equilibrium" in text
assert "earlier 20-draw perturbation exercise around the rejected vector is non-authoritative" in text

for phrase in [
    "we characterize the complete primitive parameter region",
    "we prove global strategic sign reversal",
    "we establish uniqueness of the repaired equilibrium",
]:
    assert phrase.lower() not in text.lower(), f"killed/overclaim language found: {phrase}"

print(f"PASS: Stage-11R manuscript audit ({len(used)} citation keys, {len(TEX_FILES)} tex files)")
