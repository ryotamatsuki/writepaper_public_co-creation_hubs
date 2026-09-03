from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
TEX_FILES = [ROOT/"paper"/"main.tex", *sorted((ROOT/"sections").glob("*.tex"))]
text = "\n".join(p.read_text(encoding="utf-8") for p in TEX_FILES)

# No internal workflow language in manuscript.
forbidden_internal = [
    r"\bStage 4\b", r"\bStage 6\b", r"\bStage 7(?:\.5)?\b", r"\bStage 8\b", r"\bStage 9\b",
    r"\bhard kill\b", r"\bre-kill\b", r"\bNO-GO\b", r"\bcanonical workflow\b"
]
for pat in forbidden_internal:
    assert not re.search(pat, text, re.I), f"internal workflow language in manuscript: {pat}"

# Citation keys used in LaTeX must exist in bibliography.
bib = (ROOT/"references"/"references.bib").read_text(encoding="utf-8")
bib_keys = set(re.findall(r"@\w+\{([^,]+),", bib))
cite_blocks = re.findall(r"\\cite[tp]?\{([^}]+)\}", text)
used = {k.strip() for block in cite_blocks for k in block.split(",") if k.strip()}
missing = used - bib_keys
assert not missing, f"missing bibliography keys: {sorted(missing)}"

# Headline theorem and proof-status discipline.
assert r"\BR_i^{B3\prime}>0>\BR_i^{G\prime}" in text, "headline sign reversal missing"
assert r"\begin{theorem}[Local strategic sign reversal]" in text, "analytic theorem missing"
assert r"\begin{proposition}[Exact certification at the canonical parameter vector]" in text, "exact canonical certificate proposition missing"
assert "global sign theorem" in text, "global-theorem limitation must remain explicit"
assert "not an empirical calibration" in text, "canonical vector must not be framed as empirical calibration"
assert "primitive" in text and "necessary-and-sufficient" in text, "primitive-classification limitation must remain explicit"

# Prevent regression to the superseded numerical-only proof ceiling.
legacy_proof_phrases = [
    "numerically certified existence result",
    "numerically certified existence statement",
    "computationally certified local strategic reversal",
]
for phrase in legacy_proof_phrases:
    assert phrase.lower() not in text.lower(), f"superseded numerical-only proof status found: {phrase}"

# Killed-claim safeguards in manuscript prose.
dangerous = [
    "strategic substitutes are socially desirable",
    "strategic complements are socially desirable",
    "optimal policy is to fix",
    "we prove that endogenous pricing makes public investments strategic substitutes",
    "we characterize the complete primitive parameter region",
]
for phrase in dangerous:
    assert phrase.lower() not in text.lower(), f"killed/overclaim language found: {phrase}"

print(f"PASS: manuscript audit ({len(used)} citation keys, {len(TEX_FILES)} tex files)")
