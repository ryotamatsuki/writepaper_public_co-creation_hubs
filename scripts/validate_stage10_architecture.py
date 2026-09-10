from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCH = (ROOT / "docs" / "STAGE_10_EXPOSITION_ARCHITECTURE.md").read_text(encoding="utf-8")
SECTION_MAP = (ROOT / "docs" / "STAGE_10_SECTION_MAP.md").read_text(encoding="utf-8")
MAIN = (ROOT / "paper" / "main.tex").read_text(encoding="utf-8")
RESULTS = (ROOT / "sections" / "03_main_results.tex").read_text(encoding="utf-8")
WELFARE = (ROOT / "sections" / "04_welfare.tex").read_text(encoding="utf-8")
APP = (ROOT / "sections" / "appendices.tex").read_text(encoding="utf-8")

required_architecture = [
    "T1 — first-order B3 complementarity",
    "T2 — private repricing channel",
    "T3 — local strategic sign reversal",
    "T4 — repaired all-regime numerical search",
    "T5 — remote-public route dominance at repaired vector",
    "support-side welfare cap",
    "W1 — aggregate fee transfer cancellation",
    "W2 — local coordination wedge",
    "W3 — repaired-state G/B3 welfare ranking",
    "NO REQUIRED FIGURE",
    "tab:strategic",
    "tab:welfare",
    "tab:parameters",
]
for needle in required_architecture:
    assert needle in ARCH, f"Stage-10/11R exposition architecture missing: {needle}"

required_sections = [
    "sections/01_model.tex",
    "sections/02_equilibrium.tex",
    "sections/03_main_results.tex",
    "sections/04_welfare.tex",
    "sections/05_robustness.tex",
    "sections/06_institutional_empirical.tex",
    "sections/07_related_literature.tex",
    "sections/08_introduction.tex",
    "sections/09_discussion.tex",
    "sections/10_conclusion.tex",
    "sections/appendices.tex",
]
for path in required_sections:
    assert path in SECTION_MAP, f"Stage-10 section map missing: {path}"

assert "\\label{tab:strategic}" in RESULTS
assert "../generated/tables/strategic_results.tex" in RESULTS
assert "\\label{tab:welfare}" in WELFARE
assert "../generated/tables/welfare_comparison.tex" in WELFARE
assert "\\label{tab:parameters}" in APP
assert "../generated/tables/canonical_parameters.tex" in APP

assert "\\begin{proposition}[First-order benchmark complementarity]" in RESULTS
assert "\\begin{proposition}[Private-price feedback]" in RESULTS
assert "\\begin{theorem}[Local strategic sign reversal]" in RESULTS
assert "Repaired all-regime computational search" in RESULTS
assert "search evidence" in RESULTS
assert "certified global-equilibrium" in RESULTS  # appears only in an explicit negation
assert "NO REQUIRED FIGURE" in ARCH
assert "one-vector search evidence" in ARCH

# No production figure may appear without a formal architecture reopening.
tex = "\n".join(
    p.read_text(encoding="utf-8")
    for p in [ROOT / "paper" / "main.tex", *sorted((ROOT / "sections").glob("*.tex"))]
)
assert "\\begin{figure}" not in tex, "production figure added without architecture reauthorization"

inputs = [
    "../sections/08_introduction",
    "../sections/01_model",
    "../sections/02_equilibrium",
    "../sections/03_main_results",
    "../sections/04_welfare",
    "../sections/05_robustness",
    "../sections/06_institutional_empirical",
    "../sections/07_related_literature",
    "../sections/09_discussion",
    "../sections/10_conclusion",
]
positions = [MAIN.index(f"\\input{{{x}}}") for x in inputs]
assert positions == sorted(positions), "reader-facing section order drifted from construction order"

print("PASS: Stage-10/11R section map and search-evidence exposition architecture")
