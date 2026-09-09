from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

ACTIVE = [
    ROOT / "paper/main.tex",
    ROOT / "sections/03_main_results.tex",
    ROOT / "sections/04_welfare.tex",
    ROOT / "sections/05_robustness.tex",
    ROOT / "sections/08_introduction.tex",
    ROOT / "sections/09_discussion.tex",
    ROOT / "sections/10_conclusion.tex",
]


def text(path: Path) -> str:
    if not path.exists():
        raise RuntimeError(f"missing required file: {path.relative_to(ROOT)}")
    return path.read_text(encoding="utf-8")


def require(haystack: str, needle: str, where: str) -> None:
    if needle not in haystack:
        raise RuntimeError(f"missing required scope phrase in {where}: {needle}")


def forbid(haystack: str, needle: str, where: str) -> None:
    if needle in haystack:
        raise RuntimeError(f"stale/overstated scope token in {where}: {needle}")


def main() -> None:
    active = {p.relative_to(ROOT).as_posix(): text(p) for p in ACTIVE}
    joined = "\n".join(active.values())

    # Rejected draft witness must not remain as active equilibrium authority.
    for token in [
        "0.684028",
        "0.656020",
        "0.1144",
        "-0.0241",
        "prop:exactcertificate",
        "exact canonical certificate",
        "20-draw perturbation",
    ]:
        forbid(joined, token, "active manuscript")

    abstract = active["paper/main.tex"]
    require(abstract, "regular interior stationary branches", "paper/main.tex")
    require(abstract, "all-regime computational audit", "paper/main.tex")
    require(abstract, "repaired global-equilibrium result is a computational witness rather than a general theorem", "paper/main.tex")

    main_results = active["sections/03_main_results.tex"]
    require(main_results, "Suppose there exist regular interior B3 and G stationary branches", "sections/03_main_results.tex")
    require(main_results, "does not by itself prove global public optimality", "sections/03_main_results.tex")
    require(main_results, "Repaired all-regime computational witness", "sections/03_main_results.tex")
    require(main_results, "not an exact interval proof of global optimality for all primitives", "sections/03_main_results.tex")

    welfare = active["sections/04_welfare.tex"]
    require(welfare, "This is not a solved first-best or global social-optimum comparison", "sections/04_welfare.tex")
    require(welfare, "one-witness numerical comparison only", "sections/04_welfare.tex")

    robust = active["sections/05_robustness.tex"]
    require(robust, "does not contain a certified theorem for arbitrary project-type distributions", "sections/05_robustness.tex")
    require(robust, "non-authoritative", "sections/05_robustness.tex")
    require(robust, "No broader genericity claim is made", "sections/05_robustness.tex")

    discussion = active["sections/09_discussion.tex"]
    require(discussion, "not global-equilibrium authority", "sections/09_discussion.tex")
    require(discussion, "not evidence of direct user poaching between public hubs", "sections/09_discussion.tex")
    require(discussion, "not a solved first-best problem", "sections/09_discussion.tex")

    conclusion = active["sections/10_conclusion.tex"]
    require(conclusion, "does not itself establish global public optimality", "sections/10_conclusion.tex")
    require(conclusion, "not a general theorem", "sections/10_conclusion.tex")

    appendix = text(ROOT / "sections/appendices.tex")
    require(appendix, "Archived Exact Stationary-Root Diagnostic", "sections/appendices.tex")
    require(appendix, "It is not evidence that the old root is a global public Nash equilibrium or subgame-perfect equilibrium", "sections/appendices.tex")
    require(appendix, "The earlier 20-draw perturbation exercise around the rejected vector is non-authoritative", "sections/appendices.tex")

    # The active result layer must be routed to the repaired v2.1 witness.
    results = json.loads(text(ROOT / "generated/results/canonical_results.json"))
    p = results["parameters"]
    expected = {"beta": 0.01, "gamma": 0.825, "tau": 0.35}
    for key, value in expected.items():
        if abs(float(p[key]) - value) > 1e-12:
            raise RuntimeError(f"stale repaired parameter {key}: {p[key]} != {value}")
    if not (results["computed"]["G"]["BR_slope"] < 0 < results["computed"]["B3"]["BR_slope"]):
        raise RuntimeError("repaired sign reversal missing from active result layer")
    if "stage4a_v21_repaired" not in text(ROOT / "scripts/generate_results.py"):
        raise RuntimeError("generated result pipeline is not routed to repaired Stage 4A source")
    if "public_two_sided_platform_hard_kill" in text(ROOT / "scripts/generate_results.py"):
        raise RuntimeError("old central-regime generator remains active")

    # Upstream authorities required by the quantifier gate.
    required_reviews = [
        "reviews/STAGE_04A_V21_REPAIRED_GLOBAL_CERTIFICATION_2026-09-08.md",
        "reviews/STAGE_06_V21_NOVELTY_REKILL_2026-09-08.md",
        "reviews/STAGE_07_V21_WELFARE_GENERALITY_2026-09-09.md",
        "reviews/STAGE_075_V21_FULL_THEORY_FREEZE_DECISION_2026-09-10.md",
    ]
    for rel in required_reviews:
        if not (ROOT / rel).exists():
            raise RuntimeError(f"missing upstream authority: {rel}")

    print("STAGE75A_QUANTIFIER_SCOPE_AUDIT: PASS")


if __name__ == "__main__":
    main()
