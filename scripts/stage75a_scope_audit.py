from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

ACTIVE = [
    ROOT / "paper/main.tex",
    ROOT / "sections/02_equilibrium.tex",
    ROOT / "sections/03_main_results.tex",
    ROOT / "sections/04_welfare.tex",
    ROOT / "sections/05_robustness.tex",
    ROOT / "sections/08_introduction.tex",
    ROOT / "sections/09_discussion.tex",
    ROOT / "sections/10_conclusion.tex",
    ROOT / "sections/appendices.tex",
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

    # Astra Stage-11R numerical-evidence repair: active prose must not retain the old
    # computational global-equilibrium qualification.
    for token in [
        "computationally certified all-regime witness",
        "computational global-equilibrium witness",
        "computational global-equilibrium existence witness",
        "all-regime computational certification at one baseline parameter vector",
        "Global best-response status is addressed separately",
    ]:
        forbid(joined, token, "active manuscript")

    abstract = active["paper/main.tex"]
    require(abstract, "symmetric regular beta-zero central-interior full-game stationary state", "paper/main.tex")
    require(abstract, "search evidence rather than a certified global-equilibrium or uniqueness result", "paper/main.tex")

    equilibrium = active["sections/02_equilibrium.tex"]
    require(equilibrium, "\\bar p_T(\\beta)", "sections/02_equilibrium.tex")
    require(equilibrium, "their stationary public investments generally differ", "sections/02_equilibrium.tex")
    require(equilibrium, "not a planner problem", "sections/02_equilibrium.tex")

    main_results = active["sections/03_main_results.tex"]
    require(main_results, "symmetric regular beta-zero central-interior full-game stationary state", "sections/03_main_results.tex")
    require(main_results, "support-side participation is strictly interior", "sections/03_main_results.tex")
    require(main_results, "common continuation neighborhood", "sections/03_main_results.tex")
    require(main_results, "Repaired all-regime computational search", "sections/03_main_results.tex")
    require(main_results, "search evidence", "sections/03_main_results.tex")
    require(main_results, "does not supply a rigorous upper bound on regret", "sections/03_main_results.tex")

    welfare = active["sections/04_welfare.tex"]
    require(welfare, "r_hm_h-\\frac{m_h^2}{2}", "sections/04_welfare.tex")
    require(welfare, "The all-regime numerical evaluator instead uses", "sections/04_welfare.tex")
    require(welfare, "direct national derivative", "sections/04_welfare.tex")
    require(welfare, "one-state-pair numerical comparison only", "sections/04_welfare.tex")

    robust = active["sections/05_robustness.tex"]
    require(robust, "search evidence, not a certified global-equilibrium result", "sections/05_robustness.tex")
    require(robust, "No broader genericity, equilibrium-existence, or numerical-certification claim is made", "sections/05_robustness.tex")

    discussion = active["sections/09_discussion.tex"]
    require(discussion, "search evidence, not a certified global-equilibrium or uniqueness result", "sections/09_discussion.tex")
    require(discussion, "not evidence of direct user poaching between public hubs", "sections/09_discussion.tex")
    require(discussion, "not a solved first-best problem", "sections/09_discussion.tex")

    conclusion = active["sections/10_conclusion.tex"]
    require(conclusion, "search evidence rather than a certified global-equilibrium, equilibrium-existence, or uniqueness result", "sections/10_conclusion.tex")
    require(conclusion, "symmetric regular beta-zero central-interior full-game stationary state", "sections/10_conclusion.tex")

    appendix = active["sections/appendices.tex"]
    require(appendix, "Support-side surplus with participation caps", "sections/appendices.tex")
    require(appendix, "\\varepsilon", "sections/appendices.tex")
    require(appendix, "search evidence only", "sections/appendices.tex")
    require(appendix, "It is not evidence that the old root is a global public Nash equilibrium or subgame-perfect equilibrium", "sections/appendices.tex")

    # Controlling amendment chain.
    t3_reopen = text(ROOT / "reviews/STAGE_075A_V21_T3_SYMMETRY_LIMITED_REOPEN_2026-09-10.md")
    t3_amend = text(ROOT / "theory_freeze_v21/STAGE_08_V21_T3_SYMMETRY_AMENDMENT_2026-09-10.md")
    astra_amend = text(ROOT / "theory_freeze_v21/STAGE_08_V21_ASTRA_WELFARE_EVIDENCE_AMENDMENT_2026-09-10.md")
    welfare_corr = text(ROOT / "reviews/STAGE_07_V21_SATURATED_PARTNER_SURPLUS_CORRECTION_2026-09-10.md")
    stage4_repair = text(ROOT / "reviews/STAGE_04A_V21_ASTRA_GLOBAL_EVIDENCE_REPAIR_2026-09-10.md")
    require(t3_reopen, "CERTIFICATION REGRESSION", "Stage-7.5A T3 limited reopen")
    require(t3_amend, "symmetric regular beta-zero full-game stationary state", "Stage-8 T3 amendment")
    require(astra_amend, "ALL-REGIME COMPUTATIONAL SEARCH EVIDENCE", "Stage-8 Astra amendment")
    require(astra_amend, "NO CERTIFIED GLOBAL REGRET BOUND", "Stage-8 Astra amendment")
    require(welfare_corr, "r_h*m_h - m_h^2/2", "Stage-7 welfare correction")
    require(stage4_repair, "SEARCH EVIDENCE", "Stage-4A evidence repair")

    results = json.loads(text(ROOT / "generated/results/canonical_results.json"))
    p = results["parameters"]
    expected = {"beta": 0.01, "gamma": 0.825, "tau": 0.35}
    for key, value in expected.items():
        if abs(float(p[key]) - value) > 1e-12:
            raise RuntimeError(f"stale repaired parameter {key}: {p[key]} != {value}")
    if not (results["computed"]["G"]["BR_slope"] < 0 < results["computed"]["B3"]["BR_slope"]):
        raise RuntimeError("repaired local sign reversal missing from active result layer")
    gs = results["computed"]["global_search"]
    if gs["evidence_level"] != "SEARCH EVIDENCE" or gs["certified_regret_upper_bound"] is not None:
        raise RuntimeError("Stage-11R global numerical evidence is overstated")
    if "SEARCH EVIDENCE" not in results["proof_status"]["repaired_witness"]:
        raise RuntimeError("generated proof-status layer does not reflect Astra evidence downgrade")
    if results["journal_target"] != "NOT SELECTED — STAGE 12 BLOCKED PENDING ASTRA RECHECK":
        raise RuntimeError("Stage 12 must remain blocked pending Astra recheck")

    required_reviews = [
        "reviews/STAGE_04A_V21_REPAIRED_GLOBAL_CERTIFICATION_2026-09-08.md",
        "reviews/STAGE_04A_V21_ASTRA_GLOBAL_EVIDENCE_REPAIR_2026-09-10.md",
        "reviews/STAGE_06_V21_NOVELTY_REKILL_2026-09-08.md",
        "reviews/STAGE_07_V21_WELFARE_GENERALITY_2026-09-09.md",
        "reviews/STAGE_07_V21_SATURATED_PARTNER_SURPLUS_CORRECTION_2026-09-10.md",
        "reviews/STAGE_075_V21_FULL_THEORY_FREEZE_DECISION_2026-09-10.md",
    ]
    for rel in required_reviews:
        if not (ROOT / rel).exists():
            raise RuntimeError(f"missing upstream authority: {rel}")

    print("STAGE75A_QUANTIFIER_SCOPE_AUDIT: PASS — ASTRA STAGE11R AMENDMENTS ACTIVE")


if __name__ == "__main__":
    main()
