from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from stage4a_v21_repaired.code.independent_repaired_audit import (
    P,
    XG,
    XB,
    SOLVED_EQUILIBRIUM,
    private_br,
    local_derivatives,
    welfare,
    profit,
    multistart_state,
)


def _deriv(f, x: float, h: float = 3e-4) -> float:
    return (f(x + h) - f(x - h)) / (2.0 * h)


def build_results() -> dict:
    pg = private_br(XG, XG, P, grid_n=81, full=True)
    if pg.status != SOLVED_EQUILIBRIUM or pg.p is None:
        raise RuntimeError(f"repaired private continuation unresolved: {pg.status}")
    pbar = float(pg.p)

    g_h11, g_h12, g_slope = local_derivatives("G", XG, h=5e-4)
    b_h11, b_h12, b_slope = local_derivatives("B3", XB, pbar=pbar, h=5e-4)

    g_state = multistart_state(XG, XG, pbar, P, full=True)
    b_state = multistart_state(XB, XB, pbar, P, full=True)

    g_wi = float(welfare(1, XG, XG, pbar, P, full=True))
    b_wi = float(welfare(1, XB, XB, pbar, P, full=True))
    g_pi = float(profit(XG, XG, pbar, P, full=True))
    b_pi = float(profit(XB, XB, pbar, P, full=True))
    g_wn = float(2.0 * g_wi + g_pi)
    b_wn = float(2.0 * b_wi + b_pi)

    def reduced_p(x: float) -> float:
        pr = private_br(x, XG, P, grid_n=61, full=False)
        if pr.status != SOLVED_EQUILIBRIUM or pr.p is None:
            raise RuntimeError(f"private derivative continuation unresolved: {pr.status}")
        return float(pr.p)

    def rival_welfare(x: float) -> float:
        p = reduced_p(x)
        return float(welfare(2, x, XG, p, P, full=False))

    def private_profit(x: float) -> float:
        p = reduced_p(x)
        return float(profit(x, XG, p, P, full=False))

    rival_dw = _deriv(rival_welfare, XG)
    dpi = _deriv(private_profit, XG)
    dpdx = _deriv(reduced_p, XG)

    return {
        "_generated": "DO NOT EDIT — GENERATED FILE",
        "sources": [
            "stage4a_v21_repaired/code/independent_repaired_audit.py",
            "reviews/STAGE_04A_V21_REPAIRED_GLOBAL_CERTIFICATION_2026-09-08.md",
            "reviews/STAGE_07_V21_WELFARE_GENERALITY_2026-09-09.md",
        ],
        "parameters": {k: float(v) for k, v in P.items()},
        "computed": {
            "G": {
                "x": float(XG),
                "p_T": pbar,
                "own_second": float(g_h11),
                "cross_second": float(g_h12),
                "BR_slope": float(g_slope),
            },
            "B3": {
                "x": float(XB),
                "p_T": pbar,
                "own_second": float(b_h11),
                "cross_second": float(b_h12),
                "BR_slope": float(b_slope),
            },
            "private_price_response": float(dpdx),
            "welfare": {
                "G": {
                    "x": float(XG),
                    "p_T": pbar,
                    "n_T": float(g_state["n"][2]),
                    "W_i": g_wi,
                    "Pi_T": g_pi,
                    "W_N": g_wn,
                },
                "B3": {
                    "x": float(XB),
                    "p_T": pbar,
                    "n_T": float(b_state["n"][2]),
                    "W_i": b_wi,
                    "Pi_T": b_pi,
                    "W_N": b_wn,
                },
            },
            "welfare_difference": float(g_wn - b_wn),
            "coordination": {
                "rival_welfare": float(rival_dw),
                "private_profit": float(dpi),
                "national_wedge": float(rival_dw + dpi),
            },
        },
        "proof_status": {
            "analytic_headline": "LOCAL SUFFICIENT-CONDITION THEOREM ON REGULAR STATIONARY BRANCHES",
            "repaired_witness": "ALL-REGIME COMPUTATIONAL GLOBAL-EQUILIBRIUM CERTIFICATION",
            "old_exact_certificate": "LOCAL STATIONARY-ROOT DIAGNOSTIC ONLY — NOT GLOBAL EQUILIBRIUM AUTHORITY",
            "robustness": "NO BROAD FUNCTION-CLASS ROBUSTNESS THEOREM CLAIMED",
        },
        "journal_target": "JPET — DEFENSIBLE BUT BORDERLINE",
    }


def main() -> None:
    out = ROOT / "generated/results/canonical_results.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(build_results(), indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(out)


if __name__ == "__main__":
    main()
