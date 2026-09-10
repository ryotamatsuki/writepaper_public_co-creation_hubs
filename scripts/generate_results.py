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
    GLOBAL_EVIDENCE_LEVEL,
    private_br,
    public_br,
    local_derivatives,
    welfare,
    profit,
    multistart_state,
)


def _deriv(f, x: float, h: float) -> float:
    return (f(x + h) - f(x - h)) / (2.0 * h)


def _second(f, x: float, h: float) -> float:
    return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h)


def build_results() -> dict:
    pg = private_br(XG, XG, P, grid_n=81, full=True)
    if pg.status != SOLVED_EQUILIBRIUM or pg.p is None:
        raise RuntimeError(f"repaired private continuation unresolved: {pg.status}")
    pbar = float(pg.p)

    g_h11, g_h12, g_slope = local_derivatives("G", XG, h=5e-4)
    b_h11, b_h12, b_slope = local_derivatives("B3", XB, pbar=pbar, h=5e-4)

    g_state = multistart_state(XG, XG, pbar, P, full=True)
    b_state = multistart_state(XB, XB, pbar, P, full=True)

    g_w1 = float(welfare(1, XG, XG, pbar, P, full=True))
    g_w2 = float(welfare(2, XG, XG, pbar, P, full=True))
    b_w1 = float(welfare(1, XB, XB, pbar, P, full=True))
    b_w2 = float(welfare(2, XB, XB, pbar, P, full=True))
    g_pi = float(profit(XG, XG, pbar, P, full=True))
    b_pi = float(profit(XB, XB, pbar, P, full=True))
    g_wn = float(g_w1 + g_w2 + g_pi)
    b_wn = float(b_w1 + b_w2 + b_pi)

    def reduced_p(x: float) -> float:
        pr = private_br(x, XG, P, grid_n=61, full=False)
        if pr.status != SOLVED_EQUILIBRIUM or pr.p is None:
            raise RuntimeError(f"private derivative continuation unresolved: {pr.status}")
        return float(pr.p)

    def own_welfare(x: float) -> float:
        p = reduced_p(x)
        return float(welfare(1, x, XG, p, P, full=False))

    def rival_welfare(x: float) -> float:
        p = reduced_p(x)
        return float(welfare(2, x, XG, p, P, full=False))

    def private_profit(x: float) -> float:
        p = reduced_p(x)
        return float(profit(x, XG, p, P, full=False))

    def national_welfare(x: float) -> float:
        p = reduced_p(x)
        return float(
            welfare(1, x, XG, p, P, full=False)
            + welfare(2, x, XG, p, P, full=False)
            + profit(x, XG, p, P, full=False)
        )

    welfare_h = 1e-3
    own_dw = _deriv(own_welfare, XG, welfare_h)
    rival_dw = _deriv(rival_welfare, XG, welfare_h)
    dpi = _deriv(private_profit, XG, welfare_h)
    national_direct = _deriv(national_welfare, XG, welfare_h)
    national_components = float(own_dw + rival_dw + dpi)
    dpdx = _deriv(reduced_p, XG, welfare_h)

    # Stationarity is checked numerically rather than inferred from stored candidate values.
    private_h = 5e-5
    private_foc = _deriv(lambda z: float(profit(XG, XG, z, P, full=False)), pbar, private_h)
    private_soc = _second(lambda z: float(profit(XG, XG, z, P, full=False)), pbar, private_h)
    public_h = 7e-4
    g_public_foc = _deriv(own_welfare, XG, public_h)
    b_public_foc = _deriv(
        lambda z: float(welfare(1, z, XB, pbar, P, full=False)), XB, public_h
    )

    # These are search outputs, not certified regret bounds.
    g_search = public_br("G", XG, XG, par=P, x_grid_n=51, price_grid=31)
    b_search = public_br("B3", XB, XB, pbar=pbar, par=P, x_grid_n=61, price_grid=31)

    return {
        "_generated": "DO NOT EDIT — GENERATED FILE",
        "provenance": {
            "stage8_merge_sha": "ad927ca783a6123ea4fc6f55f65598ebd6ab583b",
            "stage75a_certified_input_sha": "eeb48a3dd76ab6f43d5de175b12c3f374746db0f",
            "stage8_freeze": "theory_freeze_v21/CANONICAL_THEORY_FREEZE_2026-09-10.md",
            "stage75a_scope_ledger": "reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md",
            "stage11r_repair": "reviews/STAGE_11R_ASTRA_REPAIR_REPORT.md",
        },
        "sources": [
            "stage4a_v21_repaired/code/independent_repaired_audit.py",
            "public_two_sided_platform_welfare_generality/PARTNER_SURPLUS_DERIVATION.md",
            "reviews/STAGE_11R_ASTRA_REPAIR_REPORT.md",
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
            "stationarity": {
                "participation_residual_inf_G": float(g_state["residual"]),
                "participation_residual_inf_B3": float(b_state["residual"]),
                "private_price_foc": float(private_foc),
                "private_price_soc": float(private_soc),
                "G_public_own_foc": float(g_public_foc),
                "B3_public_own_foc": float(b_public_foc),
                "private_step": float(private_h),
                "public_step": float(public_h),
            },
            "global_search": {
                "evidence_level": GLOBAL_EVIDENCE_LEVEL,
                "public_interval": [0.0, 1.0],
                "G_public_grid_points": 51,
                "B3_public_grid_points": 61,
                "public_local_refinement": True,
                "G_deviation_private_reoptimization": True,
                "B3_fixed_matched_price": True,
                "G_best_detected_x": float(g_search.x),
                "G_best_detected_gain": float(g_search.gain),
                "B3_best_detected_x": float(b_search.x),
                "B3_best_detected_gain": float(b_search.gain),
                "certified_regret_upper_bound": None,
            },
            "private_price_response": float(dpdx),
            "welfare": {
                "G": {
                    "x": float(XG),
                    "p_T": pbar,
                    "n_T": float(g_state["n"][2]),
                    "W_i": g_w1,
                    "Pi_T": g_pi,
                    "W_N": g_wn,
                    "support_gross": [float(v) for v in g_state["r"]],
                    "support_mass": [float(v) for v in g_state["b"]],
                },
                "B3": {
                    "x": float(XB),
                    "p_T": pbar,
                    "n_T": float(b_state["n"][2]),
                    "W_i": b_w1,
                    "Pi_T": b_pi,
                    "W_N": b_wn,
                    "support_gross": [float(v) for v in b_state["r"]],
                    "support_mass": [float(v) for v in b_state["b"]],
                },
            },
            "welfare_difference": float(g_wn - b_wn),
            "coordination": {
                "step": float(welfare_h),
                "own_welfare": float(own_dw),
                "rival_welfare": float(rival_dw),
                "private_profit": float(dpi),
                "sum_components": national_components,
                "national_direct": float(national_direct),
                "decomposition_error": float(national_direct - national_components),
            },
        },
        "proof_status": {
            "analytic_headline": "LOCAL SUFFICIENT-CONDITION THEOREM ON REGULAR STATIONARY BRANCHES",
            "repaired_witness": "ALL-REGIME COMPUTATIONAL SEARCH EVIDENCE — NO CERTIFIED GLOBAL REGRET BOUND",
            "welfare_numerics": "LOCAL/ONE-STATE NUMERICAL RESULTS AT REPORTED STATIONARY SEARCH CANDIDATES",
            "old_exact_certificate": "LOCAL STATIONARY-ROOT DIAGNOSTIC ONLY — NOT GLOBAL EQUILIBRIUM AUTHORITY",
            "robustness": "NO BROAD FUNCTION-CLASS ROBUSTNESS THEOREM CLAIMED",
        },
        "journal_target": "NOT SELECTED — STAGE 12 BLOCKED PENDING ASTRA RECHECK",
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
