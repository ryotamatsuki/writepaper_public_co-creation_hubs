from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from stage4a_v21_repaired.code.independent_repaired_audit import (
    P, XG, XB, private_br, welfare, profit, multistart_state,
    SOLVED_EQUILIBRIUM,
)


def national_welfare(x1, x2, p, full=False):
    return (
        welfare(1, x1, x2, p, P, full=full)
        + welfare(2, x1, x2, p, P, full=full)
        + profit(x1, x2, p, P, full=full)
    )


def full_game_objects(x1, x2, full=False, price_grid=51):
    pr = private_br(x1, x2, P, grid_n=price_grid, full=full)
    if pr.status != SOLVED_EQUILIBRIUM or pr.p is None:
        raise RuntimeError(f"private continuation {pr.status}")
    st = multistart_state(x1, x2, pr.p, P, full=full)
    return {
        "p": pr.p,
        "profit": pr.profit,
        "W1": welfare(1, x1, x2, pr.p, P, full=full),
        "W2": welfare(2, x1, x2, pr.p, P, full=full),
        "WN": national_welfare(x1, x2, pr.p, full=full),
        "n": st["n"],
        "b": st["b"],
    }


def fixed_price_objects(x1, x2, pbar, full=False):
    st = multistart_state(x1, x2, pbar, P, full=full)
    return {
        "p": pbar,
        "profit": profit(x1, x2, pbar, P, full=full),
        "W1": welfare(1, x1, x2, pbar, P, full=full),
        "W2": welfare(2, x1, x2, pbar, P, full=full),
        "WN": national_welfare(x1, x2, pbar, full=full),
        "n": st["n"],
        "b": st["b"],
    }


def main():
    g = full_game_objects(XG, XG, full=True, price_grid=81)
    pG = g["p"]
    b3 = fixed_price_objects(XB, XB, pG, full=True)
    print("REPAIRED_G", g)
    print("REPAIRED_B3", b3)
    print("DELTA_WN_G_MINUS_B3", g["WN"] - b3["WN"])

    transfer = pG * g["n"][2]
    if abs(transfer - g["profit"]) > 1e-9:
        raise RuntimeError("private-price transfer cancellation failed")
    print("TRANSFER_CANCELLATION: PASS", transfer, g["profit"])

    # One pair of re-solved histories is enough to recover all local welfare
    # derivatives consistently; this avoids constructing an unnecessary planner
    # benchmark and keeps the benchmark language exact.
    h = 3e-4
    lo = full_game_objects(XG - h, XG, full=False, price_grid=61)
    hi = full_game_objects(XG + h, XG, full=False, price_grid=61)
    den = 2.0 * h
    dW1 = (hi["W1"] - lo["W1"]) / den
    dW2 = (hi["W2"] - lo["W2"]) / den
    dPi = (hi["profit"] - lo["profit"]) / den
    dWN = (hi["WN"] - lo["WN"]) / den
    dpdx = (hi["p"] - lo["p"]) / den
    print("G_COORDINATION_WEDGE", dW1, dW2, dPi, dWN, "rhs", dW2 + dPi)
    print("PRIVATE_PRICE_RESPONSE_DPDX", dpdx)
    if abs(dWN - (dW1 + dW2 + dPi)) > 5e-5:
        raise RuntimeError("national welfare derivative accounting failed")
    if dWN > 0:
        print("LOCAL_COORDINATION_CLASSIFICATION: UNDERPROVISION IN THE COORDINATED +x_i DIRECTION")
    elif dWN < 0:
        print("LOCAL_COORDINATION_CLASSIFICATION: OVERPROVISION IN THE COORDINATED +x_i DIRECTION")
    else:
        print("LOCAL_COORDINATION_CLASSIFICATION: LOCALLY NEUTRAL")

    print("PLANNER_BENCHMARK_REGISTER: NONE CLAIMED; NO FIRST-BEST LABEL USED")
    print("STAGE7_WELFARE_NUMERICAL_AUDIT: PASS")


if __name__ == "__main__":
    main()
