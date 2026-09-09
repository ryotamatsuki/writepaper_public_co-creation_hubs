from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

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
        "regional": st["regional"],
    }


def fixed_price_objects(x1, x2, pbar, full=False):
    st = multistart_state(x1, x2, pbar, P, full=full)
    pi = profit(x1, x2, pbar, P, full=full)
    return {
        "p": pbar,
        "profit": pi,
        "W1": welfare(1, x1, x2, pbar, P, full=full),
        "W2": welfare(2, x1, x2, pbar, P, full=full),
        "WN": national_welfare(x1, x2, pbar, full=full),
        "n": st["n"],
        "b": st["b"],
        "regional": st["regional"],
    }


def deriv(f, x, h=3e-4):
    return (f(x + h) - f(x - h)) / (2.0 * h)


def symmetric_restricted_coordination(mode, pbar=None):
    # Precisely restricted benchmark, not first best: the planner chooses only a
    # common public investment x1=x2=x and inherits either endogenous private
    # pricing (G) or the matched fixed private price (B3).
    def obj(x):
        if mode == "G":
            return full_game_objects(float(x), float(x), full=False, price_grid=31)["WN"]
        return fixed_price_objects(float(x), float(x), pbar, full=False)["WN"]

    grid = np.linspace(0.0, 1.0, 21)
    vals = np.array([obj(float(x)) for x in grid])
    candidates = [(float(grid[j]), float(vals[j])) for j in range(len(grid))]
    for j in range(1, len(grid) - 1):
        if vals[j] >= vals[j - 1] and vals[j] >= vals[j + 1]:
            opt = minimize_scalar(
                lambda z: -obj(float(z)),
                bounds=(grid[j - 1], grid[j + 1]),
                method="bounded",
                options={"xatol": 2e-5},
            )
            candidates.append((float(opt.x), float(-opt.fun)))
    xstar, valstar = max(candidates, key=lambda t: t[1])
    if mode == "G":
        chk = full_game_objects(xstar, xstar, full=True, price_grid=81)
    else:
        chk = fixed_price_objects(xstar, xstar, pbar, full=True)
    return xstar, valstar, chk["WN"]


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

    def W1_reduced(x): return full_game_objects(x, XG, full=False, price_grid=51)["W1"]
    def W2_reduced(x): return full_game_objects(x, XG, full=False, price_grid=51)["W2"]
    def Pi_reduced(x): return full_game_objects(x, XG, full=False, price_grid=51)["profit"]
    def WN_reduced(x): return full_game_objects(x, XG, full=False, price_grid=51)["WN"]
    dW1 = deriv(W1_reduced, XG)
    dW2 = deriv(W2_reduced, XG)
    dPi = deriv(Pi_reduced, XG)
    dWN = deriv(WN_reduced, XG)
    print("G_COORDINATION_WEDGE", dW1, dW2, dPi, dWN, "rhs", dW2 + dPi)
    if abs(dWN - (dW1 + dW2 + dPi)) > 5e-5:
        raise RuntimeError("national welfare derivative accounting failed")

    def p_reduced(x): return full_game_objects(x, XG, full=False, price_grid=61)["p"]
    dpdx = deriv(p_reduced, XG)
    print("PRIVATE_PRICE_RESPONSE_DPDX", dpdx)

    xcg, coarse_g, full_g = symmetric_restricted_coordination("G")
    xcb, coarse_b, full_b = symmetric_restricted_coordination("B3", pbar=pG)
    print("SYMMETRIC_RESTRICTED_COORD_G", xcg, coarse_g, full_g, "decentralized", XG)
    print("SYMMETRIC_RESTRICTED_COORD_B3", xcb, coarse_b, full_b, "decentralized", XB)

    if full_g <= g["WN"] + 1e-6:
        print("G_COORDINATION_CLASSIFICATION: NO MATERIAL SYMMETRIC WELFARE GAIN")
    elif xcg > XG:
        print("G_COORDINATION_CLASSIFICATION: UNDERINVESTMENT RELATIVE TO SYMMETRIC RESTRICTED-INSTRUMENT COORDINATION")
    else:
        print("G_COORDINATION_CLASSIFICATION: OVERINVESTMENT RELATIVE TO SYMMETRIC RESTRICTED-INSTRUMENT COORDINATION")

    print("STAGE7_WELFARE_NUMERICAL_AUDIT: PASS")


if __name__ == "__main__":
    main()
