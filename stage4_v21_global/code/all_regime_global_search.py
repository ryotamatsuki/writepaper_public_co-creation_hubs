"""Stage 4 v2.1 all-regime global-equilibrium repair search.

This construction-stage solver uses the primitive all-route evaluator preserved by
Stage 4A. It never maps branch failure to low payoff. Any unresolved continuation
raises and blocks a positive verdict.

The repair being tested is NOT a model change. It asks whether the existing remote
public-access cost tau can be placed in a nonempty primitive region where the
previously certified symmetric stationary pair is also a global public best
response. The central stationary equations do not contain tau while the rival
public route is inactive.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from stage4a_v2.code.independent_adversarial_audit import (
    P as BASE, XG, XB, PG, fixed_point, multistart_state, regional_welfare
)

SOLVED_EQUILIBRIUM = "SOLVED_EQUILIBRIUM"
MULTIPLE_EQUILIBRIA = "MULTIPLE_EQUILIBRIA"
UNRESOLVED = "UNRESOLVED"
NUMERICAL_FAILURE = "NUMERICAL_FAILURE"


@dataclass
class PriceResult:
    status: str
    p: float | None
    profit: float | None


@dataclass
class BRResult:
    status: str
    x: float | None
    welfare: float | None
    candidate_welfare: float | None
    gain_over_candidate: float | None


def par_tau(tau: float):
    p = BASE.copy()
    p["tau"] = float(tau)
    return p


def _state(x1, x2, p, par, rigorous=False):
    if rigorous:
        return multistart_state(x1, x2, p, par)
    return fixed_point(x1, x2, p, (.4, .4, .5), par)


def _profit(x1, x2, p, par, rigorous=False):
    st = _state(x1, x2, p, par, rigorous)
    return float(p * st["n"][2])


def private_best_response(x1, x2, par, grid_n=41, rigorous=False):
    upper = par["v"] + par["alpha"] - par["kT"]
    grid = np.linspace(0.0, upper, grid_n)
    vals = []
    for p in grid:
        try:
            vals.append(_profit(x1, x2, float(p), par, rigorous))
        except Exception as exc:
            raise RuntimeError(f"UNRESOLVED private continuation at p={p}: {exc}")
    vals = np.asarray(vals)
    candidates = []
    for j, value in enumerate(vals):
        left = vals[j-1] if j else -np.inf
        right = vals[j+1] if j+1 < len(vals) else -np.inf
        if value >= left and value >= right:
            candidates.append((float(grid[j]), float(value)))
            lo = grid[max(0, j-1)]
            hi = grid[min(len(grid)-1, j+1)]
            if hi > lo:
                try:
                    opt = minimize_scalar(
                        lambda z: -_profit(x1, x2, float(z), par, rigorous),
                        bounds=(lo, hi), method="bounded",
                        options={"xatol": 2e-8 if rigorous else 2e-6},
                    )
                except Exception as exc:
                    raise RuntimeError(f"UNRESOLVED private refinement: {exc}")
                candidates.append((float(opt.x), float(-opt.fun)))
    if not candidates:
        return PriceResult(UNRESOLVED, None, None)
    candidates.sort(key=lambda z: z[1], reverse=True)
    best = candidates[0]
    tied = [c for c in candidates if abs(c[1]-best[1]) <= 1e-8]
    distinct = []
    for c in tied:
        if not any(abs(c[0]-d[0]) <= 2e-5 for d in distinct):
            distinct.append(c)
    status = MULTIPLE_EQUILIBRIA if len(distinct) > 1 else SOLVED_EQUILIBRIUM
    return PriceResult(status, best[0], best[1])


def welfare_G(i, x1, x2, par, price_grid=41, rigorous=False):
    pr = private_best_response(x1, x2, par, price_grid, rigorous)
    if pr.status not in (SOLVED_EQUILIBRIUM, MULTIPLE_EQUILIBRIA) or pr.p is None:
        raise RuntimeError(f"private continuation {pr.status}")
    if pr.status == MULTIPLE_EQUILIBRIA:
        raise RuntimeError("MULTIPLE_EQUILIBRIA private continuation")
    return float(regional_welfare(i, x1, x2, pr.p, par)), pr.p


def welfare_B3(i, x1, x2, pbar, par, rigorous=False):
    _state(x1, x2, pbar, par, rigorous)
    return float(regional_welfare(i, x1, x2, pbar, par))


def _global_public_br(mode, rival_x, candidate_x, par, pbar=None,
                      x_grid_n=61, price_grid=31, rigorous=False):
    xs = np.linspace(0.0, 1.0, x_grid_n)
    vals = []
    for x in xs:
        try:
            if mode == "G":
                w, _ = welfare_G(1, float(x), rival_x, par, price_grid, rigorous)
            else:
                w = welfare_B3(1, float(x), rival_x, pbar, par, rigorous)
            vals.append(w)
        except Exception as exc:
            raise RuntimeError(f"UNRESOLVED public deviation x={x}: {exc}")
    vals = np.asarray(vals)
    local_idx = []
    for j, value in enumerate(vals):
        left = vals[j-1] if j else -np.inf
        right = vals[j+1] if j+1 < len(vals) else -np.inf
        if value >= left and value >= right:
            local_idx.append(j)
    candidates = [(float(xs[j]), float(vals[j])) for j in local_idx]
    for j in local_idx:
        lo = xs[max(0, j-1)]
        hi = xs[min(len(xs)-1, j+1)]
        if hi <= lo:
            continue
        def obj(x):
            if mode == "G":
                return -welfare_G(1, float(x), rival_x, par, price_grid, rigorous)[0]
            return -welfare_B3(1, float(x), rival_x, pbar, par, rigorous)
        try:
            opt = minimize_scalar(obj, bounds=(lo, hi), method="bounded",
                                  options={"xatol": 2e-6 if rigorous else 2e-4})
        except Exception as exc:
            raise RuntimeError(f"UNRESOLVED public refinement: {exc}")
        candidates.append((float(opt.x), float(-opt.fun)))
    if not candidates:
        return BRResult(UNRESOLVED, None, None, None, None)
    best = max(candidates, key=lambda z: z[1])
    if mode == "G":
        cand_w = welfare_G(1, candidate_x, rival_x, par, price_grid, rigorous)[0]
    else:
        cand_w = welfare_B3(1, candidate_x, rival_x, pbar, par, rigorous)
    return BRResult(SOLVED_EQUILIBRIUM, best[0], best[1], cand_w, best[1]-cand_w)


def coarse_tau_screen():
    rows = []
    for tau in (0.05, 0.20, 0.33, 0.35, 0.50):
        par = par_tau(tau)
        g = _global_public_br("G", XG, XG, par, x_grid_n=41, price_grid=25)
        pg = private_best_response(XG, XG, par, grid_n=61, rigorous=True)
        if pg.status != SOLVED_EQUILIBRIUM:
            raise RuntimeError(f"tau={tau}: private candidate status {pg.status}")
        b = _global_public_br("B3", XB, XB, par, pbar=pg.p,
                              x_grid_n=61, rigorous=False)
        rows.append((tau, g, b, pg))
        print("TAU_SCREEN", tau,
              "G_BR", g.x, "G_GAIN", g.gain_over_candidate,
              "B3_BR", b.x, "B3_GAIN", b.gain_over_candidate,
              "pG", pg.p)
    return rows


def rigorous_candidate(tau):
    par = par_tau(tau)
    rival_globally_dominated = par["kL"] + par["tau"] >= par["v"] + par["alpha"]
    print("RIVAL_GLOBAL_DOMINANCE", rival_globally_dominated,
          "remote_cost", par["kL"]+par["tau"],
          "max_quality", par["v"]+par["alpha"])

    pg = private_best_response(XG, XG, par, grid_n=121, rigorous=True)
    if pg.status != SOLVED_EQUILIBRIUM:
        raise RuntimeError(f"rigorous private candidate {pg.status}")
    print("RIGOROUS_PG", pg.p, "delta_from_old", pg.p-PG)

    g = _global_public_br("G", XG, XG, par, x_grid_n=101,
                          price_grid=41, rigorous=False)
    b = _global_public_br("B3", XB, XB, par, pbar=pg.p,
                          x_grid_n=151, rigorous=False)
    print("DENSE_G", g)
    print("DENSE_B3", b)

    alt_x = g.x
    p_alt = private_best_response(alt_x, XG, par, grid_n=121, rigorous=True)
    if p_alt.status != SOLVED_EQUILIBRIUM:
        raise RuntimeError(f"rigorous alternative private status {p_alt.status}")
    wg_c = float(regional_welfare(1, XG, XG, pg.p, par))
    wg_a = float(regional_welfare(1, alt_x, XG, p_alt.p, par))
    wb_c = float(regional_welfare(1, XB, XB, pg.p, par))
    wb_a = float(regional_welfare(1, b.x, XB, pg.p, par))
    print("MULTISTART_CONFIRM_G", "candidate", wg_c, "alt_x", alt_x,
          "alt", wg_a, "gain", wg_a-wg_c, "p_alt", p_alt.p)
    print("MULTISTART_CONFIRM_B3", "candidate", wb_c, "alt_x", b.x,
          "alt", wb_a, "gain", wb_a-wb_c)

    g_pass = (wg_a-wg_c <= 1e-4 and abs(alt_x-XG) <= 5e-3)
    b_pass = (wb_a-wb_c <= 1e-4 and abs(b.x-XB) <= 5e-3)
    return dict(tau=tau, rival_globally_dominated=rival_globally_dominated,
                pg=pg, g=g, b=b, wg_c=wg_c, wg_a=wg_a,
                wb_c=wb_c, wb_a=wb_a, g_pass=g_pass, b_pass=b_pass)


def main():
    rows = coarse_tau_screen()
    candidates = [tau for tau, g, b, pg in rows
                  if tau >= 0.33 and g.gain_over_candidate <= 5e-4
                  and b.gain_over_candidate <= 5e-4]
    if not candidates:
        print("STAGE4_GLOBAL_SEARCH: NO PASSING TAU IN SCREEN")
        raise SystemExit(2)
    tau = min(candidates)
    out = rigorous_candidate(tau)
    print("SELECTED_TAU", tau)
    print("G_GLOBAL_PASS", out["g_pass"])
    print("B3_GLOBAL_PASS", out["b_pass"])
    if not (out["rival_globally_dominated"] and out["g_pass"] and out["b_pass"]):
        print("STAGE4_GLOBAL_SEARCH: FAIL")
        raise SystemExit(3)
    print("STAGE4_GLOBAL_SEARCH: CONSTRUCTION PASS; ROUTE TO STAGE 4A")


if __name__ == "__main__":
    main()
