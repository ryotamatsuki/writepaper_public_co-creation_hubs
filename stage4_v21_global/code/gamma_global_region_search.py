"""Search existing primitive gamma for a genuine global sign-reversal equilibrium.

No model primitive is added or redefined. tau is fixed at 0.35, which implies
kL+tau=0.62 > v+alpha=0.60, so a nonresident's rival public hub is globally
dominated by the outside option at every history. We then vary only the existing
public investment cost curvature gamma and ask whether the smooth central
sign-reversal stationary pair is also a global public best response in G and B3.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from public_two_sided_platform_hard_kill.code.numerical_hard_kill import (
    BASE as OLD_BASE, derivs, pstar
)
from stage4_v21_global.code.all_regime_global_search import (
    SOLVED_EQUILIBRIUM, private_best_response, _global_public_br
)

TAU = 0.35


def old_par(gamma):
    p = OLD_BASE.copy()
    p["tau"] = TAU
    p["gamma"] = float(gamma)
    return p


def new_par(gamma):
    o = old_par(gamma)
    return dict(v=o["v"], alpha=o["alpha"], beta=o["beta"], rho=o["rho"],
                rhoT=o["rhoT"], kL=o["k"], kT=o["kT"], tau=o["tau"], gamma=o["gamma"])


def central_roots(P, mode, pbar=None):
    xs = np.linspace(.20, .98, 80)
    vals = []
    for x in xs:
        try:
            g = derivs(1, float(x), P, mode, pbar, 2e-4)[0]
            vals.append(float(g) if np.isfinite(g) else np.nan)
        except Exception:
            vals.append(np.nan)
    roots = []
    for j in range(len(xs)-1):
        if not (np.isfinite(vals[j]) and np.isfinite(vals[j+1])):
            continue
        if vals[j] == 0 or vals[j]*vals[j+1] < 0:
            try:
                f = lambda z: derivs(1, float(z), P, mode, pbar, 2e-4)[0]
                r = brentq(f, float(xs[j]), float(xs[j+1]), xtol=3e-6)
                d = derivs(1, float(r), P, mode, pbar, 1.5e-4)
                slope = -d[2]/d[1]
                if np.isfinite(d[1]) and d[1] < 0:
                    if not any(abs(r-a[0]) < 1e-4 for a in roots):
                        roots.append((float(r), tuple(float(v) for v in d), float(slope)))
            except Exception:
                pass
    return roots


def local_pair(gamma):
    P = old_par(gamma)
    grows = central_roots(P, "G", None)
    out = []
    for xg, dg, sg in grows:
        try:
            pg = float(pstar(xg, xg, P))
            brows = central_roots(P, "B3", pg)
        except Exception:
            continue
        for xb, db, sb in brows:
            if sg < 0 < sb:
                out.append(dict(xg=xg, dg=dg, sg=sg, pg=pg, xb=xb, db=db, sb=sb))
    return out


def coarse_global(pair, gamma):
    par = new_par(gamma)
    xg, xb = pair["xg"], pair["xb"]
    # Re-solve the full private continuation all-regime at the G candidate.
    pg = private_best_response(xg, xg, par, grid_n=61, rigorous=True)
    if pg.status != SOLVED_EQUILIBRIUM:
        raise RuntimeError(f"private candidate {pg.status}")
    g = _global_public_br("G", xg, xg, par, x_grid_n=41, price_grid=25)
    b = _global_public_br("B3", xb, xb, par, pbar=pg.p, x_grid_n=61, price_grid=25)
    return pg, g, b


def main():
    assert new_par(.9)["kL"] + TAU > new_par(.9)["v"] + new_par(.9)["alpha"]
    passing = []
    for gamma in (.45, .50, .55, .60, .65, .70, .75, .80, .85, .90, 1.00, 1.10):
        pairs = local_pair(gamma)
        print("GAMMA_LOCAL", gamma, "pairs", len(pairs))
        for pair in pairs:
            try:
                pg, g, b = coarse_global(pair, gamma)
            except Exception as exc:
                print("GAMMA_UNRESOLVED", gamma, repr(exc))
                continue
            print("GAMMA_SCREEN", gamma,
                  "xG", pair["xg"], "sG", pair["sg"],
                  "xB", pair["xb"], "sB", pair["sb"],
                  "pG", pg.p,
                  "G_BR", g.x, "G_GAIN", g.gain_over_candidate,
                  "B3_BR", b.x, "B3_GAIN", b.gain_over_candidate)
            if (g.gain_over_candidate <= 5e-4 and b.gain_over_candidate <= 5e-4
                    and abs(g.x-pair["xg"]) <= 8e-3
                    and abs(b.x-pair["xb"]) <= 8e-3):
                passing.append((gamma, pair, pg, g, b))
    if not passing:
        print("GAMMA_GLOBAL_REGION_SEARCH: NO CONSTRUCTION PASS")
        raise SystemExit(2)
    gamma, pair, pg, g, b = passing[0]
    print("SELECTED_GAMMA", gamma)
    print("SELECTED_PAIR", pair)
    print("SELECTED_PG", pg)
    print("SELECTED_GLOBAL_G", g)
    print("SELECTED_GLOBAL_B3", b)
    print("GAMMA_GLOBAL_REGION_SEARCH: CONSTRUCTION PASS")


if __name__ == "__main__":
    main()
