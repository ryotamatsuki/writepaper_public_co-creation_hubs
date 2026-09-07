"""Search existing (beta, gamma) primitives for a genuine global sign-reversal equilibrium.

No primitive or timing is added or redefined. We keep tau=.35 so a nonresident's
rival public hub is globally dominated by the outside option (kL+tau=.62 >
v+alpha=.60). We vary only the existing network-effect strength beta and public
investment cost curvature gamma. A construction PASS requires:

1. a regular symmetric central stationary pair with B3 slope > 0 > G slope;
2. all-regime private continuation solved at the G candidate;
3. global public best-response search over x_i in [0,1] returns the candidate in
   both G and the matched-price B3 benchmark, within construction tolerances.

Stage 4A must independently certify any selected candidate; this script is only
the construction-stage search and fail-closed screen.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from public_two_sided_platform_hard_kill.code.numerical_hard_kill import BASE as OLD_BASE, pstar
from stage4_v21_global.code.gamma_global_region_search import central_roots
from stage4_v21_global.code.all_regime_global_search import (
    SOLVED_EQUILIBRIUM,
    private_best_response,
    _global_public_br,
)

TAU = .35
BETAS = (.001, .002, .005, .01, .015, .02, .03, .04, .05)
GAMMAS = (.75, .80, .825, .85, .875, .90, .95)


def old_par(beta, gamma):
    p = OLD_BASE.copy()
    p['tau'] = TAU
    p['beta'] = float(beta)
    p['gamma'] = float(gamma)
    return p


def new_par(beta, gamma):
    o = old_par(beta, gamma)
    return dict(v=o['v'], alpha=o['alpha'], beta=o['beta'], rho=o['rho'],
                rhoT=o['rhoT'], kL=o['k'], kT=o['kT'], tau=o['tau'],
                gamma=o['gamma'])


def local_pairs(beta, gamma):
    P = old_par(beta, gamma)
    out = []
    for xg, dg, sg in central_roots(P, 'G', None):
        try:
            pg = float(pstar(xg, xg, P))
        except Exception:
            continue
        for xb, db, sb in central_roots(P, 'B3', pg):
            if sg < 0 < sb:
                out.append(dict(xg=xg, dg=dg, sg=sg, pg=pg,
                                xb=xb, db=db, sb=sb))
    return out


def global_check(beta, gamma, pair, dense=False):
    par = new_par(beta, gamma)
    xg, xb = pair['xg'], pair['xb']
    pg = private_best_response(xg, xg, par,
                               grid_n=101 if dense else 61,
                               rigorous=True)
    if pg.status != SOLVED_EQUILIBRIUM:
        raise RuntimeError(f'private candidate {pg.status}')
    g = _global_public_br('G', xg, xg, par,
                          x_grid_n=101 if dense else 41,
                          price_grid=41 if dense else 25,
                          rigorous=False)
    b = _global_public_br('B3', xb, xb, par, pbar=pg.p,
                          x_grid_n=121 if dense else 61,
                          price_grid=41 if dense else 25,
                          rigorous=False)
    return pg, g, b


def construction_pass(pair, g, b, tol_gain=5e-4, tol_x=8e-3):
    return (g.gain_over_candidate <= tol_gain
            and b.gain_over_candidate <= tol_gain
            and abs(g.x-pair['xg']) <= tol_x
            and abs(b.x-pair['xb']) <= tol_x)


def main():
    # Analytic route-dominance guard, independent of beta/gamma.
    test = new_par(.01, .85)
    assert test['kL'] + TAU > test['v'] + test['alpha']

    passing = []
    for gamma in GAMMAS:
        for beta in BETAS:
            pairs = local_pairs(beta, gamma)
            if pairs:
                print('BG_LOCAL', 'beta', beta, 'gamma', gamma, 'pairs', len(pairs))
            for pair in pairs:
                print('BG_PAIR', 'beta', beta, 'gamma', gamma,
                      'xG', pair['xg'], 'sG', pair['sg'],
                      'xB', pair['xb'], 'sB', pair['sb'], 'pG_local', pair['pg'])
                try:
                    pg, g, b = global_check(beta, gamma, pair, dense=False)
                except Exception as exc:
                    print('BG_UNRESOLVED', beta, gamma, repr(exc))
                    continue
                print('BG_SCREEN', 'beta', beta, 'gamma', gamma,
                      'pG_allregime', pg.p,
                      'G_BR', g.x, 'G_GAIN', g.gain_over_candidate,
                      'B3_BR', b.x, 'B3_GAIN', b.gain_over_candidate)
                if construction_pass(pair, g, b):
                    passing.append((beta, gamma, pair, pg, g, b))

    if not passing:
        print('BETA_GAMMA_GLOBAL_REGION_SEARCH: NO CONSTRUCTION PASS')
        raise SystemExit(2)

    # Prefer a candidate away from beta=0 and away from the original beta=.05
    # failure, then perform a denser all-regime confirmation.
    passing.sort(key=lambda row: (abs(row[0]-.01), abs(row[1]-.85)))
    beta, gamma, pair, _, _, _ = passing[0]
    pg, g, b = global_check(beta, gamma, pair, dense=True)
    print('BG_DENSE', 'beta', beta, 'gamma', gamma,
          'xG', pair['xg'], 'sG', pair['sg'],
          'xB', pair['xb'], 'sB', pair['sb'],
          'pG', pg.p,
          'G_BR', g.x, 'G_GAIN', g.gain_over_candidate,
          'B3_BR', b.x, 'B3_GAIN', b.gain_over_candidate)

    if not construction_pass(pair, g, b, tol_gain=1e-4, tol_x=5e-3):
        print('BETA_GAMMA_GLOBAL_REGION_SEARCH: DENSE CONFIRMATION FAIL')
        raise SystemExit(3)

    print('SELECTED_BETA', beta)
    print('SELECTED_GAMMA', gamma)
    print('SELECTED_PAIR', pair)
    print('SELECTED_PG', pg)
    print('SELECTED_GLOBAL_G', g)
    print('SELECTED_GLOBAL_B3', b)
    print('BETA_GAMMA_GLOBAL_REGION_SEARCH: CONSTRUCTION PASS')


if __name__ == '__main__':
    main()
