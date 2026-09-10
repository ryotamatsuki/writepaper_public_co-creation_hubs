"""Independent Stage-11 regression for the v2.1 public-hub manuscript.

This reviewer-side implementation intentionally does not import the Stage-4A audit,
production result generator, or analytic derivation scripts. It reconstructs the
beta-zero identities and the repaired witness directly from model primitives.

The numerical part is an adversarial regression, not a uniqueness proof. It checks
candidate deviations, large off-path histories, participation continuation from ten
starts, a sampled alternative-equilibrium search, welfare accounting, and local signs.
"""
from __future__ import annotations

import math
import numpy as np
import sympy as sp
from scipy.optimize import minimize_scalar

P = dict(v=.1, alpha=.5, beta=.01, rho=.15, rhoT=.05,
         kL=.27, kT=.02, tau=.35, gamma=.825)
XG = .8371022382025995
XB = .8258903860237495
EXPECTED_P = .018403679612460814
EXPECTED_WG = .5967883276987462
EXPECTED_WB = .5859284219438758
EXPECTED_WEDGE = .49044674754794976
P_MAX = P['v'] + P['alpha'] - P['kT']
STARTS = [
    (0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2), (2, 2, 0),
    (2, 0, 2), (0, 2, 2), (2, 2, 2), (.5, .5, .5), (1, 1, 1),
]


def verify_symbolic() -> None:
    # B3: beta-zero welfare is separable in D1 and D2, so the cross effect is zero.
    D1, D2, T, a, d, k, alpha, eps, rhoT = sp.symbols(
        'D1 D2 T a d k alpha eps rhoT', positive=True)
    t10, t20, s0 = d / D1, d / D2, a / T
    b10, b20, bT0 = rhoT + D1 / alpha, rhoT + D2 / alpha, rhoT
    W0 = (T * (t10**2 - s0**2) / 2 - a * (t10 - s0)
          + (T + D1) * (1 - t10**2) / 2 - k * (1 - t10)
          + (b10**2 + b20**2 + bT0**2) / 4)
    assert sp.simplify(alpha**2 * sp.diff(W0, D1, D2)) == 0

    # First-order continuation in delta=alpha*beta, derived from the primitive fixed point.
    Q = t10 + t20 - 2 * a / T
    dt1 = -t10 * ((1 - t10) - Q) / D1
    t1 = t10 + eps * dt1
    s = s0 - eps * a * Q / T**2
    qT = T + eps * Q
    q1 = T + D1 + eps * (1 - t10)
    b1 = b10 + eps * (1 - t10) / alpha
    b2 = b20 + eps * (1 - t20) / alpha
    bT = bT0 + eps * Q / alpha
    W = (qT * (t1**2 - s**2) / 2 - a * (t1 - s)
         + q1 * (1 - t1**2) / 2 - k * (1 - t1)
         + (b1**2 + b2**2 + bT**2) / 4)
    W_delta = sp.diff(W, eps).subs(eps, 0)
    cross_delta = sp.factor(alpha**2 * sp.diff(W_delta, D1, D2).subs(k, a + d))
    assert sp.simplify(cross_delta - alpha**2 * d**3 / (D1**3 * D2**2)) == 0

    # G: independently substitute the beta-zero private optimum into regional project surplus
    # and differentiate the reduced welfare directly, rather than using the production chain-rule code.
    D1g, D2g, Tg, kg, kTg, ag = sp.symbols('D1g D2g Tg kg kTg ag', positive=True)
    p = sp.symbols('p', real=True)
    aa = kTg + p
    dd = kg - aa
    tt1, tt2, ss = dd / D1g, dd / D2g, aa / Tg
    demand = tt1 + tt2 - 2 * ss
    pstar = sp.solve(sp.Eq(sp.diff(p * demand, p), 0), p)[0]
    PS1 = (Tg * (tt1**2 - ss**2) / 2 - aa * (tt1 - ss)
           + (Tg + D1g) * (1 - tt1**2) / 2 - kg * (1 - tt1))
    reduced = sp.simplify(PS1.subs(p, pstar))
    D = sp.symbols('D', positive=True)
    cross_g = sp.factor(ag**2 * sp.diff(reduced, D1g, D2g).subs({D1g: D, D2g: D}))
    target = -3 * Tg * ag**2 * kg**2 / (16 * D * (D + Tg)**3)
    assert sp.simplify(cross_g - target) == 0
    print('INDEPENDENT_SYMBOLIC_T1_T2: PASS')


def allocation(b: np.ndarray, price: float):
    q = np.array([0., P['v'] + P['alpha'] * b[0],
                  P['v'] + P['alpha'] * b[1], P['v'] + P['alpha'] * b[2]])
    aggregate = np.zeros(3)
    regions = []
    for r in (1, 2):
        costs = np.array([
            0.,
            P['kL'] if r == 1 else P['kL'] + P['tau'],
            P['kL'] if r == 2 else P['kL'] + P['tau'],
            P['kT'] + price,
        ])
        cuts = [0., 1.]
        for h in range(4):
            for g in range(h + 1, 4):
                den = q[h] - q[g]
                if abs(den) > 1e-13:
                    z = (costs[h] - costs[g]) / den
                    if 0 < z < 1:
                        cuts.append(float(z))
        cuts = sorted(set(round(z, 12) for z in cuts))
        shares = np.zeros(4)
        surplus = 0.
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            if hi - lo <= 1e-12:
                continue
            mid = (lo + hi) / 2
            route = int(np.argmax(mid * q - costs))
            shares[route] += hi - lo
            surplus += q[route] * (hi * hi - lo * lo) / 2 - costs[route] * (hi - lo)
        aggregate += shares[1:]
        regions.append((shares, surplus))
    return aggregate, regions, q


def solve_state(x1: float, x2: float, price: float, start=(.5, .5, .5),
                damping=.7, tol=2e-11, maxit=1500):
    n = np.array(start, dtype=float)
    for _ in range(maxit):
        b = np.clip([
            P['rho'] + x1 + P['beta'] * n[0],
            P['rho'] + x2 + P['beta'] * n[1],
            P['rhoT'] + P['beta'] * n[2],
        ], 0, 1)
        demand, regions, q = allocation(b, price)
        nn = (1 - damping) * n + damping * demand
        if np.max(np.abs(nn - n)) < tol:
            b = np.clip([
                P['rho'] + x1 + P['beta'] * nn[0],
                P['rho'] + x2 + P['beta'] * nn[1],
                P['rhoT'] + P['beta'] * nn[2],
            ], 0, 1)
            demand, regions, q = allocation(b, price)
            residual = float(np.max(np.abs(demand - nn)))
            if residual > 1e-8:
                raise RuntimeError(f'UNRESOLVED residual {residual}')
            return nn, b, regions, q
        n = nn
    raise RuntimeError('UNRESOLVED participation continuation')


def multistart(x1: float, x2: float, price: float):
    states = [solve_state(x1, x2, price, s) for s in STARTS]
    ns = np.array([s[0] for s in states])
    spread = float(np.max(np.ptp(ns, axis=0)))
    if spread > 1e-6:
        raise RuntimeError(f'MULTIPLE_EQUILIBRIA participation spread {spread}')
    return states[0], spread


def welfare(i: int, x1: float, x2: float, price: float, many_starts=False):
    st = multistart(x1, x2, price)[0] if many_starts else solve_state(x1, x2, price)
    n, b, regions, _ = st
    x = x1 if i == 1 else x2
    return regions[i - 1][1] + .25 * float(np.sum(b * b)) - P['gamma'] * x * x / 2


def profit(x1: float, x2: float, price: float, many_starts=False):
    st = multistart(x1, x2, price)[0] if many_starts else solve_state(x1, x2, price)
    return price * st[0][2]


def private_br(x1: float, x2: float, grid_n=41):
    # Any price above P_MAX gives zero private demand because q_T <= v+alpha.
    grid = np.linspace(0, P_MAX, grid_n)
    vals = np.array([profit(x1, x2, float(p)) for p in grid])
    j = int(np.argmax(vals))
    lo, hi = grid[max(0, j - 1)], grid[min(grid_n - 1, j + 1)]
    if hi == lo:
        return float(grid[j]), float(vals[j])
    opt = minimize_scalar(lambda z: -profit(x1, x2, float(z)), bounds=(lo, hi),
                          method='bounded', options={'xatol': 2e-8, 'maxiter': 80})
    return float(opt.x), float(-opt.fun)


def wg(x: float, rival: float, pgrid=31):
    price, _ = private_br(x, rival, pgrid)
    return welfare(1, x, rival, price), price


def wb3(x: float, rival: float, pbar: float):
    return welfare(1, x, rival, pbar)


def public_br(mode: str, rival: float, pbar=None, xgrid=31, pgrid=31):
    xs = np.linspace(0, 1, xgrid)
    vals = np.array([
        wg(float(x), rival, pgrid)[0] if mode == 'G' else wb3(float(x), rival, pbar)
        for x in xs
    ])
    local = []
    for j, val in enumerate(vals):
        left = vals[j - 1] if j else -np.inf
        right = vals[j + 1] if j + 1 < xgrid else -np.inf
        if val >= left and val >= right:
            local.append(j)
    cand = [(float(xs[j]), float(vals[j])) for j in local]
    for j in local:
        lo, hi = xs[max(0, j - 1)], xs[min(xgrid - 1, j + 1)]
        if hi <= lo:
            continue
        if mode == 'G':
            obj = lambda z: -wg(float(z), rival, pgrid)[0]
        else:
            obj = lambda z: -wb3(float(z), rival, pbar)
        opt = minimize_scalar(obj, bounds=(lo, hi), method='bounded',
                              options={'xatol': 3e-6, 'maxiter': 60})
        cand.append((float(opt.x), float(-opt.fun)))
    if not cand:
        raise RuntimeError('UNRESOLVED public best response')
    return max(cand, key=lambda z: z[1])


def numerical_attack() -> None:
    assert P['kL'] + P['tau'] > P['v'] + P['alpha']
    # Strict route dominance: z*q_remote - (kL+tau) <= v+alpha-(kL+tau) < 0.
    assert P['v'] + P['alpha'] - (P['kL'] + P['tau']) < 0
    print('REMOTE_PUBLIC_STRICT_DOMINANCE: PASS')

    pg, _ = private_br(XG, XG, 61)
    _, spread = multistart(XG, XG, pg)
    assert spread < 1e-6 and abs(pg - EXPECTED_P) < 3e-6

    gbr = public_br('G', XG, xgrid=31, pgrid=31)
    bbr = public_br('B3', XB, pbar=pg, xgrid=41, pgrid=31)
    cand_g = wg(XG, XG, 41)[0]
    cand_b = wb3(XB, XB, pg)
    assert abs(gbr[0] - XG) < 4e-3 and gbr[1] - cand_g < 3e-5
    assert abs(bbr[0] - XB) < 4e-3 and bbr[1] - cand_b < 3e-5
    print('INDEPENDENT_CANDIDATE_DEVIATION_G:', gbr, 'gain', gbr[1] - cand_g)
    print('INDEPENDENT_CANDIDATE_DEVIATION_B3:', bbr, 'gain', bbr[1] - cand_b)

    # Large deviations intentionally exit the regular branch. Reoptimize the private price
    # and demand, then verify continuation from ten dispersed starts.
    adversarial = [0., .05, .10, .18, .20, .35, .50, .75, .95, 1.0]
    max_spread = 0.
    for x in adversarial:
        price, _ = private_br(x, XG, 31)
        _, s = multistart(x, XG, price)
        max_spread = max(max_spread, s)
    assert max_spread < 1e-6
    print('OFF_PATH_10_START_CONTINUATIONS: PASS max_spread', max_spread)

    # Sample the global best-response correspondence outside the reported equilibrium.
    # This is an alternative-equilibrium attack, not a uniqueness certificate.
    rivals_g = [0., .25, .50, .75, XG, 1.0]
    rivals_b = [0., .25, .50, .75, XB, 1.0]
    br_g = [(r, public_br('G', r, xgrid=21, pgrid=21)[0]) for r in rivals_g]
    br_b = [(r, public_br('B3', r, pbar=pg, xgrid=25, pgrid=21)[0]) for r in rivals_b]
    # No sampled fixed point away from the reported neighborhood should appear.
    assert all(abs(br - r) > .02 for r, br in br_g if abs(r - XG) > .05)
    assert all(abs(br - r) > .02 for r, br in br_b if abs(r - XB) > .05)
    print('ALTERNATIVE_EQUILIBRIUM_SAMPLED_ATTACK_G:', br_g)
    print('ALTERNATIVE_EQUILIBRIUM_SAMPLED_ATTACK_B3:', br_b)

    # National welfare and its local directional decomposition.
    def comps(x1, x2):
        price, _ = private_br(x1, x2, 51)
        return np.array([welfare(1, x1, x2, price), welfare(2, x1, x2, price),
                         profit(x1, x2, price)]), price

    cg, pg2 = comps(XG, XG)
    wn_g = float(cg.sum())
    wn_b = float(welfare(1, XB, XB, pg) + welfare(2, XB, XB, pg) + profit(XB, XB, pg))
    assert abs(wn_g - EXPECTED_WG) < 2e-5
    assert abs(wn_b - EXPECTED_WB) < 2e-5
    assert wn_g > wn_b
    h = 1e-3
    cp, _ = comps(XG + h, XG)
    cm, _ = comps(XG - h, XG)
    deriv = (cp - cm) / (2 * h)
    assert abs(float(deriv.sum()) - EXPECTED_WEDGE) < 2e-4
    assert deriv[1] > 0 and deriv[2] < 0
    print('INDEPENDENT_WELFARE:', wn_g, wn_b, wn_g - wn_b)
    print('INDEPENDENT_WEDGE_COMPONENTS:', deriv.tolist(), 'sum', float(deriv.sum()))

    # Local BR slope signs, estimated from reduced welfare Hessians.
    def slope(mode, x, pbar=None, h=7e-4):
        if mode == 'G':
            f = lambda a, b: wg(a, b, 41)[0]
        else:
            f = lambda a, b: wb3(a, b, pbar)
        f0 = f(x, x)
        h11 = (f(x + h, x) - 2 * f0 + f(x - h, x)) / h**2
        h12 = (f(x + h, x + h) - f(x + h, x - h)
               - f(x - h, x + h) + f(x - h, x - h)) / (4 * h**2)
        return h11, h12, -h12 / h11

    dg = slope('G', XG)
    db = slope('B3', XB, pg)
    assert dg[0] < 0 and dg[1] < 0 and dg[2] < 0
    assert db[0] < 0 and db[1] > 0 and db[2] > 0
    print('INDEPENDENT_LOCAL_SIGNS G', dg, 'B3', db)
    print('STAGE11_INDEPENDENT_REGRESSION: PASS')


if __name__ == '__main__':
    verify_symbolic()
    numerical_attack()
