"""Stage 4A v2.1 all-regime search of the repaired sign-reversal witness.

This evaluator reconstructs primitive route choice, partner participation, regional
welfare, private pricing, public finite deviations, and local derivatives.  Stage 11R
repairs the support-side welfare accounting outside the interior region: the gross
participation benefit is not clipped, while participation mass is clipped to [0,1].

Selected existing-primitive vector: beta=.01, gamma=.825, tau=.35.
The public/private searches are deliberately reported as SEARCH EVIDENCE.  Finite
grids, local refinements, and multistart fixed-point checks do not constitute a
certified regret bound or an exact global-equilibrium certificate.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize_scalar

P = dict(v=.1, alpha=.5, beta=.01, rho=.15, rhoT=.05,
         kL=.27, kT=.02, tau=.35, gamma=.825)
XG = .8371022382025995
XB = .8258903860237495
ROUTES = ("H1", "H2", "HT")
SOLVED_EQUILIBRIUM = "SOLVED_EQUILIBRIUM"
MULTIPLE_EQUILIBRIA = "MULTIPLE_EQUILIBRIA"
UNRESOLVED = "UNRESOLVED"
GLOBAL_EVIDENCE_LEVEL = "SEARCH EVIDENCE"


@dataclass
class PriceResult:
    status: str
    p: float | None
    profit: float | None


@dataclass
class BRResult:
    x: float
    welfare: float
    candidate_welfare: float
    gain: float


def upper_envelope(m1, m2, mT, p, par=P):
    q = {
        "0": 0.,
        "H1": par['v'] + par['alpha'] * m1,
        "H2": par['v'] + par['alpha'] * m2,
        "HT": par['v'] + par['alpha'] * mT,
    }
    aggregate = {h: 0. for h in ROUTES}
    regional = []
    for r in (1, 2):
        a = {
            "0": 0.,
            "H1": par['kL'] if r == 1 else par['kL'] + par['tau'],
            "H2": par['kL'] if r == 2 else par['kL'] + par['tau'],
            "HT": par['kT'] + p,
        }
        rr = ("0",) + ROUTES
        cuts = [0., 1.]
        for j, h in enumerate(rr):
            for g in rr[j + 1:]:
                den = q[h] - q[g]
                if abs(den) > 1e-14:
                    z = (a[h] - a[g]) / den
                    if 0 < z < 1:
                        cuts.append(float(z))
        cuts = sorted(set(round(z, 13) for z in cuts))
        shares = {h: 0. for h in rr}
        surplus = 0.
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            if hi - lo <= 1e-13:
                continue
            z = (lo + hi) / 2
            vals = {h: z * q[h] - a[h] for h in rr}
            best = max(vals, key=vals.get)
            shares[best] += hi - lo
            surplus += q[best] * (hi * hi - lo * lo) / 2 - a[best] * (hi - lo)
        for h in ROUTES:
            aggregate[h] += shares[h]
        regional.append(dict(shares=shares, surplus=surplus))
    return np.array([aggregate['H1'], aggregate['H2'], aggregate['HT']]), regional, q


def partner_gross_benefit(n, x1, x2, par=P):
    """Gross support-side participation benefit r_h; this object is never clipped."""
    return np.asarray([
        par['rho'] + x1 + par['beta'] * n[0],
        par['rho'] + x2 + par['beta'] * n[1],
        par['rhoT'] + par['beta'] * n[2],
    ], dtype=float)


def partner_mass(n, x1, x2, par=P):
    """Participation mass m_h = clip(r_h,0,1), retained for compatibility."""
    return np.clip(partner_gross_benefit(n, x1, x2, par), 0., 1.)


def partner_surplus_national(gross, mass=None):
    """Primitive support surplus sum_h int_0^m (r-c) dc."""
    r = np.asarray(gross, dtype=float)
    m = np.clip(r, 0., 1.) if mass is None else np.asarray(mass, dtype=float)
    return float(np.sum(r * m - .5 * m * m))


def partner_surplus_scalar(r):
    m = float(np.clip(r, 0., 1.))
    return r * m - .5 * m * m


def support_surplus_regression():
    # Exact branches implied by c~U[0,1].
    for r in (-.25, 0.):
        assert abs(partner_surplus_scalar(r)) < 1e-14
    for r in (.2, .7):
        assert abs(partner_surplus_scalar(r) - .5 * r * r) < 1e-14
    for r in (1., 1.2, 1.7):
        assert abs(partner_surplus_scalar(r) - (r - .5)) < 1e-14
    # Continuity at both clipping boundaries.
    eps = 1e-8
    assert abs(partner_surplus_scalar(-eps) - partner_surplus_scalar(0.)) < 1e-12
    assert abs(partner_surplus_scalar(eps) - partner_surplus_scalar(0.)) < 1e-12
    assert abs(partner_surplus_scalar(1. - eps) - partner_surplus_scalar(1.)) < 2e-8
    assert abs(partner_surplus_scalar(1. + eps) - partner_surplus_scalar(1.)) < 2e-8
    print('PARTNER_SURPLUS_INTERIOR_SATURATION_BOUNDARIES: PASS')


def fixed_point(x1, x2, p, start, par=P, tol=3e-12, maxit=3000, damping=.7):
    n = np.asarray(start, dtype=float)
    for _ in range(maxit):
        r = partner_gross_benefit(n, x1, x2, par)
        m = np.clip(r, 0., 1.)
        demand, _, _ = upper_envelope(*m, p, par)
        nn = (1 - damping) * n + damping * demand
        if np.linalg.norm(nn - n, np.inf) < tol:
            r = partner_gross_benefit(nn, x1, x2, par)
            m = np.clip(r, 0., 1.)
            demand, regional, q = upper_envelope(*m, p, par)
            residual = float(np.linalg.norm(nn - demand, np.inf))
            if residual > 3e-9:
                raise RuntimeError('UNRESOLVED fixed-point residual')
            return dict(n=nn, r=r, b=m, regional=regional, q=q, residual=residual)
        n = nn
    raise RuntimeError('UNRESOLVED participation continuation')


def multistart_state(x1, x2, p, par=P, full=False):
    starts = [(0, 0, 0), (.5, .5, .5), (1, 1, 1)]
    if full:
        starts += [(1, 0, 0), (0, 1, 0), (0, 0, 2), (2, 2, 0)]
    states = [fixed_point(x1, x2, p, s, par) for s in starts]
    ref = states[0]['n']
    if any(np.linalg.norm(st['n'] - ref, np.inf) > 2e-7 for st in states[1:]):
        raise RuntimeError('MULTIPLE_EQUILIBRIA participation continuation')
    return states[0]


def welfare(i, x1, x2, p, par=P, full=False):
    st = multistart_state(x1, x2, p, par, full)
    project = st['regional'][i - 1]['surplus']
    # Each region receives half of national support surplus.  In the interior m=r,
    # this reduces to sum_h m_h^2/4; outside the interior the primitive integral is used.
    partner = .5 * partner_surplus_national(st['r'], st['b'])
    x = x1 if i == 1 else x2
    return project + partner - par['gamma'] * x * x / 2


def profit(x1, x2, p, par=P, full=False):
    return p * multistart_state(x1, x2, p, par, full)['n'][2]


def private_br(x1, x2, par=P, grid_n=41, full=False):
    upper = par['v'] + par['alpha'] - par['kT']
    grid = np.linspace(0, upper, grid_n)
    vals = np.array([profit(x1, x2, float(p), par, full) for p in grid])
    cand = []
    for j, val in enumerate(vals):
        left = vals[j - 1] if j else -np.inf
        right = vals[j + 1] if j + 1 < len(vals) else -np.inf
        if val >= left and val >= right:
            cand.append((float(grid[j]), float(val)))
            lo = grid[max(0, j - 1)]
            hi = grid[min(len(grid) - 1, j + 1)]
            if hi > lo:
                opt = minimize_scalar(
                    lambda z: -profit(x1, x2, float(z), par, full),
                    bounds=(lo, hi), method='bounded',
                    options={'xatol': 1e-8 if full else 5e-7},
                )
                cand.append((float(opt.x), float(-opt.fun)))
    if not cand:
        return PriceResult(UNRESOLVED, None, None)
    cand.sort(key=lambda z: z[1], reverse=True)
    best = cand[0]
    tied = [z for z in cand if abs(z[1] - best[1]) <= 2e-9]
    distinct = []
    for z in tied:
        if not any(abs(z[0] - d[0]) <= 2e-5 for d in distinct):
            distinct.append(z)
    status = MULTIPLE_EQUILIBRIA if len(distinct) > 1 else SOLVED_EQUILIBRIUM
    return PriceResult(status, best[0], best[1])


def wg(x1, x2, par=P, price_grid=41, full=False):
    pr = private_br(x1, x2, par, price_grid, full)
    if pr.status != SOLVED_EQUILIBRIUM or pr.p is None:
        raise RuntimeError(f'private continuation {pr.status}')
    return welfare(1, x1, x2, pr.p, par, full), pr.p


def wb3(x1, x2, pbar, par=P, full=False):
    return welfare(1, x1, x2, pbar, par, full)


def public_br(mode, rival, candidate, pbar=None, par=P, x_grid_n=51, price_grid=31):
    """Grid + local-refinement search.  This is not a certified global maximum routine."""
    xs = np.linspace(0, 1, x_grid_n)
    vals = []
    for x in xs:
        vals.append(
            wg(float(x), rival, par, price_grid, False)[0]
            if mode == 'G' else wb3(float(x), rival, pbar, par, False)
        )
    vals = np.asarray(vals)
    idx = []
    for j, value in enumerate(vals):
        left = vals[j - 1] if j else -np.inf
        right = vals[j + 1] if j + 1 < len(vals) else -np.inf
        if value >= left and value >= right:
            idx.append(j)
    cand = [(float(xs[j]), float(vals[j])) for j in idx]
    for j in idx:
        lo = xs[max(0, j - 1)]
        hi = xs[min(len(xs) - 1, j + 1)]
        if hi <= lo:
            continue
        if mode == 'G':
            obj = lambda x: -wg(float(x), rival, par, price_grid, False)[0]
        else:
            obj = lambda x: -wb3(float(x), rival, pbar, par, False)
        opt = minimize_scalar(obj, bounds=(lo, hi), method='bounded', options={'xatol': 3e-6})
        cand.append((float(opt.x), float(-opt.fun)))
    if not cand:
        raise RuntimeError('UNRESOLVED public BR search')
    best = max(cand, key=lambda z: z[1])
    cw = (
        wg(candidate, rival, par, price_grid, False)[0]
        if mode == 'G' else wb3(candidate, rival, pbar, par, False)
    )
    return BRResult(best[0], best[1], cw, best[1] - cw)


def local_derivatives(mode, x, pbar=None, h=5e-4):
    f = (
        (lambda a, b: wg(a, b, P, price_grid=51, full=False)[0])
        if mode == 'G' else
        (lambda a, b: wb3(a, b, pbar, P, False))
    )
    f0 = f(x, x)
    xp = f(x + h, x)
    xm = f(x - h, x)
    pp = f(x + h, x + h)
    pm = f(x + h, x - h)
    mp = f(x - h, x + h)
    mm = f(x - h, x - h)
    h11 = (xp - 2 * f0 + xm) / (h * h)
    h12 = (pp - pm - mp + mm) / (4 * h * h)
    return h11, h12, -h12 / h11


def central_first(f, x, h):
    return (f(x + h) - f(x - h)) / (2 * h)


def central_second(f, x, h):
    return (f(x + h) - 2 * f(x) + f(x - h)) / (h * h)


def full_point(x, rival, mode, pbar=None):
    if mode == 'G':
        pr = private_br(x, rival, P, grid_n=81, full=True)
        if pr.status != SOLVED_EQUILIBRIUM:
            raise RuntimeError(f'target private {pr.status}')
        return welfare(1, x, rival, pr.p, P, True), pr.p
    return welfare(1, x, rival, pbar, P, True), pbar


def run():
    support_surplus_regression()
    assert P['kL'] + P['tau'] > P['v'] + P['alpha']
    print('RIVAL_PUBLIC_NONRESIDENT_DOMINANCE: PASS')

    pg = private_br(XG, XG, P, grid_n=81, full=True)
    if pg.status != SOLVED_EQUILIBRIUM or pg.p is None:
        raise RuntimeError(f'on-path private {pg.status}')
    gstate = multistart_state(XG, XG, pg.p, P, full=True)
    print('ON_PATH_PARTICIPATION_RESIDUAL_INF', gstate['residual'])
    print('INDEPENDENT_PG', pg)

    # Local stationarity checks do not assume that the stored candidate is optimal.
    hp = 5e-5
    p_foc = central_first(lambda z: profit(XG, XG, z, P, False), pg.p, hp)
    p_soc = central_second(lambda z: profit(XG, XG, z, P, False), pg.p, hp)
    hx = 7e-4
    g_foc = central_first(lambda z: wg(z, XG, P, price_grid=51, full=False)[0], XG, hx)
    b_foc = central_first(lambda z: wb3(z, XB, pg.p, P, False), XB, hx)
    print('PRIVATE_FOC_SOC', p_foc, p_soc)
    print('PUBLIC_FOC_G_B3', g_foc, b_foc)
    if abs(p_foc) > 2e-5 or p_soc >= 0:
        raise RuntimeError('private stationary candidate failed')
    if abs(g_foc) > 2e-5 or abs(b_foc) > 2e-5:
        raise RuntimeError('public stationary candidate failed')

    g = public_br('G', XG, XG, par=P, x_grid_n=51, price_grid=31)
    b = public_br('B3', XB, XB, pbar=pg.p, par=P, x_grid_n=61, price_grid=31)
    print('INDEPENDENT_SEARCH_G', g)
    print('INDEPENDENT_SEARCH_B3', b)
    if g.gain > 3e-5 or abs(g.x - XG) > 4e-3:
        raise RuntimeError('PROFITABLE FINITE DEVIATION DETECTED IN G SEARCH')
    if b.gain > 3e-5 or abs(b.x - XB) > 4e-3:
        raise RuntimeError('PROFITABLE FINITE DEVIATION DETECTED IN B3 SEARCH')

    # Explicitly include boundaries, the old low-investment region, the neighborhood of
    # support saturation, the detected maxima, and the repaired candidates.
    points = sorted(set([
        0., .05, .10, .18, .20, .35, .50, .75, .84, .85, .90,
        float(b.x), XB, float(g.x), XG, .95, 1.0,
    ]))
    gvals = []
    bvals = []
    for x in points:
        vg, pgx = full_point(x, XG, 'G')
        vb, _ = full_point(x, XB, 'B3', pg.p)
        gvals.append((x, vg))
        bvals.append((x, vb))
        st = multistart_state(x, XG, pgx, P, full=True)
        print('FULL_TARGET', x, vg, vb, 'gross', st['r'].tolist(), 'mass', st['b'].tolist())
    wg_c, _ = full_point(XG, XG, 'G')
    wb_c, _ = full_point(XB, XB, 'B3', pg.p)
    if max(v for _, v in gvals) - wg_c > 3e-6:
        raise RuntimeError('TARGETED SEARCH FOUND PROFITABLE G DEVIATION')
    if max(v for _, v in bvals) - wb_c > 3e-6:
        raise RuntimeError('TARGETED SEARCH FOUND PROFITABLE B3 DEVIATION')
    print('TARGETED_FULL_MULTISTART_SEARCH: PASS')

    for h in (7e-4, 5e-4, 3.5e-4):
        dg = local_derivatives('G', XG, h=h)
        db = local_derivatives('B3', XB, pg.p, h=h)
        print('LOCAL_SIGNS', h, 'G', dg, 'B3', db)
        if not (dg[0] < 0 and dg[1] < 0 and dg[2] < 0):
            raise RuntimeError('G local substitute sign not robust')
        if not (db[0] < 0 and db[1] > 0 and db[2] > 0):
            raise RuntimeError('B3 local complement sign not robust')

    print('GLOBAL_EVIDENCE_LEVEL:', GLOBAL_EVIDENCE_LEVEL)
    print('NO_CERTIFIED_REGRET_BOUND: TRUE')
    print('STAGE4A_REPAIRED_ALL_REGIME_SEARCH: PASS')


if __name__ == '__main__':
    run()
