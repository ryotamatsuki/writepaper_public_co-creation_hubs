"""v2.0 Stage 4A clean-room adversarial audit.

This file deliberately does NOT import the production central-regime solver.
It reconstructs project route choice from the upper envelope of primitive
utilities, clips partner participation to the primitive c~U[0,1] support,
solves the induced participation fixed point directly, and globally scans the
private price over the economically relevant compact interval.

It also re-derives the two beta-zero analytic identities with SymPy.
The purpose is falsification, not reproduction of the production path.
"""
from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.optimize import minimize_scalar

P = dict(v=.1, alpha=.5, beta=.05, rho=.15, rhoT=.05,
         kL=.27, kT=.02, tau=.05, gamma=.9)
XG = 0.684028196361275
XB = 0.656020390747393
PG = 0.02274665145
XDEV = 0.1875
ROUTES = ("H1", "H2", "HT")


def _upper_envelope(b1, b2, bT, p, par=P):
    q = {"0": 0.0,
         "H1": par["v"] + par["alpha"] * b1,
         "H2": par["v"] + par["alpha"] * b2,
         "HT": par["v"] + par["alpha"] * bT}
    aggregate = {h: 0.0 for h in ROUTES}
    regional = []
    for r in (1, 2):
        a = {"0": 0.0,
             "H1": par["kL"] if r == 1 else par["kL"] + par["tau"],
             "H2": par["kL"] if r == 2 else par["kL"] + par["tau"],
             "HT": par["kT"] + p}
        rr = ("0",) + ROUTES
        cuts = [0.0, 1.0]
        for j, h in enumerate(rr):
            for g in rr[j + 1:]:
                den = q[h] - q[g]
                if abs(den) > 1e-14:
                    z = (a[h] - a[g]) / den
                    if 0.0 < z < 1.0:
                        cuts.append(float(z))
        cuts = sorted(set(round(z, 13) for z in cuts))
        shares = {h: 0.0 for h in rr}
        surplus = 0.0
        envelope = []
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            if hi - lo <= 1e-13:
                continue
            z = (lo + hi) / 2.0
            vals = {h: z * q[h] - a[h] for h in rr}
            h = max(vals, key=vals.get)
            shares[h] += hi - lo
            surplus += q[h] * (hi * hi - lo * lo) / 2.0 - a[h] * (hi - lo)
            envelope.append((lo, hi, h))
        for h in ROUTES:
            aggregate[h] += shares[h]
        regional.append(dict(shares=shares, surplus=surplus, envelope=envelope))
    return np.array([aggregate["H1"], aggregate["H2"], aggregate["HT"]]), regional, q


def _partner_mass(n, x1, x2, par=P):
    # This is the primitive c~U[0,1] participation rule, including clipping.
    return np.clip([par["rho"] + x1 + par["beta"] * n[0],
                    par["rho"] + x2 + par["beta"] * n[1],
                    par["rhoT"] + par["beta"] * n[2]], 0.0, 1.0)


def fixed_point(x1, x2, p, start=(.4, .4, .5), par=P,
                tol=1e-12, maxit=5000, damping=.7):
    n = np.asarray(start, dtype=float)
    for _ in range(maxit):
        b = _partner_mass(n, x1, x2, par)
        demand, _, _ = _upper_envelope(*b, p, par)
        nn = (1.0 - damping) * n + damping * demand
        if np.linalg.norm(nn - n, ord=np.inf) < tol:
            b = _partner_mass(nn, x1, x2, par)
            demand, regional, q = _upper_envelope(*b, p, par)
            if np.linalg.norm(nn - demand, ord=np.inf) > 1e-9:
                raise RuntimeError("fixed-point residual too large")
            return dict(n=nn, b=b, regional=regional, q=q)
        n = nn
    raise RuntimeError("UNRESOLVED participation continuation")


def multistart_state(x1, x2, p, par=P):
    starts = ((.4, .4, .5), (0., 0., 0.), (1., 1., 0.),
              (0., 0., 2.), (1., 1., 1.))
    states = [fixed_point(x1, x2, p, s, par) for s in starts]
    n0 = states[0]["n"]
    if any(np.linalg.norm(s["n"] - n0, ord=np.inf) > 1e-7 for s in states[1:]):
        raise RuntimeError("multiple continuation candidates detected")
    return states[0]


def profit(x1, x2, p, par=P):
    return p * multistart_state(x1, x2, p, par)["n"][2]


def regional_welfare(i, x1, x2, p, par=P):
    st = multistart_state(x1, x2, p, par)
    project_surplus = st["regional"][i - 1]["surplus"]
    partner_surplus = .25 * float(np.sum(st["b"] ** 2))
    x = x1 if i == 1 else x2
    return project_surplus + partner_surplus - par["gamma"] * x * x / 2.0


def global_private_price(x1, x2, par=P):
    # Since b_T<=1, q_T<=v+alpha=.6. For p>=.58, kT+p>=.6,
    # hence U_T(z)<=0 for all z<=1 and private demand is zero. Therefore
    # every positive-profit global optimum lies in [0,.58].
    upper = par["v"] + par["alpha"] - par["kT"]
    grid = np.linspace(0.0, upper, 201)
    vals = np.array([profit(x1, x2, p, par) for p in grid])
    candidates = []
    for j, value in enumerate(vals):
        left = vals[j - 1] if j else -np.inf
        right = vals[j + 1] if j + 1 < len(vals) else -np.inf
        if value >= left and value >= right:
            lo = grid[max(0, j - 1)]
            hi = grid[min(len(grid) - 1, j + 1)]
            if hi > lo:
                opt = minimize_scalar(lambda z: -profit(x1, x2, float(z), par),
                                      bounds=(lo, hi), method="bounded",
                                      options={"xatol": 1e-10})
                candidates.append((float(opt.x), float(-opt.fun)))
            candidates.append((float(grid[j]), float(value)))
    if not candidates:
        raise RuntimeError("UNRESOLVED private-price continuation")
    return max(candidates, key=lambda item: item[1])


def symbolic_local_checks():
    # Independent beta-zero derivation of the full-game cross derivative.
    alpha, T, kL, kT, D1, D2 = sp.symbols(
        "alpha T kL kT D1 D2", positive=True)
    A = 1 / D1 + 1 / D2
    B = 2 / T
    p = ((kL - kT) * A - kT * B) / (2 * (A + B))
    d = kL - kT - p
    s = (kT + p) / T
    t1 = d / D1
    q1 = T + D1
    ps1 = (T * (t1**2 - s**2) / 2
           - (kT + p) * (t1 - s)
           + q1 * (1 - t1**2) / 2 - kL * (1 - t1))
    mg = sp.factor(alpha**2 * sp.diff(sp.diff(ps1, D1), D2).subs(D2, D1))
    expected_g = -3 * T * alpha**2 * kL**2 / (16 * D1 * (D1 + T)**3)
    assert sp.simplify(mg - expected_g) == 0

    # Clean-room first-order expansion of B3 around delta=alpha*beta=0.
    d0, a, v = sp.symbols("d0 a v", positive=True)
    e = sp.symbols("e")
    t10, t20 = d0 / D1, d0 / D2
    uT = t10 + t20 - 2 * a / T
    u1 = sp.simplify(-t10 * ((1 - t10) - uT) / D1)
    u2 = sp.simplify(-t20 * ((1 - t20) - uT) / D2)
    t1e, t2e, qTe = t10 + e*u1, t20 + e*u2, T + e*uT
    q1e, q2e = T + D1 + e*(1 - t10), T + D2 + e*(1 - t20)
    se = a / qTe
    ps = (qTe * (t1e**2 - se**2) / 2 - a * (t1e - se)
          + q1e * (1 - t1e**2) / 2 - (a + d0) * (1 - t1e))
    partner = (((q1e-v)/alpha)**2 + ((q2e-v)/alpha)**2
               + ((qTe-v)/alpha)**2) / 4
    coeff = sp.simplify(sp.diff(ps + partner, e).subs(e, 0))
    dM_ddelta = sp.factor(alpha**2 * sp.diff(sp.diff(coeff, D1), D2))
    expected_b3_delta = alpha**2 * d0**3 / (D1**3 * D2**2)
    assert sp.simplify(dM_ddelta - expected_b3_delta) == 0
    return str(mg), str(dM_ddelta)


def run():
    mg, mb = symbolic_local_checks()

    # G canonical stationary configuration versus an off-regime finite deviation.
    p_cand, _ = global_private_price(XG, XG)
    w_g_cand = regional_welfare(1, XG, XG, p_cand)
    p_dev, _ = global_private_price(XDEV, XG)
    w_g_dev = regional_welfare(1, XDEV, XG, p_dev)
    st_g_dev = multistart_state(XDEV, XG, p_dev)

    # B3 uses the matched canonical private price but must still allow public
    # finite deviations over x_i in [0,1].
    w_b_cand = regional_welfare(1, XB, XB, PG)
    w_b_dev = regional_welfare(1, XDEV, XB, PG)
    st_b_dev = multistart_state(XDEV, XB, PG)

    # Exact primitive observation behind the free-riding deviation:
    # with n_1^F=0 at x=0.1875, b1=rho+x=.3375 and
    # q1=v+alpha*b1=.26875<kL=.27, so H1 is below the outside option for all z<=1.
    q1_inactive = P["v"] + P["alpha"] * (P["rho"] + XDEV)
    assert q1_inactive < P["kL"]
    assert st_g_dev["n"][0] < 1e-8 and st_b_dev["n"][0] < 1e-8
    assert w_g_dev - w_g_cand > 0.02
    assert w_b_dev - w_b_cand > 0.02

    print("INDEPENDENT LOCAL IDENTITIES: PASS")
    print("M_G(beta=0) =", mg)
    print("d M_B3 / d delta at zero =", mb)
    print("G candidate p, W1 =", p_cand, w_g_cand)
    print("G deviation x1, p, W1, gain =", XDEV, p_dev, w_g_dev, w_g_dev-w_g_cand)
    print("G deviation masses H1,H2,HT =", st_g_dev["n"])
    print("B3 candidate W1 =", w_b_cand)
    print("B3 deviation x1, W1, gain =", XDEV, w_b_dev, w_b_dev-w_b_cand)
    print("B3 deviation masses H1,H2,HT =", st_b_dev["n"])
    print("GLOBAL-EQUILIBRIUM ADVERSARIAL RESULT: COUNTEREXAMPLE FOUND")


if __name__ == "__main__":
    run()
