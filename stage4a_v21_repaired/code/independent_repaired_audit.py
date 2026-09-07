"""Stage 4A v2.1 clean-room audit of the repaired global sign-reversal witness.

This file deliberately does NOT import any Stage-4 repair solver, the old
central-regime solver, or the old exact verifier. It reconstructs primitive route
choice, partner participation, regional welfare, private pricing, public finite
deviations, and local derivatives independently.

Selected construction witness from Stage 4R:
    beta=.01, gamma=.825, tau=.35
with all other primitives unchanged.

The script is adversarial/fail-closed: unresolved participation, economically
material multiple continuations, or profitable finite public deviations raise.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize_scalar

P = dict(v=.1, alpha=.5, beta=.01, rho=.15, rhoT=.05,
         kL=.27, kT=.02, tau=.35, gamma=.825)
XG = 0.8371022382025995
XB = 0.8258903860237495
ROUTES = ("H1", "H2", "HT")

SOLVED_EQUILIBRIUM = "SOLVED_EQUILIBRIUM"
MULTIPLE_EQUILIBRIA = "MULTIPLE_EQUILIBRIA"
UNRESOLVED = "UNRESOLVED"


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


def upper_envelope(b1, b2, bT, p, par=P):
    q = {"0": 0.0,
         "H1": par["v"] + par["alpha"] * b1,
         "H2": par["v"] + par["alpha"] * b2,
         "HT": par["v"] + par["alpha"] * bT}
    aggregate = {h: 0.0 for h in ROUTES}
    regional = []
    for r in (1, 2):
        access = {"0": 0.0,
                  "H1": par["kL"] if r == 1 else par["kL"] + par["tau"],
                  "H2": par["kL"] if r == 2 else par["kL"] + par["tau"],
                  "HT": par["kT"] + p}
        rr = ("0",) + ROUTES
        cuts = [0.0, 1.0]
        for j, h in enumerate(rr):
            for g in rr[j+1:]:
                den = q[h] - q[g]
                if abs(den) > 1e-14:
                    z = (access[h] - access[g]) / den
                    if 0.0 < z < 1.0:
                        cuts.append(float(z))
        cuts = sorted(set(round(z, 13) for z in cuts))
        shares = {h: 0.0 for h in rr}
        surplus = 0.0
        envelope = []
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            if hi-lo <= 1e-13:
                continue
            zmid = (lo+hi)/2.0
            vals = {h: zmid*q[h]-access[h] for h in rr}
            best = max(vals, key=vals.get)
            shares[best] += hi-lo
            surplus += (q[best]*(hi*hi-lo*lo)/2.0
                        - access[best]*(hi-lo))
            envelope.append((lo, hi, best))
        for h in ROUTES:
            aggregate[h] += shares[h]
        regional.append(dict(shares=shares, surplus=surplus, envelope=envelope))
    return np.array([aggregate["H1"], aggregate["H2"], aggregate["HT"]]), regional, q


def partner_mass(n, x1, x2, par=P):
    return np.clip([par["rho"] + x1 + par["beta"]*n[0],
                    par["rho"] + x2 + par["beta"]*n[1],
                    par["rhoT"] + par["beta"]*n[2]], 0.0, 1.0)


def fixed_point(x1, x2, p, start, par=P, tol=2e-12, maxit=4000, damping=.7):
    n = np.asarray(start, dtype=float)
    for _ in range(maxit):
        b = partner_mass(n, x1, x2, par)
        demand, _, _ = upper_envelope(*b, p, par)
        nn = (1.0-damping)*n + damping*demand
        if np.linalg.norm(nn-n, ord=np.inf) < tol:
            b = partner_mass(nn, x1, x2, par)
            demand, regional, q = upper_envelope(*b, p, par)
            if np.linalg.norm(nn-demand, ord=np.inf) > 2e-9:
                raise RuntimeError("UNRESOLVED fixed-point residual")
            return dict(n=nn, b=b, regional=regional, q=q)
        n = nn
    raise RuntimeError("UNRESOLVED participation continuation")


def multistart_state(x1, x2, p, par=P, full=False):
    starts = [(0.,0.,0.), (.5,.5,.5), (1.,1.,1.)]
    if full:
        starts += [(1.,0.,0.), (0.,1.,0.), (0.,0.,2.), (2.,2.,0.)]
    states = [fixed_point(x1, x2, p, s, par) for s in starts]
    ref = states[0]["n"]
    for st in states[1:]:
        if np.linalg.norm(st["n"]-ref, ord=np.inf) > 2e-7:
            raise RuntimeError("MULTIPLE_EQUILIBRIA participation continuation")
    return states[0]


def welfare(i, x1, x2, p, par=P, full=False):
    st = multistart_state(x1, x2, p, par, full=full)
    project = st["regional"][i-1]["surplus"]
    # Partner welfare follows the primitive threshold-surplus integral used by the model.
    partner = .25 * float(np.sum(st["b"]**2))
    x = x1 if i == 1 else x2
    return project + partner - par["gamma"]*x*x/2.0


def profit(x1, x2, p, par=P, full=False):
    return p * multistart_state(x1, x2, p, par, full=full)["n"][2]


def private_br(x1, x2, par=P, grid_n=61, full=False):
    # bT<=1 => qT<=v+alpha=.6. At p>=.58, kT+p>=.6 and HT cannot
    # strictly beat the outside option for any z<=1, so positive profit is impossible.
    upper = par["v"] + par["alpha"] - par["kT"]
    grid = np.linspace(0.0, upper, grid_n)
    vals = np.array([profit(x1,x2,float(p),par,full=full) for p in grid])
    candidates = []
    for j, val in enumerate(vals):
        left = vals[j-1] if j else -np.inf
        right = vals[j+1] if j+1<len(vals) else -np.inf
        if val >= left and val >= right:
            candidates.append((float(grid[j]), float(val)))
            lo = grid[max(0,j-1)]; hi = grid[min(len(grid)-1,j+1)]
            if hi > lo:
                opt = minimize_scalar(lambda z: -profit(x1,x2,float(z),par,full=full),
                                      bounds=(lo,hi), method="bounded",
                                      options={"xatol": 5e-9 if full else 5e-7})
                candidates.append((float(opt.x), float(-opt.fun)))
    if not candidates:
        return PriceResult(UNRESOLVED, None, None)
    candidates.sort(key=lambda z:z[1], reverse=True)
    best = candidates[0]
    tied = [c for c in candidates if abs(c[1]-best[1]) <= 2e-9]
    distinct = []
    for c in tied:
        if not any(abs(c[0]-d[0]) <= 2e-5 for d in distinct):
            distinct.append(c)
    status = MULTIPLE_EQUILIBRIA if len(distinct)>1 else SOLVED_EQUILIBRIUM
    return PriceResult(status, best[0], best[1])


def wg(x1, x2, par=P, price_grid=61, full=False):
    pr = private_br(x1,x2,par,grid_n=price_grid,full=full)
    if pr.status != SOLVED_EQUILIBRIUM or pr.p is None:
        raise RuntimeError(f"private continuation {pr.status}")
    return welfare(1,x1,x2,pr.p,par,full=full), pr.p


def wb3(x1, x2, pbar, par=P, full=False):
    return welfare(1,x1,x2,pbar,par,full=full)


def public_br(mode, rival, candidate, pbar=None, par=P,
              x_grid_n=81, price_grid=41, full=False):
    xs = np.linspace(0.,1.,x_grid_n)
    vals = []
    for x in xs:
        if mode == "G":
            val = wg(float(x),rival,par,price_grid=price_grid,full=full)[0]
        else:
            val = wb3(float(x),rival,pbar,par,full=full)
        vals.append(val)
    vals = np.asarray(vals)
    idx = []
    for j,v in enumerate(vals):
        l = vals[j-1] if j else -np.inf
        r = vals[j+1] if j+1<len(vals) else -np.inf
        if v >= l and v >= r:
            idx.append(j)
    candidates = [(float(xs[j]),float(vals[j])) for j in idx]
    for j in idx:
        lo=xs[max(0,j-1)]; hi=xs[min(len(xs)-1,j+1)]
        if hi<=lo: continue
        def obj(x):
            if mode == "G":
                return -wg(float(x),rival,par,price_grid=price_grid,full=full)[0]
            return -wb3(float(x),rival,pbar,par,full=full)
        opt=minimize_scalar(obj,bounds=(lo,hi),method="bounded",
                            options={"xatol":2e-6 if full else 1e-5})
        candidates.append((float(opt.x),float(-opt.fun)))
    if not candidates:
        raise RuntimeError("UNRESOLVED public BR")
    best=max(candidates,key=lambda z:z[1])
    if mode=="G": cand=wg(candidate,rival,par,price_grid=price_grid,full=full)[0]
    else: cand=wb3(candidate,rival,pbar,par,full=full)
    return BRResult(best[0],best[1],cand,best[1]-cand)


def local_derivatives(mode, x, pbar=None, h=4e-4):
    if mode=="G":
        f=lambda a,b: wg(a,b,P,price_grid=81,full=False)[0]
    else:
        f=lambda a,b: wb3(a,b,pbar,P,full=False)
    f0=f(x,x)
    xp=f(x+h,x); xm=f(x-h,x)
    pp=f(x+h,x+h); pm=f(x+h,x-h); mp=f(x-h,x+h); mm=f(x-h,x-h)
    h11=(xp-2*f0+xm)/(h*h)
    h12=(pp-pm-mp+mm)/(4*h*h)
    slope=-h12/h11
    return h11,h12,slope


def targeted_full_multistart_checks(pG):
    # Histories chosen to attack boundaries, the old free-riding region, central
    # candidates, and high-investment corners. Each price continuation is globally
    # re-optimized with the fuller start set.
    points=(0.0,.10,.18,.20,.35,.50,XB,XG,.95,1.0)
    for x in points:
        pr=private_br(x,XG,P,grid_n=81,full=True)
        if pr.status!=SOLVED_EQUILIBRIUM:
            raise RuntimeError(f"target x={x}: private {pr.status}")
        _=welfare(1,x,XG,pr.p,P,full=True)
        _=welfare(1,x,XB,pG,P,full=True)
    print("TARGETED_FULL_MULTISTART: PASS")


def run():
    # Primitive analytic dominance of rival public route for nonresidents.
    assert P["kL"]+P["tau"] > P["v"]+P["alpha"]
    print("RIVAL_PUBLIC_NONRESIDENT_DOMINANCE: PASS",
          P["kL"]+P["tau"], ">", P["v"]+P["alpha"])

    pg=private_br(XG,XG,P,grid_n=121,full=True)
    if pg.status!=SOLVED_EQUILIBRIUM:
        raise RuntimeError(f"on-path private {pg.status}")
    print("INDEPENDENT_PG",pg)

    # Independent global public BR attack. The coarse global sweep uses three
    # participation starts at every evaluated history; the selected maxima and
    # targeted dangerous histories are then rechecked with seven starts.
    g=public_br("G",XG,XG,par=P,x_grid_n=81,price_grid=51,full=False)
    b=public_br("B3",XB,XB,pbar=pg.p,par=P,x_grid_n=101,price_grid=51,full=False)
    print("INDEPENDENT_GLOBAL_G",g)
    print("INDEPENDENT_GLOBAL_B3",b)

    if g.gain > 2e-5 or abs(g.x-XG)>3e-3:
        raise RuntimeError("PROFITABLE FINITE DEVIATION IN G")
    if b.gain > 2e-5 or abs(b.x-XB)>3e-3:
        raise RuntimeError("PROFITABLE FINITE DEVIATION IN B3")

    # Recheck detected best replies with fuller continuation multistart and dense
    # private-price search.
    pr_alt=private_br(g.x,XG,P,grid_n=121,full=True)
    if pr_alt.status!=SOLVED_EQUILIBRIUM:
        raise RuntimeError(f"G best-reply private {pr_alt.status}")
    wg_c=welfare(1,XG,XG,pg.p,P,full=True)
    wg_a=welfare(1,g.x,XG,pr_alt.p,P,full=True)
    wb_c=welfare(1,XB,XB,pg.p,P,full=True)
    wb_a=welfare(1,b.x,XB,pg.p,P,full=True)
    print("FULL_CONFIRM_G",wg_c,wg_a,wg_a-wg_c,pr_alt.p)
    print("FULL_CONFIRM_B3",wb_c,wb_a,wb_a-wb_c)
    if wg_a-wg_c>2e-6 or wb_a-wb_c>2e-6:
        raise RuntimeError("FULL MULTISTART PROFITABLE DEVIATION")

    targeted_full_multistart_checks(pg.p)

    # Local strategic signs, independently finite-differenced from primitive payoff.
    for h in (6e-4,4e-4,3e-4):
        dg=local_derivatives("G",XG,h=h)
        db=local_derivatives("B3",XB,pbar=pg.p,h=h)
        print("LOCAL_SIGNS",h,"G",dg,"B3",db)
        if not (dg[0]<0 and dg[1]<0 and dg[2]<0):
            raise RuntimeError("G local substitute sign not robust")
        if not (db[0]<0 and db[1]>0 and db[2]>0):
            raise RuntimeError("B3 local complement sign not robust")

    print("STAGE4A_REPAIRED_GLOBAL_CERTIFICATION: PASS")


if __name__=="__main__":
    run()
