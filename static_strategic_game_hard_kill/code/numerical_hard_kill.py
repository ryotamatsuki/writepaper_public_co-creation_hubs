"""Stage 4R-MC-G: numerical hard-kill diagnostics.

Provenance:
- Python 3.13.5
- random seed: 20260902
- target draws per experiment: 100_000
- sign tolerance: 1e-10

These simulations are counterexample/sanity searches, not proofs.
"""
import random

SEED = 20260902
N = 100_000
TOL = 1e-10


def experiment_a():
    """Broad admissible central-regime state draws."""
    rng = random.Random(SEED)
    neg_dp = neg_g = pos_g = neg_b2 = pos_b2 = 0
    b2_pos_g_neg = b2_neg_g_pos = 0
    accepted = attempts = 0

    while accepted < N:
        attempts += 1
        c = rng.uniform(0.02, 0.45)
        u = rng.uniform(0.02, 0.95-c)
        a = 2*c+u
        d = c+u
        x1 = rng.uniform(0.05, 1.0)
        x2 = rng.uniform(0.05, 1.0)

        if not (x1*d > x2*c and x2*d > x1*c):
            continue
        K = a*x1*x2/(x1+x2)
        upper = min(x1*d, x2*d)
        if upper <= K + 1e-5:
            continue
        m_lo = K + 1e-5
        m_hi = min(1.0, 2*upper-K-1e-5)
        if m_hi <= m_lo:
            continue
        m = rng.uniform(m_lo, m_hi)
        qstar = (m+K)/2
        kappa = rng.uniform(0.0, max(1e-6, qstar-1e-5))
        if not (qstar > kappa and qstar < upper):
            continue

        accepted += 1
        dp = -a*x2*x2/(2*(x1+x2)**2)
        cg = -a*(a*x1*x2*(x1+7*x2)+2*m*(x2*x2-x1*x1))/(4*u*(x1+x2)**4)
        cb2 = x1*x2*a*a*(x1-5*x2)/(u*(x1+x2)**4)

        neg_dp += dp < -TOL
        neg_g += cg < -TOL
        pos_g += cg > TOL
        neg_b2 += cb2 < -TOL
        pos_b2 += cb2 > TOL
        b2_pos_g_neg += cb2 > TOL and cg < -TOL
        b2_neg_g_pos += cb2 < -TOL and cg > TOL

    return {
        "attempts": attempts, "accepted": accepted, "dp_negative": neg_dp,
        "g_negative": neg_g, "g_positive": pos_g,
        "b2_negative": neg_b2, "b2_positive": pos_b2,
        "b2_complement_g_substitute": b2_pos_g_neg,
        "b2_substitute_g_complement": b2_neg_g_pos,
    }


def experiment_b():
    """Construct symmetric interior G stationary points and test planner/B3 diagnostics."""
    rng = random.Random(SEED)
    accepted = attempts = 0
    planner_pos = planner_neg = b3_positive = 0
    b2_subset = b2_over = b2_under = 0

    while accepted < N:
        attempts += 1
        c = rng.uniform(0.02, 0.45)
        u = rng.uniform(0.02, 0.95-c)
        a = 2*c+u
        x = rng.uniform(0.05, 1.0)
        m_lo = a*x/2 + 1e-5
        m_hi = min(1.0, x*(a/2+u)-1e-5)
        if m_hi <= m_lo:
            continue
        m = rng.uniform(m_lo, m_hi)
        qstar = (m+a*x/2)/2
        if qstar <= 1e-5:
            continue
        kappa = rng.uniform(0.0, qstar-1e-5)

        B = 12*c*c+28*c*u+15*u*u
        gamma = (3*B*x*x+4*a*m*x+4*m*m)/(32*u*x**3)
        soc = ((a*x-2*m)*(a*x+m)-8*gamma*u*x**3)/(8*u*x**3)
        if gamma <= 0 or soc >= -TOL:
            continue

        accepted += 1
        b3_marginal = a*a/(8*u)
        b3_positive += b3_marginal > TOL

        planner_marginal = (B*x*x-4*a*m*x-20*m*m)/(32*u*x*x)
        planner_pos += planner_marginal > TOL
        planner_neg += planner_marginal < -TOL

        xN_b2 = (8*c*c+20*c*u+11*u*u)/(8*gamma*u)
        xSP_b2 = (4*c+3*u)/(2*gamma)
        if (0 < xN_b2 < 1 and 0 < xSP_b2 < 1
                and kappa <= a*min(xN_b2, xSP_b2)/2):
            b2_subset += 1
            b2_over += xN_b2 > xSP_b2 + TOL
            b2_under += xN_b2 < xSP_b2 - TOL

    return {
        "attempts": attempts, "accepted": accepted,
        "b3_marginal_positive": b3_positive,
        "planner_positive": planner_pos, "planner_negative": planner_neg,
        "b2_subset": b2_subset, "b2_over": b2_over, "b2_under": b2_under,
    }


if __name__ == "__main__":
    print("seed", SEED, "N", N, "tol", TOL)
    print("Experiment A", experiment_a())
    print("Experiment B", experiment_b())
