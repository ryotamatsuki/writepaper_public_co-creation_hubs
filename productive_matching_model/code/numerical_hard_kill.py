"""Numerical falsification for Stage 3M-4M productive-matching model."""

import numpy as np

SEED = 20260902
N = 500_000
TOL = 1e-12
rng = np.random.default_rng(SEED)

# Admissible draws: A>Delta>0, p in (0,1), r in [0,1].
A = rng.uniform(0.2, 4.0, N)
Delta = A * rng.uniform(0.02, 0.95, N)
p = rng.uniform(0.005, 0.95, N)
r = rng.uniform(0.0, 1.0, N)

h = 2*p - p*p
k = (1-p)**2 * r

Gamma0 = -Delta*p*(8*A + Delta*(11 + 7*p**3 - 14*p**2))/18
Gamma1 = Gamma0 + Delta*k*(8*A + Delta*(11 - 14*h - 7*k))/18
Gamma_ns = k-p  # positive scaling V omitted

print("M0 positive / negative / near-zero:",
      np.sum(Gamma0>TOL), np.sum(Gamma0<-TOL), np.sum(np.abs(Gamma0)<=TOL))
print("M1 positive / negative / near-zero:",
      np.sum(Gamma1>TOL), np.sum(Gamma1<-TOL), np.sum(np.abs(Gamma1)<=TOL))
print("Nonstrategic positive / negative:", np.sum(Gamma_ns>0), np.sum(Gamma_ns<0))
print("NS positive but Cournot negative:", np.sum((Gamma_ns>0) & (Gamma1<0)))
print("NS negative but Cournot positive:", np.sum((Gamma_ns<0) & (Gamma1>0)))

# Closed-form expected regional welfare.
def R_closed(A, D, u, v):
    return (
        8*A*A + 20*A*D*u - 4*A*D*v - 14*D*D*u*v + 17*D*D*u + 5*D*D*v
    ) / 36

# Independent enumeration of four Bernoulli matching states.
def R_enum(A, D, u, v):
    total = 0.0
    for mi in (0,1):
        for mj in (0,1):
            prob = (u if mi else 1-u) * (v if mj else 1-v)
            qi = (A + 2*D*mi - D*mj)/3
            qj = (A + 2*D*mj - D*mi)/3
            Q = qi+qj
            total += prob * (qi*qi + Q*Q/4)
    return total

max_err = 0.0
for _ in range(10_000):
    aa = rng.uniform(0.2,3.0)
    dd = rng.uniform(0.001,0.99*aa)
    uu = rng.uniform(0,1)
    vv = rng.uniform(0,1)
    max_err = max(max_err, abs(R_closed(aa,dd,uu,vv)-R_enum(aa,dd,uu,vv)))
print("max expected-welfare enumeration error:", max_err)

# Local open-set perturbation checks.
def gamma(A,D,p,r):
    h = 2*p-p*p
    k = (1-p)**2*r
    g0 = -D*p*(8*A+D*(11+7*p**3-14*p**2))/18
    return g0 + D*k*(8*A+D*(11-14*h-7*k))/18

def perturb(base, eps=0.05, n=100_000):
    factors = rng.uniform(1-eps,1+eps,(n,4))
    vals = np.asarray(base)*factors
    aa,dd,pp,rr = vals.T
    valid = (aa>dd) & (dd>0) & (pp>0) & (pp<1) & (rr>=0) & (rr<=1)
    gg = gamma(aa[valid],dd[valid],pp[valid],rr[valid])
    return valid.sum(), np.sum(gg>0), np.sum(gg<0), gg.min(), gg.max()

print("negative-witness perturbation:", perturb((1.0,0.5,0.1,0.1)))
print("positive-witness perturbation:", perturb((1.0,0.5,0.1,0.2)))

# Expected reference output with this seed and sampling design:
# M0: 0 / 500000 / 0
# M1: 106678 / 393322 / 0
# nonstrategic: 127377 / 372623
# NS positive but Cournot negative: 20699
# NS negative but Cournot positive: 0
# enumeration max error should be approximately machine precision (~1e-15).
