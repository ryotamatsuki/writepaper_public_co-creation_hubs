"""Numerical falsification for Stage 3R-4R Gateway Reboot."""
import numpy as np

SEED = 420260902
N_DRAWS = 300_000
rng = np.random.default_rng(SEED)

beta = rng.uniform(0.001,0.48,N_DRAWS)
q = 1-2*beta
rho = q*rng.uniform(0.005,0.95,N_DRAWS)
tau = rho*rng.uniform(0.001,0.999,N_DRAWS)
B = (q-rho)*rng.uniform(0.001,0.999,N_DRAWS)

r = B/q
h = (B+rho-beta*tau)/q
x = h-tau
s = (B+rho)/q

assert np.all((r>0)&(r<x)&(x<h)&(h<s)&(s<1))

D0 = (h*h-r*r)/2
D1 = (s*s-x*x)/2
Gamma = D1-D0
Gamma_closed = (
    2*beta*(1-beta)*tau*tau
    -(rho-tau)*(2*B+rho-tau)
)/(2*q*q)
assert np.max(np.abs(Gamma-Gamma_closed)) < 1e-12

X0 = (x*x-r*r)/2
X1 = (s*s-h*h)/2
G0 = D0+X0
G1 = D1+X1
assert np.all(X0>0)
assert np.all(X1>0)
assert np.max(np.abs((G1-G0)-2*Gamma)) < 1e-12

Bstar = beta*(1-beta)*tau*tau/(rho-tau)-(rho-tau)/2
assert np.sum(np.sign(Bstar-B) != np.sign(Gamma)) == 0

avgG=(G0+G1)/2
assert np.all(avgG>D0)
assert np.all(avgG>D1)

print("draws =", N_DRAWS)
print("Gamma>0 =", int(np.sum(Gamma>0)))
print("Gamma<0 =", int(np.sum(Gamma<0)))
print("B* inside canonical B-domain =", int(np.sum((Bstar>0)&(Bstar<q-rho))))
print("minimum X0 =", float(np.min(X0)))
print("minimum X1 =", float(np.min(X1)))
print("max |G1-G0-2Gamma| =", float(np.max(np.abs((G1-G0)-2*Gamma))))

# Independent fixed-point verification on a random subset.
def fixed_point(Bv,bv,rhov,tauv,state,tol=1e-13,maxit=10000):
    e1,e2=state
    bonuses=[]
    for eown,eother in ((e1,e2),(e2,e1)):
        opts=[0.0]
        if eown:
            opts.append(rhov)
        if eother:
            opts.append(rhov-tauv)
        bonuses.append(max(opts))
    Nv=0.0
    for _ in range(maxit):
        z=np.array([min(1.0,max(0.0,Bv+bv*Nv+a)) for a in bonuses])
        Nnew=float(z.sum())
        if abs(Nnew-Nv)<tol:
            return z
        Nv=Nnew
    raise RuntimeError("fixed point did not converge")

idx=rng.choice(N_DRAWS,size=10_000,replace=False)
maxerr=0.0
for k in idx:
    targets={
        (0,0):np.array([r[k],r[k]]),
        (1,0):np.array([h[k],x[k]]),
        (1,1):np.array([s[k],s[k]]),
    }
    for st,target in targets.items():
        z=fixed_point(B[k],beta[k],rho[k],tau[k],st)
        maxerr=max(maxerr,float(np.max(np.abs(z-target))))
print("max fixed-point cutoff error =", maxerr)
