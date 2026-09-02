import numpy as np

def dist_single(p,h,m,l):
    q=1-p
    return {0.0:q*q, l:p*q, m:p/2, h:p/2}

def dist_both(p,h,m,l,referral=True):
    q=1-p
    if referral:
        return {0.0:q**3, l:p*q, m:p*q*q, h:p}
    return {0.0:q*q, l:p*q, h:p}

def moments(d):
    mu=sum(x*pr for x,pr in d.items())
    nu=sum(x*x*pr for x,pr in d.items())
    return mu,nu

def market_state(a,xi,xj):
    qi=(a+2*xi-xj)/3
    qj=(a+2*xj-xi)/3
    pi=qi*qi
    Q=qi+qj
    cs=Q*Q/2
    return pi+cs/2

def enumerate_K(a,d):
    return sum(pi*pj*market_state(a,xi,xj)
               for xi,pi in d.items() for xj,pj in d.items())

def moment_K(a,d):
    mu,nu=moments(d)
    return 2*a*a/9+4*a*mu/9+11*nu/18-7*mu*mu/18

rng=np.random.default_rng(20260902)
maxerr=0.0
for _ in range(10_000):
    p=rng.uniform(.01,.99)
    l=rng.uniform(.02,1)
    m=l+rng.uniform(.02,1.5)
    h=m+rng.uniform(.02,1.5)
    a=h+rng.uniform(.02,5)
    for d in (dist_single(p,h,m,l),dist_both(p,h,m,l,True),dist_both(p,h,m,l,False)):
        maxerr=max(maxerr,abs(enumerate_K(a,d)-moment_K(a,d)))
print("max absolute state-enumeration error:",maxerr)
assert maxerr < 1e-11
