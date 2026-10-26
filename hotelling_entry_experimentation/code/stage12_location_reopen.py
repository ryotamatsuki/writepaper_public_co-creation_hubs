"""Stage 12: endogenous-location preliminary analytics and Stage 11 counterexample.

Exact rational mixed-duopoly fee/location response for a fixed public a,
and a candidate symmetric-side cubic. No global location SPNE is claimed.
Standard-library-only.
"""
from fractions import Fraction as F
from math import sqrt

M=R=tI=tR=F(1)
mu=F(1,2)
u,tau=F(1,8),F(1,4)
S=u+tau


def interior_price_equilibrium(a,b):
    assert 0 <= a < b <= 1
    d=b-a
    mid=(a+b)/2
    hi=tI*d
    hr=tR*d
    gI=mu*R*S
    gR=mu*M*S
    det=hi*hr-gI*gR
    assert det>0
    shift=F(1,2)-mid
    x=F(1,2)+shift*hr*(hi+gI)/det
    y=F(1,2)+shift*hi*(hr+gR)/det
    assert 0<x<1 and 0<y<1
    pPI=2*(hi*x-mu*R*tau*y)
    pPR=2*(hr*y-mu*M*u*x)
    pGI=pPI+2*mu*R*tau*(y-F(1,2))
    pGR=pPR+2*mu*M*u*(x-F(1,2))
    assert min(pPI,pPR,pGI,pGR)>0
    profit=M*x*pPI+R*y*pPR
    return dict(x=x,y=y,pPI=pPI,pPR=pPR,pGI=pGI,pGR=pGR,
                profit=profit)


def interior_location_cubic(d,r):
    return -2*d**3+(1+11*r)*d**2-(6*r+8*r*r)*d+6*r*r


def bisect(func,lo,hi,n=100):
    flo=func(lo)
    assert flo*func(hi)<0
    for _ in range(n):
        mid=(lo+hi)/2
        if flo*func(mid)<=0:
            hi=mid
        else:
            lo=mid
            flo=func(mid)
    return (lo+hi)/2


if __name__=="__main__":
    original=interior_price_equilibrium(F(1,4),F(3,4))
    deviation=interior_price_equilibrium(F(1,4),F(37,50))
    assert original["profit"]==F(13,32)
    assert deviation["profit"]>original["profit"]
    print("Stage11 a=.25,b=.75:",{k:float(v) for k,v in original.items()})
    print("Deviation a=.25,b=.74:",{k:float(v) for k,v in deviation.items()})
    print("Profitable deviation:",float(deviation["profit"]-original["profit"]))
    # For symmetric side rewards (beta=2/3) and mu=1/10:
    # u=tau=2/9, r=mu*u/t=1/45.
    r=1/45
    d=bisect(lambda D:interior_location_cubic(D,r),.45,.50)
    x=(d*(1-d)-2*r)/(d-4*r)
    a=(1-x)/2
    b=a+d
    assert 0<a<1 and 0<b<1 and 0<x<1
    assert abs(interior_location_cubic(d,r))<1e-12
    assert abs((d-r)*(1-2*d)-r*x)<1e-12
    print("Interior cubic local candidate: d,x,a,b",d,x,a,b)
    print("PASS: Stage11 locations are NOT location-equilibrium positions")
