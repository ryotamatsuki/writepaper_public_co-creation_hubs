"""Exact-fraction sufficient-condition verifier for Stage 15.

The symbolic proof is in STAGE_15_PURE_SPNE_NONEXISTENCE_OPEN_SET.md.
These tests check strict inequalities, NOT fee strategy space exhaustively.
IMPORTANT: theorem assumes mandatory full subscription, not finite-v opt-out.
"""
from fractions import Fraction as Q


def certificate(r, d, m):
    assert 0 < r < d < 2*r and Q(1,2) < m < 1
    a=m-d/2; b=m+d/2
    assert 0<a<b<1
    z=2*m-1
    bound=(2*r-d)*(d-r)/(d*(3*d+2*r))
    assert 0 < z < bound
    A0=2*d*(1-m)-r
    sbar=d*z/(2*r-d)
    szero=A0/(4*d)
    sint=A0/(4*(d-r))
    pmin=2*(d-r)
    welfare_G_minus_P=2*d*z
    assert A0>0 and 0<sbar < szero < sint
    assert sbar<Q(1,2)
    assert pmin>A0
    assert welfare_G_minus_P>0
    return dict(a=a,b=b,r=r,d=d,m=m,A0=A0,sbar=sbar,
                szero=szero,sint=sint,pmin=pmin,
                welfare_G_minus_P=welfare_G_minus_P,
                open_parameter_bound=bound)


if __name__=="__main__":
    example=certificate(Q(1,45),Q(3,100),Q(51,100))
    assert example["A0"]==Q(323,45000)
    assert example["sbar"]==Q(27,650)
    assert example["szero"]==Q(323,5400)
    assert example["sint"]==Q(323,1400)
    assert example["pmin"]==Q(7,450)
    assert example["welfare_G_minus_P"]==Q(3,2500)
    print("Exact rational witness:",example)
    for r in (Q(1,100000),Q(1,1000),Q(1,45),Q(1,10),Q(2,9)):
        certificate(r,Q(3,2)*r,Q(101,200))
    print("PASS: strict OPEN sufficient conditions and arbitrarily weak network strengths")
