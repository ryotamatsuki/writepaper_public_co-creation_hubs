"""Stage 13: endogenous-location interior candidate and two globally
valid but payoff-distinct off-path Nash fee equilibria. Standard library only.
This does NOT certify a global location SPNE.
"""
from math import isclose

r=1/45  # mu=0.1, beta=2/3, A=c=t=M=R=1; u=tau=2/9
poly=lambda d:-2*d**3+(1+11*r)*d*d-(6*r+8*r*r)*d+6*r*r

lo,hi=.45,.55
assert poly(lo)*poly(hi)<0
for _ in range(100):
    mid=(lo+hi)/2
    if poly(lo)*poly(mid)<=0: hi=mid
    else:lo=mid
d=(lo+hi)/2
x=(d*(1-d)-2*r)/(d-4*r)
a=(1-x)/2
b=a+d
pP=2*(d-r)*x
pG=2*d*x-r
profit=2*x*pP

def intsq(loc,l,h):return ((h-loc)**3-(l-loc)**3)/3

W2=8-2*(intsq(a,0,1-x)+intsq(b,1-x,1))+2*r*(x*x+(1-x)**2)
W1=8-2/12+2*r
C=W2-W1
H=8-2/4+2*r

# A far-from-interior PRIVATE LOCATION deviation.
dd=.01
bnew=a+dd
m=(a+bnew)/2
Delta=r-2*dd*m
assert dd<r and Delta>0 and abs(bnew-.5)<abs(a-.5)

def revenue_upper(X,Y,q):
    return ((q+2*dd*(1-m)-r)*(X+Y)
            -2*dd*(X*X+Y*Y)+4*r*X*Y)
for q in (0,1):
    # Proof, not the grid: set s=(X+Y)/2,h=(X-Y)/2.
    # Q=2(q+2dd(1-m)-r)s+4(r-dd)s²-4(r+dd)h².
    # Since dd<r, this is convex in s and maximized at s=0 or 1.
    # P-all yields 2(q+Delta)>0; hence global private fee BR.
    # Public globally maximizes W by serving everyone at P:
    # socially efficient W is convex in s when dd<2r.
    assert 4-(q+Delta)-max(bnew*bnew,(1-bnew)**2)>0
    for i in range(101):
        for j in range(101):
            assert revenue_upper(i/100,j/100,q)<=2*(q+Delta)+1e-10
    print("off-path fee Nash public q:",q,
          "private both-side p:",q+Delta,"private gross:",2*(q+Delta))

assert 2*Delta<profit<2*(1+Delta)
assert d>2*r and 0<a<b<1 and 0<x<1 and pG>0 and pP>0
print("stationary interior: D,x,a,b:",d,x,a,b)
print("public/private joining fees:",pG,pP)
print("private gross P:",profit)
print("gross W mixed/single, marginal C:",W2,W1,C)
print("P-only monopoly gross H:",H)
print("CONDITIONAL stage-0: FP<=P -> build iff FG<=C;")
print("P<FP<=H -> build iff FG<=FP; FP>H -> build iff FG<=W1.")
print("PASS: off-path high-price best-response lemma; full SPNE OPEN")
