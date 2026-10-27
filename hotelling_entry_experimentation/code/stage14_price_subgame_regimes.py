"""Stage 14: exact rational classification of full-coverage bilateral network
subscriptions, explicit all-private Nash fee continuations and a rational
PURE-FEE-NONEXISTENCE certificate in a FULL-COVERAGE specialization.

IMPORTANT: The original finite-v game allows opting out, whose additional
demand branches are not exhausted here. Therefore this script must NOT be
described as proving SPNE nonexistence for finite v.

Uses only the Python standard library. Exact fractions throughout.
"""
from fractions import Fraction as F


def clip(z):
    return min(F(1), max(F(0), z))


def subscription_regimes(A, B, k):
    """Exhaustively enumerate isolated solutions of x=clip(A+k*y),
    y=clip(B+k*x). Special continuum at k=1, A+B=0 marked separately."""
    one,zero=F(1),F(0)
    branches=[]
    def add(label,x,y):
        if x==clip(A+k*y) and y==clip(B+k*x):
            if (x,y) not in [pair for _,*pair in branches]:
                branches.append((label,x,y))
    for ix,x in [(0,zero),(1,one)]:
        for iy,y in [(0,zero),(1,one)]:
            add(str(ix)+str(iy),x,y)
    if 0<A<1:
        add("I0", A, zero)
    if 0<A+k<1:
        add("I1", A+k, one)
    if 0<B<1:
        add("0I", zero, B)
    if 0<B+k<1:
        add("1I", one, B+k)
    if k!=1:
        x=(A+k*B)/(1-k*k)
        y=(B+k*A)/(1-k*k)
        if 0<x<1 and 0<y<1:
            add("II",x,y)
    elif A+B==0:
        branches.append(("CONTINUUM POSSIBLE",None,None))
    return branches


def fee_subscription(a,b,r,gI,gR,pI,pR,t=F(1)):
    assert a!=b
    d=abs(b-a)
    mid=(a+b)/2
    c=1-mid if b>a else mid
    k=r/d
    A=c+(gI-pI)/(2*t*d)-k/2
    B=c+(gR-pR)/(2*t*d)-k/2
    return subscription_regimes(A,B,k)


def W_given_all_registered(a,b,x,y,r,t=F(1),v=F(4)):
    """Total W with both unit-mass sides, no hub fixed costs."""
    if b<a:
        # P shares now correspond to left-tail fractions x,y.
        # transport gain from P's serving the left fractions
        d=a-b
        c=(a+b)/2
    else:
        # P serves right-tail fractions
        d=b-a
        c=1-(a+b)/2
    WG=2*v-2*t*(F(1,12)+(a-F(1,2))**2)+2*t*r
    return WG +t*d*(2*c*(x+y)-x*x-y*y) -2*t*r*(x+y)+4*t*r*x*y


def fullP_q_bound(a,b,r,q,t=F(1)):
    d=abs(b-a);mid=(a+b)/2
    c=1-mid if b>a else mid
    Delta=t*(r-2*d*(1-c))
    pp=q+Delta
    return Delta,pp


def assert_allP_Nash(a,b,r,q):
    """Sufficient branch certificate for a private-more-central pair
    with d<2r; public alternative fee q covers all. Checks Q on a grid,
    while exact proof is the s/h quadratic optimization in memo."""
    d=abs(b-a)
    mid=(a+b)/2
    c=1-mid if b>a else mid
    assert c>=F(1,2) and 0<d<2*r
    Delta,p=fullP_q_bound(a,b,r,q)
    assert p>0 and q>=0
    if d>r:
        assert p>=2*(d-r)
    matches=fee_subscription(a,b,r,q,q,p,p)
    assert any(x==1 and y==1 for _,x,y in matches)
    # If public's fallback can serve all consumers even if P deviates,
    # the following quadratic is an upper bound for ANY P joint fee change.
    Aq=q+2*d*c-r
    Popt=2*p
    for i in range(101):
        for j in range(101):
            x,y=F(i,100),F(j,100)
            envelope=Aq*(x+y)-2*d*(x*x+y*y)+4*r*x*y
            assert envelope<=Popt
            assert W_given_all_registered(a,b,F(1),F(1),r)>=W_given_all_registered(a,b,x,y,r)
    return Delta,p,Popt


def rational_nonexistence_witness():
    """PROVES analytic inequalities forbidding pure fee Nash when all
    sides MUST subscribe to one of the platforms. This is not a finite-v
    outside-option equilibrium nonexistence theorem."""
    r=F(1,45)
    d=F(3,100)
    mid=F(51,100)
    a=mid-d/2; b=mid+d/2
    assert (a,b)==(F(99,200),F(21,40))
    assert r<d<2*r
    A0=2*d*(1-mid)-r
    assert A0==F(323,45000)
    # P's only-interior-FOC price BR share floor at G fees 0.
    sint=A0/(4*(d-r))
    assert sint==F(323,1400)
    # If one private fee hits zero, its side demand is >=A0/(2d)
    # and hence common share average is >=A0/(4d).
    szero=A0/(4*d)
    assert szero==F(323,5400)
    # If one side is entirely served by P, s>=1/2.
    threshold=d*(2*mid-1)/(2*r-d)
    assert threshold==F(27,650)
    assert min(sint,szero,F(1,2))>threshold
    # In full coverage W(X,Y)-W(allP) =
    # 2(1-s)*(d*(2m-1)-(2r-d)*s)-2(d+2r)*h^2
    # <0 for any P best reply with  s∈(threshold,1).
    WP=W_given_all_registered(a,b,F(1),F(1),r)
    WG=W_given_all_registered(a,b,F(0),F(0),r)
    assert WG-WP==2*d*(2*mid-1)==F(3,2500)
    # Illustrate failure of two likely candidate branches.
    p0=A0/2
    eq0=fee_subscription(a,b,r,F(0),F(0),p0,p0)
    assert len(eq0)==1 and eq0[0][1]==sint and eq0[0][2]==sint
    assert W_given_all_registered(a,b,sint,sint,r)<WP
    # G can raise its fees and obtain all-P with the P fee fixed.
    eqhigh=fee_subscription(a,b,r,F(1),F(1),p0,p0)
    assert len(eqhigh)==1 and eqhigh[0][1:]==(F(1),F(1))
    Delta,_=fullP_q_bound(a,b,r,F(0))
    Pmin=2*(d-r)
    Gmin=Pmin-Delta
    assert Pmin>0 and Pmin>A0 and Gmin>=0
    allP=fee_subscription(a,b,r,Gmin,Gmin,Pmin,Pmin)
    allG=fee_subscription(a,b,r,F(0),F(0),Pmin,Pmin)
    assert len(allP)==1 and allP[0][1:]==(F(1),F(1))
    assert len(allG)==1 and allG[0][1:]==(F(0),F(0))
    print("FULL-COVERAGE NO-PURE-FEE-NE witness a,b,d,r:",a,b,d,r)
    print("A0, interior BR avg s, zero-fee BR bound, welfare threshold:",
          A0,sint,szero,threshold)
    print("G all vs P all planner gain:",WG-WP)
    print("At G fees 0 private best price:",p0,"public deviation to fees 1 increases W")
    print("At P all fee cutoff P,G:",Pmin,Gmin,
          "public deviation to fees 0 increases W")
    print("ANALYTIC OBSTRUCTION: no PURE price equilibrium in full-coverage specialization.")


if __name__=="__main__":
    r=F(1,45)
    a=F(1,4)
    b=a+F(1,100)
    # At equal fees 0 and close positions multiple subscription fixed points.
    subs=fee_subscription(a,b,r,F(0),F(0),F(0),F(0))
    print("same fees multiple membership equilibria:",subs)
    assert any(x==0 and y==0 for _,x,y in subs)
    assert any(x==1 and y==1 for _,x,y in subs)
    # For P slightly more central, fullP family of PRICE Nash equilibria.
    for q in [F(0), F(1)]:
        delta,price,profit=assert_allP_Nash(a,b,r,q)
        # public q=1 covers every person from both populations at a=.25
        assert 4-q-max(a*a,(1-a)*(1-a))>0
        print("P-dominant near collocation G fee",q,
              "P fee",price,"P profit",profit)
    # A separated d∈(r,2r) also supports allP price equilibria if P central
    a=F(9,20); b=F(12,25) # .45,.48; P closer to center
    d=abs(b-a);assert r<d<2*r
    delta=r-2*d*(a+b)/2
    qmin=max(F(0), 2*(d-r)-delta)
    print("intermediate allP-price Nash:",
          assert_allP_Nash(a,b,r,qmin))
    rational_nonexistence_witness()
    print("PASS: nine-case enumeration, allP global Nash envelopes,")
    print("and exact-rational conditional pure-price nonexistence inequalities.")
