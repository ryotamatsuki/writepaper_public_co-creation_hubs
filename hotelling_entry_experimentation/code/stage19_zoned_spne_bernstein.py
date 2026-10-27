"""Stage 19: exact rational Bernstein certificates and solved equilibrium
for a new ZONED game: G a∈[0,2/5], P b∈[3/5,1], even in monopoly.

Requires sympy and scipy. Full model/proof: STAGE_19_*_SPNE.md.
IMPORTANT: not an unrestricted-location game.
"""
import sympy as s
from scipy.optimize import root
from fractions import Fraction

a,b=s.symbols("a b", real=True)
k=s.Rational(1,450)
h=b-a
D=s.expand(h*(h-2*k)-4*k*k)
C=(1-a-b)/2
nX=s.expand(D/2+h*h*C)
nY=s.expand(D/2+h*(h+2*k)*C)
def ex(z): return s.expand(z)
pI=ex(2*h*nX-k*nY)
pR=ex(-2*k*nX+(2*h-k)*nY)
gI=ex(pI+k*(2*nY-D))
gR=ex(pR+3*k*(nX+nY-D))
piNum=ex(2*h*(nX*nX+nY*nY)-3*k*nX*nY-k*nY*nY)
# W_a = envelope derivative; pi_b = full derivative of equilibrium profit
waNum=ex((D-nX)**2+(D-nY)**2-2*a*D*(2*D-nX-nY))
pbNum=ex(s.diff(piNum,b)*D-2*piNum*s.diff(D,b))
waaNum=ex(s.diff(waNum,a)*D-2*waNum*s.diff(D,a))
wbaNum=ex(s.diff(waNum,b)*D-2*waNum*s.diff(D,b))
pbbNum=ex(s.diff(pbNum,b)*D-3*pbNum*s.diff(D,b))
pbaNum=ex(s.diff(pbNum,a)*D-3*pbNum*s.diff(D,a))

def bern_bounds(poly, a0, a1, b0, b1):
    """Exact min/max of Bernstein CONTROL COEFFICIENTS on rectangle.
    If all coefficients same strict sign -> polynomial has strict sign."""
    u,v=s.symbols("u v")
    poly=s.Poly(ex(poly.subs({a:a0+(a1-a0)*u,b:b0+(b1-b0)*v}, simultaneous=True)),u,v)
    n,m=poly.degree_list()
    coeff=dict(poly.terms())
    bounds=[]
    for i in range(n+1):
        for j in range(m+1):
            z=sum(coeff.get((l,q),0)*s.binomial(i,l)/s.binomial(n,l)*s.binomial(j,q)/s.binomial(m,q) for l in range(i+1) for q in range(j+1))
            bounds.append(z)
    return min(bounds), max(bounds)

def pos(poly, rect, name):
    mn,mx=bern_bounds(poly,*rect)
    assert mn>0, (name,"not Bernstein-positive",mn,mx)
    return mn,mx

def neg(poly, rect, name):
    mn,mx=bern_bounds(poly,*rect)
    assert mx<0, (name,"not Bernstein-negative",mn,mx)
    return mn,mx

Z=(s.Rational(0),s.Rational(2,5),s.Rational(3,5),s.Rational(1))
for name,pol in [
    ("D",D),("X*D",nX),("(1-X)*D",D-nX),
    ("Y*D",nY),("(1-Y)*D",D-nY),
    ("private I fee*D",pI),("private R fee*D",pR),
    ("public I fee*D",gI),("public R fee*D",gR),
    ("upper GI",s.Rational(5,2)*D-gI),("upper GR",s.Rational(5,2)*D-gR),
]:
    pos(pol,Z,name)
neg(waaNum,Z,"Waa<0")
neg(pbbNum,Z,"Pbb<0")

Q=lambda x:s.Rational(x,100)
pos(waNum,(Q(0),Q(0),Q(60),Q(100)),"G derivative at a=0")
neg(waNum,(Q(40),Q(40),Q(60),Q(100)),"G derivative at a=.4")
pos(pbNum,(Q(0),Q(40),Q(60),Q(60)),"P derivative at b=.6")
neg(pbNum,(Q(0),Q(40),Q(100),Q(100)),"P derivative at b=1")

# First best-reply image.
pos(waNum,(Q(19),Q(19),Q(60),Q(100)),"G derivative at .19")
neg(waNum,(Q(35),Q(35),Q(60),Q(100)),"G derivative at .35")
pos(pbNum,(Q(0),Q(40),Q(65),Q(65)),"P derivative at .65")
neg(pbNum,(Q(0),Q(40),Q(85),Q(85)),"P derivative at .85")
# Second image, given first image of opponent.
pos(waNum,(Q(21),Q(21),Q(65),Q(85)),"G derivative at .21")
neg(waNum,(Q(29),Q(29),Q(65),Q(85)),"G derivative at .29")
pos(pbNum,(Q(19),Q(35),Q(72),Q(72)),"P derivative at .72")
neg(pbNum,(Q(19),Q(35),Q(80),Q(80)),"P derivative at .80")

Z2=(Q(21),Q(29),Q(72),Q(80))
pos(wbaNum,Z2,"W_ab>0")
pos(pbaNum,Z2,"P_ba>0")
# D^3 (or D^4) denominators are common for each curvature pair.
neg(wbaNum+s.Rational(1,2)*waaNum,Z2,"W_ab<-W_aa/2")
neg(pbaNum+s.Rational(1,2)*pbbNum,Z2,"P_ba<-P_bb/2")

waf=s.lambdify((a,b),waNum/D**2,"numpy")
pbf=s.lambdify((a,b),pbNum/D**3,"numpy")
sol=root(lambda v:[waf(*v),pbf(*v)], [.25,.75],tol=1e-11)
assert sol.success or max(abs(q) for q in sol.fun)<1e-9
aa,bb=sol.x
xf=s.lambdify((a,b),nX/D,"numpy")
yf=s.lambdify((a,b),nY/D,"numpy")
fee_fns=[s.lambdify((a,b),z/D,"numpy") for z in [gI,gR,pI,pR]]
X,Y=float(xf(aa,bb)),float(yf(aa,bb))
g_I,g_R,p_I,p_R=(float(fn(aa,bb)) for fn in fee_fns)
pi=X*p_I+Y*p_R
def mismatch(lo,hi,loc):
    return ((hi-loc)**3-(lo-loc)**3)/3
trans=(mismatch(0,1-X,aa)+mismatch(1-X,1,bb)+
       mismatch(0,1-Y,aa)+mismatch(1-Y,1,bb))
grossW=8-trans+2*float(k)*(X*Y*Y+(1-X)*(1-Y)**2)
N=s.Rational(1759,225)
H=s.Rational(1639,225)
Delta=grossW-float(N)
assert 0<Delta<pi<float(H)<float(N)
print("PASS exact polynomial Bernstein global curvature, boundary and contraction certificates")
print("GLOBAL ZONED LOCATION EQUILIBRIUM: a,b",aa,bb)
print("Shares P I,R",X,Y)
print("Fees G I,R",g_I,g_R,"P I,R",p_I,p_R)
print("Private gross fee revenue",pi)
print("Gross welfare duopoly",grossW, "single",float(N))
print("Single-private gross revenue",float(H),"marginal G welfare",Delta)
print("Matched links P,G",.1*X*Y,.1*(1-X)*(1-Y))
print("Pre-match effort P,G",(2/3)*.1*Y,(2/3)*.1*(1-Y))
for fg,fp in ((.05,.05),(.1,.3),(8,7.5),(.4,1)):
    W00=0;W10=float(N)-fg;W01=float(N)-fp;W11=grossW-fg-fp
    ne=[]
    if fg>=float(N) and fp>=float(H):ne.append("none")
    if fg<=float(N) and fp>=pi:ne.append("public")
    if fg>=Delta and fp<=float(H):ne.append("private")
    if fg<=Delta and fp<=pi:ne.append("both")
    allW={"none":W00,"public":W10,"private":W01,"both":W11}
    print("costs",fg,fp,"NE",ne,"planner",max(allW,key=allW.get),
          "Welfare",allW)
