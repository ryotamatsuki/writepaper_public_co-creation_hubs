from __future__ import annotations
import itertools
from fractions import Fraction
import sympy as sp
from exact_model import build

M=build(); t,q,p=M['t'],M['q'],M['p']; H,R,A,B=M['H'],M['Rred'],M['AGnum'],M['Bnum']
assert sp.factor(H)==50*p+1000*q**2-50*q*t-125*q+1
psol=sp.solve(H,p)[0]
R2=sp.cancel(sp.primitive(sp.Poly(sp.expand(R.subs(p,psol)),t,q))[1].as_expr()/q**2)
A2=sp.cancel(sp.primitive(sp.Poly(sp.expand(A.subs(p,psol)),t,q))[1].as_expr()/q**2)
resT=sp.factor(sp.resultant(R2,A2,q)); degs=sorted((sp.degree(f,t),e) for f,e in sp.factor_list(resT)[1])
assert degs==[(1,3),(1,8),(2,3),(31,1)]
resB=sp.factor(sp.resultant(H,B,q)); bf=sp.factor_list(resB)[1]
assert sorted((sp.degree(f,t),sp.degree(f,p),e) for f,e in bf)==[(0,1,2),(14,7,1)]
print('elimination: G core degree 31; B3 matched-price core degree 14: PASS')

class I:
    def __init__(self,lo,hi=None):
        self.lo=Fraction(lo); self.hi=Fraction(lo if hi is None else hi)
        if self.lo>self.hi:self.lo,self.hi=self.hi,self.lo
    def __add__(self,o):o=o if isinstance(o,I) else I(o);return I(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,o):return self+(-(o if isinstance(o,I) else I(o)))
    def __rsub__(self,o):return I(o)-self
    def __mul__(self,o):
        o=o if isinstance(o,I) else I(o); z=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi];return I(min(z),max(z))
    __rmul__=__mul__
    def inv(self):
        assert not (self.lo<=0<=self.hi)
        return I(1/self.hi,1/self.lo) if self.lo>0 else I(1/self.lo,1/self.hi)
    def __truediv__(self,o):return self*(o if isinstance(o,I) else I(o)).inv()
    def __rtruediv__(self,o):return I(o)/self
    def __pow__(self,n):
        if n<0:return self.inv()**(-n)
        if n==0:return I(1)
        if n%2==0 and self.lo<=0<=self.hi:return I(0,max(abs(self.lo),abs(self.hi))**n)
        z=[self.lo**n,self.hi**n];return I(min(z),max(z))
    def sign(self):return 1 if self.lo>0 else (-1 if self.hi<0 else 0)
    def mid(self):return (self.lo+self.hi)/2
    def rad(self):return (self.hi-self.lo)/2
    def subset_interior(self,o):return self.lo>o.lo and self.hi<o.hi
    def fp(self):return (float(self.lo),float(self.hi))

def polyI(expr,vars,boxes):
    P=sp.Poly(expr,*vars,domain=sp.QQ); out=I(0)
    for mon,c in P.terms():
        term=I(Fraction(int(c.p),int(c.q)))
        for e,X in zip(mon,boxes):term=term*(X**e)
        out=out+term
    return out

def ev(expr,env,cache):
    if expr in cache:return cache[expr]
    if expr.is_Symbol:r=env[expr]
    elif expr.is_Integer:r=I(int(expr))
    elif expr.is_Rational:r=I(Fraction(int(expr.p),int(expr.q)))
    elif expr.is_Add:
        r=I(0)
        for a in expr.args:r=r+ev(a,env,cache)
    elif expr.is_Mul:
        r=I(1)
        for a in expr.args:r=r*ev(a,env,cache)
    elif expr.is_Pow and expr.exp.is_Integer:r=ev(expr.base,env,cache)**int(expr.exp)
    else:raise TypeError(expr)
    cache[expr]=r;return r

def det(A):
    n=len(A); out=I(0)
    for perm in itertools.permutations(range(n)):
        inv=sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n)); term=I(-1 if inv%2 else 1)
        for i,j in enumerate(perm):term=term*A[i][j]
        out=out+term
    return out

def invmat(A):
    n=len(A); D=det(A); assert D.sign()!=0; C=[[None]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            minor=[[A[r][c] for c in range(n) if c!=i] for r in range(n) if r!=j]
            C[i][j]=(((-1)**(i+j))*det(minor))/D
    return C

def mv(A,v):return [sum((A[i][j]*v[j] for j in range(len(v))),I(0)) for i in range(len(A))]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),I(0))

TB,QB=sp.symbols('tb qb'); vars=[t,q,p,TB,QB]
FB=[H,R,A,H.subs({t:TB,q:QB}),B.subs({t:TB,q:QB})]
centers=['0.5848900730456509736','0.13885157350636662587','0.02274665143268023016','0.61009086253699917757','0.14026693513637059270']
boxes=[]; den=10**10
for s in centers:
    y=Fraction(s); n=y.numerator*den//y.denominator; boxes.append(I(Fraction(n,den),Fraction(n+1,den)))
mid=[X.mid() for X in boxes]; sub={z:sp.Rational(m.numerator,m.denominator) for z,m in zip(vars,mid)}
J=sp.Matrix([[sp.diff(f,z).subs(sub) for z in vars] for f in FB]); C=J.inv(); Cf=[[Fraction(int(C[i,j].p),int(C[i,j].q)) for j in range(5)] for i in range(5)]
fm=[]
for f in FB:
    z=sp.factor(f.subs(sub)); fm.append(Fraction(int(z.p),int(z.q)))
JX=[[polyI(sp.diff(f,z),vars,boxes) for z in vars] for f in FB]
MM=[]
for i in range(5):
    row=[]
    for j in range(5):
        s=I(0)
        for k0 in range(5):s=s+Cf[i][k0]*JX[k0][j]
        row.append(I(1 if i==j else 0)-s)
    MM.append(row)
K=[]
for i in range(5):
    c=mid[i]-sum((Cf[i][k0]*fm[k0] for k0 in range(5)),Fraction(0)); Y=I(c)
    for j in range(5):Y=Y+MM[i][j]*I(-boxes[j].rad(),boxes[j].rad())
    K.append(Y)
assert all(K[i].subset_interior(boxes[i]) for i in range(5))
TG,QG,PG,TBB,QBB=K
print('coupled G/B3 exact rational Krawczyk certificate: PASS')

def reduced_hessian(E,zvars,xvars,W,env):
    cache={}; AA=[[ev(sp.diff(E[r],zvars[c]),env,cache) for c in range(len(zvars))] for r in range(len(E))]; Ainv=invmat(AA)
    zr=[]
    for xx in xvars:zr.append([-u for u in mv(Ainv,[ev(sp.diff(E[r],xx),env,cache) for r in range(len(E))])])
    Wz=[ev(sp.diff(W,z),env,cache) for z in zvars]; HH=[]
    for ii,xi in enumerate(xvars):
        row=[]; vi=zr[ii]
        for jj,xj in enumerate(xvars):
            vj=zr[jj]; cc=[]
            for r in range(len(E)):
                zc=ev(sp.diff(E[r],xi,xj),env,cache)
                for a,za in enumerate(zvars):zc=zc+ev(sp.diff(E[r],xi,za),env,cache)*vj[a]+ev(sp.diff(E[r],za,xj),env,cache)*vi[a]
                for a,za in enumerate(zvars):
                    for b,zb in enumerate(zvars):zc=zc+ev(sp.diff(E[r],za,zb),env,cache)*vi[a]*vj[b]
                cc.append(zc)
            zij=[-u for u in mv(Ainv,cc)]
            zc=ev(sp.diff(W,xi,xj),env,cache)
            for a,za in enumerate(zvars):zc=zc+ev(sp.diff(W,xi,za),env,cache)*vj[a]+ev(sp.diff(W,za,xj),env,cache)*vi[a]
            for a,za in enumerate(zvars):
                for b,zb in enumerate(zvars):zc=zc+ev(sp.diff(W,za,zb),env,cache)*vi[a]*vj[b]
            row.append(zc+dot(Wz,zij))
        HH.append(row)
    return HH,zr

v,alpha,beta,rho,rhoT,k,kT,gamma=[M[z] for z in ['v','alpha','beta','rho','rhoT','k','kT','gamma']]
t1,t2,x1,x2=M['t1'],M['t2'],M['x1'],M['x2']; G,Efull,W,pi=M['G'],M['Efull'],M['W'],M['pi']; y=M['y']
Pq={v:Fraction(1,10),alpha:Fraction(1,2),beta:Fraction(1,20),rho:Fraction(3,20),rhoT:Fraction(1,20),k:Fraction(27,100),kT:Fraction(1,50),gamma:Fraction(9,10)}
def xrecI(T,Q,Pp):
    return (Q+(I(Pq[k])-I(Pq[kT])-Pp)/T-I(Pq[alpha]*Pq[beta])*(I(1)-T)-I(Pq[v]))/I(Pq[alpha])-I(Pq[rho])
XG=xrecI(TG,QG,PG); XB=xrecI(TBB,QBB,PG)
base={z:I(val) for z,val in Pq.items()}
envB={**base,t1:TBB,t2:TBB,q:QBB,p:PG,x1:XB,x2:XB}; HB,_=reduced_hessian(G,y,[x1,x2],W,envB)
envG={**base,t1:TG,t2:TG,q:QG,p:PG,x1:XG,x2:XG}; HG,zrg=reduced_hessian(Efull,[t1,t2,q,p],[x1,x2],W,envG)
HP,_=reduced_hessian(G,y,[p],pi,envG)
assert HB[0][0].sign()<0 and HB[0][1].sign()>0
assert HG[0][0].sign()<0 and HG[0][1].sign()<0
assert HP[0][0].sign()<0
assert zrg[0][3].sign()<0 and zrg[1][3].sign()<0
print('B3 H11,H12',HB[0][0].fp(),HB[0][1].fp())
print('G H11,H12',HG[0][0].fp(),HG[0][1].fp())
print('private price SOC',HP[0][0].fp(),'p_x',zrg[0][3].fp())
print('exact canonical sign reversal certificate: PASS')
