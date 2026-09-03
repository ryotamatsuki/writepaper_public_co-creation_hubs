from __future__ import annotations
import sympy as sp

R=sp.Rational
v,alpha,beta,rho,rhoT,k,kT,gamma=R(1,10),R(1,2),R(1,20),R(3,20),R(1,20),R(27,100),R(1,50),R(9,10)
delta=alpha*beta
x1,x2,p=sp.symbols('x1 x2 p')
t1,t2,q=sp.symbols('t1 t2 q')
x,t=sp.symbols('x t')
a=kT+p; d=k-a
A1=v+alpha*(rho+x1); A2=v+alpha*(rho+x2); AT=v+alpha*rhoT
q1=A1+delta*(1-t1); q2=A2+delta*(1-t2)
F=sp.Matrix([
    t1*(q1-q)-d,
    t2*(q2-q)-d,
    q-AT-delta*(t1+t2-2*a/q),
])
y=sp.Matrix([t1,t2,q])
J=F.jacobian(y)
yp=-J.inv()*F.diff(p)

Pi=p*(t1+t2-2*a/q)
Pi_p=sp.diff(Pi,p)+(sp.Matrix([sp.diff(Pi,z) for z in y]).T*yp)[0]
Rnum=sp.fraction(sp.together(Pi_p))[0]

s=a/q
n1=1-t1; n2=1-t2; nT=t1+t2-2*s
b1=rho+x1+beta*n1; b2=rho+x2+beta*n2; bT=rhoT+beta*nT
PS=q*(t1**2-s**2)/2-a*(t1-s)+q1*(1-t1**2)/2-k*(1-t1)
W=sp.expand(PS+(b1**2+b2**2+bT**2)/4-gamma*x1**2/2)

E=sp.Matrix([F[0],F[1],F[2],Rnum])
z=sp.Matrix([t1,t2,q,p])
sym={x1:x,x2:x,t1:t,t2:t}
Ez=E.jacobian(z).subs(sym)
Ex1=E.diff(x1).subs(sym)
zx1=-Ez.inv()*Ex1
Wz=sp.Matrix([sp.diff(W,zz) for zz in z]).subs(sym)
Wx=sp.diff(W,x1).subs(sym)
FOCG=sp.together(Wx+(Wz.T*zx1)[0])
Gnum=sp.fraction(FOCG)[0]

xrec=sp.simplify((q+d/t-delta*(1-t)-v)/alpha-rho)
Pq=sp.expand(q**2-(AT+2*delta*t)*q+2*delta*a)

Rs=sp.factor(sp.together((Rnum.subs(sym)).subs(x,xrec)))
Rq=sp.factor(sp.rem(sp.fraction(Rs)[0],Pq,q))
Priv=sp.factor(sp.resultant(Pq,Rq,q))
priv_factors=sp.factor_list(Priv)[1]
Priv_core=[f for f,e in priv_factors if sp.degree(f,p)==4 and sp.degree(f,t)==7][0]

Gsub=sp.together(Gnum.subs(x,xrec))
Gqrem=sp.rem(sp.fraction(Gsub)[0],Pq,q)
assert sp.degree(Rq,q)==1 and sp.degree(Gqrem,q)==1
Gtp=sp.factor(sp.resultant(Rq,Gqrem,q))
Gred=sp.rem(sp.Poly(sp.expand(Gtp),p),sp.Poly(sp.expand(Priv_core),p)).as_expr()
if Gred==0:
    raise AssertionError('unexpected zero G remainder')

Rest=sp.resultant(sp.expand(Priv_core),sp.expand(Gred),p)
tfacs=sp.factor_list(Rest)[1]
t_degrees=sorted((sp.degree(f,t),e) for f,e in tfacs)
T31=[f for f,e in tfacs if sp.degree(f,t)==31][0]
PT=sp.Poly(T31,t,domain=sp.QQ)
lo,hi=R(58489,100000),R(58490,100000)
assert PT.count_roots(lo,hi)==1
assert sp.sign(PT.eval(lo))*sp.sign(PT.eval(hi))==-1

Resp=sp.resultant(sp.expand(Priv_core),sp.expand(Gred),t)
pfacs=sp.factor_list(Resp)[1]
P31=[f for f,e in pfacs if sp.degree(f,p)==31][0]
PP=sp.Poly(P31,p,domain=sp.QQ)
plo,phi=R(22746,1000000),R(22747,1000000)
assert PP.count_roots(plo,phi)==1
assert sp.sign(PP.eval(plo))*sp.sign(PP.eval(phi))==-1

print('canonical private core degrees (t,p)=',(sp.degree(Priv_core,t),sp.degree(Priv_core,p)))
print('G t-resultant factor degrees=',t_degrees)
print('G witness t interval=',(lo,hi),'unique')
print('G witness p interval=',(plo,phi),'unique')
print('canonical exact elimination/root geometry: PASS')
