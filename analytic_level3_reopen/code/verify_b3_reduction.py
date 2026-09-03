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
yx1=-J.inv()*F.diff(x1)
yx2=-J.inv()*F.diff(x2)

s=a/q
n1=1-t1;n2=1-t2;nT=t1+t2-2*s
b1=rho+x1+beta*n1;b2=rho+x2+beta*n2;bT=rhoT+beta*nT
PS=q*(t1**2-s**2)/2-a*(t1-s)+q1*(1-t1**2)/2-k*(1-t1)
W=sp.expand(PS+(b1**2+b2**2+bT**2)/4-gamma*x1**2/2)

def D(expr,var):
    r=yx1 if var==x1 else yx2
    return sp.diff(expr,var)+(sp.Matrix([sp.diff(expr,z) for z in y]).T*r)[0]

FOC=sp.together(D(W,x1))
CROSS=sp.together(D(FOC,x2))
OWN=sp.together(D(FOC,x1))

sym={x1:x,x2:x,t1:t,t2:t}
xrec=sp.simplify((q+d/t-delta*(1-t)-v)/alpha-rho)
Pq=sp.expand(q**2-(AT+2*delta*t)*q+2*delta*a)

def reduced_num(expr):
    e=sp.together(expr.subs(sym))
    n=sp.fraction(e)[0]
    n=sp.fraction(sp.together(n.subs(x,xrec)))[0]
    return sp.factor(sp.rem(n,Pq,q))

Fnum=reduced_num(FOC)
Cnum=reduced_num(CROSS)
Onum=reduced_num(OWN)

assert sp.degree(Fnum,q)==1 and sp.degree(Fnum,t)==8 and sp.degree(Fnum,p)==4
assert sp.degree(Cnum,q)==1 and sp.degree(Cnum,t)==18 and sp.degree(Cnum,p)==8
assert sp.degree(Onum,q)==1 and sp.degree(Onum,t)==18 and sp.degree(Onum,p)==9

Btp=sp.factor(sp.resultant(Pq,Fnum,q))
factors=sp.factor_list(Btp)[1]
core=[f for f,e in factors if sp.degree(f,t)==14][0]
assert sp.degree(core,t)==14 and sp.degree(core,p)==7

print('B3 FOC reduced degrees (q,t,p)=',(1,8,4))
print('B3 cross numerator reduced degrees (q,t,p)=',(1,18,8))
print('B3 own numerator reduced degrees (q,t,p)=',(1,18,9))
print('B3 equilibrium core degrees (t,p)=',(14,7))
print('B3 equilibrium-manifold reduction: PASS')
