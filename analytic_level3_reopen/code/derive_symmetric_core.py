from __future__ import annotations
import sympy as sp

alpha,beta,v,rhoT,k,kT=sp.symbols('alpha beta v rhoT k kT', positive=True)
t,q,p=sp.symbols('t q p', finite=True)
delta=alpha*beta
a=kT+p
d=k-a
AT=v+alpha*rhoT

Pq=sp.expand(q**2-(AT+2*delta*t)*q+2*delta*a)
delta_inv=sp.factor(q*(q-AT)/(2*(t*q-a)))
assert sp.simplify(Pq.subs(delta,delta_inv)) == 0

u=d/t-delta*t
w=1-2*delta*a/q**2
H=sp.factor(u*w-2*delta*t)
J=sp.Matrix([[u,0,-t],[0,u,-t],[-delta,-delta,w]])
assert sp.factor(J.det()-u*H) == 0

q_p=sp.factor(-2*delta*(q+u)/(q*H))
t_p=sp.factor((t*q_p-1)/u)
s=a/q
s_p=sp.factor(1/q-a*q_p/q**2)
R=sp.together(2*(t-s)+2*p*(t_p-s_p))
Rnum=sp.fraction(R)[0]
Rrem=sp.factor(sp.rem(sp.expand(Rnum),Pq,q))
assert sp.degree(Rrem,q) == 1
assert sp.degree(Rrem,t) == 5
assert sp.degree(Rrem,p) == 3

Res=sp.factor(sp.resultant(Pq,Rrem,q))
fac=sp.factor_list(Res)[1]
degrees=sorted((sp.degree(f,t),sp.degree(f,p),e) for f,e in fac)
assert any(dp==4 and dt==7 for dt,dp,e in degrees)
assert any(sp.factor(f-(kT+p))==0 and e==2 for f,e in fac)

print('symmetric participation degree(q)=',sp.degree(Pq,q))
print('delta inversion=',delta_inv)
print('det J = u*H with H=',H)
print('private FOC remainder degrees (q,t,p)=',
      (sp.degree(Rrem,q),sp.degree(Rrem,t),sp.degree(Rrem,p)))
print('private resultant factor degrees=',degrees)
print('L3-1 symmetric core: PASS')
