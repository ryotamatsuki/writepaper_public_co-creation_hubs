from __future__ import annotations
import sympy as sp


def build():
    t1,t2,q,p,x1,x2=sp.symbols('t1 t2 q p x1 x2')
    t,x=sp.symbols('t x')
    v,alpha,beta,rho,rhoT,k,kT,gamma=sp.symbols('v alpha beta rho rhoT k kT gamma')
    delta=alpha*beta; a=kT+p; d=k-a; AT=v+alpha*rhoT
    q1=v+alpha*(rho+x1)+delta*(1-t1); q2=v+alpha*(rho+x2)+delta*(1-t2)
    G1=sp.expand(t1*(q1-q)-d); G2=sp.expand(t2*(q2-q)-d)
    G3=sp.expand(q**2-(AT+delta*(t1+t2))*q+2*delta*a)
    G=sp.Matrix([G1,G2,G3]); y=[t1,t2,q]; J=G.jacobian(y)
    nT=t1+t2-2*a/q; pi=p*nT
    Gp=G.diff(p); piy=sp.Matrix([sp.diff(pi,z) for z in y]); pip=sp.diff(pi,p)
    Rnum=sp.together(pip-(piy.T*J.adjugate()*Gp)[0]/J.det()).as_numer_denom()[0]
    s=a/q; b1=rho+x1+beta*(1-t1); b2=rho+x2+beta*(1-t2); bT=rhoT+beta*nT
    PS=q*(t1**2-s**2)/2-a*(t1-s)+q1*(1-t1**2)/2-k*(1-t1)
    W=PS+(b1**2+b2**2+bT**2)/4-gamma*x1**2/2
    P={v:sp.Rational(1,10),alpha:sp.Rational(1,2),beta:sp.Rational(1,20),rho:sp.Rational(3,20),rhoT:sp.Rational(1,20),k:sp.Rational(27,100),kT:sp.Rational(1,50),gamma:sp.Rational(9,10)}
    Rsy=sp.factor(Rnum.subs({t1:t,t2:t,x1:x,x2:x}))
    Rcore=max((f for f,e in sp.factor_list(Rsy)[1]),key=sp.count_ops)
    G1s=G1.subs({t1:t,x1:x}); G3s=G3.subs({t1:t,t2:t})
    Es=sp.Matrix([G1s,G3s,Rcore]); z=[t,q,p]; A=Es.jacobian(z); bx=Es.diff(x)
    Ac=A.subs(P); bxc=bx.subs(P); D=sp.expand(Ac.det()); zp=[]
    for j in range(3):
        Aj=Ac.copy(); Aj[:,j]=-bxc; zp.append(sp.cancel(Aj.det()/D))
    u=sp.diff(G1,t1).subs({t1:t,x1:x,**P}); ta=sp.cancel(-P[alpha]*t/u)
    zt1=sp.cancel((zp[0]+ta)/2); zt2=sp.cancel((zp[0]-ta)/2); zq=sp.cancel(zp[1]/2); zpp=sp.cancel(zp[2]/2)
    sub={t1:t,t2:t,x1:x,x2:x,**P}
    AG=sp.cancel(sp.diff(W,x1).subs(sub)+sp.diff(W,t1).subs(sub)*zt1+sp.diff(W,t2).subs(sub)*zt2+sp.diff(W,q).subs(sub)*zq+sp.diff(W,p).subs(sub)*zpp)
    xrec=sp.cancel((q+(P[k]-P[kT]-p)/t-P[alpha]*P[beta]*(1-t)-P[v])/P[alpha]-P[rho])
    AGnum=sp.together(AG.subs(x,xrec)).as_numer_denom()[0]
    AGnum=sp.primitive(sp.Poly(AGnum,t,q,p))[1].as_expr()
    Rred=sp.together(Rcore.subs(P).subs(x,xrec)).as_numer_denom()[0]
    Rred=sp.primitive(sp.Poly(Rred,t,q,p))[1].as_expr()
    H=sp.primitive(sp.Poly(G3s.subs(P),t,q,p))[1].as_expr()
    E2=sp.Matrix([G1s.subs(P),G3s.subs(P)]); A2=E2.jacobian([t,q]); bx2=E2.diff(x); D2=sp.expand(A2.det()); zp2=[]
    for j in range(2):
        Aj=A2.copy(); Aj[:,j]=-bx2; zp2.append(sp.cancel(Aj.det()/D2))
    zt1b=sp.cancel((zp2[0]+ta)/2); zt2b=sp.cancel((zp2[0]-ta)/2); zqb=sp.cancel(zp2[1]/2)
    AB=sp.cancel(sp.diff(W,x1).subs(sub)+sp.diff(W,t1).subs(sub)*zt1b+sp.diff(W,t2).subs(sub)*zt2b+sp.diff(W,q).subs(sub)*zqb)
    Bnum=sp.together(AB.subs(x,xrec)).as_numer_denom()[0]
    Bnum=sp.primitive(sp.Poly(Bnum,t,q,p))[1].as_expr()
    Efull=sp.Matrix([G1,G2,G3,Rnum])
    return locals()
