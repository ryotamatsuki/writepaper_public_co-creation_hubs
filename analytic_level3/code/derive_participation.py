from __future__ import annotations
import sympy as sp


def symbolic_objects():
    alpha,beta,v,rho,rhoT,k,kT,p = sp.symbols(
        'alpha beta v rho rhoT k kT p', positive=True, finite=True
    )
    t1,t2,qT,x1,x2 = sp.symbols('t1 t2 qT x1 x2', positive=True, finite=True)
    delta=alpha*beta
    a=kT+p
    d=k-a
    AT=v+alpha*rhoT
    A1=v+alpha*(rho+x1)
    A2=v+alpha*(rho+x2)
    q1=A1+delta*(1-t1)
    q2=A2+delta*(1-t2)
    F1=sp.expand(t1*(q1-qT)-d)
    F2=sp.expand(t2*(q2-qT)-d)
    FT=sp.expand(qT-AT-delta*(t1+t2-2*a/qT))
    q_poly=sp.expand(qT**2-(AT+delta*(t1+t2))*qT+2*delta*a)
    x1_rec=sp.simplify((qT+d/t1-delta*(1-t1)-v)/alpha-rho)
    x2_rec=sp.simplify((qT+d/t2-delta*(1-t2)-v)/alpha-rho)
    return locals()


def verify():
    o=symbolic_objects()
    assert sp.simplify(o['qT']*o['FT']-o['q_poly']) == 0
    assert sp.simplify(o['F1'].subs(o['x1'],o['x1_rec'])) == 0
    assert sp.simplify(o['F2'].subs(o['x2'],o['x2_rec'])) == 0

    t=sp.symbols('t', positive=True, finite=True)
    u,w=sp.symbols('u w', finite=True)
    J=sp.Matrix([[u,0,-t],[0,u,-t],[-o['delta'],-o['delta'],w]])
    H=sp.expand(u*w-2*o['delta']*t)
    assert sp.factor(J.det()-u*H) == 0
    return True


if __name__ == '__main__':
    verify()
    print('participation symbolic identities: PASS')
