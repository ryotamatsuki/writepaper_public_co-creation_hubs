from __future__ import annotations
import sympy as sp


def verify_b3_first_order():
    D1,D2,T,a,d,k,alpha,rhoT,eps=sp.symbols(
        'Delta1 Delta2 T a d k alpha rhoT eps', positive=True
    )
    t10=d/D1; t20=d/D2
    Q1=t10+t20-2*a/T
    qT=T+eps*Q1
    s=a/T-eps*a*Q1/T**2
    q1=T+D1+eps*(1-t10)
    t11=-d*((1-t10)-Q1)/D1**2
    t1=t10+eps*t11
    b1=rhoT+D1/alpha+eps/alpha*(1-t10)
    b2=rhoT+D2/alpha+eps/alpha*(1-t20)
    bT=rhoT+eps/alpha*Q1
    PS=qT*(t1**2-s**2)/2-a*(t1-s)+q1*(1-t1**2)/2-k*(1-t1)
    partner=(b1**2+b2**2+bT**2)/4
    W=PS+partner
    W_eps=sp.diff(W,eps).subs(eps,0)
    coeff=sp.factor(alpha**2*sp.diff(W_eps,D1,D2))
    target=alpha**2*d**3/(D1**3*D2**2)
    assert sp.simplify(coeff.subs(k,a+d)-target)==0
    return target


def verify_beta0_g_cross():
    alpha,Delta,T,k,kT=sp.symbols('alpha Delta T k kT', positive=True)
    p=(-Delta*kT+T*k-T*kT)/(2*(Delta+T))
    a=kT+p
    d=k-a
    p_i=-alpha*k*T/(4*(Delta+T)**2)
    p_ij=-alpha**2*k*T**2/(4*Delta*(Delta+T)**3)
    m=sp.factor(d/Delta-a/T)
    Wp=-m
    Wxp=d*alpha/Delta**2
    Wpp=(Delta+T)/(Delta*T)
    cross=sp.factor(Wxp*p_i+Wpp*p_i**2+Wp*p_ij)
    target=-3*T*alpha**2*k**2/(16*Delta*(Delta+T)**3)
    assert sp.simplify(cross-target)==0
    return target


if __name__=='__main__':
    print('d M_B3/d(delta)|0 =',verify_b3_first_order())
    print('M_G(beta=0) =',verify_beta0_g_cross())
    print('small-beta analytic identities: PASS')
