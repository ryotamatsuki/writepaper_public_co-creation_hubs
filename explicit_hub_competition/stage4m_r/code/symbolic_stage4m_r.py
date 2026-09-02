import sympy as sp
from math import comb
import itertools

p,h,m,l,a,sigma=sp.symbols("p h m l a sigma", positive=True)
q=1-p

muS=sp.simplify(p*(h+m)/2+p*q*l)
nuS=sp.simplify(p*(h**2+m**2)/2+p*q*l**2)
muB=sp.simplify(p*h+p*q*l+p*q**2*m)
nuB=sp.simplify(p*h**2+p*q*l**2+p*q**2*m**2)
muB0=sp.simplify(p*h+p*q*l)
nuB0=sp.simplify(p*h**2+p*q*l**2)

def K(mu,nu):
    return sp.Rational(2,9)*a**2+sp.Rational(4,9)*a*mu+sp.Rational(11,18)*nu-sp.Rational(7,18)*mu**2

Gamma=sp.factor(K(muB,nuB)-2*K(muS,nuS)+K(0,0)+sigma*(muB-2*muS))
Gamma0=sp.factor(K(muB0,nuB0)-2*K(muS,nuS)+K(0,0)+sigma*(muB0-2*muS))
dmu=sp.factor(muB-2*muS)
dnu=sp.factor(nuB-2*nuS)
dmu0=sp.factor(muB0-2*muS)
dnu0=sp.factor(nuB0-2*nuS)

assert sp.simplify(dmu + p*(q*l+p*(2-p)*m)) == 0
assert sp.simplify(dnu + p*(q*l**2+p*(2-p)*m**2)) == 0
assert sp.simplify(dmu0 + p*(q*l+m)) == 0
assert sp.simplify(dnu0 + p*(q*l**2+m**2)) == 0
assert sp.factor(sp.diff(Gamma,a)-sp.Rational(4,9)*dmu)==0
assert sp.factor(sp.diff(Gamma,sigma)-dmu)==0

Q,X,Y,R=sp.symbols("Q X Y R", nonnegative=True)
Gnorm=sp.factor(Gamma.subs({h:1,m:X,l:Y,p:1-Q,a:1,sigma:0}))
P=sp.factor(-36*Gnorm/(1-Q))
Pcube=sp.expand(P.subs(X,Y+(1-Y)*R))

G0norm=sp.factor(Gamma0.subs({h:1,m:X,l:Y,p:1-Q,a:1,sigma:0}))
P0=sp.factor(-36*G0norm/(1-Q))
P0cube=sp.expand(P0.subs(X,Y+(1-Y)*R))

def bernstein_coeffs(poly_expr, vars_):
    poly=sp.Poly(poly_expr,*vars_)
    degs=tuple(poly.degree(v) for v in vars_)
    coeff={mon:sp.Rational(c) for mon,c in poly.terms()}
    out={}
    for ks in itertools.product(*[range(d+1) for d in degs]):
        s=sp.Rational(0)
        for inds in itertools.product(*[range(k+1) for k in ks]):
            ac=coeff.get(tuple(inds),0)
            if ac:
                fac=sp.Rational(1)
                for k,i,n in zip(ks,inds,degs):
                    fac*=sp.Rational(comb(k,i),comb(n,i))
                s += ac*fac
        out[ks]=sp.simplify(s)
    return degs,out

deg,bc=bernstein_coeffs(Pcube,(Q,R,Y))
deg0,bc0=bernstein_coeffs(P0cube,(Q,R,Y))
assert all(v>=0 for v in bc.values())
assert all(v>=0 for v in bc0.values())
assert any(v>0 for v in bc.values())
assert any(v>0 for v in bc0.values())

print("muS =", muS)
print("muB =", muB)
print("Gamma =", Gamma)
print("Gamma0 =", Gamma0)
print("Bernstein referral degree",deg,"zeros",sum(v==0 for v in bc.values()),"min positive",min(v for v in bc.values() if v>0))
print("Bernstein no-ref degree",deg0,"zeros",sum(v==0 for v in bc0.values()),"min positive",min(v for v in bc0.values() if v>0))
