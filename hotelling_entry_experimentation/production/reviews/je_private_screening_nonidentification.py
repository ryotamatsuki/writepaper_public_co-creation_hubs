#!/usr/bin/env python3
"""Constructive model-internal non-identification by private matching observations.

A new independent exact-Fraction proof, NOT a general prior-art novelty
certificate or a global four-party equilibrium result.
"""
from fractions import Fraction as F

chi=F(9,10)
n_B=2
y_B=F(20,73)
x_B_standard=F(55,219)
x_B_adverse=F(2,5)

def product_market(x):
    # Maintain the SAME new-venture Cournot output and profit while adjusting
    # incumbent output through market net-demand intercepts.
    A=(n_B+1)*x+chi*y_B
    D=2*y_B+chi*n_B*x
    den=2*(n_B+1)-n_B*chi**2
    assert (2*A-chi*D)/den==x
    assert ((n_B+1)*D-n_B*chi*A)/den==y_B
    x0=A/(n_B+1)
    cs1=((n_B*x)**2+2*chi*n_B*x*y_B+y_B*y_B)/2
    cs0=(n_B*x0)**2/2
    m=cs1+n_B*x*x-cs0-n_B*x0*x0
    assert m==-4*(219*x-14)/5329
    return A,D,m

den_A=2*(5+1)-5*chi**2
y_A=((5+1)-5*chi)/den_A
PI=y_A*y_A+y_B*y_B
zeta=F(1,200)
k=F(1,2500)
p0=F(1,20)
pb=p0+zeta*zeta*PI/k
R0=p0*PI+zeta*zeta*PI*PI/(2*k)
Ni,Nj=F(1),F(3)
gamma=F(1,10)
omega=F(1,5)
u=F(7,25)
d=F(1)
C=gamma*Ni*Nj
U=u+gamma*Nj*omega*R0
eta=gamma*Nj*omega*PI/2

def ell(t):
    DD=d+eta*t
    return (DD-U)/(2*DD)

def ellprime(t):
    DD=d+eta*t
    return U*eta/(2*DD*DD)

l0,l1=ell(F(0)),ell(F(1,2))
e0,e1=ellprime(F(0)),ellprime(F(1,2))
assert F(0)<l0<l1<F(1)
assert F(0)<e1<e0

aA,dA,M1=product_market(x_B_standard)
aB,dB,M2=product_market(x_B_adverse)
assert aA==dA==F(1)
assert (aB,dB)==(F(528,365),F(2314,1825))
assert M2<M1<0

def incidence(M):
    return M+omega*PI,pb*M+omega*R0

K1,S1=incidence(M1)
K2,S2=incidence(M2)
assert K2<K1<0 and S2<S1<0
# Verify two SCOPED derivative sign bounds across all t in [0, 1/2]
# via exact rational endpoint monotonicities. The bracket is positive.
Bmin=u+l0+C*S1+C*F(1,2)*K1*l1/2
assert Bmin>0
lower=e1*Bmin+C*K1*l1*l1/4
# For adverse market, S2<0, K2<0, so bracket <= u+l1
upper=e0*(u+l1)+C*K2*l0*l0/4
assert lower>F(1,10000),lower
assert upper<-F(1,100000),upper

def slope(M,t):
    K,S=incidence(M)
    z=ell(t)
    return (u+z+C*(S+t*K*z/2))*ellprime(t)+C*K*z*z/4

def ceiling(M,t):
    K,S=incidence(M)
    z=ell(t)
    return u*z-d/8+d*z*z/2+C*(S*z+t*K*z*z/4)

for t in (F(0),F(1,4),F(1,2)):
    assert slope(M1,t)>=lower
    assert slope(M2,t)<=upper
    # same P-only posted fee and admission cutoff in both economies
    private_fee=(U+d+eta*t)/2
    assert private_fee>0 and ell(t)>0

print("PASS: exact same venture output/profit, baseline research and all private fee/participation primitives")
print("PASS: A_B,D_B worlds",str(aA),str(dA),"vs",str(aB),str(dB))
print("PASS: regional incidence M_B",str(M1),"vs",str(M2))
print("PASS: all-lambda favorable lower",str(lower),float(lower))
print("PASS: all-lambda adverse upper",str(upper),float(upper))
print("PASS: screening-only statistics do not identify conditional public-value slope sign")
print("Scope: conditional region-B comparison, partner A fully participating; NOT the entire stage-1 equilibrium")
