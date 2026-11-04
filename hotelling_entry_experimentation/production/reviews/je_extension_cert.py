#!/usr/bin/env python3
"""Exact conditional JE extension verifier. Not a four-player equilibrium proof."""
from fractions import Fraction as Q
import sympy as sp

u,d,ni,nj,gamma,omega,pi,market,return0,lamb=sp.symbols(
    'u d ni nj gamma omega pi market return0 lamb', positive=True)
s0=sp.symbols('s0',real=True)
mu=sp.Rational(1,2)
C=gamma*ni*nj
U=u+gamma*nj*omega*return0
eta=gamma*nj*omega*pi*mu
D=d+eta*lamb
ell=(D-U)/(2*D)
ell_frozen=(d-U)/(2*D)
T=lambda t: ni*(u*t-d/8+d*t*t/2)+C*(s0*t+lamb*(market+omega*pi)*t*t/4)
weak_opt=sp.factor(sp.diff(sp.diff(T(ell),lamb),gamma).subs(gamma,0))
weak_fix=sp.factor(sp.diff(sp.diff(T(ell_frozen),lamb),gamma).subs(gamma,0))
assert sp.simplify(weak_opt-ni*nj/(16*d*d)*(market*(d-u)**2+omega*pi*(d*d+3*u*u)))==0
assert sp.simplify(weak_opt-weak_fix-ni*nj*omega*pi*(d+u)/(8*d))==0

# General distribution reduces exactly to the original uniform equation.
f,fp,J,t,mu_j=sp.symbols('f fp J t mu_j',real=True)
B=2*d*f+(u+d*t)*fp
A=u*(u+d*t)*f*f/(d*B)
general=ni*nj*mu_j*(market*J+omega*pi*(J+A))
special=general.subs({fp:0,f:1,J:t*t/2,mu_j:sp.Rational(1,2)}).subs(t,(d-u)/(2*d))
assert sp.simplify(special-weak_opt)==0

def cournot(n):
    chi=Q(9,10)
    den=2*(n+1)-n*chi**2
    x=(2-chi)/den
    y=((n+1)-n*chi)/den
    x0=Q(1,n+1)
    cs1=((n*x)**2+2*chi*n*x*y+y*y)/2
    cs0=(n*x0)**2/2
    return cs1+n*x*x-cs0-n*x0*x0,y*y

ma,pa=cournot(5)
mb,pb=cournot(2)
PI=pa+pb
p0=Q(1,20);zeta=Q(1,200);k=Q(1,2500)
R0=p0*PI+zeta*zeta*PI*PI/(2*k)
PBASE=p0+zeta*zeta*PI/k
S0=PBASE*mb+Q(1,5)*R0
NI=Q(1);NJ=Q(3);GAMMA=Q(1,10);OMEGA=Q(1,5)
BASIC=Q(7,25);PREMIUM=Q(1);MEAN=Q(1,2)
CC=GAMMA*NI*NJ
UU=BASIC+GAMMA*NJ*OMEGA*R0
ETA=GAMMA*NJ*OMEGA*PI*MEAN
KK=mb+OMEGA*PI
assert KK<0 and S0<0

def uniform_ceiling(L,theta):
    return NI*(BASIC*L-PREMIUM/8+PREMIUM*L*L/2)+CC*(S0*L+theta*KK*L*L/4)
def opt(t):
    dd=PREMIUM+ETA*t
    return (dd-UU)/(2*dd)
def fixed(t):
    dd=PREMIUM+ETA*t
    return (PREMIUM-UU)/(2*dd)

optimal_0=uniform_ceiling(opt(Q(0)),Q(0))
optimal_half=uniform_ceiling(opt(Q(1,2)),Q(1,2))
frozen_half=uniform_ceiling(fixed(Q(1,2)),Q(1,2))
assert optimal_0==uniform_ceiling(fixed(Q(0)),Q(0))
assert optimal_half>optimal_0>frozen_half
assert abs(float(optimal_half-optimal_0)-0.00010659909532825)<1e-12
assert abs(float(frozen_half-optimal_0)+0.00042306040973083)<1e-12

# Exact-rational full lambda-domain Beta(2,2) conditional slope certificate.
def m(t):return 1-3*t*t+2*t**3
def g(t):return 6*t*(1-t)
def gp(t):return 6-12*t
def Jbeta(t):return 2*t**3-Q(3,2)*t**4
left,right=Q(3,10),Q(33,100)
Dmax=PREMIUM+ETA*Q(1,2)
assert m(left)-left*g(left)>0 and m(right)-right*g(right)>0
assert (PREMIUM*(m(left)-left*g(left))-UU*g(left))>0
assert (Dmax*(m(right)-right*g(right))-UU*g(right))<0
# Beta hazard derivative is 6(2t^2+1)/[(1-t)^2(1+2t)^2]>0.
assert gp(right)>0
upperden=Dmax*(2*Dmax*g(right)+(UU+Dmax*right)*gp(left))
lprime_lower=ETA*UU*g(left)/upperden
bracket_lower=NI*(BASIC+PREMIUM*left)+CC*S0+CC*Q(1,2)*KK*MEAN*right
assert bracket_lower>0
slope_lower=bracket_lower*g(left)*lprime_lower+CC*KK*MEAN*Jbeta(right)
assert slope_lower>Q(1,10000)

print('PASS: frozen-fee vs repriced sign algebra; exact original uniform reduction')
print('PASS: original fixed vs optimised delta:',float(frozen_half-optimal_0),float(optimal_half-optimal_0))
print('PASS: conditional Beta(2,2) Tprime >',float(slope_lower),'over [0,0.5]')
print('NOT CLAIMED: general F finite matching, global Beta off-path price equilibria, full SPNE, or Lean formalisation')
