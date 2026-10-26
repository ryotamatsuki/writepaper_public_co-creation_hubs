"""Stage 9 regression: exact ex-post payoff accounting, and two distinct
endogenous participation equilibria under the same one-hub pricing history.
NOT a complete SPNE solver. Requires Python standard library only.
"""
from math import exp, sqrt, expm1

A=c=1.
beta=rho=.5
d2=1.8
kappa=.02
lam=4.
v_R=.02
p_R=.10
t_R=.10

def bargaining(m):
    scale=1 if m==1 else d2
    rm=A*(m**rho-m*(1-beta)*(m**rho-(m-1)**rho))
    effort=rm/(c*scale)
    private=rm*rm/(2*c*scale)
    each_partner=(1-beta)*A*(m**rho-(m-1)**rho)*effort
    social=A*m**rho*effort-c*scale*effort*effort/2
    assert abs(social-private-m*each_partner)<1e-12
    return effort,private,each_partner,social

e1,u1,tau1,s1=bargaining(1)
e2,u2,tau2,s2=bargaining(2)
assert abs(u1-.125)<1e-12 and abs(tau1-.25)<1e-12
assert abs(s1-.375)<1e-12

def q(r):
    if r<=0:return 0.
    h=-expm1(-lam*r)
    return -r*expm1(-h/r)

def partner_receipt(r):
    return tau1*q(r)/r if r>0 else 0.

def participation_best_response(r):
    surplus=partner_receipt(r)+v_R-p_R
    return min(1.,sqrt(max(0.,surplus/t_R)))

def f(r):return participation_best_response(r)-r

# Positive equilibrium by IVT: f(.90)>0, f(.95)<0.
lo,hi=.90,.95
assert f(lo)>0 and f(hi)<0
for _ in range(80):
    mid=(lo+hi)/2
    if f(mid)>0:lo=mid
    else:hi=mid
r=(lo+hi)/2
assert abs(f(r))<1e-12
assert q(r)*u1-kappa>0
assert abs(v_R-p_R-t_R*r*r+partner_receipt(r))<1e-12

# Zero equilibrium: all innovators select k=0, no proposals; no
# positive partner profit exists from a unilateral participation move.
assert v_R-p_R<0 and kappa>0

print("m=1 effort, innovator net income, per-partner payment, net social value:",
      e1,u1,tau1,s1)
print("m=2:",e2,u2,tau2,s2)
print("EqA: N_R=0; k=0; completed matches=0")
print("EqB: N_R=%.12f; q=%.12f; matches=%.12f; partner receipt=%.12f; innovator net contact gain=%.12f"%(
  r,q(r),q(r),partner_receipt(r),q(r)*u1-kappa))
print("PASS: bargaining accounting and two distinct equilibrium fixed points")
