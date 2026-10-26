"""Analytic-identity and grid audit for restricted two-sided bilinear SPNE.

Restrictions and proof: STAGE_11_BILINEAR_GLOBAL_SPNE.md.
Finite grid checks supplement, and DO NOT replace, analytic global proof.
Standard library only.
"""
from fractions import Fraction as F

M=R=tI=tR=F(1)
mu=F(1,2)
beta=F(1,2)
A=c=F(1)
vI=vR=F(4)
FP=F(1,10)
FG=F(1,50)

u=beta*beta*A*A/(2*c)
tau=beta*(1-beta)*A*A/c
S=u+tau
hI=tI/2
hR=tR/2
T=M*tI+R*tR
B=M*vI+R*vR
pI=hI-mu*R*tau
pR=hR-mu*M*u
Z=mu*M*R*S
P_duo=(M*pI+R*pR)/2
P_mon=B-T/4+Z
W_mon_use=B-T/12+Z
W_duo_use=B-T/48+Z/2
C=T/16-Z/2

assert mu*M<=1 and mu*R<=1
assert hI*hR>mu*mu*M*R*S*S
assert pI>0 and pR>0
assert vI>=max(F(3,4)*tI,pI+F(9,16)*tI)
assert vR>=max(F(3,4)*tR,pR+F(9,16)*tR)
assert (u,tau,S)==(F(1,8),F(1,4),F(3,8))
assert (pI,pR)==(F(3,8),F(7,16))
assert P_duo==F(13,32)
assert P_mon==F(123,16)
assert C==F(1,32)
assert W_mon_use==B+F(1,48)
assert W_duo_use==B+F(5,96)
assert P_duo-C==3*T/16

def private_upper_bound(x,y):
    return (2*M*hI*x*(1-x)+2*R*hR*y*(1-y)
            -Z*(x+y)+2*Z*x*y)

def private_gap(x,y):
    d,e=x-F(1,2), y-F(1,2)
    return 2*M*hI*d*d+2*R*hR*e*e-2*Z*d*e

def public_gap(x,y):
    d,e=x-F(1,2), y-F(1,2)
    return M*hI*d*d+R*hR*e*e-2*Z*d*e

assert private_upper_bound(F(1,2),F(1,2))==P_duo
for ii in range(51):
    for jj in range(51):
        x,y=F(ii,50),F(jj,50)
        assert P_duo-private_upper_bound(x,y)==private_gap(x,y)
        assert private_gap(x,y)>=0
        assert public_gap(x,y)>=0

def clip(a):
    return max(0., min(1., a))

def shares_under_P_deviation(pPI,pPR):
    x=y=.5
    wI=float(mu*R*u/hI)
    wR=float(mu*M*tau/hR)
    for _ in range(100):
        xx=clip(.5-(pPI-float(pI))/(2*float(hI))+wI*(y-.5))
        yy=clip(.5-(pPR-float(pR))/(2*float(hR))+wR*(xx-.5))
        if max(abs(xx-x),abs(yy-y))<1e-13:
            return xx,yy
        x,y=xx,yy
    raise AssertionError('fixed point did not converge')

best=(-1,None)
for i in range(161):
    for j in range(161):
        pi,pj=i/40,j/40
        x,y=shares_under_P_deviation(pi,pj)
        profit=float(M)*pi*x+float(R)*pj*y
        assert profit<=float(P_duo)+1e-9,(pi,pj,x,y,profit)
        if profit>best[0]:
            best=profit,(pi,pj,x,y)

def monopoly_revenue_envelope(x,y):
    return (M*x*(vI-tI*x*x/4+u*mu*R*y)
           +R*y*(vR-tR*y*y/4+tau*mu*M*x))
for ii in range(51):
    for jj in range(51):
        assert monopoly_revenue_envelope(F(ii,50),F(jj,50))<=P_mon

def policy(fp,fg):
    without=(W_mon_use-fp) if P_mon>=fp else F(0)
    enters_with=(P_duo>=fp)
    with_public=(W_duo_use-fp-fg) if enters_with else (W_mon_use-fg)
    return {'private_without_G':P_mon>=fp,
            'private_with_G':enters_with,
            'W_without_G':without,'W_with_G':with_public,
            'build':with_public>=without}

out=policy(FP,FG)
assert out['private_without_G'] and out['private_with_G'] and out['build']
assert out['W_with_G']-out['W_without_G']==C-FG==F(9,800)
assert policy(F(1,2),FG)['build'] and not policy(F(1,2),FG)['private_with_G']
assert policy(F(8),F(1,2))['build'] and not policy(F(8),F(1,2))['private_without_G']

print('u, tau, S:',u,tau,S)
print('duopoly fees (I,R):',pI,pR)
print('private gross with/without public:',P_duo,P_mon)
print('gross user welfare mixed/single:',W_duo_use,W_mon_use)
print('private entry social increment:',C)
print('private two-price 161x161 grid best:',best)
print('welfare and entry:',out)
print('PASS: exact rational identities, corners and 25,921 price deviations')
