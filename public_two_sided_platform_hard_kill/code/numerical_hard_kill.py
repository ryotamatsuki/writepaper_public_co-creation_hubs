"""Stage 4R-TP deterministic sign-reversal certificate.
Requires numpy/scipy. Seed for local perturbations: 20260903.
The script solves the smooth 0->HT->Hi regime only and rejects invalid draws.
"""
import numpy as np
from scipy.optimize import root, minimize_scalar, brentq

BASE=dict(v=.1,alpha=.5,beta=.05,rho=.15,rhoT=.05,k=.27,kT=.02,tau=.05,gamma=.9)

def eq(x1,x2,p,P):
    de=P['alpha']*P['beta']; a=P['kT']+p; d=P['k']-a
    if d<=0:return None
    A1=P['v']+P['alpha']*(P['rho']+x1); A2=P['v']+P['alpha']*(P['rho']+x2); AT=P['v']+P['alpha']*P['rhoT']
    def F(y):
        t1,t2,qT=y; q1=A1+de*(1-t1); q2=A2+de*(1-t2); s=a/qT
        return [t1*(q1-qT)-d,t2*(q2-qT)-d,qT-AT-de*(t1+t2-2*s)]
    s=root(F,[.6,.6,.14])
    if not s.success:return None
    t1,t2,qT=s.x; z0=a/qT; q1=A1+de*(1-t1);q2=A2+de*(1-t2)
    n1=1-t1;n2=1-t2;nT=t1+t2-2*z0
    b1=P['rho']+x1+P['beta']*n1;b2=P['rho']+x2+P['beta']*n2;bT=P['rhoT']+P['beta']*nT
    if not (0<z0<t1<1 and 0<z0<t2<1 and q1>qT and q2>qT and min(n1,n2,nT)>0 and all(0<b<1 for b in [b1,b2,bT])):return None
    return dict(t1=t1,t2=t2,qT=qT,z0=z0,q1=q1,q2=q2,n1=n1,n2=n2,nT=nT,b1=b1,b2=b2,bT=bT)

def profit(p,x1,x2,P):
    e=eq(x1,x2,p,P); return -1e9 if e is None else p*e['nT']

def pstar(x1,x2,P):
    hi=P['k']-P['kT']-1e-5; grid=np.linspace(1e-6,hi,30); vals=[profit(p,x1,x2,P) for p in grid]; j=int(np.argmax(vals)); lo=grid[max(0,j-1)]; hh=grid[min(len(grid)-1,j+1)]
    r=minimize_scalar(lambda p:-profit(p,x1,x2,P),bounds=(lo,hh),method='bounded',options={'xatol':1e-9})
    return r.x

def W(i,x1,x2,p,P):
    e=eq(x1,x2,p,P)
    if e is None:return np.nan
    t=e['t1'] if i==1 else e['t2']; q=e['q1'] if i==1 else e['q2']; a=P['kT']+p; s=e['z0']; qT=e['qT']
    ps=qT*(t*t-s*s)/2-a*(t-s)+q*(1-t*t)/2-P['k']*(1-t)
    partner=.25*(e['b1']**2+e['b2']**2+e['bT']**2); x=x1 if i==1 else x2
    return ps+partner-P['gamma']*x*x/2

def WG(i,x1,x2,P): return W(i,x1,x2,pstar(x1,x2,P),P)

def derivs(i,x,P,mode='G',pbar=None,h=1.5e-4):
    f=(lambda a,b:WG(i,a,b,P)) if mode=='G' else (lambda a,b:W(i,a,b,pbar,P))
    f0=f(x,x); xp=f(x+h,x); xm=f(x-h,x); pp=f(x+h,x+h);pm=f(x+h,x-h);mp=f(x-h,x+h);mm=f(x-h,x-h)
    gx=(xp-xm)/(2*h); hxx=(xp-2*f0+xm)/h**2; hxy=(pp-pm-mp+mm)/(4*h*h)
    return gx,hxx,hxy

def stationary(P,mode,pbar,bracket):
    f=lambda x:derivs(1,x,P,mode,pbar,2e-4)[0]
    x=brentq(f,*bracket,xtol=2e-5); d=derivs(1,x,P,mode,pbar)
    return x,d,-d[2]/d[1]

if __name__=='__main__':
    xg,dg,sg=stationary(BASE,'G',None,(.62,.75)); pg=pstar(xg,xg,BASE)
    xb,db,sb=stationary(BASE,'B3',pg,(.60,.72))
    print('G',xg,pg,dg,sg); print('B3',xb,db,sb)
    assert dg[1]<0 and db[1]<0 and sg<0<sb
    rng=np.random.default_rng(20260903); ok=0
    for _ in range(20):
        P=BASE.copy()
        for key in ['v','alpha','beta','rho','rhoT','k','kT','gamma']: P[key]*=1+rng.uniform(-.005,.005)
        try:
            xg,dg,sg=stationary(P,'G',None,(.62,.75)); pg=pstar(xg,xg,P); xb,db,sb=stationary(P,'B3',pg,(.60,.72)); ok+=int(dg[1]<0 and db[1]<0 and sg<0<sb)
        except Exception: pass
    print('local reversal successes',ok,'/20'); assert ok==20
