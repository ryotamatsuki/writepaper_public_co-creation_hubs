"""Stage 7R-TP deterministic B3-vs-G welfare reproduction."""
import numpy as np
from scipy.optimize import root,minimize_scalar,brentq
BASE=dict(v=.1,alpha=.5,beta=.05,rho=.15,rhoT=.05,k=.27,kT=.02,tau=.05,gamma=.9)
def eq(x1,x2,p,P):
    de=P['alpha']*P['beta']; a=P['kT']+p; d=P['k']-a
    A1=P['v']+P['alpha']*(P['rho']+x1); A2=P['v']+P['alpha']*(P['rho']+x2); AT=P['v']+P['alpha']*P['rhoT']
    def F(y):
        t1,t2,qT=y; q1=A1+de*(1-t1); q2=A2+de*(1-t2); s=a/qT
        return [t1*(q1-qT)-d,t2*(q2-qT)-d,qT-AT-de*(t1+t2-2*s)]
    r=root(F,[.6,.6,.14])
    if not r.success:return None
    t1,t2,qT=r.x; s=a/qT; q1=A1+de*(1-t1); q2=A2+de*(1-t2)
    n1=1-t1;n2=1-t2;nT=t1+t2-2*s
    b1=P['rho']+x1+P['beta']*n1;b2=P['rho']+x2+P['beta']*n2;bT=P['rhoT']+P['beta']*nT
    if not(0<s<t1<1 and 0<s<t2<1 and min(n1,n2,nT)>0):return None
    return dict(t1=t1,t2=t2,qT=qT,s=s,q1=q1,q2=q2,n1=n1,n2=n2,nT=nT,b1=b1,b2=b2,bT=bT)
def profit(p,x1,x2,P):
    e=eq(x1,x2,p,P); return -1e9 if e is None else p*e['nT']
def pstar(x1,x2,P):
    hi=P['k']-P['kT']-1e-5; grid=np.linspace(1e-6,hi,30); vals=[profit(p,x1,x2,P) for p in grid]
    j=int(np.argmax(vals)); lo=grid[max(0,j-1)]; hh=grid[min(len(grid)-1,j+1)]
    return minimize_scalar(lambda p:-profit(p,x1,x2,P),bounds=(lo,hh),method='bounded',options={'xatol':1e-9}).x
def W(i,x1,x2,p,P):
    e=eq(x1,x2,p,P); t=e['t1'] if i==1 else e['t2']; q=e['q1'] if i==1 else e['q2']
    a=P['kT']+p; s=e['s']; qT=e['qT']
    ps=qT*(t*t-s*s)/2-a*(t-s)+q*(1-t*t)/2-P['k']*(1-t)
    partner=.25*(e['b1']**2+e['b2']**2+e['bT']**2); x=x1 if i==1 else x2
    return ps+partner-P['gamma']*x*x/2
def gx(x,mode,pbar=None,h=2e-4):
    f=lambda a,b: W(1,a,b,pstar(a,b,BASE),BASE) if mode=='G' else W(1,a,b,pbar,BASE)
    return (f(x+h,x)-f(x-h,x))/(2*h)
def solve_symmetric():
    xg=brentq(lambda x:gx(x,'G'),.62,.75); pg=pstar(xg,xg,BASE); xb=brentq(lambda x:gx(x,'B3',pg),.60,.72)
    return xg,pg,xb
if __name__=='__main__':
    xg,pg,xb=solve_symmetric()
    for lab,x in [('G',xg),('B3',xb)]:
        e=eq(x,x,pg,BASE); wi=W(1,x,x,pg,BASE); pi=profit(pg,x,x,BASE)
        print(lab,dict(x=x,p=pg,nT=e['nT'],W_i=wi,Pi=pi,W_N=2*wi+pi))
