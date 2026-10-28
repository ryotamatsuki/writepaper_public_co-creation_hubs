#!/usr/bin/env python3
"""Exact-rational P/P and first-stage baseline certificate (review-only).

Frozen STAGE33-2026-10-28-THEORY-FREEZE-v1. No production-solver import.
The mathematical argument for existence/uniqueness is in the companion
STAGE11R_P0_PP_GLOBAL_CERTIFICATE.md. This script checks its sufficient
rational inequalities, 2D Poincare-Miranda faces, all first-stage signs.
"""
from fractions import Fraction as F
from decimal import Decimal
from itertools import product
import json

class I:
    def __init__(self,lo,hi=None):
        self.lo=F(lo);self.hi=F(lo if hi is None else hi)
        assert self.lo<=self.hi
    @staticmethod
    def wrap(v):return v if isinstance(v,I) else I(v)
    def __add__(self,v):
        v=I.wrap(v);return I(self.lo+v.lo,self.hi+v.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,v):return self+-I.wrap(v)
    def __rsub__(self,v):return I.wrap(v)+-self
    def __mul__(self,v):
        v=I.wrap(v);t=[self.lo*v.lo,self.lo*v.hi,self.hi*v.lo,self.hi*v.hi]
        return I(min(t),max(t))
    __rmul__=__mul__
    def __truediv__(self,v):
        v=I.wrap(v)
        assert v.lo>0 or v.hi<0
        return self*I(1/v.hi,1/v.lo)
    def __rtruediv__(self,v):return I.wrap(v)/self
    def sign(self):return '+' if self.lo>0 else '-' if self.hi<0 else '?'
    def approx(self):return [float(self.lo),float(self.hi)]

def primitives(lam,pooled):
    n=(5,2);chi=F(9,10)
    m=[];pi=F(0);base=[]
    for ni in n:
        den=2*(ni+1)-ni*chi*chi
        x=(2-chi)/den;y=(ni+1-ni*chi)/den;zero=F(1,ni+1)
        cs1=((ni*x)*(ni*x)+2*chi*ni*x*y+y*y)/2
        cs0=(ni*zero)*(ni*zero)/2
        m.append(cs1+ni*x*x-cs0-ni*zero*zero)
        pi+=y*y;base.append(cs0+ni*zero*zero)
    r0=F(1,20)*pi+F(5,1000)**2*pi*pi/(2*F(4,10000))
    pb=F(1,20)+F(5,1000)**2*pi/F(4,10000)
    nn=(F(3),F(1));omega=(F(4,5),F(1,5))
    aa=[F(2,5)*F(1,4)*nn[1-i]*omega[i]*r0 for i in range(2)]
    bb=[F(2,5)*F(1,4)*nn[1-i]*omega[i]*lam*pi for i in range(2)]
    if pooled:aa=[aa[i]+bb[i]/4 for i in range(2)];bb=[F(0),F(0)]
    return dict(m=m,pi=pi,base=base,r0=r0,pb=pb,a=aa,b=bb,
      u=(F(119,500),F(7,25)),fg=(F(9,100),F(81,2000)),
      fp=(F(13,20),F(6,25)),n=nn,omega=omega,lam=lam,pooled=pooled)

def bounds(p):
    u,a,b=p['u'],p['a'],p['b'];k=[a[i]+b[i] for i in range(2)]
    out=[]
    for i in range(2):
        z=k[0]*k[1]
        lo=(1-z-u[i]-a[i])/(2+b[i]/2-z)
        hi=(1+b[i]/2-u[i])/(2+b[i]/2)
        price_lo=u[i]+lo;price_hi=u[i]+a[i]+(1+b[i]/2)*hi
        h=b[i]+b[i]*k[1-i]+b[1-i]*k[i]+k[i]
        lips=(1+b[i]/2)*h/(2*(1-z))+k[i]
        assert 0<lo<hi<1 and z<1
        assert u[i]+a[i]<price_lo<price_hi<u[i]+1
        assert 0<=lips<F(1,100)
        out.append(dict(region='AB'[i],cutoff=[float(lo),float(hi)],
          fee=[float(price_lo),float(price_hi)],br_contraction_bound=float(lips)))
    return out

def stationary(x,p):
    u,a,b=p['u'],p['a'],p['b']
    d=[1+b[i]*(1-x[1-i]*x[1-i])/2 for i in range(2)]
    c=[a[i]+b[i]*x[0]*x[1] for i in range(2)]
    fee=[u[i]+a[i]*(1-x[1-i])+d[i]*x[i] for i in range(2)]
    return [(d[i]-c[0]*c[1]/d[1-i])*(1-x[i])-fee[i] for i in range(2)]

def faces(p,center,radius=F(1,10**9)):
    box=[[v-radius,v+radius] for v in center]
    evidence=[]
    for i in range(2):
        for side,idx,want in [('lower',0,'+'),('upper',1,'-')]:
            x=[I(box[j][0],box[j][1]) for j in range(2)]
            x[i]=I(box[i][idx])
            z=stationary(x,p)[i]
            assert z.sign()==want,('failed face',i,side,z.approx())
            evidence.append(dict(region='AB'[i],face=side,FOC=z.approx(),sign=z.sign()))
    return box,evidence

STATES=('0','G','P','B')
FLIP={'G':{'0':'G','G':'0','P':'B','B':'P'},'P':{'0':'P','P':'0','G':'B','B':'G'}}
def basic(i,state):
    return [F(15,100),F(1,10)][i]+F(2,5)*[F(22,100),F(45,100)][i]*int(state!='0')

def profile(s,p,box):
    a,b,n=p['a'],p['b'],p['n']
    fee=[I(0),I(0)];mq=[None,None]
    if s==('P','P'):
        x=[I(box[i][0],box[i][1]) for i in range(2)]
        for i in range(2):
            j=1-i
            fee[i]=I(p['u'][i])+a[i]*(1-x[j])+(1+b[i]*(1-x[j]*x[j])/2)*x[i]
            mq[i]=(1-x[i],(1-x[i]*x[i])/2)
    else:
        for i in range(2):
            if s[i]=='0':mq[i]=(I(0),I(0))
            if s[i] in ('G','B'):mq[i]=(I(1),I(F(1,2)))
        for i in range(2):
            if s[i]!='P':continue
            U=basic(i,s[1-i])+a[i]*mq[1-i][0]
            T=1+b[i]*mq[1-i][1]
            assert T.lo>U.hi
            fee[i]=(U+T)/2
            cut=(T-U)/(2*T)
            mq[i]=(1-cut,(1-cut*cut)/2)
        for i in range(2):
            if s[i]=='B':fee[i]=I(F(1,2))
    wealth=[];profit=[]
    for i in range(2):
        d=basic(i,s[1-i])
        s0=p['pb']*p['m'][i]+p['omega'][i]*p['r0']
        s1=p['lam']*(p['m'][i]+p['omega'][i]*p['pi'])
        if p['pooled']:s0+=s1/4;s1=F(0)
        z={'0':I(0),'G':n[i]*d,'P':n[i]*(d*mq[i][0]+mq[i][1]),
          'B':n[i]*(d+F(3,8))}[s[i]]
        w=I(p['base'][i])+z-p['fg'][i]*int(s[i] in ('G','B'))-p['fp'][i]*int(s[i] in ('P','B'))
        w+=F(3,10)*mq[0][0]*mq[1][0]*s0+F(3,10)*mq[0][1]*mq[1][1]*s1
        v=I(0) if s[i] in ('0','G') else n[i]*fee[i]*(mq[i][0] if s[i]=='P' else F(1,2))-p['fp'][i]
        wealth.append(w);profit.append(v)
    return dict(welfare=wealth,private=profit,fees=fee,mq=mq)

def profile_certificate(p,box,expected):
    prof={s:profile(s,p,box) for s in product(STATES,repeat=2)}
    losers=[];winner=[]
    for s in prof:
        slacks=[]
        for i in range(2):
            for act in ('G','P'):
                t=list(s);t[i]=FLIP[act][t[i]];t=tuple(t)
                key='welfare' if act=='G' else 'private'
                slacks.append((prof[s][key][i]-prof[t][key][i],i,act))
        if s==expected:
            assert all(g.lo>0 for g,_,_ in slacks)
            winner=[dict(actor='AB'[i]+'_'+act,slack=g.approx()) for g,i,act in slacks]
        else:
            attack=[(g,i,act) for g,i,act in slacks if g.hi<0]
            assert len(attack)>0,('unrefuted configuration',s)
            g,i,act=min(attack,key=lambda it:it[0].hi)
            losers.append(dict(config=''.join(s),actor='AB'[i]+'_'+act,slack=g.approx()))
    planner=('B','B');w=prof[planner]['welfare'][0]+prof[planner]['welfare'][1]
    gaps=[]
    for s in prof:
        if s==planner:continue
        v=w-prof[s]['welfare'][0]-prof[s]['welfare'][1]
        assert v.lo>0,('planner alternative not eliminated',s,v.approx())
        gaps.append(dict(config=''.join(s),planner_gap=v.approx()))
    return dict(NE=''.join(expected),planner='BB',selected_flips=winner,
      all_15_other_profiles_refuted=losers,planner_dominates_other_15=gaps,
      min_planner_gap_lower=min(t['planner_gap'][0] for t in gaps),
      P_P_fee_intervals=[v.approx() for v in prof[('P','P')]['fees']])

def main():
    data=[]
    scenarios=[(F(0),False,('B','P'),('0.380848479444378','0.359890064348462')),
      (F(1,2),False,('B','B'),('0.381077381587935','0.360088394043391')),
      (F(1,2),True,('B','P'),('0.380493572505174','0.359632374104520'))]
    for lam,pooled,expected,xx in scenarios:
        p=primitives(lam,pooled);cen=[F(Decimal(t)) for t in xx]
        bb=bounds(p);box,f=faces(p,cen);result=profile_certificate(p,box,expected)
        data.append(dict(lambda_=str(lam),pooled=pooled,price_contraction=bb,
          root_box=[list(map(str,bi)) for bi in box],rational_positivity_faces=f,first_stage=result))
        print('PASS',lam,'pooled',pooled,'P/P global uniqueness + exact rational root enclosure; all 16 profiles, 64 flips, planner')
    with open('stage11r_pp_global_certificate.json','w') as out:
        json.dump({'theory':'STAGE33-2026-10-28-THEORY-FREEZE-v1',
          'proof_scope':'Three baseline scenarios only. Exact rational inequalities + analytic BR contraction lemma in report; no claim of proof assistant verifying full lemma or arbitrary parameters.',
          'scenarios':data},out,indent=2)
if __name__=='__main__':main()
