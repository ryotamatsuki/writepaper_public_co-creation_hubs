#!/usr/bin/env python3
"""Minimal copula extension; all frozen primitives imported without mutation.

Moment identities/tests use Fraction. Equilibrium enumeration is floating-point
numerics, not a proof of global price Nash existence or uniqueness for kappa<1.
Run from any directory: python exploration_kappa_20261028/numerics.py
"""
from pathlib import Path
import sys, json, math
import numpy as np
from fractions import Fraction as F
from itertools import product
from scipy.optimize import minimize_scalar, root

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent))
sys.path.insert(0,str(ROOT/'hotelling_entry_experimentation/production/reviews'))
import stage11r_pp_global_certificate as frozen

from moments import moments

def basic(i,state): return float(frozen.basic(i,state))

def primitives(lam,pooled):
    p=frozen.primitives(F(str(lam)),pooled)
    return {k:([float(x) for x in v] if isinstance(v,(tuple,list)) else float(v) if isinstance(v,F) else v) for k,v in p.items()}

def pp_demand(fees,p,kappa,start=0.5):
    mq=[(start,start/2,start/2)]*2
    for it in range(200):
        new=[moments(p['u'][i]+p['a'][i]*mq[1-i][0],p['b'][i]*mq[1-i][1],1.0,fees[i],kappa) for i in (0,1)]
        if max(abs(new[i][j]-mq[i][j]) for i in (0,1) for j in (0,1))<2e-14:return new
        mq=new
    raise AssertionError('participation iteration failed')

def global_max(fun,upper,extra=(),grid=80):
    """Whole compact feasible fee interval, local-peak refinement + all kinks.

    An exhaustive grid is a diagnostic, NOT a mathematical global certificate.
    Revenue is zero above upper by the bound U+h+D <= u+a+b/2+1.
    """
    xx=sorted(set([upper*i/grid for i in range(grid+1)]+[float(x) for x in extra if 0<=x<=upper]))
    yy=[fun(x) for x in xx]; cand=list(zip(xx,yy))
    for j in range(1,len(xx)-1):
        if yy[j]>=yy[j-1] and yy[j]>=yy[j+1]:
            sol=minimize_scalar(lambda f:-fun(f),bounds=(xx[j-1],xx[j+1]),method='bounded',options={'xatol':2e-13})
            cand.append((sol.x,-sol.fun))
    return max(cand,key=lambda t:t[1])

def pp_prices(p,kappa,start=(0.65,0.65)):
    def revenue(i,f,other):
        ff=[other,other];ff[i]=f
        return p['n'][i]*f*pp_demand(ff,p,kappa)[i][0]
    def br(i,other,grid=80):
        upper=p['u'][i]+p['a'][i]+1+p['b'][i]/2
        # Fee kinks depend on participation; scan and refine all local peaks.
        return global_max(lambda f:revenue(i,f,other),upper,grid=grid)
    prices=list(start)
    for it in range(80):
        new=[br(i,prices[1-i])[0] for i in (0,1)]
        if max(abs(prices[i]-new[i]) for i in (0,1))<8e-8:
            prices=new;break
        prices=new
    else:raise AssertionError('price BR iteration failed')
    # Sharpen the numerical candidate using the exact central-strip demand
    # derivatives, then recheck every root against the whole fee interval.
    def stationary(ff):
        qq=pp_demand(ff,p,kappa);jac=np.eye(4);direct=np.zeros((4,2))
        for i in (0,1):
            j=1-i;U=p['u'][i]+p['a'][i]*qq[j][0];h=p['b'][i]*qq[j][1];x=ff[i]-U
            assert h<x<1,(i,x,h)
            mf=kappa/(1+h)+(1-kappa)
            mh=kappa*x/(1+h)**2+(1-kappa)/2
            qf=kappa*x/(1+h)**2+(1-kappa)/2
            qh=kappa*x*x/(1+h)**3+(1-kappa)/3
            jac[2*i,2*j]=-mf*p['a'][i];jac[2*i,2*j+1]=-mh*p['b'][i]
            jac[2*i+1,2*j]=-qf*p['a'][i];jac[2*i+1,2*j+1]=-qh*p['b'][i]
            direct[2*i,i]=-mf;direct[2*i+1,i]=-qf
        deriv=np.linalg.solve(jac,direct)
        return [qq[i][0]+ff[i]*deriv[2*i,i] for i in (0,1)]
    sharp=root(stationary,prices,tol=1e-11)
    assert sharp.success or max(map(abs,stationary(sharp.x)))<1e-12
    prices=list(sharp.x)
    mq=pp_demand(prices,p,kappa)
    regret=[br(i,prices[1-i],grid=400)[1]-revenue(i,prices[i],prices[1-i]) for i in (0,1)]
    demand_seeds=max(abs(pp_demand(prices,p,kappa,0)[i][j]-pp_demand(prices,p,kappa,1)[i][j]) for i in (0,1) for j in (0,1))
    assert max(regret)<2e-10,(prices,regret)
    return prices,mq,{'numerical_max_price_regret':max(regret),'participation_seed_discrepancy':demand_seeds,'iterations':it+1}

def profile(s,p,kappa,pp_cache=None):
    fee=[0.,0.];mq=[None,None];diag={}
    if s==('P','P'):
        fee,mq,diag=pp_cache or pp_prices(p,kappa)
    else:
        for i in (0,1):
            if s[i]=='0':mq[i]=(0.,0.,0.)
            if s[i] in ('G','B'):mq[i]=(1.,.5,.5)
        for i in (0,1):
            j=1-i
            if s[i]=='P':
                U=basic(i,s[j])+p['a'][i]*mq[j][0];h=p['b'][i]*mq[j][1]
                fee[i],_=global_max(lambda f:p['n'][i]*f*moments(U,h,1.,f,kappa)[0],U+h+1,(U,U+h,U+1,U+h+1))
                A=kappa/(1+h)+(1-kappa);c=(1-kappa)*h/2
                candidate=(1+c+A*U)/(2*A)
                if h<candidate-U<1:
                    old_rev=fee[i]*moments(U,h,1.,fee[i],kappa)[0]
                    assert candidate*moments(U,h,1.,candidate,kappa)[0]>=old_rev-1e-12
                    fee[i]=candidate
                mq[i]=moments(U,h,1.,fee[i],kappa)
            if s[i]=='B':fee[i]=.5
    wealth=[];profit=[]
    for i in (0,1):
        u=basic(i,s[1-i]); m,q,r=mq[i]
        s0=p['pb']*p['m'][i]+p['omega'][i]*p['r0']
        s1=p['lam']*(p['m'][i]+p['omega'][i]*p['pi'])
        if p['pooled']:s0+=s1/4;s1=0
        z={'0':0.,'G':p['n'][i]*u,'P':p['n'][i]*(u*m+r),'B':p['n'][i]*(u+3/8)}[s[i]]
        w=p['base'][i]+z-p['fg'][i]*(s[i] in ('G','B'))-p['fp'][i]*(s[i] in ('P','B'))
        w+=.3*mq[0][0]*mq[1][0]*s0+.3*mq[0][1]*mq[1][1]*s1
        v=p['n'][i]*fee[i]*(m if s[i]=='P' else .5)-p['fp'][i] if s[i] in ('P','B') else 0
        wealth.append(w);profit.append(v)
    return dict(welfare=wealth,private=profit,fees=fee,mqr=mq,diagnostics=diag)

def all_profiles(p,kappa):
    pp=pp_prices(p,kappa)
    prof={s:profile(s,p,kappa,pp) for s in product(frozen.STATES,repeat=2)}
    ne=[];margins={}
    for s in prof:
        slacks=[]
        for i in (0,1):
            for act in ('G','P'):
                t=list(s);t[i]=frozen.FLIP[act][s[i]];t=tuple(t)
                key='welfare' if act=='G' else 'private'
                slacks.append(prof[s][key][i]-prof[t][key][i])
        margins[''.join(s)]=slacks
        if min(slacks)>=-1e-10:ne.append(''.join(s))
    planner=max(prof,key=lambda s:sum(prof[s]['welfare']))
    threshold=prof[('B','B')]['welfare'][1]-prof[('B','P')]['welfare'][1]+p['fg'][1]
    return {'kappa':kappa,'NE':ne,'planner':''.join(planner),'public_B_threshold_given_B':threshold,'profiles':{''.join(s):v for s,v in prof.items()},'four_entry_slacks':margins}

def exact_tests():
    # Empty, full, diagonal reduction, square triangles, and independence theta-selection.
    U,h,D=F(1,5),F(2,5),F(1)
    for k in [F(0),F(1,3),F(1)]:
        assert moments(U,h,D,U,k)==(1,F(1,2),F(1,2))
        assert moments(U,h,D,U+h+D,k)==(0,0,0)
    assert moments(F(0),F(1),F(1),F(1),F(0))==(F(1,2),F(1,3),F(1,3))
    assert moments(F(0),F(0),F(1),F(1,2),F(0))==(F(1,2),F(1,4),F(3,8))
    for f in [F(0),F(1,5),F(7,10),F(6,5),F(2)]:
        m,q,r=moments(U,h,D,f,F(1));assert q==r==m-m*m/2
    # Corr(theta,z)=k exactly, uniform margins.
    for k in [F(0),F(1,3),F(1)]:assert (k/F(3)+(1-k)/4-F(1,4))*12==k
    # At k=0, h>0, selection in theta remains strict.
    m,q,r=moments(U,h,D,F(4,5),F(0));assert q>m/2 and r>m/2
    # At h=0, fees and m,r invariant in k, but q changes.
    assert moments(U,F(0),D,F(4,5),F(0))[0::2]==moments(U,F(0),D,F(4,5),F(1))[0::2]
    for k in [F(0),F(1,3),F(1)]:assert moments(F(0),F(0),D,F(1,2),k)[2]==F(3,8)
    # Right derivative at h=0 of the reoptimised conditional central-strip fee.
    for k in [F(0),F(1)]:
        eps=F(1,10**8);c=(D-U)/(2*D)
        A=k/(D+eps)+(1-k)/D;cc=(1-k)*eps/(2*D)
        f=(1+cc+A*U)/(2*A)
        old_f=(D+U)/2;old_m=(D+U)/(2*D)
        new_m=moments(U,eps,D,f,k)[0]
        assert abs((f-old_f)/eps-(1+k)/4)<F(1,10**7)
        assert abs((new_m-old_m)/eps-(F(1,4)+k*(c-F(3,4)))/D)<F(1,10**7)
    return 'PASS exact Fraction moments, corners, kappa endpoints, corr, residual skill selection'

def seed_checks():
    rows=[]
    for lam,pooled in [(0,False),(.5,False),(.5,True)]:
        for kappa in [0,.5,1]:
            p=primitives(lam,pooled)
            solutions=[pp_prices(p,kappa,start) for start in [(0,0),(.65,.65),(1.3,1.3)]]
            spread=max(abs(x[0][i]-solutions[0][0][i]) for x in solutions for i in (0,1))
            assert spread<1e-11
            rows.append({'lambda':lam,'pooled':pooled,'kappa':kappa,'three_price_seed_spread':spread})
    (Path(__file__).parent/'price_seed_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('PASS PP three-seed diagnostics',max(x['three_price_seed_spread'] for x in rows))

def participation_certificate():
    rows=[]
    for lam,pooled in [(F(0),False),(F(1,2),False),(F(1,2),True)]:
        pp=frozen.primitives(lam,pooled)
        bounds=[pp['a'][i]+pp['b'][i] for i in (0,1)]
        assert max(bounds)<1
        rows.append({'lambda':str(lam),'pooled':pooled,'Delta_A':'1','Delta_B':'1',
          'regional_lipschitz_bounds_exact':list(map(str,bounds)),
          'regional_lipschitz_bounds_float':list(map(float,bounds)),
          'joint_contraction_upper_bound':float(max(bounds))})
    out={'scope':'all fees, all kappa in [0,1], participation fixed point only; no uniqueness of prices asserted',
      'proof':'Each m,q is separately 1/Delta-Lipschitz in U and h under diagonal and independent components, preserved by mixture and clipping. Therefore joint sup-norm map has Lipschitz <= max_i(a_i+b_i)/Delta_i < 1. Banach gives a unique participation fixed point for every posted fee pair.',
      'scenarios':rows}
    (Path(__file__).parent/'participation_contraction.json').write_text(json.dumps(out,indent=2)+'\n')

def baseline_regression(report,lam,pooled):
    xx={ (0,False):('0.380848479444378','0.359890064348462'),(.5,False):('0.381077381587935','0.360088394043391'),(.5,True):('0.380493572505174','0.359632374104520')}[(lam,pooled)]
    p=frozen.primitives(F(str(lam)),pooled);box,_=frozen.faces(p,list(map(F,xx)))
    errors={'fees':0.,'m':0.,'q':0.,'welfare':0.,'private':0.}
    for s in product(frozen.STATES,repeat=2):
        old=frozen.profile(s,p,box);new=report['profiles'][''.join(s)]
        for key in ('fees','welfare','private'):
            for i in (0,1):errors[key]=max(errors[key],abs(new[key][i]-float((old[key][i].lo+old[key][i].hi)/2)))
        for i in (0,1):
            for j,key in enumerate(('m','q')):errors[key]=max(errors[key],abs(new['mqr'][i][j]-float((old['mq'][i][j].lo+old['mq'][i][j].hi)/2)))
    assert max(errors.values())<1e-7,errors
    return errors

def main():
    participation_certificate()
    output={'scope':'Exploratory numerical equilibria for kappa<1; exact frozen certificates retained at kappa=1; no global existence/uniqueness theorem asserted.','exact_tests':exact_tests(),'scenarios':[]}
    for lam,pooled,name in [(0,False,'skill_zero'),(.5,False,'skill_half'),(.5,True,'pooled_half')]:
        p=primitives(lam,pooled)
        rows=[all_profiles(p,k/10) for k in range(11)]
        errs=baseline_regression(rows[-1],lam,pooled)
        output['scenarios'].append({'name':name,'lambda':lam,'pooled':pooled,'kappa1_frozen_max_abs_errors':errs,'rows':rows})
        print(name,[(r['kappa'],r['NE'],r['public_B_threshold_given_B']) for r in rows],flush=True)
    (Path(__file__).parent/'results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(output['exact_tests'])
if __name__=='__main__':
    if '--seed-checks' in sys.argv:seed_checks()
    else:main()
