#!/usr/bin/env python3
"""Stage11R P1: all 16 entry configurations over open parameter rectangles.

Exact rational endpoint arithmetic; no scanning/grid interpolation.
Imports only the already-certificated frozen baseline economic primitives and
P/P global price lemma from stage11r_pp_global_certificate.py. The novel
independent check here is the entire interval of every four-player deviation,
the restricted planner, binding cost thresholds, and an explicit non-cost
omega_B neighborhood with omega_A=1-omega_B.

The 'interval' class is an inclusion arithmetic, not a Lean proof of economics.
Any failed/ambiguous strict inequality raises AssertionError (fail closed).
"""
import json
from itertools import product
from fractions import Fraction as F
from decimal import Decimal
import stage11r_pp_global_certificate as c
I=c.I
BASE_COST_B=F(81,2000)
COST_INNER=(F('0.0404252'),F('0.0405315'))
SCENARIOS=(
 (F(0),False,('B','P'),('0.380848479444378','0.359890064348462')),
 (F(1,2),False,('B','B'),('0.381077381587935','0.360088394043391')),
 (F(1,2),True,('B','P'),('0.380493572505174','0.359632374104520')),
)
ACTORS=((0,'G'),(0,'P'),(1,'G'),(1,'P'))
STATES=tuple(product(c.STATES,repeat=2))

def payoff_slack(profile,s,i,actor):
    t=list(s);t[i]=c.FLIP[actor][t[i]];t=tuple(t)
    k='welfare' if actor=='G' else 'private'
    return profile[s][k][i]-profile[t][k][i]

def certify_all(p,box,expected,tag):
    prof={s:c.profile(s,p,box) for s in STATES}
    stable_slack=[]
    losers=[]
    for s in STATES:
        margins=[(payoff_slack(prof,s,i,act),i,act) for i,act in ACTORS]
        if s==expected:
            assert all(v.lo>0 for v,_,_ in margins),(tag,'NE not strictly secured',s,[v.approx() for v,_,_ in margins])
            stable_slack=[{'actor':'AB'[i]+'_'+act,'lower':float(v.lo),'upper':float(v.hi)} for v,i,act in margins]
        else:
            witnesses=[(v,i,act) for v,i,act in margins if v.hi<0]
            assert witnesses,(tag,'alternative may be NE',s,[(v.approx(),i,act) for v,i,act in margins])
            v,i,act=min(witnesses,key=lambda t:t[0].hi)
            losers.append({'state':''.join(s),'actor':'AB'[i]+'_'+act,'profitable_lower':float(-v.hi)})
    target=('B','B');sw=prof[target]['welfare'][0]+prof[target]['welfare'][1]
    losers_planner=[]
    for s in STATES:
        if s==target:continue
        gap=sw-prof[s]['welfare'][0]-prof[s]['welfare'][1]
        assert gap.lo>0,(tag,'planner uncertain',s,gap.approx())
        losers_planner.append({'state':''.join(s),'gap_lower':float(gap.lo)})
    return {'NE':''.join(expected),'planner':'BB',
      'weakest_selected_flip_lower':min(x['lower'] for x in stable_slack),
      'weakest_planner_gap_lower':min(x['gap_lower'] for x in losers_planner),
      'entry_flip_witnesses':losers,'selected_flips':stable_slack,
      'planner_15_witnesses':losers_planner}

def pp_box_and_sufficient_conditions(p,center,radius):
    # Every inequality below is exact rational inclusion arithmetic.
    a,b,u=p['a'],p['b'],p['u']
    hh=[I.wrap(a[i]+b[i]) for i in (0,1)]
    z=hh[0]*hh[1]
    assert max(h.hi for h in hh)<1 and z.hi<1
    bounds=[]
    for i in (0,1):
        lo=(1-z-u[i]-a[i])/(2+b[i]/2-z)
        hi=(1+b[i]/2-u[i])/(2+b[i]/2)
        assert lo.lo>0 and hi.hi<1 and lo.hi<hi.lo
        fee_lo=u[i]+lo
        fee_hi=u[i]+a[i]+(1+b[i]/2)*hi
        assert (u[i]+a[i]).hi<fee_lo.lo and fee_hi.hi<(u[i]+1).lo
        h=(b[i]+b[i]*hh[1-i]+b[1-i]*hh[i]+hh[i])
        B=hh[i]+(1+b[i]/2)*h/(2*(1-z))
        assert B.hi<F(1,100),(i,B.approx())
        bounds.append({'region':'AB'[i],'cutoff_global_BR_lower':lo.approx(),
           'cutoff_global_BR_upper':hi.approx(),'fee_BR_lower':fee_lo.approx(),
           'fee_BR_upper':fee_hi.approx(),'BR_lipschitz_upper':float(B.hi)})
    box,faces=c.faces(p,center,radius)
    return box,faces,bounds

def omega_cost_box(p,omega_eps,fg_a_eps,fg_b_eps):
    # Linked omega_A=1-omega_B; interval over-approximation may lose
    # correlations, never assert spurious dependence/independence.
    wB=I(F(1,5)-omega_eps,F(1,5)+omega_eps)
    wA=1-wB
    p['omega']=(wA,wB)
    aa=[F(1,10)*p['n'][1-i]*p['omega'][i]*p['r0'] for i in (0,1)]
    bb=[F(1,10)*p['n'][1-i]*p['omega'][i]*p['lam']*p['pi'] for i in (0,1)]
    if p['pooled']:
        aa=[aa[i]+bb[i]/4 for i in (0,1)]
        bb=[I(0),I(0)]
    p['a'],p['b']=aa,bb
    p['fg']=(I(F(9,100)-fg_a_eps,F(9,100)+fg_a_eps),
             I(BASE_COST_B-fg_b_eps,BASE_COST_B+fg_b_eps))
    return p

def public_B_threshold(p,box):
    prof={s:c.profile(s,p,box) for s in (('B','P'),('B','B'))}
    # W_B(BB)-W_B(BP) = threshold-F_GB (exact affine cost dependence).
    gain=prof[('B','B')]['welfare'][1]-prof[('B','P')]['welfare'][1]
    return I(BASE_COST_B)+gain

def main():
    baseline=[];robust=[];cost_only=[]
    thresholds={}
    for lam,pooled,expected,xx in SCENARIOS:
        name=('pooled_' if pooled else 'skill_')+('half' if lam else 'zero')
        cen=[F(Decimal(x)) for x in xx]
        pb=c.primitives(lam,pooled)
        base_box,faces,br=pp_box_and_sufficient_conditions(pb,cen,F(1,10**9))
        thresh=public_B_threshold(pb,base_box)
        thresholds[name]=thresh
        baseline.append({'scenario':name,'threshold':thresh.approx(),
                         'face_signs':[t['sign'] for t in faces],
                         'price_BR_bounds':br})
        print('P/P threshold',name,thresh.approx(),flush=True)
        # Uniform inner interval: combines all three NE outcomes, not only
        # the two half-lambda technology comparisons.
        ip=c.primitives(lam,pooled)
        ip['fg']=(F(9,100),I(*COST_INNER))
        one=certify_all(ip,base_box,expected,name+'_cost_only')
        cost_only.append({'scenario':name,**one})
        print('PASS cost-only all sixteen',name,flush=True)
        # Explicit noncost parameter neighborhood of omega_B, preserving
        # omega_A+omega_B=1, together with both government cost perturbations.
        rp=c.primitives(lam,pooled)
        rp=omega_cost_box(rp,F(1,10**6),F(1,10**6),F(1,10**6))
        rbox,rfaces,rbr=pp_box_and_sufficient_conditions(rp,cen,F(1,10**7))
        diag=certify_all(rp,rbox,expected,name+'_3d')
        robust.append({'scenario':name,'omega_B':[0.2-1e-6,0.2+1e-6],
             'F_GA':[0.09-1e-6,0.09+1e-6],
             'F_GB':[float(BASE_COST_B-F(1,10**6)),float(BASE_COST_B+F(1,10**6))],
             'cutoff_box_radius':1e-7,'faces':[x['sign'] for x in rfaces],
             'fee_BR_bounds':rbr,**diag})
        print('PASS omega/cost open-rectangle',name,flush=True)
    # Distinguish 2-way hetero vs pooled at lambda=.5 from all 3-way
    # configurations (which additionally impose the skill lambda=0 BP choice).
    t0=thresholds['skill_zero'];th=thresholds['skill_half'];tp=thresholds['pooled_half']
    assert t0.lo>tp.hi and t0.hi<th.lo,(t0.approx(),tp.approx(),th.approx())
    lo,hi=COST_INNER
    assert lo>t0.hi and hi<th.lo
    report={'theory':'STAGE33-2026-10-28-THEORY-FREEZE-v1',
       'proof':'All strict parameter-box checks use endpoint rational arithmetic plus analytic global price best-response lemma. No grids; independent human/Lean full-game proof not asserted.',
       'conditional_lambda_half_skill_pooled_band':[tp.approx(),th.approx()],
       'three_scenario_binding_public_B_cost_thresholds':{'lower_B_G_skill_lambda_zero':t0.approx(),
             'upper_B_G_skill_lambda_half':th.approx(),
             'pooled_B_G_lambda_half_not_binding':tp.approx()},
       'three_scenario_certified_cost_only_interior':[float(lo),float(hi)],
       'cost_width_3way_inner':float(hi-lo),
       'all_three_cost_only':cost_only,
       'explicit_open_parameter_neighborhoods':robust,
       'base_thresholds_and_pp_sufficient_bounds':baseline}
    with open('stage11r_p1_open_region_certificate.json','w') as fp:json.dump(report,fp,indent=2)
    print('PASS ALL: three-way cost band; all 16 x 4 entry deviations; restricted planner; 3D noncost neighborhood')
if __name__=='__main__':main()
