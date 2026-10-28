#!/usr/bin/env python3
"""Exact decomposition of existing frozen P/B comparisons; no new primitives.

This checks algebra and the original conditional comparisons, NOT priority,
all-lambda SPNE, or new Lean coverage. Standard library only.
"""
from fractions import Fraction as F
from pathlib import Path
import sys, json, hashlib, subprocess

HERE=Path(__file__).resolve().parent
ROOT=next((p for p in HERE.parents if (p/'formal/Stage33Formal.lean').exists()),HERE)
REVIEW=ROOT/'hotelling_entry_experimentation/production/reviews'
sys.path.insert(0,str(REVIEW if REVIEW.exists() else HERE))
import stage11r_pp_global_certificate as b
OUT=Path(__file__).resolve().parent
C=F(3,10)

def val(x): return {'exact':str(x),'approx':float(x)}

def comparison(lam,pooled):
    p=b.primitives(lam,pooled);i=1;d=F(1)
    u=b.basic(i,'B');U=u+p['a'][i];T=d+p['b'][i]/2
    ell=max(F(0),(T-U)/(2*T))
    S0=p['pb']*p['m'][i]+p['omega'][i]*p['r0']
    K=p['m'][i]+p['omega'][i]*p['pi']
    direct=p['n'][i]*(u*ell-d/8+d*ell*ell/2)
    baseline=C*S0*ell
    technological=C*lam*K*(ell/4 if pooled else ell*ell/4)
    threshold=direct+baseline+technological
    bp=b.profile(('B','P'),p,[None,None]);bb=b.profile(('B','B'),p,[None,None])
    exact_diff=bb['welfare'][i]-bp['welfare'][i]
    assert exact_diff.lo==exact_diff.hi==threshold-p['fg'][i]
    return dict(lam=lam,pooled=pooled,p=p,u=u,U=U,T=T,ell=ell,
                direct=direct,baseline=baseline,technological=technological,
                threshold=threshold,K=K,S0=S0)

def main():
    zero=comparison(F(0),False)
    skill=comparison(F(1,2),False)
    pool=comparison(F(1,2),True)
    lH,lP=skill['ell'],pool['ell'];delta=lH-lP
    p=skill['p'];u=skill['u'];U=skill['U'];h=p['b'][1]/2
    assert delta==h*(1+h+2*U)/(4*(1+h))>0
    access=delta*(p['n'][1]*(u+(lH+lP)/2)+C*skill['S0'])
    incidence=C*F(1,2)*skill['K']*(lH*lH-lP)/4
    assert access+incidence==skill['threshold']-pool['threshold']
    # The baseline's incremental successful venture is locally harmful.
    assert skill['K']<0 and skill['S0']<0
    assert skill['direct']>skill['threshold']>0
    # The three-way cost band is small because the lambda change alters the
    # private cutoff only slightly and its service gain is partly offset.
    dz=skill['direct']-zero['direct']
    dm=skill['baseline']-zero['baseline']
    dt=skill['technological']
    assert dz+dm+dt==skill['threshold']-zero['threshold']
    rows=[]
    for name,r in [('skill_zero',zero),('skill_half',skill),('pooled_half',pool)]:
        rows.append({'scenario':name,**{k:val(r[k]) for k in
            ('ell','direct','baseline','technological','threshold','K','S0')}})
    band=skill['threshold']-zero['threshold'];cost=p['fg'][1]
    result={
      'scope':'Exact frozen-model diagnostic only; no all-lambda equilibrium or novelty certification',
      'base_commit':'d2520c9afbf9b2b742a5f04fbba6159eb489db61',
      'scenarios':rows,
      'skill_vs_pooled':{'cutoff_difference':val(delta),'access_and_baseline':val(access),
                        'technological_incidence':val(incidence),
                        'threshold_difference':val(access+incidence)},
      'skill_half_minus_zero':{'cutoff_change':val(lH-zero['ell']),
            'direct_service':val(dz),'baseline_matching':val(dm),
            'skill_success_incidence':val(dt),'cost_band_width':val(band),
            'band_percent_of_frozen_cost':float(100*band/cost),
            'upper_slack_percent_of_frozen_cost':float(100*(skill['threshold']-cost)/cost)},
      'unchanged_formal_source_sha256':hashlib.sha256((ROOT/'formal/Stage33Formal.lean').read_bytes()).hexdigest() if (ROOT/'formal/Stage33Formal.lean').exists() else None,
    }
    (OUT/'diagnostic_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS exact original thresholds, cutoff comparison, reoptimised decomposition, cost width')
    print(json.dumps(result['skill_half_minus_zero'],indent=2))

if __name__=='__main__':main()
