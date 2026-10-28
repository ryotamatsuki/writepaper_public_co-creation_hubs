#!/usr/bin/env python3
"""P1b: consistency of empirical implications with frozen model primitives.

Checks identities and strict inequalities using fractions.Fraction. These are
*conditional theoretical moments*, NOT a regression, empirical validation,
causal identification, proof of generic entry or full Lean formalization.
Reuses the independent mathematical continuation module from the prior P0a
audited branch and the same three fixed Stage8 environments.
"""
from fractions import Fraction as F
from itertools import product
import json
import stage11r_pp_global_certificate as b

SCENARIOS=((F(0),False,"skill_lambda_zero"),
           (F(1,2),False,"skill_lambda_half"),
           (F(1,2),True,"pooled_lambda_half"))
DATA=[]

def verify_screening(p):
    """Fixed rival moments and fixed private-only U,T on interior tariff arm."""
    U=F(2,5);T=F(5,4);N=F(3)
    assert 0<U<T
    l1,l2=F(1,5),F(3,5)
    outcome=[]
    for l in (l1,l2):
        fee=U+T*l
        mass=1-l
        moment=(1-l*l)/2
        mean=moment/mass
        assert mean==(1+l)/2
        outcome.append({"f":float(fee),"m":float(mass),"q":float(moment),"mean":float(mean)})
    assert outcome[1]['f']>outcome[0]['f']
    assert outcome[1]['m']<outcome[0]['m']
    assert outcome[1]['q']<outcome[0]['q']
    assert outcome[1]['mean']>outcome[0]['mean']
    assert -1/T<0 and -l1/T<0 and 1/(2*T)>0
    return {"comparison":"fixed U,T and rival moments; two admissible fee/cutoff points",
            "analytical_derivatives":{"dm_df":str(-1/T),"dq_df_at_l1":str(-l1/T),
                                     "dmean_df":str(1/(2*T))},
            "points":outcome}

def verify_coprovision(p):
    """Within opponent's full participation: P-only versus public+premium B."""
    N=p['n'][1]
    U=p['u'][1]+p['a'][1]  # other region A: m_A=1
    T=1+p['b'][1]/2       # other region A: q_A=1/2, Delta_B=1
    assert 0<U<T
    cutoff=(T-U)/(2*T)
    paid_participation=1-cutoff
    RP=N*(U+T)*(U+T)/(4*T)
    RB=N/4
    assert 0<cutoff<1 and RP>RB and paid_participation>F(1,2)
    assert RB-p['fp'][1]>0  # B provider can survive with public basic service
    assert p['fp'][1]==F(6,25)
    return {"opposite_region_full_coverage":True,"P_fee":float((U+T)/2),
            "P_only_take_up_share":float(paid_participation),
            "P_gross_revenue":float(RP),
            "B_premium_fee":0.5,"B_premium_share":0.5,
            "B_gross_revenue":float(RB),"B_net_private_profit":float(RB-p['fp'][1]),
            "conditional_revenue_difference_positive":True}

def verify_foreign_welfare(p):
    """Identity for local B addition at fixed region A status B."""
    a=('B','P');z=('B','B')
    pa=b.profile(a,p,[None,None])  # P/P fee box unused for these statuses
    pb=b.profile(z,p,[None,None])
    ell=1-pa['mq'][1][0]      # exact Fraction-wrapped interval
    assert ell.lo==ell.hi and 0<ell.lo<1
    wA=pb['welfare'][0]-pa['welfare'][0]
    S0=p['pb']*p['m'][0]+p['omega'][0]*p['r0']
    S1=p['lam']*(p['m'][0]+p['omega'][0]*p['pi'])
    if p['pooled']:
        S0+=S1/4
        S1=F(0)
    C=F(2,5)*F(1,4)*p['n'][0]*p['n'][1]
    RHS=C*(S0*ell.lo+S1*ell.lo*ell.lo/4)
    assert wA.lo==wA.hi==RHS,("foreign-welfare identity failed",wA.approx(),float(RHS))
    # No unconditional sign claim. Observables must identify consumers,
    # incumbents, venture profit division, and counterfactual market outcomes.
    return {"conditional_foreign_welfare_exact":str(RHS),
            "conditional_foreign_welfare_approx":float(RHS),
            "source":"W_A(B,B)-W_A(B,P) equals foreign incremental matching welfare",
            "pooled":p['pooled']}

for lam,pooled,name in SCENARIOS:
    p=b.primitives(lam,pooled)
    DATA.append({"scenario":name,
                 "screening":verify_screening(p),
                 "coprovision":verify_coprovision(p),
                 "regional_wedge":verify_foreign_welfare(p)})
    print("PASS",name,"conditional fee screening, co-provision revenue and regional welfare accounting")

with open("stage11r_p1b_empirical_implications_checks.json","w") as file:
    json.dump({"version":"STAGE33-2026-10-28-THEORY-FREEZE-v1",
               "status":"algebraic model-implied identities and exact baseline consistency; NOT observed data or causal estimates",
               "scenarios":DATA},file,indent=2)
print("PASS all 3 model-to-observable mappings across all 3 frozen baseline environments")
