#!/usr/bin/env python3
"""Independent Beta(2,2) portability attack on the frozen model's *technology*.

This is a separate numerical counter-model, not a replacement for the uniform
Stage8 baseline and not a complete proof of P/P equilibrium existence.
All other real and strategic primitives are held at their Stage9 baseline.
"""
import importlib.util
import json
from itertools import product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("freeze",ROOT/"verification/stage33_reproduce.py")
s=importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)

def mass(c):return 1-3*c*c+2*c**3
def skill(c):return .5-2*c**3+1.5*c**4
def moment(c):return mass(c),skill(c)

def clip(c):return max(0.,min(1.,c))

def best(fun,lo,hi):
    z,v=s.golden_max(fun,lo,hi)
    candidates=[(lo,fun(lo)),(z,v),(hi,fun(hi))]
    return max(candidates,key=lambda a:a[1])[0]

def private_fee(U,T):
    # Revenue f when f <= U (full coverage); none for f>=U+T.
    f0=U
    f1=best(lambda f:f*mass(clip((f-U)/T)),U,U+T)
    fee=max((f0,f1),key=lambda f:f*mass(clip((f-U)/T)))
    return fee,clip((fee-U)/T)

def premium_fee(i):
    delta=s.DELTA[i]
    f=best(lambda x:x*mass(clip(x/delta)),0.,delta)
    return f,moment(f/delta)

def pp_membership(fees,lam):
    cutoffs=[.5,.5]
    for iteration in range(2000):
        new=[]
        for i in (0,1):
            j=1-i
            a,b=s.matching_coefficients(i,lam,False)
            m,q=moment(cutoffs[j])
            t=s.DELTA[i]+b*q
            new.append(clip((fees[i]-s.direct_basic(i,"P")-a*m)/t))
        err=max(abs(new[i]-cutoffs[i]) for i in (0,1))
        cutoffs=[(new[i]+cutoffs[i])/2 for i in (0,1)]
        if err<1e-13:break
    else:raise AssertionError("Beta P/P membership continuation unresolved")
    return tuple(moment(t) for t in cutoffs)

def pp_equilibrium(lam):
    fees=[.65,.65]
    def revenue(i,own,other):
        ff=(own,other) if i==0 else (other,own)
        return s.N[i]*own*pp_membership(ff,lam)[i][0]
    for _ in range(160):
        replies=[]
        for i in (0,1):
            grid=[(z/20,revenue(i,z/20,fees[1-i])) for z in range(61)]
            z=max(range(61),key=lambda k:grid[k][1])
            left,right=max(0,(z-1)/20),min(3,(z+1)/20)
            replies.append(best(lambda x:revenue(i,x,fees[1-i]),left,right))
        err=max(abs(replies[i]-fees[i]) for i in (0,1))
        fees=[(fees[i]+replies[i])/2 for i in (0,1)]
        if err<1e-6:break
    else:raise AssertionError("Beta P/P fee reply did not converge")
    for i in (0,1):
        eq=revenue(i,fees[i],fees[1-i])
        bound=max(revenue(i,z/40,fees[1-i]) for z in range(121))
        if bound>eq+1e-6:raise AssertionError("Beta P/P profitable grid deviation")
    return tuple(fees),pp_membership(fees,lam)

def profile(state,lam):
    if state==("P","P"):
        fees,mq=pp_equilibrium(lam)
    else:
        fees=[0.,0.]
        mq=[None,None]
        for i in (0,1):
            if state[i]=="0":mq[i]=(0.,0.)
            if state[i] in ("G","B"):mq[i]=(1.,.5)
        for i in (0,1):
            if state[i]!="P":continue
            j=1-i
            a,b=s.matching_coefficients(i,lam,False)
            U=s.direct_basic(i,state[j])+a*mq[j][0]
            T=s.DELTA[i]+b*mq[j][1]
            fees[i],cut=private_fee(U,T)
            mq[i]=moment(cut)
        for i in (0,1):
            if state[i]=="B":fees[i],_=premium_fee(i)
    W=[];P=[]
    for i in (0,1):
        j=1-i
        basic=s.direct_basic(i,state[j])
        s0=s.PB*s.M[i]+s.OMEGA[i]*s.R0
        s1=lam*(s.M[i]+s.OMEGA[i]*s.PI)
        if state[i]=="0":z=0.
        elif state[i]=="G":z=s.N[i]*basic
        elif state[i]=="P":z=s.N[i]*(basic*mq[i][0]+s.DELTA[i]*mq[i][1])
        else:
            premium,up_mq=premium_fee(i)
            z=s.N[i]*(basic+s.DELTA[i]*up_mq[1])
        W.append(s.BASE_W[i]+z-s.FG[i]*(state[i] in ("G","B"))
                -s.FP[i]*(state[i] in ("P","B"))
                +s.C*mq[0][0]*mq[1][0]*s0
                +s.C*mq[0][1]*mq[1][1]*s1)
        prem_mass= mq[i][0] if state[i]=="P" else premium_fee(i)[1][0]
        P.append(0. if state[i] in ("0","G") else s.N[i]*fees[i]*prem_mass-s.FP[i])
    return {"welfare":tuple(W),"private":tuple(P),"fees":fees,"mq":mq}

def solve(lam):
    states={z:profile(z,lam) for z in product(s.STATES,repeat=2)}
    best=max(sum(x["welfare"]) for x in states.values())
    planner=[z for z in states if abs(sum(states[z]["welfare"])-best)<1e-8]
    NE=[]
    for z in states:
        if all(states[z]["welfare" if a=="G" else "private"][i]+1e-8>=
               states[tuple(s.flip(z[k],a) if k==i else z[k] for k in (0,1))]["welfare" if a=="G" else "private"][i]
               for i in (0,1) for a in ("G","P")):
            NE.append(z)
    return {"lambda":lam,"NE":NE,"planner":planner,
       "region_B_addition":states[("B","B")]["welfare"][1]-states[("B","P")]["welfare"][1],
       "PP_fees":states[("P","P")]["fees"]}

if __name__=="__main__":
    assert abs(mass(.5)-.5)<1e-15
    assert abs(skill(0)-.5)<1e-15
    out=[solve(0.),solve(.5)]
    print(json.dumps({"distribution":"Beta(2,2)","results":out,
        "scope":"Selected candidate numerical P/P continuations only, not a new frozen theorem"},indent=2))
    if out[0]["NE"]==[("B","P")] and out[1]["NE"]==[("B","B")]:
        raise AssertionError("Beta alternative accidentally reproduces target ownership switch")
    print("PASS: adverse Beta(2,2) ownership portability test at frozen costs")
