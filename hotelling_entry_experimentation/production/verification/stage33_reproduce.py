#!/usr/bin/env python3
"""Stage 9 independent baseline reproduction of the frozen Stage 33 example.

Standard-library only. All four independently chosen binary entry actions are
evaluated in all 16 profiles. The P/P price continuation is COMPUTATIONAL:
best-reply iteration, bounded global grid search, and convergence diagnostics;
this is NOT a proof of global Nash existence or uniqueness for all primitives.
Theory authority: STAGE33-2026-10-28-THEORY-FREEZE-v1.
"""
import csv
import json
import math
import sys
from pathlib import Path
from itertools import product

REGIONS = ("A", "B")
STATES = ("0", "G", "P", "B")  # none, public, private, both
N = (3.0, 1.0)
INCUMBENTS = (5, 2)
V = (0.15, 0.10)
W = (0.22, 0.45)
DELTA = (1.0, 1.0)
FG = (0.09, 0.0405)
FP = (0.65, 0.24)
OMEGA = (0.8, 0.2)
RHO, MU, P0, ZETA, K, CHI = 0.4, 0.25, 0.05, 0.005, 0.0004, 0.9
C = RHO * MU * N[0] * N[1]

def market(n):
    """Frozen linear-quadratic two-good Cournot market, A=D=b=d=1."""
    den = 2 * (n + 1) - n * CHI**2
    x, y, x0 = (2 - CHI)/den, ((n + 1) - n * CHI)/den, 1/(n + 1)
    cs1 = ((n*x)**2 + 2*CHI*n*x*y + y*y)/2
    cs0 = (n*x0)**2/2
    m = (cs1 + n*x*x) - (cs0 + n*x0*x0)
    return m, y*y, cs0 + n*x0*x0

MARKETS = [market(n) for n in INCUMBENTS]
M = tuple(t[0] for t in MARKETS)
PI = sum(t[1] for t in MARKETS)
BASE_W = tuple(t[2] for t in MARKETS)
EFFORT = ZETA * PI/K
PB = P0 + ZETA * EFFORT
R0 = P0*PI + ZETA**2*PI**2/(2*K)

def clip(x):
    return min(1.0, max(0.0, x))

def direct_basic(i, other_state):
    return V[i] + RHO*W[i]*(other_state != "0")

def matching_coefficients(i, strength, pooled):
    j = 1-i
    a = RHO*MU*N[j]*OMEGA[i]*R0
    b = RHO*MU*N[j]*OMEGA[i]*strength*PI
    if pooled:
        a += b/4
        b = 0.0
    return a,b

def pp_membership(fees, strength, pooled):
    """Demand fixed point at simultaneous posted prices; damped iteration."""
    cutoffs = [0.5,0.5]
    for _ in range(1000):
        updated = []
        for i in range(2):
            j = 1-i
            a,b = matching_coefficients(i,strength,pooled)
            mj = 1-cutoffs[j]
            qj = (1-cutoffs[j]**2)/2
            t = DELTA[i]+b*qj
            updated.append(clip((fees[i]-direct_basic(i,"P")-a*mj)/t))
        err = max(abs(updated[i]-cutoffs[i]) for i in range(2))
        cutoffs = [(updated[i]+cutoffs[i])/2 for i in range(2)]
        if err < 1e-14:
            break
    else:
        raise AssertionError("P/P participation fixed point did not converge")
    return tuple((1-t,(1-t*t)/2) for t in cutoffs)

def golden_max(fun,left,right):
    ratio = (math.sqrt(5)-1)/2
    c,d = right-ratio*(right-left),left+ratio*(right-left)
    fc,fd = fun(c),fun(d)
    for _ in range(80):
        if fc < fd:
            left,c,fc = c,d,fd
            d = left+ratio*(right-left)
            fd = fun(d)
        else:
            right,d,fd = d,c,fc
            c = right-ratio*(right-left)
            fc = fun(c)
    z = (left+right)/2
    return z,fun(z)

def pp_price_equilibrium(strength, pooled):
    """Numerical baseline continuation; not an unrestricted theorem."""
    prices = [0.65,0.65]
    def revenue(i,own,other):
        ff = (own,other) if i == 0 else (other,own)
        return N[i]*own*pp_membership(ff,strength,pooled)[i][0]
    for iteration in range(300):
        proposed = []
        for i in range(2):
            # Scan all prices in [0,3]; refine around the best grid cell.
            grid = [(3*z/60,revenue(i,3*z/60,prices[1-i]))
                    for z in range(61)]
            imax = max(range(61),key=lambda z:grid[z][1])
            left,right = max(0,3*(imax-1)/60),min(3,3*(imax+1)/60)
            own,_ = golden_max(lambda f:revenue(i,f,prices[1-i]),left,right)
            proposed.append(own)
        err = max(abs(proposed[i]-prices[i]) for i in range(2))
        prices = [(prices[i]+proposed[i])/2 for i in range(2)]
        if err < 1e-12:
            break
    else:
        raise AssertionError("P/P price best replies did not converge")
    # Additional numeric global-deviation GRID diagnostic on [0,3].
    for i in range(2):
        at_eq = revenue(i,prices[i],prices[1-i])
        maxgrid = max(revenue(i,3*z/240,prices[1-i]) for z in range(241))
        assert at_eq + 1e-7 >= maxgrid
    return tuple(prices),pp_membership(prices,strength,pooled)

def profile(state, strength, pooled=False):
    assert len(state)==2 and all(x in STATES for x in state)
    if state == ("P","P"):
        fees, mq = pp_price_equilibrium(strength,pooled)
    else:
        fees = [0.0,0.0]
        mq = [None,None]
        for i in range(2):
            if state[i] == "0":
                mq[i]=(0.0,0.0)
            elif state[i] in ("G","B"):
                mq[i]=(1.0,0.5)
        for i in range(2):
            if state[i] != "P":
                continue
            j = 1-i
            a,b = matching_coefficients(i,strength,pooled)
            U = direct_basic(i,state[j])+a*mq[j][0]
            T = DELTA[i]+b*mq[j][1]
            if U<T:
                fees[i]=(U+T)/2
                cutoff=(T-U)/(2*T)
            else:
                fees[i]=U
                cutoff=0.0
            mq[i]=(1-cutoff,(1-cutoff**2)/2)
        for i in range(2):
            if state[i]=="B":
                fees[i]=DELTA[i]/2
    welfare,private = [],[]
    for i in range(2):
        j = 1-i
        s0 = PB*M[i]+OMEGA[i]*R0
        s1 = strength*(M[i]+OMEGA[i]*PI)
        if pooled:
            s0 += s1/4
            s1 = 0.0
        if state[i]=="0":
            z = 0.0
        elif state[i]=="G":
            z = N[i]*direct_basic(i,state[j])
        elif state[i]=="P":
            z = N[i]*(direct_basic(i,state[j])*mq[i][0]+DELTA[i]*mq[i][1])
        else:
            z = N[i]*(direct_basic(i,state[j])+3*DELTA[i]/8)
        wi = (BASE_W[i]+z-FG[i]*(state[i] in ("G","B"))
              -FP[i]*(state[i] in ("P","B"))
              +C*mq[0][0]*mq[1][0]*s0+C*mq[0][1]*mq[1][1]*s1)
        priv = (0.0 if state[i] in ("0","G") else
                N[i]*fees[i]*(mq[i][0] if state[i]=="P" else 0.5)-FP[i])
        welfare.append(wi)
        private.append(priv)
    return dict(fees=tuple(fees),mq=tuple(mq),
                welfare=tuple(welfare),private=tuple(private))

def flip(current,actor):
    if actor=="G":
        return {"0":"G","G":"0","P":"B","B":"P"}[current]
    return {"0":"P","P":"0","G":"B","B":"G"}[current]

def solve(strength,pooled=False):
    assert PB+strength <= 1+1e-12
    profiles = {s:profile(s,strength,pooled)
                for s in product(STATES,repeat=2)}
    equilibria = []
    for s,here in profiles.items():
        stable = True
        for i in range(2):
            for actor in ("G","P"):
                alt = list(s)
                alt[i] = flip(s[i],actor)
                candidate = profiles[tuple(alt)]
                key = "welfare" if actor=="G" else "private"
                if candidate[key][i] > here[key][i]+1e-9:
                    stable = False
        if stable:equilibria.append(s)
    max_welfare = max(sum(x["welfare"]) for x in profiles.values())
    planner = [s for s,x in profiles.items()
               if abs(sum(x["welfare"])-max_welfare)<1e-9]
    return profiles,equilibria,planner,max_welfare

def local_b_entry(strength,pooled=False):
    d=profile(("B","B"),strength,pooled)["welfare"][1]
    p=profile(("B","P"),strength,pooled)["welfare"][1]
    return d-p

def check():
    scenarios=((0.0,False,("B","P")),
               (0.5,False,("B","B")),
               (0.5,True,("B","P")))
    for strength,pooled,expected in scenarios:
        profiles,equilibria,planner,_=solve(strength,pooled)
        assert len(profiles)==16
        assert equilibria==[expected],(strength,pooled,equilibria)
        assert planner==[("B","B")],(strength,pooled,planner)
    heterogeneous=solve(0.5,False)[3]
    homogeneous=solve(0.5,True)[3]
    assert abs(heterogeneous-2.407438368575148)<1e-11
    assert abs(heterogeneous-homogeneous)<1e-11
    p=profile(("B","P"),0.5,False)
    assert abs(p["fees"][1]-0.6410074280605)<1e-10
    assert abs(p["mq"][1][0]-0.639945174957)<1e-10
    assert abs(local_b_entry(0.5,False)-0.00003167010245)<1e-10
    assert abs(local_b_entry(0.5,True)+0.00045677076684)<1e-10
    lo,hi=0.0,0.5
    assert local_b_entry(lo)<0<local_b_entry(hi)
    for _ in range(60):
        mid=(lo+hi)/2
        if local_b_entry(mid)>0:hi=mid
        else:lo=mid
    threshold=(lo+hi)/2
    assert abs(threshold-0.3511840318401)<1e-9,threshold
    return threshold

def main():
    threshold=check()
    rows=[]
    for strength,pooled in ((0.0,False),(0.5,False),(0.5,True)):
        all_profiles,equilibria,planner,sw=solve(strength,pooled)
        scenario=("pooled" if pooled else "skill")
        for s,record in all_profiles.items():
            rows.append({
                "lambda":str(strength),"scenario":scenario,
                "entry_A":s[0],"entry_B":s[1],
                "fee_A":format(record["fees"][0],".12f"),
                "fee_B":format(record["fees"][1],".12f"),
                "m_A":format(record["mq"][0][0],".12f"),
                "m_B":format(record["mq"][1][0],".12f"),
                "q_A":format(record["mq"][0][1],".12f"),
                "q_B":format(record["mq"][1][1],".12f"),
                "W_A":format(record["welfare"][0],".12f"),
                "W_B":format(record["welfare"][1],".12f"),
                "private_A":format(record["private"][0],".12f"),
                "private_B":format(record["private"][1],".12f"),
                "is_NE":int(s in equilibria),
                "is_planner":int(s in planner)})
        print(json.dumps({"lambda":strength,"scenario":scenario,
            "NE":equilibria,"planner":planner,
            "restricted_planner_value":round(sw,12)},sort_keys=True))
    print("local-public-entry cutoff lambda =",format(threshold,".13f"))
    if "--csv" in sys.argv:
        index = sys.argv.index("--csv")
        output=Path(sys.argv[index+1])
        output.parent.mkdir(parents=True,exist_ok=True)
        with output.open("w",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print("wrote",output,"rows",len(rows))
    print("PASS: baseline regression, all 16 profiles/scenario, fee, welfare and cutoff")

if __name__ == "__main__":
    main()
