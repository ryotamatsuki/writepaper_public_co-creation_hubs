"""Stage 4: endogenous-location interior diagnostic (no SPNE certification).

Quadratic Hotelling, public welfare and private profit.  Full coverage
in mixed duopoly; single private entrant; public program n in {1,2}.
The equations assume BOTH firms are active and all market users are
served.  Location results are solutions to interior FOCs only.
"""
from math import isclose

T,V,ALPHA,SIGMA,GAMMA,B,K=1.0,2.0,.09,1.4,1.0,.5,.001
FP,FG=.55,.4

def pilot(n):
    assert n in (1,2)
    r=ALPHA*n/(1+SIGMA*(n-1))
    e=r
    return dict(n=n,e=e,delta=GAMMA*n*e,
                H=r*e-.5*e*e+B*n*e-K*n)

def mixed_interior(m,n):
    p=pilot(n)
    delta=p["delta"]
    assert 0<=delta<T/4, "Interior-location candidate invalid."
    a=.25-delta/T
    b=.75-delta/T
    z=.5-2*delta/T
    price=T/2+2*delta
    profit=m*(T+4*delta)**2/(4*T)
    gross_social_entry=m*(T+4*delta)**2/(16*T)
    assert 0<a<z<b<1
    assert isclose(profit,4*gross_social_entry)
    enters=profit>=FP
    Wpubliconly=m*(V-T/12)+p["H"]-FG
    Wbuild=Wpubliconly+(gross_social_entry-FP if enters else 0)
    return {**p,"M":m,"a":a,"b":b,"cutoff":z,
            "price":price,"gross_private_profit":profit,
            "gross_social_entry":gross_social_entry,
            "private_enters":enters,"W_build":Wbuild}

def private_only(m):
    # For the declared V=2,T=1 parameters, the private monopoly
    # optimally locates at 1/2 and charges V-T/4 with full coverage.
    assert isclose(V,2*T)
    price=V-T/4
    entered=m*price>=FP
    W=m*(V-T/12)-FP if entered else 0.0
    return dict(enters=entered,price=price,W=W)

def policy(m):
    base=private_only(m)
    alternatives=[mixed_interior(m,n) for n in (1,2)]
    best=max(alternatives,key=lambda x:x["W_build"])
    choice=f"BUILD n={best['n']}" if best["W_build"]>=base["W"] else "DO NOT BUILD"
    return base,alternatives,choice

def efficient_welfare(a,b,delta,t=T):
    """Social welfare excluding common v: exact efficient allocation.

    This is a global *allocation* diagnostic, NOT a price Nash game.
    """
    def sq_integral(c,lo,hi):
        return ((hi-c)**3-(lo-c)**3)/3
    if abs(a-b)<1e-14:
        return delta-t*sq_integral(b,0,1) if delta>=0 else -t*sq_integral(a,0,1)
    z=(a+b)/2 + (delta/(2*t*(a-b)))
    z=max(0.0,min(1.0,z))
    if a<b:
        return delta*(1-z)-t*(sq_integral(a,0,z)+sq_integral(b,z,1))
    return delta*z-t*(sq_integral(b,0,z)+sq_integral(a,z,1))

def check_public_global_allocation(n,grid=1000):
    """Sample all possible public locations for fixed private b."""
    delta=pilot(n)["delta"]
    a=.25-delta/T
    b=.75-delta/T
    current=efficient_welfare(a,b,delta)
    vals=[(efficient_welfare(i/grid,b,delta),i/grid) for i in range(grid+1)]
    best=max(vals)
    return best[0]-current,best[1],a

if __name__=="__main__":
    print("Pilot / endogenous location candidates")
    for n in (1,2):
        p=mixed_interior(1.0,n)
        print("n",n,"e",p["e"],"delta",p["delta"],"H",p["H"],
              "a",p["a"],"b",p["b"],"cutoff",p["cutoff"],
              "private_profit",p["gross_private_profit"],
              "entry_social_gain",p["gross_social_entry"],
              "entry_mass_threshold",FP/p["gross_private_profit"])
        err,besta,current=check_public_global_allocation(n)
        print("public global efficient allocation grid gain",err,
              "grid argmax",besta,"FOC candidate",current)
        assert err<1e-6
    for m in (.10,.20,.50,.85,.90,1.0,1.18,1.20,1.50,2.10,3.0):
        base,opts,choice=policy(m)
        print("M",m,"policy",choice,"no_public",round(base["W"],6),
              "W1",round(opts[0]["W_build"],6),
              "W2",round(opts[1]["W_build"],6))
    assert policy(1.0)[2]=="BUILD n=1"
    assert policy(1.2)[2]=="DO NOT BUILD"
    assert policy(2.1)[2]=="BUILD n=2"
    print("All restricted-benchmark assertions passed.")
