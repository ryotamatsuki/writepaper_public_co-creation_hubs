"""Stage 6: constructive SPNE certificate for a restricted mixed Hotelling game.

ASSUMPTIONS:
- Public hub maximizes total surplus; private hub maximizes profit.
- Quadratic taste mismatch, mass M of uniform users, no variable costs.
- Nonnegative prices and outside option zero.
- Baseline value v >= 4*t (chosen here: v=4, t=1).
- Spillover delta arises from optimal pilot effort e=r_n, with
  r_n=alpha*n/(1+sigma*(n-1)), delta=gamma*n*r_n.
- d=delta/t in [(sqrt(5)-2)/4,1/4), so the particular price
  CONTINUATION RULE supports the interior location equilibrium.
- Exactly one prospective private hub. Binary program scope {1,2}.
- This is a constructive equilibrium existence example; NOT uniqueness,
  an equilibrium refinement, or a novelty-certified theorem.
"""
from math import sqrt, isclose

T=1.0
V=4.0
ALPHA=0.09
SIGMA=1.4
GAMMA=1.0
B=0.5
K=0.001
FP=0.55
FG=0.4
D_MIN=(sqrt(5.0)-2.0)/4.0
assert V>=4*T


def pilot(n):
    r=ALPHA*n/(1+SIGMA*(n-1))
    e=r
    return e,GAMMA*n*e,e*e/2+B*n*e-K*n


def selected_prices(a,b,delta,t=T):
    """A Nash price continuation (pG,pP,P share) at EVERY a,b in [0,1].

    Nonnegative prices. Full market coverage follows when v>=4t,
    given delta/t<1/4. In all-private regimes, SELECT the smallest
    public price among equilibria. Other price equilibria exist.
    """
    assert 0<=a<=1 and 0<=b<=1 and 0<delta<t/4
    D=abs(a-b)
    if D<1e-13:
        return 0.0,delta,1.0,"same-location"
    if b>a:
        z=(a+b)/2-delta/(2*t*D)
        if z>0:
            p=2*t*D*(1-z)
            return p,p,1-z,"interior"
        advantage=delta-t*D*(a+b)
    else:
        z=(a+b)/2+delta/(2*t*D)
        if z<1:
            p=2*t*D*z
            return p,p,z,"interior"
        advantage=delta-t*D*(2-a-b)
    assert advantage>=-1e-12
    q=max(0.0,2*t*D-advantage)
    return q,q+advantage,1.0,"all-private"


def private_share(a,b,delta,pG,pP,t=T):
    if abs(b-a)<1e-13:
        return float(delta+pG>=pP-1e-12)
    if b>a:
        z=(a+b)/2+(pP-pG-delta)/(2*t*(b-a))
        return min(1.0,max(0.0,1-z))
    z=(a+b)/2+(delta+pG-pP)/(2*t*(a-b))
    return min(1.0,max(0.0,z))


def efficient_gross_welfare(a,b,delta,v=V,t=T):
    """Exact optimal assignment gross user surplus; no fixed costs."""
    def intsq(c,lo,hi):
        return ((hi-c)**3-(lo-c)**3)/3
    if abs(a-b)<1e-13:
        return v+delta-t*intsq(b,0,1)
    if b>a:
        z=max(0.0,min(1.0,(a+b)/2-delta/(2*t*(b-a))))
        return (v-t*(intsq(a,0,z)+intsq(b,z,1))
                +delta*(1-z))
    z=max(0.0,min(1.0,(a+b)/2+delta/(2*t*(a-b))))
    return v-t*(intsq(b,0,z)+intsq(a,z,1))+delta*z


def location_solution(delta,t=T):
    d=delta/t
    assert D_MIN<=d<.25
    return .25-d,.75-d


def location_checks(delta,grid=10000):
    a,b=location_solution(delta)
    q,p,share,branch=selected_prices(a,b,delta)
    expected_profit=(T+4*delta)**2/(4*T)
    assert branch=="interior" and isclose(p*share,expected_profit,rel_tol=1e-10)

    for i in range(grid+1):
        other=i/grid
        q1,p1,s1,_=selected_prices(a,other,delta)
        assert p1*s1<=expected_profit+2e-10, (delta,"P deviation",other)
        welfare=efficient_gross_welfare(other,b,delta)
        w_star=efficient_gross_welfare(a,b,delta)
        assert welfare<=w_star+2e-10,(delta,"G deviation",other)
    return dict(delta=delta,a=a,b=b,p=p,private_share=share,
                profit=expected_profit,entry_welfare=(T+4*delta)**2/(16*T))


def verify_price_best_responses(delta,step=20):
    # Grid audit supplements, but does not replace, the price-game proof.
    for i in range(step+1):
        a=i/step
        for j in range(step+1):
            b=j/step
            q,p,share,regime=selected_prices(a,b,delta)
            # Public alone at q can serve every user; outside is irrelevant
            # for ANY private price deviation.
            assert min(V-q-T*a*a,V-q-T*(1-a)**2)>=0
            incumbent_profit=p*share
            for h in range(251):
                pp=h/50
                alt=pp*private_share(a,b,delta,q,pp)
                assert alt<=incumbent_profit+1e-8,(a,b,q,p,pp,alt)
    return True


def policy(mass):
    # Without public hub, private monopoly chooses characteristic 1/2.
    # v>=4t => full-cover monopoly price v-t/4.
    no_entry_private_gross=mass*(V-T/4)
    no_g_w=(mass*(V-T/12)-FP
            if no_entry_private_gross>=FP else 0.0)
    policies=[]
    for n in (1,2):
        e,delta,H=pilot(n)
        a,b=location_solution(delta)
        private_gross=mass*(T+4*delta)**2/(4*T)
        gain=mass*(T+4*delta)**2/(16*T)
        entering=private_gross>=FP  # enter on indifference
        w=mass*(V-T/12)+H-FG
        if entering:
            w+=gain-FP
        policies.append(dict(n=n,delta=delta,e=e,H=H,
                             entry=entering,a=a,b=b,W=w))
    choice=max(policies,key=lambda x:x["W"])
    return dict(M=mass,choice=("BUILD n="+str(choice["n"])
                               if choice["W"]>=no_g_w else "NO BUILD"),
                no_public_welfare=no_g_w,alternatives=policies)


if __name__=="__main__":
    for n in (1,2):
        e,d,H=pilot(n)
        print("LOCATION",location_checks(d,grid=10000))
        print("PRICE",n,verify_price_best_responses(d,step=20))
    for m in (.05,.08,.10,.14,.20,.50,.85,.90,1.0,1.18,1.20,1.50,2.10):
        out=policy(m)
        print("POLICY",m,out["choice"],
              "no-public",round(out["no_public_welfare"],6),
              "by n",[(x["n"],round(x["W"],6),x["entry"])
                      for x in out["alternatives"]])
    assert policy(0.05)["choice"]=="NO BUILD"
    assert policy(.20)["choice"]=="BUILD n=2"
    assert policy(1.0)["choice"]=="BUILD n=1"
    assert policy(1.2)["choice"]=="NO BUILD"
    assert policy(2.10)["choice"]=="BUILD n=2"
    print("PASS: tested constructive continuation and stage-0 policy")
