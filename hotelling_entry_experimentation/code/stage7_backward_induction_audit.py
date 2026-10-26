"""Independent backward-induction audit of Stage 6 (Stage 7 checkpoint).

Exact rational arithmetic, separate realized entry paths, welfare decomposition,
and broad (finite-grid) unilateral location checks. Python standard library only.
NOT a proof of unique SPNE or novelty. See STAGE_6_CONSTRUCTIVE_SPNE.md
for the analytic restricted existence proof and off-path price selection.
"""
from fractions import Fraction as F
from math import sqrt

t, v, alpha, sigma, gamma = F(1), F(4), F(9,100), F(7,5), F(1)
B, k, Fp, Fg = F(1,2), F(1,1000), F(11,20), F(2,5)
MIN_D = (sqrt(5)-2)/4


def pilot(n):
    r=alpha*n/(1+sigma*(n-1))
    delta=gamma*n*r
    H=r*r/2+B*n*r-k*n
    assert MIN_D <= float(delta/t) < .25
    return {"n":n,"e":r,"pilot_private_net":r*r/2,"delta":delta,"H":H}


def mixed(n, mass):
    p=pilot(n)
    delta=p["delta"]
    d=delta/t
    a=F(1,4)-d; b=F(3,4)-d; z=F(1,2)-2*d
    price=t/2+2*delta
    private_gross=mass*price*(1-z)
    gov_fee_revenue=mass*price*z
    consumer_surplus=mass*(v-t/F(48)+delta/2+delta*delta/t-price)
    user_welfare=consumer_surplus+private_gross+gov_fee_revenue
    expected_user_welfare=mass*(v-t/F(48)+delta/2+delta*delta/t)
    assert user_welfare==expected_user_welfare
    gain=mass*(t+4*delta)**2/(16*t)
    assert private_gross==4*gain
    return {**p,"a":a,"b":b,"z":z,"price":price,
            "gov_fee_revenue":gov_fee_revenue,"private_gross":private_gross,
            "private_net":private_gross-Fp,"consumer_surplus":consumer_surplus,
            "social_use":user_welfare,"entry_social_gain":gain}


def private_no_public(mass):
    gross=mass*(v-t/4)
    entered=gross>=Fp
    return {
        "entry":entered,
        "location":F(1,2) if entered else None,
        "price":v-t/4 if entered else None,
        "private_gross":gross if entered else F(0),
        "private_net":gross-Fp if entered else F(0),
        "welfare":mass*(v-t/12)-Fp if entered else F(0),
    }


def public_n(mass,n):
    benchmark=mixed(n,mass)
    entry=benchmark["private_gross"]>=Fp
    W=mass*(v-t/12)+benchmark["H"]-Fg
    if entry:
        W+=benchmark["entry_social_gain"]-Fp
    return {"entry":entry,
            "public_location":benchmark["a"] if entry else F(1,2),
            "private_location":benchmark["b"] if entry else None,
            "public_price":benchmark["price"] if entry else F(0),
            "private_price":benchmark["price"] if entry else None,
            "private_net":benchmark["private_net"] if entry else F(0),
            "welfare":W,"mixed_hypothetical":benchmark}


def policy(mass):
    m=F(str(mass))
    options={"no_public":private_no_public(m)["welfare"],
             "public_n1":public_n(m,1)["welfare"],
             "public_n2":public_n(m,2)["welfare"]}
    return max(options,key=options.get),options


def selected_price(a,b,delta):
    # Selected off-path Nash price equilibrium from Stage 6.
    a,b,delta=float(a),float(b),float(delta)
    D=abs(a-b)
    if D<1e-13:return 0.0,delta,1.0
    if b>a:
        z=(a+b)/2-delta/(2*D)
        if z>0:
            price=2*D*(1-z)
            return price,price,1-z
        gap=delta-D*(a+b)
    else:
        z=(a+b)/2+delta/(2*D)
        if z<1:
            price=2*D*z
            return price,price,z
        gap=delta-D*(2-a-b)
    public_price=max(0.,2*D-gap)
    return public_price,public_price+gap,1.0


def efficient_use(a,b,delta):
    def int_square(location,left,right):
        return ((right-location)**3-(left-location)**3)/3
    if abs(a-b)<1e-13:
        return 4+delta-int_square(a,0,1)
    D=abs(b-a)
    if b>a:
        cutoff=max(0.,min(1.,(a+b)/2-delta/(2*D)))
        return 4-(int_square(a,0,cutoff)+int_square(b,cutoff,1))+delta*(1-cutoff)
    cutoff=max(0.,min(1.,(a+b)/2+delta/(2*D)))
    return 4-(int_square(b,0,cutoff)+int_square(a,cutoff,1))+delta*cutoff


def global_location_grid_check(n,steps=10000):
    base=mixed(n,F(1))
    a,b,d=map(float,(base["a"],base["b"],base["delta"]))
    reference_profit=float(base["private_gross"])
    reference_welfare=efficient_use(a,b,d)
    for i in range(steps+1):
        alt=i/steps
        q,p,s=selected_price(a,alt,d)
        assert p*s<=reference_profit+2e-10,(n,"private",alt,p*s,reference_profit)
        w=efficient_use(alt,b,d)
        assert w<=reference_welfare+2e-10,(n,"government",alt,w,reference_welfare)


if __name__=="__main__":
    for n in (1,2):
        p=pilot(n)
        print("pilot",n,{k:float(val) if isinstance(val,F) else val
                         for k,val in p.items()})
        global_location_grid_check(n)
        print("location best-response grid passed",n)
        m=mixed(n,F(1))
        print("mixed prices/shares/profits",float(m["price"]),float(m["z"]),
              float(m["private_gross"]),float(m["entry_social_gain"]))
    assert policy(F(1))[0]=="public_n1"
    assert not public_n(F(1),1)["entry"]
    assert public_n(F(1),1)["public_location"]==F(1,2)
    assert public_n(F(1),2)["entry"]
    assert policy(F(6,5))[0]=="no_public"
    assert policy(F(21,10))[0]=="public_n2"
    for M in ("0.05","0.10","0.20","0.50","0.90","1.00","1.20","1.50","2.10"):
        chosen, vals=policy(M)
        print(f"M={M} chosen={chosen} W:",
              {label:round(float(value),8) for label,value in vals.items()})
    print("entry cutoffs",float(Fp/(v-t/4)),
          float(Fp/mixed(1,F(1))["private_gross"]),
          float(Fp/mixed(2,F(1))["private_gross"]))
    print("PASS: exact rational accounting and independent location-grid audit")
