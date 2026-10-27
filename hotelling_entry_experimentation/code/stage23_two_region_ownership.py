"""Stage 23: Exact rational simultaneous-entry enumeration for two asymmetric regions.

Three stage-1 decision makers:
  G chooses site vector in {(0,0),(1,0),(0,1),(1,1)} jointly;
  P_U, P_R independently choose binary entry.
16 profiles checked, INCLUDING 3- and 4-site off-path ones.

Downstream sites are anchored to urban/rural regions, participants fully
enroll LOCALLY if a hub exists, and split equally between G/P when both
operate locally. NO ENDOGENOUS TWO-SIDED FEE CHOICE: the private broker
collects a FIXED contractual percentage of ex-post partner bargaining proceeds.
Innovators invest before match formation, effort best response is exact.
This is an SPNE of the restricted game, NOT the unrestricted Hotelling model.

Run with Python stdlib only. All payoff/Nash operations use fractions.Fraction.
"""
from fractions import Fraction as F
from itertools import product
from csv import DictWriter

ACT=tuple(product((0,1), repeat=2))
HIST=tuple(product(ACT,ACT))
beta=F(2,3); mu=F(2,5); sigma=F(7,10)
v=F(1,20); t=F(1,50); A=c=F(1)
kappa=(beta-beta**2/2)*A*A*mu*mu/c
phi=sigma*beta*(1-beta)*A*A*mu*mu/c
fp=(F(18,1000),F(12,1000))
fg=(F(110,1000),F(37,1000))
alpha=F(1,2)

def name(g,p):
    return ''.join('B' if g[i] and p[i] else
                   'G' if g[i] else 'P' if p[i] else '-' for i in range(2))

def payoff(g,p,r,rho,fg=fg,fp=fp,alpha=alpha):
    m=(F(1),r); n=(g[0]+p[0],g[1]+p[1])
    hubs=[[],[]]; weighted=[F(0),F(0)]
    for i in range(2):
        if not n[i]: continue
        for o,on in (('G',g[i]),('P',p[i])):
            if not on:continue
            gamma=F(1) if o=='G' else alpha
            hub=dict(o=o,I=m[i]/n[i],R=m[i]/n[i],gamma=gamma)
            hubs[i].append(hub)
            weighted[i]+=hub['R']*gamma
    gross=F(0)
    for i in range(2):
        if n[i]:gross+=2*m[i]*(v-t/(12*n[i]**2))
        for hub in hubs[i]:
            z=hub['R']+rho*hub['gamma']*weighted[1-i]
            q=mu*z
            assert 0<=q<=1
            hub['e']=beta*A*q/c
            gross+=kappa*hub['I']*z*z
    priv={}
    for i in range(2):
        for hub in hubs[i]:
            if hub['o']!='P':continue
            cross=sum(other['gamma']*other['I']*other['e']
                      for other in hubs[1-i])
            priv[i]=sigma*(1-beta)*A*mu*hub['R']*(
                hub['I']*hub['e']+rho*hub['gamma']*cross)
    social=gross-sum(fg[i]*g[i]+fp[i]*p[i] for i in range(2))
    return social,priv,gross

def solve(r,rho,fg=fg,fp=fp,alpha=alpha):
    data={(g,p):payoff(g,p,r,rho,fg,fp,alpha) for g,p in HIST}
    ne=[]
    for g,p in HIST:
        w,pi,_=data[g,p]
        if any(data[gd,p][0]>w for gd in ACT):continue
        profitable_deviation=False
        for i in range(2):
            alt=tuple(1-p[j] if i==j else p[j] for j in range(2))
            currently=pi.get(i,F(0))-fp[i] if p[i] else F(0)
            newly=data[g,alt][1].get(i,F(0))-fp[i] if alt[i] else F(0)
            if newly>currently:
                profitable_deviation=True
                break
        if not profitable_deviation:ne.append((g,p))
    best=max(data[s][0] for s in HIST)
    planner=[s for s in HIST if data[s][0]==best]
    return ne,planner,data

def root_crossings(r,fg=fg,fp=fp,alpha=alpha):
    """Exact pairwise two-site PP/PG and PG/GG threshold coefficients
    followed by numerical quadratic roots; global optimum separately audited."""
    K=float(kappa*r*(1+r)); h=float(alpha)
    from math import sqrt
    def positive_root(quadratic,linear,cost):
        if cost<=0:return 0.0
        return 2*cost/(linear+sqrt(linear**2+4*quadratic*cost))
    first=positive_root(K*h*h*(1-h*h),2*K*h*(1-h),
                        float(fg[1]-fp[1]))
    second=positive_root(K*(1-h*h),2*K*(1-h),
                         float(fg[0]-fp[0]))
    return first,second

def run():
    assert kappa==F(16,225) and phi==F(28,1125)
    r=F(7,10)
    cutoffs=root_crossings(r)
    assert .49789<cutoffs[0]<.49791
    assert .70956<cutoffs[1]<.70958
    print("rho cutoffs",cutoffs)
    for rho,expected in ((F(2,5),'PP'),(F(3,5),'PG'),(F(4,5),'GG')):
        ne,planner,data=solve(r,rho)
        assert len(planner)==1 and name(*planner[0])==expected
        print("r",r,"rho",rho,"NE",[name(*s) for s in ne],
              "planner",[name(*s) for s in planner],
              "W*",float(data[planner[0]][0]))
    ne,planner,data=solve(r,F(3,5))
    assert set(name(*s) for s in ne)=={"PP","PG","GP","GG"}
    for s in ne:
        print(name(*s),"W",float(data[s][0]),
              "Private gross",{i:float(z) for i,z in data[s][1].items()})
    for rho in (F(0),F(1,4),F(1,2),F(3,4),F(1)):
        ne,planner,_=solve(r,rho)
        print("rho",float(rho),"pure NE",','.join(name(*s) for s in ne),
              "social optimum",','.join(name(*s) for s in planner))
    with open("stage23_two_region_equilibria.csv","w",newline="",encoding="utf-8-sig") as f:
        fields=("rural_ratio","cross_efficiency","pure_NE_all",
                "social_optima_all","n_pure_NE","best_social_welfare",
                "W_PP","W_PG","W_GP","W_GG")
        w=DictWriter(f,fieldnames=fields);w.writeheader()
        for ir in range(1,21):
            rr=F(ir,20)
            for j in range(21):
                rho=F(j,20)
                ne,planner,data=solve(rr,rho)
                rec={"rural_ratio":float(rr),"cross_efficiency":float(rho),
                     "pure_NE_all":"|".join(name(*x) for x in ne),
                     "social_optima_all":"|".join(name(*x) for x in planner),
                     "n_pure_NE":len(ne),
                     "best_social_welfare":float(data[planner[0]][0])}
                for x in HIST:
                    label=name(*x)
                    if label in ("PP","PG","GP","GG"):
                        rec["W_"+label]=float(data[x][0])
                w.writerow(rec)
    print("All 16 profiles and 420 exact grid points checked.")

if __name__=="__main__":
    run()
