#!/usr/bin/env python3
"""Post-freeze referee sensitivity study, original Stage9 primitives unchanged.
Rebuilds all 16 strategic profiles for a small preregistered factor grid.
Costs shift ONLY first-stage government welfare, not downstream solutions.
The P/P candidate remains numerical and is NOT a global continuation theorem.
"""
import csv
import importlib.util
import json
import math
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"verification/stage33_reproduce.py"
spec=importlib.util.spec_from_file_location("stage33_frozen",SOURCE)
s=importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)
OUT=ROOT/"output"
OUT.mkdir(exist_ok=True)
BASE={"N":s.N,"DELTA":s.DELTA,"CHI":s.CHI,
      "INCUMBENTS":s.INCUMBENTS,"FG":s.FG}
STATES=s.STATES
assert abs(s.local_b_entry(.5,False)-.00003167010245)<1e-10

def reset(n_b=1.0,delta_b=1.0,chi=.9,n_inc_b=2):
    s.N=(3.0,float(n_b))
    s.DELTA=(1.0,float(delta_b))
    s.CHI=float(chi)
    s.INCUMBENTS=(5,int(n_inc_b))
    s.C=s.RHO*s.MU*s.N[0]*s.N[1]
    mkt=[s.market(n) for n in s.INCUMBENTS]
    s.M=tuple(x[0] for x in mkt)
    s.PI=sum(x[1] for x in mkt)
    s.BASE_W=tuple(x[2] for x in mkt)
    s.EFFORT=s.ZETA*s.PI/s.K
    s.PB=s.P0+s.ZETA*s.EFFORT
    s.R0=s.P0*s.PI+s.ZETA**2*s.PI**2/(2*s.K)
    assert 0<s.PB and s.PB+.5<=1

def classify(profiles,government_B_cost):
    """Reclassify 16 profiles from underlying payoffs, without rerunning prices."""
    def pay(state,index,actor):
        v=profiles[state]["welfare" if actor=="G" else "private"][index]
        if actor=="G" and index==1 and state[1] in ("G","B"):
            v-=government_B_cost-BASE["FG"][1]
        return v
    ne=[]
    best=-float("inf")
    pl=[]
    for state in profiles:
        if all(pay(state,i,actor)+1e-9>=
               pay(tuple(s.flip(state[k],actor) if k==i else state[k] for k in (0,1)),i,actor)
               for i in (0,1) for actor in ("G","P")):
            ne.append("".join(state))
        val=sum(pay(state,i,"G") for i in (0,1))
        if val>best+1e-9:
            best,pl=val,["".join(state)]
        elif abs(val-best)<=1e-9:pl.append("".join(state))
    return "+".join(sorted(ne)) or "none","+".join(sorted(pl))

def wage(state,profiles,i):
    return profiles[state]["welfare"][i]

def scan(label,lambda_,delta_b=1.,n_b=1.,chi=.9,n_inc_b=2,pooled=False):
    reset(n_b,delta_b,chi,n_inc_b)
    profiles,_,_,_=s.solve(lambda_,pooled)
    local=wage(("B","B"),profiles,1)-wage(("B","P"),profiles,1)
    foreign=wage(("B","B"),profiles,0)-wage(("B","P"),profiles,0)
    cost_star=BASE["FG"][1]+local
    records=[]
    for cost_int in range(3980,4191,10):  # 0.03980 to 0.04190 steps 0.00010
        cost=cost_int/100000.0
        ne,plan=classify(profiles,cost)
        records.append({"family":label,"lambda":lambda_,"Delta_B":delta_b,
            "N_B":n_b,"chi":chi,"n_B":n_inc_b,"pooled":int(pooled),
            "F_GB":f"{cost:.5f}","NE":ne,"planner":plan,
            "local_critical_cost":f"{cost_star:.12f}",
            "foreign_gain":f"{foreign:.12f}",
            "national_critical_cost":f"{cost_star+foreign:.12f}"})
    return records

def draw_phase_svg(rows,filename):
    """Black-and-white category/contrast phase plot, no visual inflation."""
    xs=sorted({float(r["F_GB"]) for r in rows})
    ys=sorted({float(r["lambda"]) for r in rows})
    w,h=780,395; left=72;right=22;top=35;bottom=60
    dx=(w-left-right)/len(xs);dy=(h-top-bottom)/len(ys)
    xidx={x:i for i,x in enumerate(xs)}
    yidx={y:i for i,y in enumerate(ys)}
    palette={"BB": "#d0d0d0", "BP":"#6f6f6f", "none":"#ffffff"}
    arr=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="white"/>',
         f'<text x="{w//2}" y="20" text-anchor="middle" font-size="15">Equilibrium map: government fixed cost vs skill complementarity</text>']
    for r in rows:
        x=float(r["F_GB"]);y=float(r["lambda"]);z=r["NE"]
        px=left+xidx[x]*dx;py=top+(len(ys)-yidx[y]-1)*dy
        color=palette.get(z,"#f4f4f4")
        arr.append(f'<rect x="{px:.1f}" y="{py:.1f}" width="{dx:.1f}" height="{dy:.1f}" fill="{color}" stroke="white" stroke-width=".4"/>')
        if z not in palette:
            arr.append(f'<text x="{px+dx/2:.1f}" y="{py+dy/2:.1f}" font-size="8" text-anchor="middle">{z}</text>')
    arr.append(f'<rect x="{left}" y="{top}" width="{w-left-right}" height="{h-top-bottom}" fill="none" stroke="black"/>')
    arr.append(f'<text x="{w//2}" y="{h-12}" text-anchor="middle" font-size="13">F_GB (original baseline = 0.04050)</text>')
    arr.append(f'<text x="18" y="{h//2}" transform="rotate(-90 18 {h//2})" text-anchor="middle" font-size="13">lambda</text>')
    for v in (xs[0],.0405,xs[-1]):
        ix=min(range(len(xs)),key=lambda i:abs(xs[i]-v))
        arr.append(f'<text x="{left+(ix+.5)*dx:.1f}" y="{h-bottom+16}" font-size="11" text-anchor="middle">{xs[ix]:.5f}</text>')
    for v in ys:
        arr.append(f'<text x="{left-9}" y="{top+(len(ys)-yidx[v]-.5)*dy+4:.1f}" font-size="11" text-anchor="end">{v:.2f}</text>')
    arr.append(f'<text x="{left+300}" y="{h-34}" font-size="11">Dark: BP; light: BB; others labelled</text>')
    arr.append("</svg>")
    filename.write_text("\n".join(arr),encoding="utf-8")

def main():
    rows=[]
    # Every scan fully recomputes sixteen profiles for each predeclared parameter value.
    for v in (0.,.10,.20,.30,.40,.50):
        rows+=scan("lambda_cost",v)
    for v in (.8,.9,1.,1.1,1.2):
        rows+=scan("delta_cost",.5,delta_b=v)
    for v in (.8,1.,1.2):
        rows+=scan("population_cost",.5,n_b=v)
    for v in (.85,.9,.95):
        rows+=scan("substitution_cost",.5,chi=v)
    # The pooled comparison at lambda=.5 must be treated as a distinct model scenario.
    rows+=scan("pooled_cost",.5,pooled=True)
    output=OUT/"referee_2d_entry_sensitivity.csv"
    with output.open("w",newline="",encoding="utf-8") as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    fig=OUT/"referee_lambda_fixed_cost.svg"
    draw_phase_svg([r for r in rows if r["family"]=="lambda_cost"],fig)
    reset()
    for lab,lam,pool,exp in (("baseline_zero",0.,False,"BP"),
                              ("baseline_half",.5,False,"BB"),
                              ("pooled_half",.5,True,"BP")):
        profiles,nes,plans,sw=s.solve(lam,pool)
        assert "".join(nes[0])==exp and plans==[("B","B")],(lab,nes,plans)
    print(json.dumps({"rows":len(rows),"grid_cells_per_scenario":22,
       "lambda_values":6,"delta_B_values":5,"N_B_values":3,
       "chi_values":3,"pool_comparison":1,"original_baseline_restored":True,
       "interpretation":"Selected baseline numerical continuations only; not parameter-global price-Nash proof"},indent=2))
    print("PASS: 2D sensitivity phase and unchanged frozen baseline")
if __name__=="__main__":main()
