#!/usr/bin/env python3
"""PR61-frozen model: exact four-entrant reduction and government best-response cycle.

No primitives, manuscripts, Lean files, original baselines, or fee algorithms are altered.
This reuses Stage11R's already proved global P/P price-equilibrium certificate.
The exact rational boxes below are STRICT subrectangles of the government-cost
threshold intervals; claims apply only to the original skill lambda=1/2 case.
"""
from fractions import Fraction as F
from decimal import Decimal
from itertools import product
from pathlib import Path
import json, sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "reviews"))
import stage11r_pp_global_certificate as b

CENTERS = {
  "skill_zero": (F(0), False, ("0.380848479444378", "0.359890064348462")),
  "skill_half": (F(1,2), False, ("0.381077381587935", "0.360088394043391")),
  "pooled_half": (F(1,2), True, ("0.380493572505174", "0.359632374104520")),
}
COST_BOX = {
  "A": (F(11556,100000),F(11563,100000)), # (.11556,.11563)
  "B": (F(4054,100000), F(4057,100000)),   # (.04054,.04057)
}
def iv(x):
    return {"lower":str(x.lo),"upper":str(x.hi),
            "lower_approx":float(x.lo),"upper_approx":float(x.hi)}
def strictly_pos(x):return x.lo>0
def strictly_neg(x):return x.hi<0
def thresholds(p,price_box):
    states = [("P","P"),("B","P"),("P","B"),("B","B")]
    profiles = {s:b.profile(s,p,price_box) for s in states}
    def w(s,i):return profiles[s]["welfare"][i]
    A_P=w(("B","P"),0)-w(("P","P"),0)+p["fg"][0]
    A_B=w(("B","B"),0)-w(("P","B"),0)+p["fg"][0]
    B_P=w(("P","B"),1)-w(("P","P"),1)+p["fg"][1]
    B_B=w(("B","B"),1)-w(("B","P"),1)+p["fg"][1]
    return {"A_P":A_P,"A_B":A_B,"B_P":B_P,"B_B":B_B},profiles

def uniform_price_contraction():
    """Conservative rational envelope across all skill lambda in [0, 1/2]."""
    p=b.primitives(F(1,2),False)
    a,pb,u=p["a"],p["b"],p["u"]
    h=[a[i]+pb[i] for i in range(2)]
    z=h[0]*h[1]
    assert 0<z<1
    rows=[]
    for i in range(2):
        j=1-i
        lower=(1-z-u[i]-a[i])/(2+pb[i]/2)
        upper=(1+pb[i]/2-u[i])/(2+pb[i]/2)
        Lips=h[i]+(1+pb[i]/2)*(pb[i]+pb[i]*h[j]+pb[j]*h[i]+h[i])/(2*(1-z))
        assert 0<a[i]<lower<upper<1
        assert a[i]+(1+pb[i]/2)*upper<1
        assert 0<Lips<F(1,100)
        rows.append({"region":"AB"[i],"uniform_cutoff_lower":str(lower),
                    "uniform_cutoff_upper":str(upper),
                    "uniform_price_BR_lipschitz_upper":str(Lips),
                    "upper_approx":float(Lips)})
    return rows

def main():
    results={}
    for label,(lam,pooled,center) in CENTERS.items():
        p=b.primitives(lam,pooled)
        # Independent analytical global-price-BR sufficient inequalities and
        # exact root-face enclosure previously checked at the frozen calibration.
        b.bounds(p)
        box,faces=b.faces(p,[F(Decimal(x)) for x in center])
        t,profiles=thresholds(p,box)
        # Strong private-entry dominance for ANY other first-stage actions:
        # with own G present, fee receipts=N_i*Delta_i/4.
        # without own G, fee u_i+Delta_i/2 guarantees mass>=1/2
        # since u_i>0, opposite expected benefits and success increments >=0.
        for i in (0,1):
            assert p["fp"][i]<p["n"][i]/4
            assert b.basic(i,"0")>0
        assert strictly_neg(t["B_B"]-t["B_P"])
        if label=="skill_zero":
            assert strictly_neg(t["A_B"]-t["A_P"])
            # Both governments are strategic substitutes in this scenario:
            # on this open, exactly rational cost box there are TWO strict
            # asymmetric pure equilibria (BP and PB) and one mixed equilibrium.
            aa=(F(1147741,10000000),F(1147745,10000000))
            bb=(F(4043,100000),F(4048,100000))
            assert t["A_B"].hi<aa[0]<aa[1]<t["A_P"].lo
            assert t["B_B"].hi<bb[0]<bb[1]<t["B_P"].lo
            results[label]=results.get(label,{})
            results[label]["two_pure_equilibria_cost_rectangle"]={
                "A":[str(v) for v in aa],"B":[str(v) for v in bb],
                "strict_pure_equilibria":["BP","PB"],
                "additional_mixed_equilibrium":True,
            }
        else:
            assert strictly_pos(t["A_B"]-t["A_P"])
        results[label].update({
          "thresholds":{k:iv(x) for k,x in t.items()},
          "A_cross_effect":iv(t["A_B"]-t["A_P"]),
          "B_cross_effect":iv(t["B_B"]-t["B_P"]),
          "P_P_root_boxes":[[str(z) for z in v] for v in box],
          "private_entry_costs":[str(v) for v in p["fp"]],
          "private_entry_dominance_bounds":[str(v/4) for v in p["n"]],
        })
        if label=="skill_half":
            # A best response to B: entry only when B enters.
            # B best response to A: entry only when A does not enter.
            assert t["A_P"].hi<COST_BOX["A"][0]<COST_BOX["A"][1]<t["A_B"].lo
            assert t["B_B"].hi<COST_BOX["B"][0]<COST_BOX["B"][1]<t["B_P"].lo
            fa=(COST_BOX["A"][0]+COST_BOX["A"][1])/2
            fb=(COST_BOX["B"][0]+COST_BOX["B"][1])/2
            prob_B=(b.I(fa)-t["A_P"])/(t["A_B"]-t["A_P"])
            prob_A=(b.I(fb)-t["B_P"])/(t["B_B"]-t["B_P"])
            assert 0<prob_A.lo<prob_A.hi<1
            assert 0<prob_B.lo<prob_B.hi<1
            results[label]["certified_no_pure_equilibrium_cost_rectangle"]={
              k:[str(z) for z in v] for k,v in COST_BOX.items()
            }
            results[label]["cost_midpoint"]={"A":str(fa),"B":str(fb)}
            results[label]["unique_mixed_government_equilibrium_midpoint"]={
              "prob_A_public":iv(prob_A),"prob_B_public":iv(prob_B)
            }
            # The restricted social planner's sixteen-profile comparison is
            # also strict throughout the entire analytic cost rectangle.
            # Build each welfare DIFFERENCE before interval evaluation so
            # shared public-cost terms cancel exactly (dependency-aware).
            planner=("P","B")
            planner_old=b.profile(planner,p,box)
            gaps={}
            for st in product(b.STATES,repeat=2):
                if st==planner:continue
                old=b.profile(st,p,box)
                gap=(planner_old["welfare"][0]+planner_old["welfare"][1]
                    -old["welfare"][0]-old["welfare"][1])
                for i,k in enumerate(("A","B")):
                    govt_diff=int(planner[i] in ("G","B"))-int(st[i] in ("G","B"))
                    if govt_diff:
                        lo,hi=COST_BOX[k]
                        gap=gap+govt_diff*(p["fg"][i]-b.I(lo,hi))
                assert gap.lo>0,("restricted planner not P/B",st,gap.approx())
                gaps["".join(st)]=iv(gap)
            results[label]["restricted_planner"]={
                "selected":"PB",
                "strictly_beats_all_15_alternatives_in_cost_rectangle":True,
                "minimum_welfare_gap_lower_bound":min(x["lower_approx"] for x in gaps.values()),
            }
            # Strict welfare inefficiency of the unique mixed equilibrium:
            # its support contains non-planner profiles with positive probability.
            one=b.I(1)
            costs_mid={"A":fa,"B":fb}
            def welfare_mid(st):
                r=b.profile(st,p,box)
                v=r["welfare"][0]+r["welfare"][1]
                for i,k in enumerate(("A","B")):
                    if st[i] in ("G","B"):
                        v=v+p["fg"][i]-costs_mid[k]
                return v
            weights={
              ("P","P"):(one-prob_A)*(one-prob_B),
              ("P","B"):(one-prob_A)*prob_B,
              ("B","P"):prob_A*(one-prob_B),
              ("B","B"):prob_A*prob_B,
            }
            expected=sum((welfare_mid(st)*pr for st,pr in weights.items()),b.I(0))
            loss=welfare_mid(planner)-expected
            assert loss.lo>0
            results[label]["midpoint_mixed_vs_planner_welfare_loss"]=iv(loss)
    results["uniform_P_P_price_contraction_skill_lambda_0_to_half"]=uniform_price_contraction()
    # Strict A cross-sign reversal at lambda=0 and lambda=1/2, plus
    # uniform continuity/uniqueness, implies >=1 zero by intermediate value
    # theorem (not uniqueness or characterization of all lambda equilibria).
    results["sign_change_statement"]="A government cross-effect has at least one zero in (0,1/2) under the frozen skill technology; B cross-effect negative at the certified endpoints, not proved negative at all intermediate lambda."
    out=Path(__file__).with_name("government_game_certificate.json")
    out.write_text(json.dumps(results,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("PASS strict private-entry dominance, A cross-sign reversal, B substitutes at all certified endpoints, and open no-pure/mixed-SPNE cost rectangle")
    for label in CENTERS:
        print(label,results[label]["A_cross_effect"]["lower_approx"],results[label]["B_cross_effect"]["upper_approx"])
    print("private costs are original; cost rectangle is an analytical counterfactual, NOT the original cost vector")

if __name__=="__main__":main()
