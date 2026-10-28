#!/usr/bin/env python3
"""Independent Stage 11 audit of Stage 10's exported 48-case payoff panel.

Rechecks all four unilateral binary entry deviations and the restricted
16-profile planner without consuming the production script's NE/planner flags.
Also computes exact-from-CSV local/national public-addition margins and a
cost-only perturbation. Requires Python 3.12+ standard library; no external IO.
It does NOT independently construct private/private tariff equilibria.
"""
import csv
import json
import sys
from pathlib import Path

DEFAULT=Path(__file__).resolve().parents[1] / "output" / "stage33_profiles.csv"
CHANGE={
    "G":{"0":"G","G":"0","P":"B","B":"P"},
    "P":{"0":"P","P":"0","G":"B","B":"G"}
}
COST_BASE=0.0405

def main(filename):
    with open(filename,newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==48,("Expected 48 profile rows",len(rows))
    result=[]
    for lam,mode,expected in [
        (0.0,"skill",("B","P")),
        (0.5,"skill",("B","B")),
        (0.5,"pooled",("B","P"))
    ]:
        dat={(r["entry_A"],r["entry_B"]):r for r in rows
             if float(r["lambda"])==lam and r["scenario"]==mode}
        assert len(dat)==16
        def wealth(s,i,cost=COST_BASE):
            r=dat[s]
            return float(r["W_"+("A" if i==0 else "B")]) - (
                cost-COST_BASE if i==1 and s[1] in ("G","B") else 0.0)
        def private(s,i):
            return float(dat[s]["private_"+("A" if i==0 else "B")])
        def status(cost=COST_BASE):
            ne=[]
            for s in dat:
                stable=True
                for i in (0,1):
                    for a in ("G","P"):
                        t=list(s)
                        t[i]=CHANGE[a][s[i]]
                        t=tuple(t)
                        own=wealth(s,i,cost) if a=="G" else private(s,i)
                        devi=wealth(t,i,cost) if a=="G" else private(t,i)
                        if devi-own>1e-9:
                            stable=False
                if stable:
                    ne.append(s)
            planner_score=max(wealth(s,0,cost)+wealth(s,1,cost) for s in dat)
            planner=[s for s in dat if abs(wealth(s,0,cost)+wealth(s,1,cost)-planner_score)<1e-9]
            return sorted(ne),sorted(planner)
        ne,planner=status()
        assert ne==[expected],(mode,lam,ne)
        assert planner==[("B","B")],(mode,lam,planner)
        assert ne==sorted(s for s in dat if int(dat[s]["is_NE"]))
        assert planner==sorted(s for s in dat if int(dat[s]["is_planner"]))
        s_bp=("B","P")
        s_bb=("B","B")
        local=wealth(s_bb,1)-wealth(s_bp,1)
        foreign=wealth(s_bb,0)-wealth(s_bp,0)
        local_star=COST_BASE+local
        national_star=COST_BASE+local+foreign
        def min_slack(s):
            slacks=[]
            for i in (0,1):
                for a in ("G","P"):
                    t=list(s);t[i]=CHANGE[a][s[i]];t=tuple(t)
                    own=wealth(s,i) if a=="G" else private(s,i)
                    dev=wealth(t,i) if a=="G" else private(t,i)
                    slacks.append((own-dev,i,a))
            return min(slacks)
        result.append(dict(lambda_=lam,scenario=mode,NE=ne,restricted_planner=planner,
            weakest_one_flip_margin=min_slack(expected)[0],
            weakest_one_flip_actor="AB"[min_slack(expected)[1]]+"_"+min_slack(expected)[2],
            local_B_addition=local,foreign_A_incidence=foreign,
            national_addition=local+foreign,
            local_B_FG_threshold=local_star,national_B_FG_threshold=national_star,
            full_BB_aggregate_welfare=wealth(s_bb,0)+wealth(s_bb,1),
            NE_at_FG_B_plus_0001=status(COST_BASE+0.0001)[0]))
    skill=result[1]
    pooled=result[2]
    interval=(pooled["local_B_FG_threshold"],skill["local_B_FG_threshold"])
    assert interval[0]<COST_BASE<interval[1]
    assert abs(skill["local_B_addition"]-0.00003167010245)<1e-9
    assert abs(pooled["local_B_addition"]+0.00045677076684)<1e-9
    print(json.dumps({"scenarios":result,"distinct_ownership_cost_interval":interval,
          "interval_width":interval[1]-interval[0],
          "baseline_cost":COST_BASE,
          "upper_distance_as_percent_of_cost":100*(interval[1]-COST_BASE)/COST_BASE,
          "interpretation":"CSV-level independent entry/planner audit only; no independent P/P equilibrium proof"},
         ensure_ascii=False,indent=2))
    print("PASS: 48-row independent entry/planner and cost fragility audit")

if __name__=="__main__":
    main(Path(sys.argv[1]) if len(sys.argv)>1 else DEFAULT)
