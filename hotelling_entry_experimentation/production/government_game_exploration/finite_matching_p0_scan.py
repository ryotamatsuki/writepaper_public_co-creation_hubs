#!/usr/bin/env python3
"""Pointwise exact rational finite-gamma SPNE sensitivity in existing p0.
Global P/P BR + rational face + all 16 profiles at each preselected p0.
No continuous p0 box is asserted; no retuning of cost comparisons.
"""
from fractions import Fraction as F
from decimal import Decimal
from pathlib import Path
from itertools import product
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import finite_matching_spne as f
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"reviews"))
import stage11r_pp_global_certificate as b

POINTS=(F(1,1000),F(5,1000),F(10,1000),F(20,1000),F(24,1000),F(30,1000),F(50,1000))
def center(p):
 x=[.381,.360]
 for it in range(15):
  z=b.stationary(x,p);q=[float(v) for v in z]
  if max(map(abs,q))<1e-14:break
  eps=1e-6
  r=b.stationary([x[0]+eps,x[1]],p)
  s=b.stationary([x[0],x[1]+eps],p)
  J00=(float(r[0])-q[0])/eps;J10=(float(r[1])-q[1])/eps
  J01=(float(s[0])-q[0])/eps;J11=(float(s[1])-q[1])/eps
  det=J00*J11-J01*J10
  assert abs(det)>1e-6
  x[0]-=(q[0]*J11-q[1]*J01)/det
  x[1]-=(J00*q[1]-J10*q[0])/det
  assert 0<x[0]<1 and 0<x[1]<1
 assert max(abs(float(q)) for q in b.stationary(x,p))<1e-11
 return [F(Decimal(f"{k:.15f}")) for k in x]

def main():
 results=[]
 for p0 in POINTS:
  p=f.modp(p0); b.bounds(p)
  root,_=b.faces(p,center(p),radius=F(1,10**9))
  pp=b.profile(("P","P"),p,root);bp=b.profile(("B","P"),p,root)
  fee=bp["fees"][1]-pp["fees"][1]
  cutoff=pp["mq"][1][0]-bp["mq"][1][0]
  assert fee.lo>0
  assert cutoff.lo>0 if p0<=F(20,1000) else cutoff.hi<0
  chosen={}
  for tag,(fa,fb,expected) in f.SCENARIOS.items():
   p["fg"]=(fa,fb)
   winner,losers,prof=f.strict_ne(p,root,expected)
   assert len(winner)==4 and len(losers)==15
   chosen[tag]="".join(expected)
  results.append({"p0":str(p0),"fee_increase":[float(fee.lo),float(fee.hi)],
                  "B_participation_cutoff_increase":[float(cutoff.lo),float(cutoff.hi)],
                  "resulting_B_participation_decreases":cutoff.lo>0,
                  "price_equilibrium_root_box":[[str(x) for x in interval] for interval in root],
                  "strict_unique_pure_spne_at_fixed_two_costs":chosen})
  print("EXACT p0",p0,"neighbor fee",results[-1]["fee_increase"],"neighbor cutoff",results[-1]["B_participation_cutoff_increase"],"SPNE",chosen)
 assert len(results)==7
 out=Path(__file__).with_name("finite_matching_p0_scan.json")
 out.write_text(json.dumps({"scope":"Seven exact pointwise p0 certificates, not a continuous p0 range or uniqueness of crossing. Two fixed government cost vectors apply at all p0; original model maintained.",
 "p0_points":results},indent=2)+"\n")
 print("PASS 7 distinct global P/P prices, 16-state SPNE pairs, fee rise at every point, participation sign reversal without new model or tuned government costs")
if __name__=="__main__":main()
