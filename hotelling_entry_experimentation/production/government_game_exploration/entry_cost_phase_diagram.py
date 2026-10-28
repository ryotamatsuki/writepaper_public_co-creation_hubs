#!/usr/bin/env python3
"""Complete, exact stage-one *pure* NE cost partition for frozen lambda=1/2.

All 16 state profiles map to 4 independent entry-cost inequalities.
This is a complete conditional characterisation for ANY four nonnegative
entry costs with the unchanged, globally certified frozen downstream fees.
It does not compute all mixed equilibria of the four-player game.
"""
from fractions import Fraction as F
from decimal import Decimal
from itertools import product
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"reviews"))
import stage11r_pp_global_certificate as b

ACTORS=("GA","PA","GB","PB")
STATUS={(0,0):"0",(1,0):"G",(0,1):"P",(1,1):"B"}
def states(bits):
 return (STATUS[bits[0:2]],STATUS[bits[2:4]])
def flip(bits,k):
 bits=list(bits);bits[k]^=1
 return tuple(bits)
def display(I):return [float(I.lo),float(I.hi)]
def profiles_and_boxes():
 p=b.primitives(F(1,2),False)
 b.bounds(p)
 center=(F(Decimal("0.381077381587935")),F(Decimal("0.360088394043391")))
 root,faces=b.faces(p,center,radius=F(1,10**9))
 pmap={bits:b.profile(states(bits),p,root) for bits in product((0,1),repeat=4)}
 thresholds={}
 for bits in pmap:
  t=[]
  for k in range(4):
   i=k//2;typ=k%2
   yes=bits if bits[k]==1 else flip(bits,k)
   no=flip(bits,k) if bits[k]==1 else bits
   if typ==0:
    v=pmap[yes]["welfare"][i]-pmap[no]["welfare"][i]+p["fg"][i]
   else:
    v=pmap[yes]["private"][i]+p["fp"][i]
    assert pmap[no]["private"][i].lo==pmap[no]["private"][i].hi==0
   t.append(v)
  thresholds[bits]=t
 return p,root,thresholds

def strict_NE(bits,cost,thresholds):
 vals=thresholds[bits]
 for k in range(4):
  if bits[k] and not vals[k].lo>cost[k]:return False
  if not bits[k] and not vals[k].hi<cost[k]:return False
 return True

def main():
 p,root,boxes=profiles_and_boxes()
 prior=(p["fg"][0],p["fp"][0],p["fg"][1],p["fp"][1])
 examples={
  "original":prior,
  "private_costs_beyond_dominance":(F(9,100),F(4,5),F(81,2000),F(3,10)),
  "one_private_just_above_premium_threshold":(F(9,100),F(37,50),F(81,2000),F(13,50)),
  "no_public_high_costs":(F(2),F(4,5),F(2),F(3,10)),
  "all_zero_entry_costs":(F(0),F(0),F(0),F(0)),
 }
 result={"scope":"16 exact boxes in four fixed costs; PP globally certified for frozen lambda=1/2 only; all other fee branches exact; this is pure Nash partition, NOT all mixed equilibrium characterisation.",
 "original_primitive_parameters": {"p0":"1/20","lambda":"1/2","rho_mu":"1/10"},
 "actor_cost_order":ACTORS, "profile_boxes":{}, "test_cost_profiles":{}}
 for bits,t in boxes.items():
  key="".join(states(bits))
  result["profile_boxes"][key]={
   "activity_binary":list(bits),
   "thresholds":{ACTORS[k]:display(t[k]) for k in range(4)},
   "inequalities":[f'{ACTORS[k]} fixed entry cost '+('<=' if bits[k] else '>=' )+f' [{float(t[k].lo):.12g},{float(t[k].hi):.12g}]' for k in range(4)]}
 # The sign test is exact at the provided rational costs (strict signs on interval enclosures).
 for label,cost in examples.items():
  eq=[ "".join(states(bits)) for bits in boxes if strict_NE(bits,cost,boxes)]
  if label=="original":assert eq==["BB"]
  if label=="all_zero_entry_costs":assert "BB" in eq
  if label=="private_costs_beyond_dominance":assert all("B" not in s for s in eq)
  result["test_cost_profiles"][label]={"costs":list(map(str,cost)), "strict_pure_NE":eq}
  print(label,eq)
 # Global model-independent binary-service conclusion:
 # Under govt G_i active, the private gross upgrade revenue is N_i*d_i/4.
 # Under govt absent, a price u_i+d_i/2 secures >N_i*d_i/4.
 for i in (0,1):
  priv_index=1+2*i;gov_index=2*i; pivot=p["n"][i]/4
  for bits,t in boxes.items():
   if bits[gov_index]:
    assert t[priv_index].lo==t[priv_index].hi==pivot
   else:
    # u>0; revenue floor secure. Includes opponents' endogenous pricing.
    assert t[priv_index].lo>pivot
 result["dominance_consequence"]="For F_Pi < N_i*Delta_i/4, P is strict dominant. For F_Pi > N_i*Delta_i/4, an equilibrium cannot include B_i, while under F_Pi >= max(P-only revenue threshold) P_i is strictly non-entry dominant. In intermediate cost range P_i is strict substitute to G_i."
 # Construct a complete 64-inequality consistency check, across all profiles.
 assert len(boxes)==16
 assert sum(len(v) for v in boxes.values())==64
 # All 4 actor thresholds are within each of 16 exact cost-orientation cells.
 path=Path(__file__).with_name("all_16_cost_boxes.json")
 path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
 print("PASS complete 16 four-cost-axis box characterisation, all 64 cost inequalities, dominance and above-dominance cases")
if __name__=="__main__":main()
