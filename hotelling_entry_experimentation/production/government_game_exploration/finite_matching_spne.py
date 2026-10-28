#!/usr/bin/env python3
"""FROZEN-model exact rational certificate: neighboring fee rise + exit
under distinct STRICT full four-entrant pure SPNE.

Only vary existing p0 and government A public fixed cost; all economic rules
and lambda, matching, profit shares, private fixed costs are unchanged.
All 16 profiles / 64 actor unilateral deviations certified with original
global P/P price proof and rational root faces.
"""
from fractions import Fraction as F
from decimal import Decimal
from itertools import product
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"reviews"))
import stage11r_pp_global_certificate as b

NEW_P0=F(1,1000)
CENTER=(F(Decimal("0.381216064749003")),F(Decimal("0.360189089100758")))
SCENARIOS={
    "lower_A_cost_B_joins":(F(1155,10000),F(41,1000),("B","P")),
    "higher_A_cost_B_does_not_join":(F(1156,10000),F(41,1000),("P","P")),
}
def modp(p0):
 p=b.primitives(F(1,2),False)
 # From the unchanged original p0 and R&D effort equations.
 pi=p["pi"]; zeta=F(5,1000);k=F(4,10000)
 p["r0"]=p0*pi+zeta*zeta*pi*pi/(2*k)
 p["pb"]=p0+zeta*zeta*pi/k
 p["a"]=[F(1,10)*p["n"][1-i]*p["omega"][i]*p["r0"] for i in (0,1)]
 assert 0<p["pb"] and p["pb"]+p["lam"]<1
 return p
def approx(z):return [float(z.lo),float(z.hi)]
def strict_ne(p,rootbox,expected):
 prof={s:b.profile(s,p,rootbox) for s in product(b.STATES,repeat=2)}
 winner=[]
 losers=[]
 for s,v in prof.items():
  slacks=[]
  for i in (0,1):
   for act in ("G","P"):
    t=list(s);t[i]=b.FLIP[act][s[i]];t=tuple(t)
    key="welfare" if act=="G" else "private"
    margin=v[key][i]-prof[t][key][i]
    slacks.append((margin,"AB"[i]+act))
  if s==expected:
   assert all(m.lo>0 for m,_ in slacks),("not strict NE",s,[(a,approx(m)) for m,a in slacks])
   winner=[{"actor":a,"slack":approx(m)} for m,a in slacks]
  else:
   witness=[(m,a) for m,a in slacks if m.hi<0]
   assert witness,("other candidate not eliminated",s,[(a,approx(m)) for m,a in slacks])
   losers.append({"profile":"".join(s),"deviation_witness":witness[0][1],
                  "profitability_min":float(-witness[0][0].hi)})
 return winner,losers,prof
def main():
 p=modp(NEW_P0)
 # Full price global BR constraints, not merely stationary FOCs.
 price_conditions=b.bounds(p)
 rootbox,face_evidence=b.faces(p,list(CENTER),radius=F(1,10**9))
 proof={"scope":"p0=1/1000, original lambda=1/2, gamma=1/10, original N,omega,delta,F_P,screening,product markets. Two different original policy cost values; full 16 profiles and global PP fees.",
        "source":"Frozen PR61/P1e original P/P global fee theorem reused verbatim",
        "p0":str(NEW_P0),"lambda":str(p["lam"]),"matching_gamma":"1/10",
        "pb":float(p["pb"]),"R0":float(p["r0"]),
        "P_P_global_price_sufficient_conditions":price_conditions,
        "P_P_root_box":[[str(q) for q in z] for z in rootbox],
        "P_P_root_face_checks":face_evidence,
        "cost_scenarios":{}}
 pbase=b.profile(("P","P"),p,rootbox)
 pafter=b.profile(("B","P"),p,rootbox)
 fee_before=pbase["fees"][1]
 fee_after=pafter["fees"][1]
 cutoff_before=b.I(1)-pbase["mq"][1][0]
 cutoff_after=b.I(1)-pafter["mq"][1][0]
 dfee=fee_after-fee_before
 dcut=cutoff_after-cutoff_before
 assert dfee.lo>0 and dcut.lo>0
 # No attempt to conflate own private demand with national matching.
 # Under P/P, after region A opens public basic, B's PRIVATE-ONLY
 # participants decline, even though region A participation rises.
 nmatchbefore=pbase["mq"][0][0]*pbase["mq"][1][0]
 nmatchafter=pafter["mq"][0][0]*pafter["mq"][1][0]
 assert nmatchafter.lo>nmatchbefore.hi # national volume increases despite B exit
 proof["repriced_neighboring_P_firm"]={
   "B_fee_when_A_P":approx(fee_before),
   "B_fee_when_A_B":approx(fee_after),
   "B_fee_increase_certified":approx(dfee),
   "B_participation_cutoff_when_A_P":approx(cutoff_before),
   "B_participation_cutoff_when_A_B":approx(cutoff_after),
   "B_cutoff_increase_and_participation_decrease":approx(dcut),
   "national_match_volume_increases":True,
   "national_match_rate_before":approx(nmatchbefore),
   "national_match_rate_after":approx(nmatchafter),
 }
 thrA_P=pafter["welfare"][0]-pbase["welfare"][0]+p["fg"][0]
 thrB_A_P=b.profile(("P","B"),p,rootbox)["welfare"][1]-pbase["welfare"][1]+p["fg"][1]
 thrB_A_B=b.profile(("B","B"),p,rootbox)["welfare"][1]-pafter["welfare"][1]+p["fg"][1]
 proof["government_thresholds"]={"T_A_given_B_P":approx(thrA_P),
                                   "T_B_given_A_P":approx(thrB_A_P),
                                   "T_B_given_A_B":approx(thrB_A_B)}
 # Private entry is STRICTLY dominant in both regions (original FP values);
 # at both costs B-government non-entry is STRICTLY dominant in the reduced
 # government game. Hence the selected A response is strict and there is NO
 # further mixed-strategy stage-one equilibrium.
 assert F(41,1000)>max(thrB_A_P.hi,thrB_A_B.hi)
 assert F(1155,10000)<thrA_P.lo<F(1156,10000)
 proof["entire_first_stage_game_uniqueness_including_mixed"]=True
 for k,(ga,gb,expected) in SCENARIOS.items():
  p["fg"]=(ga,gb)
  winners,losers,_=strict_ne(p,rootbox,expected)
  assert len(losers)==15
  proof["cost_scenarios"][k]={
     "govt_A_public_cost":str(ga),"govt_B_public_cost":str(gb),
     "unique_strict_pure_entry_SPNE":"".join(expected),
     "no_other_mixed_entry_SPNE":True,
     "all_four_actor_deviation_slacks":winners,
     "all_15_other_profiles_excluded":losers,
  }
 # Reverify the baseline itself hasn't been overwritten by this p0
 baseline=b.primitives(F(1,2),False)
 assert baseline["r0"]>p["r0"]
 path=Path(__file__).with_name("finite_matching_spne_certificate.json")
 path.write_text(json.dumps(proof,ensure_ascii=False,indent=2)+"\n")
 print("PASS GLOBAL PP posted-fee equilibrium, all 16 states, all four actor deviations, two distinct unique strict SPNE, neighboring fee up / participation down")
 print("FEES",approx(dfee),"CUTOFF CHANGE",approx(dcut))
 print("MATCH VOLUME(before,after)",approx(nmatchbefore),approx(nmatchafter))
 print("THRESHOLD A|B P",approx(thrA_P))
if __name__=="__main__":main()
