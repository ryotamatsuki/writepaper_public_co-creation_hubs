#!/usr/bin/env python3
"""P1g audit: independent exact rational checks of three frozen-model deductions.
No new state variables, technologies, entrant objective, policy or price rule.
Checks analytical identities and signs, not a full-price-game new theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"reviews"))
import stage11r_pp_global_certificate as core

def H(u,d,ell,N=1):return N*(u*ell-d/F(8)+d*ell*ell/2)

def public_service_diagnostics():
 d=F(1); u=d
 assert H(u,d,F(0)) == -d/8
 assert H(F(1,4),d,F(0))<0
 assert H(F(1,4),d,F(1,2))>0
 # Strictly increasing in ell: derivative N(u+d ell)>0.
 # Unique root ell* = sqrt((u/d)^2+1/4)-u/d.
 # Exact rational bracketing (u/d=1/4): root between 3/10 and 1/3.
 assert H(F(1,4),d,F(3,10))<0
 assert H(F(1,4),d,F(1,3))>0
 return dict(full_private_access_public_option_service_loss_exact=str(H(F(1),F(1),F(0))),
   root_interval_rational=["3/10","1/3"],
   strict_dH_dell="N*(u+d*ell)>0 for u,d,N positive",
   scope="Partner fully covered; at cutoff zero even partner private endogenously same because m,q fixed. Gross direct service only. Matching terms vanish if local participation unchanged.")

def private_receipts():
 p=core.primitives(F(1,2),False)
 out=[]
 for i in (0,1):
  d=F(1);n=p['n'][i];u=core.basic(i,"0")
  floor=n*(u+d/2)/2
  B=n*d/4
  assert 0<B<floor
  assert p['fp'][i]<B
  # Without a free own public service, fee f=u+d/2 guarantees
  # participation theta>=1/2 for every nonnegative m_j,q_j,
  # including when opponent users re-equilibrate after a deviation.
  out.append(dict(region="AB"[i],d=str(d),N=str(n),u_min=str(u),
    revenue_floor_without_own_public=str(floor),
    premium_revenue_with_public=str(B),
    private_fixed_cost=str(p['fp'][i]),
    strictly_less_revenue_with_free_public=True))
 return out

def cross_border_crowding_sign(p0):
 p=core.primitives(F(1,2),False)
 pi=p['pi']
 # Existing same Cournot and R&D primitives. No change in model form.
 # zeta=1/200; k=1/2500 -> zeta^2/(2k)=1/32.
 r0=p0*pi+pi*pi/F(32)
 pb=p0+pi/F(16)
 lam=F(1,2);ui=p['u'][0];uj=p['u'][1];d=F(1)
 ell_i=(d-ui)/(2*d)
 numerator=lam*pi*uj*ell_i-2*r0*d
 assert 0<pb+lam<1
 # Cross-region P-only region j: after region i adds free G,
 # ell_j(BP)-ell_j(PP)= gamma*N_i*omega_j*ell_i/(4*d^2)
 #  * numerator + O(gamma^2). Positive => FEWER j participants.
 return dict(p0=str(p0),R0=str(r0),pb=str(pb),
  ell_i_zero_gamma=str(ell_i),
  first_order_sign_numerator=str(numerator),
  first_order_sign="j-participation-falls" if numerator>0 else "j-participation-rises",
  feasibility_probability_success="p_b+lambda<1")

def main():
 direct=public_service_diagnostics()
 receipts=private_receipts()
 baseline=cross_border_crowding_sign(F(1,20))
 possibility=cross_border_crowding_sign(F(1,1000))
 assert F(baseline["first_order_sign_numerator"])<0
 assert F(possibility["first_order_sign_numerator"])>0
 data=dict(scope="Exact identity and weak matching first-order criterion, NOT any new globally certified PP equilibrium for alternative p0",
   direct_service_loss=direct,private_fee_crowding=receipts,
   induced_other_region_participation={"baseline":baseline,
     "valid_existing_primitive_counterexample":possibility,
     "asymptotic_formula":"ell_j(BP)-ell_j(PP) = gamma*N_i*omega_j*ell_i/(4*d_j^2) * (lambda*Pi*u_j*ell_i-2*R0*d_j) + O(gamma^2), where i switches P to B and both remained P before",
     "assumptions":"0<u_i<d_i;0<u_j<d_j; omega_j>0; unique smooth inner PP Nash near gamma=0; rho held fixed as mu tends to zero; all successful probability bounds; no new primitive parameter."})
 out=Path(__file__).with_name("additional_identities.json")
 out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
 print("PASS strict public service-quality loss at full private coverage")
 print("PASS private-receipts floor despite any endogenous other-region participation")
 print("PASS cross-jurisdiction participation response sign can flip in feasible existing model by changing p0")
 print("scope: weak-gamma and conditional service results; full game not re-certified at alternative p0")
if __name__=="__main__":main()
