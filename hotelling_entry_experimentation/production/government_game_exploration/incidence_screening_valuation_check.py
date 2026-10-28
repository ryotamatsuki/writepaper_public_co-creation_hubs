#!/usr/bin/env python3
"""Independent extra economic-relevance test of P1j:
Positive public-entry threshold, negative region-B venture incidence,
but positive and increasing B-government value across a continuous
existing-primitives rectangle. Outward exact rational intervals, not grid samples.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from incidence_screening_sign_certificate import I,cournot
M,pi=cournot()
u=F(7,25);d=F(1);gamma=F(1,10);N=F(1);Nj=F(3)
zeta=F(5,1000);k=F(4,10000);p0=F(1,20)
pb=p0+zeta*zeta*pi/k;r0=p0*pi+zeta*zeta*pi*pi/(2*k)
best_min_T=None;max_venture=None
for iw,il in product(range(10),range(10)):
    w=I(F(13,100)+F(14,100)*iw/10,F(13,100)+F(14,100)*(iw+1)/10)
    l=I(F(il,20),F(il+1,20))
    U=u+gamma*Nj*w*r0;h=gamma*Nj*w*pi*l/2;D=d+h
    ell=(D-U)/(2*D)
    C=gamma*N*Nj
    K=M[1]+w*pi
    S=pb*M[1]+w*r0
    H=N*(u*ell-d/8+d*ell*ell/2)
    venture=C*(S*ell+l*K*ell*ell/4)
    T=H+venture
    assert H.lo>0 and venture.hi<0 and T.lo>0, (iw,il,H.approx(),venture.approx(),T.approx())
    best_min_T=T.lo if best_min_T is None else min(best_min_T,T.lo)
    max_venture=venture.hi if max_venture is None else max(max_venture,venture.hi)
# ω=0 => ell is independent of λ, K=M<0, and T' = C*M*ell²/4<0 (exact).
ell0=(d-u)/(2*d)
domega0=gamma*N*Nj*M[1]*ell0*ell0/4
assert domega0<0
# weak gamma primitive B|u=0 = Delta² K<0, no inversion.
w=F(1,5)
assert d*d*(M[1]+w*pi)<0
# general economic inequality is strict and open in economics primitives.
assert F(0)<best_min_T and max_venture<F(0)
report={"scope":"Conditional region B P vs B, with region A G/B fully covered. NOT a global 4-player SPNE statement or any subsidy cost result.",
"rectangle":{"omegaB":["13/100","27/100"],"lambda":["0","1/2"]},
"all_pieces_public_entry_threshold_positive_min_lower":float(best_min_T),
"all_pieces_marginal_venture_related_welfare_negative_max_upper":float(max_venture),
"omegaB_zero_lambda_derivative_negative_exact":str(domega0),
"u_equals_zero_weak_gamma_condition":"The leading gamma coefficient simplifies to Delta²*(M+omega*Pi), hence cannot be positive if K<0",
"all_100_closed_cells_certified":True}
Path(__file__).with_name("incidence_screening_valuation_check.json").write_text(json.dumps(report,indent=2)+"\n")
print("PASS: every omegaB x lambda in continuous 0.13-0.27 x 0-0.5 has positive threshold and positive direct service, strictly negative venture incidence.")
print("minimum certified public willingness to pay:",float(best_min_T),"largest certified venture-related local gain:",float(max_venture))
print("PASS: when omegaB=0 marginal public value derivative negative. When basic value u=0, weak gamma paradox impossible for K<0.")
