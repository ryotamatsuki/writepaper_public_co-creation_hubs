#!/usr/bin/env python3
"""P1k: independent direct-vs-repriced local policy-comparative-static proof.

No model/primitive changes. Proves for a WHOLE rectangular range
omegaB in [0.13,0.27], lambda in [0,0.5] at ORIGINAL gamma=0.1:
  K_B<0, S_B^0<0, d T_reoptimised/dlambda > 0,
  d T_fixed_fee/dlambda < 0, and their difference >0.
Fees 'fixed' means infinitesimal technology perturbation evaluated at
the current private monopoly fee f*(lambda0); NOT alternative equilibrium.
All signs certified with rational endpoint inclusion over 100 rectangles,
not point grids or floating solvers.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
import incidence_screening_sign_certificate as b

def q(i):return [float(i.lo),float(i.hi)]
def calc():
 m,pi=b.cournot();n_i,n_j=F(1),F(3)
 g=F(1,10);C=g*n_i*n_j;u=F(7,25);d=F(1)
 zeta=F(5,1000);k=F(4,10000);p0=F(1,20)
 pb=p0+zeta*zeta*pi/k;R0=p0*pi+zeta*zeta*pi*pi/(2*k)
 assert m[1]<0 and m[0]+m[1]+pi>0
 minima={name:None for name in ('reprice','negative_fixed','increment','quality','basic_access')}
 extrema={'K_max':None,'S_max':None}
 for w_idx,l_idx in product(range(10),repeat=2):
  w=b.I(F(13,100)+F(14,100)*w_idx/10,F(13,100)+F(14,100)*(w_idx+1)/10)
  lam=b.I(F(l_idx,20),F(l_idx+1,20))
  U=u+g*n_j*w*R0
  eta=g*n_j*w*pi/2
  D=d+eta*lam
  ell=(D-U)/(2*D)
  K=m[1]+w*pi
  S=pb*m[1]+w*R0
  reprice=U*eta/(2*D*D)
  fixed=-ell*eta/D
  assert ell.lo>0 and ell.hi<F(1,2) and K.hi<0 and S.hi<0
  def dt(ep):
   return (n_i*(u+d*ell)+C*S)*ep+(C*K/4)*(ell*ell+2*lam*ell*ep)
  Dreprice=dt(reprice)
  Dfixed=dt(fixed)
  difference=Dreprice-Dfixed
  # Exact decomposition of price-induced cutoff:
  # d ell*/d lambda - (d ell/d lambda)|f=f*(lambda0)
  # = eta/(2D), always >0.
  repricing_cutoff=eta/(2*D)
  # derivative of the value T as a function of current cutoff
  sensitivity=n_i*(u+d*ell)+C*S+C*lam*K*ell/2
  assert repricing_cutoff.lo>0 and sensitivity.lo>0
  assert difference.lo>0 and Dreprice.lo>0 and Dfixed.hi<0
  identity=repricing_cutoff*sensitivity
  # The interval difference via dependency may be wider than isolated identity
  # but must have the same positivity; prove exactly at rational corner points.
  if w_idx in [0,5,9] and l_idx in [0,5,9]:
   for w0,lam0 in product([w.lo,w.hi],[lam.lo,lam.hi]):
    u0=u+g*n_j*w0*R0;eta0=g*n_j*w0*pi/2;D0=d+eta0*lam0
    ell0=(D0-u0)/(2*D0);S0=pb*m[1]+w0*R0;K0=m[1]+w0*pi
    E=u0*eta0/(2*D0*D0);F0=-ell0*eta0/D0
    T=lambda e:(n_i*(u+d*ell0)+C*S0)*e + (C*K0/4)*(ell0*ell0+2*lam0*ell0*e)
    assert T(E)-T(F0)==eta0/(2*D0)*(n_i*(u+d*ell0)+C*S0+C*lam0*K0*ell0/2)
  metrics={'reprice':Dreprice.lo,
           'negative_fixed':-Dfixed.hi,
           'increment':difference.lo,
           'quality':sensitivity.lo,
           'basic_access':repricing_cutoff.lo}
  for key,val in metrics.items():
   minima[key]=val if minima[key] is None else min(minima[key],val)
  for key,val in [('K_max',K.hi),('S_max',S.hi)]:
   extrema[key]=val if extrema[key] is None else max(extrema[key],val)
 # Exact weak-matching coefficient identities verified symbolically by hand:
 # U=u+O(gamma), eta=O(gamma), ell0=(d-u)/(2d).
 # endogenous: M(d-u)^2 + omega Pi(d²+3u²)
 # fixed: M(d-u)^2 - omega Pi(d-u)(d+3u)
 # difference: 2 omega Pi d(d+u)
 # Nonzero basic service value is indispensable for K<0 open share band.
 assert F(13,100)<F(1,5)<F(27,100)
 res={'scope':'Finite gamma original 1/10, conditional partner A fully covered; fixed fee is a local non-equilibrium partial derivative, NOT a new first-stage game. Theorem is NOT full four-player global equilibrium result.',
  'parameter_rectangle':{'omegaB':['13/100','27/100'],'lambda':['0','1/2'],'cells_covering_continuum':100},
  'proved_signs':{'endogenous_fee_policy_threshold_derivative':'>0','holding_private_fee_fixed_policy_threshold_derivative':'<0','endogenous_minus_fixed':'>0','K_B':'<0','S_B0':'<0'},
  'min_rational_lower_bounds':{name:float(val) for name,val in minima.items()},
  'max_rational_upper_bounds':{name:float(val) for name,val in extrema.items()},
  'weak_gamma_exact_leading_numerators':{
  'fee_reoptimised':'M_i*(Delta_i-u_i)^2 + omega_i*Pi*(Delta_i^2+3u_i^2)',
  'fee_fixed':'M_i*(Delta_i-u_i)^2 - omega_i*Pi*(Delta_i-u_i)*(Delta_i+3u_i)',
  'difference':'2*omega_i*Pi*Delta_i*(Delta_i+u_i)'},
  'identity':'dT_reprice/dlambda-dT_fixed/dlambda = eta/(2D)*[N(u+d*ell)+C*S0+C*lambda*K*ell/2], evaluated at fee-optimal ell(lambda).'}
 path=Path(__file__).with_name('repricing_counterfactual_sign_results.json')
 path.write_text(json.dumps(res,ensure_ascii=False,indent=2)+'\n')
 print('PASS 100 exact rational coverage rectangles, sign flip: endogenous T_lambda>0 but fixed-fee T_lambda<0 with K<0 and S0<0')
 print('minima',res['min_rational_lower_bounds'])
 print('maxima',res['max_rational_upper_bounds'])
if __name__=='__main__':calc()
