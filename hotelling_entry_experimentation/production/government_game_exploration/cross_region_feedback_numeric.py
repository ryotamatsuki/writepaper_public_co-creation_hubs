#!/usr/bin/env python3
"""Independent no-dependency Newton continuation: P/P fees fully repriced.
Numerical diagnostics only; exact local comparative-static coefficient is
checked separately by additional_identities.py and its derivation report.
"""
from fractions import Fraction as F
import math
PI=F(1656500,14969161)
N=(3.,1.);OMEGA=(.8,.2);U=(.238,.28)
def test(p0,gamma):
 r0=float(F(p0)*PI+PI*PI/F(32)); PIF=float(PI); lam=.5
 a=[gamma*N[1-i]*OMEGA[i]*r0 for i in (0,1)]
 b=[gamma*N[1-i]*OMEGA[i]*lam*PIF for i in (0,1)]
 def foc(l):
  d=[1+b[i]*(1-l[1-i]**2)/2 for i in (0,1)]
  c=[a[i]+b[i]*l[0]*l[1] for i in (0,1)]
  fee=[U[i]+a[i]*(1-l[1-i])+d[i]*l[i] for i in (0,1)]
  return [(1-l[i])*(d[i]-c[0]*c[1]/d[1-i])-fee[i] for i in (0,1)]
 x=[.381,.360]
 for _ in range(12):
  q=foc(x)
  if max(abs(z) for z in q)<1e-13:break
  eps=1e-6
  s=foc([x[0]+eps,x[1]]); t=foc([x[0],x[1]+eps])
  a00=(s[0]-q[0])/eps;a10=(s[1]-q[1])/eps
  a01=(t[0]-q[0])/eps;a11=(t[1]-q[1])/eps
  det=a00*a11-a01*a10
  assert abs(det)>1e-5
  x[0]-=(q[0]*a11-q[1]*a01)/det
  x[1]-=(a00*q[1]-a10*q[0])/det
  assert 0<x[0]<1 and 0<x[1]<1
 assert max(abs(z) for z in foc(x))<1e-12
 # After A switches P->B, its m_A=1, q_A=1/2.
 # Private B's global conditional monopoly optimum is this cutoff.
 U_B=U[1]+a[1];T_B=1+b[1]/2
 cutoff_B=(T_B-U_B)/(2*T_B)
 observed=(cutoff_B-x[1])/gamma
 ell_A=(1-U[0])/2
 predicted=N[0]*OMEGA[1]*ell_A/4*(lam*PIF*U[1]*ell_A-2*r0)
 return observed,predicted

def main():
 for p0,expected in [(F(1,20),-1),(F(1,1000),1)]:
  errors=[]
  for gamma in (1e-3,1e-4,1e-5):
   obs,pred=test(p0,gamma)
   assert obs*expected>0 and pred*expected>0
   errors.append(abs(obs-pred))
   print(f"p0={p0} gamma={gamma} actual_coefficient={obs:.14g} predicted={pred:.14g} abs_error={errors[-1]:.4g}")
  assert errors[2]<errors[1]<errors[0]
  assert errors[-1]<1e-8
 print("PASS independent 2x2 endogenous P/P posted-fee response agrees with derived weak-matching cross-region participation effect")
if __name__=="__main__":main()
