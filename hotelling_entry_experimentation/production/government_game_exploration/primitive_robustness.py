#!/usr/bin/env python3
"""Certified discrete sensitivity in original model; weak-coupling exact formula."""
from fractions import Fraction as F
from decimal import Decimal
from itertools import product
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reviews'))
import stage11r_pp_global_certificate as b

def primitive(lam,wb):
 p=b.primitives(lam,False);p['omega']=(1-wb,wb)
 p['a']=[F(1,10)*p['n'][1-i]*p['omega'][i]*p['r0'] for i in (0,1)]
 p['b']=[F(1,10)*p['n'][1-i]*p['omega'][i]*lam*p['pi'] for i in (0,1)]
 return p

def weak(lam,wb):
 p=primitive(lam,wb);u=p['u'];ell=[(1-v)/2 for v in u];out=[]
 for i in (0,1):
  j=1-i;li,lj=ell[i],ell[j]
  market=p['pb']*li*lj+lam*li*li*lj*lj/4
  owner=p['r0']*lj*(1-3*u[i])/4+lam*p['pi']*(li*li*lj*lj/4+(u[i]+li)*u[i]*lj*lj/4)
  coeff=p['n'][i]*p['n'][j]*(p['m'][i]*market+p['omega'][i]*owner)
  assert market>0
  out.append(dict(coefficient=str(coeff),approx=float(coeff),
                  market_weight=str(market),owner_weight=str(owner),
                  critical_market_effect=str(-p['omega'][i]*owner/market)))
 return out

def fresh_root(p):
 x=[.381,.360]
 for it in range(10):
  f=[float(t) for t in b.stationary(x,p)]
  if max(map(abs,f))<1e-14:break
  eps=1e-6
  v=[[float(t) for t in b.stationary([x[0]+eps,x[1]],p)],
     [float(t) for t in b.stationary([x[0],x[1]+eps],p)]]
  a=(v[0][0]-f[0])/eps;c=(v[0][1]-f[1])/eps
  d=(v[1][0]-f[0])/eps;e=(v[1][1]-f[1])/eps
  det=a*e-c*d
  assert abs(det)>1e-5
  x[0]-=(e*f[0]-d*f[1])/det
  x[1]-=(a*f[1]-c*f[0])/det
  assert 0<x[0]<1 and 0<x[1]<1
 assert max(abs(float(z)) for z in b.stationary(x,p))<1e-12
 return [F(Decimal(f'{v:.15f}')) for v in x]

def main():
 q0=weak(F(0),F(1,5));qh=weak(F(1,2),F(1,5))
 assert F(q0[0]['coefficient'])<0<F(qh[0]['coefficient'])
 assert F(q0[1]['coefficient'])<0 and F(qh[1]['coefficient'])<0
 roots=[]
 for i in (0,1):
  a=F(q0[i]['coefficient']);d=2*(F(qh[i]['coefficient'])-a)
  assert d>0
  roots.append(dict(region='AB'[i],exact=str(-a/d),approx=float(-a/d)))
 assert 0<F(roots[0]['exact'])<F(1,2)<F(roots[1]['exact'])
 rows=[]
 for lam,wb in product((F(2,5),F(1,2),F(3,5)),(F(9,50),F(1,5),F(11,50))):
  p=primitive(lam,wb);assert p['pb']+lam<1
  assert all(p['fp'][i]<p['n'][i]/4 for i in (0,1))
  sufficient=b.bounds(p)
  box,_=b.faces(p,fresh_root(p))
  pp=b.profile(('P','P'),p,box);bp=b.profile(('B','P'),p,box)
  pb=b.profile(('P','B'),p,box);bb=b.profile(('B','B'),p,box)
  ga=bb['welfare'][0]+pp['welfare'][0]-bp['welfare'][0]-pb['welfare'][0]
  gb=bb['welfare'][1]+pp['welfare'][1]-bp['welfare'][1]-pb['welfare'][1]
  assert ga.lo>0 and gb.hi<0,('fail',lam,wb,ga.approx(),gb.approx())
  rows.append(dict(lam=str(lam),omega_B=str(wb),A=[float(ga.lo),float(ga.hi)],
                   B=[float(gb.lo),float(gb.hi)],cost_area_lower=float(ga.lo*(-gb.hi)),
                   global_fee_BRL_upper=max(t['br_contraction_bound'] for t in sufficient)))
 data=dict(scope='Nine individually certified rational primitive points, not a continuous primitive box; weak-coupling coefficient exact, not full-gamma sign theorem',
    weak_gamma={'lambda_zero':q0,'lambda_half':qh,'crossings':roots},points=rows,
    summary={'number':len(rows),'min_cost_area_lower':min(r['cost_area_lower'] for r in rows),
             'max_fee_BRL_upper':max(r['global_fee_BRL_upper'] for r in rows)})
 Path(__file__).with_name('primitive_robustness_results.json').write_text(json.dumps(data,indent=2)+'\n')
 print('PASS weak-gamma primitive thresholds; 9 original-model P/P global price and asymmetric government BR signatures')
 print('lambda roots:',[z['approx'] for z in roots])
if __name__=='__main__': main()
