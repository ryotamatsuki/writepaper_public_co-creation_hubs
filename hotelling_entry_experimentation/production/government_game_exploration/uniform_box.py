#!/usr/bin/env python3
"""True all-parameter interval certificate for small lambda/omega rectangles.

Checks uniform global-price sufficient inequalities, Poincare-Miranda
root-face signs, and strategic cross-effect signs for EACH parameter
in the certified (lambda, omega_B) rectangle, not merely grid points.
Uses only Fraction and rigorously outward rational interval b.I.
"""
from fractions import Fraction as F
from decimal import Decimal
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reviews'))
import stage11r_pp_global_certificate as b

def enclose(rad):
 lam=b.I(F(1,2)-rad,F(1,2)+rad)
 wb=b.I(F(1,5)-rad,F(1,5)+rad)
 p=b.primitives(F(1,2),False)
 p['omega']=(b.I(1)-wb,wb);p['lam']=lam
 p['a']=[F(1,10)*p['n'][1-i]*p['omega'][i]*p['r0'] for i in (0,1)]
 p['b']=[F(1,10)*p['n'][1-i]*p['omega'][i]*lam*p['pi'] for i in (0,1)]
 u,a,bb=p['u'],p['a'],p['b'];h=[a[i]+bb[i] for i in (0,1)]
 z=h[0]*h[1];assert z.lo>0 and z.hi<1
 for i in (0,1):
  j=1-i
  L=(1-z-u[i]-a[i])/(2+bb[i]/2-z)
  U=(1+bb[i]/2-u[i])/(2+bb[i]/2)
  PL=u[i]+L;PH=u[i]+a[i]+(1+bb[i]/2)*U
  Lip=h[i]+(1+bb[i]/2)*(bb[i]+bb[i]*h[j]+bb[j]*h[i]+h[i])/(2*(1-z))
  assert L.lo>0 and U.hi<1 and L.hi<U.lo
  assert (u[i]+a[i]).hi<PL.lo and PL.hi<PH.lo and PH.hi<u[i]+1
  assert Lip.lo>=0 and Lip.hi<1
 # Choose one common cutoff root box for every real primitive in the box.
 # Radius scales mildly with the primitive radius and checked face signs,
 # not just a single observed numerical root.
 r=F(1,10**8)+rad/F(10)
 center=[F(Decimal('0.381077381587935')),F(Decimal('0.360088394043391'))]
 box=[[v-r,v+r] for v in center]
 for i in (0,1):
  for side,expect in ((0,'+'),(1,'-')):
   x=[b.I(k[0],k[1]) for k in box]
   x[i]=b.I(box[i][side])
   s=b.stationary(x,p)[i].sign()
   assert s==expect,("fee face",i,side,s,rad)
 pp=b.profile(('P','P'),p,box);bp=b.profile(('B','P'),p,box)
 pb=b.profile(('P','B'),p,box);both=b.profile(('B','B'),p,box)
 ds=[]
 for i in (0,1):
  d=both['welfare'][i]+pp['welfare'][i]-bp['welfare'][i]-pb['welfare'][i]
  ds.append(d)
 assert ds[0].lo>0 and ds[1].hi<0,("cross",rad,ds[0].approx(),ds[1].approx())
 return {'radius':str(rad),'lambda':[str(lam.lo),str(lam.hi)],
         'omega_B':[str(wb.lo),str(wb.hi)],
         'root_cutoff_boxes':[[str(q) for q in v] for v in box],
         'A_cross_rational':[str(ds[0].lo),str(ds[0].hi)],
         'B_cross_rational':[str(ds[1].lo),str(ds[1].hi)],
         'A_cross_float':ds[0].approx(),'B_cross_float':ds[1].approx()}

def main():
 passed=[];failed=[]
 for rad in [F(1,10**8),F(1,10**7),F(1,10**6),F(1,10**5),F(1,10**4),F(1,10**3),F(1,500),F(1,100)]:
  try:
   result=enclose(rad);passed.append(result)
  except AssertionError as e:
   failed.append({'radius':str(rad),'reason':repr(e)[:500]})
 assert passed
 # Only certified rectangle is claimed; an unsuccessful larger test is not
 # a mathematical counterexample to the underlying property.
 data={'scope':'Uniform continuous primitive rectangle, rational interval enclosure; neither all-lambda nor population-distribution robustness',
       'largest_certified_by_this_enclosure':passed[-1],
       'successful':passed,'uncertified_larger_attempts':failed}
 Path(__file__).with_name('uniform_box_results.json').write_text(json.dumps(data,indent=2)+'\n')
 print('PASS uniform price/entry sign for all lambda and omega_B in a certified continuous rectangle')
 print('largest radius',passed[-1]['radius'])
 print('A interval',passed[-1]['A_cross_float'],'B interval',passed[-1]['B_cross_float'])
 print('failed tests do not disprove larger robustness:',failed)
if __name__=='__main__':main()
