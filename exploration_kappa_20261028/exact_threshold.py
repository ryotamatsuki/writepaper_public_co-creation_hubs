#!/usr/bin/env python3
"""Exact rational certificate of the fixed opposite-full BP/BB threshold.

This proves the conditional threshold, not all sixteen entry equilibria for
every kappa, and does not certify the P/P price game for kappa<1.
"""
from fractions import Fraction as F
from pathlib import Path
import sys,json,math
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'hotelling_entry_experimentation/production/reviews'))
import stage11r_pp_global_certificate as c
from moments import moments

def add(a,b):return [(a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0)) for i in range(max(len(a),len(b)))]
def scale(a,s):return [x*s for x in a]
def mul(a,b):
    q=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):q[i+j]+=x*y
    return q
def deriv(a):return [i*a[i] for i in range(1,len(a))]
def ev(a,x):return sum(t*x**i for i,t in enumerate(a))
def bernstein(a):
    n=len(a)-1
    return [sum(a[j]*F(math.comb(k,j),math.comb(n,j)) for j in range(k+1)) for k in range(n+1)]

def main():
    p=c.primitives(F(1,2),False);i=1
    U=c.basic(i,'B')+p['a'][i];h=p['b'][i]/2
    assert 0<U<1 and 0<h<1
    # Exclude both noncentral branches uniformly in kappa by crude bounds.
    # Lower strip: revenue<=U+h. Upper strip: revenue<=max fee * max tail mass.
    common_lower=(1+U)**2/4
    assert common_lower>U+h
    assert common_lower>(U+h+1)*h/(1+h)
    B=[1+h,-h] # A=B/(1+h), D=1.
    A=[x/(1+h) for x in B]
    numer=add([1+h/2,-h/2],scale(A,-U)) # x=(1+c-AU)/(2A)
    # h<x<1 follows from the two affine numerators' endpoint signs.
    assert all(ev(add(numer,scale(A,-2*h)),k)>0 for k in (0,1))
    assert all(ev(add(scale(A,2),scale(numer,-1)),k)>0 for k in (0,1))
    def threshold(k):
        aa=ev(A,k);f=(1+(1-k)*h/2+aa*U)/(2*aa)
        m,q,r=moments(U,h,F(1),f,k)
        S0=p['pb']*p['m'][i]+p['omega'][i]*p['r0']
        S1=p['lam']*(p['m'][i]+p['omega'][i]*p['pi'])
        return p['n'][i]*(c.basic(i,'B')*(1-m)+F(3,8)-r)+F(3,10)*((1-m)*S0+(F(1,2)-q)*S1/2)
    # Direct polynomial construction, not interpolation or fitted regression.
    # x=X/B; B and X are affine in kappa. M is the affine mass m.
    X=scale(numer,(1+h)/2)
    M=add([1+h/2,-h/2],scale(X,-1/(1+h)))
    den=mul(B,B);X2=mul(X,X);XB=mul(X,B)
    diagonal=scale(add(den,scale(X2,-1/(1+h)**2)),F(1,2))
    indep_q=add(scale(den,F(1,2)+h/3),scale(XB,-F(1,2)))
    indep_r=add(add(scale(den,F(1,2)-h*h/6),scale(X2,-F(1,2))),scale(XB,h/2))
    Q=add(mul([0,1],diagonal),mul([1,-1],indep_q))
    R=add(mul([0,1],diagonal),mul([1,-1],indep_r))
    missing_mass=mul(add([1],scale(M,-1)),den)
    S0=p['pb']*p['m'][i]+p['omega'][i]*p['r0']
    S1=p['lam']*(p['m'][i]+p['omega'][i]*p['pi'])
    direct=add(add(scale(missing_mass,c.basic(i,'B')),scale(den,F(3,8))),scale(R,-1))
    research=add(scale(missing_mass,S0),scale(add(scale(den,F(1,2)),scale(Q,-1)),S1/2))
    num=add(scale(direct,p['n'][i]),scale(research,F(3,10)))
    assert len(num)<=4 # every term is at most affine times affine squared
    # Independent exact rational evaluations catch implementation mistakes.
    for k in (F(0),F(1,7),F(1,3),F(3,5),F(2,3),F(9,10),F(1)):
        assert ev(num,k)/ev(den,k)==threshold(k)
    derivative_num=add(mul(deriv(num),den),scale(mul(num,deriv(den)),-1))
    assert all(x>0 for x in bernstein(derivative_num))
    lo,hi=F(0),F(1);target=p['fg'][1]
    assert threshold(lo)<target<threshold(hi)
    for _ in range(55):
        mid=(lo+hi)/2
        if threshold(mid)<target:lo=mid
        else:hi=mid
    out={'scope':'conditional B-government BP-to-BB comparison; opposite A already B; exact global conditional fee optimum; not full-game uniform-kappa certificate',
         'strictly_increasing_on_0_1':True,
         'proof':'Threshold rational numerator degree <=3, positive denominator. Derivative numerator has strictly positive Bernstein coefficients on [0,1].',
         'threshold_at_kappa_0':float(threshold(F(0))),'threshold_at_kappa_1':float(threshold(F(1))),
         'kappa_crossing_exact_bracket':[str(lo),str(hi)],'kappa_crossing_float_bracket':[float(lo),float(hi)],
         'F_GB':float(target),'derivative_numerator_bernstein_coefficients':[str(x) for x in bernstein(derivative_num)],
         'threshold_numerator_coefficients':[str(x) for x in num],'threshold_denominator_coefficients':[str(x) for x in den]}
    (Path(__file__).parent/'exact_threshold.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS exact conditional global fee, strictly increasing threshold, unique crossing',out['kappa_crossing_float_bracket'])
if __name__=='__main__':main()
