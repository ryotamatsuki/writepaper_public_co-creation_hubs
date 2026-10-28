#!/usr/bin/env python3
"""Stage11R P1j. Exact rational certificate for an incidence-screening policy paradox.

Frozen PR #61 economic primitives and conditional public-addition policy, not a
new model and NOT a first-stage SPNE theorem. Partner jurisdiction is G or B
and therefore full participation; technology lambda in [0,1/2].

Analytical identity:
 T_i(lambda) = N_i(u_i ell_i-Delta_i/8+Delta_i ell_i²/2)
               + C(S_i0 ell_i + lambda K_i ell_i²/4),
 ell_i=(D_i-U_i)/(2D_i), U_i=u_i+gamma N_j omega_i R0,
 D_i=Delta_i + eta_i lambda, eta_i=gamma N_j omega_i Pi/2.

 (16 D_i^3/C)*dT_i/dlambda =
   omega_i Pi U_i [4 D_i u_i + 2 Delta_i (D_i-U_i)]
   + 8 D_i U_i eta_i S_i0
   + K_i (D_i-U_i)[D_i(D_i+U_i)-2 Delta_i U_i].
All identities are algebraic. Interval arithmetic uses exact Fraction endpoints.
No grids of sampled equilibria are used as proofs. Each interval cell covers
EVERY point of its rectangle.

Weak gamma (at fixed rho, u_i constant), uniformly over compact lambda:
 dT_i/dlambda = gamma*N_i*N_j/(16 Delta_i²) *
  [ M_i(Delta_i-u_i)^2 + omega_i Pi(Delta_i²+3u_i²) ] + O(gamma²).
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

class I:
    def __init__(self,lo,hi=None):
        if isinstance(lo,I):
            self.lo,self.hi=lo.lo,lo.hi
            return
        self.lo=F(lo)
        self.hi=F(lo if hi is None else hi)
        assert self.lo<=self.hi
    def __add__(self,z):
        z=I(z);return I(self.lo+z.lo,self.hi+z.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,z):return self+-I(z)
    def __rsub__(self,z):return I(z)+-self
    def __mul__(self,z):
        z=I(z);v=[self.lo*z.lo,self.lo*z.hi,self.hi*z.lo,self.hi*z.hi]
        return I(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,z):
        z=I(z)
        assert z.lo>0 or z.hi<0
        return self*I(1/z.hi,1/z.lo)
    def __rtruediv__(self,z):return I(z)/self
    def approx(self):return [float(self.lo),float(self.hi)]

def cournot():
    n=(5,2)
    chi=F(9,10)
    m=[]
    pi=F(0)
    for ni in n:
        den=2*(ni+1)-ni*chi*chi
        x=(2-chi)/den
        y=F(ni+1-ni*chi,1)/den
        x0=F(1,ni+1)
        cs1=((ni*x)**2+2*chi*ni*x*y+y*y)/2
        cs0=(ni*x0)**2/2
        m.append(cs1+ni*x*x-cs0-ni*x0*x0)
        pi+=y*y
    return m,pi

def test():
    m,pi=cournot()
    assert m[0]<0 and m[1]<0 and pi>0
    assert m[0]+m[1]+pi>0  # national incremental product-market surplus
    gamma=F(1,10)
    rho=F(2,5)
    u=F(7,25) # v_B+rho w_B = 0.10+0.4*0.45
    delta=F(1)
    n_i=F(1)
    n_j=F(3)
    zeta=F(5,1000)
    k=F(4,10000)
    p0=F(1,20)
    pb=p0+zeta*zeta*pi/k
    r0=p0*pi+zeta*zeta*pi*pi/(2*k)
    assert pb+F(1,2)<1
    wcrit=(-m[1])*(delta-u)**2/(pi*(delta*delta+3*u*u))
    w_upper=(-m[1])/pi
    assert F(13,100)>wcrit and F(27,100)<w_upper
    assert wcrit<w_upper and wcrit>0
    # A genuine nondegenerate region: for ALL omegaB and lambda in
    # [0.13,0.27] x [0,0.5], derivative positive while S0,K_B negative.
    wlo,whi=F(13,100),F(27,100)
    nw,nl=10,10
    minP=None; maxK=None; maxS0=None
    minell=None; maxell=None; mindT=None
    result=[]
    for idx_w,idx_l in product(range(nw),range(nl)):
        w=I(wlo+(whi-wlo)*F(idx_w,nw),
            wlo+(whi-wlo)*F(idx_w+1,nw))
        lam=I(F(idx_l,2*nl),F(idx_l+1,2*nl))
        eta=gamma*n_j*w*pi/2
        U=u+gamma*n_j*w*r0
        D=delta+eta*lam
        Si0=pb*m[1]+w*r0
        Ki=m[1]+w*pi
        ell=(D-U)/(2*D)
        P=w*pi*U*(4*D*u+2*delta*(D-U))+8*D*U*eta*Si0+Ki*(D-U)*(D*(D+U)-2*delta*U)
        dT=gamma*n_i*n_j*P/(16*D*D*D)
        assert P.lo>0,(idx_w,idx_l,'polynomial sign',P.approx())
        assert dT.lo>0,(idx_w,idx_l,'derivative',dT.approx())
        assert Ki.hi<0 and Si0.hi<0,(idx_w,idx_l,'innovation harms local region',Ki.approx(),Si0.approx())
        assert ell.lo>0 and ell.hi<F(1,2),(idx_w,idx_l,'interior participation',ell.approx())
        minP=P.lo if minP is None else min(minP,P.lo)
        maxK=Ki.hi if maxK is None else max(maxK,Ki.hi)
        maxS0=Si0.hi if maxS0 is None else max(maxS0,Si0.hi)
        mindT=dT.lo if mindT is None else min(mindT,dT.lo)
        minell=ell.lo if minell is None else min(minell,ell.lo)
        maxell=ell.hi if maxell is None else max(maxell,ell.hi)
    # At the ORIGINAL PROFIT INCIDENCE wB=0.2, finite matching gamma 0.1,
    # lambda in the whole [0,0.5] interval, not only at individual grids.
    def point(lam):
        w=F(1,5)
        eta=gamma*n_j*w*pi/2
        U=u+gamma*n_j*w*r0
        D=delta+eta*lam
        ell=(D-U)/(2*D)
        C=gamma*n_i*n_j
        Ki=m[1]+w*pi
        Si0=pb*m[1]+w*r0
        ellprime=U*eta/(2*D*D)
        T=n_i*(u*ell-delta/8+delta*ell*ell/2)+C*(Si0*ell+lam*Ki*ell*ell/4)
        dT=(n_i*(u+delta*ell)+C*Si0)*ellprime+C*Ki*(ell*ell+2*lam*ell*ellprime)/4
        P=w*pi*U*(4*D*u+2*delta*(D-U))+8*D*U*eta*Si0+Ki*(D-U)*(D*(D+U)-2*delta*U)
        assert 16*D**3*dT==C*P
        return {'lambda':str(lam),'public_entry_threshold':float(T),
                'marginal_lambda_effect':float(dT),
                'cutoff':float(ell)}
    vals=[point(F(0)),point(F(1,4)),point(F(1,2))]
    # Weak-coupling criterion: nonempty omega range has local K negative
    # AND the marginal value of access strictly rises with complementarity.
    lower=F(13,100)
    upper=F(27,100)
    assert wcrit<lower and upper<w_upper
    info={
      'scope':'Frozen model, region B P-to-B conditional on region A already full G/B. NOT a complete four-player equilibrium theorem or an all-parameter global pricing certificate.',
      'joint_parameter_rectangle':{'omega_B':[str(wlo),str(whi)], 'lambda':['0','1/2'],
                 'omega_A':'1 - omega_B', 'gamma':'1/10', 'partner_status':'G or B (full participation)', 'exact_fraction_cells':nw*nl},
      'proof':'Exact rational inclusion arithmetic for entire 2D rectangle, no interpolation of point values; 16D^3*T_lambda=C*P identity (analytical).',
      'uniform_certificates':{'min_derivative_numerator':float(minP),
            'min_dT_dlambda':float(mindT),
            'max_local_market_K':float(maxK),
            'max_baseline_market_S0':float(maxS0),
            'cutoff_interior_enclosure':[float(minell),float(maxell)]},
      'weak_gamma_general_sign': 'M_i*(Delta_i-u_i)^2+omega_i*Pi*(Delta_i^2+3*u_i^2)',
      'weak_gamma_critical_omega':float(wcrit),
      'weak_gamma_local_market_negative_omega_upper':float(w_upper),
      'national_marginal_venture_incidence_positive':float(m[0]+m[1]+pi),
      'M_B':float(m[1]),'Pi':float(pi),'baseline':vals
    }
    out=Path(__file__).with_name('incidence_screening_sign_results.json')
    out.write_text(json.dumps(info,ensure_ascii=False,indent=2)+'\n')
    print('PASS exactly certified full rectangle omegaB in [0.13,0.27] and lambda in [0,0.5] with local K<0,S0<0 and dT/dlambda>0')
    print('WEAK GAMMA omega share interval',float(wcrit),float(w_upper))
    print('CONTINUOUS 2D PROOF min rational derivative lower:',float(mindT),'max K:',float(maxK),'max S0:',float(maxS0))
    print('BASELINE',vals)

if __name__=='__main__':test()
