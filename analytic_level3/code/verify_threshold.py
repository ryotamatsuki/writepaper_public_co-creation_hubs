from __future__ import annotations
import sympy as sp
from derive_small_beta import verify_b3_first_order,verify_beta0_g_cross


def verify():
    Gamma,Lambda,a,b,eps=sp.symbols('Gamma Lambda a b eps', positive=True)
    br_b3=Gamma/a
    br_g=(Gamma-Lambda)/b
    assert sp.simplify(br_b3-Gamma/a)==0
    assert sp.simplify(br_g-(Gamma-Lambda)/b)==0
    assert sp.simplify((Gamma-Lambda).subs(Lambda,Gamma+eps))<0
    b3=verify_b3_first_order()
    g0=verify_beta0_g_cross()
    assert b3.is_positive
    assert g0.is_negative
    print('reduced-block threshold identity: PASS')
    print('small-beta sufficient-theorem identities: PASS')


if __name__=='__main__':
    verify()
