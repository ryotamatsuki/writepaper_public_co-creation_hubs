from __future__ import annotations
import sympy as sp


def br_slope(Fxx,Fxy):
    return -Fxy/Fxx


def verify_sign_logic():
    a,b=sp.symbols('a b', positive=True)
    assert sp.simplify(br_slope(-a,b)-b/a)==0
    assert sp.simplify(br_slope(-a,-b)+b/a)==0
    return True


if __name__=='__main__':
    verify_sign_logic()
    print('best-response sign logic under SOC: PASS')
