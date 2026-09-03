from __future__ import annotations
import sympy as sp


def correction_blocks():
    M,Wxip,Wxjp,Wpp,Wp,pi,pj,pij=sp.symbols(
        'M W_xip W_xjp W_pp W_p p_i p_j p_ij', finite=True
    )
    Pprice=Wxip*pj+Wxjp*pi+Wpp*pi*pj+Wp*pij
    return M,Pprice,sp.expand(M+Pprice)


def verify_polynomial_chain_rule():
    x,y=sp.symbols('x y')
    c=sp.symbols('c0:10')
    p=c[0]+c[1]*x+c[2]*y+c[3]*x*y+c[4]*x**2+c[5]*y**2
    z=sp.symbols('z')
    W=c[6]*x*y+c[7]*x*z+c[8]*y*z+c[9]*z**2+x**2*y*z
    reduced=sp.expand(W.subs(z,p))
    lhs=sp.diff(reduced,x,y)
    Wxy=sp.diff(W,x,y).subs(z,p)
    Wxp=sp.diff(W,x,z).subs(z,p)
    Wyp=sp.diff(W,y,z).subs(z,p)
    Wpp=sp.diff(W,z,z).subs(z,p)
    Wp=sp.diff(W,z).subs(z,p)
    rhs=Wxy+Wxp*sp.diff(p,y)+Wyp*sp.diff(p,x)+Wpp*sp.diff(p,x)*sp.diff(p,y)+Wp*sp.diff(p,x,y)
    assert sp.simplify(lhs-rhs)==0
    return True


if __name__=='__main__':
    verify_polynomial_chain_rule()
    print('full-game chain-rule decomposition: PASS')
