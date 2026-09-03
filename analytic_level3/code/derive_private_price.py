from __future__ import annotations
import sympy as sp


def price_policy_blocks():
    Rp,Ri,Rj,Rij,Rip,Rjp,Rpp = sp.symbols(
        'R_p R_i R_j R_ij R_ip R_jp R_pp', nonzero=True
    )
    pi=-Ri/Rp
    pj=-Rj/Rp
    pij=-(Rij+Rip*pj+Rjp*pi+Rpp*pi*pj)/Rp
    return {'p_i':sp.factor(pi),'p_j':sp.factor(pj),'p_ij':sp.factor(pij)}


def verify():
    b=price_policy_blocks()
    Rp,Ri,Rj,Rij,Rip,Rjp,Rpp=sp.symbols(
        'R_p R_i R_j R_ij R_ip R_jp R_pp', nonzero=True
    )
    residual=sp.expand(Rij+Rip*b['p_j']+Rjp*b['p_i']+Rpp*b['p_i']*b['p_j']+Rp*b['p_ij'])
    assert sp.simplify(residual)==0
    return True


if __name__=='__main__':
    verify()
    print('private-price IFT blocks: PASS')
