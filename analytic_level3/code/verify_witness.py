from __future__ import annotations
import importlib.util,math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'public_two_sided_platform_hard_kill/code/numerical_hard_kill.py'


def load_stage4():
    spec=importlib.util.spec_from_file_location('stage4_l3',SRC)
    mod=importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


def verify():
    s4=load_stage4(); P=s4.BASE
    xg,dg,brg=s4.stationary(P,'G',None,(.62,.75))
    pg=s4.pstar(xg,xg,P)
    xb,db,brb=s4.stationary(P,'B3',pg,(.60,.72))
    e=s4.eq(xg,xg,pg,P)
    assert e is not None
    assert dg[1]<-0.1 and dg[2]<-0.001 and brg<-0.01
    assert db[1]<-0.05 and db[2]>0.005 and brb>0.05

    h=1e-5
    px=(s4.pstar(xg+h,xg,P)-s4.pstar(xg-h,xg,P))/(2*h)
    assert px<-0.01 and abs(px-(-0.016505))<8e-4

    delta=P['alpha']*P['beta']; a=P['kT']+pg
    AT=P['v']+P['alpha']*P['rhoT']
    B=AT+delta*(e['t1']+e['t2'])
    disc=B*B-8*delta*a
    assert disc>0
    upper=(B+math.sqrt(disc))/2
    lower=(B-math.sqrt(disc))/2
    assert abs(upper-e['qT'])<5e-8
    assert a/lower>1.0

    t=e['t1']; qT=e['qT']; qi=e['q1']
    u=qi-qT-delta*t
    w=1-2*delta*a/(qT*qT)
    H=u*w-2*delta*t
    det=u*H
    assert abs(det)>1e-2

    print(f'G x={xg:.12f} p={pg:.12f} own={dg[1]:.12f} cross={dg[2]:.12f} BR={brg:.12f}')
    print(f'B3 x={xb:.12f} own={db[1]:.12f} cross={db[2]:.12f} BR={brb:.12f}')
    print(f'p_x={px:.12f} disc={disc:.12f} detJ={det:.12f}')
    print('canonical witness reconciliation: PASS')


if __name__=='__main__':
    verify()
