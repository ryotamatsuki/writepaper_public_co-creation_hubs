"""Search the existing public baseline-partner primitive rho.

The Stage-4A counterexample exploits positive partner surplus at a nearly inactive
own public hub. Lower rho directly reduces that outside/free-riding payoff without
changing the model. We keep gamma=.9 and set tau=.35 so the rival public route is
globally dominated for nonresidents; then search rho for an interior local
B3-complements/G-substitutes pair that is also a global public best response.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from public_two_sided_platform_hard_kill.code.numerical_hard_kill import BASE as OLD_BASE
from stage4_v21_global.code.gamma_global_region_search import central_roots
from public_two_sided_platform_hard_kill.code.numerical_hard_kill import pstar
from stage4_v21_global.code.all_regime_global_search import (
    SOLVED_EQUILIBRIUM, private_best_response, _global_public_br
)

TAU=.35
GAMMA=.9


def old_par(rho):
    p=OLD_BASE.copy(); p['tau']=TAU; p['gamma']=GAMMA; p['rho']=float(rho); return p


def new_par(rho):
    o=old_par(rho)
    return dict(v=o['v'],alpha=o['alpha'],beta=o['beta'],rho=o['rho'],rhoT=o['rhoT'],
                kL=o['k'],kT=o['kT'],tau=o['tau'],gamma=o['gamma'])


def local_pairs(rho):
    P=old_par(rho); out=[]
    for xg,dg,sg in central_roots(P,'G',None):
        try: pg=float(pstar(xg,xg,P))
        except Exception: continue
        for xb,db,sb in central_roots(P,'B3',pg):
            if sg<0<sb:
                out.append(dict(xg=xg,dg=dg,sg=sg,pg=pg,xb=xb,db=db,sb=sb))
    return out


def global_check(rho,pair):
    par=new_par(rho); xg=pair['xg']; xb=pair['xb']
    pg=private_best_response(xg,xg,par,grid_n=61,rigorous=True)
    if pg.status!=SOLVED_EQUILIBRIUM: raise RuntimeError(pg.status)
    g=_global_public_br('G',xg,xg,par,x_grid_n=41,price_grid=25)
    b=_global_public_br('B3',xb,xb,par,pbar=pg.p,x_grid_n=61,price_grid=25)
    return pg,g,b


def main():
    assert new_par(.15)['kL']+TAU > new_par(.15)['v']+new_par(.15)['alpha']
    passing=[]
    for rho in (0.0,.01,.02,.03,.04,.05,.06,.075,.09,.105,.12,.135,.15,.18,.21):
        pairs=local_pairs(rho)
        print('RHO_LOCAL',rho,'pairs',len(pairs))
        for pair in pairs:
            try: pg,g,b=global_check(rho,pair)
            except Exception as exc:
                print('RHO_UNRESOLVED',rho,repr(exc)); continue
            print('RHO_SCREEN',rho,'xG',pair['xg'],'sG',pair['sg'],
                  'xB',pair['xb'],'sB',pair['sb'],'pG',pg.p,
                  'G_BR',g.x,'G_GAIN',g.gain_over_candidate,
                  'B3_BR',b.x,'B3_GAIN',b.gain_over_candidate)
            if (g.gain_over_candidate<=5e-4 and b.gain_over_candidate<=5e-4
                and abs(g.x-pair['xg'])<=8e-3 and abs(b.x-pair['xb'])<=8e-3):
                passing.append((rho,pair,pg,g,b))
    if not passing:
        print('RHO_GLOBAL_REGION_SEARCH: NO CONSTRUCTION PASS'); raise SystemExit(2)
    rho,pair,pg,g,b=passing[0]
    print('SELECTED_RHO',rho); print('SELECTED_PAIR',pair); print('SELECTED_PG',pg)
    print('SELECTED_GLOBAL_G',g); print('SELECTED_GLOBAL_B3',b)
    print('RHO_GLOBAL_REGION_SEARCH: CONSTRUCTION PASS')

if __name__=='__main__': main()
