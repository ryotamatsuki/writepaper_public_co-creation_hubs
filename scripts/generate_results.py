from __future__ import annotations
import json, re
import numpy as np
from source_loader import ROOT, load
S4_PATH='public_two_sided_platform_hard_kill/code/numerical_hard_kill.py'
W7_PATH='public_two_sided_platform_welfare_generality/scripts/welfare_b3_vs_g.py'
WITNESS=ROOT/'public_two_sided_platform_stage8_theory_freeze/NUMERICAL_WITNESS_REGISTER.md'
def _frozen_display():
    t=WITNESS.read_text(encoding='utf-8')
    def one(pattern):
        m=re.search(pattern,t)
        if not m: raise RuntimeError(f'missing freeze value: {pattern}')
        return float(m.group(1))
    return {'G_x':one(r'x_1=x_2=([.0-9-]+)'),'G_p':one(r'p_T=([.0-9-]+)'),'G_br':one(r"BR'_G=([.0-9-]+)"),'B3_strategic_x':one(r"Stage-4 strategic solver reports `x_B3=([.0-9-]+)"),'B3_welfare_x':one(r"Stage-7 welfare recomputation reports `x_B3=([.0-9-]+)"),'B3_br':one(r"BR'_B3=([.0-9-]+)"),'p_x':one(r"p\^\*_{x_i}=([.0-9-]+)")}
def build_results():
    s4=load('stage4_numeric',S4_PATH); w7=load('stage7_welfare',W7_PATH)
    xg,dg,brg=s4.stationary(s4.BASE,'G',None,(.62,.75)); pg=s4.pstar(xg,xg,s4.BASE)
    xb_s4,db,brb=s4.stationary(s4.BASE,'B3',pg,(.60,.72)); xg7,pg7,xb_w7=w7.solve_symmetric(); welfare={}
    for lab,x in [('G',xg7),('B3',xb_w7)]:
        e=w7.eq(x,x,pg7,w7.BASE); wi=w7.W(1,x,x,pg7,w7.BASE); pi=w7.profit(pg7,x,x,w7.BASE)
        welfare[lab]={'x':float(x),'p_T':float(pg7),'n_T':float(e['nT']),'W_i':float(wi),'Pi_T':float(pi),'W_N':float(2*wi+pi)}
    h=1e-5; px=(s4.pstar(xg+h,xg,s4.BASE)-s4.pstar(xg-h,xg,s4.BASE))/(2*h)
    def d1(f,a,b): return (f(a+h,b)-f(a-h,b))/(2*h)
    fw2=lambda a,b:w7.W(2,a,b,w7.pstar(a,b,w7.BASE),w7.BASE); fpi=lambda a,b:w7.profit(w7.pstar(a,b,w7.BASE),a,b,w7.BASE)
    dw=d1(fw2,xg7,xg7); dpi=d1(fpi,xg7,xg7); rng=np.random.default_rng(20260903); ok=0
    for _ in range(20):
        P=s4.BASE.copy()
        for key in ['v','alpha','beta','rho','rhoT','k','kT','gamma']: P[key]*=1+rng.uniform(-.005,.005)
        try:
            x1,d1s,s1=s4.stationary(P,'G',None,(.62,.75)); p1=s4.pstar(x1,x1,P); x2,d2s,s2=s4.stationary(P,'B3',p1,(.60,.72)); ok+=int(d1s[1]<0 and d2s[1]<0 and s1<0<s2)
        except Exception: pass
    return {'_generated':'DO NOT EDIT — GENERATED FILE','sources':[S4_PATH,W7_PATH,'public_two_sided_platform_stage8_theory_freeze/NUMERICAL_WITNESS_REGISTER.md'],'parameters':{k:float(v) for k,v in s4.BASE.items()},'frozen_display':_frozen_display(),'computed':{'G':{'x':float(xg),'p_T':float(pg),'own_second':float(dg[1]),'cross_second':float(dg[2]),'BR_slope':float(brg)},'B3':{'x_stage4_script':float(xb_s4),'own_second':float(db[1]),'cross_second':float(db[2]),'BR_slope':float(brb)},'B3_welfare_script_x':float(xb_w7),'private_price_response':float(px),'welfare':welfare,'welfare_difference':float(welfare['G']['W_N']-welfare['B3']['W_N']),'coordination':{'rival_welfare':float(dw),'private_profit':float(dpi),'national_wedge':float(dw+dpi)},'local_perturbation_successes':ok,'local_perturbation_total':20},'proof_status':{'headline':'NUMERICALLY SUPPORTED ONLY','robustness':'LOCAL/CONDITIONAL'},'journal_target':'NOT SELECTED'}
def main():
    out=ROOT/'generated/results/canonical_results.json'; out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(build_results(),indent=2,sort_keys=True)+'\n',encoding='utf-8'); print(out)
if __name__=='__main__': main()
