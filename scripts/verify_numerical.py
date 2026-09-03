from __future__ import annotations
import re
from generate_results import build_results, ROOT
def close(a,b,t): assert abs(a-b)<=t, f'{a} != {b} within {t}'
def welfare_freeze():
    text=(ROOT/'public_two_sided_platform_stage8_theory_freeze/WELFARE_REGISTER.md').read_text(encoding='utf-8')
    def find(label,field):
        line=next(x for x in text.splitlines() if x.startswith(f'- {label}:')); m=re.search(rf'{re.escape(field)}≈([+.0-9-]+)',line); assert m,(label,field); return float(m.group(1))
    return {'G_Wi':find('G','W_i'),'G_Pi':find('G','Pi_T'),'G_WN':find('G','W^N'),'B3_Wi':find('B3','W_i'),'B3_Pi':find('B3','Pi_T'),'B3_WN':find('B3','W^N'),'diff':find('G−B3','W^N')}
r=build_results(); c=r['computed']; f=r['frozen_display']; wf=welfare_freeze(); close(c['G']['x'],f['G_x'],5e-5); close(c['G']['p_T'],f['G_p'],5e-6); close(c['G']['BR_slope'],f['G_br'],5e-4); assert c['G']['BR_slope']<0<c['B3']['BR_slope']; close(c['B3']['BR_slope'],f['B3_br'],5e-4)
for x in (c['B3']['x_stage4_script'],c['B3_welfare_script_x']): assert min(abs(x-f['B3_strategic_x']),abs(x-f['B3_welfare_x']))<=5e-5
close(c['private_price_response'],f['p_x'],6e-5); w=c['welfare']
for got,exp in [(w['G']['W_i'],wf['G_Wi']),(w['G']['Pi_T'],wf['G_Pi']),(w['G']['W_N'],wf['G_WN']),(w['B3']['W_i'],wf['B3_Wi']),(w['B3']['Pi_T'],wf['B3_Pi']),(w['B3']['W_N'],wf['B3_WN']),(c['welfare_difference'],wf['diff'])]: close(got,exp,6e-5)
assert abs(c['coordination']['national_wedge']-.42277)<=6e-5; assert c['local_perturbation_successes']==20; print('PASS: numerical witness, welfare, price response, signs, and 20/20 local certificate')
