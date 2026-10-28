#!/usr/bin/env python3
"""Regression for the authorised p_b>0 hypothesis correction, not a new theorem.
Run with Python 3 and SymPy 1.14+. No frozen scientific file is modified.
"""
from pathlib import Path
from collections import Counter
import subprocess,re,json,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parents[3]
BASE='3d2eaf2'
PROD='hotelling_entry_experimentation/production/'
def original(path):return subprocess.check_output(['git','show',BASE+':'+path],cwd=ROOT).decode()
def flatten(path,reader):
    return re.sub(r'\\input\{([^}]+)\}',lambda m:flatten(PROD+m[1]+('' if m[1].endswith('.tex') else '.tex'),reader),reader(path))
def norm(text):return re.sub(r'\s+',' ',text).strip()
def displays(text):return Counter(norm(m[0]) for m in re.finditer(r'\\begin\{(equation\*?|align\*?|gather\*?|multline\*?)\}.*?\\end\{\1\}|\\\[.*?\\\]',text,re.S))
def labels(text):return Counter(re.findall(r'\\label\{([^}]+)\}',text))
# Exact algebra in the diagnosed proof bridge.
p0,zeta,w=s.symbols('p0 zeta omega',nonnegative=True)
pi,k,loss,pb=s.symbols('Pi k loss pb',positive=True)
M=s.symbols('M',real=True)
p_b=p0+zeta**2*pi/k
R0=p0*pi+zeta**2*pi**2/(2*k)
K=M+w*pi
S0=p_b*M+w*R0
identity=s.simplify(S0-(p_b*K-w*zeta**2*pi**2/(2*k)))
assert identity==0
zero_boundary=s.simplify(S0.subs({p0:0,zeta:0}))
assert zero_boundary==0  # The former strict conclusion fails at this old-domain edge.
positive_baseline_form=-pb*loss-w*zeta**2*pi**2/(2*k)
assert positive_baseline_form.is_negative is True
# Both ways to obtain a positive baseline under maintained nonnegative primitives.
p0_positive=s.symbols('p0_positive',positive=True)
zeta_positive=s.symbols('zeta_positive',positive=True)
assert s.expand(S0.subs(M,-loss-w*pi)).subs(p0,p0_positive).is_negative is True
assert s.expand(S0.subs(M,-loss-w*pi)).subs(zeta,zeta_positive).is_negative is True
report={'baseline':BASE,'sympy_version':s.__version__,'symbolic_regression':{'S0_identity_residual':str(identity),'zero_baseline_S0':str(zero_boundary),'positive_baseline_K_negative_implies_S0_negative':True,'p0_positive_branch':True,'zeta_positive_branch':True},'active_targets':{},'protected_changed':[]}
for target in ('main','supplement'):
    path=PROD+'manuscript/'+target+'.tex'
    old=flatten(path,original);new=flatten(path,lambda p:(ROOT/p).read_text())
    a,b=displays(old),displays(new)
    assert a==b,(target,a-b,b-a)
    assert labels(old)==labels(new),target
    oldtables=Counter(norm(m[0]) for m in re.finditer(r'\\begin\{tabular\}.*?\\end\{tabular\}',old,re.S))
    newtables=Counter(norm(m[0]) for m in re.finditer(r'\\begin\{tabular\}.*?\\end\{tabular\}',new,re.S))
    assert oldtables==newtables,target
    report['active_targets'][target]={'display_instances':sum(a.values()),'display_multiset_equal':True,'label_multiset_equal':True,'tabular_multiset_equal':True}
# The formal statement changes only through the authorised hypothesis insertion.
p=PROD+'manuscript/sections/04_propositions.tex'
def prop(text):return text.split('\\paragraph{Proposition 1')[1].split('The critical share')[0]
before=prop(original(p));after=prop((ROOT/p).read_text())
oldphrase='and nonnegative $R_0$.'
newphrase='nonnegative $R_0$, and a positive baseline success intercept $p_b>0$ as defined in \\eqref{eq:effort}.'
assert before.replace(oldphrase,newphrase)==after
report['proposition_statement_change']='EXACTLY authorised p_b>0 hypothesis insertion; all other statement text unchanged'
tracked=subprocess.check_output(['git','ls-tree','-r','--name-only',BASE],cwd=ROOT,text=True).splitlines()
protected=[p for p in tracked if (p.startswith('formal/') or p.startswith(PROD)) and Path(p).suffix in ('.py','.lean','.csv','.json')]
archived=[p for p in tracked if p.startswith(PROD+'revision_kawasaki_results/') or p.startswith(PROD+'reviews/KAWASAKI_')]
for p in set(protected+archived):
    now=ROOT/p
    if not now.exists() or subprocess.check_output(['git','show',BASE+':'+p],cwd=ROOT)!=now.read_bytes():report['protected_changed'].append(p)
assert not report['protected_changed'],report['protected_changed']
report['protected_scientific_files_compared']=len(protected)
report['archived_PR73_files_compared']=len(archived)
report['lean']='NOT EXECUTED BY THIS CHECK; existing Stage33Formal does not encode Proposition 1; source byte-identical'
out=Path(__file__).with_suffix('.json')
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
