from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
F=ROOT/'theory_freeze_v21'/'CANONICAL_THEORY_FREEZE_2026-09-10.md'
text=F.read_text(encoding='utf-8')
needles=[
    'THEORY FROZEN — GO TO REPRODUCIBILITY SETUP',
    'beta=.01',
    'gamma=.825',
    'tau=.35',
    '0.8371022382025995',
    '0.8258903860237495',
    'LOCAL SUFFICIENT-CONDITION THEOREM',
    'all-regime computational global-equilibrium existence witness',
    'old vector `(beta=.05, gamma=.9, tau=.05)` is rejected',
    'reviews/STAGE_075A_V21_GENERALITY_QUANTIFIER_RED_TEAM_2026-09-10.md',
]
for n in needles:
    assert n in text, n
print('PASS: v2.1 repaired theory freeze identity and claim ceiling')
