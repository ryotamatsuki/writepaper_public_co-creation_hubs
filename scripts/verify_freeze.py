from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'public_two_sided_platform_stage8_theory_freeze'
checks={
    'FREEZE_MANIFEST.md':['53405fa904dd31817639a734de2063158ec69321','BR_i^{B3′}>0>BR_i^{G′}','NOT SELECTED'],
    'PROOF_STATUS_TABLE.md':['NUMERICALLY SUPPORTED ONLY'],
    'KILLED_CLAIMS_LOCKFILE.md':['P1-R','P3–P5','global G welfare dominance over B3','Status: **LOCKED**'],
    'JOURNAL_POSITIONING_DEFERRED.md':['NO PRIMARY JOURNAL IS SELECTED AT STAGE 8','Stage 12 — Journal Positioning','submission fee = 0','mandatory publication/APC = 0'],
}
for fn,needles in checks.items():
    text=(D/fn).read_text(encoding='utf-8')
    for n in needles:
        assert n in text,(fn,n)
print('PASS: Stage-8 freeze identifiers, claim ceiling, killed claims, journal deferral')
