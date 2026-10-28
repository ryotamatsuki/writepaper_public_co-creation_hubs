#!/usr/bin/env python3
"""Ensure the frozen-model audit is a reproducible, gap-visible coverage ledger."""
import csv,collections
from pathlib import Path

ROOT=Path(__file__).resolve().parent
p=ROOT/'EQUILIBRIUM_COVERAGE_MATRIX.csv'
with p.open(encoding='utf-8',newline='') as f:
 rows=list(csv.DictReader(f))
assert len(rows)>=80,('unexpectedly incomplete',len(rows))
ids=[r['ID'] for r in rows]
assert len(ids)==len(set(ids)),'duplicate audit identifier'
required={'S','C','P','E','W','R','L'}
counts=collections.Counter(r['ID'][0] for r in rows)
assert required <= set(counts),('missing category',required-set(counts))
for r in rows:
 assert all(str(v).strip() for v in r.values()),('missing description / status / evidence / next step',r)
 assert r['ID'][0] in required,('unknown audit category',r['ID'])
 assert '完了' not in r['監査ステータス'] or r['現時点の検証資産']!='該当なし'
 # Prevent conflating original finite point checks with global general results.
 if r['ID'] in {'P05','P06','E08','R03','R04'}:
  assert ('基準' in r['監査ステータス'] or '連続' in r['監査ステータス']
          or '各点' in r['監査ステータス'] or '矩形' in r['監査ステータス'])
for group,n in sorted(counts.items()):
 print(group,n)
print('PASS coverage ledger',len(rows),'unique tracked questions, 7 domains, all gaps explicitly labelled')
