from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
files=[ROOT/'generated/results/canonical_results.json',*sorted((ROOT/'generated/tables').glob('*.tex')),ROOT/'generated/figures/README.md']
data={'_generated':'DO NOT EDIT — GENERATED FILE','files':{}}
for p in files:
 b=p.read_bytes(); data['files'][str(p.relative_to(ROOT))]={'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)}
out=ROOT/'generated/results/manifest.json'; out.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8'); print(out)
