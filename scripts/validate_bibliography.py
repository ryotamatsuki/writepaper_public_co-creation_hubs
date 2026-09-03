from pathlib import Path
import re
p=Path(__file__).resolve().parents[1]/'references/references.bib'; t=p.read_text(encoding='utf-8'); keys=re.findall(r'@\w+\{([^,]+),',t); assert keys and len(keys)==len(set(keys)),'duplicate/empty keys'
for entry in re.split(r'\n@',t):
    if not entry.strip(): continue
    e='@'+entry if not entry.startswith('@') else entry
    for field in ('author','title','year'): assert re.search(rf'\b{field}\s*=',e,re.I),f'missing {field}'
print('PASS: bibliography parse-level validation',len(keys),'entries')
