from __future__ import annotations
import subprocess,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
SCRIPTS=['derive_participation.py','derive_private_price.py','derive_b3_br.py','derive_full_game_br.py','derive_small_beta.py']
for name in SCRIPTS:
    r=subprocess.run([sys.executable,str(HERE/name)],capture_output=True,text=True,timeout=45)
    if r.returncode:
        raise SystemExit(f'{name} failed:\n{r.stdout}\n{r.stderr}')
    print(r.stdout.strip())
print('analytic symbolic identity suite: PASS')
