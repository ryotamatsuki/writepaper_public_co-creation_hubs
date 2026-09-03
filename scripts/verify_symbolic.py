from __future__ import annotations
import subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; required=['public_two_sided_platform_hard_kill/code/derive_project_cutoffs.py','public_two_sided_platform_hard_kill/code/derive_private_price.py','public_two_sided_platform_hard_kill/code/verify_identities.py','public_two_sided_platform_welfare_generality/scripts/verify_welfare_identities.py']
for rel in required: subprocess.run([sys.executable,str(ROOT/rel)],cwd=str((ROOT/rel).parent),check=True)
print('PASS: frozen symbolic/identity source scripts')
