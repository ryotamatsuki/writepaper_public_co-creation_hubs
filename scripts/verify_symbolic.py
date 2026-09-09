from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
checks = [
    ROOT / "analytic_level3" / "code" / "verify_symbolic_identities.py",
    ROOT / "public_two_sided_platform_welfare_generality" / "scripts" / "verify_welfare_identities.py",
]

for script in checks:
    assert script.exists(), f"missing symbolic verification source: {script.relative_to(ROOT)}"
    subprocess.run([sys.executable, str(script)], cwd=str(script.parent), check=True)

print("PASS: v2.1 small-beta theorem identities and fee-transfer identity")
