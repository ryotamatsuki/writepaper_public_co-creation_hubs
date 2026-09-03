from pathlib import Path
import re
log=(Path(__file__).resolve().parents[1]/"paper"/"main.log").read_text(encoding="utf-8",errors="replace")
for pat in ("There were undefined citations","There were undefined references"):
    assert pat not in log, pat
assert not re.search(r"(Citation|Reference) `[^']+' .* undefined", log), "undefined citation/reference remains"
print("PASS: final LaTeX log has no unresolved citations/references")
