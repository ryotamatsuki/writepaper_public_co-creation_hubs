from pathlib import Path
import json, re, zipfile

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "paper/main.tex"
SUB = ROOT / "submission/jpet"
TITLE_PAGE = ROOT / "Title_Page.docx"
OUT = ROOT / "generated/stage14/qa_summary.json"
TITLE = "Strategic Interaction among Public Innovation Hubs under Private Repricing"

main = MAIN.read_text(encoding="utf-8")
for label in ["Introduction.", "Methods.", "Results.", "Conclusion."]:
    assert f"\\textbf{{{label}}}" in main
kw = re.search(r"\\textbf\{Keywords:\}\s*(.+?)\\\\", main, re.S)
assert kw
keywords = [x.strip() for x in kw.group(1).split(";") if x.strip()]
assert 1 <= len(keywords) <= 7
jel = re.search(r"\\textbf\{JEL Classification:\}\s*([^\n]+)", main)
assert jel
jel_codes = re.findall(r"\b[A-Z][0-9]{2}\b", jel.group(1))
assert len(jel_codes) >= 2
assert "Data Availability Statement" in main
assert "OpenAI ChatGPT was used" in main

text = "\n".join([main] + [p.read_text(encoding="utf-8") for p in sorted((ROOT / "sections").glob("*.tex"))])
for pattern in [r"\bTODO\b", r"\bFIXME\b", r"\bTK\b", r"citation needed", r"\bplaceholder\b", r"\bStage 13\b", r"\bStage 14\b", r"\bL3-1\b", r"\breviewer attack\b"]:
    assert not re.search(pattern, text, re.I), pattern

required = ["README.md", "metadata.md", "cover_letter.md", "declarations.md", "source_inventory.md", "reviewer_candidates.md"]
for name in required:
    assert (SUB / name).exists(), name

with zipfile.ZipFile(TITLE_PAGE) as z:
    xml = z.read("word/document.xml").decode("utf-8")
plain = re.sub(r"<[^>]+>", " ", xml)
flags = {
    "current_title_present": TITLE.lower() in plain.lower(),
    "email_present": bool(re.search(r"[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}", plain, re.I)),
    "orcid_present": bool(re.search(r"\b\d{4}-\d{4}-\d{4}-[\dX]{4}\b", plain, re.I)),
    "funding_language_present": bool(re.search(r"funding|financial support|grant|no external funding", plain, re.I)),
    "conflict_language_present": bool(re.search(r"conflict of interest|competing interest|no competing", plain, re.I)),
    "data_availability_language_present": bool(re.search(r"data availability", plain, re.I)),
}
assert flags["current_title_present"]
assert flags["email_present"]

warnings = []
if not flags["orcid_present"]:
    warnings.append("ORCID must be entered in Wiley Authors/Research Exchange before submission freeze")
if not flags["funding_language_present"]:
    warnings.append("Funding declaration needs factual author confirmation before submission freeze")
if not flags["conflict_language_present"]:
    warnings.append("Conflict-of-interest declaration needs factual author confirmation before submission freeze")

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps({
    "structured_abstract": True,
    "keyword_count": len(keywords),
    "jel_codes": jel_codes,
    "title_page_redacted_flags": flags,
    "stage15_metadata_warnings": warnings,
    "personal_values_logged": False,
}, indent=2) + "\n", encoding="utf-8")
print("PASS: Stage-14 submission structure audit")
print("title_page_required_fields=REDACTED")
print(f"stage15_metadata_warning_count={len(warnings)}")
