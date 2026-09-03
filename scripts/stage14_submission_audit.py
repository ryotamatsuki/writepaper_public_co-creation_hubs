from __future__ import annotations

from pathlib import Path
import json
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "paper" / "main.tex"
SUB = ROOT / "submission" / "jpet"
TITLE_PAGE = ROOT / "Title_Page.docx"
OUT = ROOT / "generated" / "stage14" / "qa_summary.json"

EXPECTED_TITLE = "Strategic Interaction among Public Innovation Hubs under Private Repricing"
REQUIRED_PACKAGE = [
    SUB / "README.md",
    SUB / "metadata.md",
    SUB / "cover_letter.md",
    SUB / "declarations.md",
    SUB / "source_inventory.md",
    SUB / "reviewer_candidates.md",
]


def plain_docx_text(path: Path) -> str:
    with zipfile.ZipFile(path) as zf:
        xml = zf.read("word/document.xml").decode("utf-8")
    xml = re.sub(r"<w:tab[^>]*/>", "\t", xml)
    xml = re.sub(r"</w:p>", "\n", xml)
    return re.sub(r"<[^>]+>", "", xml)


main = MAIN.read_text(encoding="utf-8")

# Current JPET initial-submission metadata gates.
for label in ["Introduction.", "Methods.", "Results.", "Conclusion."]:
    assert f"\\textbf{{{label}}}" in main, f"structured abstract label missing: {label}"

kw_match = re.search(r"\\textbf\{Keywords:\}\s*(.+?)\\\\", main, re.S)
assert kw_match, "keywords line missing"
keywords = [x.strip().rstrip(".") for x in kw_match.group(1).split(";") if x.strip()]
assert 1 <= len(keywords) <= 7, f"JPET permits up to 7 keywords, found {len(keywords)}"

jel_match = re.search(r"\\textbf\{JEL Classification:\}\s*([^\n]+)", main)
assert jel_match, "JEL classification line missing"
jel_codes = re.findall(r"\b[A-Z][0-9]{2}\b", jel_match.group(1))
assert len(jel_codes) >= 2, "at least primary/secondary JEL candidates should be prepared"

assert "Data Availability Statement" in main, "Data Availability Statement missing"
assert "Generative AI tools were used" in main, "AIGC disclosure missing"

# No silent return of workflow/meta language or unresolved manuscript placeholders.
manuscript_files = [MAIN, *sorted((ROOT / "sections").glob("*.tex"))]
manuscript_text = "\n".join(p.read_text(encoding="utf-8") for p in manuscript_files)
for pat in [
    r"\bTODO\b", r"\bFIXME\b", r"\bTK\b", r"citation needed", r"\bplaceholder\b",
    r"\bStage 13\b", r"\bStage 14\b", r"\bL3-1\b", r"\breviewer attack\b",
]:
    assert not re.search(pat, manuscript_text, re.I), f"submission-source placeholder/internal marker: {pat}"

for p in REQUIRED_PACKAGE:
    assert p.exists(), f"required Stage-14 package file missing: {p.relative_to(ROOT)}"

# Redacted title-page audit: report presence only; never print names, email addresses,
# phone numbers, ORCID values, or other contact details.
assert TITLE_PAGE.exists(), "Title_Page.docx missing"
title_page_text = plain_docx_text(TITLE_PAGE)
title_page_flags = {
    "current_title_present": EXPECTED_TITLE.lower() in title_page_text.lower(),
    "email_present": bool(re.search(r"[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}", title_page_text, re.I)),
    "orcid_present": bool(re.search(r"\b\d{4}-\d{4}-\d{4}-[\dX]{4}\b", title_page_text, re.I)),
    "funding_language_present": bool(re.search(r"funding|financial support|grant|no external funding", title_page_text, re.I)),
    "conflict_language_present": bool(re.search(r"conflict of interest|competing interest|no competing", title_page_text, re.I)),
    "data_availability_language_present": bool(re.search(r"data availability", title_page_text, re.I)),
}

# A stale title page is a package blocker; other declaration flags are surfaced for
# human factual review rather than guessed by automation.
assert title_page_flags["current_title_present"], "Title_Page.docx does not contain current manuscript title"
assert title_page_flags["email_present"], "Title_Page.docx lacks an identifiable author email field"
assert title_page_flags["orcid_present"], "Title_Page.docx lacks an ORCID-shaped field"

summary = {
    "expected_title": EXPECTED_TITLE,
    "structured_abstract": True,
    "keyword_count": len(keywords),
    "jel_codes": jel_codes,
    "required_package_files": [str(p.relative_to(ROOT)) for p in REQUIRED_PACKAGE],
    "title_page_redacted_flags": title_page_flags,
    "pii_emitted": False,
}
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

print("PASS: Stage-14 submission structure audit")
print(f"keywords={len(keywords)}; JEL={','.join(jel_codes)}")
for key, value in title_page_flags.items():
    print(f"title_page.{key}={value}")
print("title_page.PII=REDACTED")
