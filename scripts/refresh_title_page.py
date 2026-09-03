from __future__ import annotations

from pathlib import Path
import re
import shutil
import tempfile
import zipfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "Title_Page.docx"
CURRENT_TITLE = "Strategic Interaction among Public Innovation Hubs under Private Repricing"

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def paragraph_text(p: ET.Element) -> str:
    return "".join(t.text or "" for t in p.iter(W + "t")).strip()


assert DOCX.exists(), "Title_Page.docx missing"

with zipfile.ZipFile(DOCX, "r") as zin:
    xml = zin.read("word/document.xml")
    root = ET.fromstring(xml)
    paragraphs = list(root.iter(W + "p"))
    texts = [paragraph_text(p) for p in paragraphs]

    if any(CURRENT_TITLE == t for t in texts):
        print("PASS: title page already uses current manuscript title")
        raise SystemExit(0)

    # Identify a title-like paragraph without logging its text. The existing title page
    # predates Stage 14, so we look for a substantial early paragraph containing several
    # manuscript-title concepts. This preserves all author/contact/declaration fields.
    concepts = re.compile(r"strategic|public|hub|platform|private|pricing|repric|innovation|competition", re.I)
    candidates: list[tuple[int, ET.Element, int]] = []
    for idx, (p, text) in enumerate(zip(paragraphs, texts)):
        if not text or idx > 20:
            continue
        score = len(concepts.findall(text))
        if len(text) >= 25 and score >= 2:
            candidates.append((idx, p, score))

    assert candidates, "unable to identify a title-like paragraph in Title_Page.docx"
    # Prefer the concept-richest candidate; break ties in favor of the earliest paragraph.
    candidates.sort(key=lambda item: (-item[2], item[0]))
    idx, p, score = candidates[0]

    runs = list(p.iter(W + "t"))
    assert runs, "identified title paragraph has no Word text runs"
    runs[0].text = CURRENT_TITLE
    for node in runs[1:]:
        node.text = ""

    new_xml = ET.tostring(root, encoding="utf-8", xml_declaration=True)

    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "Title_Page.docx"
        with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zout:
            for info in zin.infolist():
                data = new_xml if info.filename == "word/document.xml" else zin.read(info.filename)
                zout.writestr(info, data)
        shutil.copy2(out, DOCX)

print(f"PASS: refreshed title-page manuscript title in paragraph_index={idx}; author metadata preserved")
