#!/usr/bin/env python3
"""Build JE submission *preview* from frozen source, without claiming author consent.

Run from any directory. Needs latexmk, pdfinfo, pdflatex, BibTeX and TeX packages.
Does not alter manuscript sources, model inputs, archived PDFs, or evidence.
Creates a minimal preview archive; never a submission-authorised archive.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

PRODUCTION = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[3]
OUT = PRODUCTION / "_stage14_preview"
MANUSCRIPT = PRODUCTION / "manuscript"
COVER = Path(__file__).with_name("cover_letter.md")
PREVIEW = "JE_STAGE14_PREVIEW_NEEDS_AUTHOR_CONFIRMATION.zip"
INPUT_RE = re.compile(r"\\input\{([^}]+)\}")

def flatten_source(path: Path, seen: tuple[Path, ...] = ()) -> str:
    full = path.resolve()
    if full in seen:
        raise RuntimeError("Recursive input: " + str(path))
    txt = full.read_text(encoding="utf-8")
    def expand(match: re.Match[str]) -> str:
        name = match.group(1)
        child = (PRODUCTION / name).with_suffix(".tex")
        if not child.is_file():
            raise FileNotFoundError("Missing TeX input: " + str(child))
        return flatten_source(child, seen + (full,))
    return INPUT_RE.sub(expand, txt)

def run(cmd: list[str], cwd: Path) -> str:
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    output = proc.stdout + "\n" + proc.stderr
    if proc.returncode:
        raise RuntimeError(f"Failed ({proc.returncode}): {' '.join(cmd)}\n" + output[-5000:])
    return output

def hash_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def pdf_pages(path: Path) -> int:
    info = run(["pdfinfo", str(path)], path.parent)
    m = re.search(r"^Pages:\s*(\d+)\s*$", info, flags=re.M)
    if not m:
        raise RuntimeError("Cannot verify PDF page count: " + str(path))
    return int(m.group(1))

def main() -> None:
    for binary in ("latexmk", "pdflatex", "bibtex", "pdfinfo"):
        if not shutil.which(binary):
            raise RuntimeError(f"Required build dependency absent: {binary}")
    if not COVER.exists():
        raise FileNotFoundError(COVER)

    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="stage14-je-clean-") as t:
        build = Path(t)
        for name in ("main", "supplement"):
            text = flatten_source(MANUSCRIPT / f"{name}.tex")
            text = text.replace(r"\bibliography{manuscript/references}", r"\bibliography{references}")
            if INPUT_RE.search(text):
                raise RuntimeError("An unexpanded input remains in " + name)
            (build / (name + ".tex")).write_text(text, encoding="utf-8")
        shutil.copy2(MANUSCRIPT / "references.bib", build / "references.bib")

        for name in ("main", "supplement"):
            run(["latexmk", "-pdf", "-bibtex", "-interaction=nonstopmode",
                 "-halt-on-error", "-file-line-error", f"{name}.tex"], build)
        errors = []
        warnings = []
        patterns = (
            r"Citation .* undefined",
            r"Reference .* undefined",
            r"There were undefined references",
            r"There were undefined citations",
            r"multiply defined",
            r"Overfull \\hbox",
            r"Overfull \\vbox",
        )
        for name in ("main", "supplement"):
            log = (build / (name + ".log")).read_text(encoding="utf-8", errors="replace")
            for pattern in patterns:
                matches = list(re.finditer(pattern, log, flags=re.I))
                for match in matches[:8]:
                    if pattern.startswith('Overfull'):
                        errors.append(f"{name}: {pattern}: " + log[max(0,match.start()-50):match.start()+300].replace('\\n',' | '))
                    else:
                        errors.append(f"{name}: {pattern}")
            if re.search(r"Underfull \\[hv]box", log):
                warnings.append(f"{name}: Underfull box present (review manually)")
        if errors:
            raise RuntimeError("LaTeX quality gate failed: " + "; ".join(errors))

        counts = {name: pdf_pages(build / (name + ".pdf")) for name in ("main", "supplement")}
        if counts["main"] < 1 or counts["supplement"] < 1:
            raise RuntimeError("A PDF has no pages")

        submission = {
            "manuscript.pdf": build / "main.pdf",
            "main.tex": build / "main.tex",
            "ESM_1.pdf": build / "supplement.pdf",
            "supplement.tex": build / "supplement.tex",
            "references.bib": build / "references.bib",
            "cover_letter_DRAFT.txt": COVER,
        }
        zip_path = OUT / PREVIEW
        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
            for target_name, src in submission.items():
                z.write(src, arcname=target_name)

        try:
            commit = run(["git", "rev-parse", "HEAD"], REPO).strip()
        except (FileNotFoundError, RuntimeError):
            commit = "UNVERIFIED"
        report = {
            "status": "CLEAN BUILD COMPLETE; PREVIEW ONLY; AUTHOR CONFIRMATION REQUIRED",
            "source_commit": commit,
            "source_pdfs_from_current_head": True,
            "pdf_pages": counts,
            "pdf_sha256": {
                "manuscript.pdf": hash_file(build / "main.pdf"),
                "ESM_1.pdf": hash_file(build / "supplement.pdf"),
            },
            "archive": PREVIEW,
            "archive_sha256": hash_file(zip_path),
            "latex_warnings": warnings,
            "manual_visual_qa": "NOT DONE BY THIS SCRIPT",
            "author_confirmation": "REQUIRED",
            "references_metadata_primary_audit": "INCOMPLETE",
            "journal_portal_file_prompt_audit": "INCOMPLETE",
        }
        (OUT / "build_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
