from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "generated/results/manifest.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if not MANIFEST.exists():
        raise RuntimeError("manifest missing; run make manifest")
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    files = data.get("files", {})
    if not files:
        raise RuntimeError("manifest contains no files")

    for rel, meta in sorted(files.items()):
        path = ROOT / rel
        if not path.exists():
            raise RuntimeError(f"manifest target missing: {rel}")
        actual_bytes = len(path.read_bytes())
        actual_hash = sha256(path)
        if int(meta["bytes"]) != actual_bytes:
            raise RuntimeError(
                f"manifest byte-count mismatch for {rel}: {meta['bytes']} != {actual_bytes}"
            )
        if meta["sha256"] != actual_hash:
            raise RuntimeError(
                f"manifest hash mismatch for {rel}: {meta['sha256']} != {actual_hash}"
            )

    print(f"MANIFEST_INTEGRITY: PASS ({len(files)} files)")


if __name__ == "__main__":
    main()
