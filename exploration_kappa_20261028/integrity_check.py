#!/usr/bin/env python3
"""Verify old tracked assets are byte-identical to the PR59 source head."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BASE = 'd2520c9afbf9b2b742a5f04fbba6159eb489db61'
MANIFEST = Path(__file__).with_name('freeze_integrity.json')


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def main():
    # Check EVERY pre-existing tracked file, not only the named proof assets.
    paths = git('ls-tree', '-rz', '--name-only', BASE).decode().split('\0')
    aggregate = hashlib.sha256()
    named = {}
    count = 0
    for rel in paths:
        if not rel:
            continue
        expected = git('show', f'{BASE}:{rel}')
        actual = (ROOT/rel).read_bytes()
        assert actual == expected, f'Frozen base file changed: {rel}'
        digest = hashlib.sha256(actual).hexdigest()
        aggregate.update(rel.encode()+b'\0'+digest.encode()+b'\n')
        count += 1
        if rel.startswith(('formal/', 'hotelling_entry_experimentation/production/manuscript/')) or rel in (
            'hotelling_entry_experimentation/STAGE8_CANONICAL_THEORY_FREEZE_V1.md',
            'hotelling_entry_experimentation/STAGE8_THEORY_FREEZE_MANIFEST.md',
            'hotelling_entry_experimentation/production/reviews/stage11r_pp_global_certificate.py',
        ):
            named[rel] = digest
    out = {'base_commit': BASE, 'all_base_tracked_files_checked': count,
           'aggregate_path_and_sha256_digest': aggregate.hexdigest(),
           'named_frozen_asset_sha256': named,
           'scope': 'Byte integrity only; not a new Lean compilation or mathematical certificate'}
    if '--write' in sys.argv:
        MANIFEST.write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
    else:
        assert json.loads(MANIFEST.read_text()) == out
    print(f'PASS unchanged base files: {count}; digest={aggregate.hexdigest()}')


if __name__ == '__main__':
    main()
