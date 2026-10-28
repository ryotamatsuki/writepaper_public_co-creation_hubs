#!/usr/bin/env python3
"""Run the original checks in isolation, check freeze integrity, and log scope."""
from pathlib import Path
import subprocess, tempfile, shutil, json, hashlib, re, sys

HERE=Path(__file__).resolve().parent
PROD=HERE.parent
ROOT=PROD.parents[1]
BASE='d2520c9afbf9b2b742a5f04fbba6159eb489db61'
NAMES=['stage11r_pp_global_certificate.py','stage11r_p1_open_region_certificate.py',
       'stage11r_p1b_empirical_implications_check.py','stage11r_p1c_reduction_check.py']

def main():
    protected=['formal/Stage33Formal.lean',
      'hotelling_entry_experimentation/STAGE8_CANONICAL_THEORY_FREEZE_V1.md',
      'hotelling_entry_experimentation/STAGE8_THEORY_FREEZE_MANIFEST.md',
      'hotelling_entry_experimentation/production/verification/stage33_reproduce.py']
    protected += ['hotelling_entry_experimentation/production/reviews/'+n for n in NAMES]
    hashes={}
    for path in protected:
        original=subprocess.check_output(['git','show',BASE+':'+path],cwd=ROOT)
        now=(ROOT/path).read_bytes()
        assert original==now,path
        hashes[path]=hashlib.sha256(now).hexdigest()
    logs=[]
    with tempfile.TemporaryDirectory(prefix='jict_checks_') as td:
        td=Path(td)
        for name in NAMES:shutil.copy2(PROD/'reviews'/name,td/name)
        shutil.copy2(PROD/'verification/stage33_reproduce.py',td/'stage33_reproduce.py')
        for name in NAMES+['stage33_reproduce.py']:
            p=subprocess.run([sys.executable,name],cwd=td,capture_output=True,text=True)
            logs.append('CHECK '+name+'\n'+p.stdout+p.stderr)
            assert p.returncode==0,name
    p=subprocess.run([sys.executable,str(HERE/'diagnostic_check.py')],cwd=ROOT,capture_output=True,text=True)
    logs.append('CHECK diagnostic_check.py\n'+p.stdout+p.stderr)
    assert p.returncode==0
    # Isolated flat replication must run without the surrounding repository.
    with tempfile.TemporaryDirectory(prefix='jict_flat_') as td:
        td=Path(td)
        for name in NAMES+['diagnostic_check.py']:
            shutil.copy2(HERE/'submission_files'/name,td/name)
        p=subprocess.run([sys.executable,'diagnostic_check.py'],cwd=td,capture_output=True,text=True)
        assert p.returncode==0,p.stdout+p.stderr
        logs.append('PASS isolated flat diagnostic replication\n')
    (HERE/'validation.txt').write_text('\n'.join(logs))
    result={'base_commit':BASE,'unchanged_scientific_sources_sha256':hashes,
            'original_checks_passed':NAMES+['stage33_reproduce.py'],
            'exact_diagnostic_passed':True,'isolated_flat_diagnostic_passed':True,
            'new_lean_or_full_lambda_certificate':False}
    (HERE/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS original five checks, exact diagnostic, flat reproduction, unchanged frozen sources')

if __name__=='__main__':main()
