#!/usr/bin/env python3
"""Build anonymous, flat LaTeX review sources and two inspected PDF targets."""
from pathlib import Path
import subprocess, tempfile, shutil, re, json, os, zipfile
from datetime import datetime, timezone

HERE=Path(__file__).resolve().parent
PRODUCTION=HERE.parent
SOURCE=HERE/'submission_files'
OUT=HERE/'output/pdf'

def flatten(path):
    text=path.read_text()
    def inc(m):
        target=PRODUCTION/m.group(1)
        if not target.suffix:target=target.with_suffix('.tex')
        return flatten(target)
    return re.sub(r'\\input\{([^}]+)\}',inc,text)

def main():
    SOURCE.mkdir(exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    for name in ('main','supplement'):
        text=flatten(PRODUCTION/f'manuscript/{name}.tex')
        text=text.replace('\\bibliography{manuscript/references}','\\bibliography{references}')
        (SOURCE/f'{name}.tex').write_text(text)
    shutil.copy2(PRODUCTION/'manuscript/references.bib',SOURCE/'references.bib')
    for name in ('stage11r_pp_global_certificate.py','stage11r_p1_open_region_certificate.py',
                 'stage11r_p1b_empirical_implications_check.py','stage11r_p1c_reduction_check.py'):
        shutil.copy2(PRODUCTION/'reviews'/name,SOURCE/name)
    shutil.copy2(PRODUCTION/'verification/stage33_reproduce.py',SOURCE/'stage33_reproduce.py')
    shutil.copy2(HERE/'diagnostic_check.py',SOURCE/'diagnostic_check.py')
    for name in ('incidence_screening_sign_certificate.py','incidence_screening_valuation_check.py'):
        shutil.copy2(PRODUCTION/'government_game_exploration'/name,SOURCE/name)
    result={'scope':'Anonymous review PDFs; external review and author declarations still pending',
            'abstract_words':len((SOURCE/'main.tex').read_text().split('\\begin{abstract}')[1].split('\\end{abstract}')[0].split()),
            'files':[]}
    with tempfile.TemporaryDirectory(prefix='jict_p1d_') as tmp:
        tmp=Path(tmp)
        for p in SOURCE.iterdir():
            if p.suffix in ('.tex','.bib'):shutil.copy2(p,tmp/p.name)
        for name,dest in [('main','jict_revised_manuscript.pdf'),('supplement','jict_online_supplement.pdf')]:
            env=dict(os.environ,SOURCE_DATE_EPOCH=str(int(datetime(2026,10,28,tzinfo=timezone.utc).timestamp())),FORCE_SOURCE_DATE='1')
            run=subprocess.run(['latexmk','-pdf','-bibtex','-interaction=nonstopmode','-halt-on-error',f'{name}.tex'],cwd=tmp,capture_output=True,text=True,env=env)
            (HERE/f'{name}_build.log').write_text(run.stdout+run.stderr)
            if run.returncode:raise RuntimeError(f'{name} build failed; see {HERE/name}_build.log')
            log=(tmp/f'{name}.log').read_text(errors='replace')
            bad=re.findall(r'.*(?:undefined|multiply defined|Overfull).*',log)
            # Underfull boxes may be deliberate ragged tables or references;
            # they are retained for visual review, never silently classified.
            under=re.findall(r'.*Underfull.*',log)
            shutil.copy2(tmp/f'{name}.pdf',OUT/dest)
            info=subprocess.check_output(['pdfinfo',str(OUT/dest)],text=True)
            pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M).group(1))
            result['files'].append({'target':dest,'pages':pages,'unresolved_or_overfull':bad,'underfull_boxes':under})
            if bad:raise AssertionError(f'{name}: reference/layout warnings: {bad}')
    assert 150<=result['abstract_words']<=250
    assert result['files'][0]['pages']<=40
    blocked=re.compile(r'ryotamatsuki|ryota\.matsuki|Ryota Matsuki',re.I)
    for p in SOURCE.iterdir():
        if p.suffix in ('.tex','.bib','.py','.md'):
            assert not blocked.search(p.read_text()),f'Identifying source: {p.name}'
    for p in OUT.glob('*.pdf'):
        extracted=subprocess.check_output(['pdftotext',str(p),'-'],text=True)
        assert not blocked.search(extracted),f'Identifying PDF: {p.name}'
    result['anonymous_source_and_pdf_scan']='PASS'
    package=HERE/'output/jict_anonymous_review_package.zip'
    items=[p for p in SOURCE.iterdir() if p.suffix in ('.tex','.bib','.py','.md')]+list(OUT.glob('*.pdf'))
    with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(items,key=lambda p:p.name):
            meta=zipfile.ZipInfo(p.name,date_time=(2026,10,28,0,0,0))
            meta.compress_type=zipfile.ZIP_DEFLATED
            z.writestr(meta,p.read_bytes())
    result['flat_anonymous_package_files']=sorted(p.name for p in items)
    (HERE/'build_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
