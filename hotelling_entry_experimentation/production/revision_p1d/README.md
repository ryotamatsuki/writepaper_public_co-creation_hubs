# JICT frozen-model revision responding to the independent final review

Base: PR #59 head `d2520c9afbf9b2b742a5f04fbba6159eb489db61`.
The kappa extension in PR #60 is not included.
No primitive, player, policy instrument, timing, research technology, product game,
original numerical setting or Lean source is changed.

Read [RESPONSE_TO_REVIEW.md](RESPONSE_TO_REVIEW.md) for the point-by-point response
and remaining academic/publication limits.

The revised main result is a model-specific diagnostic: after repricing under
skill-dependent and equal-full-participation-mean success, the public-access
threshold gap separates private exclusion and regional technological incidence.
The benchmark's public entry is driven by restored direct service value despite
negative local commercialisation returns. The narrow cost band is quantified,
not hidden or enlarged by recalibration.

Outputs:

- `output/pdf/jict_revised_manuscript.pdf`: anonymous, double-spaced manuscript.
- `output/pdf/jict_online_supplement.pdf`: original proofs and the diagnostic proof.
- `submission_files/`: flat LaTeX and standalone standard-library replication.
- `diagnostic_results.json`: exact frozen-setting decomposition.
- `validation.json` and `validation.txt`: original checks and protected-source hashes.
- `build_checks.json`: page count, abstract length and LaTeX warning audit.
- `visual_qa.json`: actual final PDFs and all-page visual inspection record.
- `AUTHOR_INFORMATION_REVIEW.md`: author-specific facts requiring confirmation,
  kept outside the anonymous review materials.

From repository root:

```bash
python hotelling_entry_experimentation/production/revision_p1d/build_review.py
python hotelling_entry_experimentation/production/revision_p1d/validate_revision.py
```

The build uses latexmk, a standard LaTeX installation with natbib/endfloat/xr-hyper,
and pdfinfo. Python checks need only the standard library. Compiled PDF metadata
uses a fixed source date for reproducibility.

Do not equate successful checks with historical priority, an external-human
review, all-lambda SPNE, new Lean coverage, or permission to submit/merge.
The main text remains above JICT's typical 20-page preference (but below its
40-page maximum), and the scholarly novelty gate is not automatically cleared.
