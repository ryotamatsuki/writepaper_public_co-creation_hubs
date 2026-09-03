# Reproducibility Decision Log
- 2026-09-03: use Stage7.5 branch after Stage8 merge as canonical base.
- Reuse current Stage4/Stage7 verification sources through wrappers; do not reuse legacy JRS/RoRR model code.
- Use generic `article` LaTeX; journal positioning remains deferred.
- Treat `.656020/.656016` as the Stage8 frozen reporting reconciliation. Current source scripts produce values within that frozen tolerance; production records both computed values and both frozen display values rather than silently overwriting either.
- No substantive production figure is required at Stage9.
