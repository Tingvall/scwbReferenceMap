# scwbReferenceMap 0.2.0

- Added an optional continuous-acceptance policy that can require both
  reference support and confident membership in selected reference annotation
  values before publishing accepted continuous outputs.
- Added per-query, per-continuous-key acceptance decisions and reasons while
  retaining raw continuous estimates for every query.
- Preserved the previous support-only acceptance behavior as the default,
  including acceptance of both supported and borderline queries.
- Gated transfer now rejects calibrations with incomplete reference annotations
  so missing neighbor labels cannot inflate annotation confidence.
- Raised the minimum compatible scWorkbench version to 0.1.2.

# scwbReferenceMap 0.1.0

- Added generic leave-one-category-out reference-map calibration with
  configurable holdout metadata, k values, primary k, support quantiles,
  optional categorical annotations, and optional continuous vector families.
- Added reusable retained calibration artifacts tied to authenticated
  reference-space identities.
- Added authenticated query-to-reference transfer with supported, borderline,
  and out-of-reference calls; caller threshold overrides; categorical
  confidence diagnostics; weighted continuous summaries; raw audit estimates;
  and out-of-reference withholding.
- Added exact synthetic qualification for non-default holdout fields, category
  coverage gaps, annotation-only and continuous-only transfer, leakage
  prevention, cross-k diagnostics, identity mismatches, deterministic replay,
  promotion closure, finalized read-only reopen, and missing dependencies.
- Uses the approved ABI-v3 reference-only effective-cell declaration for the
  authenticated neighbor ledger and self-linked calibration artifact.
