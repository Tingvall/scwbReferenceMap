# scwbReferenceMap

`scwbReferenceMap` is the official, dataset-agnostic scWorkbench operation
pack for calibrating a frozen reference and mapping other cells into it.
Installation is inert; activate the pack through the explicit scWorkbench
operation-pack boundary.

## Operations

- `scwbreferencemap::reference_map_calibrate` performs leave-one-category-out
  calibration using a configurable categorical metadata field (default
  `orig.ident`). It accepts an optional categorical annotation and optional
  named continuous vectors. A single supplied continuous vector is selected as
  primary automatically; for a family, set `primary_continuous_key` to the
  chosen trajectory/pseudotime or other continuous axis.
- `scwbreferencemap::reference_map_transfer` authenticates an externally
  produced query-to-reference neighbor ledger against the frozen reference,
  applies calibrated or caller-overridden distance thresholds, and emits
  `supported`, `borderline`, or `out_of_reference` calls. It retains raw
  categorical and continuous estimates but withholds accepted transferred
  values for out-of-reference queries.

The provisional defaults are `k = 15, 30, 50`, primary `k = 30`, median
first-30 neighbor distance, and q95/q99 support thresholds. Transfer may
override the two numeric thresholds without changing the retained calibration.
Categorical transfer reports confidence, top-two margin, normalized entropy,
and cross-k agreement. Continuous transfer uses
`exp(-(distance / kth_distance)^2)` weights and reports weighted median, q25,
q75, IQR, and cross-k sensitivity while preserving missing reference values.
Accepted continuous values are published as a sparse typed vector family under
caller-supplied keys. Accepted categorical values are published as a typed
table with `query_cell_id`, caller-supplied `output_key`, and `value`; withheld
out-of-reference cells are absent rather than assigned a synthetic label.

## Reference identity and leakage protection

Calibration consumes the authenticated reference space plus a reference matrix
used only to authorize the selected metadata field. Query information is not a
calibration input. The retained calibration model is linked to the exact
reference-space identity and contains the calibrated annotation/continuous
values. Transfer reauthenticates the model linkage, reference axis, query and
reference coordinates, neighbor ranks, candidate tie order, and every ledger
distance before producing outputs.

The pack never modifies source bundles directly. Publication remains a
scWorkbench promotion step and all outputs from each operation form one atomic
promotion closure.

## Core support

ABI v3 now admits an explicit `effective_cell_coverage: "reference_only"`
declaration for the authenticated neighbor ledger and frozen self-linked
calibration model. They remain fully identity- and provenance-bound but do not
filter the disjoint query-cell axis. See
[CORE_LIMITATION.md](CORE_LIMITATION.md) for the resolved preflight limitation.
PyBench is unchanged.
