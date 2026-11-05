# Resolved core limitation: non-filtering authenticated ledger inputs

## Observed blocker

`reference_map_transfer` needs four independently authenticated inputs:

1. a complete query space;
2. a complete, disjoint reference space;
3. the retained calibration model linked to the reference space; and
4. the `knn_neighbors` observation table whose `query_cell_id` is on the query
   axis and whose `cell_id` is on the reference/candidate axis.

ABI v3 can correctly declare both spaces as complete typed-space reads and can
declare the calibration as a contained-model read. It has no declaration for a
non-filtering authenticated ordinary input. During the v3-to-v2 catalog
projection, the neighbor `obs_table` therefore participates in effective-cell
intersection. Core primitive-axis extraction does not treat `query_cell_id` as
the table's cell axis; the candidate `cell_id` values filter out every disjoint
query cell. Execution is rejected before the pack runner receives its readers:

```text
requested cells: 3
effective cells: 0
neighbors (obs_table): 0/3
```

Changing the table so `cell_id` contains query cells would violate the existing
authenticated `knn_neighbors` ledger contract. Adding mapping behavior to
`knn_neighbors`, weakening the query/reference separation, or encoding the
ledger as an unrelated model would also violate the approved design.

## Implemented primitive

The operation-pack API now lets a required, alignment-none `obs_table` declare
itself **non-filtering/reference-only for effective-cell resolution** while
remaining fully identity-bound and available to the runner. The same explicit
declaration is admitted for a required `model_ref` only when a self-linked
contained-model reader independently verifies its durable source and fitted
cell linkage.

## Approved core change

The `neighbors` input now declares:

```json
"effective_cell_coverage": "reference_only"
```

Admission permits it only for a required `obs_table` whose manifest alignment
is `none`, or a required, authenticated self-linked `model_ref`. The v3 catalog
projection sets that input's coverage to `reference_only`; request construction,
identity verification, reader delivery, reopen validation, and provenance
recording remain unchanged.

The broader inference—automatically projecting every ABI-v3 input with
`alignment: none` as `reference_only`—was not used. The explicit declaration is
safer because it limits compatibility impact. PyBench remains unchanged.
