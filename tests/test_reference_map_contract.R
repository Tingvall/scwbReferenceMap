test_reference_map_end_to_end_contract <- function() {
  x <- reference_map_fixture(); on.exit(reference_map_close(x), add = TRUE)
  space_before <- wb_space(x$bundle, "reference")
  metadata_before <- wb_meta(x$bundle, cells = x$reference_cells,
    cols = c("cell_id", "donor"))
  calibration <- reference_map_calibrate(x)
  stopifnot(setequal(names(calibration$run$outputs), c("calibration",
    "calibration_metrics", "coverage_limitations", "heldout_predictions",
    "thresholds")),
    identical(calibration$run$outputs$heldout_predictions$table$cell_id,
      x$reference_cells),
    identical(calibration$run$outputs$coverage_limitations$table$
      holdout_category, c("donor_a", "donor_b", "donor_c")),
    identical(calibration$run$params$holdout_group_field, "donor"),
    identical(calibration$run$params$primary_k, 3L),
    identical(wb_space(x$bundle, "reference"), space_before),
    identical(wb_meta(x$bundle, cells = x$reference_cells,
      cols = c("cell_id", "donor")), metadata_before))
  retained <- reference_map_promote_calibration(x, calibration)
  transfer <- reference_map_transfer(x, retained$model_id)
  summary <- transfer$run$outputs$per_query_summary$table
  stopifnot(identical(summary$query_cell_id, x$query_cells),
    identical(summary$support_call,
      c("supported", "borderline", "out_of_reference")),
    identical(transfer$run$outputs$support_call$cell_id, x$query_cells),
    identical(transfer$run$outputs$support_call$value,
      c("supported", "borderline", "out_of_reference")),
    identical(transfer$run$outputs$accepted_annotation$table$query_cell_id,
      x$query_cells[1:2]),
    identical(unique(transfer$run$outputs$accepted_annotation$table$output_key),
      "mapped_state"),
    !"q_out" %in% transfer$run$outputs$accepted_annotation$table$query_cell_id,
    all(c("diagnostic.distance_statistic", "mapped_trajectory") %in%
      transfer$run$outputs$numeric_outputs$key),
    !any(transfer$run$outputs$numeric_outputs$key == "mapped_trajectory" &
      transfer$run$outputs$numeric_outputs$cell_id == "q_out"),
    any(transfer$run$outputs$continuous_estimates$table$query_cell_id ==
      "q_out"),
    identical(transfer$run$params$supported_threshold, .2),
    identical(transfer$run$params$out_of_reference_threshold, 20))
  promoted <- wb_session_promote_current(x$bundle, transfer$session,
    outputs = names(transfer$run$outputs), created_by = "synthetic-contract")
  stopifnot(all(promoted$promoted$status == "ok"))
  TRUE
}

test_reference_map_optional_modes_and_negatives_contract <- function() {
  for (mode in c("annotation", "continuous")) {
    x <- reference_map_fixture(annotation = mode == "annotation",
      continuous = mode == "continuous")
    calibration <- reference_map_calibrate(x)
    stopifnot(nrow(calibration$run$outputs$heldout_predictions$table) == 12L)
    retained <- reference_map_promote_calibration(x, calibration)
    transfer <- reference_map_transfer(x, retained$model_id)
    if (mode == "annotation") stopifnot(
      nrow(transfer$run$outputs$accepted_annotation$table) == 2L,
      all(is.na(transfer$run$outputs$continuous_estimates$table$reference_key)))
    if (mode == "continuous") stopifnot(
      nrow(transfer$run$outputs$accepted_annotation$table) == 0L,
      "mapped_trajectory" %in% transfer$run$outputs$numeric_outputs$key)
    reference_map_close(x)
  }
  x <- reference_map_fixture(); on.exit(reference_map_close(x), add = TRUE)
  calibration <- reference_map_calibrate(x)
  retained <- reference_map_promote_calibration(x, calibration)
  bad <- reference_map_ledger(x)
  bad$table$distance[[1L]] <- bad$table$distance[[1L]] + 1
  rejected <- try(reference_map_transfer(x, retained$model_id, bad),
    silent = TRUE)
  reversed <- try(reference_map_transfer(x, retained$model_id,
    supported_threshold = 5, out_of_reference_threshold = 2), silent = TRUE)
  mismatched <- try(reference_map_transfer(x, retained$model_id,
    reference_space = "query"), silent = TRUE)
  stopifnot(inherits(rejected, "try-error"), inherits(reversed, "try-error"),
    inherits(mismatched, "try-error"))
  TRUE
}

test_reference_map_coverage_leakage_and_stability_contract <- function() {
  labels <- paste0("state_", rep(seq_len(4L), times = 3L))
  labels[[1L]] <- "rare_state"
  x <- reference_map_fixture(annotation_values = labels)
  on.exit(reference_map_close(x), add = TRUE)
  calibration <- reference_map_calibrate(x)
  coverage <- calibration$run$outputs$coverage_limitations$table
  stopifnot(identical(coverage$holdout_category,
      c("donor_a", "donor_b", "donor_c")),
    grepl("rare_state", coverage$absent_annotation_categories[[1L]],
      fixed = TRUE),
    identical(calibration$run$outputs$heldout_predictions$table$cell_id,
      x$reference_cells))
  retained <- reference_map_promote_calibration(x, calibration)
  transfer <- reference_map_transfer(x, retained$model_id)
  categorical <- transfer$run$outputs$categorical_estimates$table
  continuous <- transfer$run$outputs$continuous_estimates$table
  stopifnot(identical(unique(categorical$k), c(2L, 3L)),
    all(is.finite(categorical$confidence)),
    all(categorical$normalized_entropy >= 0 &
      categorical$normalized_entropy <= 1),
    all(is.finite(continuous$cross_k_sensitivity)),
    any(continuous$n_non_missing < continuous$k))

  shifted <- x$query
  shifted[,] <- shifted + 1000
  y <- reference_map_fixture(annotation_values = labels,
    query_values = shifted)
  on.exit(reference_map_close(y), add = TRUE)
  other <- reference_map_calibrate(y)
  stopifnot(identical(calibration$run$outputs$thresholds$table,
      other$run$outputs$thresholds$table),
    identical(calibration$run$outputs$calibration_metrics$table,
      other$run$outputs$calibration_metrics$table),
    identical(calibration$run$outputs$coverage_limitations$table,
      other$run$outputs$coverage_limitations$table))
  TRUE
}

test_reference_map_replay_promotion_reopen_and_dependency_contract <- function() {
  x <- reference_map_fixture()
  root <- x$root
  on.exit({
    try(wb_close_bundle(x$bundle), silent = TRUE)
    unlink(root, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  calibration <- reference_map_calibrate(x)
  retained <- reference_map_promote_calibration(x, calibration)
  first <- reference_map_transfer(x, retained$model_id)
  second <- reference_map_transfer(x, retained$model_id)
  stopifnot(identical(first$run$outputs, second$run$outputs),
    identical(first$run$outputs$per_query_summary$table$query_cell_id,
      x$query_cells))

  manifest <- jsonlite::read_json(system.file("scworkbench",
    "operation-pack.json", package = "scwbReferenceMap"),
    simplifyVector = FALSE)
  operation <- Filter(function(value) identical(value$name,
    "reference_map_transfer"), manifest$operations)[[1L]]
  stopifnot(all(vapply(operation$promotion_closure, function(value)
    setequal(unlist(value, use.names = FALSE), names(first$run$outputs)),
    logical(1))))
  promoted <- wb_session_promote_current(x$bundle, first$session,
    outputs = names(first$run$outputs), created_by = "synthetic-closure")
  stopifnot(all(promoted$promoted$status == "ok"),
    setequal(promoted$promoted$output, names(first$run$outputs)),
    isTRUE(wb_finalize_recipient_bundle(x$bundle)$ready))
  artifact_ids <- promoted$promoted$artifact_id
  wb_close_bundle(x$bundle)
  x$bundle <- wb_open_bundle(root, read_only = TRUE)
  reopened <- lapply(artifact_ids, function(id) wb_art_load(x$bundle, id))
  stopifnot(length(reopened) == length(first$run$outputs),
    all(vapply(reopened, function(value) inherits(value, "wb_primitive"),
      logical(1))))

  operation$dependencies <- c(list(list(requirement_id = "missing_fixture",
    kind = "r_package", name = "scwbDefinitelyMissingFixture",
    version = list(scheme = "r_package_version_constraint_v1",
      specifier = ">=1.0.0,<2.0.0"), disposition = "required",
    import_probe = "missing_export")), operation$dependencies)
  missing <- tryCatch(getFromNamespace(".wb_operation_pack_v2_preflight",
    "scworkbench")(operation, list(), list(), registration = list(
      identity = list(backend = "R"))), error = identity)
  stopifnot(inherits(missing, "wb_operation_pack_error"),
    grepl("unavailable", conditionMessage(missing), fixed = TRUE))
  TRUE
}
