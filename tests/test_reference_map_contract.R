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
    identical(transfer$run$outputs$continuous_acceptance$table$
      acceptance_reason,
      c("accepted", "accepted", "out_of_reference")),
    identical(transfer$run$params$supported_threshold, .2),
    identical(transfer$run$params$out_of_reference_threshold, 20),
    identical(transfer$run$params$continuous_acceptance$mode,
      "support_only"))
  promoted <- wb_session_promote_current(x$bundle, transfer$session,
    outputs = names(transfer$run$outputs), created_by = "synthetic-contract")
  stopifnot(all(promoted$promoted$status == "ok"))
  TRUE
}

test_reference_map_continuous_acceptance_contract <- function() {
  x <- reference_map_fixture()
  on.exit(reference_map_close(x), add = TRUE)
  calibration <- reference_map_calibrate(x)
  retained <- reference_map_promote_calibration(x, calibration)
  policy <- list(mode = "support_and_annotation",
    accepted_annotation_values = "state_1",
    minimum_annotation_confidence = .5)
  transfer <- reference_map_transfer(x, retained$model_id,
    continuous_acceptance = policy)
  acceptance <- transfer$run$outputs$continuous_acceptance$table
  stopifnot(identical(acceptance$query_cell_id, x$query_cells),
    identical(acceptance$acceptance_reason,
      c("accepted", "annotation_not_accepted", "out_of_reference")),
    identical(acceptance$accepted, c(TRUE, FALSE, FALSE)),
    all(is.finite(acceptance$raw_estimate)),
    identical(transfer$run$outputs$numeric_outputs$cell_id[
      transfer$run$outputs$numeric_outputs$key == "mapped_trajectory"],
      "q_supported"),
    all(x$query_cells %in%
      transfer$run$outputs$continuous_estimates$table$query_cell_id))

  labels <- rep("out_domain", length(x$reference_cells))
  labels[c(1L, 5L)] <- "in_domain"
  y <- reference_map_fixture(annotation_values = labels)
  on.exit(reference_map_close(y), add = TRUE)
  calibration_y <- reference_map_calibrate(y)
  retained_y <- reference_map_promote_calibration(y, calibration_y)
  low <- reference_map_transfer(y, retained_y$model_id,
    continuous_acceptance = list(mode = "support_and_annotation",
      accepted_annotation_values = "in_domain",
      minimum_annotation_confidence = .9))
  low_acceptance <- low$run$outputs$continuous_acceptance$table
  stopifnot(identical(low_acceptance$acceptance_reason[[1L]],
      "annotation_confidence_too_low"),
    !low_acceptance$accepted[[1L]], is.finite(low_acceptance$raw_estimate[[1L]]))

  z <- reference_map_fixture(continuous_values = rep(NA_real_, 12L))
  on.exit(reference_map_close(z), add = TRUE)
  calibration_z <- suppressWarnings(reference_map_calibrate(z))
  retained_z <- reference_map_promote_calibration(z, calibration_z)
  unavailable <- reference_map_transfer(z, retained_z$model_id)
  unavailable_acceptance <-
    unavailable$run$outputs$continuous_acceptance$table
  stopifnot(identical(unavailable_acceptance$acceptance_reason,
      c("continuous_value_unavailable", "continuous_value_unavailable",
        "out_of_reference")),
    all(is.na(unavailable_acceptance$raw_estimate)))

  invalid <- try(reference_map_transfer(x, retained$model_id,
    continuous_acceptance = list(mode = "support_and_annotation",
      accepted_annotation_values = "unknown", minimum_annotation_confidence = .5)),
    silent = TRUE)
  incomplete_labels <- paste0("state_", rep(seq_len(4L), times = 3L))
  incomplete_labels[c(5L, 9L)] <- NA_character_
  incomplete <- reference_map_fixture(annotation_values = incomplete_labels)
  on.exit(reference_map_close(incomplete), add = TRUE)
  incomplete_calibration <- reference_map_calibrate(incomplete)
  incomplete_retained <- reference_map_promote_calibration(incomplete,
    incomplete_calibration)
  incomplete_transfer <- try(reference_map_transfer(incomplete,
    incomplete_retained$model_id,
    continuous_acceptance = list(mode = "support_and_annotation",
      accepted_annotation_values = "state_1",
      minimum_annotation_confidence = .5)), silent = TRUE)
  stopifnot(inherits(invalid, "try-error"),
    inherits(incomplete_transfer, "try-error"),
    grepl("requires complete reference annotations",
      as.character(incomplete_transfer), fixed = TRUE))
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
