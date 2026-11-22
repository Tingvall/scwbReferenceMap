.run_reference_map_transfer <- function(readers, writers, params, context) {
  query <- .refmap_read_space(readers$query_coordinates,
    context$typed_space_reads$query_coordinates)
  reference <- .refmap_read_space(readers$reference_coordinates,
    context$typed_space_reads$reference_coordinates)
  .refmap_assert(ncol(query) == ncol(reference),
    "Query and reference spaces have different dimensions.")
  verified <- readers$calibration_payload()
  calibration <- .refmap_read_calibration(verified)
  reference_cells <- as.character(unlist(calibration$reference$cell_id,
    use.names = FALSE))
  .refmap_assert(identical(reference_cells, rownames(reference)),
    "Calibration and reference-space cell identities differ.")

  neighbors <- readers$neighbors()
  .refmap_assert(inherits(neighbors, "wb_primitive") &&
      identical(as.character(neighbors$kind)[1L], "obs_table"),
    "neighbors must be an obs_table primitive.")
  ledger <- neighbors$table
  required <- c("query_cell_id", "neighbor_rank", "cell_id", "distance")
  .refmap_assert(is.data.frame(ledger) && all(required %in% names(ledger)) &&
      nrow(ledger) > 0L, "Neighbor ledger is malformed.")
  ledger$query_cell_id <- as.character(ledger$query_cell_id)
  ledger$cell_id <- as.character(ledger$cell_id)
  ledger$neighbor_rank <- as.integer(ledger$neighbor_rank)
  ledger$distance <- as.numeric(ledger$distance)
  query_cells <- rownames(query)
  .refmap_assert(identical(unique(ledger$query_cell_id), query_cells) &&
      all(ledger$cell_id %in% reference_cells) &&
      all(is.finite(ledger$distance) & ledger$distance >= 0),
    "Neighbor ledger axes do not match the authenticated spaces.")
  config <- calibration$configuration
  k_values <- as.integer(unlist(config$k_values, use.names = FALSE))
  primary_k <- as.integer(config$primary_k)
  for (cell in query_cells) {
    rows <- ledger[ledger$query_cell_id == cell, , drop = FALSE]
    .refmap_assert(identical(rows$neighbor_rank, seq_len(nrow(rows))) &&
        nrow(rows) >= max(k_values) && !anyDuplicated(rows$cell_id),
      "Neighbor ranks are incomplete, duplicated, or out of order.")
    observed <- .refmap_distance(query[cell, ],
      reference[rows$cell_id, , drop = FALSE])
    tolerance <- 1e-12 * pmax(1, observed, rows$distance)
    .refmap_assert(all(abs(observed - rows$distance) <= tolerance),
      "Neighbor distances do not authenticate against the supplied spaces.")
    expected <- order(observed, match(rows$cell_id, reference_cells),
      method = "radix")
    .refmap_assert(identical(expected, seq_len(nrow(rows))),
      "Neighbor ledger rank or candidate tie policy is invalid.")
  }

  calibrated_thresholds <- as.numeric(unlist(config$thresholds,
    use.names = FALSE))
  supported_threshold <- as.numeric(params$supported_threshold %||%
    calibrated_thresholds[[1L]])
  out_threshold <- as.numeric(params$out_of_reference_threshold %||%
    calibrated_thresholds[[2L]])
  .refmap_assert(length(supported_threshold) == 1L &&
      length(out_threshold) == 1L && is.finite(supported_threshold) &&
      is.finite(out_threshold) && supported_threshold >= 0 &&
      out_threshold >= supported_threshold,
    "Transfer thresholds must be finite with supported <= out_of_reference.")

  annotation_raw <- calibration$reference$annotation
  annotation <- if (is.null(annotation_raw)) NULL else stats::setNames(
    vapply(annotation_raw, function(x) if (is.null(x)) NA_character_ else
      as.character(x)[1L], character(1)), reference_cells)
  annotation_universe <- as.character(unlist(
    calibration$reference$annotation_categories, use.names = FALSE))
  acceptance_policy <- params$continuous_acceptance %||%
    list(mode = "support_only")
  acceptance_mode <- as.character(acceptance_policy$mode %||%
    "support_only")
  .refmap_assert(length(acceptance_mode) == 1L &&
      acceptance_mode %in% c("support_only", "support_and_annotation"),
    "continuous_acceptance mode is unsupported.")
  accepted_annotation_values <- as.character(unlist(
    acceptance_policy$accepted_annotation_values %||% character(),
    use.names = FALSE))
  minimum_annotation_confidence <- as.numeric(
    acceptance_policy$minimum_annotation_confidence %||% .5)
  if (identical(acceptance_mode, "support_and_annotation")) {
    .refmap_assert(!is.null(annotation), paste(
      "support_and_annotation continuous acceptance requires a calibrated",
      "reference annotation."))
    .refmap_assert(all(!is.na(annotation) & nzchar(annotation)), paste(
      "support_and_annotation continuous acceptance requires complete",
      "reference annotations."))
    .refmap_assert(length(accepted_annotation_values) > 0L &&
        !anyNA(accepted_annotation_values) &&
        all(nzchar(accepted_annotation_values)) &&
        !anyDuplicated(accepted_annotation_values),
      paste("accepted_annotation_values must contain unique, non-empty",
        "annotation values."))
    .refmap_assert(all(accepted_annotation_values %in% annotation_universe),
      paste("accepted_annotation_values contains values absent from the",
        "calibrated annotation universe."))
    .refmap_assert(length(minimum_annotation_confidence) == 1L &&
        is.finite(minimum_annotation_confidence) &&
        minimum_annotation_confidence >= 0 &&
        minimum_annotation_confidence <= 1,
      "minimum_annotation_confidence must be finite and in [0, 1].")
  }
  continuous_raw <- calibration$reference$continuous %||% list()
  continuous <- lapply(continuous_raw, function(values) stats::setNames(
    vapply(values, function(x) if (is.null(x)) NA_real_ else as.numeric(x)[1L],
      numeric(1)), reference_cells))
  continuous_names <- names(continuous)
  output_keys <- as.character(unlist(params$continuous_output_keys %||%
    continuous_names, use.names = FALSE))
  .refmap_assert(length(output_keys) == length(continuous_names) &&
      !anyNA(output_keys) && all(nzchar(output_keys)) &&
      !anyDuplicated(output_keys),
    "continuous_output_keys must uniquely name every calibrated vector.")
  names(output_keys) <- continuous_names

  per_query <- list(); categorical_rows <- list(); continuous_rows <- list()
  continuous_acceptance_rows <- list()
  accepted_annotation <- stats::setNames(rep(NA_character_, length(query_cells)),
    query_cells)
  accepted_continuous <- stats::setNames(lapply(continuous_names,
    function(x) stats::setNames(rep(NA_real_, length(query_cells)), query_cells)),
    continuous_names)
  for (i in seq_along(query_cells)) {
    cell <- query_cells[[i]]
    rows <- ledger[ledger$query_cell_id == cell, , drop = FALSE]
    statistic <- .refmap_distance_statistic(rows$distance, primary_k)
    status <- if (statistic <= supported_threshold) "supported" else if (
      statistic <= out_threshold) "borderline" else "out_of_reference"
    predictions <- list(); votes <- list()
    if (!is.null(annotation)) for (k in k_values) {
      take <- seq_len(k)
      vote <- .refmap_vote(annotation[rows$cell_id[take]],
        rows$distance[take], annotation_universe)
      votes[[as.character(k)]] <- vote
      predictions[[as.character(k)]] <- vote$label
      categorical_rows[[length(categorical_rows) + 1L]] <- data.frame(
        query_cell_id = cell, k = k, raw_prediction = vote$label,
        confidence = vote$confidence, top_two_margin = vote$margin,
        normalized_entropy = vote$entropy, stringsAsFactors = FALSE)
    }
    primary_vote <- votes[[as.character(primary_k)]] %||% list(
      label = NA_character_, confidence = NA_real_, margin = NA_real_,
      entropy = NA_real_)
    cross_k <- if (length(predictions)) mean(unlist(predictions,
      use.names = FALSE) == primary_vote$label, na.rm = TRUE) else NA_real_
    if (!is.null(annotation) && !identical(status, "out_of_reference"))
      accepted_annotation[[cell]] <- primary_vote$label
    primary_continuous_value <- NA_real_
    for (name in continuous_names) {
      estimates <- numeric(length(k_values)); names(estimates) <- k_values
      for (k in k_values) {
        take <- seq_len(k); weight <- .refmap_weight(rows$distance[take])
        quantile <- .refmap_weighted_quantile(
          continuous[[name]][rows$cell_id[take]], weight, c(.25, .5, .75))
        estimates[[as.character(k)]] <- quantile[[2L]]
        continuous_rows[[length(continuous_rows) + 1L]] <- data.frame(
          query_cell_id = cell, reference_key = name, k = k,
          weighted_q25 = quantile[[1L]], raw_estimate = quantile[[2L]],
          weighted_q75 = quantile[[3L]], iqr = quantile[[3L]] - quantile[[1L]],
          n_non_missing = sum(!is.na(continuous[[name]][rows$cell_id[take]])),
          stringsAsFactors = FALSE)
      }
      primary_estimate <- estimates[[as.character(primary_k)]]
      sensitivity <- if (all(is.na(estimates))) NA_real_ else
        max(abs(estimates - primary_estimate), na.rm = TRUE)
      indexes <- vapply(continuous_rows, function(x)
        identical(x$query_cell_id, cell) && identical(x$reference_key, name),
        logical(1))
      for (j in which(indexes)) continuous_rows[[j]]$cross_k_sensitivity <-
        sensitivity
      acceptance_reason <- if (identical(status, "out_of_reference")) {
        "out_of_reference"
      } else if (identical(acceptance_mode, "support_and_annotation") &&
          !primary_vote$label %in% accepted_annotation_values) {
        "annotation_not_accepted"
      } else if (identical(acceptance_mode, "support_and_annotation") &&
          (is.na(primary_vote$confidence) ||
            primary_vote$confidence < minimum_annotation_confidence)) {
        "annotation_confidence_too_low"
      } else if (!is.finite(primary_estimate)) {
        "continuous_value_unavailable"
      } else {
        "accepted"
      }
      accepted <- identical(acceptance_reason, "accepted")
      if (accepted)
        accepted_continuous[[name]][[cell]] <- primary_estimate
      continuous_acceptance_rows[[length(continuous_acceptance_rows) + 1L]] <-
        data.frame(query_cell_id = cell, reference_key = name,
          output_key = unname(output_keys[[name]]), accepted = accepted,
          acceptance_reason = acceptance_reason,
          raw_estimate = primary_estimate,
          accepted_value = if (accepted) primary_estimate else NA_real_,
          stringsAsFactors = FALSE)
      if (identical(name, as.character(config$primary_continuous_key %||% "")))
        primary_continuous_value <- if (accepted) primary_estimate else NA_real_
    }
    per_query[[i]] <- data.frame(query_cell_id = cell, support_call = status,
      distance_statistic = statistic, supported_threshold = supported_threshold,
      out_of_reference_threshold = out_threshold,
      accepted_annotation = accepted_annotation[[cell]],
      confidence = primary_vote$confidence, top_two_margin = primary_vote$margin,
      normalized_entropy = primary_vote$entropy,
      cross_k_agreement = cross_k,
      primary_continuous_key = as.character(config$primary_continuous_key %||% NA),
      primary_continuous_value = primary_continuous_value,
      stringsAsFactors = FALSE)
  }
  per_query_table <- do.call(rbind, per_query); rownames(per_query_table) <- NULL
  if (!length(categorical_rows)) categorical_rows <- list(data.frame(
    query_cell_id = query_cells, k = NA_integer_, raw_prediction = NA_character_,
    confidence = NA_real_, top_two_margin = NA_real_,
    normalized_entropy = NA_real_, stringsAsFactors = FALSE))
  categorical_table <- do.call(rbind, categorical_rows)
  if (!length(continuous_rows)) continuous_rows <- list(data.frame(
    query_cell_id = query_cells, reference_key = NA_character_, k = NA_integer_,
    weighted_q25 = NA_real_, raw_estimate = NA_real_, weighted_q75 = NA_real_,
    iqr = NA_real_, n_non_missing = NA_integer_,
    cross_k_sensitivity = NA_real_, stringsAsFactors = FALSE))
  continuous_table <- do.call(rbind, continuous_rows)
  if (length(continuous_acceptance_rows)) {
    continuous_acceptance_table <- do.call(rbind, continuous_acceptance_rows)
    rownames(continuous_acceptance_table) <- NULL
  } else continuous_acceptance_table <- data.frame(
    query_cell_id = character(), reference_key = character(),
    output_key = character(), accepted = logical(),
    acceptance_reason = character(), raw_estimate = numeric(),
    accepted_value = numeric(), stringsAsFactors = FALSE)
  grouped <- as.data.frame(table(factor(per_query_table$support_call,
    levels = c("supported", "borderline", "out_of_reference"))),
    stringsAsFactors = FALSE)
  names(grouped) <- c("support_call", "n_query")
  grouped$proportion <- grouped$n_query / length(query_cells)

  diagnostics <- list(distance_statistic = stats::setNames(
      per_query_table$distance_statistic, query_cells),
    confidence = stats::setNames(per_query_table$confidence, query_cells),
    top_two_margin = stats::setNames(per_query_table$top_two_margin, query_cells),
    normalized_entropy = stats::setNames(per_query_table$normalized_entropy,
      query_cells), cross_k_agreement = stats::setNames(
      per_query_table$cross_k_agreement, query_cells))
  numeric <- c(diagnostics, accepted_continuous)
  names(numeric)[seq_along(diagnostics)] <- paste0("diagnostic.",
    names(diagnostics))
  if (length(continuous_names))
    names(numeric)[length(diagnostics) + seq_along(continuous_names)] <-
      unname(output_keys[continuous_names])
  .refmap_assert(!anyDuplicated(names(numeric)),
    "Continuous output keys collide with diagnostic keys.")
  numeric <- Filter(function(value) any(is.finite(value)), numeric)
  numeric_cells <- lapply(numeric, function(value)
    query_cells[is.finite(value)])
  writers$declare_sparse_cell_family_layout("numeric_outputs",
    names(numeric), unname(numeric_cells))
  numeric_family <- scworkbench::wb_prim_cell_vector_family(
    key = unlist(Map(function(key, cells) rep(key, length(cells)),
      names(numeric), numeric_cells), use.names = FALSE),
    cell_id = unlist(numeric_cells, use.names = FALSE),
    value = unlist(Map(function(value, cells) value[cells], numeric,
      numeric_cells), use.names = FALSE), value_type = "continuous",
    meta = list(primary_continuous_key = config$primary_continuous_key,
      caller_output_keys = as.list(output_keys)))
  support <- stats::setNames(per_query_table$support_call, query_cells)
  writers$write_terminal_primitive("support_call",
    scworkbench::wb_prim_cell_vector(support, value_type = "categorical",
      meta = list(output_key = as.character(params$support_output_key)[1L])))
  annotation_key <- as.character(params$annotation_output_key)[1L]
  annotation_cells <- query_cells[!is.na(accepted_annotation)]
  writers$write_terminal_primitive("accepted_annotation",
    scworkbench::wb_prim_obs_table(data.frame(
      query_cell_id = annotation_cells,
      output_key = rep(annotation_key, length(annotation_cells)),
      value = unname(accepted_annotation[annotation_cells]),
      stringsAsFactors = FALSE),
      table_type = "reference_map_accepted_annotation",
      entity_kind = "query_cell"))
  writers$write_terminal_primitive("numeric_outputs", numeric_family)
  writers$write_terminal_primitive("per_query_summary",
    scworkbench::wb_prim_obs_table(per_query_table,
      table_type = "reference_map_per_query", entity_kind = "query_cell"))
  writers$write_terminal_primitive("categorical_estimates",
    scworkbench::wb_prim_obs_table(categorical_table,
      table_type = "reference_map_categorical_estimates",
      entity_kind = "query_cell_k"))
  writers$write_terminal_primitive("continuous_estimates",
    scworkbench::wb_prim_obs_table(continuous_table,
      table_type = "reference_map_continuous_estimates",
      entity_kind = "query_cell_vector_k"))
  writers$write_terminal_primitive("continuous_acceptance",
    scworkbench::wb_prim_obs_table(continuous_acceptance_table,
      table_type = "reference_map_continuous_acceptance",
      entity_kind = "query_cell_vector"))
  writers$write_terminal_primitive("grouped_summary",
    scworkbench::wb_prim_obs_table(grouped,
      table_type = "reference_map_grouped_summary", entity_kind = "support_call"))
  invisible(list(status = "ok"))
}
