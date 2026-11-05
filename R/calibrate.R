.run_reference_map_calibrate <- function(readers, writers, params, context) {
  space <- .refmap_read_space(readers$reference_coordinates,
    context$typed_space_reads$reference_coordinates)
  cells <- rownames(space)
  metadata <- readers$holdout_metadata()
  .refmap_assert(identical(as.character(unlist(metadata$cell_id,
    use.names = FALSE)), cells),
    "Holdout metadata and reference-space cell axes differ.")
  fields <- metadata$fields
  .refmap_assert(length(fields) == 1L &&
      identical(as.character(fields[[1L]]$name)[1L],
        as.character(params$holdout_group_field)[1L]),
    "The requested holdout metadata field was not authenticated.")
  groups <- vapply(fields[[1L]]$values, function(x)
    if (is.null(x)) NA_character_ else as.character(x)[1L], character(1))
  categories <- unique(groups[!is.na(groups) & nzchar(groups)])
  .refmap_assert(length(categories) >= 2L,
    "Leave-one-category-out calibration requires at least two categories.")

  annotation <- NULL
  if (is.function(readers$reference_annotation))
    annotation <- .refmap_cell_vector(readers$reference_annotation(), cells,
      "categorical")
  continuous <- list()
  if (is.function(readers$reference_continuous))
    continuous <- .refmap_family(readers$reference_continuous(), cells)
  primary_continuous <- as.character(params$primary_continuous_key %||%
    character())
  if (length(continuous)) {
    if (!length(primary_continuous) && length(continuous) == 1L)
      primary_continuous <- names(continuous)
    .refmap_assert(length(primary_continuous) == 1L &&
        primary_continuous %in% names(continuous),
      paste("primary_continuous_key must select one supplied continuous",
        "reference vector when a family is supplied."))
  } else .refmap_assert(!length(primary_continuous),
    "primary_continuous_key was supplied without continuous reference values.")

  k_values <- sort(unique(as.integer(unlist(params$k_values,
    use.names = FALSE))))
  primary_k <- as.integer(params$primary_k)
  .refmap_assert(length(k_values) > 0L && all(k_values >= 1L) &&
      primary_k %in% k_values,
    "k_values must be positive and contain primary_k.")
  quantiles <- sort(unique(as.numeric(unlist(params$support_quantiles,
    use.names = FALSE))))
  .refmap_assert(length(quantiles) == 2L && all(quantiles > 0 & quantiles < 1),
    "support_quantiles must contain two distinct probabilities in (0, 1).")

  predictions <- list(); coverage <- list(); cursor <- 0L
  annotation_universe <- if (is.null(annotation)) character() else
    unique(annotation[!is.na(annotation) & nzchar(annotation)])
  for (category in categories) {
    test_index <- which(groups == category)
    train_index <- which(!is.na(groups) & groups != category)
    .refmap_assert(length(train_index) >= max(k_values), paste0(
      "Holdout category '", category,
      "' leaves fewer training cells than max(k_values)."))
    ledger <- .refmap_exact_knn(space[test_index, , drop = FALSE],
      space[train_index, , drop = FALSE], max(k_values),
      as.integer(params$query_chunk_size))
    for (query in cells[test_index]) {
      query_rows <- ledger[ledger$query_cell_id == query, , drop = FALSE]
      row <- data.frame(cell_id = query, holdout_category = category,
        distance_statistic = .refmap_distance_statistic(query_rows$distance,
          primary_k), stringsAsFactors = FALSE)
      if (!is.null(annotation)) {
        for (k in k_values) {
          take <- seq_len(k)
          vote <- .refmap_vote(annotation[query_rows$cell_id[take]],
            query_rows$distance[take], annotation_universe)
          row[[paste0("annotation_k", k)]] <- vote$label
          row[[paste0("annotation_correct_k", k)]] <-
            !is.na(annotation[[query]]) && identical(vote$label,
              annotation[[query]])
        }
      }
      for (name in names(continuous)) for (k in k_values) {
        take <- seq_len(k)
        estimate <- .refmap_weighted_quantile(
          continuous[[name]][query_rows$cell_id[take]],
          .refmap_weight(query_rows$distance[take]), .5)
        row[[paste0(name, "_estimate_k", k)]] <- estimate
        row[[paste0(name, "_absolute_error_k", k)]] <-
          if (is.na(continuous[[name]][[query]])) NA_real_ else
            abs(estimate - continuous[[name]][[query]])
      }
      cursor <- cursor + 1L; predictions[[cursor]] <- row
    }
    absent <- if (is.null(annotation)) character() else
      setdiff(annotation_universe, unique(annotation[train_index]))
    coverage[[length(coverage) + 1L]] <- data.frame(
      holdout_category = category, n_held_out = length(test_index),
      n_training = length(train_index),
      absent_annotation_categories = paste(absent, collapse = "|"),
      stringsAsFactors = FALSE)
  }
  heldout <- do.call(rbind, predictions); rownames(heldout) <- NULL
  coverage_table <- do.call(rbind, coverage); rownames(coverage_table) <- NULL
  if (anyNA(groups) || any(!nzchar(groups[!is.na(groups)])))
    coverage_table <- rbind(coverage_table, data.frame(
      holdout_category = "<missing>",
      n_held_out = sum(is.na(groups) | !nzchar(groups)), n_training = NA_integer_,
      absent_annotation_categories = "not_calibrated", stringsAsFactors = FALSE))
  thresholds <- stats::quantile(heldout$distance_statistic,
    probs = quantiles, names = FALSE, type = 7)
  threshold_table <- data.frame(
    role = c("supported_upper", "borderline_upper"),
    quantile = quantiles, threshold = as.numeric(thresholds),
    statistic = "median_first_primary_k_distance", primary_k = primary_k,
    stringsAsFactors = FALSE)

  metric_rows <- list(data.frame(metric = "n_calibrated",
    key = "all", k = NA_integer_, value = nrow(heldout),
    stringsAsFactors = FALSE))
  if (!is.null(annotation)) for (k in k_values)
    metric_rows[[length(metric_rows) + 1L]] <- data.frame(
      metric = "annotation_accuracy", key = "reference_annotation", k = k,
      value = mean(heldout[[paste0("annotation_correct_k", k)]], na.rm = TRUE),
      stringsAsFactors = FALSE)
  for (name in names(continuous)) for (k in k_values)
    metric_rows[[length(metric_rows) + 1L]] <- data.frame(
      metric = "continuous_median_absolute_error", key = name, k = k,
      value = stats::median(heldout[[paste0(name, "_absolute_error_k", k)]],
        na.rm = TRUE), stringsAsFactors = FALSE)
  metrics <- do.call(rbind, metric_rows); rownames(metrics) <- NULL

  artifact <- list(
    schema = "scwbreferencemap.reference_map_calibration.v1",
    configuration = list(holdout_group_field = params$holdout_group_field,
      holdout_categories = as.list(categories), k_values = as.list(k_values),
      primary_k = primary_k,
      distance_statistic = "median_first_primary_k_distance",
      support_quantiles = as.list(quantiles),
      thresholds = as.list(as.numeric(thresholds)),
      primary_continuous_key = if (length(primary_continuous))
        primary_continuous else NULL),
    reference = list(cell_id = as.list(cells),
      annotation = if (is.null(annotation)) NULL else as.list(annotation),
      annotation_categories = as.list(annotation_universe),
      continuous = lapply(continuous, as.list)),
    coverage = lapply(seq_len(nrow(coverage_table)), function(i)
      as.list(coverage_table[i, , drop = FALSE])))
  tree <- writers$open_model_tree("calibration")
  writeBin(.refmap_json_bytes(artifact), tree$path(tree$primary_path))
  reference_space_id <- as.character(context$input_ids$reference_space)[1L]
  model <- scworkbench::wb_prim_model_ref(reference_space_id,
    reference_space_id, engine = "scwbreferencemap")
  model$model_type <- "reference_map_calibration"
  model$holdout_group_field <- as.character(params$holdout_group_field)[1L]
  model$primary_k <- primary_k
  tree$finish(model)
  writers$write_terminal_primitive("heldout_predictions",
    scworkbench::wb_prim_obs_table(heldout,
      table_type = "reference_map_heldout_predictions", entity_kind = "cell"))
  writers$write_terminal_primitive("calibration_metrics",
    scworkbench::wb_prim_obs_table(metrics,
      table_type = "reference_map_calibration_metrics",
      entity_kind = "calibration_metric"))
  writers$write_terminal_primitive("thresholds",
    scworkbench::wb_prim_obs_table(threshold_table,
      table_type = "reference_map_thresholds", entity_kind = "threshold"))
  writers$write_terminal_primitive("coverage_limitations",
    scworkbench::wb_prim_obs_table(coverage_table,
      table_type = "reference_map_coverage_limitations",
      entity_kind = "holdout_category"))
  invisible(list(status = "ok"))
}
