reference_map_activate_pack <- function() {
  library(scworkbench)
  qualified <- "scwbreferencemap::reference_map_calibrate"
  if (qualified %in% wb_list_ops()) return(invisible(TRUE))
  path <- system.file("scworkbench", "operation-pack.json",
    package = "scwbReferenceMap")
  if (!nzchar(path)) stop("scwbReferenceMap is not installed")
  wb_operation_pack_activate(path, lib.loc = .libPaths())
  invisible(TRUE)
}

reference_map_fixture <- function(annotation = TRUE, continuous = TRUE,
    annotation_values = NULL, continuous_values = NULL, query_values = NULL) {
  reference_map_activate_pack()
  root <- tempfile("scwbReferenceMap-test-")
  bundle <- wb_init_bundle(root, overwrite = TRUE)
  reference_cells <- paste0("r", seq_len(12L))
  query_cells <- c("q_supported", "q_borderline", "q_out")
  wb_define_cells(bundle, c(reference_cells, query_cells))
  state_index <- rep(seq_len(4L), times = 3L)
  donor_index <- rep(0:2, each = 4L)
  reference <- cbind(x = state_index + donor_index * .04,
    y = (state_index %% 2L) + donor_index * .03)
  rownames(reference) <- reference_cells
  query <- rbind(q_supported = reference[1L, ] + c(.01, .01),
    q_borderline = c(10, 10), q_out = c(50, 50))
  if (!is.null(query_values)) {
    query <- as.matrix(query_values)
    stopifnot(identical(rownames(query), query_cells), ncol(query) == 2L)
  }
  wb_add_space(bundle, "reference", reference)
  wb_add_space(bundle, "query", query)
  counts <- Matrix::Matrix(matrix(seq_len(24L), nrow = 2L,
    dimnames = list(c("g1", "g2"), reference_cells)), sparse = TRUE)
  wb_attach_matrix(bundle, "RNA", counts, feature_type = "gene")
  wb_attach_metadata(bundle, data.frame(cell_id = reference_cells,
    donor = rep(c("donor_a", "donor_b", "donor_c"), each = 4L),
    stringsAsFactors = FALSE))
  labels <- if (is.null(annotation_values)) paste0("state_", state_index) else
    annotation_values
  stopifnot(length(labels) == length(reference_cells))
  annotation_value <- if (annotation) scworkbench::wb_prim_cell_vector(
    stats::setNames(labels, reference_cells),
    value_type = "categorical") else NULL
  trajectory <- if (is.null(continuous_values))
    rep(seq(0, 1, length.out = 4L), times = 3L) else continuous_values
  stopifnot(length(trajectory) == length(reference_cells))
  if (is.null(continuous_values)) trajectory[[12L]] <- NA_real_
  continuous_value <- if (continuous)
    scworkbench::wb_prim_cell_vector_family(rep("trajectory", 12L),
      reference_cells, trajectory, value_type = "continuous") else NULL
  list(root = root, bundle = bundle, reference_cells = reference_cells,
    query_cells = query_cells, reference = reference, query = query,
    annotation = annotation_value, continuous = continuous_value)
}

reference_map_close <- function(x) {
  try(wb_close_bundle(x$bundle), silent = TRUE)
  unlink(x$root, recursive = TRUE, force = TRUE)
}

reference_map_calibrate <- function(x) {
  inputs <- list(reference_matrix = "RNA", reference_space = "reference",
    holdout_group_field = "donor", k_values = c(2L, 3L), primary_k = 3L,
    support_quantiles = list(.8, .95), query_chunk_size = 2L,
    space_chunk_size = 3L, min_cells = 1L)
  if (!is.null(x$annotation)) inputs$reference_annotation <- x$annotation
  if (!is.null(x$continuous)) {
    inputs$reference_continuous <- x$continuous
    inputs$primary_continuous_key <- "trajectory"
  }
  wb_session_run_op(wb_session_new(), x$bundle,
    "scwbreferencemap::reference_map_calibrate", inputs = inputs)
}

reference_map_ledger <- function(x) {
  rows <- lapply(seq_along(x$query_cells), function(i) {
    distance <- sqrt(rowSums((x$reference - matrix(x$query[i, ],
      nrow = nrow(x$reference), ncol = ncol(x$reference), byrow = TRUE))^2))
    order_index <- order(distance, seq_along(distance), method = "radix")[1:3]
    data.frame(query_cell_id = x$query_cells[[i]], neighbor_rank = 1:3,
      cell_id = x$reference_cells[order_index], distance = distance[order_index],
      stringsAsFactors = FALSE)
  })
  scworkbench::wb_prim_obs_table(do.call(rbind, rows),
    table_type = "knn_neighbors", entity_kind = "cell",
    meta = list(space_id = "reference", k = 3L, include_query = FALSE))
}

reference_map_promote_calibration <- function(x, calibration) {
  promoted <- wb_session_promote_current(x$bundle, calibration$session,
    outputs = names(calibration$run$outputs), created_by = "synthetic-contract")
  stopifnot(all(promoted$promoted$status == "ok"))
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "calibration"][[1L]]
  list(promotion = promoted, model_id = model_id)
}

reference_map_transfer <- function(x, model_id, ledger = reference_map_ledger(x),
    supported_threshold = .2, out_of_reference_threshold = 20,
    reference_space = "reference", continuous_acceptance = NULL) {
  inputs <- list(
      in_space = "query", reference_space = reference_space,
      neighbors = ledger,
      calibration = wb_ref_artifact(model_id, "model"),
      supported_threshold = supported_threshold,
      out_of_reference_threshold = out_of_reference_threshold,
      support_output_key = "map_support",
      annotation_output_key = "mapped_state",
      cell_ids = x$query_cells, min_cells = 1L, space_chunk_size = 2L)
  if (!is.null(x$continuous))
    inputs$continuous_output_keys <- "mapped_trajectory"
  if (!is.null(continuous_acceptance))
    inputs$continuous_acceptance <- continuous_acceptance
  wb_session_run_op(wb_session_new(), x$bundle,
    "scwbreferencemap::reference_map_transfer", inputs = inputs)
}
