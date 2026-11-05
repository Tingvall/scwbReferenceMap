`%||%` <- function(x, y) if (is.null(x)) y else x

.refmap_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

.refmap_require <- function(package, purpose) {
  if (!requireNamespace(package, quietly = TRUE))
    stop(sprintf("Package '%s' is required for %s.", package, purpose),
      call. = FALSE)
}

.refmap_read_space <- function(reader, capability) {
  chunks <- scworkbench::wb_partition_axis(as.character(capability$cells),
    as.integer(capability$requested_chunk_size))
  blocks <- lapply(chunks, function(cells) reader(cells = cells)$values)
  values <- do.call(rbind, blocks)
  rownames(values) <- as.character(capability$cells)
  .refmap_assert(is.matrix(values) && nrow(values) > 0L && ncol(values) > 0L &&
      all(is.finite(values)), "Reference-map space coordinates are malformed.")
  values
}

.refmap_distance <- function(query, candidate) {
  sqrt(rowSums((candidate - matrix(query, nrow = nrow(candidate),
    ncol = ncol(candidate), byrow = TRUE))^2))
}

.refmap_exact_knn <- function(query, candidate, k, chunk_size) {
  .refmap_require("RANN", "exact reference-map neighbor search")
  .refmap_assert(ncol(query) == ncol(candidate) && nrow(candidate) > 0L,
    "Query and candidate reference-map spaces must share dimensions.")
  keep <- min(as.integer(k), nrow(candidate))
  search <- min(nrow(candidate), keep + 1L)
  chunks <- scworkbench::wb_partition_axis(rownames(query),
    as.integer(chunk_size))
  rows <- list(); cursor <- 0L
  for (cells in chunks) {
    qi <- match(cells, rownames(query))
    nearest <- RANN::nn2(candidate, query[qi, , drop = FALSE], k = search,
      eps = 0)
    index <- nearest$nn.idx
    if (is.null(dim(index))) index <- matrix(index, nrow = length(cells),
      byrow = TRUE)
    index <- as.matrix(index)
    .refmap_assert(identical(as.integer(dim(index)),
      c(as.integer(length(cells)), as.integer(search))),
      "RANN returned an invalid exact-neighbor index matrix.")
    for (j in seq_along(cells)) {
      selected <- as.integer(index[j, ])
      distance <- .refmap_distance(query[qi[[j]], ],
        candidate[selected, , drop = FALSE])
      order_index <- order(distance, selected, method = "radix")
      selected <- selected[order_index]; distance <- distance[order_index]
      if (search > keep && identical(distance[[keep]], distance[[keep + 1L]])) {
        distance <- .refmap_distance(query[qi[[j]], ], candidate)
        selected <- order(distance, seq_along(distance), method = "radix")
        selected <- selected[seq_len(keep)]
        distance <- distance[selected]
      } else {
        selected <- selected[seq_len(keep)]
        distance <- distance[seq_len(keep)]
      }
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(query_cell_id = cells[[j]],
        neighbor_rank = seq_len(keep), cell_id = rownames(candidate)[selected],
        distance = as.numeric(distance), stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

.refmap_cell_vector <- function(value, cells, value_type) {
  .refmap_assert(inherits(value, "wb_primitive") &&
      identical(as.character(value$kind)[1L], "cell_vector"),
    "Expected a cell_vector primitive.")
  index <- match(cells, as.character(value$cell_id))
  .refmap_assert(length(value$cell_id) == length(cells) && !anyNA(index) &&
      !anyDuplicated(value$cell_id),
    "Reference vector must exactly cover the authenticated reference axis.")
  raw <- value$value[index]
  if (identical(value_type, "continuous")) {
    result <- suppressWarnings(as.numeric(raw))
    .refmap_assert(all(is.na(raw) | !is.na(result)) &&
        all(is.finite(result[!is.na(result)])),
      "Continuous reference values must be finite numeric values or missing.")
  } else result <- as.character(raw)
  stats::setNames(result, cells)
}

.refmap_family <- function(value, cells) {
  .refmap_assert(inherits(value, "wb_primitive") &&
      identical(as.character(value$kind)[1L], "cell_vector_family"),
    "Expected a cell_vector_family primitive.")
  key <- as.character(value$key); cell <- as.character(value$cell_id)
  keys <- unique(key)
  .refmap_assert(length(keys) > 0L && !anyNA(keys) && all(nzchar(keys)) &&
      !anyDuplicated(paste(key, cell, sep = "\r")),
    "Continuous reference family keys are empty or duplicated.")
  out <- lapply(keys, function(name) {
    take <- key == name
    .refmap_assert(sum(take) == length(cells) &&
        setequal(cell[take], cells),
      "Each continuous reference vector must exactly cover the reference axis.")
    index <- match(cells, cell[take])
    result <- as.numeric(value$value[take][index])
    .refmap_assert(all(is.na(result) | is.finite(result)),
      "Continuous reference values must be finite or missing.")
    stats::setNames(result, cells)
  })
  stats::setNames(out, keys)
}

.refmap_weight <- function(distance) {
  kth <- distance[[length(distance)]]
  if (kth == 0) as.numeric(distance == 0) else exp(-(distance / kth)^2)
}

.refmap_weighted_quantile <- function(value, weight, probabilities) {
  keep <- !is.na(value) & is.finite(value) & is.finite(weight) & weight > 0
  if (!any(keep)) return(rep(NA_real_, length(probabilities)))
  value <- value[keep]; weight <- weight[keep]
  order_index <- order(value, seq_along(value), method = "radix")
  value <- value[order_index]; weight <- weight[order_index]
  cumulative <- cumsum(weight) / sum(weight)
  vapply(probabilities, function(probability)
    value[[which(cumulative >= probability)[[1L]]]], numeric(1))
}

.refmap_vote <- function(label, distance, universe) {
  keep <- !is.na(label) & nzchar(label)
  if (!any(keep)) return(list(label = NA_character_, confidence = NA_real_,
    margin = NA_real_, entropy = NA_real_))
  label <- label[keep]; weight <- .refmap_weight(distance)[keep]
  score <- stats::setNames(rep(0, length(universe)), universe)
  for (i in seq_along(label)) score[[label[[i]]]] <-
    score[[label[[i]]]] + weight[[i]]
  ranked <- order(-score, seq_along(score), method = "radix")
  total <- sum(score)
  probability <- if (total > 0) score / total else score
  positive <- probability[probability > 0]
  entropy <- if (length(universe) <= 1L || !length(positive)) 0 else
    -sum(positive * log(positive)) / log(length(universe))
  list(label = names(score)[ranked[[1L]]],
    confidence = unname(probability[ranked[[1L]]]),
    margin = unname(probability[ranked[[1L]]] -
      if (length(ranked) > 1L) probability[ranked[[2L]]] else 0),
    entropy = unname(entropy))
}

.refmap_distance_statistic <- function(distance, primary_k) {
  stats::median(distance[seq_len(min(as.integer(primary_k),
    length(distance)))])
}

.refmap_json_bytes <- function(value) charToRaw(jsonlite::toJSON(value,
  auto_unbox = TRUE, null = "null", na = "null", digits = 17,
  dataframe = "rows", pretty = FALSE))

.refmap_read_calibration <- function(reader) {
  .refmap_assert(identical(as.character(reader$primary_path)[1L],
    "calibration.json"), "Calibration artifact primary path changed.")
  value <- jsonlite::fromJSON(rawToChar(reader$read_bytes(reader$primary_path)),
    simplifyVector = FALSE)
  .refmap_assert(identical(value$schema,
    "scwbreferencemap.reference_map_calibration.v1"),
    "Calibration artifact schema is unsupported.")
  value
}
