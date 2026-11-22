args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE)

parameter <- function(type, default = NULL, minimum = NULL, maximum = NULL,
                      nullable = is.null(default)) {
  out <- list(type = if (nullable) c(type, "null") else type)
  if (!is.null(default)) out$default <- default
  if (!is.null(minimum)) out$minimum <- minimum
  if (!is.null(maximum)) out$maximum <- maximum
  out
}
identity <- function(name, probe, spec, package = sub("_runtime$", "", name),
                     disposition = "required") list(requirement_id = name,
  kind = "r_package", name = package,
  version = list(scheme = "r_package_version_constraint_v1",
    specifier = spec), disposition = disposition, import_probe = probe)
input <- function(name, kind, required = TRUE, reference_only = FALSE) {
  out <- list(name = name, kind = kind, required = required,
    alignment = switch(kind, matrix = "cell_and_feature_ordered", space =
      "cell_ordered", cell_vector = "cell_ordered", cell_vector_family =
      "cell_ordered", model_ref = "model_linked", obs_table = "none"))
  if (isTRUE(reference_only)) out$effective_cell_coverage <- "reference_only"
  out
}
output <- function(name, kind, group, table_type = NULL,
                   value_kind = NULL, model_policy = NULL) {
  out <- list(name = name, category = "terminal_artifact", kind = kind,
    storage_profile = if (identical(kind, "model_ref"))
      "model_ref_contained_tree_v1" else "terminal_primitive_v1",
    allocated_id_key = paste0(name, "_id"), required = TRUE,
    promotion_group = group)
  if (!is.null(table_type)) out$table_type <- table_type
  if (!is.null(value_kind)) out$value_kind <- value_kind
  if (!is.null(model_policy)) out$contained_model_sibling_closure_policy <-
    model_policy
  out
}
retained <- function(params) lapply(names(params), function(name) {
  schema <- params[[name]]
  nullable_default <- "default" %in% names(schema) &&
    is.null(schema$default) && "null" %in% as.character(unlist(schema$type,
      use.names = FALSE))
  list(name = name, source = "effective_parameter_value", parameter = name,
    required = !nullable_default)
})
full_space_read <- function(name, input, role, chunk = "space_chunk_size") list(
  name = name, input = input, reader = "read_typed_space", kind = "space",
  role = role, selection = list(policy = "full_input_axis_v1"),
  chunk_size_parameter = chunk,
  call_discipline = "ordered_partition_exactly_once")
closure <- function(names) stats::setNames(lapply(names, function(x)
  unname(as.list(sort(names, method = "radix")))), names)
replay <- function(names, model = character()) stats::setNames(lapply(names,
  function(name) if (name %in% model) "loadable_reference" else
    "exact_canonical_value"), names)
sparse_family_policy <- function(output, input) {
  scheme <- "scworkbench.operation_pack.standalone_sparse_cell_family_policy.v1"
  value <- list(schema = scheme, output = output,
    policy = "standalone_sparse_cell_family_v1",
    axis_source = list(source = "typed_space_input_v1", input = input),
    key_order = "declared_verbatim_v1",
    row_order = "key_major_then_source_relative_cell_v1",
    per_key_coverage = "per_key_source_relative_subset_v1",
    missing_values = "forbidden",
    declaration_discipline = "exactly_once_before_write",
    maximum_axis_cells = 1000000L, maximum_keys = 64L,
    maximum_rows = 1000000L, maximum_id_bytes = 1024L,
    maximum_control_bytes = 67108864L)
  value$policy_identity <- getFromNamespace(
    ".wb_operation_pack_v2_identity_scoped", "scworkbench")(
      scheme, "policy excluding policy_identity", value)
  value
}

calibration_params <- list(
  holdout_group_field = list(type = c("string", "null"), default = "orig.ident",
    minLength = 1L),
  k_values = list(type = "array", default = list(15L, 30L, 50L),
    minItems = 1L, uniqueItems = TRUE,
    items = list(type = "integer", minimum = 1L, maximum = 10000L)),
  primary_k = parameter("integer", 30L, 1L, 10000L, FALSE),
  distance_statistic = list(type = "string",
    default = "median_first_primary_k_distance",
    enum = list("median_first_primary_k_distance")),
  support_quantiles = list(type = "array", default = list(.95, .99),
    minItems = 2L, maxItems = 2L, uniqueItems = TRUE,
    items = list(type = "number", exclusiveMinimum = 0,
      exclusiveMaximum = 1)),
  primary_continuous_key = list(type = c("string", "null"), default = NULL,
    minLength = 1L),
  space_chunk_size = parameter("integer", 1024L, 1L, 1000000L, FALSE),
  query_chunk_size = parameter("integer", 1024L, 1L, 1000000L, FALSE))
calibration_outputs <- c("calibration", "heldout_predictions",
  "calibration_metrics", "thresholds", "coverage_limitations")
model_policy <- list(policy = "staged_contained_model_sibling_closure_v1",
  engine = "scwbreferencemap", model_type = "reference_map_calibration",
  primary_path = "calibration.json", source_space_input = "reference_space",
  dimensions_parameter = NULL, required_sibling_terminals = unname(as.list(
    sort(setdiff(calibration_outputs, "calibration"), method = "radix"))),
  file_order = "utf8_path_lexicographic_v1", maximum_file_count = 4L,
  maximum_total_bytes = 67108864L)
calibrate <- list(
  name = "reference_map_calibrate",
  qualified_name = "scwbreferencemap::reference_map_calibrate",
  title = "Calibrate a reusable reference map",
  description = paste("Leave one configured metadata category out at a time",
    "to calibrate reference-space support and optional transfer behavior."),
  runner = ".run_reference_map_calibrate",
  inputs = list(input("reference_matrix", "matrix"),
    input("reference_space", "space"),
    input("reference_annotation", "cell_vector", FALSE),
    input("reference_continuous", "cell_vector_family", FALSE)),
  params_schema = list(type = "object", properties = calibration_params,
    additionalProperties = FALSE),
  registry_reads = list(list(name = "holdout_metadata",
    category = "authenticated_core_registry_snapshot",
    kind = "ordered_cell_metadata", reader = "read_ordered_cell_metadata",
    matrix_input = "reference_matrix", field_parameters = list(list(
      parameter = "holdout_group_field", role = "categorical",
      cardinality = "zero_or_one")),
    row_order = "matrix_effective_cell_order_v1",
    missing_value_policy = "preserve_null_v1", call_discipline = "single_use")),
  typed_space_reads = list(full_space_read("reference_coordinates",
    "reference_space", "candidate")),
  dimension_policies = lapply(setdiff(calibration_outputs, "calibration"),
    function(name) list(output = name, policy = "data_dependent_table_v1")),
  outputs = list(output("calibration", "model_ref", "calibration",
      model_policy = model_policy),
    output("heldout_predictions", "obs_table", "calibration",
      "reference_map_heldout_predictions"),
    output("calibration_metrics", "obs_table", "calibration",
      "reference_map_calibration_metrics"),
    output("thresholds", "obs_table", "calibration",
      "reference_map_thresholds"),
    output("coverage_limitations", "obs_table", "calibration",
      "reference_map_coverage_limitations")),
  dependencies = list(identity("scworkbench_runtime", "wb_prim_model_ref",
      ">=0.1.0,<0.2.0"), identity("jsonlite_runtime", "toJSON",
      ">=1.8.0,<3.0.0", "jsonlite"), identity("rann_runtime", "nn2",
      ">=2.6.0,<3.0.0", "RANN")),
  promotion_closure = closure(calibration_outputs),
  replay_equivalence = replay(calibration_outputs, "calibration"),
  retained_parameters = retained(calibration_params),
  public_interface = list(intent = list("reference", "mapping", "calibration"),
    inputs = list(
      list(name = "reference_matrix", label = "reference_matrix",
        help = "Matrix supplying the authenticated reference metadata axis."),
      list(name = "reference_space", label = "reference_space",
        help = "Authenticated reference coordinates."),
      list(name = "reference_annotation", label = "reference_annotation",
        help = "Optional categorical annotation to calibrate."),
      list(name = "reference_continuous", label = "reference_continuous",
        help = "Optional continuous vectors, including selected pseudotime.")),
    parameters = list(
      list(name = "holdout_group_field", type = "string",
        help = "Categorical metadata field used for held-out categories."),
      list(name = "k_values", type = "list", help = "Neighbor k values."),
      list(name = "primary_k", type = "int", help = "Primary neighbor k."),
      list(name = "distance_statistic", type = "choice",
        help = "Support distance statistic."),
      list(name = "support_quantiles", type = "list",
        help = "Supported and borderline calibration quantiles."),
      list(name = "primary_continuous_key", type = "string",
        help = "Selected primary trajectory or other continuous key."),
      list(name = "space_chunk_size", type = "int",
        help = "Bounded coordinate read size."),
      list(name = "query_chunk_size", type = "int",
        help = "Bounded exact-neighbor query size.")),
    outputs = lapply(calibration_outputs, function(name) list(name = name,
      description = paste("Reference-map calibration", name, "artifact.")))))

transfer_params <- list(
  dims = list(type = c("array", "null"), default = NULL,
    uniqueItems = TRUE, items = list(type = "integer", minimum = 1L)),
  supported_threshold = parameter("number", nullable = TRUE, minimum = 0),
  out_of_reference_threshold = parameter("number", nullable = TRUE,
    minimum = 0),
  support_output_key = list(type = "string", default = "reference_support",
    minLength = 1L),
  annotation_output_key = list(type = "string",
    default = "reference_annotation", minLength = 1L),
  continuous_output_keys = list(type = c("array", "null"), default = NULL,
    uniqueItems = TRUE, items = list(type = "string", minLength = 1L)),
  continuous_acceptance = list(type = "object",
    default = list(mode = "support_only"), additionalProperties = FALSE,
    properties = list(
      mode = list(type = "string", default = "support_only",
        enum = list("support_only", "support_and_annotation")),
      accepted_annotation_values = list(type = c("array", "string", "null"),
        default = NULL, minItems = 1L, uniqueItems = TRUE,
        items = list(type = "string", minLength = 1L)),
      minimum_annotation_confidence = list(type = "number", default = .5,
        minimum = 0, maximum = 1))),
  space_chunk_size = parameter("integer", 1024L, 1L, 1000000L, FALSE))
transfer_outputs <- c("support_call", "accepted_annotation",
  "numeric_outputs", "per_query_summary", "categorical_estimates",
  "continuous_estimates", "continuous_acceptance", "grouped_summary")
transfer <- list(
  name = "reference_map_transfer",
  qualified_name = "scwbreferencemap::reference_map_transfer",
  title = "Transfer a frozen reference map",
  description = paste("Authenticate a query-to-reference ledger against a",
    "frozen calibration, classify support, and transfer optional values."),
  runner = ".run_reference_map_transfer",
  inputs = list(input("in_space", "space"),
    input("reference_space", "space"),
    input("neighbors", "obs_table", reference_only = TRUE),
    input("calibration", "model_ref", reference_only = TRUE)),
  params_schema = list(type = "object", properties = transfer_params,
    additionalProperties = FALSE),
  typed_space_reads = list(full_space_read("query_coordinates", "in_space",
      "query"), full_space_read("reference_coordinates", "reference_space",
      "candidate")),
  self_linked_contained_model_reads = list(list(
    name = "calibration_payload", input = "calibration",
    reader = "read_self_linked_contained_model_payload",
    allowed_models = list(list(engine = "scwbreferencemap",
      model_type = "reference_map_calibration")),
    linkage_policy = "authenticated_existing_source_and_cell_linkage_v1",
    maximum_file_count = 4L, maximum_total_bytes = 67108864L,
    call_discipline = "exactly_once")),
  dimension_policies = c(
    list(list(output = "support_call", policy = "authenticated_output_axes_v1",
      shape = "one_per_cell", cell_axis = list(source =
        "effective_cells_v1", input = "in_space"),
      row_order = "input_cell_order_v1", missing_values = "forbidden"),
      list(output = "accepted_annotation", policy = "data_dependent_table_v1"),
      sparse_family_policy("numeric_outputs", "in_space")),
    lapply(setdiff(transfer_outputs,
      c("support_call", "accepted_annotation", "numeric_outputs")),
      function(name) list(output = name, policy = "data_dependent_table_v1"))),
  outputs = list(
    output("support_call", "cell_vector", "transfer", value_kind = "categorical"),
    output("accepted_annotation", "obs_table", "transfer",
      "reference_map_accepted_annotation"),
    output("numeric_outputs", "cell_vector_family", "transfer",
      value_kind = "continuous"),
    output("per_query_summary", "obs_table", "transfer",
      "reference_map_per_query"),
    output("categorical_estimates", "obs_table", "transfer",
      "reference_map_categorical_estimates"),
    output("continuous_estimates", "obs_table", "transfer",
      "reference_map_continuous_estimates"),
    output("continuous_acceptance", "obs_table", "transfer",
      "reference_map_continuous_acceptance"),
    output("grouped_summary", "obs_table", "transfer",
      "reference_map_grouped_summary")),
  dependencies = list(identity("scworkbench_runtime", "wb_prim_obs_table",
      ">=0.1.0,<0.2.0"), identity("jsonlite_runtime", "fromJSON",
      ">=1.8.0,<3.0.0", "jsonlite")),
  promotion_closure = closure(transfer_outputs),
  replay_equivalence = replay(transfer_outputs),
  retained_parameters = retained(transfer_params),
  public_interface = list(intent = list("reference", "mapping", "transfer"),
    inputs = list(
      list(name = "in_space", label = "query_space",
        help = "Authenticated query coordinates."),
      list(name = "reference_space", label = "reference_space",
        help = "Authenticated calibrated reference coordinates."),
      list(name = "neighbors", label = "neighbors",
        help = "Authenticated query-to-reference neighbor ledger."),
      list(name = "calibration", label = "calibration",
        help = "Frozen retained reference-map calibration.")),
    parameters = list(
      list(name = "dims", type = "list", help = "Calibrated dimensions."),
      list(name = "supported_threshold", type = "dbl",
        help = "Optional transfer-only supported threshold override."),
      list(name = "out_of_reference_threshold", type = "dbl",
        help = "Optional transfer-only out-of-reference threshold override."),
      list(name = "support_output_key", type = "string",
        help = "Caller-supplied support output key."),
      list(name = "annotation_output_key", type = "string",
        help = "Caller-supplied annotation output key."),
      list(name = "continuous_output_keys", type = "string_vec",
        help = "Caller-supplied keys for continuous outputs."),
      list(name = "continuous_acceptance", type = "list",
        help = paste("Policy for accepting continuous outputs after reference",
          "support and optional annotation-membership checks.")),
      list(name = "space_chunk_size", type = "int",
        help = "Bounded coordinate read size.")),
    outputs = lapply(transfer_outputs, function(name) list(name = name,
      description = paste("Reference-map transfer", name, "artifact.")))))

manifest <- list(
  schema_version = "scworkbench.operation_pack.v3", pack_abi_version = "3",
  identity_scheme = "scworkbench.operation_pack.canonical_json.v1",
  namespace = "scwbreferencemap",
  package = list(name = "scwbReferenceMap", version = "0.2.0"),
  backend = "R", compatibility = list(
    scworkbench = list(minimum = "0.1.2", maximum_exclusive = "0.2.0"),
    pack_abi = list(minimum = "3", maximum = "3")),
  entry_point = "scwbReferenceMap::register_operation_pack",
  operations = list(calibrate, transfer))
jsonlite::write_json(manifest,
  file.path(root, "inst", "scworkbench", "operation-pack.json"),
  auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
