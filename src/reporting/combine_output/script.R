## VIASH START
processed_dir <- "resources_test/openproblems/task_results_v4/processed"

par <- list(
  # Inputs
  input_task_info = paste0(processed_dir, "/task_info.json"),
  input_quality_control = paste0(processed_dir, "/quality_control.json"),
  input_metric_info = paste0(processed_dir, "/metric_info.json"),
  input_method_info = paste0(processed_dir, "/method_info.json"),
  input_dataset_info = paste0(processed_dir, "/dataset_info.json"),
  input_results = paste0(processed_dir, "/results.json")
  # Outputs
  output = "task_results.json"
)
## VIASH END

################################################################################
#                              MAIN SCRIPT
################################################################################

cat("====== Combine output ======\n")

cat("\n>>> Reading input files...\n")
cat("Reading task info from '", par$input_task_info, "'...\n", sep = "")
task_info <- jsonlite::read_json(par$input_task_info)

cat("Reading quality control from '", par$input_quality_control, "'...\n", sep = "")
quality_control <- jsonlite::read_json(par$input_quality_control)

cat("Reading metric info from '", par$input_metric_info, "'...\n", sep = "")
metric_info <- jsonlite::read_json(par$input_metric_info)

cat("Reading method info from '", par$input_method_info, "'...\n", sep = "")
method_info <- jsonlite::read_json(par$input_method_info)

cat("Reading dataset info from '", par$input_dataset_info, "'...\n", sep = "")
dataset_info <- jsonlite::read_json(par$input_dataset_info)

cat("Reading results from '", par$input_results, "'...\n", sep = "")
results <- jsonlite::read_json(par$input_results)

cat("\n>>> Combining outputs...\n")
# Create combined output according to task_results_schema.json
combined_output <- list(
  task_info = task_info,
  dataset_info = dataset_info,
  method_info = method_info,
  metric_info = metric_info,
  results = results,
  quality_control = quality_control
)

cat("\n>>> Writing output file...\n")
cat("Writing combined output to '", par$output, "'...\n", sep = "")
jsonlite::write_json(
  combined_output,
  par$output,
  pretty = TRUE,
  null = "null",
  na = "null",
  auto_unbox = TRUE
)

cat("\n>>> Done!\n")
