### VIASH START
par <- list(
  input_task_info = "test-results/label_projection/processed/task_info.json",
  input_quality_control = "test-results/label_projection/processed/quality_control.json",
  input_metric_info = "test-results/label_projection/processed/metric_info.json",
  input_method_info = "test-results/label_projection/processed/method_info.json",
  input_dataset_info = "test-results/label_projection/processed/dataset_info.json",
  input_results = "test-results/label_projection/processed/results.json",
  input_metric_executions = "test-results/label_projection/processed/metric_execution_info.json",
  output = "tasks.json"
)
### VIASH END

task_info <- jsonlite::read_json(par$input_task_info)
quality_control <- jsonlite::read_json(par$input_quality_control)
metric_info <- jsonlite::read_json(par$input_metric_info)
method_info <- jsonlite::read_json(par$input_method_info)
dataset_info <- jsonlite::read_json(par$input_dataset_info)
results <- jsonlite::read_json(par$input_results)
metric_executions <- jsonlite::read_json(par$input_metric_executions)

task_id_clean <- stringr::str_remove(task_info$task_id, "^task_")

task <- list(
  task_id = task_info$task_id,
  commit_sha = task_info$commit_sha, # Missing from new output
  task_name = task_info$task_name,
  task_summary = task_info$task_summary,
  task_description = task_info$task_description,
  repo_url = task_info$repo,
  version = task_info$version, # build_main in new output
  license = task_info$license,
  authors = task_info$authors,
  quality_control = quality_control,
  metrics = metric_info,
  methods = method_info,
  datasets = dataset_info,
  results = results,
  metric_executions = metric_executions,
  image_url = paste0("/images/benchmarks/", task_id_clean, ".svg"),
  is_prerelease = FALSE
)

jsonlite::write_json(
  list(task),
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)
