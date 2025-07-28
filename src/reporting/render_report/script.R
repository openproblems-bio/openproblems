## VIASH START
processed_dir <- "resources_test/openproblems/task_results_v4/processed"

par <- list(
  # Inputs
  input_task_results = paste0(processed_dir, "/task_info.json"),
  # Outputs
  output = "report.html"
)
## VIASH END

################################################################################
#                              MAIN SCRIPT
################################################################################

cat("====== Render report ======\n")

cat("\n>>> Copying input file...\n")
file.copy(
  par$input_task_results,
  file.path(meta$resources_dir, "task_results.json"),
  overwrite = TRUE
)

cat("\n>>> Rendering ...\n")
cat("Quarto version: ", as.character(quarto::quarto_version()), sep = "")
xfun::in_dir(
  meta$resources_dir,
  quarto::quarto_render(
    input = "report.qmd",
    output_file = "report.html",
    execute_params = list(
      task_results_json = "task_results.json",
      logo = "logo.svg",
      functions = "functions.R"
    )
  )
)

cat("\n>>> Copying output file...\n")
file.copy(
  file.path(meta$resources_dir, "report.html"),
  par$output,
  overwrite = TRUE
)

cat("\n>>> Done!\n")
