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

cat("\n>>> Copying input file to temporary directory...\n")
tmp_dir <- file.path(tempdir(), "render-report")
dir.create(tmp_dir, recursive = TRUE)
cat("Temporary directory: ", tmp_dir, "\n", sep = "")
file.copy(
  par$input_task_results,
  file.path(tmp_dir, "task_results.json"),
  overwrite = TRUE
)

cat("\n>>> Copying resources to temporary directory...\n")
cat("Copying 'report.qmd'...\n")
file.copy(
  file.path(meta$resources_dir, "report.qmd"),
  tmp_dir,
  overwrite = TRUE
)
cat("Copying 'logo.svg'...\n")
file.copy(
  file.path(meta$resources_dir, "logo.svg"),
  tmp_dir,
  overwrite = TRUE
)
cat("Copying 'functions.R'...\n")
file.copy(
  file.path(meta$resources_dir, "functions.R"),
  tmp_dir,
  overwrite = TRUE
)

cat("\n>>> Rendering report...\n")
cat("Quarto version: ", as.character(quarto::quarto_version()), sep = "")
xfun::in_dir(
  tmp_dir,
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
  file.path(tmp_dir, "report.html"),
  par$output,
  overwrite = TRUE
)

cat("\n>>> Done!\n")
