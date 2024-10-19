requireNamespace("openproblems.docs", quietly = TRUE)
requireNamespace("processx", quietly = TRUE)

## VIASH START
par <- list(
  "input" = "path/to/input",
  "output" = "path/to/input/README.md"
)
## VIASH END

cat("Read task metadata\n")
metadata <- openproblems.docs::read_task_metadata(par$input)

cat("Render README.qmd content\n")
qmd_content <- openproblems.docs::render_task_readme_qmd(metadata)

cat("Write README.qmd to file\n")
if (!dir.exists(meta$temp_dir)) {
  dir.create(meta$temp_dir, recursive = TRUE)
}
qmd_file <- tempfile(
  pattern = "README_",
  fileext = ".qmd",
  tmpdir = meta$temp_dir
)
writeLines(qmd_content, qmd_file)

cat("Render README.qmd to README.md\n")
out <- processx::run(
  command = "quarto",
  args = c("render", qmd_file, "--output", "-"),
  echo = TRUE
)

writeLines(out$stdout, par$output)
