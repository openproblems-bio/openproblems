requireNamespace("assertthat", quietly = TRUE)

## VIASH START
## VIASH END

input <- system.file("extdata", "example_project", "api", package = "openproblems.docs")

output_path <- "output.md"

cat(">> Running the script as test\n")
out <- processx::run(
  meta[["executable"]],
  args = c("--input", input, "--output", output_path)
)

cat(">> Checking whether output files exist\n")
assertthat::assert_that(file.exists(output_path))

cat(">> Checking file contents\n")
lines <- readLines(output_path)
assertthat::assert_that(any(grepl("# Template", lines)))
assertthat::assert_that(any(grepl("## Description", lines)))
# assertthat::assert_that(any(grepl("## Authors", lines)))
assertthat::assert_that(any(grepl("flowchart TB", lines)))
assertthat::assert_that(any(grepl("## File format:", lines)))

cat("All checks succeeded!\n")
