requireNamespace("assertthat", quietly = TRUE)

## VIASH START
## VIASH END

opv2 <- paste0(meta$resources_dir, "/openproblems-v2")
output_path <- "output.md"

cat(">> Running the script as test\n")
system(paste(
  meta["executable"],
  "--task", "label_projection",
  "--output", output_path,
  "--task_dir", paste0(opv2, "/src/tasks/label_projection"),
  "--viash_yaml", paste0(opv2, "/_viash.yaml")
))

cat(">> Checking whether output files exist\n")
assertthat::assert_that(file.exists(output_path))

cat(">> Checking file contents\n")
lines <- readLines(output_path)
assertthat::assert_that(any(grepl("# Label projection", lines)))
assertthat::assert_that(any(grepl("# Description", lines)))
assertthat::assert_that(any(grepl("# Motivation", lines)))
assertthat::assert_that(any(grepl("# Authors", lines)))
assertthat::assert_that(any(grepl("flowchart LR", lines)))
assertthat::assert_that(any(grepl("# File format:", lines)))

cat("All checks succeeded!\n")
