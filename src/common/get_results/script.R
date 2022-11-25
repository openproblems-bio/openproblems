library(tidyverse)
library(rlang)

## VIASH START
par <- list(
  input_scores = "resources/label_projection/benchmarks/openproblems_v1/combined.extract_scores.output.tsv",
  input_execution = "resources/label_projection/benchmarks/openproblems_v1/nextflow_log.tsv",
  output = "resources/label_projection/output/results.yaml"
)
## VIASH END

# read scores
raw_scores <- read_tsv(par$input_scores) %>%
  spread(metric_ids, metric_values)

# read nxf log
id_regex <- "^(.*)\\.([^\\.]*)"
nxf_log <- read_tsv(par$input_execution) %>%
  mutate(
    id = tag,
    dataset_id = gsub("^([^\\.]*)\\.([^\\.]*).*", "\\1/\\2", id),
    method_id = gsub(".*\\.", "", id)
  )

# process execution info
execution_info <- nxf_log %>%
  transmute(
    method_id,
    dataset_id,
    status,
    realtime = lubridate::duration(toupper(realtime)),
    pcpu = as.numeric(gsub("%", "", pcpu)),
    vmem_gb = as.numeric(gsub(" GB", "", vmem)),
    peak_vmem_gb = as.numeric(gsub(" GB", "", peak_vmem)),
    read_bytes_mb = as.numeric(gsub(" MB", "", read_bytes)),
    write_bytes_mb = as.numeric(gsub(" MB", "", write_bytes))
  )

df <- full_join(raw_scores, execution_info, by = c("method_id", "dataset_id"))

yaml::write_yaml(purrr::transpose(df), par$output)

