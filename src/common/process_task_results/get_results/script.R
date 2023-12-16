library(tidyverse)
library(rlang)

## VIASH START
par <- list(
  input_scores = "resources/batch_integration/results/scores.yaml",
  input_execution = "resources/batch_integration/results/trace.txt",
  output = "test.json"
)
## VIASH END

# read scores
raw_scores <- yaml::yaml.load_file(par$input_scores)
score_df <- as_tibble(map_df(raw_scores, as.data.frame))

scores <- score_df %>%
  complete(dataset_id, method_id, metric_ids, fill = list(metric_values = NA_real_), normalization_id) %>%
  group_by(metric_ids, dataset_id) %>%
  mutate(
    scaled_score = dynutils::scale_minmax(metric_values) %|% 0
  ) %>%
  group_by(dataset_id, method_id) %>%
  summarise(
    metric_values = list(as.list(setNames(metric_values, metric_ids))),
    scaled_scores = list(as.list(setNames(scaled_score, metric_ids))),
    mean_score = mean(scaled_score),
    .groups = "drop"
  )

# read nxf log
nxf_log <- read_tsv(par$input_execution) %>%
    mutate(
    id = name,
    process_id = gsub(".*:(.*)_process.*", "\\1", id),
    method_id = gsub(".*\\.([^)]*)\\)", "\\1", id)
  ) %>%
  filter(process_id == method_id)


# process execution info
execution_info <- nxf_log %>%
  rowwise() %>%
  transmute(
    dataset_id = gsub(".*\\(([^/]*)\\/.*", "\\1", id),
    normalization_id = gsub(".*\\/([^.]*)\\..*", "\\1", id),
    method_id,
    resources = list(list(
      exit_code = exit,
      duration_sec = as.numeric(lubridate::duration(toupper(realtime))),
      cpu_pct = as.numeric(gsub("%", "", `%cpu`)),
      peak_memory_mb = as.numeric(gsub(" *GB", "", peak_vmem)) * 1024,
      disk_read_mb = as.numeric(gsub(" *MB", "", rchar)),
      disk_write_mb = as.numeric(gsub(" *MB", "", wchar))
    ))
  ) %>%
  ungroup()

df <- full_join(scores, execution_info, by = c("method_id", "dataset_id")) %>% 
  filter(!is.na(mean_score))

jsonlite::write_json(
  purrr::transpose(df),
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)