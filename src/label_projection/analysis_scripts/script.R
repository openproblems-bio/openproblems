library(tidyverse)

out_dir <- "resources/label_projection/benchmarks/openproblems_v1/"

# read scores
scores <- read_tsv(paste0(out_dir, "combined.extract_scores.output.tsv")) %>%
  rename(metric_id = metric_ids, metric_value = metric_values)

# read nxf log
output_regex <- "^(.*)\\.([^\\.]*)\\.output[^\\.]*\\.h5ad"
nxf_log <- read_tsv(paste0(out_dir, "nextflow_log.tsv")) %>%
  mutate(
    output_file = gsub(".*output=([^,]*).*", "\\1", params),
    id = gsub(output_regex, "\\1", output_file),
    dataset_id = gsub("^([^\\.]*)\\.([^\\.]*).*", "\\1/\\2", id),
    component_id = gsub(output_regex, "\\2", output_file)
  )
nxf_log %>% select(id:component_id)

# process execution info
execution_info <- nxf_log %>%
  filter(component_id %in% method_info$id) %>%
  transmute(
    method_id = component_id,
    dataset_id,
    status,
    realtime = lubridate::duration(toupper(realtime)),
    pcpu = as.numeric(gsub("%", "", pcpu)),
    vmem_gb = as.numeric(gsub(" GB", "", vmem)),
    peak_vmem_gb = as.numeric(gsub(" GB", "", peak_vmem)),
    read_bytes_mb = as.numeric(gsub(" MB", "", read_bytes)),
    write_bytes_mb = as.numeric(gsub(" MB", "", write_bytes))
  )

# get method info
ns_list_methods <- yaml::yaml.load(processx::run("viash", c("ns", "list", "-q", "label_projection.*methods"))$stdout)

method_info <- map_df(ns_list_methods, function(conf) {
  tryCatch({
    info <- c(
      list(
        id = conf$functionality$name,
        namespace = conf$functionality$namespace,
        description = conf$functionality$description
      ),
      conf$functionality$info
    )
    as.data.frame(info)
  }, error = function(err) {
    cat(err$message, "\n", sep = "")
    data.frame(id = conf$functionality$name)
  })
})

# get metric info
ns_list_metrics <- yaml::yaml.load(processx::run("viash", c("ns", "list", "-q", "label_projection.*metrics"))$stdout)

metric_info <- map_df(ns_list_metrics, function(conf) {
  tryCatch({
    map_df(conf$functionality$info$metrics, as.data.frame)
  }, error = function(err) {
    cat(err$message, "\n", sep = "")
    data.frame(id = conf$functionality$name)
  })
})

# get data table
ranking <- scores %>%
  left_join(metric_info %>% select(metric_id = id, maximise), by = "metric_id") %>%
  inner_join(method_info %>% select(method_id = id, method_label = label), by = "method_id") %>%
  group_by(metric_id, dataset_id) %>%
  mutate(rank = rank(ifelse(maximise, -metric_value, metric_value))) %>%
  ungroup() %>%
  group_by(method_id, method_label) %>%
  summarise(mean_rank = mean(rank), .groups = "drop") %>%
  arrange(mean_rank)

df <-
  method_info %>%
  select(id, type, label) %>%
  rename_all(function(x) paste0("method_", x)) %>%
  inner_join(scores %>% spread(metric_id, metric_value), by = "method_id") %>%
  left_join(execution_info, by = c("dataset_id", "method_id")) %>%
  mutate(method_label = factor(method_label, levels = rev(ranking$method_label))) %>%
  arrange(method_label)


# get feature info
feature_info_exec <- tribble(
  ~id, ~label, ~log_x,
  "realtime", "Duration (s)", FALSE,
  "pcpu", "CPU usage (%)", FALSE,
  "vmem_gb", "Memory usage (GB)", FALSE,
  "peak_vmem_gb", "Peak memory (GB)", FALSE,
  "read_bytes_mb", "Read disk (MB)", FALSE,
  "write_bytes_mb", "Write disk (MB)", FALSE
)
feature_info <- bind_rows(
  metric_info %>% transmute(id, label, log_x = FALSE),
  feature_info_exec
)

plots <- pmap(feature_info, function(id, label, log_x) {
  g <- ggplot(df) +
    geom_path(aes_string(id, "method_label", group = "dataset_id"), alpha = .2) +
    geom_point(aes_string(id, "method_label", colour = "method_type")) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = label, y = NULL) +
    expand_limits(x = 0)
  if (log_x) {
    g <- g + scale_x_log10()
  }
  g
})
g <- patchwork::wrap_plots(plots, ncol = 2, byrow = FALSE)


ggsave(paste0(out_dir, "plot.pdf"), g, width = 8, height = 12)
ggsave(paste0(out_dir, "plot.png"), g, width = 8, height = 12)


