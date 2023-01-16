library(tidyverse)

out_dir <- "resources/label_projection/benchmarks/openproblems_v1/"

# read results
results <- map_df(
  yaml::read_yaml("resources/label_projection/output/results.yaml"),
  as_tibble
)

# get method info
method_info <- map_df(
  yaml::read_yaml("resources/label_projection/output/method_info.yaml"),
  as_tibble
)

# get metric info
metric_info <- map_df(
  yaml::read_yaml("resources/label_projection/output/metric_info.yaml"),
  as_tibble
)



# get data table
ranking <- scores %>%
  left_join(metric_info %>% select(metric_id = id, maximize), by = "metric_id") %>%
  inner_join(method_info %>% select(method_id = id, method_label = label), by = "method_id") %>%
  group_by(metric_id, dataset_id) %>%
  mutate(rank = rank(ifelse(maximize, -metric_value, metric_value))) %>%
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


