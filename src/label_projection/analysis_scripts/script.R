library(tidyverse)

out_dir <- "resources/label_projection/benchmarks/openproblems_v1/"

scores <- read_tsv(paste0(out_dir, "combined.extract_scores.output.tsv")) %>%
  rename(metric_id = metric_ids, metric_value = metric_values)

output_regex <- "^(.*)\\.([^\\.]*)\\.output[^\\.]*\\.h5ad"
nxf_log <- read_tsv(paste0(out_dir, "nextflow_log.tsv")) %>%
  mutate(
    output_file = gsub(".*output=([^,]*).*", "\\1", params),
    id = gsub(output_regex, "\\1", output_file),
    dataset_id = gsub("^([^\\.]*)\\.([^\\.]*).*", "\\1/\\2", id),
    component_id = gsub(output_regex, "\\2", output_file)
  )
nxf_log %>% select(id:component_id)

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

ns_list_metrics <- yaml::yaml.load(processx::run("viash", c("ns", "list", "-q", "label_projection.*metrics"))$stdout)

metric_info <- map_df(ns_list_metrics, function(conf) {
  tryCatch({
    map_df(conf$functionality$info$metrics, as.data.frame)
  }, error = function(err) {
    cat(err$message, "\n", sep = "")
    data.frame(id = conf$functionality$name)
  })
})

# make scores plot
df <- scores %>%
  left_join(method_info %>% select(id, type, label) %>% rename_all(function(x) paste0("method_", x)), by = "method_id") %>%
  left_join(metric_info %>% select(id, label, min, max, maximise) %>% rename_all(function(x) paste0("metric_", x)), by = "metric_id")

ordering <- df %>%
  group_by(metric_id, dataset_id) %>%
  mutate(rank = rank(ifelse(metric_maximise, -metric_value, metric_value))) %>%
  ungroup() %>%
  group_by(method_id, method_label) %>%
  summarise(mean_rank = mean(rank)) %>%
  arrange(mean_rank)

df$method_label <- factor(df$method_label, levels = rev(ordering$method_label))

g1 <- ggplot(df %>% arrange(method_label)) +
  geom_path(aes(metric_value, method_label, group = dataset_id), alpha = .2) +
  geom_point(aes(metric_value, method_label, colour = method_type)) +
  facet_wrap(~metric_label, ncol = 1) +
  theme_bw() +
  labs(title = "OpenProblems v2 - Label projection v0.1")

# make execution plot
# todo: capture walltime, memory, disk, cpu
execution_info <- nxf_log %>%
  filter(component_id %in% method_info$id) %>%
  transmute(method_id = component_id, dataset_id, exit, duration = lubridate::duration(toupper(duration)), realtime = lubridate::duration(toupper(realtime))) %>%
  left_join(method_info %>% select(id, type, label) %>% rename_all(function(x) paste0("method_", x)), by = "method_id")

execution_info$method_label <- factor(execution_info$method_label, levels = rev(ordering$method_label))

g2 <-
  ggplot(execution_info %>% arrange(method_label)) +
  geom_path(aes(duration, method_label, group = dataset_id), alpha = .2) +
  geom_point(aes(duration, method_label, colour = method_type)) +
  theme_bw() +
  labs(title = "Execution time", x = "Log duration (s)") +
  scale_x_log10()

g <- patchwork::wrap_plots(g1, g2, ncol = 1, heights = c(4, 1))

ggsave(paste0(out_dir, "plot.pdf"), g, width = 6, height = 10)
ggsave(paste0(out_dir, "plot.png"), g, width = 6, height = 10)


