library(tidyverse)

scores <- read_tsv("output/label_projection/combined.extract_scores.output.tsv") %>%
  rename(metric_id = metric_ids, metric_value = metric_values)

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



df <- scores %>%
  left_join(method_info %>% select(id, type, label) %>% rename_all(function(x) paste0("method_", x)), by = "method_id") %>%
  left_join(metric_info %>% select(id, label, min, max, maximise) %>% rename_all(function(x) paste0("metric_", x)), by = "metric_id") %>%
  mutate(method_label = factor(method_label, levels = c("True labels", "Multilayer perceptron")))

ordering <- df %>%
  group_by(metric_id, dataset_id) %>%
  mutate(rank = rank(ifelse(metric_maximise, -metric_value, metric_value))) %>%
  ungroup() %>%
  group_by(method_id, method_label) %>%
  summarise(mean_rank = mean(rank)) %>%
  arrange(mean_rank)

df$method_label <- factor(df$method_label, levels = rev(ordering$method_label))

g <- ggplot(df %>% arrange(method_label)) +
  geom_path(aes(metric_value, method_label, group = dataset_id), alpha = .2) +
  geom_point(aes(metric_value, method_label, colour = method_type)) +
  facet_wrap(~metric_label, ncol = 1) +
  theme_bw() +
  labs(title = "OpenProblems v2 - Label projection v0.1")

ggsave("output/label_projection/plot.pdf", g, width = 6, height = 8)
ggsave("output/label_projection/plot.png", g, width = 6, height = 8)
