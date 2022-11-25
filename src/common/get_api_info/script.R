library(tidyverse)
library(rlang)

## VIASH START
par <- list(
  input = "src/label_projection",
  output = "resources/label_projection/output/api.yaml"
)
## VIASH END

comp_yamls <- list.files(paste0(par$input, "/api"), pattern = "comp_", full.names = TRUE)
file_yamls <- list.files(paste0(par$input, "/api"), pattern = "anndata_", full.names = TRUE)

# list component - file args links
comp_file <- map_df(comp_yamls, function(yaml_file) {
  conf <- yaml::read_yaml(yaml_file)

  map_df(conf$functionality$arguments, function(arg) {
    tibble(
      comp_name = basename(yaml_file) %>% gsub("\\.yaml", "", .),
      arg_name = str_replace_all(arg$name, "^-*", ""),
      direction = arg$direction %||% "input",
      file_name = basename(arg$`__inherits__`) %>% gsub("\\.yaml", "", .)
    )
  })
})

# get component info
comp_info <- map_df(comp_yamls, function(yaml_file) {
  conf <- yaml::read_yaml(yaml_file)

  tibble(
    name = basename(yaml_file) %>% gsub("\\.yaml", "", .),
    label = name %>% gsub("comp_", "", .) %>% gsub("_", " ", .)
  )
})

# get file info
file_info <- map_df(file_yamls, function(yaml_file) {
  arg <- yaml::read_yaml(yaml_file)
  
  tibble(
    name = basename(yaml_file) %>% gsub("\\.yaml", "", .),
    description = arg$description,
    short_description = arg$info$short_description,
    example = arg$example,
    label = name %>% gsub("anndata_", "", .) %>% gsub("_", " ", .)
  )
})

# get file - slot args
file_slot <- map_df(file_yamls, function(yaml_file) {
  arg <- yaml::read_yaml(yaml_file)

  map2_df(names(arg$info$slots), arg$info$slots, function(group_name, slot) {
    df <- map_df(slot, as.data.frame)
    df$struct <- group_name
    df$file_name <- basename(yaml_file) %>% gsub("\\.yaml", "", .)
    as_tibble(df)
  })
}) %>% 
  mutate(multiple = multiple %|% FALSE)

out <- list(
  comp_info = purrr::transpose(comp_info),
  file_info = purrr::transpose(file_info),
  comp_file_io = purrr::transpose(comp_file),
  file_schema = purrr::transpose(file_slot)
)

yaml::write_yaml(out, par$output)
