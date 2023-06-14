
anndata_struct_names <- c("obs", "var", "obsm", "obsp", "varm", "varp", "layers", "uns")

read_anndata_spec <- function(path) {
  spec <- read_and_merge_yaml(path)
  list(
    info = read_anndata_info(spec, path),
    slots = read_anndata_slots(spec, path)
  )
}
read_anndata_info <- function(spec, path) {
  # TEMP: make it readable
  spec$info$slots <- NULL
  df <- list_as_tibble(spec)
  if (list_contains_tibble(spec$info)) {
    df <- dplyr::bind_cols(df, list_as_tibble(spec$info))
  }
  df$file_name <- basename(path) %>% gsub("\\.yaml", "", .)
  as_tibble(df)
}
read_anndata_slots <- function(spec, path) {
  map_df(
    anndata_struct_names,
    function(struct_name, slot) {
      slot <- spec$info$slots[[struct_name]]
      if (is.null(slot)) return(NULL)
      df <- map_df(slot, as.data.frame)
      df$struct <- struct_name
      df$file_name <- basename(path) %>% gsub("\\.yaml", "", .)
      df$required <- df$required %||% TRUE %|% TRUE
      df$multiple <- df$multiple %||% FALSE %|% FALSE
      as_tibble(df)
    }
  )
}

format_slots <- function(spec) {
  example <- spec$slots %>%
    group_by(struct) %>%
    summarise(
      str = paste0(unique(struct), ": ", paste0("'", name, "'", collapse = ", "))
    ) %>%
    arrange(match(struct, anndata_struct_names))

  c("    AnnData object", paste0("     ", example$str))
}

format_slots_as_kable <- function(spec) {
  if (nrow(spec$slots) == 0) return("")
  spec$slots %>%
    mutate(
      tag_str = pmap_chr(lst(required), function(required) {
        out <- c()
        if (!required) {
          out <- c(out, "Optional")
        }
        if (length(out) == 0) {
          ""
        } else {
          paste0("(_", paste(out, collapse = ", "), "_) ")
        }
      })
    ) %>%
    transmute(
      Slot = paste0("`", struct, "[\"", name, "\"]`"),
      Type = paste0("`", type, "`"),
      Description = paste0(
        tag_str,
        description %>% gsub(" *\n *", " ", .) %>% gsub("\\. *$", "", .), 
        "."
      )
    ) %>%
    knitr::kable()
}

list_contains_tibble <- function(li) {
  is.list(li) && any(sapply(li, is.atomic))
}

list_as_tibble <- function(li) {
  as.data.frame(li[sapply(li, is.atomic)], check.names = FALSE)
}

read_comp_spec <- function(path) {
  spec_yaml <- read_and_merge_yaml(path)
  list(
    info = read_comp_info(spec_yaml, path),
    args = read_comp_args(spec_yaml, path)
  )
}

read_comp_info <- function(spec_yaml, path) {
  # TEMP: make it readable
  spec_yaml$functionality$arguments <- NULL
  spec_yaml$functionality$argument_groups <- NULL
  
  df <- list_as_tibble(spec_yaml$functionality)
  if (list_contains_tibble(spec_yaml$functionality$info)) {
    df <- dplyr::bind_cols(df, list_as_tibble(spec_yaml$functionality$info))
  }
  if (list_contains_tibble(spec_yaml$functionality$info$type_info)) {
    df <- dplyr::bind_cols(df, list_as_tibble(spec_yaml$functionality$info$type_info))
  }
  df$file_name <- basename(path) %>% gsub("\\.yaml", "", .)
  as_tibble(df)
}
read_comp_args <- function(spec_yaml, path) {
  arguments <- spec_yaml$functionality$arguments
  for (arg_group in spec_yaml$functionality$argument_groups) {
    arguments <- c(arguments, arg_group$arguments)
  }
  map_df(arguments, function(arg) {
    df <- list_as_tibble(arg)
    if (list_contains_tibble(arg$info)) {
      df <- dplyr::bind_cols(df, list_as_tibble(arg$info))
    }
    df$file_name <- basename(path) %>% gsub("\\.yaml", "", .)
    df$arg_name <- stringr::str_replace_all(arg$name, "^-*", "")
    df$direction <- df$direction %||% "input" %|% "input"
    df$parent <- df$`__merge__` %||% NA_character_ %>% basename() %>% gsub("\\.yaml", "", .)
    df$required <- df$required %||% FALSE %|% FALSE
    df$default <- df$default %||% NA_character_ %>% as.character
    df$example <- df$example %||% NA_character_ %>% as.character
    df$description <- df$description %||% NA_character_ %>% as.character
    df$summary <- df$summary %||% NA_character_ %>% as.character
    df
  })
}

format_comp_args_as_tibble <- function(spec) {
  if (nrow(spec$args) == 0) return("")
  spec$args %>%
    mutate(
      tag_str = pmap_chr(lst(required, direction), function(required, direction) {
        out <- c()
        if (!required) {
          out <- c(out, "Optional")
        }
        if (direction == "output") {
          out <- c(out, "Output")
        }
        if (length(out) == 0) {
          ""
        } else {
          paste0("(_", paste(out, collapse = ", "), "_) ")
        }
      })
    ) %>%
    transmute(
      Name = paste0("`--", arg_name, "`"),
      Type = paste0("`", type, "`"),
      Description = paste0(
        tag_str,
        (summary %|% description) %>% gsub(" *\n *", " ", .) %>% gsub("\\. *$", "", .), 
        ".",
        ifelse(!is.na(default), paste0(" Default: `", default, "`."), "")
      )
    ) %>%
    knitr::kable()
}

# path <- "src/datasets/api/comp_processor_knn.yaml"
render_component <- function(path) {
  spec <- read_comp_spec(path)

  cat(strip_margin(glue::glue("
    §# Component type: {spec$info$label}
    §
    §Path: [`src/{spec$info$namespace}`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/{spec$info$namespace})
    §
    §{spec$info$summary}
    §
    §Arguments:
    §
    §:::{{.small}}
    §{paste(format_comp_args_as_tibble(spec), collapse = '\n')}
    §:::
    §
    §"), symbol = "§"))
}

# path <- "src/datasets/api/file_pca.yaml"
render_file <- function(path) {
  spec <- read_anndata_spec(path)

  cat(strip_margin(glue::glue("
    §# File format: {spec$info$label}
    §
    §Example file: `{spec$info$example %|% '<Missing>'}`
    §
    §{spec$info$summary}
    §
    §Format:
    §
    §:::{{.small}}
    §{paste(format_slots(spec), collapse = '\n')}
    §:::
    §
    §Slot description:
    §
    §:::{{.small}}
    §{paste(format_slots_as_kable(spec), collapse = '\n')}
    §:::
    §
    §"), symbol = "§"))
}