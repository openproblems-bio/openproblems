
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
  df$description <- df$description %||% NA_character_ %>% as.character
  df$summary <- df$summary %||% NA_character_ %>% as.character
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
    df$arg_name <- gsub("^-*", "", arg$name)
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
render_component <- function(spec) {
  if (is.character(spec)) {
    spec <- read_comp_spec(spec)
  }

  strip_margin(glue::glue("
    §## Component type: {spec$info$label}
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
    §"), symbol = "§")
}

# path <- "src/datasets/api/file_pca.yaml"
render_file <- function(spec) {
  if (is.character(spec)) {
    spec <- read_anndata_spec(spec)
  }

  if (!"label" %in% names(spec$info)) {
    spec$info$label <- basename(spec$info$example)
  }

  example <-
    if (is.null(spec$info$example) || is.na(spec$info$example)) {
      ""
    } else {
      paste0("Example file: `", spec$info$example, "`")
    }

  description <-
    if (is.null(spec$info$description) || is.na(spec$info$description)) {
      ""
    } else {
      paste0("Description:\n\n", spec$info$description)
    }

  strip_margin(glue::glue("
    §## File format: {spec$info$label}
    §
    §{spec$info$summary %||% ''}
    §
    §{example}
    §
    §{description}
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
    §"), symbol = "§")
}

# path <- "src/tasks/denoising"
read_task_api <- function(path) {
  cli::cli_inform("Looking for project root")
  project_path <- .ram_find_project(path)
  api_dir <- paste0(path, "/api")

  cli::cli_inform("Reading task info")
  task_info_yaml <- list.files(api_dir, pattern = "task_info.ya?ml", full.names = TRUE)
  assertthat::assert_that(length(task_info_yaml) == 1)
  task_info <- read_and_merge_yaml(task_info_yaml, project_path)

  cli::cli_inform("Reading task authors")
  authors <- map_df(task_info$authors, function(aut) {
    aut$roles <- paste(aut$roles, collapse = ", ")
    list_as_tibble(aut)
  })

  cli::cli_inform("Reading component yamls")
  comp_yamls <- list.files(api_dir, pattern = "comp_.*\\.ya?ml", full.names = TRUE)
  comps <- map(comp_yamls, read_comp_spec)
  comp_info <- map_df(comps, "info")
  comp_args <- map_df(comps, "args")
  names(comps) <- basename(comp_yamls) %>% gsub("\\..*$", "", .)

  cli::cli_inform("Reading file yamls")
  file_yamls <- .ram_resolve_path(
    path = na.omit(unique(comp_args$`__merge__`)),
    project_path = project_path,
    parent_path = api_dir
  )
  files <- map(file_yamls, read_anndata_spec)
  names(files) <- basename(file_yamls) %>% gsub("\\..*$", "", .)
  file_info <- map_df(files, "info")
  file_slots <- map_df(files, "slots")

  cli::cli_inform("Generating task graph")
  task_graph <- create_task_graph(file_info, comp_info, comp_args)

  list(
    task_info = task_info,
    file_specs = files,
    file_info = file_info,
    file_slots = file_slots,
    comp_specs = comps,
    comp_info = comp_info,
    comp_args = comp_args,
    task_graph = task_graph,
    authors = authors
  )
}


create_task_graph <- function(file_info, comp_info, comp_args) {
  clean_id <- function(id) {
    gsub("graph", "graaf", id)
  }
  nodes <-
    bind_rows(
      file_info %>%
        mutate(id = file_name, label = label, is_comp = FALSE),
      comp_info %>%
        mutate(id = file_name, label = label, is_comp = TRUE)
    ) %>%
      select(id, label, everything()) %>%
      mutate(str = paste0(
        "  ",
        clean_id(id),
        ifelse(is_comp, "[/\"", "(\""),
        label,
        ifelse(is_comp, "\"/]", "\")")
      ))
  edges <- bind_rows(
    comp_args %>%
      filter(type == "file", direction == "input") %>%
      mutate(
        from = parent,
        to = file_name,
        arrow = "---"
      ),
    comp_args %>%
      filter(type == "file", direction == "output") %>%
      mutate(
        from = file_name,
        to = parent,
        arrow = "-->"
      )
  ) %>%
    select(from, to, everything()) %>%
    mutate(str = paste0("  ", clean_id(from), arrow, clean_id(to)))

  igraph::graph_from_data_frame(
    edges,
    vertices = nodes,
    directed = TRUE
  )
}

.task_graph_get_root <- function(task_api) {
  root <- names(which(igraph::degree(task_api$task_graph, mode = "in") == 0))
  if (length(root) > 1) {
    warning(
      "There should probably only be one node with in-degree equal to 0.\n",
      "  Nodes with in-degree == 0: ", paste(root, collapse = ", ")
    )
  }
  root[[1]]
}

render_task_graph <- function(task_api, root = .task_graph_get_root(task_api)) {
  order <- names(igraph::bfs(task_api$task_graph, root)$order)

  vdf <- igraph::as_data_frame(task_api$task_graph, "vertices") %>%
    arrange(match(name, order))
  edf <- igraph::as_data_frame(task_api$task_graph, "edges") %>%
    arrange(match(from, order), match(to, order))

  strip_margin(glue::glue("
    §```mermaid
    §flowchart LR
    §{paste(vdf$str, collapse = '\n')}
    §{paste(edf$str, collapse = '\n')}
    §```
    §"), symbol = "§")
}



# Recursive function to process each property with indentation
.render_example_process_property <- function(prop, prop_name = NULL, indent_level = 0) {
  if (is.null(prop_name)) {
    prop_name <- ""
  }

  out <- c()

  # define helper variables
  indent_spaces <- strrep(" ", indent_level)
  next_indent_spaces <- strrep(" ", indent_level + 2)

  # add comment if available
  if ("description" %in% names(prop)) {
    comment <- gsub("\n", paste0("\n", indent_spaces, "# "), stringr::str_trim(prop$description))
    out <- c(out, indent_spaces, "# ", comment, "\n")
  }

  # add variable
  out <- c(out, indent_spaces, prop_name, ": ")

  if (prop$type == "object" && "properties" %in% names(prop)) {
    # Handle object with properties
    prop_names <- setdiff(names(prop$properties), "additionalProperties")
    sub_props <- unlist(lapply(prop_names, function(sub_prop_name) {
      prop_out <- .render_example_process_property(
        prop$properties[[sub_prop_name]],
        sub_prop_name,
        indent_level + 2
      )
      c(prop_out, "\n")
    }))
    c(out, "\n", sub_props[-length(sub_props)])
  } else if (prop$type == "array") {
    if (is.list(prop$items) && "properties" %in% names(prop$items)) {
      # Handle array of objects
      array_items_yaml <- unlist(lapply(names(prop$items$properties), function(item_prop_name) {
        prop_out <- .render_example_process_property(
          prop$items$properties[[item_prop_name]],
          item_prop_name,
          indent_level + 4
        )
        c(prop_out, "\n")
      }))
      c(out, "\n", next_indent_spaces, "- ", array_items_yaml[-1])
    } else {
      # Handle simple array
      c(out, "[ ... ]")
    }
  } else {
    c(out, "...")
  }
}

# Function for rendering an example yaml based on a JSON schema
render_example <- function(json_schema) {
  if (!"properties" %in% names(json_schema)) {
    return("")
  }
  text <-
    unlist(lapply(names(json_schema$properties), function(prop_name) {
      out <- .render_example_process_property(
        json_schema$properties[[prop_name]],
        prop_name,
        0
      )
      c(out, "\n")
    }))

  paste(text, collapse = "")
}