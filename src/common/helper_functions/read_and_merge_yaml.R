#' Read a Viash YAML
#'
#' If the YAML contains a "__merge__" key anywhere in the yaml,
#' the path specified in that YAML will be read and the two
#' lists will be merged. This is a recursive procedure.
#'
#' @param path Path to Viash YAML
read_and_merge_yaml <- function(path, project_path = .ram_find_project(path)) {
  path <- normalizePath(path, mustWork = FALSE)
  data <- tryCatch({
    suppressWarnings(yaml::read_yaml(path))
  }, error = function(e) {
    stop("Could not read ", path, ". Error: ", e)
  })
  .ram_process_merge(data, data, path, project_path)
}

.ram_find_project <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  check <- paste0(dirname(path), "/_viash.yaml")
  if (file.exists(check)) {
    dirname(check)
  } else if (check == "//_viash.yaml") {
    NULL
  } else {
    .ram_find_project(dirname(check))
  }
}

.ram_is_named_list <- function(obj) {
  is.null(obj) || (is.list(obj) && (length(obj) == 0 || !is.null(names(obj))))
}

.ram_process_merge <- function(data, root_data, path, project_path) {
  if (.ram_is_named_list(data)) {
    # check whether children have `__merge__` entries
    processed_data <- lapply(data, function(dat) {
      .ram_process_merge(dat, root_data, path, project_path)
    })
    processed_data <- lapply(names(data), function(nm) {
      dat <- data[[nm]]
      .ram_process_merge(dat, root_data, path, project_path)
    })
    names(processed_data) <- names(data)

    # if current element has __merge__, read list2 yaml and combine with data
    new_data <-
      if ("__merge__" %in% names(processed_data) && !.ram_is_named_list(processed_data$`__merge__`)) {
        new_data_path <- .ram_resolve_path(
          path = processed_data$`__merge__`,
          project_path = project_path,
          parent_path = dirname(path)
        )
        read_and_merge_yaml(new_data_path, project_path)
      } else if ("$ref" %in% names(processed_data) && !.ram_is_named_list(processed_data$`$ref`)) {
        ref_parts <- strsplit(processed_data$`$ref`, "#")[[1]]

        # resolve the path in $ref
        x <-
          if (ref_parts[[1]] == "") {
            root_data
          } else {
            new_data_path <- .ram_resolve_path(
              path = ref_parts[[1]],
              project_path = project_path,
              parent_path = dirname(path)
            )
            new_data_path <- normalizePath(new_data_path, mustWork = FALSE)

            # read in the new data
            tryCatch({
              suppressWarnings(yaml::read_yaml(new_data_path))
            }, error = function(e) {
              stop("Could not read ", new_data_path, ". Error: ", e)
            })
          }
        x_root <- x
        

        # Navigate the path and retrieve the referenced data
        ref_path_parts <- unlist(strsplit(ref_parts[[2]], "/"))
        for (part in ref_path_parts) {
          if (part == "") {
            next
          } else if (part %in% names(x)) {
            x <- x[[part]]
          } else {
            stop("Could not find ", processed_data$`$ref`, " in ", path)
          }
        }

        # postprocess the new data
        if (ref_parts[[1]] == "") {
          x
        } else {
          .ram_process_merge(x, x_root, new_data_path, project_path)
        }
      } else {
        list()
      }

    .ram_deep_merge(new_data, processed_data)
  } else if (is.list(data)) {
    lapply(data, function(dat) {
      .ram_process_merge(dat, root_data, path, project_path)
    })
  } else {
    data
  }
}

.ram_resolve_path <- function(path, project_path, parent_path) {
  ifelse(
    grepl("^/", path),
    paste0(project_path, "/", path),
    fs::path_abs(path, parent_path)
  )
}

.ram_deep_merge <- function(list1, list2) {
  if (.ram_is_named_list(list1) && .ram_is_named_list(list2)) {
    # if list1 and list2 are objects, recursively merge
    keys <- unique(c(names(list1), names(list2)))
    out <- lapply(keys, function(key) {
      if (key %in% names(list1)) {
        if (key %in% names(list2)) {
          .ram_deep_merge(list1[[key]], list2[[key]])
        } else {
          list1[[key]]
        }
      } else {
        list2[[key]]
      }
    })
    names(out) <- keys
    out
  } else if (is.list(list1) && is.list(list2)) {
    # if list1 and list2 are both lists, append
    c(list1, list2)
  } else {
    # else override list1 with list2
    list2
  }
}