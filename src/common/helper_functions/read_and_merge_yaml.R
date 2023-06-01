#' Read a Viash YAML
#'
#' If the YAML contains a "__merge__" key anywhere in the yaml,
#' the path specified in that YAML will be read and the two
#' lists will be merged. This is a recursive procedure.
#'
#' @param path Path to Viash YAML
read_and_merge_yaml <- function(path) {
  data <- suppressWarnings(yaml::read_yaml(path))
  .ram_process_merge(data, path)
}

.ram_is_named_list <- function(obj) {
  is.null(obj) || (is.list(obj) && (length(obj) == 0 || !is.null(names(obj))))
}

.ram_process_merge <- function(data, path) {
  if (.ram_is_named_list(data)) {
    # check whether children have `__merge__` entries
    processed_data <- lapply(data, function(dat) {
      .ram_process_merge(dat, path)
    })
    names(processed_data) <- names(data)

    # if current element has __merge__, read list2 yaml and combine with data
    new_data <- 
      if ("__merge__" %in% names(processed_data)) {
        new_data_path <- paste0(dirname(path), "/", processed_data$`__merge__`)
        read_and_merge_yaml(new_data_path)
      } else {
        list()
      }

    .ram_deep_merge(new_data, processed_data)
  } else if (is.list(data)) {
    lapply(data, function(dat) {
      .ram_process_merge(dat, path)
    })
  } else {
    data
  }
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