if (nchar(Sys.getenv("BIOCVERSION")) > 0) {
  renv::settings$bioconductor.version(Sys.getenv("BIOCVERSION"))
}

as_integer_version <- function(v) {
  class(v) <- "list"
  v[[1]]
}

compare_version <- function(v1, v2) {
  v1 <- as_integer_version(v1)
  v2 <- as_integer_version(v2)
  for (i in seq_len(min(length(v1), length(v2)))) {
    if (v1[i] != v2[i]) {
      return(FALSE)
    }
  }
  return(TRUE)
}

check_available <- function(remote) {
  remote <- with_retries(renv:::renv_remotes_resolve, spec = remote)
  tryCatch(
    {
      version <- packageVersion(remote$Package)
      if (!is.null(remote$Version)) {
        compare_version(version, numeric_version(remote$Version))
      } else {
        TRUE
      }
    },
    error = function(...) {
      FALSE
    }
  )
}

strip_comments <- function(remote) {
  gsub("\\s*#.*", "", remote)
}

with_retries <- function(func,
                         attempts = 5,
                         sleep_once = 3,
                         sleep_multiple = 60,
                         backoff = 2,
                         ...) {
  result <- NULL
  attempt <- 1
  sleep <- sleep_once
  while (is.null(result) && attempt < attempts) {
    attempt <- attempt + 1
    try(
      result <- func(...)
    )
    closeAllConnections()
    Sys.sleep(sleep)
    if (sleep == sleep_once) {
      sleep <- sleep_multiple
    } else {
      sleep <- sleep * backoff
    }
  }
  if (is.null(result)) {
    # last attempt
    result <- func(...)
  }
  result
}

patch_renv <- function() {
  if (!requireNamespace("memoise", quietly = TRUE)) install.packages("memoise")
  renv_env <- parent.env(getNamespace("renv"))
  parent_env <- parent.env(renv_env)
  hack <- "__hack__"
  # set the new env between renv imports and base env, only if not already done
  if (!exists(hack, envir = parent.env(renv_env))) {
    # make a new env, with memoized renv_remotes_resolve
    renv_remotes_resolve_memoised <- memoise::memoise(
      renv:::renv_remotes_resolve
    )
    env <- new.env(parent = parent_env)
    assign(hack, TRUE, envir = env)
    assign("renv_remotes_resolve", renv_remotes_resolve_memoised, envir = env)
    parent.env(renv_env) <- env ### insert our custom env
  }
}

install_renv <- function(requirements_file, ...) {
  patch_renv()
  remotes <- scan(requirements_file, what = character(), sep = "\n")
  remotes <- sapply(remotes, strip_comments)
  remotes_installed <- sapply(remotes, check_available)
  remotes_to_install <- remotes[!remotes_installed]
  message(paste0("Installing ", length(remotes_to_install), " packages"))
  if (length(remotes_to_install) > 0) {
    with_retries(renv::install, packages = remotes_to_install, ...)
  }
}
