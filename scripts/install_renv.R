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
  remote <- renv:::renv_remotes_resolve(remote)
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

install_with_retries <- function(remotes,
                                 attempts = 3,
                                 sleep = 3,
                                 backoff = 2,
                                 ...) {
  result <- NULL
  attempt <- 1
  while (is.null(result) && attempt <= attempts - 1) {
    attempt <- attempt + 1
    try(
      result <- renv::install(remotes, ...)
    )
    Sys.sleep(sleep)
    sleep <- sleep * backoff
  }
  if (is.null(result)) {
    # last attempt
    renv::install(remotes, ...)
  }
}

install_renv <- function(requirements_file, ...) {
  remotes <- scan(requirements_file, what = character(), sep = "\n")
  remotes <- sapply(remotes, strip_comments)
  remotes_installed <- sapply(remotes, check_available)
  remotes_to_install <- remotes[!remotes_installed]
  if (length(remotes_to_install) > 0) {
    install_with_retries(remotes_to_install, ...)
  }
}
