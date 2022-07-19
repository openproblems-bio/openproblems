upgraded_remote_version <- function(remote) {
  if (remote$Source == "Repository") {
    out <- paste0(remote$Package, "@", remote$Version)
    if (remote$Repository == "BioCsoft") {
      out <- paste0("bioc::", out)
    }
  } else if (remote$Source == "GitHub") {
    out <- paste0(
      remote$RemoteUsername, "/", remote$RemoteRepo, "@", remote$RemoteSha[1:7]
    )
    if (remote$RemoteRef != remote$RemoteSha) {
      out <- paste0(out, " # ", remote$RemoteRef)
    }
  }
  return(out)
}

upgrade_first_available <- function(remotes) {
  for (remote in remotes) {
    cat(paste0(remote, "\n"))
    parsed_spec <- renv:::renv_remotes_parse(remote)
    result <- renv::update(parsed_spec$package, prompt = FALSE)
    if (class(result) != "logical") {
      upgraded_remotes <- sapply(result, upgraded_remote_version)
      return(upgraded_remotes)
    }
  }
}

check_pinned <- function(remote) {
  length(grep("#.*pinned", remote)) == 0
}

drop_pinned <- function(remotes) {
  remotes[sapply(remotes, check_pinned)]
}

strip_comments <- function(remote) {
  gsub("\\s*#.*", "", remote)
}

write_updates_to_file <- function(curr_remotes, new_remotes, filename) {
  cat("Running write_updates\n")
  current_names <- gsub("@.*", "", curr_remotes)
  new_names <- gsub("@.*", "", new_remotes)
  keep_current <- !(current_names %in% new_names)
  write_remotes <- sort(c(curr_remotes[keep_current], new_remotes))
  cat(format(write_remotes))
  writeLines(write_remotes, filename, sep = "\n")
  cat(format(scan(filename, what = character(), sep = "\n")))
}

upgrade_renv <- function(requirements_file) {
  remotes <- scan(requirements_file, what = character(), sep = "\n")
  if (length(remotes) > 0) {
    remotes_parsed <- drop_pinned(remotes)
    remotes_parsed <- sapply(remotes_parsed, strip_comments)
    capture.output(suppressWarnings(suppressMessages(
      upgraded_remotes <- upgrade_first_available(remotes_parsed)
    )), file = nullfile())
    if (!is.null(upgraded_remotes)) {
      cat("Upgrades are available:\n")
      cat(paste(upgraded_remotes, collapse = "\n"))
      cat("\n")
      write_updates_to_file(remotes, upgraded_remotes, requirements_file)
    } else {
      cat("No upgrades available\n")
    }
  }
}

upgrade_renv(commandArgs(trailingOnly = TRUE)[1])
