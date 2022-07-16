upgraded_remote_version <- function(remote) {
  if (remote$Source == "Repository") {
    return(paste0(remote$Package, "@", remote$Version))
  } else if (remote$Source == "GitHub") {
    out <- paste0(
      remote$RemoteUsername, "/", remote$RemoteRepo, "@", remote$RemoteSha[1:7]
    )
    if (remote$RemoteRef != remote$RemoteSha) {
      out <- paste0(out, " # ", remote$RemoteRef)
    }
    return(out)
  }
}

upgrade_first_available <- function(remotes) {
  for (remote in remotes) {
    cat(paste0(remote, "\n"))
    parsed_spec <- renv:::renv_remotes_parse(remote)
    result <- renv::update(parsed_spec$package, prompt = FALSE)
    if (class(result) != "logical") {
      upgraded_remotes <- sapply(result, upgraded_remote_version)
      return(paste(upgraded_remotes, collapse = "\n"))
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

upgrade_renv <- function(requirements_file) {
  remotes <- scan(requirements_file, what = character(), sep = "\n")
  if (length(remotes) > 0) {
    remotes <- drop_pinned(remotes)
    remotes <- sapply(remotes, strip_comments)
    capture.output(suppressWarnings(suppressMessages(
      output <- upgrade_first_available(remotes)
    )), file = nullfile())
    if (!is.null(output)) {
      cat("Upgrades are available:\n")
      cat(output)
      cat("\n")
    } else {
      cat("No upgrades available\n")
    }
  }
}

upgrade_renv(commandArgs(trailingOnly = TRUE)[1])
