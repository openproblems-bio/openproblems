#' Read bibliography
#'
#' @param bib_file Path to a bibliography BibTex file
#'
#' @returns A list with two elements `doi` and `bibtex` where the names are
#'   reference keys and the values are the corresponding DOIs or BibTeX entries
read_bibliography <- function(bib_file) {
  bibentries <- bibtex::read.bib(bib_file)

  dois <- lapply(bibentries, function(.entry) {
    if (!is.null(.entry$doi)) {
      .entry$doi
    } else if (!is.null(.entry$DOI)) {
      .entry$DOI
    } else {
      NULL
    }
  }) |>
    purrr::compact()

  bibtex <- lapply(bibentries, function(.entry) {
    format(.entry, "bibtex")
  })

  list(
    doi = dois,
    bibtex = bibtex
  )
}

#' Get references list
#'
#' Convert a reference field from a config file into a references list
#' resolving any legacy reference keys
#'
#' @param reference The reference field value
#' @param bibliography A bibliography list as returned by `read_bibliography()`
#'
#' @returns A list with two elements `doi` and `bibtex`, where each is a character vector
#'   containing corresponding DOIs or BibTeX entries
get_references_list <- function(reference, bibliography) {
  # If null, return empty references
  if (
    is.null(reference) ||
      length(reference) == 0 ||
      (length(reference) == 1 && reference == "")
  ) {
    return(list(doi = character(0), bibtex = character(0)))
  }

  # If reference is a list, assume it is in the current format
  if (is.list(reference)) {
    return(
      list(
        doi = reference$doi %||% character(0),
        bibtex = reference$bibtex %||% character(0)
      )
    )
  }

  # If not a list, check if it is a DOI or BibTeX entry
  if (startsWith(reference, "@")) {
    return(list(doi = character(0), bibtex = reference))
  } else if (startsWith(reference, "1")) {
    return(list(doi = reference, bibtex = character(0)))
  }

  # Otherwise, assume it is a bibliography key
  if (reference %in% names(bibliography$doi)) {
    return(list(doi = bibliography$doi[[reference]], bibtex = character(0)))
  } else if (reference %in% names(bibliography$bibtex)) {
    return(list(doi = character(0), bibtex = bibliography$bibtex[[reference]]))
  } else {
    stop("Reference key '", reference, "' not found in bibliography")
  }
}

#' Get authors list
#'
#' Convert a list of authors from a config file into a structured list
#'
#' @param authors The authors field from a config file
#'
#' @returns An authors list in the expected format
get_authors_list <- function(authors) {
  `%||%` <- rlang::`%||%`

  purrr::map(authors, function(.author) {
    other_fields <- setdiff(names(.author$info), c("github", "orcid"))

    list(
      name = jsonlite::unbox(.author$name),
      roles = .author$roles %||% character(0),
      github = jsonlite::unbox(.author$info$github),
      orcid = jsonlite::unbox(.author$info$orcid),
      info = .author$info[other_fields] %||% setNames(list(), character(0))
    )
  })
}
