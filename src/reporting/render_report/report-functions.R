# Tables ----

#' Get references table
#'
#' @param references References list from results JSON
#'
#' @returns A `reactable` table containing the references
#'
#' @details
#' Information for DOI references is retrieved from CrossRef. BibTeX references
#' are formatted as code and ID references are shown as IDs.
get_references_table <- function(references) {
  if (all(c("doi", "bibtex") %in% names(references))) {
    references_df <- data.frame(
      reference_type = character(0),
      reference = character(0)
    )

    dois <- references$doi
    if (!(is.null(dois) || length(dois) == 0)) {
      doi_strs <- unlist(rcrossref::cr_cn(references$doi, format = "text"))
      references_df <- dplyr::bind_rows(
        references_df,
        data.frame(
          reference_type = "DOI",
          reference = doi_strs
        )
      )
    }

    bibtex <- references$bibtex
    if (!(is.null(bibtex) || length(bibtex) == 0)) {
      bibtex_strs <- purrr::map_chr(bibtex, function(.bibtex) {
        prettify_bibtex(.bibtex, output = "html")
      })
      references_df <- dplyr::bind_rows(
        references_df,
        data.frame(
          reference_type = "BibTeX",
          reference = bibtex_strs
        )
      )
    }
  } else {
    references_df <- data.frame(
      reference_type = "ID",
      reference = unlist(references)
    )
  }

  reactable::reactable(
    references_df,
    columns = list(
      reference = reactable::colDef(
        name = "References",
        cell = function(value, index, column) {
          reference_type <- references_df$reference_type[[index]]

          if (reference_type == "ID") {
            paste("ID:", value)
          } else {
            value
          }
        },
        style = function(value, row) {
          reference_type <- references_df$reference_type[[row]]

          if (reference_type == "BibTeX") {
            list("font-family" = "monospace")
          } else if (reference_type == "ID") {
            list("font-family" = "monospace")
          }
        },
        html = TRUE
      ),
      reference_type = reactable::colDef(show = FALSE)
    )
  )
}

#' Get source table
#'
#' @param details_df A data frame containing details
#' @param source_columns A character vector of column names to include in the
#'   source table
#'
#' @returns A `reactable` table containing the source information
#'
#' @details
#' The source columns are formatted as monospace text
get_source_table <- function(details_df, source_columns) {
  source_df <- details_df[, source_columns, drop = FALSE]

  reactable::reactable(
    source_df,
    columns = purrr::map(names(source_columns), function(.label) {
      reactable::colDef(
        name = .label,
        style = list("font-family" = "monospace")
      )
    }) |>
      purrr::set_names(source_columns),
    sortable = FALSE
  )
}


#' Get description table
#'
#' @param description_df A data frame containing the description information
#'
#' @returns A `reactable` table containing the description
#'
#' @details
#' The description Markdown is rendered as HTML
get_description_table <- function(description_df) {
  reactable::reactable(
    description_df,
    columns = list(
      description = reactable::colDef(
        name = "Description",
        cell = function(value) {
          commonmark::markdown_html(value)
        },
        html = TRUE
      )
    ),
    sortable = FALSE
  )
}

#' Get links table
#'
#' @param details_df A data frame containing details
#' @param link_columns A character vector of column names to include in the
#'   links table
#'
#' @returns A `reactable` table containing the links
#'
#' @details
#' The link columns are formatted as HTML links
get_links_table <- function(details_df, link_columns) {
  links_df <- details_df[, link_columns, drop = FALSE]

  reactable::reactable(
    links_df,
    columns = purrr::map(names(link_columns), function(.label) {
      reactable::colDef(
        name = .label,
        cell = format_html_link,
        html = TRUE
      )
    }) |>
      purrr::set_names(link_columns),
    sortable = FALSE
  )
}

#' Get additional information table
#'
#' @param additional_info A list containing additional information to display in
#'   a table
#'
#' @returns A `reactable` table containing the additional information or a HTML
#'   div with a message
#'
#' @details
#' Nicer heading are created from the column names, otherwise values are shown
#' as given. The additional information can contain any fields so we cannot
#' handle them specifically.
#'
#' If there are is no additional information, a div containing a messsage is
#' returned. A message is also returned if the additional information fails to
#' render.
get_additional_info_table <- function(additional_info) {
  if (is.null(additional_info) || length(additional_info) == 0) {
    return(htmltools::div(
      "No additional information found",
      style = "padding: 0.5rem"
    ))
  }

  tryCatch(
    {
      additional_data <- additional_info |>
        purrr::map(\(.x) {paste(.x, collapse = ", ")}) |>
        as.data.frame()

      colnames(additional_data) <- colnames(additional_data) |>
        stringr::str_replace_all("_", " ") |>
        stringr::str_to_sentence()

      reactable::reactable(
        additional_data
      )
    }, error = function(e) {
      htmltools::div(
        paste(
          "Additional information failed to render with error: ",
          e$message
        ),
        style = "padding: 0.5rem",
      )
    }
  )
}

#' Get quality control table
#'
#' @param qc_df A data frame containing quality control information
#'
#' @returns A `reactable` table containing the quality control checks
get_qc_table <- function(qc_df) {
  reactable::reactable(
    qc_df[, c("label", "severity")],

    columns = list(
      label = reactable::colDef(name = "Check"),
      severity = reactable::colDef(
        name = "Severity",
        cell = function(value) {
          switch(value,
                 "1" = "❌",
                 "2" = "❌❌",
                 "3" = "❌❌❌",
          )
        }
      )
    ),

    details = function(index, column) {
      details_df <- qc_df[index, , drop = FALSE]

      details_table <- reactable::reactable(
        details_df[, c("value", "condition", "severity_value")],
        columns = list(
          value = reactable::colDef(
            name = "Value",
            format = reactable::colFormat(digits = 2),
            width = 100
          ),
          condition = reactable::colDef(
            name = "Condition",
            style = list("font-family" = "monospace")
          ),
          severity_value = reactable::colDef(
            name = "Severity value",
            format = reactable::colFormat(digits = 2)
          )
        ),
        sortable = FALSE
      )

      message_table <- reactable::reactable(
        details_df[, "message", drop = FALSE],
        columns = list(
          message = reactable::colDef(
            name = "Message",
            cell = function(value) {
              stringr::str_replace_all(
                value,
                "\n",
                "<br>"
              )
            },
            html = TRUE
          )
        ),
        sortable = FALSE
      )

      htmltools::div(
        style = "padding: 1rem",
        details_table,
        message_table
      )
    },

    striped = TRUE,
    highlight = TRUE,
    defaultSorted = "severity",
    defaultSortOrder = "desc",
    defaultPageSize = 25,
    showPageSizeOptions = TRUE,

    rowStyle = reactable::JS(
      "function(rowInfo) {
        return {
          borderLeft: '2px solid #104E8B',
          fontWeight: 400
        }
      }"
    )
  )
}

# Plotting ----

#' Plot scaling
#'
#' @param complete_scores A long data frame containing all scaled metric scores
#' @param sel_metric The metric to plot
#' @param method_details A data frame containing method details
#' @param metric_details A data frame containing metric details
#'
#' @returns A `ggplot` object showing the scaling of the selected metric
#'
#' @details
#' Creates a normalization plot showing scaling of metric values, highlighting
#' values outside the [0, 1] range. A main panel shows all datasets and a
#' secondary panel is faceted by dataset.
plot_scaling <- function(
    complete_scores,
    sel_metric,
    method_details,
    metric_details
  ) {
  plot_data <- complete_scores |>
    dplyr::filter(metric == sel_metric) |>
    dplyr::mutate(
      method = factor(
        method,
        levels = method_details$method,
        labels = method_details$method_label
      ),
      method_type = factor(
        method_type,
        levels = sort(unique(method_type)),
        labels = sort(unique(method_type)) |>
          stringr::str_replace_all("_", " ") |>
          stringr::str_to_sentence()
      ),
    )

  norm_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = scaled_value, y = method)
  ) +
    ggplot2::annotate(
      geom = "rect",
      xmin = -Inf,
      xmax = 0,
      ymin = -Inf,
      ymax = Inf,
      fill = "red",
      alpha = 0.1
    ) +
    ggplot2::annotate(
      geom = "rect",
      xmin = 1,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      fill = "red",
      alpha = 0.1
    ) +
    ggplot2::geom_vline(
      xintercept = c(0, 1),
      linetype = "dashed",
      colour = "red"
    ) +
    ggplot2::geom_path(ggplot2::aes(group = dataset)) +
    ggplot2::geom_point(ggplot2::aes(colour = method_type)) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::labs(x = "Scaled value") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )

  norm_facets <- norm_plot +
    ggplot2::facet_wrap(
      ~ dataset,
      scales = "free_x",
      labeller = ggplot2::as_labeller(
        \(.x) {stringr::str_wrap(.x, width = 10, whitespace_only = FALSE)}
      )
    )

  norm_panel <- patchwork::wrap_plots(
    norm_plot + ggplot2::labs(title = "Overall"),
    norm_facets + ggplot2::labs(title = "By dataset"),
    ncol = 1,
    guides = "collect"
  ) &
    ggplot2::theme(legend.position = "bottom")

  norm_panel + patchwork::plot_annotation(
    title = metric_details$metric_label[metric_details$metric == sel_metric],
  )
}

# Formatting ----

#' Prettify BibTeX
#'
#' @param bibtex BibTeX string to prettify
#' @param output Output format, either "md" for Markdown or "html" for HTML
#'
#' @returns A prettified BibTeX string formatted for the specified output
prettify_bibtex <- function(bibtex, output = c("md", "html")) {
  output <- match.arg(output)

  newline <- switch(
    output,
    md = "\n",
    html = "<br>"
  )

  bibtex_str <- bibtex |>
    stringr::str_squish() |>
    stringr::str_replace(", ", paste0(",", newline, "    ")) |>
    stringr::str_replace_all("\\}, ", paste0("\\},", newline, "    ")) |>
    stringr::str_replace("\\s?\\}$", paste0(newline, "\\}"))

  if (output == "html") {
    bibtex_str <- paste0("<pre>", bibtex_str, "</pre>")
  }

  bibtex_str
}

#' Format HTML link
#'
#' @param value The URL to format as an HTML link
#'
#' @returns A string containing the HTML link
format_html_link <- function(value) {
  paste0("<a href='", value, "'>", value, "</a>")
}

#' Label memory
#'
#' @param x_mb A numeric vector of memory sizes in megabytes (MB)
#' @param include_mb A logical indicating whether to include label values less
#'   than 1 GB
#'
#' @returns A character vector with memory labels
label_memory <- function(x_mb, include_mb = TRUE) {
  dplyr::case_when(
    is.na(x_mb) | x_mb < 0 ~ "NA",
    x_mb < 1 ~ "<1M",
    x_mb < 1e3 & !include_mb ~ "<1G",
    x_mb < 1e3 ~ paste0(round(x_mb), "M"),
    x_mb < 1e6 ~ paste0(round(x_mb / 1e3), "G"),
    x_mb < 1e9 ~ paste0(round(x_mb / 1e6), "T"),
    TRUE ~ ">1P"
  )
}

#' Label time
#'
#' @param time A numeric vector of time values in seconds
#'
#' @returns A character vector with time labels
label_time <- function(time) {
  dplyr::case_when(
    is.na(time) | time < 0 ~ "NA",
    time < 1e-5 ~ "0s",
    time < 1 ~ "<1s",
    time < 60 ~ paste0(floor(time), "s"),
    time < 3600 ~ paste0(floor(time / 60), "m"),
    time < 3600 * 24 ~ paste0(floor(time / 3600), "h"),
    time < 3600 * 24 * 7 ~ paste0(floor(time / 3600 / 24), "d"),
    TRUE ~ ">7d"
  )
}

# Helpers ----

#' Aggregate scores
#'
#' @param scores A vector of scores to aggregate
#'
#' @returns An aggregated mean score
#'
#' @details
#' Values are restricted to between 0 and 1 and missing values are replaced by
#' 0. For use in creating the summary FunkyHeatmap
aggregate_scores <- function(scores) {
  scores[is.na(scores)] <- 0
  scores[scores < 0] <- 0
  scores[scores > 1] <- 1

  mean(scores, na.rm = TRUE)
}
