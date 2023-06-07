strip_margin <- function(text, symbol = "\\|") {
  gsub(paste0("(^|\n)[ \t]*", symbol), "\\1", text)
}