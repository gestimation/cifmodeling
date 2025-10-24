.onLoad <- function(...) {
  op <- options()
  op.set <- list(
    cifplot.palette.default = NULL
  )
  toset <- !(names(op.set) %in% names(op))
  if (any(toset)) options(op.set[toset])
  invisible()
}

get_default_palette <- function(k) {
  pal <- getOption("cifplot.palette.default", NULL)
  if (is.null(pal)) {
    return(.default_fallback(k))
  }
  pal <- as.character(pal)
  if (!length(pal)) {
    return(.default_fallback(k))
  }
  .validate_color_names(pal)
  rep_len(unname(pal), k)
}
