#' @export
print.cifplot <- function(x, ...) {
  p <- x$plot %||% x$patchwork
  if (!is.null(p)) print(p)
  invisible(x)
}

#' @export
print.cifpanel <- function(x, ...) {
  p <- x$patchwork
  if (!is.null(p)) print(p)
  invisible(x)
}
