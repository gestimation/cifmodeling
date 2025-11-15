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

#' Methods for \code{polyreg} objects
#'
#' S3 methods to extract coefficients, variance–covariance matrices,
#' sample size, and formatted summaries from objects returned by
#' \code{polyreg()}.
#'
#' @name polyreg-methods
#'
#' @param object A \code{"polyreg"} object returned by \code{polyreg()}.
#' @param type Character string; one of \code{"default"}, \code{"sandwich"},
#'   or \code{"bootstrap"}. When \code{"default"}, the function chooses between
#'   sandwich and bootstrap variance based on the original \code{polyreg()}
#'   settings, using \code{outcome.type}, \code{report.sandwich.conf}, and
#'   \code{report.boot.conf}. (Used only by \code{vcov.polyreg()}.)
#' @param x A \code{"summary.polyreg"} object, as returned by
#'   \code{summary.polyreg()}.
#' @param digits Number of digits to print for parameter estimates.
#'   (Used only by \code{print.summary.polyreg()}.)
#' @param ... Further arguments passed to or from methods.
#'
#' @return
#' \itemize{
#'   \item \code{coef.polyreg()} returns a numeric vector of regression
#'     coefficients.
#'   \item \code{vcov.polyreg()} returns a variance–covariance matrix.
#'   \item \code{nobs.polyreg()} returns the number of observations.
#'   \item \code{summary.polyreg()} returns a list of tidy and glance
#'     summaries by event.
#'   \item \code{print.summary.polyreg()} is called for its side effect
#'     of printing a formatted, \code{modelsummary}-like table to the
#'     console and returns \code{x} invisibly.
#' }
#'
#' @seealso \code{\link{polyreg}}
#'
#' @export
#' @rdname polyreg-methods
coef.polyreg <- function(object, ...) {
  object$coef
}

#' @export
#' @rdname polyreg-methods
vcov.polyreg <- function(object,
                         type = c("default", "sandwich", "bootstrap"),
                         ...) {

  type <- match.arg(type)

  outcome.type         <- object$outcome.type
  report.sandwich.conf <- object$optim.method$report.sandwich.conf
  report.boot.conf     <- object$optim.method$report.boot.conf

  V_sand <- object$vcov
  V_boot <- object$cov_bootstrap

  ok_mat <- function(V) {
    !is.null(V) && is.matrix(V) && length(V) > 0 && any(is.finite(V))
  }

  if (type == "default") {
    if (outcome.type %in% c("competing-risk", "survival", "binomial")) {
      type <- if (isFALSE(report.sandwich.conf)) "bootstrap" else "sandwich"
    } else if (outcome.type %in% c("proportional-survival", "proportional-competing-risk")) {
      type <- if (isFALSE(report.boot.conf)) "sandwich" else "bootstrap"
    } else {
      type <- "sandwich"
    }
  }

  if (type == "sandwich") {
    if (ok_mat(V_sand)) return(V_sand)
    if (ok_mat(V_boot)) {
      warning(
        "Sandwich variance is not available; ",
        "returning bootstrap variance instead.",
        call. = FALSE
      )
      return(V_boot)
    }
    stop("No valid sandwich or bootstrap variance–covariance matrix.")
  }

  if (type == "bootstrap") {
    if (ok_mat(V_boot)) return(V_boot)
    if (ok_mat(V_sand)) {
      warning(
        "Bootstrap variance is not available; ",
        "returning sandwich variance instead.",
        call. = FALSE
      )
      return(V_sand)
    }
    stop("No valid bootstrap or sandwich variance–covariance matrix.")
  }

  stop("Internal error in vcov.polyreg: unsupported 'type'.")
}

#' @export
#' @rdname polyreg-methods
nobs.polyreg <- function(object, ...) {
  nrow(object$diagnostics)
}

#' @export
#' @rdname polyreg-methods
summary.polyreg <- function(object, ...) {
  s <- object$summary
  class(s) <- c("summary.polyreg", class(s))
  s
}

#' @export
#' @rdname polyreg-methods
print.summary.polyreg <- function(summary,
                                  digits = 3,
                                  ...) {
  event_names <- names(summary)

  cat("\n")
  cat(sprintf("%-20s", ""))
  for (ev in event_names) {
    cat(sprintf("  %-12s", ev))
  }
  cat("\n")

  cat(strrep("-", 20 + 2 + 12 * length(event_names)), "\n")

  terms_all <- unique(unlist(lapply(summary, function(ev) ev$tidy$term)))

  for (term in terms_all) {

    cat(sprintf("%-20s", term), "\n")

    est_row  <- character(length(event_names))
    ci_row   <- character(length(event_names))
    p_row    <- character(length(event_names))

    for (j in seq_along(event_names)) {
      ev <- summary[[event_names[j]]]
      td <- ev$tidy[ev$tidy$term == term, , drop = FALSE]
      if (nrow(td) == 0) {
        est_row[j] <- ""
        ci_row[j]  <- ""
        p_row[j]   <- ""
      } else {
        est_row[j] <- sprintf("%.*f", digits, td$estimate[1])
        ci_row[j]  <- sprintf("[%.*f, %.*f]",
                              digits, td$conf.low[1],
                              digits, td$conf.high[1])
        p_row[j]   <- sprintf("(p=%0.3f)", td$p.value[1])
      }
    }

    cat(sprintf("%-20s", ""))
    for (val in est_row) cat(sprintf("  %-12s", val))
    cat("\n")

    cat(sprintf("%-20s", ""))
    for (val in ci_row) cat(sprintf("  %-12s", val))
    cat("\n")

    cat(sprintf("%-20s", ""))
    for (val in p_row) cat(sprintf("  %-12s", val))
    cat("\n\n")
  }

  cat(strrep("-", 20 + 2 + 12 * length(event_names)), "\n\n")

  glance_rows <- c(
    "effect.measure",
    "n.events",
    "median.follow.up",
    "range.follow.up",
    "n.parameters",
    "converged.by",
    "nleqslv.message"
  )

  for (rowname in glance_rows) {
    cat(sprintf("%-20s", rowname))

    for (evnm in event_names) {
      gl  <- summary[[evnm]]$glance
      val <- if (rowname %in% names(gl)) gl[[rowname]][1] else ""
      cat(sprintf("  %-12s", as.character(val)))
    }
    cat("\n")
  }

  cat("\n")
  invisible(summary)
}
