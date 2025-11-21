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

#' Methods for polyreg objects
#'
#' S3 methods to extract coefficients, variance-covariance matrix,
#' sample size, formatted summaries, and tidy/glance/augment
#' from objects returned by `polyreg()`.
#'
#' @name polyreg-methods
#'
#' @param object A `"polyreg"` object returned by `polyreg()`.
#' @param type Character string; one of `"default"`, `"sandwich"`,
#'   or `"bootstrap"`. When `"default"`, the function chooses between
#'   sandwich and bootstrap variance based on the original `polyreg()`
#'   settings, using `outcome.type`, `report.sandwich.conf`, and
#'   `report.boot.conf`. (Used only by `vcov.polyreg()`.)
#' @param x A `"summary.polyreg"`` object, as returned by
#'   `summary.polyreg()`.
#' @param digits Number of digits to print for parameter estimates.
#'   (Used only by `print.summary.polyreg()`.)
#' @param ... Further arguments passed to or from methods.
#'
#' @return
#' \itemize{
#'   \item `coef.polyreg()` returns a numeric vector of regression
#'     coefficients.
#'   \item `vcov.polyreg()` returns a variance-covariance matrix.
#'   \item `nobs.polyreg()` returns the number of observations.
#'   \item `summary.polyreg()` returns a list of tidy and glance
#'     summaries by event.
#'   \item `print.summary.polyreg()` is called for its side effect
#'     of printing a formatted, modelsummary-like table to the
#'     console and returns `x` invisibly.
#'   \item `tidy.polyreg()` returns a list of tidy by event.
#'   \item `glance.polyreg()` returns a list of glance by event.
#'   \item `augment.polyreg()` returns an augmented data frame.
#' }
#'
#' @seealso [polyreg()] for log-odds product modeling of CIFs
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
  report.sandwich.conf <- object$boot.method$report.sandwich.conf
  report.boot.conf     <- object$boot.method$report.boot.conf

  V_sand <- object$vcov
  V_boot <- object$vcov_bootstrap

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
    stop("No valid sandwich or bootstrap variance-covariance matrix.")
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
    stop("No valid bootstrap or sandwich variance-covariance matrix.")
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
print.summary.polyreg <- function(x,
                                  digits = 3,
                                  ...) {
  event_names <- names(x)

  cat("\n")
  cat(sprintf("%-20s", ""))
  for (ev in event_names) {
    cat(sprintf("  %-12s", ev))
  }
  cat("\n")

  cat(strrep("-", 20 + 2 + 12 * length(event_names)), "\n")

  terms_all <- unique(unlist(lapply(x, function(ev) ev$tidy$term)))

  for (term in terms_all) {

    cat(sprintf("%-20s", term), "\n")

    est_row  <- character(length(event_names))
    ci_row   <- character(length(event_names))
    p_row    <- character(length(event_names))

    for (j in seq_along(event_names)) {
      ev <- x[[event_names[j]]]
      td <- ev$tidy[ev$tidy$term == term, , drop = FALSE]
      if (nrow(td) == 0) {
        est_row[j] <- ""
        ci_row[j]  <- ""
        p_row[j]   <- ""
      } else {
        a <- exp(td$estimate[1])
        l <- exp(td$conf.low[1])
        h <- exp(td$conf.high[1])
        est_row[j] <- sprintf("%.*f", digits, a)
        ci_row[j]  <- sprintf("[%.*f, %.*f]",
                              digits, l,
                              digits, h)
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
      gl  <- x[[evnm]]$glance
      val <- if (rowname %in% names(gl)) gl[[rowname]][1] else ""
      cat(sprintf("  %-12s", as.character(val)))
    }
    cat("\n")
  }

  cat("\n")
  invisible(x)
}

#' @export
#' @rdname polyreg-methods
effect_label.polyreg <- function(x,
                                 event               = c("event1", "event2"),
                                 add.time.point      = TRUE,
                                 add.outcome         = TRUE,
                                 add.exposure.levels = TRUE,
                                 add.conf            = TRUE,
                                 add.p               = TRUE,
                                 value.time          = NULL,
                                 unit.time           = NULL,
                                 digits              = 2,
                                 p_digits            = 2,
                                 p_cut               = 0.05,
                                 ...) {

  event <- match.arg(event)

  if (event == "event1") {
    td <- generics::tidy(x, event = "event1")
    eff <- x$estimand$effect.measure1
    ev_id <- 1L
  } else {
    td <- generics::tidy(x, event = if (x$outcome.type %in% c("competing-risk",
                                                              "proportional-competing-risk")) {
      "event2"
    } else {
      "event1"
    })
    eff  <- x$estimand$effect.measure2
    ev_id <- 2L
  }

  exposure <- x$exposure
  is_exp <- grepl(paste0("^", exposure, ","), td$term)
  td_exp <- td[is_exp, , drop = FALSE]

  if (nrow(td_exp) == 0L) {
    stop("No exposure effect terms were found in tidy(x) for this polyreg object.",
         call. = FALSE)
  }

  if (is.null(value.time)) {
    tp <- x$estimand$time.point
    if (length(tp) == 1L && is.finite(tp)) {
      value.time <- tp
    }
  }

  if (is.null(unit.time)) {
    unit.time <- if (!is.null(value.time)) "years" else ""
  }

  labels <- character(nrow(td_exp))

  for (i in seq_len(nrow(td_exp))) {
    est <- exp(td_exp$estimate[i])
    lcl <- exp(td_exp$conf.low[i])
    ucl <- exp(td_exp$conf.high[i])
    p   <- td_exp$p.value[i]

    level_str <- ""
    if (add.exposure.levels) {
      s <- td_exp$term[i]
      s <- sub(paste0("^", exposure, ",\\s*"), "", s)
      level_str <- s
    }

    eff_part <- eff
    if (add.exposure.levels && nzchar(level_str)) {
      eff_part <- sprintf("%s (%s)", eff, level_str)
    }

    prefix <- eff_part
    if (add.outcome) {
      prefix <- paste(prefix, sprintf("of event %d", ev_id))
    }

    if (add.time.point && !is.null(value.time)) {
      if (nzchar(unit.time)) {
        tp_txt <- sprintf("at %g %s", value.time, unit.time)
      } else {
        tp_txt <- sprintf("at %g", value.time)
      }
      prefix <- paste(prefix, tp_txt)
    }

    est_txt <- sprintf("%.*f", digits, est)
    parts   <- character(0)

    if (add.conf) {
      ci_txt <- sprintf("%.*f to %.*f", digits, lcl, digits, ucl)
      parts  <- c(parts, ci_txt)
    }

    if (add.p && !is.na(p)) {
      if (p < p_cut) {
        p_txt <- sprintf("p < %.*f", p_digits, p_cut)
      } else {
        p_txt <- sprintf("p = %.*f", p_digits, p)
      }
      parts <- c(parts, p_txt)
    }

    if (length(parts) > 0L) {
      inner <- sprintf("%s (%s)", est_txt, paste(parts, collapse = ", "))
    } else {
      inner <- est_txt
    }
    labels[i] <- sprintf("%s = %s", prefix, inner)
  }
  paste(labels, collapse = "\n")
}

#' @export
#' @rdname polyreg-methods
tidy.polyreg <- function(x,
                         event = c("event1", "event2", "both"),
                         ...) {
  event <- match.arg(event)
  s     <- x$summary
  ot    <- x$outcome.type

  if (ot %in% c("competing-risk", "proportional-competing-risk")) {
    s1 <- s$event1
    s2 <- s$event2
  } else {
    s1 <- s[[1L]]
    s2 <- NULL
  }

  make_df <- function(part, label) {
    if (is.null(part)) return(NULL)
    df <- part$tidy
    df$event <- label
    df
  }

  if (event == "event1") {
    return(make_df(s1, "event1"))
  }

  if (event == "event2") {
    if (is.null(s2)) {
      stop("No competing event in this polyreg model.", call. = FALSE)
    }
    return(make_df(s2, "event2"))
  }

  out <- rbind(
    make_df(s1, "event1"),
    make_df(s2, "event2")
  )
  rownames(out) <- NULL
  return(out)
}

#' @export
#' @rdname polyreg-methods
glance.polyreg <- function(x,
                           event = c("event1", "event2"),
                           ...) {
  event <- match.arg(event)
  s     <- x$summary
  ot    <- x$outcome.type

  if (ot %in% c("competing-risk", "proportional-competing-risk")) {
    s1 <- s$event1
    s2 <- s$event2
  } else {
    s1 <- s[[1L]]
    s2 <- NULL
  }

  make_df <- function(part, label) {
    if (is.null(part)) return(NULL)
    df <- part$glance
    df$event <- label
    df
  }

  if (event == "event1") {
    return(make_df(s1, "event1"))
  }

  if (event == "event2") {
    if (is.null(s2)) {
      stop("No competing event in this polyreg model.", call. = FALSE)
    }
    return(make_df(s2, "event2"))
  }

  out <- rbind(
    make_df(s1, "event1"),
    make_df(s2, "event2")
  )
  rownames(out) <- NULL
  out
}

#' @export
#' @rdname polyreg-methods
augment.polyreg <- function(x, ...) {
  df <- x$diagnostics

  if ("ip.weight" %in% names(df)) {
    names(df)[names(df) == "ip.weight"] <- ".weights"
  }

  if ("potential.CIFs" %in% names(df)) {
    cif <- df$potential.CIFs
    df$potential.CIFs <- NULL

    if (is.matrix(cif) || is.data.frame(cif)) {
      cif_df <- as.data.frame(cif)
      if (is.null(colnames(cif_df))) {
        colnames(cif_df) <- paste0("event", seq_len(ncol(cif_df)))
      }
      colnames(cif_df) <- paste0(".fitted_", colnames(cif_df))
      df <- cbind(df, cif_df)
    } else {
      df$.fitted <- cif
    }
  }

  if ("influence.function" %in% names(df)) {
    names(df)[names(df) == "influence.function"] <- ".influence"
  }
  df
}
