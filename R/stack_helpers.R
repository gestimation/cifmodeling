#' @keywords internal
#' @noRd
cifplot_build_stacked_survfit <- function(df_long, cause_order = NULL, add_t0 = TRUE) {
  if (!all(c("time", "cause", "cif") %in% names(df_long))) {
    stop("df_long must contain columns 'time', 'cause', and 'cif'.", call. = FALSE)
  }

  df <- df_long
  df$cause <- as.character(df$cause)
  df$time <- as.numeric(df$time)
  df$cif <- as.numeric(df$cif)

  if (!length(df$time)) {
    stop("df_long must contain at least one observation.", call. = FALSE)
  }

  if (isTRUE(add_t0)) {
    t0 <- data.frame(time = 0, cause = df$cause[1], cif = 0, stringsAsFactors = FALSE)
    df <- rbind(t0, df)
  }

  df <- df[order(df$cause, df$time), , drop = FALSE]
  df$cif[is.na(df$cif)] <- 0

  causes <- if (!is.null(cause_order) && length(cause_order)) {
    as.character(cause_order)
  } else {
    sort(unique(df$cause))
  }

  grid <- sort(unique(df$time))
  n_time <- length(grid)
  n_cause <- length(causes)

  if (!n_cause) stop("No causes detected for stacked CIF construction.", call. = FALSE)

  mat <- matrix(0, nrow = n_time, ncol = n_cause,
                dimnames = list(time = grid, cause = causes))

  for (j in seq_along(causes)) {
    cj <- causes[j]
    idx <- which(df$cause == cj)
    vals <- rep(NA_real_, n_time)
    if (length(idx)) {
      pos <- match(df$time[idx], grid)
      vals[pos] <- df$cif[idx]
    }
    last <- 0
    for (k in seq_along(vals)) {
      if (!is.na(vals[k])) last <- vals[k]
      vals[k] <- last
    }
    vals[!is.finite(vals)] <- 0
    vals <- pmin(pmax(vals, 0), 1)
    vals <- cummax(vals)
    mat[, j] <- vals
  }

  upper <- t(apply(mat, 1, cumsum))
  lower <- upper - mat
  upper <- pmin(upper, 1)
  lower <- pmax(lower, 0)

  time_v <- rep(grid, times = n_cause)
  lower_v <- as.numeric(lower)
  upper_v <- as.numeric(upper)
  counts <- rep.int(n_time, n_cause)
  names(counts) <- causes

  sf <- list(
    time      = time_v,
    n.risk    = rep(NA_integer_, length(time_v)),
    n.event   = rep(0L, length(time_v)),
    n.censor  = rep(0L, length(time_v)),
    surv      = upper_v,
    std.err   = rep(0, length(time_v)),
    lower     = lower_v,
    upper     = upper_v,
    conf.type = "plain",
    conf.int  = TRUE,
    type      = "right",
    strata    = counts,
    n         = counts
  )
  class(sf) <- "survfit"
  attr(sf, "stack_bounds") <- list(
    data = data.frame(
      time  = time_v,
      cause = factor(rep(causes, each = n_time), levels = causes),
      ymin  = lower_v,
      ymax  = upper_v
    )
  )
  sf
}

#' @keywords internal
#' @noRd
cifplot_draw_stacked_from_survfit <- function(sf, engine = c("ggsurvfit", "ggplot")) {
  engine <- match.arg(engine)
  if (engine == "ggsurvfit" && requireNamespace("ggsurvfit", quietly = TRUE)) {
    p <- ggsurvfit::ggsurvfit(sf) + ggsurvfit::add_confidence_interval()
  } else {
    bounds <- attr(sf, "stack_bounds")
    data_bounds <- bounds$data
    p <- ggplot2::ggplot(
      data_bounds,
      ggplot2::aes(x = time, ymin = ymin, ymax = ymax, fill = cause)
    ) + ggplot2::geom_ribbon()
  }
  p
}

#' @keywords internal
#' @noRd
cifplot_draw_stacked_from_survfit_list <- function(sflist, engine = c("ggsurvfit", "ggplot")) {
  engine <- match.arg(engine)
  nm <- names(sflist)
  if (is.null(nm) || !length(nm)) {
    nm <- paste0("Stratum ", seq_along(sflist))
  }
  plist <- lapply(seq_along(sflist), function(i) {
    cifplot_draw_stacked_from_survfit(sflist[[i]], engine = engine) +
      ggplot2::ggtitle(nm[[i]])
  })
  if (length(plist) == 1L) {
    return(plist[[1]])
  }
  if (requireNamespace("patchwork", quietly = TRUE)) {
    Reduce(`+`, plist) + patchwork::plot_layout(nrow = 1)
  } else {
    warning(
      "Multiple strata detected; showing first panel only. Install 'patchwork' for multipanel support.",
      call. = FALSE
    )
    plist[[1]]
  }
}

#' @keywords internal
#' @noRd
cifplot_survfit_to_stacked_survfits <- function(survfit_object, add_t0 = TRUE) {
  fit_type <- tolower(as.character(survfit_object$type %||% ""))
  if (!identical(fit_type, "aalen-johansen")) {
    stop("type.y = 'stacked' requires a competing-risk (Aalen-Johansen) survfit object.", call. = FALSE)
  }

  time_all <- as.numeric(survfit_object$time %||% numeric())
  surv_all <- as.numeric(survfit_object$surv %||% numeric())
  aj_all   <- as.numeric(survfit_object$aj %||% numeric())

  if (!length(time_all)) {
    stop("Survfit object has no time points for stacked plotting.", call. = FALSE)
  }

  if (!length(aj_all)) {
    aj_all <- pmax(0, 1 - surv_all)
  }

  aj_all <- aj_all[seq_along(time_all)]
  surv_all <- surv_all[seq_along(time_all)]

  if (!length(surv_all)) {
    stop("Survfit object lacks survival estimates required for stacked plotting.", call. = FALSE)
  }

  cause_names <- c("Event of interest", "Competing events")

  if (!is.null(survfit_object$strata) && length(survfit_object$strata)) {
    counts <- as.integer(survfit_object$strata)
    strata_names <- names(survfit_object$strata)
  } else {
    counts <- length(time_all)
    strata_names <- NULL
  }

  split_idx <- split(seq_along(time_all), rep.int(seq_along(counts), counts))
  if (!is.null(strata_names) && length(strata_names) == length(split_idx)) {
    names(split_idx) <- strata_names
  } else if (is.null(strata_names)) {
    if (length(split_idx) == 1L) {
      names(split_idx) <- ""
    } else {
      names(split_idx) <- paste0("Stratum ", seq_along(split_idx))
    }
  }

  sflist <- vector("list", length(split_idx))
  names(sflist) <- names(split_idx)

  for (i in seq_along(split_idx)) {
    idx <- split_idx[[i]]
    ti <- time_all[idx]
    si <- surv_all[idx]
    ci_main <- pmin(pmax(aj_all[idx], 0), 1)
    ci_other <- pmax(0, pmin(1, 1 - si - ci_main))

    df_i <- rbind(
      data.frame(time = ti, cause = cause_names[1], cif = ci_main, stringsAsFactors = FALSE),
      data.frame(time = ti, cause = cause_names[2], cif = ci_other, stringsAsFactors = FALSE)
    )

    sflist[[i]] <- cifplot_build_stacked_survfit(df_i, cause_order = cause_names, add_t0 = add_t0)
  }

  list(survfits = sflist, causes = cause_names)
}
