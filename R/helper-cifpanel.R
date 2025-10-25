#' Prepare curve fits and plotting arguments for cifpanel().
#'
#' This internal helper extracts per-panel arguments, evaluates
#' [cifcurve()] for each panel, and prepares the argument lists needed by
#' [cifplot()]. It keeps the logic centralized so that `cifpanel()` only
#' focuses on layout/patchwork composition.
#'
#' @param K Number of panels.
#' @param formulas List of formulas (one per panel).
#' @param data Data frame passed to [cifcurve()].
#' @param code.events List of numeric codes describing events per panel.
#' @param outcome.flags Character vector of normalized outcome flags
#'   (`"S"` or `"C"`).
#' @param outcome.list Optional list of outcome types (passed through to
#'   [cifcurve()]).
#' @param typey.list,labely.list,typex.list,labelx.list Panel-wise plot
#'   controls forwarded to [cifplot()].
#' @param limsx.list,limsy.list Panel-wise axis limits for [cifplot()].
#' @param breakx.list,breaky.list Panel-wise axis breaks for [cifplot()].
#' @param addCI.list,addCen.list,addCR.list,addIC.list,addQ.list Panel-wise
#'   logical toggles forwarded to [cifplot()].
#' @param strata.list Optional legend labels per panel.
#' @param legend.position Legend position forwarded to [cifplot()].
#' @param dots List of additional arguments (e.g., style, fonts).
#' @param fonts Optional list with `family` and `size`; if `NULL`, extracted
#'   from `dots` via [cifpanel()] helpers.
#'
#' @return A list with elements `curves`, `plot_args`, and `K`.
#' @keywords internal
#' @noRd
panel_prepare <- function(
    K,
    formulas,
    data,
    code.events,
    outcome.flags,
    outcome.list = NULL,
    typey.list = NULL,
    labely.list = NULL,
    typex.list = NULL,
    labelx.list = NULL,
    limsx.list = NULL,
    limsy.list = NULL,
    breakx.list = NULL,
    breaky.list = NULL,
    addCI.list = NULL,
    addCen.list = NULL,
    addCR.list = NULL,
    addIC.list = NULL,
    addQ.list = NULL,
    strata.list = NULL,
    legend.position = "top",
    dots = list(),
    fonts = NULL
) {
  if (is.null(fonts)) {
    fonts <- panel_extract_fonts(dots)
  }

  curves <- vector("list", length = K)
  plot_args <- vector("list", length = K)

  for (i in seq_len(K)) {
    pair <- code.events[[i]]
    if (outcome.flags[i] == "S") {
      ce1 <- pair[1]
      ce2 <- NULL
      cc  <- pair[2]
    } else {
      ce1 <- pair[1]
      ce2 <- pair[2]
      cc  <- pair[3]
    }

    norm_inputs <- cifplot_normalize_formula_data(formulas[[i]], data)
    data_i <- norm_inputs$data

    args_est <- panel_drop_nulls(list(
      formula        = formulas[[i]],
      data           = data_i,
      outcome.type   = if (!is.null(outcome.list)) outcome.list[[i]] else NULL,
      code.event1    = ce1,
      code.event2    = ce2,
      code.censoring = cc
    ))
    fit_i <- do.call(cifcurve, args_est)
    curves[[i]] <- fit_i

    args_plot <- panel_drop_nulls(list(
      x                       = fit_i,
      type.y                  = if (!is.null(typey.list))   typey.list[[i]]   else NULL,
      label.y                 = if (!is.null(labely.list))  labely.list[[i]]  else NULL,
      label.x                 = if (!is.null(labelx.list))  labelx.list[[i]]  else NULL,
      limits.y                = if (!is.null(limsy.list))   limsy.list[[i]]   else NULL,
      limits.x                = if (!is.null(limsx.list))   limsx.list[[i]]   else NULL,
      breaks.x                = if (!is.null(breakx.list))  breakx.list[[i]]  else NULL,
      breaks.y                = if (!is.null(breaky.list))  breaky.list[[i]]  else NULL,
      addConfidenceInterval    = if (!is.null(addCI.list))  addCI.list[[i]]   else TRUE,
      addCensorMark            = if (!is.null(addCen.list)) addCen.list[[i]]  else TRUE,
      addCompetingRiskMark     = if (!is.null(addCR.list))  addCR.list[[i]]   else FALSE,
      addIntercurrentEventMark = if (!is.null(addIC.list))  addIC.list[[i]]   else FALSE,
      addQuantileLine          = if (!is.null(addQ.list))   addQ.list[[i]]    else FALSE,
      label.strata             = if (!is.null(strata.list)) strata.list[[i]]  else NULL,
      style                    = dots$style %||% "CLASSIC",
      font.family              = fonts$family,
      font.size                = fonts$size,
      legend.position          = legend.position
    ))
    plot_args[[i]] <- args_plot
  }

  list(curves = curves, plot_args = plot_args, K = K)
}

panel_as_formula_global <- function(f) {
  if (is.character(f))   return(stats::as.formula(f, env = .GlobalEnv))
  if (inherits(f, "formula")) return(stats::as.formula(f, env = .GlobalEnv))
  .err("formula_must_be")
}

panel_recycle_to <- function(x, n) {
  if (length(x) == n) return(x)
  if (length(x) == 0L) .err("recycle_empty")
  rep_len(x, n)
}
panel_to_list <- function(x) if (is.null(x)) NULL else if (is.list(x)) x else as.list(x)
panel_drop_nulls <- function(lst) lst[!vapply(lst, is.null, logical(1))]
panel_strip_overrides_from_dots <- function(dots, override_names) {
  if (length(override_names) == 0L) return(dots)
  dots[setdiff(names(dots), override_names)]
}
panel_is_surv <- function(x) {
  x <- toupper(as.character(x %||% ""))
  x %in% c("S", "SURVIVAL")
}
panel_is_comp <- function(x) {
  x <- toupper(as.character(x %||% ""))
  x %in% c("C", "COMPETING-RISK", "COMPETING_RISK", "COMPETINGRISK")
}
panel_norm_outcome <- function(x) {
  if (panel_is_surv(x)) return("S")
  if (panel_is_comp(x)) return("C")
  stop("Unknown outcome.type: ", x, " (use 'S'/'SURVIVAL' or 'C'/'COMPETING-RISK').")
}
panel_validate_code_events <- function(code_events_list, outcome_flags) {
  stopifnot(length(code_events_list) == length(outcome_flags))
  for (i in seq_along(code_events_list)) {
    pair <- code_events_list[[i]]
    if (outcome_flags[i] == "S") {
      if (!(is.numeric(pair) && length(pair) == 2L))
        .err("code_events_len_surv", i = i)
    } else {
      if (!(is.numeric(pair) && length(pair) == 3L))
        .err("code_events_len_cr", i = i)
    }
  }
  invisible(TRUE)
}

panel_extract_fonts <- function(dots) {
  list(
    family = dots$font.family %||% "sans",
    size   = dots$font.size   %||% 12
  )
}
panel_build_theme <- function(font.family = "sans", font.size = 12) {
  base  <- font.size
  big   <- base * 1.20
  mid   <- base * 1.00
  small <- base * 0.85
  ggplot2::theme(
    text          = ggplot2::element_text(family = font.family, size = base),
    plot.title    = ggplot2::element_text(family = font.family, size = big, face = "bold"),
    plot.subtitle = ggplot2::element_text(family = font.family, size = mid),
    plot.caption  = ggplot2::element_text(family = font.family, size = small),
    axis.title.x  = ggplot2::element_text(family = font.family, size = mid),
    axis.title.y  = ggplot2::element_text(family = font.family, size = mid),
    axis.text.x   = ggplot2::element_text(family = font.family, size = small),
    axis.text.y   = ggplot2::element_text(family = font.family, size = small),
    legend.title  = ggplot2::element_text(family = font.family, size = mid),
    legend.text   = ggplot2::element_text(family = font.family, size = small),
    strip.text    = ggplot2::element_text(family = font.family, size = mid)
  )
}
