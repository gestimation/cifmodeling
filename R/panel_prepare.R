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
    fonts <- .panel_extract_fonts(dots)
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

    args_est <- .panel_drop_nulls(list(
      formula        = formulas[[i]],
      data           = data,
      outcome.type   = if (!is.null(outcome.list)) outcome.list[[i]] else NULL,
      code.event1    = ce1,
      code.event2    = ce2,
      code.censoring = cc
    ))
    fit_i <- do.call(cifcurve, args_est)
    curves[[i]] <- fit_i

    args_plot <- .panel_drop_nulls(list(
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
