#' @title Visualize time-to-event outcomes and intercurrent events
#' @description
#' Estimate and plot survival curves using the Kaplan–Meier estimator or
#' cumulative incidence curves under competing risks using the Aalen–Johansen estimator.
#' Returns a \code{survfit}-compatible object and, by default, draws a publication-ready plot via \pkg{ggsurvfit}.
#'
#' @details
#' \strong{Estimation:}
#' \itemize{
#'   \item \code{outcome.type = "SURVIVAL"}: Kaplan–Meier estimator with Greenwood-type variance.
#'   \item \code{outcome.type = "COMPETING-RISK"}: Aalen–Johansen estimator for CIF of \code{code.event1}
#'         using IPCW for the censoring distribution. The returned \code{surv} corresponds to \code{1 - CIF}.
#' }
#' \strong{Confidence intervals:}
#' Constructed on the probability scale with the specified \code{conf.type}.
#' If \code{conf.type \%in\% c("none","n")}, the plot suppresses CI bands.
#'
#' \strong{Plotting:}
#' By default, the function calls an internal helper \code{call_ggsurvfit()} which adds
#' confidence bands, risk table, censoring marks, and optional intercurrent-event marks.
#' For CIF display, set \code{ggsurvfit.type = "risk"}.
#'
#' @param formula A model formula specifying the outcome and (optionally) \code{strata()}.
#' @param data A data frame containing variables in \code{formula}.
#' @param weights Optional name of the weight variable in \code{data}. Weights must be nonnegative; strictly positive is recommended.
#' @param subset.condition Optional character expression to subset \code{data} before analysis.
#' @param na.action Function to handle missing values (default: \code{\link[stats]{na.omit}}).
#' @param code.event1 Integer code of the event of interest (default \code{1}).
#' @param code.event2 Integer code of the competing event (default \code{2}).
#' @param code.censoring Integer code of censoring (default \code{0}).
#' @param outcome.type \code{"SURVIVAL"} (KM) or \code{"COMPETING-RISK"} (AJ).
#' @param conf.int Two-sided confidence level (default \code{0.95}).
#' @param error Character specifying variance type used internally.
#'   For \code{"SURVIVAL"} typically \code{"greenwood"}; for \code{"COMPETING-RISK"} pass options supported by \code{calculateAalenDeltaSE()} (e.g., \code{"aalen"}, \code{"delta"}, \code{"none"}).
#' @param conf.type Transformation for CI on the probability scale (default \code{"arcsine-square root"}).
#' @param report.survfit.std.err If \code{TRUE}, report SE on the log-survival scale (survfit's convention). Otherwise SE is on the probability scale.
#' @param report.ggsurvfit If \code{TRUE} (default), draw a \pkg{ggsurvfit} plot.
#' @param ggsurvfit.type \code{NULL} (survival) or \code{"risk"} (display \code{1 - survival} i.e. CIF).
#' @param addConfidenceInterval Logical add \code{add_confidence_interval()} to plot (default \code{TRUE}).
#' @param addRiskTable Logical add \code{add_risktable(risktable_stats="n.risk")} to plot (default \code{TRUE}).
#' @param addCensorMark Logical add \code{add_censor_mark()} to plot (default \code{TRUE}).
#' @param shape.censor.mark Integer point shape for censor marks (default \code{3}).
#' @param size.censor.mark Numeric point size for censor marks (default \code{2}).
#' @param addCompetingRiskMark Logical add time marks to describe event2 specified by Event(), usually the competing events (default \code{TRUE}).
#' @param shape.competing.risk.mark Integer point shape for competing-risk marks (default \code{16}).
#' @param size.competing.risk.mark Numeric point size for competing-risk marks (default \code{2}).
#' @param addIntercurrentEventMark Logical overlay user-specified time marks per strata (default \code{TRUE}).
#' @param intercurrent.event.time Named list of numeric vectors (names must match or be mapped to strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks (default \code{1}).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks (default \code{2}).
#' @param label.strata Character vector of labels for strata.
#' @param label.x Axis labels (defaults: \code{"Time"}).
#' @param label.y Axis labels (defaults: \code{"Survival probability"}).
#'   If \code{ggsurvfit.type="risk"} and \code{label.y} is unchanged, it is internally set to \code{"1 - survival probability"}.
#' @param lims.x Numeric length-2 vectors for axis limits (defaults: \code{NULL}, \code{c(0,1)}).
#' @param lims.y Numeric length-2 vectors for axis limits (defaults: \code{NULL}, \code{c(0,1)}).
#' @param font.family Character plot theme controls (defaults: \code{"sans"}).
#' @param font.size Integer plot theme controls (defaults: \code{14}).
#' @param legend.position Legend position: \code{"top"}, \code{"right"}, \code{"bottom"}, \code{"left"}, or \code{"none"} (default \code{"top"}).
#' @param filename Character specify the name of the graph output.
#'
#' @returns A \code{survfit} object. For \code{outcome.type="SURVIVAL"}, \code{$surv} is the survival function.
#' For \code{outcome.type="COMPETING-RISK"}, \code{$surv} equals \code{1 - CIF} for \code{code.event1}.
#' Standard error and CIs are provided per \code{conf.type}. Note that some methods for \code{survfit} (e.g., \code{residuals.survfit}) may not be supported.
#'
#' @seealso \code{\link{polyreg}} for log-odds product models of CIFs; \pkg{ggsurvfit} for plotting helpers.
#'
#' @examples
#' data(diabetes.complications)
#' survfit_by_group <- cifcurve(Event(t,epsilon) ~ fruitq, data = diabetes.complications,
#'                     outcome.type='COMPETING-RISK', error='delta', ggsurvfit.type = 'risk',
#'                     label.y = 'CIF of diabetic retinopathy', label.x = 'Years from registration')
#'
#' @importFrom ggsurvfit ggsurvfit add_confidence_interval add_risktable add_censor_mark
#' @importFrom ggplot2 theme_classic theme element_text labs lims geom_point aes ggsave
#' @importFrom Rcpp sourceCpp
#' @useDynLib cifmodeling, .registration = TRUE
#' @export
cifcurve <- function(formula,
                     data,
                     weights = NULL,
                     subset.condition = NULL,
                     na.action = na.omit,
                     code.event1 = 1,
                     code.event2 = 2,
                     code.censoring = 0,
                     outcome.type = "SURVIVAL",
                     conf.int = 0.95,
                     error = "greenwood",
                     conf.type = "arcsine-square root",
                     report.survfit.std.err = FALSE,
                     report.ggsurvfit = TRUE,
                     ggsurvfit.type = NULL,
                     addConfidenceInterval = FALSE,
                     addRiskTable = TRUE,
                     addCensorMark = TRUE,
                     addCompetingRiskMark = FALSE,
                     addIntercurrentEventMark = FALSE,
                     label.x = "Time",
                     label.y = "Survival probability",
                     label.strata = NULL,
                     lims.x = NULL,
                     lims.y = c(0, 1),
                     shape.censor.mark = 3,
                     size.censor.mark = 2,
                     shape.competing.risk.mark = 16,
                     size.competing.risk.mark = 2,
                     intercurrent.event.time = NULL,
                     shape.intercurrent.event.mark = 1,
                     size.intercurrent.event.mark = 2,
                     font.family = "sans",
                     font.size = 14,
                     legend.position = "top",
                     filename = NULL) {

  outcome.type  <- check_outcome.type(outcome.type)
  out_readSurv  <- readSurv(formula, data, weights, code.event1, code.event2, code.censoring, subset.condition, na.action)
  out_mcm       <- make_competingrisk_marks(out_readSurv, event_code = code.event2)

  error <- check_error(error, outcome.type)
  check_label.strata(out_readSurv, label.strata)
  lev_old <- levels(as.factor(out_readSurv$strata))
  lab_map <- setNames(if (is.null(label.strata)) lev_old else label.strata, lev_old)

  if (identical(outcome.type, "SURVIVAL")) {
    out_km <- calculateKM(out_readSurv$t, out_readSurv$d, out_readSurv$w, as.integer(out_readSurv$strata), error)
    out_km$std.err <- out_km$surv * out_km$std.err
    out_ci <- calculateCI(out_km, conf.int, conf.type, conf.lower = NULL)
    if (isTRUE(report.survfit.std.err)) out_km$std.err <- out_km$std.err / out_km$surv

    survfit_object <- list(
      time = out_km$time,
      surv = out_km$surv,
      n = out_km$n,
      n.risk = out_km$n.risk,
      n.event = out_km$n.event,
      n.censor = out_km$n.censor,
      std.err = out_km$std.err,
      upper = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$upper,
      lower = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$lower,
      conf.type = conf.type,
      call = match.call(),
      type = "kaplan-meier",
      method = "Kaplan-Meier"
    )
    names(out_km$strata) <- levels(as.factor(out_readSurv$strata))
    if (any(as.integer(out_readSurv$strata) != 1)) survfit_object$strata <- out_km$strata
    class(survfit_object) <- "survfit"

  } else {
    out_aj <- calculateAJ(out_readSurv)
    names(out_aj$strata1) <- levels(as.factor(out_readSurv$strata))

    if (any(as.integer(out_readSurv$strata) != 1)) {
      n <- table(as.integer(out_readSurv$strata))
      rep_list <- mapply(rep, n, out_aj$strata1, SIMPLIFY = FALSE)
      n.risk <- do.call(c, rep_list) - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
    } else {
      n <- length(out_readSurv$strata)
      n.risk <- n - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
    }
    out_aj$std.err <- calculateAalenDeltaSE(out_aj$time1, out_aj$aj1, out_aj$n.event1, out_aj$n.event2,
                                            n.risk, out_aj$time0, out_aj$km0, out_aj$strata1, error)
    out_aj$surv <- 1 - out_aj$aj1
    out_ci <- calculateCI(out_aj, conf.int, conf.type, conf.lower = NULL)
    if (isTRUE(report.survfit.std.err)) out_aj$std.err <- out_aj$std.err / out_aj$surv

    survfit_object <- list(
      time = out_aj$time1,
      surv = out_aj$surv,
      n = n,
      n.risk = n.risk,
      n.event = out_aj$n.event1,
      n.censor = out_aj$n.censor,
      std.err = out_aj$std.err,
      upper = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$upper,
      lower = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$lower,
      conf.type = conf.type,
      call = match.call(),
      type = "Aalen-Johansen",
      method = "Aalen-Johansen"
    )
    if (any(as.integer(out_readSurv$strata) != 1)) survfit_object$strata <- out_aj$strata1
    class(survfit_object) <- "survfit"
  }

  if (isTRUE(report.ggsurvfit) || !is.null(filename)) {
    out_ggsurvfit <- call_ggsurvfit(
      survfit_object = survfit_object,
      out_readSurv   = out_readSurv,
      ggsurvfit.type = ggsurvfit.type,
      conf.type      = conf.type,
      addConfidenceInterval = addConfidenceInterval,
      addRiskTable   = addRiskTable,
      addCensorMark  = addCensorMark,
      shape.censor.mark = shape.censor.mark,
      size.censor.mark  = size.censor.mark,
      addCompetingRiskMark = addCompetingRiskMark,
      competing.risk.time  = out_mcm,
      shape.competing.risk.mark = shape.competing.risk.mark,
      size.competing.risk.mark  = size.competing.risk.mark,
      addIntercurrentEventMark = addIntercurrentEventMark,
      intercurrent.event.time  = intercurrent.event.time,
      shape.intercurrent.event.mark = shape.intercurrent.event.mark,
      size.intercurrent.event.mark  = size.intercurrent.event.mark,
      label.x = label.x, label.y = label.y, label.strata = lab_map,
      lims.x = lims.x, lims.y = lims.y,
      font.family = font.family, font.size = font.size,
      legend.position = legend.position
    )
    if (isTRUE(report.ggsurvfit)) print(out_ggsurvfit)
    if (!is.null(filename)) ggsave(filename, plot = out_ggsurvfit, width = 6, height = 6, dpi = 300)
  }
  survfit_object
}

#' Plot survival or cumulative incidence curves with ggsurvfit
#'
#' @param survfit_object A \code{survfit} object.
#' @param out_readSurv (optional) List returned by your \code{readSurv()} to auto-set x limits.
#' @param ggsurvfit.type Character; NULL (survival) or "risk" for cumulative incidence display.
#' @param conf.type Character; same as used when constructing CI (e.g., "none", "n", "arcsine-square root").
#' @param addConfidenceInterval Logical; add CI via \code{add_confidence_interval()}.
#' @param addRiskTable Logical; add risk table via \code{add_risktable(risktable_stats="n.risk")}.
#' @param addCensorMark Logical add \code{add_censor_mark()} to plot (default \code{TRUE}).
#' @param shape.censor.mark Integer point shape for censor marks (default \code{3}).
#' @param size.censor.mark Numeric point size for censor marks (default \code{2}).
#' @param addCompetingRiskMark Logical add time marks to descrive incidents of competing risks (default \code{TRUE}).
#' @param competing.risk.time Named list of numeric vectors (names must match or be mapped to strata labels).
#' @param shape.competing.risk.mark Integer point shape for competing-risk marks (default \code{16}).
#' @param size.competing.risk.mark Numeric point size for competing-risk marks (default \code{2}).
#' @param addIntercurrentEventMark Logical overlay user-specified time marks per strata (default \code{TRUE}).
#' @param intercurrent.event.time Named list of numeric vectors (names must match or be mapped to strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks (default \code{1}).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks (default \code{2}).
#' @param label.x,label.y Axis labels. If \code{label.y = "Survival probability"} (default) and
#'   \code{ggsurvfit.type == "risk"}, it is automatically replaced with \code{"1 - survival probability"}.
#' @param label.strata Character vector of labels for strata.
#' @param lims.x Numeric length-2; x limits. If NULL and \code{out_readSurv} given, uses \code{c(0,max(out_readSurv$t))}.
#' @param lims.y Numeric length-2; y limits.
#' @param font.family,font.size Theme controls.
#' @param legend.position "top","right","bottom","left" or "none".
#' @return A \code{ggplot} object.
#' Plot survival or cumulative incidence curves with ggsurvfit
#' ...
#' ...
call_ggsurvfit <- function(
    survfit_object,
    out_readSurv = NULL,
    ggsurvfit.type = NULL,
    conf.type = "arcsine-square root",
    addConfidenceInterval = TRUE,
    addRiskTable = TRUE,
    addCensorMark = TRUE,
    shape.censor.mark = 3,
    size.censor.mark = 2,
    addCompetingRiskMark = TRUE,
    competing.risk.time = list(),
    shape.competing.risk.mark = 16,
    size.competing.risk.mark = 2,
    addIntercurrentEventMark = TRUE,
    intercurrent.event.time = list(),
    shape.intercurrent.event.mark = 1,
    size.intercurrent.event.mark = 2,
    label.x = "Time",
    label.y = "Survival probability",
    label.strata = NULL,
    lims.x = NULL,
    lims.y = c(0, 1),
    font.family = "sans",
    font.size = 14,
    legend.position = "top"
){
  if (is.null(lims.x) && !is.null(out_readSurv) && !is.null(out_readSurv$t)) {
    lims.x <- c(0, max(out_readSurv$t))
  }

  label.y.corrected <- if (identical(ggsurvfit.type, "risk") && identical(label.y, "Survival probability")) {
    "1 - survival probability"
  } else label.y

  check_ggsurvfit(
    survfit_object = survfit_object,
    lims.x = lims.x, lims.y = lims.y,
    ggsurvfit.type = ggsurvfit.type,
    addConfidenceInterval = addConfidenceInterval,
    addCensorMark = addCensorMark,
    addCompetingRiskMark = addCompetingRiskMark,
    addIntercurrentEventMark = addIntercurrentEventMark,
    shape.censor.mark = shape.censor.mark,
    shape.competing.risk.mark = shape.competing.risk.mark,
    shape.intercurrent.event.mark = shape.intercurrent.event.mark
  )

  select_y_axis <- function(sfobj) {
    if (identical(ggsurvfit.type, "risk")) ggsurvfit(sfobj, type = "risk") else ggsurvfit(sfobj)
  }

  sf_for_plot <- coerce_ci_if_needed(survfit_object, conf.type)

  p <- select_y_axis(sf_for_plot) +
    base_surv_theme(font.family, font.size, legend.position) +
    labs(x = label.x, y = label.y.corrected) +
    lims(x = lims.x, y = lims.y)
  p <- p +
    ggplot2::scale_color_discrete(
      breaks = names(label.strata),
      labels = unname(label.strata)
    ) +
    ggplot2::scale_fill_discrete(
      breaks = names(label.strata),
      labels = unname(label.strata)
    )

  if (isTRUE(addConfidenceInterval)) p <- p + add_confidence_interval()
  if (isTRUE(addRiskTable))          p <- p + add_risktable(risktable_stats = c("n.risk"))
  if (isTRUE(addCensorMark))         p <- p + add_censor_mark(shape = shape.censor.mark, size = size.censor.mark)
  if (isTRUE(addCompetingRiskMark) && length(competing.risk.time)) {
    p <- draw_marks_if_any(p, survfit_object, competing.risk.time, ggsurvfit.type,
                           shape = shape.competing.risk.mark, size = size.competing.risk.mark)
  }
  if (isTRUE(addIntercurrentEventMark) && length(intercurrent.event.time)) {
    p <- draw_marks_if_any(p, survfit_object, intercurrent.event.time, ggsurvfit.type,
                           shape = shape.intercurrent.event.mark, size = size.intercurrent.event.mark)
  }
  return(p)
}

base_surv_theme <- function(font.family, font.size, legend.position) {
  theme_classic() +
    theme(
      legend.position = legend.position,
      axis.title = element_text(size = (font.size + 2), family = font.family),
      axis.text  = element_text(size = font.size, family = font.family),
      legend.text = element_text(size = font.size, family = font.family)
    )
}

coerce_ci_if_needed <- function(survfit_object, conf.type) {
  if (!is.null(survfit_object$lower) && !is.null(survfit_object$upper)) return(survfit_object)
  if (conf.type %in% c("none","n") || length(survfit_object$strata) > 2) {
    x <- survfit_object
    x$lower <- x$surv
    x$upper <- x$surv
    return(x)
  }
  survfit_object
}

draw_marks_if_any <- function(p, survfit_object, marks, ggsurvfit.type, shape, size) {
  if (is.null(marks) || !length(marks)) return(p)
  mark_df <- make_mark_df(survfit_object, marks, ggsurvfit.type, extend = TRUE)
  if (is.null(mark_df) || !nrow(mark_df)) return(p)

  if (is.null(survfit_object$strata)) {
    p + geom_point(
      data = mark_df,
      aes(x = time, y = y),
      inherit.aes = FALSE,
      shape = shape,
      size  = size,
      show.legend = FALSE
    )
  } else {
    p + geom_point(
      data = mark_df,
      aes(x = time, y = y, group = strata, colour = strata),
      inherit.aes = FALSE,
      shape = shape,
      size  = size,
      show.legend = FALSE
    )
  }
}

align_marks_keys <- function(fit, marks) {
  canon_str <- function(x) sub("^.*=", "", as.character(x))
  if (is.null(marks) || length(marks) == 0) return(marks)
  fit_names <- if (is.null(fit$strata)) "(all)" else names(fit$strata)
  fit_key   <- canon_str(fit_names)
  out <- list()
  for (k in names(marks)) {
    k_can <- canon_str(k)
    j <- match(k_can, fit_key)
    if (is.na(j)) next
    out[[ fit_names[j] ]] <- marks[[k]]
  }
  out
}

make_mark_df <- function(fit, marks, ggsurvfit.type, extend = TRUE) {
  if (is.null(marks) || length(marks) == 0) return(NULL)
  marks2 <- align_marks_keys(fit, marks)
  strata_names <- if (is.null(fit$strata)) "(all)" else names(fit$strata)

  out <- lapply(names(marks2), function(st_lab) {
    tt <- marks2[[st_lab]]
    if (length(tt) == 0) return(NULL)
    fit_st <- if (is.null(fit$strata)) {
      if (!identical(st_lab, "(all)")) return(NULL)
      fit
    } else {
      i <- match(st_lab, strata_names); if (is.na(i)) return(NULL)
      fit[i]
    }
    sm <- summary(fit_st, times = tt, extend = extend)
    y  <- if (identical(ggsurvfit.type, "risk")) 1 - sm$surv else sm$surv
    data.frame(strata = st_lab, time = sm$time, y = y, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

calculateAJ <- function(data) {
  out_km0 <- calculateKM(data$t, data$d0, data$w, as.integer(data$strata), "none")
  km0 <- get_surv(data$t, out_km0$surv, out_km0$time, as.integer(data$strata), out_km0$strata, out_km0$strata.levels)
  ip.weight <- (data$d0 == 0) * ifelse(km0 > 0, 1 / km0, 0)
  d1_ipw <- as.matrix(data$w * data$d1 * ip.weight)

  aj1 <- time1 <- integer(0)
  n.cum.event1 <- n.cum.event2 <- n.cum.censor <- numeric(0)
  n.event1 <- n.event2 <- n.censor <- numeric(0)
  strata1 <- integer(0)
  strata_vec <- as.integer(data$strata)

  for (level in sort(unique(strata_vec))) {
    idx <- which(strata_vec == level)

    sub_t  <- data$t[idx]
    sub_d0 <- data$d0[idx]
    sub_d1 <- data$d1[idx]
    sub_d2 <- data$d2[idx]
    sub_d1_ipw <- d1_ipw[idx, , drop = FALSE]

    o <- order(sub_t)
    sub_t  <- sub_t[o]
    sub_d0 <- sub_d0[o]
    sub_d1 <- sub_d1[o]
    sub_d2 <- sub_d2[o]
    sub_d1_ipw <- sub_d1_ipw[o, , drop = FALSE]

    not_atrisk <- outer(sub_t, sub_t, FUN = ">=")
    sub_aj1 <- as.vector(not_atrisk %*% sub_d1_ipw) / length(sub_t)
    sub_n.censor <- as.vector(not_atrisk %*% as.matrix(sub_d0))
    sub_n.event1 <- as.vector(not_atrisk %*% as.matrix(sub_d1))
    sub_n.event2 <- as.vector(not_atrisk %*% as.matrix(sub_d2))

    keep <- !duplicated(rev(sub_t))
    keep <- rev(keep)

    u_t  <- sub_t[keep]
    u_aj1 <- sub_aj1[keep]
    u_nc  <- sub_n.censor[keep]
    u_ne1 <- sub_n.event1[keep]
    u_ne2 <- sub_n.event2[keep]

    oo <- order(u_t)
    u_t  <- u_t[oo]
    u_aj1 <- u_aj1[oo]
    u_nc  <- u_nc[oo]
    u_ne1 <- u_ne1[oo]
    u_ne2 <- u_ne2[oo]

    inc_nc  <- c(u_nc[1],  diff(u_nc))
    inc_ne1 <- c(u_ne1[1], diff(u_ne1))
    inc_ne2 <- c(u_ne2[1], diff(u_ne2))

    time1 <- c(time1, u_t)
    aj1   <- c(aj1, u_aj1)

    n.cum.censor <- c(n.cum.censor, u_nc)
    n.cum.event1 <- c(n.cum.event1, u_ne1)
    n.cum.event2 <- c(n.cum.event2, u_ne2)

    n.censor <- c(n.censor, inc_nc)
    n.event1 <- c(n.event1, inc_ne1)
    n.event2 <- c(n.event2, inc_ne2)

    strata1 <- c(strata1, length(u_t))
  }

  list(
    time1 = time1,
    aj1 = aj1,
    n.event1 = n.event1,
    n.event2 = n.event2,
    n.censor = n.censor,
    n.cum.event1 = n.cum.event1,
    n.cum.event2 = n.cum.event2,
    n.cum.censor = n.cum.censor,
    strata1 = strata1,
    time0 = out_km0$time,
    km0 = out_km0$surv
  )
}


get_surv <- function(
    predicted.time,
    estimated.surv,
    estimated.time,
    predicted.strata = NULL,
    estimated.strata = NULL,
    strata.levels = NULL
){
  if (anyNA(predicted.time)) stop("Invalid predicted.time: contains NA.")
  if (length(estimated.surv) != length(estimated.time))
    stop("estimated.surv and estimated.time must have the same length.")

  prepareSeries <- function(time_vec, surv_vec) {
    ok <- !(is.na(time_vec) | is.na(surv_vec))
    time_vec <- time_vec[ok]; surv_vec <- surv_vec[ok]
    if (!length(time_vec)) return(list(t = numeric(0), s = numeric(0)))
    o <- order(time_vec)
    t2 <- time_vec[o]; s2 <- surv_vec[o]
    keep <- !duplicated(t2, fromLast = TRUE)
    list(t = t2[keep], s = s2[keep])
  }

  n_pred <- length(predicted.time)
  predicted.surv <- numeric(n_pred)

  strata_mode <- !(
    is.null(predicted.strata) || is.null(estimated.strata) || is.null(strata.levels) ||
      length(estimated.strata) == 0L || length(strata.levels) == 0L
  )

  if (!strata_mode) {
    ser <- prepareSeries(estimated.time, estimated.surv)
    if (!length(ser$t)) return(rep(1.0, n_pred))
    for (i in seq_len(n_pred)) {
      idx <- findInterval(predicted.time[i], ser$t, left.open = TRUE)
      predicted.surv[i] <- if (idx > 0L) ser$s[idx] else 1.0
    }
    return(predicted.surv)
  }

  if (!is.numeric(estimated.strata) || any(estimated.strata < 0))
    stop("'estimated.strata' must be a non-negative integer vector.")
  if (sum(estimated.strata) != length(estimated.time))
    stop("sum(estimated.strata) must equal length(estimated.time).")
  K <- length(estimated.strata)
  if (length(strata.levels) != K)
    stop("'strata.levels' must have length K = length(estimated.strata).")

  if (length(predicted.strata) == 1L) {
    predicted.strata <- rep(predicted.strata, n_pred)
  } else if (length(predicted.strata) != n_pred) {
    stop("Length of predicted.strata must be 1 or match length(predicted.time).")
  }

  mapped <- if (is.factor(predicted.strata)) {
    match(as.character(predicted.strata), as.character(strata.levels))
  } else {
    match(predicted.strata, strata.levels)
  }
  if (any(is.na(mapped))) {
    bad <- unique(predicted.strata[is.na(mapped)])
    stop("Some values in predicted.strata are not found in 'strata.levels': ",
         paste(bad, collapse = ", "))
  }

  strata_start <- c(1L, head(cumsum(estimated.strata), -1L) + 1L)
  strata_end   <- cumsum(estimated.strata)

  for (i in seq_len(n_pred)) {
    s <- mapped[i]
    if (estimated.strata[s] == 0L) { out[i] <- NA_real_; next }
    idx <- strata_start[s]:strata_end[s]
    ser <- prepareSeries(estimated.time[idx], estimated.surv[idx])
    if (!length(ser$t)) { out[i] <- 1.0; next }
    j <- findInterval(predicted.time[i], ser$t, left.open = TRUE)
    predicted.surv[i] <- if (j > 0L) ser$s[j] else 1.0
  }
  return(predicted.surv)
}
