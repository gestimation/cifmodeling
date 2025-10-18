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
#' For CIF display, set \code{type.y = "risk"}.
#'
#' @param formula A model formula specifying the outcome and (optionally) \code{strata()}.
#' @param data A data frame containing variables in \code{formula}.
#' @param weights Optional name of the weight variable in \code{data}. Weights must be nonnegative; strictly positive is recommended.
#' @param subset.condition Optional character expression to subset \code{data} before analysis.
#' @param na.action Function to handle missing values (default: \code{\link[stats]{na.omit}}).

#' @param outcome.type \code{"SURVIVAL"} (KM) or \code{"COMPETING-RISK"} (AJ).
#' @param code.event1 Integer code of the event of interest (default \code{1}).
#' @param code.event2 Integer code of the competing event (default \code{2}).
#' @param code.censoring Integer code of censoring (default \code{0}).

#' @param error Character specifying variance type used internally. For \code{"SURVIVAL"} typically \code{"greenwood"}.
#'   For \code{"COMPETING-RISK"} pass options supported by \code{calculateAalenDeltaSE()} (\code{"aalen"}, \code{"delta"}, \code{"none"}).
#' @param conf.type Character transformation for CI on the probability scale (default \code{"arcsine-square root"}).
#' @param conf.int numeric two-sided confidence level (default \code{0.95}).

#' @param addConfidenceInterval Logical add \code{add_confidence_interval()} to plot. It calls geom_ribbon() (default \code{TRUE}).
#' @param addRiskTable Logical add \code{add_risktable(risktable_stats="n.risk")} to plot (default \code{TRUE}).
#' @param addEstimateTable Logical add \code{add_risktable(risktable_stats="estimate (conf.low, conf.high)")} to plot (default \code{FALSE}).
#' @param addCensorMark Logical add \code{add_censor_mark()} to plot. It calls geom_point() (default \code{TRUE}).
#' @param shape.censor.mark Integer point shape for censor marks (default \code{3}).
#' @param size.censor.mark Numeric point size for censor marks (default \code{2}).
#' @param addCompetingRiskMark Logical add time marks to describe event2 specified by Event(), usually the competing events. It calls geom_point() (default \code{TRUE}).
#' @param shape.competing.risk.mark Integer point shape for competing-risk marks (default \code{16}).
#' @param size.competing.risk.mark Numeric point size for competing-risk marks (default \code{2}).
#' @param addIntercurrentEventMark Logical overlay user-specified time marks per strata calls geom_point() (default \code{TRUE}).
#' @param intercurrent.event.time Named list of numeric vectors (names must be mapped to strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks (default \code{1}).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks (default \code{2}).
#' @param addQuantileLine Logical add \code{add_quantile()} to plot. It calls geom_segment() (default \code{TRUE}).
#' @param quantile Numeric specify quantile for \code{add_quantile()} (default \code{0.5}).

#' @param type.y \code{NULL} (survival) or \code{"risk"} (display \code{1 - survival} i.e. CIF).
#' @param label.x Character x-axis labels (default \code{"Time"}).
#' @param label.y Character y-axis labels (default internally set to \code{"Survival" or \code{"Cumulative incidence"}).
#' @param label.strata Character vector of labels for strata.
#' @param limits.x Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,max(out_readSurv$t))}.
#' @param limits.y Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,1)}.
#' @param breaks.x Numeric vectors for axis breaks (default \code{NULL}).
#' @param breaks.y Numeric vectors for axis breaks (default \code{NULL}).
#' @param use_coord_cartesian Logical specify use of coord_cartesian() (default \code{FALSE}).

#' @param style Character plot theme controls (default \code{"CLASSIC"}).
#' @param font.family Character plot theme controls (e.g. "sans", "serif", and "mono". default \code{"sans"}).
#' @param font.size Integer plot theme controls (default \code{12}).
#' @param legend.position Character specify position of legend: \code{"top"}, \code{"right"}, \code{"bottom"}, \code{"left"}, or \code{"none"} (default \code{"top"}).
#' @param report.survfit.std.err If \code{TRUE}, report SE on the log-survival scale (survfit's convention). Otherwise SE is on the probability scale.
#' @param return.ggsurvfit If \code{TRUE} (default), return a \pkg{ggsurvfit} object.
#' @param print.ggsurvfit If \code{TRUE} (default), draw a \pkg{ggsurvfit} graph
#' @param filename.ggsurvfit Character save the \pkg{ggsurvfit} graph with the name specified.

#' @returns A \code{survfit} object. For \code{outcome.type="SURVIVAL"}, \code{$surv} is the survival function.
#' For \code{outcome.type="COMPETING-RISK"}, \code{$surv} equals \code{1 - CIF} for \code{code.event1}.
#' Standard error and CIs are provided per \code{conf.type}. Note that some methods for \code{survfit} (e.g., \code{residuals.survfit}) may not be supported.
#'
#' @seealso \code{\link{polyreg}} for log-odds product models of CIFs; \pkg{ggsurvfit} for plotting helpers.

#' @examples
#' data(diabetes.complications)
#' survfit_by_group <- cifcurve(Event(t,epsilon) ~ fruitq, data = diabetes.complications,
#'                     outcome.type='COMPETING-RISK', error='delta', type.y = 'risk',
#'                     label.y = 'CIF of diabetic retinopathy', label.x = 'Years from registration')
#'
#' @importFrom ggsurvfit ggsurvfit add_confidence_interval add_risktable add_risktable_strata_symbol add_censor_mark add_quantile
#' @importFrom ggplot2 theme_classic theme_bw element_text labs lims geom_point aes ggsave scale_color_discrete scale_fill_discrete element_text element_rect element_blank scale_color_manual scale_fill_manual scale_linetype_manual scale_shape_manual
#' @importFrom Rcpp sourceCpp
#' @useDynLib cifmodeling, .registration = TRUE
#' @export
cifcurve <- function(formula,
                     data,
                     weights = NULL,
                     subset.condition = NULL,
                     na.action = na.omit,
                     outcome.type = "SURVIVAL",
                     code.event1 = 1,
                     code.event2 = 2,
                     code.censoring = 0,
                     error = NULL,
                     conf.type = "arcsine-square root",
                     conf.int = 0.95,
                     addConfidenceInterval = FALSE,
                     addRiskTable = TRUE,
                     addEstimateTable = FALSE,
                     addCensorMark = TRUE,
                     shape.censor.mark = 3,
                     size.censor.mark = 2,
                     addCompetingRiskMark = FALSE,
                     shape.competing.risk.mark = 16,
                     size.competing.risk.mark = 1,
                     addIntercurrentEventMark = FALSE,
                     intercurrent.event.time = NULL,
                     shape.intercurrent.event.mark = 1,
                     size.intercurrent.event.mark = 2,
                     addQuantileLine = FALSE,
                     quantile = 0.5,
                     type.y = NULL,
                     label.x = "Time",
                     label.y = NULL,
                     label.strata = NULL,
                     limits.x = NULL,
                     limits.y = NULL,
                     breaks.x = NULL,
                     breaks.y = NULL,
                     use_coord_cartesian = FALSE,
                     style = "CLASSIC",
                     font.family = "sans",
                     font.size = 12,
                     legend.position = "top",
                     report.survfit.std.err = FALSE,
                     return.ggsurvfit = FALSE,
                     print.ggsurvfit = TRUE,
                     filename.ggsurvfit = NULL) {

  style         <- check_style(style)
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
      time      = out_km$time,
      surv      = out_km$surv,
      n         = out_km$n,
      n.risk    = out_km$n.risk,
      n.event   = out_km$n.event,
      n.censor  = out_km$n.censor,
      std.err   = out_km$std.err,
      upper     = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$upper,
      lower     = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$lower,
      conf.type = conf.type,
      call      = match.call(),
      type      = "kaplan-meier",
      method    = "Kaplan-Meier"
    )
    if (any(as.integer(out_readSurv$strata) != 1)) {
      names(out_km$strata) <- levels(as.factor(out_readSurv$strata))
      survfit_object$strata <- out_km$strata
    }
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
      time        = out_aj$time1,
      surv        = out_aj$surv,
      n           = n,
      n.risk      = n.risk,
      n.event     = out_aj$n.event1,
      n.censor    = out_aj$n.censor,
      std.err     = out_aj$std.err,
      upper       = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$upper,
      lower       = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$lower,
      conf.type   = conf.type,
      call        = match.call(),
      type        = "aalen-johansen",
      method      = "aalen-johansen"
    )
    if (any(as.integer(out_readSurv$strata) != 1)) survfit_object$strata <- out_aj$strata1
    class(survfit_object) <- "survfit"
  }

  if (isTRUE(print.ggsurvfit) || isTRUE(return.ggsurvfit) || !is.null(filename.ggsurvfit)) {
    out_ggsurvfit <- call_ggsurvfit(
      survfit_object                = survfit_object,
      out_readSurv                  = out_readSurv,
      type.y                        = type.y,
      conf.type                     = conf.type,
      addConfidenceInterval         = addConfidenceInterval,
      addRiskTable                  = addRiskTable,
      addEstimateTable              = addEstimateTable,
      addCensorMark                 = addCensorMark,
      shape.censor.mark             = shape.censor.mark,
      size.censor.mark              = size.censor.mark,
      addCompetingRiskMark          = addCompetingRiskMark,
      competing.risk.time           = out_mcm,
      shape.competing.risk.mark     = shape.competing.risk.mark,
      size.competing.risk.mark      = size.competing.risk.mark,
      addIntercurrentEventMark      = addIntercurrentEventMark,
      intercurrent.event.time       = intercurrent.event.time,
      shape.intercurrent.event.mark = shape.intercurrent.event.mark,
      size.intercurrent.event.mark  = size.intercurrent.event.mark,
      addQuantileLine               = addQuantileLine,
      quantile                      = quantile,
      label.x                       = label.x,
      label.y                       = label.y,
      label.strata                  = lab_map,
      limits.x = limits.x, limits.y = limits.y,
      breaks.x = breaks.x, breaks.y = breaks.y,
      use_coord_cartesian           = use_coord_cartesian,
      style                         = style,
      font.family                   = font.family,
      font.size                     = font.size,
      legend.position               = legend.position
    )
    if (isTRUE(print.ggsurvfit)) print(out_ggsurvfit)
    if (isTRUE(return.ggsurvfit)) return(out_ggsurvfit)
    if (!is.null(filename.ggsurvfit)) ggsave(filename.ggsurvfit, plot = out_ggsurvfit, width = 6, height = 6, dpi = 300)
  }
  return(survfit_object)
}

#' Plot survival or cumulative incidence curves with ggsurvfit
#'
#' @param survfit_object A \code{survfit} object.
#' @param out_readSurv (optional) List returned by your \code{readSurv()} to auto-set x limits.
#' @param conf.type Character transformation for CI on the probability scale (default \code{"arcsine-square root"}).

#' @param addConfidenceInterval Logical add \code{add_confidence_interval()} to plot. It calls geom_ribbon() (default \code{TRUE}).
#' @param addRiskTable Logical add \code{add_risktable(risktable_stats="n.risk")} to plot (default \code{TRUE}).
#' @param addEstimateTable Logical add \code{add_risktable(risktable_stats="estimate (conf.low, conf.high)")} to plot (default \code{FALSE}).
#' @param addCensorMark Logical add \code{add_censor_mark()} to plot. It calls geom_point() (default \code{TRUE}).
#' @param shape.censor.mark Integer point shape for censor marks (default \code{3}).
#' @param size.censor.mark Numeric point size for censor marks (default \code{2}).
#' @param addCompetingRiskMark Logical add time marks to describe event2 specified by Event(), usually the competing events. It calls geom_point() (default \code{TRUE}).
#' @param shape.competing.risk.mark Integer point shape for competing-risk marks (default \code{16}).
#' @param size.competing.risk.mark Numeric point size for competing-risk marks (default \code{2}).
#' @param addIntercurrentEventMark Logical overlay user-specified time marks per strata calls geom_point() (default \code{TRUE}).
#' @param intercurrent.event.time Named list of numeric vectors (names must be mapped to strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks (default \code{1}).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks (default \code{2}).
#' @param addQuantileLine Logical add \code{add_quantile()} to plot. It calls geom_segment() (default \code{TRUE}).
#' @param quantile Numeric specify quantile for \code{add_quantile()} (default \code{0.5}).

#' @param type.y \code{NULL} (survival) or \code{"risk"} (display \code{1 - survival} i.e. CIF).
#' @param label.x Character x-axis labels (default \code{"Time"}).
#' @param label.y Character y-axis labels (default internally set to \code{"Survival" or \code{"Cumulative incidence"}).
#' @param label.strata Character vector of labels for strata.

#' @param limits.x Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,max(out_readSurv$t))}.
#' @param limits.y Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,1)}.
#' @param breaks.x Numeric vectors for axis breaks (default \code{NULL}).
#' @param breaks.y Numeric vectors for axis breaks (default \code{NULL}).
#' @param use_coord_cartesian Logical specify use of coord_cartesian() (default \code{FALSE}).

#' @param style Character plot theme controls (default \code{"ggsurvfit"}).
#' @param font.family Character plot theme controls (default \code{"sans"}).
#' @param font.size Integer plot theme controls (default \code{14}).
#' @param legend.position Character specify position of legend: \code{"top"}, \code{"right"}, \code{"bottom"}, \code{"left"}, or \code{"none"} (default \code{"top"}).

#' @return A \code{ggplot} object.
#' Plot survival or cumulative incidence curves with ggsurvfit
#' ...
#' ...
call_ggsurvfit <- function(
    survfit_object,
    out_readSurv = NULL,
    conf.type = NULL,
    addConfidenceInterval = TRUE,
    addRiskTable = TRUE,
    addEstimateTable = FALSE,
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
    addQuantileLine = FALSE,
    quantile = 0.5,
    type.y = NULL,
    label.x = "Time",
    label.y = NULL,
    label.strata = NULL,
    limits.x = NULL,
    limits.y = NULL,
    breaks.x = NULL,
    breaks.y = NULL,
    use_coord_cartesian = FALSE,
    style = "CLASSIC",
    font.family = "sans",
    font.size = 14,
    legend.position = "top"
){
  out_cg <- check_ggsurvfit(
    survfit_object                = survfit_object,
    type.y                        = type.y,
    conf.type                     = conf.type,
    label.y                       = label.y,
    limits.x                      = limits.x,
    limits.y                      = limits.y,
    breaks.x                      = breaks.x,
    breaks.y                      = breaks.y,
    addConfidenceInterval         = addConfidenceInterval,
    addCensorMark                 = addCensorMark,
    addCompetingRiskMark          = addCompetingRiskMark,
    addIntercurrentEventMark      = addIntercurrentEventMark,
    shape.censor.mark             = shape.censor.mark,
    shape.competing.risk.mark     = shape.competing.risk.mark,
    shape.intercurrent.event.mark = shape.intercurrent.event.mark,
    out_readSurv                  = out_readSurv,
    use_coord_cartesian           = use_coord_cartesian,
    style                         = style
  )

  p <- out_cg$out_survfit_object +
    ggplot2::labs(x = label.x, y = out_cg$label.y)

  if (!is.null(label.strata)) {
    p <- p +
      ggplot2::scale_color_discrete(breaks = names(label.strata), labels = unname(label.strata)) +
      ggplot2::scale_fill_discrete(breaks = names(label.strata), labels = unname(label.strata))
  }

  if (isTRUE(addConfidenceInterval)) p <- p + add_confidence_interval()
  if (isTRUE(addCensorMark))         p <- p + add_censor_mark(shape = shape.censor.mark, size = size.censor.mark)
  if (isTRUE(addQuantileLine))       p <- p + ggsurvfit::add_quantile(y_value=quantile)

  if (isTRUE(addEstimateTable) && isTRUE(addRiskTable)) {
    p <- p + add_risktable(risktable_stats = c("n.risk", "{round(estimate, digit=2)} ({round(conf.low, digit=2)}, {round(conf.high, digit=2)})"), stats_label=c("No. at risk", "Estimate (95% CI)"), theme=theme_risktable_font(font.family=font.family, plot.title.size=font.size), risktable_group = "risktable_stats")
    if (!identical(style, "MONOCHROME")) {
      p <- p + ggsurvfit::add_risktable_strata_symbol()
    }
  } else if (isTRUE(addEstimateTable)) {
    p <- p + add_risktable(risktable_stats = c("{round(estimate, digit=2)} ({round(conf.low, digit=2)}, {round(conf.high, digit=2)})"), stats_label=c("Estimate (95% CI)"), theme=theme_risktable_font(font.family=font.family, plot.title.size=font.size), )
    if (!identical(style, "MONOCHROME")) {
      p <- p + ggsurvfit::add_risktable_strata_symbol()
    }
  } else if (isTRUE(addRiskTable)) {
    p <- p + add_risktable(risktable_stats = c("n.risk"), stats_label=c("No. at risk"), theme=theme_risktable_font(font.family=font.family, plot.title.size=font.size))
    if (!identical(style, "MONOCHROME")) {
      p <- p + ggsurvfit::add_risktable_strata_symbol()
    }
  }
  if (isTRUE(addCompetingRiskMark) && length(competing.risk.time)) {
    p <- drawMarks(p, survfit_object, competing.risk.time, out_cg$type.y, shape = shape.competing.risk.mark, size = size.competing.risk.mark)
  }
  if (isTRUE(addIntercurrentEventMark) && length(intercurrent.event.time)) {
    p <- drawMarks(p, survfit_object, intercurrent.event.time, out_cg$type.y, shape = shape.intercurrent.event.mark, size = size.intercurrent.event.mark)
  }

  if (isTRUE(use_coord_cartesian)) {
    if (!is.null(breaks.x)) p <- p + ggplot2::scale_x_continuous(breaks = breaks.x)
    if (!is.null(breaks.y)) p <- p + ggplot2::scale_y_continuous(breaks = breaks.y)
    if (!is.null(limits.x) || !is.null(limits.y)) {
      p <- p + ggplot2::coord_cartesian(xlim = limits.x, ylim = limits.y, expand = FALSE)
    } else if (!is.null(limits.x)) {
      p <- p + ggplot2::coord_cartesian(xlim = c(0, max(out_readSurv$t)), ylim = limits.y, expand = FALSE)
    } else if (!is.null(limits.y)) {
      p <- p + ggplot2::coord_cartesian(xlim = limits.x, ylim = c(0, 1), expand = FALSE)
    } else {
      p <- p + ggplot2::coord_cartesian(xlim = c(0, max(out_readSurv$t)), ylim = c(0, 1), expand = FALSE)
    }
  } else {
    if (!is.null(breaks.y)) {
      p <- p + ggplot2::scale_y_continuous(breaks = breaks.y)
    } else if (!is.null(limits.y)) {
      p <- p + ggplot2::lims(y = limits.y)
    } else {
      p <- p + ggplot2::lims(y = c(0, 1))
    }
    if (!is.null(breaks.x)) {
      p <- p + ggplot2::scale_x_continuous(breaks = breaks.x)
    } else if (!is.null(limits.x)) {
      p <- p + ggplot2::lims(x = limits.x)
    } else {
      p <- p + ggplot2::lims(x = c(0, max(out_readSurv$t)))
    }
  }

  if (!identical(style, "GGSURVFIT")) {
    p <- applyStyle(p, style = style, font.family, font.size, legend.position)
  }

  use.polyreg <- FALSE
  time.point.polyreg <- NULL
  if (use.polyreg==TRUE & !is.null(time.point.polyreg) & length(levels(out_readSurv$strata))==2) {
    risk_ratio <- vector("numeric", length = length(time.point.polyreg))
    for (i in 1:length(time.point.polyreg)) {
      out_polyreg <- polyreg(nuisance.model = update.formula(formula, ~1), exposure = out_readSurv$strata_name, data = data, time.point = time.point.polyreg[i], outcome.type='SURVIVAL')
      risk_ratio[i] <- create_rr_text(out_polyreg$coefficient, out_polyreg$cov, 2)
      p <- p + geom_vline(xintercept = time.point.polyreg[i], linetype = "dashed", color = "black", size = 1)
    }
    if (text.position.polyreg=="bottom") {
      annotations <- data.frame(x=time.point.polyreg*1.02, y=rep(0, length(time.point.polyreg)), label = risk_ratio)
      p <- p + geom_text(data=annotations, aes(x=x, y=y, label=label), size = 4, color = "black", hjust = 0, vjust = 0)
    }
    if (text.position.polyreg=="top") {
      annotations <- data.frame(x=time.point.polyreg*1.02, y=rep(0.9, length(time.point.polyreg)), label = risk_ratio)
      p <- p + geom_text(data=annotations, aes(x=x, y=y, label=label), size = 4, color = "black", hjust = 0, vjust = 0)
    }
  }
  return(p)
}

style_theme_classic <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_classic(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 4, family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size + 4, family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_blank()
    )
}

style_theme_bold <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_classic(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 3, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size,     family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_blank()
    )
}

style_theme_framed <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_bw(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 3, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size,     family = font.family),
      legend.background = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank(),
      strip.background  = ggplot2::element_rect(fill = "grey90", color = "black"),
      panel.border      = ggplot2::element_rect(color = "black", linewidth = 0.8)
    )
}

style_theme_monochrome <- function(font.family = "sans", font.size = 14, legend.position = "top") {
  ggplot2::theme_classic(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 1, face = "bold"),
      axis.text         = ggplot2::element_text(size = font.size, color = "grey20"),
      legend.text       = ggplot2::element_text(size = font.size - 1, color = "grey20"),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_blank()
    )
}

scale_monochrome <- function(n_strata = 6) {
  ltys_all   <- c("dashed","solid","dotted","longdash","dotdash","twodash","dashed","solid","dotted","longdash","dotdash","twodash")
  shapes_all <- c(16, 1, 3, 4, 15, 17, 16, 1, 3, 4, 15, 17)
  n_use <- min(n_strata, length(ltys_all), length(shapes_all))
  list(
    ggplot2::scale_color_manual(values = rep("black", n_use)),
    ggplot2::scale_fill_manual(values = gray(seq(0.85, 0.30, length.out = n_use))),
    ggplot2::scale_linetype_manual(values = ltys_all[seq_len(n_use)]),
    ggplot2::scale_shape_manual(values = shapes_all[seq_len(n_use)])
  )
}

applyStyle <- function(
    p,
    style = c("CLASSIC", "BOLD", "FRAMED", "MONOCHROME"),
    font.family = "sans",
    font.size = 14,
    legend.position = "top",
    n_strata = 6
) {
  style <- match.arg(style)
  style_theme <- switch(
    style,
    CLASSIC    = style_theme_classic(font.family, font.size, legend.position),
    BOLD       = style_theme_bold(font.family, font.size, legend.position),
    FRAMED     = style_theme_framed(font.family, font.size, legend.position),
    MONOCHROME = style_theme_monochrome(font.family, font.size, legend.position)
  )
  p <- p + style_theme
  if (identical(style, "MONOCHROME")) {
    p <- p + scale_monochrome(n_strata = n_strata)
  }
  return(p)
}

drawMarks <- function(p, survfit_object, marks, type.y, shape, size) {
  if (is.null(marks) || !length(marks)) return(p)
  mark_df <- makeMarkDataFrame(survfit_object, marks, type.y, extend = TRUE)
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

theme_risktable_font <- function(
    axis.text.y.size = 10,
    plot.title.size = 10.75,
    font.family = "sans"
) {
  list(
    ggplot2::theme_bw(base_family = font.family),
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(
        size = 9, vjust = 1, hjust = 1, family = font.family
      ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_text(
        size = axis.text.y.size, colour = "black", face = "plain", family = font.family
      ),
      plot.margin = ggplot2::unit(c(0, 5.5, 0, 5.5), "points"),
      plot.title  = ggplot2::element_text(
        hjust = 0, vjust = 0, size = plot.title.size, family = font.family
      ),
      legend.position = "none"
    ),
    ggplot2::xlab(NULL),
    ggplot2::ylab(NULL)
  )
}

alignMarkKeys <- function(fit, marks) {
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

makeMarkDataFrame <- function(fit, marks, type.y, extend = TRUE) {
  if (is.null(marks) || length(marks) == 0) return(NULL)
  out_alignMarkKeys <- alignMarkKeys(fit, marks)
  strata_names <- if (is.null(fit$strata)) "(all)" else names(fit$strata)

  out <- lapply(names(out_alignMarkKeys), function(st_lab) {
    tt <- out_alignMarkKeys[[st_lab]]
    if (length(tt) == 0) return(NULL)
    fit_st <- if (is.null(fit$strata)) {
      if (!identical(st_lab, "(all)")) return(NULL)
      fit
    } else {
      i <- match(st_lab, strata_names); if (is.na(i)) return(NULL)
      fit[i]
    }
    sm <- summary(fit_st, times = tt, extend = extend)
    y  <- if (identical(type.y, "risk")) 1 - sm$surv else sm$surv
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
