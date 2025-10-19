#' @title Visualize time-to-event outcomes and intercurrent events
#' @description
#' Draw a publication-ready plot. Accepts a \code{survfit} object or a \code{formula+data}
#' (in which case it computes a \code{survfit} via \code{cifcurve()} first).
#' @details
#' \strong{Plotting:}
#' This function calls an internal helper \code{call_ggsurvfit()} which adds confidence bands,
#' risk table, censoring marks, and optional competing-risk and intercurrent-event marks.
#' For CIF display, set \code{type.y = "risk"}.
#'
#' @param x A model formula or a survfit object.
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
#' @param type.y \code{NULL} (survival) or \code{"risk"} (display \code{1 - survival} i.e. CIF).
#' @param label.x Character x-axis labels (default \code{"Time"}).
#' @param label.y Character y-axis labels (default internally set to \code{"Survival"} or \code{"Cumulative incidence"}).
#' @param label.strata Character vector of labels for strata.
#' @param limits.x Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,max(out_readSurv$t))}.
#' @param limits.y Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,1)}.
#' @param breaks.x Numeric vectors for axis breaks (default \code{NULL}).
#' @param breaks.y Numeric vectors for axis breaks (default \code{NULL}).
#' @param use_coord_cartesian Logical specify use of coord_cartesian() (default \code{FALSE}).
#'
#' @param addConfidenceInterval Logical add \code{add_confidence_interval()} to plot. It calls geom_ribbon() (default \code{TRUE}).
#' @param addRiskTable Logical add \code{add_risktable(risktable_stats="n.risk")} to plot (default \code{TRUE}).
#' @param addEstimateTable Logical add \code{add_risktable(risktable_stats="estimate (conf.low, conf.high)")} to plot (default \code{FALSE}).
#' @param addCensorMark Logical add \code{add_censor_mark()} to plot. It calls geom_point() (default \code{TRUE}).
#' @param shape.censor.mark Integer point shape for censor marks (default \code{3}).
#' @param size.censor.mark Numeric point size for censor marks (default \code{2}).
#' @param addCompetingRiskMark Logical add time marks to describe event2 specified by Event(), usually the competing events.
#' It calls geom_point() (default \code{TRUE}).
#' @param shape.competing.risk.mark Integer point shape for competing-risk marks (default \code{16}).
#' @param size.competing.risk.mark Numeric point size for competing-risk marks (default \code{2}).
#' @param addIntercurrentEventMark Logical overlay user-specified time marks per strata calls geom_point() (default \code{TRUE}).
#' @param intercurrent.event.time Named list of numeric vectors (names must be mapped to strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks (default \code{1}).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks (default \code{2}).
#' @param addQuantileLine Logical add \code{add_quantile()} to plot. It calls geom_segment() (default \code{TRUE}).
#' @param quantile Numeric specify quantile for \code{add_quantile()} (default \code{0.5}).

#' @param style Character plot theme controls (default \code{"CLASSIC"}).
#' @param font.family Character plot theme controls (e.g. "sans", "serif", and "mono". default \code{"sans"}).
#' @param font.size Integer plot theme controls (default \code{12}).
#' @param legend.position Character specify position of legend: \code{"top"}, \code{"right"}, \code{"bottom"}, \code{"left"}, or \code{"none"} (default \code{"top"}).
#' @param filename.ggsave Character save the \pkg{ggsurvfit} plot with the path and name specified.
#' @param width.ggsave Numeric specify width of the \pkg{ggsurvfit} plot.
#' @param height.ggsave Numeric specify height of the \pkg{ggsurvfit} plot.
#' @param dpi.ggsave Numeric specify dpi of the \pkg{ggsurvfit} plot.

#' @return A \code{ggplot} object.
#' For \code{outcome.type="COMPETING-RISK"}, \code{$surv} equals \code{1 - CIF} for \code{code.event1}.
#' Standard error and CIs are provided per \code{conf.type}. Note that some methods for \code{survfit} (e.g., \code{residuals.survfit}) may not be supported.

#' @examples
#' data(diabetes.complications)
#' survfit_by_group <- cifcurve(Event(t,epsilon) ~ fruitq, data = diabetes.complications,
#'                     outcome.type='COMPETING-RISK', error='delta')
#' cifplot(survfit_by_group, type.y = 'risk', label.y = 'CIF of diabetic retinopathy', label.x = 'Years from registration', addConfidenceInterval=FALSE, style="MONOCHROME")
#' out_readEventTime <- readEventTime(Event(t,epsilon) ~ fruitq, data = diabetes.complications, which_event = "event2")
#' cifplot(survfit_by_group, type.y = 'risk', label.y = 'CIF of diabetic retinopathy', label.x = 'Years from registration',
#' addCompetingRiskMark=TRUE, competing.risk.time=out_readEventTime)

#' @importFrom ggsurvfit ggsurvfit add_confidence_interval add_risktable add_risktable_strata_symbol add_censor_mark add_quantile
#' @importFrom ggplot2 theme_classic theme_bw element_text labs lims geom_point aes ggsave scale_color_discrete scale_fill_discrete element_text element_rect element_blank scale_color_manual scale_fill_manual scale_linetype_manual scale_shape_manual
#' @seealso \code{\link{cifcurve}} for the estimators; \code{\link{cifpanel}} for display of multiple plots; \pkg{ggsurvfit} for plotting helpers; \code{\link{polyreg}} for log-odds product models of CIFs.

#' @export
cifplot <- function(
    x,
    data = NULL,
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
    type.y = NULL,
    label.x = "Time",
    label.y = NULL,
    label.strata = NULL,
    limits.x = NULL,
    limits.y = NULL,
    breaks.x = NULL,
    breaks.y = NULL,
    use_coord_cartesian = FALSE,
    addConfidenceInterval = TRUE,
    addRiskTable = TRUE,
    addEstimateTable = FALSE,
    addCensorMark = TRUE,
    shape.censor.mark = 3,
    size.censor.mark = 2,
    addCompetingRiskMark = FALSE,
    competing.risk.time = list(),
    shape.competing.risk.mark = 16,
    size.competing.risk.mark = 2,
    addIntercurrentEventMark = FALSE,
    intercurrent.event.time = list(),
    shape.intercurrent.event.mark = 1,
    size.intercurrent.event.mark = 2,
    addQuantileLine = FALSE,
    quantile = 0.5,
    style = "CLASSIC",
    font.family = "sans",
    font.size = 12,
    legend.position = "top",
    filename.ggsave = NULL,
    width.ggsave = 6,
    height.ggsave = 6,
    dpi.ggsave = 300
) {
  if (!inherits(x, "survfit")) {
    if (is.null(data)) stop("When `x` is a formula, `data` must be provided.")
    if (addCompetingRiskMark==TRUE && length(competing.risk.time)==0) {
      competing.risk.time <- readEventTime(x, data=data, which_event = "event2", code.event1=code.event1, code.event2=code.event2, code.censoring=code.censoring)
    }
    x <- cifcurve(x, data = data, weights=weights, subset.condition=subset.condition, na.action=na.action,
                  outcome.type=outcome.type, code.event1=code.event1, code.event2=code.event2, code.censoring=code.censoring,
                  error=error, conf.type=conf.type, conf.int=conf.int)
  }

  p <- call_ggsurvfit(
    survfit_object                = x,
    out_readSurv                  = NULL,
    conf.type                     = conf.type,
    addConfidenceInterval         = addConfidenceInterval,
    addRiskTable                  = addRiskTable,
    addEstimateTable              = addEstimateTable,
    addCensorMark                 = addCensorMark,
    shape.censor.mark             = shape.censor.mark,
    size.censor.mark              = size.censor.mark,
    addCompetingRiskMark          = addCompetingRiskMark,
    competing.risk.time           = competing.risk.time,
    shape.competing.risk.mark     = shape.competing.risk.mark,
    size.competing.risk.mark      = size.competing.risk.mark,
    addIntercurrentEventMark      = addIntercurrentEventMark,
    intercurrent.event.time       = intercurrent.event.time,
    shape.intercurrent.event.mark = shape.intercurrent.event.mark,
    size.intercurrent.event.mark  = size.intercurrent.event.mark,
    addQuantileLine               = addQuantileLine,
    quantile                      = quantile,
    type.y                        = type.y,
    label.x                       = label.x,
    label.y                       = label.y,
    label.strata                  = label.strata,
    limits.x                      = limits.x,
    limits.y                      = limits.y,
    breaks.x                      = breaks.x,
    breaks.y                      = breaks.y,
    use_coord_cartesian           = use_coord_cartesian,
    style                         = style,
    font.family                   = font.family,
    font.size                     = font.size,
    legend.position               = legend.position
  )
  if (!is.null(filename.ggsave)) ggplot2::ggsave(filename.ggsave, plot = p, width = width.ggsave, height = height.ggsave, dpi.ggsave = dpi)
  return(p)
}

#' Extract per-stratum event times from (formula, data)
#'
#' @description
#' Build per-stratum event-time lists that you can pass to `competing.risk.time` or
#' `intercurrent.event.time` in `cifplot()` / `cifpanel()`.
#'
#' @param formula Event(time, status) ~ strata(...)
#' @param data Data frame.
#' @param which One of "event2", "event1", "censor", "code".
#' @param event.code Used when which="code": status value to pick (e.g., 3 for an intercurrent event).
#' @param code.event1,code.event2,code.censoring Integers describing your coding (used to build d1/d2/d0 inside readSurv).
#' @param subset.condition Optional character expression to subset `data`.
#' @param na.action Function to handle NA (default na.omit).
#' @param unique_times,drop_empty See `cif_event_times_from_frame()`.
#' @return Named list of numeric vectors (times per stratum).
#' @export
readEventTime <- function(
    formula, data, which_event = c("event2", "event1", "censor", "censoring", "user_specified"),
    code.event1 = 1, code.event2 = 2, code.censoring = 0, user_specified_code = NULL,
    subset.condition = NULL, na.action = stats::na.omit,
    unique_times = TRUE, drop_empty = TRUE
){
  which_event <- match.arg(which_event)
  out_readSurv <- readSurv(
    formula = formula, data = data, weights = NULL,
    code.event1 = code.event1, code.event2 = code.event2, code.censoring = code.censoring,
    subset.condition = subset.condition, na.action = na.action
  )
  getEventTime(
    out_readSurv = out_readSurv,
    which_event = which_event, user_specified_code = user_specified_code,
    unique_times = unique_times, drop_empty = drop_empty
  )
}

getEventTime <- function(
    out_readSurv,
    which_event = c("event2", "event1", "censor", "censoring", "user_specified"),
    user_specified_code = NULL,
    unique_times = TRUE,
    drop_empty = TRUE
){
  which_event <- match.arg(which_event)

  if (is.null(out_readSurv$t) || is.null(out_readSurv$strata)) {
    stop("out_readSurv must contain $t and $strata.")
  }
  tvec   <- out_readSurv$t
  strata <- as.character(out_readSurv$strata)

  pick <- switch(
    which_event,
    event1      = out_readSurv$d1,
    event2      = out_readSurv$d2,
    censor      = out_readSurv$d0,
    censoring   = out_readSurv$d0,
    user_specified = {
      epsilon <- out_readSurv$epsilon
      if (is.null(epsilon)) stop("out_readSurv$epsilon is required when which_event='user_specified'.")
      if (is.null(user_specified_code)) stop("Provide user_specified_code when which_event='user_specified'.")
      as.integer(epsilon == user_specified_code)
    }
  )

  if (is.null(pick)) stop("Requested 'which_event' requires a component not found in out_readSurv.")

  labs <- unique(strata)
  out  <- setNames(vector("list", length(labs)), labs)
  for (s in labs) {
    idx <- (strata == s) & (pick > 0L)
    tt  <- tvec[idx]
    if (unique_times) tt <- sort(unique(tt))
    if (length(tt) > 0L || !drop_empty) out[[s]] <- tt
  }
  if (drop_empty) out <- out[vapply(out, length, integer(1)) > 0L]
  return(out)
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
#' @param label.y Character y-axis labels (default internally set to \code{"Survival"} or \code{"Cumulative incidence"}).
#' @param label.strata Character vector of labels for strata.

#' @param limits.x Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,max(out_readSurv$t))}.
#' @param limits.y Numeric length-2 vectors for axis limits. If NULL it is internally set to \code{c(0,1)}.
#' @param breaks.x Numeric vectors for axis breaks (default \code{NULL}).
#' @param breaks.y Numeric vectors for axis breaks (default \code{NULL}).
#' @param use_coord_cartesian Logical specify use of coord_cartesian() (default \code{FALSE}).

#' @param style Character plot theme controls (default \code{"CLASSIC"}).
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

  .xmax_from_fit <- function(sf) {
    if (!is.null(sf$time)) {
      xm <- suppressWarnings(max(sf$time, na.rm = TRUE))
      if (is.finite(xm)) return(xm)
    }
    sm <- try(suppressWarnings(summary(sf)), silent = TRUE)
    if (!inherits(sm, "try-error") && !is.null(sm$time)) {
      xm <- suppressWarnings(max(sm$time, na.rm = TRUE))
      if (is.finite(xm)) return(xm)
    }
    return(1)
  }
  x_max <- .xmax_from_fit(survfit_object)

  out_cg <- check_ggsurvfit(
    survfit_object = survfit_object,
    type.y = type.y,
    conf.type = conf.type,
    label.y = label.y,
    limits.x = limits.x,
    limits.y = limits.y,
    breaks.x = breaks.x,
    breaks.y = breaks.y,
    addConfidenceInterval = addConfidenceInterval,
    addCensorMark = addCensorMark,
    addCompetingRiskMark = addCompetingRiskMark,
    addIntercurrentEventMark = addIntercurrentEventMark,
    shape.censor.mark = shape.censor.mark,
    shape.competing.risk.mark = shape.competing.risk.mark,
    shape.intercurrent.event.mark = shape.intercurrent.event.mark,
    out_readSurv = out_readSurv,
    use_coord_cartesian = use_coord_cartesian,
    style = style
  )

  lab_map <- .normalize_label_strata(survfit_object, label.strata)

  p <- out_cg$out_survfit_object +
    ggplot2::labs(x = label.x, y = out_cg$label.y)

  if (!is.null(lab_map)) {
    p <- p +
      ggplot2::scale_color_discrete(breaks = names(lab_map), labels = unname(lab_map)) +
      ggplot2::scale_fill_discrete(breaks = names(lab_map), labels = unname(lab_map))
  }

  if (isTRUE(addConfidenceInterval)) p <- p + add_confidence_interval()
  if (isTRUE(addCensorMark))         p <- p + add_censor_mark(shape = shape.censor.mark, size = size.censor.mark)
  if (isTRUE(addQuantileLine))       p <- p + ggsurvfit::add_quantile(y_value=quantile)

  if (isTRUE(addEstimateTable) && isTRUE(addRiskTable)) {
    p <- p + add_risktable(risktable_stats = c("n.risk", "{round(estimate, digits=2)} ({round(conf.low, digits=2)}, {round(conf.high, digits=2)})"), stats_label=c("No. at risk", "Estimate (95% CI)"), theme=theme_risktable_font(font.family=font.family, plot.title.size=font.size), risktable_group = "risktable_stats")
    if (!identical(style, "MONOCHROME")) {
      p <- p + ggsurvfit::add_risktable_strata_symbol()
    }
  } else if (isTRUE(addEstimateTable)) {
    p <- p + add_risktable(risktable_stats = c("{round(estimate, digits=2)} ({round(conf.low, digits=2)}, {round(conf.high, digits=2)})"), stats_label=c("Estimate (95% CI)"), theme=theme_risktable_font(font.family=font.family, plot.title.size=font.size), )
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

  x_max <- .xmax_from_fit(survfit_object)
  if (isTRUE(use_coord_cartesian)) {
    if (!is.null(breaks.x)) p <- p + ggplot2::scale_x_continuous(breaks = breaks.x)
    if (!is.null(breaks.y)) p <- p + ggplot2::scale_y_continuous(breaks = breaks.y)
    if (!is.null(limits.x) || !is.null(limits.y)) {
      p <- p + ggplot2::coord_cartesian(xlim = limits.x, ylim = limits.y, expand = FALSE)
    } else if (!is.null(limits.x)) {
      p <- p + ggplot2::coord_cartesian(xlim = limits.x, ylim = limits.y, expand = FALSE)
    } else if (!is.null(limits.y)) {
      p <- p + ggplot2::coord_cartesian(xlim = c(0, x_max), ylim = limits.y, expand = FALSE)
    } else {
      p <- p + ggplot2::coord_cartesian(xlim = c(0, x_max), ylim = c(0, 1), expand = FALSE)
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
      p <- p + ggplot2::lims(x = c(0, x_max))
    }
  }

  if (!identical(style, "GGSURVFIT")) {
    p <- applyStyle(p, style = style, font.family, font.size, legend.position)
  }

  if (!is.null(lab_map)) {
    lvls   <- names(lab_map)
    labs   <- unname(lab_map)
    n      <- length(lvls)

    if (identical(style, "MONOCHROME")) {
      ltys_all   <- c("dashed","solid","dotted","longdash","dotdash","twodash",
                      "dashed","solid","dotted","longdash","dotdash","twodash")
      shapes_all <- c(16, 1, 3, 4, 15, 17, 16, 1, 3, 4, 15, 17)
      ltys   <- ltys_all[seq_len(n)]
      shps   <- shapes_all[seq_len(n)]
      fills  <- gray(seq(0.85, 0.30, length.out = n))

      p <- p +
        ggplot2::scale_color_manual(values = rep("black", n), limits = lvls, labels = labs) +
        ggplot2::scale_fill_manual (values = fills,                limits = lvls, labels = labs) +
        ggplot2::scale_linetype_manual(values = ltys,              limits = lvls, labels = labs) +
        ggplot2::scale_shape_manual   (values = shps,              limits = lvls, labels = labs)

    } else {
      p <- p +
        ggplot2::scale_color_discrete  (limits = lvls, labels = labs) +
        ggplot2::scale_fill_discrete   (limits = lvls, labels = labs) +
        ggplot2::scale_linetype_discrete(limits = lvls, labels = labs) +
        ggplot2::scale_shape_discrete  (limits = lvls, labels = labs)
    }
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

check_ggsurvfit <- function(
    survfit_object,
    type.y,
    conf.type,
    label.y,
    limits.x, limits.y,
    breaks.x,  breaks.y,
    addConfidenceInterval,
    addCensorMark,
    addCompetingRiskMark,
    addIntercurrentEventMark,
    shape.censor.mark,
    shape.competing.risk.mark,
    shape.intercurrent.event.mark,
    out_readSurv,
    use_coord_cartesian,
    style
){
  if (isTRUE(addCensorMark) && isTRUE(addIntercurrentEventMark) &&
      identical(shape.censor.mark, shape.intercurrent.event.mark)) {
    warning("shape.censor.mark and shape.intercurrent.event.mark specify an identical type of symbol", call. = FALSE)
  }
  if (isTRUE(addCensorMark) && isTRUE(addCompetingRiskMark) &&
      identical(shape.censor.mark, shape.competing.risk.mark)) {
    warning("shape.censor.mark and shape.competing.risk.mark specify an identical type of symbol", call. = FALSE)
  }
  if (isTRUE(addCompetingRiskMark) && isTRUE(addIntercurrentEventMark) &&
      identical(shape.competing.risk.mark, shape.intercurrent.event.mark)) {
    warning("shape.competing.risk.mark and shape.intercurrent.event.mark specify an identical type of symbol", call. = FALSE)
  }

  is_len2_num <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x))
  is_increasing <- function(x) all(diff(x) > 0, na.rm = TRUE)
  is_nondec <- function(x) all(diff(x) >= 0, na.rm = TRUE)

  if (!is.null(limits.x)) {
    if (!is_len2_num(limits.x)) {
      warning("`limits.x` must be a numeric length-2 vector.", call. = FALSE)
    } else if (!is_increasing(limits.x)) {
      warning("`limits.x` must be strictly increasing, e.g., c(xmin, xmax).", call. = FALSE)
    }
    tmax <- suppressWarnings(max(survfit_object$time, na.rm = TRUE))
    if (is.finite(tmax) && is_len2_num(limits.x)) {
      if (tmax < limits.x[1] || tmax > limits.x[2]) {
        warning(sprintf("Max of `survfit_object$time` (%.4g) lies outside `limits.x` = [%.4g, %.4g].",
                        tmax, limits.x[1], limits.x[2]), call. = FALSE)
      }
    }
  } else if (!is.null(out_readSurv) && !is.null(out_readSurv$t)) {
    tmax <- suppressWarnings(max(out_readSurv$t, na.rm = TRUE))
    if (!is.finite(tmax) || tmax <= 0) {
      warning("`out_readSurv$t` provided but has no finite positive max; x-limits fallback may fail.", call. = FALSE)
    }
  }

  if (!is.null(limits.y)) {
    if (!is_len2_num(limits.y)) {
      warning("`limits.y` must be a numeric length-2 vector.", call. = FALSE)
    } else if (!is_increasing(limits.y)) {
      warning("`limits.y` must be strictly increasing, e.g., c(ymin, ymax).", call. = FALSE)
    }
    surv  <- survfit_object$surv
    upper <- survfit_object$upper
    lower <- survfit_object$lower
    if (identical(type.y, "risk")) {
      surv  <- 1 - surv
      if (!is.null(upper)) upper <- 1 - upper
      if (!is.null(lower)) lower <- 1 - lower
    }
    if (any(surv < limits.y[1] | surv > limits.y[2], na.rm = TRUE)) {
      warning(sprintf("Some point estimates fall outside `limits.y` = [%.4g, %.4g].",
                      limits.y[1], limits.y[2]), call. = FALSE)
    }
    if (isTRUE(addConfidenceInterval)) {
      if (!is.null(upper) && any(upper < limits.y[1] | upper > limits.y[2], na.rm = TRUE)) {
        warning(sprintf("Some upper CI values fall outside `limits.y` = [%.4g, %.4g].",
                        limits.y[1], limits.y[2]), call. = FALSE)
      }
      if (!is.null(lower) && any(lower < limits.y[1] | lower > limits.y[2], na.rm = TRUE)) {
        warning(sprintf("Some lower CI values fall outside `limits.y` = [%.4g, %.4g].",
                        limits.y[1], limits.y[2]), call. = FALSE)
      }
    }
  }

  check_breaks <- function(bk, nm, lims) {
    if (is.null(bk) || is.function(bk)) return(invisible())
    if (!is.numeric(bk)) {
      warning(sprintf("`%s` should be numeric (or a function).", nm), call. = FALSE); return(invisible())
    }
    if (!is_nondec(bk)) {
      warning(sprintf("`%s` must be non-decreasing.", nm), call. = FALSE)
    }
    if (!is.null(lims) && is_len2_num(lims)) {
      if (any(bk < lims[1] | bk > lims[2], na.rm = TRUE)) {
        warning(sprintf("Some `%s` are outside plotting range [%g, %g].", nm, lims[1], lims[2]), call. = FALSE)
      }
    }
    invisible()
  }
  check_breaks(breaks.x, "breaks.x", limits.x)
  check_breaks(breaks.y, "breaks.y", limits.y)

  type.y <- if (identical(tolower(type.y), "risk")) "risk"
  else if (identical(tolower(type.y), "r")) "risk"
  else if (identical(tolower(type.y), "survival")) "survival"
  else if (identical(tolower(type.y), "s")) "survival"
  else type.y

  label.y <- if (is.null(label.y) && !identical(type.y, "risk") && identical(survfit_object$type, "kaplan-meier")) "Survival"
  else if (is.null(label.y) && identical(type.y, "risk") && identical(survfit_object$type, "kaplan-meier")) "Risk"
  else if (is.null(label.y) && !identical(type.y, "survival") && identical(survfit_object$type, "aalen-johansen")) "Cumulative incidence"
  else if (is.null(label.y) && identical(type.y, "survival") && identical(survfit_object$type, "aalen-johansen")) "1 - cumulative incidence"
  else label.y

  coerce_conf <- function(survfit_object, conf.type) {
    if (!is.null(survfit_object$lower) && !is.null(survfit_object$upper)) return(survfit_object)
    if (conf.type %in% c("none","n") || length(survfit_object$strata) > 2) {
      x <- survfit_object
      x$lower <- x$surv
      x$upper <- x$surv
      return(x)
    }
    return(survfit_object)
  }
  survfit_object <- coerce_conf(survfit_object, conf.type)

  if (identical(style, "MONOCHROME")) options("ggsurvfit.switch-color-linetype" = TRUE)
  if (       !identical(type.y,     "risk") && identical(survfit_object$type, "kaplan-meier")) {
    survfit_object <- ggsurvfit(survfit_object, type = "survival")
    type.y <- "survival"
  } else if ( identical(type.y,     "risk") && identical(survfit_object$type, "kaplan-meier")) {
    survfit_object <- ggsurvfit(survfit_object, type = "risk")
    type.y <- "risk"
  } else if (!identical(type.y, "survival") && identical(survfit_object$type, "aalen-johansen")) {
    survfit_object <- ggsurvfit(survfit_object, type = "risk")
    type.y <- "risk"
  } else if ( identical(type.y, "survival") && identical(survfit_object$type, "aalen-johansen")) {
    survfit_object <- ggsurvfit(survfit_object, type = "survival")
    type.y <- "survival"
  }
  if (identical(style, "MONOCHROME")) options("ggsurvfit.switch-color-linetype" = FALSE)
  return(list(out_survfit_object=survfit_object, label.y=label.y, type.y=type.y))
}

.xmax_from_fit <- function(sf) {
  if (!is.null(sf$time)) {
    xm <- suppressWarnings(max(sf$time, na.rm = TRUE))
    if (is.finite(xm)) return(xm)
  }
  sm <- try(suppressWarnings(summary(sf)), silent = TRUE)
  if (!inherits(sm, "try-error") && !is.null(sm$time)) {
    xm <- suppressWarnings(max(sm$time, na.rm = TRUE))
    if (is.finite(xm)) return(xm)
  }
  return(1)
}

.sf_strata_names <- function(fit) {
  if (is.null(fit$strata)) return(NULL)
  nm <- names(fit$strata)
  if (is.null(nm)) return(NULL)
  nm
}

.canon_str <- function(x) sub("^.*?=", "", as.character(x))  # "strata=level" → "level"

.normalize_label_strata <- function(fit, label.strata) {
  if (is.null(label.strata)) return(NULL)

  fit_names <- .sf_strata_names(fit)
  if (is.null(fit_names)) return(NULL)
  if (!is.null(names(label.strata)) && any(nzchar(names(label.strata)))) {
    key_in  <- .canon_str(names(label.strata))
    key_fit <- .canon_str(fit_names)
    idx <- match(key_fit, key_in)
    if (all(!is.na(idx))) {
      out <- unname(label.strata[idx])
      names(out) <- fit_names
      return(out)
    } else {
      ok <- which(!is.na(idx))
      if (length(ok) > 0L) {
        out <- unname(label.strata[idx[ok]])
        names(out) <- fit_names[ok]
        warning("Some label.strata names did not match strata and were ignored: ",
                paste(setdiff(key_in, key_fit), collapse = ", "), call. = FALSE)
        return(out)
      } else {
        warning("No names in label.strata matched strata; falling back to order.", call. = FALSE)
      }
    }
  }

  if (length(label.strata) == length(fit_names)) {
    out <- label.strata
    names(out) <- fit_names
    return(out)
  }
  warning(sprintf(
    "Length of label.strata (%d) does not match number of strata (%d); labels ignored.",
    length(label.strata), length(fit_names)
  ), call. = FALSE)
  NULL
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
      plot.margin = ggplot2::margin(t = 0, r = 5.5, b = 0, l = 5.5, unit = "pt"),
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
  return(out)
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
