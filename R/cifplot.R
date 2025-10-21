#' @title Generate a survival or cumulative incidence curve with marks that represent
#' censoring, competing risks and intermediate events
#' @description
#' This function produces the Kaplan–Meier survival or Aalen–Johansen cumulative
#' incidence curve from a unified formula + data interface (\code{Event()} or \code{Surv()} on
#' the left-hand side). It auto-labels axes based on `\code{outcome.type} and \code{type.y}, can
#' add censoring/competing-risk/intercurrent-event marks, and returns a regular \code{ggplot}
#' object (compatible with \code{+} and \code{%+%}). You may also pass a survfit-compatible object directly.
#'
#' @param x A model formula or a survfit object.
#' @param data A data frame containing variables in \code{formula}.
#' @param weights Optional name of the weight variable in \code{data}. Weights must be nonnegative; strictly positive is recommended.
#' @param subset.condition Optional character expression to subset \code{data} before analysis.
#' @param na.action Function to handle missing values (default: \code{\link[stats]{na.omit}}).
#' @param outcome.type
#' Character string specifying the type of time-to-event outcome.
#' One of \code{"SURVIVAL"} (Kaplan–Meier type) or \code{"COMPETING-RISK"} (Aalen–Johansen type).
#' If \code{NULL} (default), the function automatically infers the outcome type
#' from the data: if the event variable has more than two unique levels,
#' \code{"COMPETING-RISK"} is assumed; otherwise, \code{"SURVIVAL"} is used.
#' You can also use abbreviations such as \code{"S"} or \code{"C"}.
#' Mixed or ambiguous inputs (e.g., \code{c("S", "C")}) trigger automatic
#' detection based on the event coding in \code{data}.
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
#' @param competing.risk.time Named list of numeric vectors (names must be mapped to strata labels).
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

#' @details
#' This function calls an internal helper \code{call_ggsurvfit()} which adds confidence intervals,
#' risk table, censoring marks, and optional competing-risk and intercurrent-event marks.
#'
#' ### Advanced control not required for typical use
#'
#' The arguments below fine-tune internal estimation and figure appearance.
#' **Most users do not need to change these defaults.**
#'
#' ### Standard error and confidence intervals
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `error` | Standard error for KM: `"greenwood"`, `"tsiatis"`. For CIF: `"aalen"`, `"delta"`, `"none"`. | Automatically chosen `"greenwood"` or `"delta"` |
#' | `conf.type` | Transformation for confidence intervals: `"plain"`, `"log"`, `"log-log"`, `"arcsin"`, `"logit"`, or `"none"`. | `"arcsin"` |
#' | `conf.int` | Two-sided confidence intervals level. | `0.95` |
#'
#' #### Graphical layers
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `addConfidenceInterval` | Add confidence interval ribbon. | `TRUE` |
#' | `addRiskTable` | Add numbers-at-risk table below the plot. | `TRUE` |
#' | `addEstimateTable` | Add estimates & CIs table. | `FALSE` |
#' | `addCensorMark` | Add censoring marks. | `TRUE` |
#' | `addCompetingRiskMark` | Add marks for event2 of "COMPETING-RISK" outcome. | `FALSE` |
#' | `addIntercurrentEventMark` | Add intercurrent event marks at user-specified times. | `FALSE` |
#' | `addQuantileLine` | Add quantile lines. | `FALSE` |
#' | `quantile` | Quantile for `addQuantileLine`. | `0.5` |
#'
#' #### Time for marks
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `competing.risk.time` | **Named list** of numeric vectors that contains times of competing risks. Names must match strata labels. Typically created internally. | `list()` |
#' | `intercurrent.event.time` | **Named list** of numeric vectors that contains times of intercurrent events. Names must match strata labels. Typically created by `extract_time_to_event()`. | `list()` |
#'
#' #### Appearance of marks
#'
#' | Argument | Applies to | Default |
#' |---|---|---|
#' | `shape.censor.mark` | Censoring marks | `3` (cross) |
#' | `size.censor.mark` | Censoring marks | `2` |
#' | `shape.competing.risk.mark` | Competing-risk marks | `16` (filled circle) |
#' | `size.competing.risk.mark` | Competing-risk marks | `2` |
#' | `shape.intercurrent.event.mark` | Intercurrent marks | `1` (circle) |
#' | `size.intercurrent.event.mark` | Intercurrent marks | `2` |
#'
#' #### Axes and legend
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `limits.x`, `limits.y` | Axis limits (`c(min, max)`). | Auto |
#' | `breaks.x`, `breaks.y` | Tick breaks for x and y axes. | Auto |
#' | `use_coord_cartesian` | For zooming use `coord_cartesian()`. | `FALSE` |
#' | `legend.position` |`"top"`, `"right"`, `"bottom"`, `"left"`, `"none"`. | `"top"` |
#'
#' #### Export (optional convenience)
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `filename.ggsave` | If non-`NULL`, save the plot using `ggsave()`. | `NULL` |
#' | `width.ggsave` | Size passed to `ggsave()`. | `6`|
#' | `height.ggsave` | Size passed to `ggsave()`. | `6`|
#' | `dpi.ggsave` | DPI passed to `ggsave()`. | `300` |
#'
#' **Notes.**
#' - For CIF displays, set `type.y = "risk"`. For survival scale, use `type.y = NULL` or `= "surv"`.
#' - Event coding can be controlled via `code.event1`, `code.event2`, `code.censoring`.
#'   For ADaM-style data, use `code.event1 = 0`, `code.censoring = 1`.
#' - Per-stratum time lists should have names identical to plotted strata labels.

#' @return A \code{ggplot} object.
#'
#' @examples
#' data(diabetes.complications)
#' cifplot(Event(t,epsilon) ~ fruitq,
#'         data = diabetes.complications,
#'         outcome.type="COMPETING-RISK",
#'         addRiskTable = FALSE,
#'         label.y='CIF of diabetic retinopathy',
#'         label.x='Years from registration')

#' @importFrom ggsurvfit ggsurvfit add_confidence_interval add_risktable add_risktable_strata_symbol add_censor_mark add_quantile
#' @importFrom ggplot2 theme_classic theme_bw element_text labs lims geom_point aes ggsave scale_color_discrete scale_fill_discrete element_text element_rect element_blank scale_color_manual scale_fill_manual scale_linetype_manual scale_shape_manual
#' @importFrom grDevices gray
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
      competing.risk.time <- extract_time_to_event(x, data=data, which_event = "event2", code.event1=code.event1, code.event2=code.event2, code.censoring=code.censoring)
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
  if (!is.null(filename.ggsave)) ggplot2::ggsave(filename.ggsave, plot = p, width = width.ggsave, height = height.ggsave, dpi = dpi.ggsave)
  return(p)
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
#' @param competing.risk.time Named list of numeric vectors (names must be mapped to strata labels).
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
#' @keywords internal
#' @noRd
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

  label.strata.map <- .plot_make_label.strata.map(survfit_object, label.strata)

  p <- out_cg$out_survfit_object +
    ggplot2::labs(x = label.x, y = out_cg$label.y)

  if (!is.null(label.strata.map)) {
    p <- p +
      ggplot2::scale_color_discrete(breaks = names(label.strata.map), labels = unname(label.strata.map)) +
      ggplot2::scale_fill_discrete(breaks = names(label.strata.map), labels = unname(label.strata.map))
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

  x_max <- .plot_make_x_max(survfit_object)
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

  if (!is.null(label.strata.map)) {
    lvls   <- names(label.strata.map)
    labs   <- unname(label.strata.map)
    n      <- length(lvls)

    if (identical(style, "MONOCHROME")) {
      ltys_all   <- c("dashed","solid","dotted","longdash","dotdash","twodash",
                      "dashed","solid","dotted","longdash","dotdash","twodash")
      shapes_all <- c(16, 1, 3, 4, 15, 17, 16, 1, 3, 4, 15, 17)
      ltys   <- ltys_all[seq_len(n)]
      shps   <- shapes_all[seq_len(n)]
      fills  <- gray(seq(0.85, 0.30, length.out = n))

      p <- p +
        ggplot2::scale_color_manual   (values = rep("black", n), limits = lvls, labels = labs) +
        ggplot2::scale_fill_manual    (values = fills,           limits = lvls, labels = labs) +
        ggplot2::scale_linetype_manual(values = ltys,            limits = lvls, labels = labs) +
        ggplot2::scale_shape_manual   (values = shps,            limits = lvls, labels = labs)

    } else {
      p <- p +
        ggplot2::scale_color_discrete  (limits =  lvls, labels = labs) +
        ggplot2::scale_fill_discrete   (limits =  lvls, labels = labs) +
        ggplot2::scale_linetype_discrete(limits = lvls, labels = labs) +
        ggplot2::scale_shape_discrete  (limits =  lvls, labels = labs)
    }
  }

#  use.polyreg <- FALSE
#  time.point.polyreg <- NULL
#  if (use.polyreg==TRUE & !is.null(time.point.polyreg) & length(levels(out_readSurv$strata))==2) {
#    risk_ratio <- vector("numeric", length = length(time.point.polyreg))
#    for (i in 1:length(time.point.polyreg)) {
#      out_polyreg <- polyreg(nuisance.model = update.formula(formula, ~1), exposure = out_readSurv$strata_name, data = data, time.point = time.point.polyreg[i], outcome.type='SURVIVAL')
#      risk_ratio[i] <- create_rr_text(out_polyreg$coefficient, out_polyreg$cov, 2)
#      p <- p + geom_vline(xintercept = time.point.polyreg[i], linetype = "dashed", color = "black", size = 1)
#    }
#    if (text.position.polyreg=="bottom") {
#      annotations <- data.frame(x=time.point.polyreg*1.02, y=rep(0, length(time.point.polyreg)), label = risk_ratio)
#      p <- p + geom_text(data=annotations, aes(x=x, y=y, label=label), size = 4, color = "black", hjust = 0, vjust = 0)
#    }
#    if (text.position.polyreg=="top") {
#      annotations <- data.frame(x=time.point.polyreg*1.02, y=rep(0.9, length(time.point.polyreg)), label = risk_ratio)
#      p <- p + geom_text(data=annotations, aes(x=x, y=y, label=label), size = 4, color = "black", hjust = 0, vjust = 0)
#    }
#  }
  return(p)
}

create_rr_text <- function(coefficient, cov, index, omit.conf.int=TRUE, conf.int=0.95) {
  alpha <- 1 - conf.int
  critical_value <- qnorm(1 - alpha / 2)
  coef <- coefficient[index]
  coef_se <- sqrt(diag(cov)[index])
  conf_low <- coef - critical_value * coef_se
  conf_high <- coef + critical_value * coef_se
  p_value <- floor(2 * (1 - pnorm(abs(coef) / coef_se)))
  if (omit.conf.int==TRUE) {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digits=2), ", p<0.01")
    else text <- paste0("RR=", round(exp(coef), digits=2), ", p=", p_value)
  } else {
    if (p_value<0.01) text <- paste0("RR=", round(exp(coef), digits=2), " (", round(exp(conf_low), digits=2), " to ", round(exp(conf_high), digits=2), ", p<0.01", ")")
    else text <- paste0("RR=", round(exp(coef), digits=2), " (", round(exp(conf_low), digits=2), " to ", round(exp(conf_high), digits=2), ", p=", p_value, ")")
  }
  return(text)
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
    .warn("shape_identical", a = "shape.censor.mark", b = "shape.intercurrent.event.mark")
  }
  if (isTRUE(addCensorMark) && isTRUE((addCompetingRiskMark)) &&
      identical(shape.censor.mark, shape.intercurrent.event.mark)) {
    .warn("shape_identical", a = "shape.censor.mark", b = "shape.intercurrent.event.mark")
  }
  if (isTRUE(addCompetingRiskMark) && isTRUE(addIntercurrentEventMark) &&
      identical(shape.competing.risk.mark, shape.intercurrent.event.mark)) {
    .warn("shape_identical", a = "shape.competing.risk.mark", b = "shape.intercurrent.event.mark")
  }

  is_len2_num <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x))
  is_increasing <- function(x) all(diff(x) > 0, na.rm = TRUE)
  is_nondec <- function(x) all(diff(x) >= 0, na.rm = TRUE)

  if (!is.null(limits.x)) {
    if (!(is.numeric(limits.x) && length(limits.x) == 2L && all(is.finite(limits.x)))) {
      .warn("limits_len2", arg = "limits.x")
    } else if (!all(diff(limits.x) > 0)) {
      .warn("limits_increasing", arg = "limits.x")
    }
    tmax <- suppressWarnings(max(survfit_object$time, na.rm = TRUE))
    if (is.finite(tmax)) {
      if (tmax < limits.x[1] || tmax > limits.x[2]) {
        .warn("limits_x_outside", tmax = signif(tmax, 6), arg = "limits.x",
              a = signif(limits.x[1], 6), b = signif(limits.x[2], 6))
      }
    }
  } else if (!is.null(out_readSurv) && !is.null(out_readSurv$t)) {
    tmax <- suppressWarnings(max(out_readSurv$t, na.rm = TRUE))
    if (!is.finite(tmax) || tmax <= 0) .warn("ors_tmax_bad")
  }

  if (!is.null(limits.y)) {
    if (!(is.numeric(limits.y) && length(limits.y) == 2L && all(is.finite(limits.y)))) {
      .warn("limits_len2", arg = "limits.y")
    } else if (!all(diff(limits.y) > 0)) {
      .warn("limits_increasing", arg = "limits.y")
    }
    surv  <- survfit_object$surv
    upper <- survfit_object$upper
    lower <- survfit_object$lower
    if (identical(type.y, "risk")) {
      surv  <- 1 - surv
      if (!is.null(upper)) upper <- 1 - upper
      if (!is.null(lower)) lower <- 1 - lower
    }
    if (any(surv < limits.y[1] | surv > limits.y[2], na.rm = TRUE))
      .warn("est_outside_limits_y", arg = "limits.y", a = limits.y[1], b = limits.y[2])
    if (isTRUE(addConfidenceInterval)) {
      if (!is.null(upper) && any(upper < limits.y[1] | upper > limits.y[2], na.rm = TRUE))
        .warn("upper_outside_limits_y", arg = "limits.y", a = limits.y[1], b = limits.y[2])
      if (!is.null(lower) && any(lower < limits.y[1] | lower > limits.y[2], na.rm = TRUE))
        .warn("lower_outside_limits_y", arg = "limits.y", a = limits.y[1], b = limits.y[2])
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
