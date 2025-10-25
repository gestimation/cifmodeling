#' @title Generate a survival or cumulative incidence curve with marks that represent
#' censoring, competing risks and intermediate events
#' @description
#' This function produces the Kaplan–Meier survival or Aalen–Johansen cumulative
#' incidence curve from a unified formula + data interface (\code{Event()} or \code{Surv()} on
#' the left-hand side). It auto-labels axes based on `\code{outcome.type} and \code{type.y}, can
#' add censoring/competing-risk/intercurrent-event marks, and returns a regular \code{ggplot}
#' object (compatible with \code{+} and \code{%+%}). You may also pass a survfit-compatible object directly.
#'
#' @param formula A model formula or a survfit object.
#' @param data A data frame containing variables in \code{formula}.
#' @param weights Optional name of the weight variable in \code{data}. Weights must be nonnegative; strictly positive is recommended.
#' @param subset.condition Optional character expression to subset \code{data} before analysis.
#' @param na.action Function to handle missing values (default: \code{na.omit} in \pkg{stats}).
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
#' @param code.events Optional numeric length-3 vector \code{c(event1, event2, censoring)}.
#'   When supplied, it overrides \code{code.event1}, \code{code.event2}, and \code{code.censoring}
#'   (primarily used when \code{printEachEvent = TRUE}).
#' @param error Character specifying variance type used internally. For \code{"SURVIVAL"} typically \code{"greenwood"}.
#'   For \code{"COMPETING-RISK"} pass options supported by \code{calculateAalenDeltaSE()} (\code{"aalen"}, \code{"delta"}, \code{"none"}).
#' @param conf.type Character transformation for CI on the probability scale (default \code{"arcsine-square root"}).
#' @param conf.int numeric two-sided confidence level (default \code{0.95}).
#' @param type.y \code{NULL} (survival) or \code{"risk"} (display \code{1 - survival} i.e. CIF).
#' @param label.x Character x-axis labels (default \code{"Time"}).
#' @param label.y Character y-axis labels (default internally set to \code{"Survival"} or \code{"Cumulative incidence"}).
#' @param label.strata Character vector of labels for strata. When \code{printEachVar = TRUE}, you may
#'   also supply a named list \code{list(var = c("L1","L2", ...))}.
#' @param order.strata Optional ordering of strata levels.
#'   - When \code{printEachVar = TRUE}, supply a named list
#'     \code{list(var = c("L1","L2",...))} for each RHS variable; unmatched levels are dropped.
#'   - When \code{printEachVar = FALSE}, supply a character vector \code{c("L1","L2",...)}
#'     that specifies the display order (legend/risktable) of the single stratification factor.
#'     Levels not listed are dropped.
#'   If \code{label.strata} is a named vector, its names must match the (re-ordered) levels.
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

#' @param printEachEvent Logical. If \code{TRUE} and \code{outcome.type == "COMPETING-RISK"},
#'   \code{cifplot()} internally calls \code{\link{cifpanel}} to display both event-specific
#'   cumulative incidence curves side-by-side (event 1 and event 2). Defaults to \code{FALSE}.
#'   Ignored for non-competing-risk outcomes.
#' @param printEachVar Logical. If \code{TRUE}, when multiple covariates are listed
#'   on the right-hand side (e.g., \code{~ a + b + c}), the function produces
#'   a panel of CIF plots, each stratified by one variable at a time.
#' @param rows.columns.panel Optional integer vector \code{c(nrow, ncol)} controlling
#'   the panel layout. If \code{NULL}, an automatic layout is used.

#' @param style Character plot theme controls (default \code{"CLASSIC"}).
#' @param font.family Character plot theme controls (e.g. "sans", "serif", and "mono". default \code{"sans"}).
#' @param font.size Integer plot theme controls (default \code{12}).
#' @param legend.position Character specify position of legend: \code{"top"}, \code{"right"}, \code{"bottom"}, \code{"left"}, or \code{"none"} (default \code{"top"}).
#' @param filename.ggsave Character save the \pkg{ggsurvfit} plot with the path and name specified.
#' @param width.ggsave Numeric specify width of the \pkg{ggsurvfit} plot.
#' @param height.ggsave Numeric specify height of the \pkg{ggsurvfit} plot.
#' @param dpi.ggsave Numeric specify dpi of the \pkg{ggsurvfit} plot.
#' @param ... Additional arguments forwarded to \code{cifpanel()} when \code{printEachEvent = TRUE}.

#' @details
#' This function calls an internal helper \code{call_ggsurvfit()} which adds confidence intervals,
#' risk table, censoring marks, and optional competing-risk and intercurrent-event marks.
#'
#' When \code{printEachEvent = TRUE}, two panels are created with
#' \code{code.events = list(c(e1, e2, c), c(e2, e1, c))}, where
#' \code{code.events = c(e1, e2, c)} is the input coding for event1, event2, and censoring.
#' Common legend is collected by default (\code{legend.collect = TRUE}).
#'
#' Numeric stratification variables are normalized automatically. Columns with
#' fewer than nine distinct numeric values are coerced to factors; columns with
#' nine or more distinct numeric values are split at the median into
#' \dQuote{Below median} and \dQuote{Above median} strata.
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

#' @return A \code{ggplot} object. When \code{printEachVar = TRUE}, a \pkg{patchwork}
#'   object is returned with an attribute \code{attr(x, "plots")} containing the
#'   individual \code{ggplot} panels.
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
#' @importFrom ggplot2 theme_classic theme_bw element_text labs lims geom_point aes ggsave guides scale_color_discrete scale_fill_discrete element_text element_rect element_blank scale_color_manual scale_fill_manual scale_linetype_manual scale_shape_manual scale_linetype_discrete scale_shape_discrete
#' @importFrom grDevices gray
#' @importFrom patchwork wrap_plots

#' @name cifplot
#' @seealso [polyreg()] for log-odds product modeling of CIFs; [cifcurve()] for KM/AJ estimators; [cifpanel()] for display of multiple CIFs; [ggsurvfit][ggsurvfit], [patchwork][patchwork] and [modelsummary][modelsummary] for display helpers.
#' @export
cifplot <- function(
    formula,
    data = NULL,
    weights = NULL,
    subset.condition = NULL,
    na.action = na.omit,
    outcome.type = c("COMPETING-RISK", "SURVIVAL"),
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    code.events = NULL,
    error = NULL,
    conf.type = "arcsine-square root",
    conf.int = 0.95,
    type.y = NULL,
    label.x = "Time",
    label.y = NULL,
    label.strata = NULL,
    order.strata = NULL,
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
    printEachEvent = FALSE,
    printEachVar = FALSE,
    rows.columns.panel = NULL,
    style = "CLASSIC",
    palette = NULL,
    font.family = "sans",
    font.size = 12,
    legend.position = "top",
    filename.ggsave = NULL,
    width.ggsave = 6,
    height.ggsave = 6,
    dpi.ggsave = 300,
    ...
) {
  dots <- list(...)

  if (is.null(font.family) || !nzchar(font.family)) font.family <- "sans"
  if (is.null(font.size)   || !is.finite(font.size)) font.size   <- 12

  outcome.type <- if (is.null(outcome.type)) NULL else
    match.arg(outcome.type, c("COMPETING-RISK", "SURVIVAL"))

  if (!is.null(code.events)) {
    ce_vec <- normalize_code_events(code.events)
    code.event1    <- ce_vec[1L]; code.event2 <- ce_vec[2L]; code.censoring <- ce_vec[3L]
  }


  if (!isTRUE(printEachVar)) {
    dots1 <- plot_drop_panel_only_args(dots)
    allowed <- setdiff(names(formals(cifplot_single)), "...")
    drop_extra <- c("printEachVar")
    dots2 <- dots1[setdiff(names(dots1), drop_extra)]

    dots_clean <- if (!is.null(names(dots2))) {
      dots2[intersect(names(dots2), allowed)]
    } else dots2

    args_single <- c(
      list(
        formula_or_fit           = formula,
        data                     = data,
        weights                  = weights,
        subset.condition         = subset.condition,
        na.action                = na.action,
        outcome.type             = outcome.type,
        code.event1              = code.event1,
        code.event2              = code.event2,
        code.censoring           = code.censoring,
        code.events              = code.events,
        error                    = error,
        conf.type                = conf.type,
        conf.int                 = conf.int,
        type.y                   = type.y,
        label.x                  = label.x,
        label.y                  = label.y,
        label.strata             = label.strata,
        order.strata             = order.strata,
        limits.x                 = limits.x,
        limits.y                 = limits.y,
        breaks.x                 = breaks.x,
        breaks.y                 = breaks.y,
        use_coord_cartesian      = use_coord_cartesian,
        addConfidenceInterval    = addConfidenceInterval,
        addRiskTable             = addRiskTable,
        addEstimateTable         = addEstimateTable,
        addCensorMark            = addCensorMark,
        shape.censor.mark        = shape.censor.mark,
        size.censor.mark         = size.censor.mark,
        addCompetingRiskMark     = addCompetingRiskMark,
        competing.risk.time      = competing.risk.time,
        shape.competing.risk.mark= shape.competing.risk.mark,
        size.competing.risk.mark = size.competing.risk.mark,
        addIntercurrentEventMark = addIntercurrentEventMark,
        intercurrent.event.time  = intercurrent.event.time,
        shape.intercurrent.event.mark = shape.intercurrent.event.mark,
        size.intercurrent.event.mark  = size.intercurrent.event.mark,
        addQuantileLine          = addQuantileLine,
        quantile                 = quantile,
        printEachEvent           = printEachEvent,
        style                    = style,
        palette                  = palette,
        font.family              = if (is.null(font.family) || !nzchar(font.family)) "sans" else font.family,
        font.size                = if (is.null(font.size)   || !is.finite(font.size)) 12 else font.size,
        legend.position          = legend.position,
        filename.ggsave          = filename.ggsave,
        width.ggsave             = width.ggsave,
        height.ggsave            = height.ggsave,
        dpi.ggsave               = dpi.ggsave
      ),
      dots_clean
    )
    return(do.call(cifplot_single, args_single))
  }

  if (inherits(formula, "survfit")) {
    stop("printEachVar=TRUE requires a formula interface.")
  }
  if (isTRUE(printEachEvent)) {
    stop("printEachVar=TRUE cannot be combined with printEachEvent.")
  }
  if (!inherits(formula, "formula")) {
    stop("printEachVar=TRUE requires a model formula on the left-hand side.")
  }
  if (is.null(data)) {
    stop("When `formula` is a formula, `data` must be provided.")
  }

  Terms <- stats::terms(formula, data = data)
  rhs_vars <- attr(Terms, "term.labels")
  if (length(rhs_vars) == 0L) {
    stop("printEachVar=TRUE requires >=1 variable on the right-hand side.")
  }
  invalid_terms <- rhs_vars[grepl("[:*()]", rhs_vars)]
  if (length(invalid_terms) > 0L) {
    stop(
      sprintf(
        "printEachVar=TRUE supports simple covariate terms only; remove transformations/interactions: %s",
        paste(invalid_terms, collapse = ", ")
      )
    )
  }
  response_str <- deparse(formula[[2L]])

  lab_is_list <- is.list(label.strata)
  ord_is_list <- is.list(order.strata)

  plots <- lapply(rhs_vars, function(var_name) {
    if (!var_name %in% names(data)) {
      stop(sprintf("Variable '%s' not found in data.", var_name))
    }
    dv <- data
    original_vec <- dv[[var_name]]
    strata_var_name <- paste0(var_name, "_strata")

    norm_res <- cifplot_normalize_strata_var(original_vec)
    strata_vec <- norm_res$values

    dv[[strata_var_name]] <- strata_vec

    ord_vec <- NULL
    if (!is.null(order.strata)) {
      if (ord_is_list) {
        ord_vec <- order.strata[[var_name]]
      } else if (length(rhs_vars) == 1L) {
        ord_vec <- order.strata
      }
    }
    if (!is.null(ord_vec)) {
      if (!is.character(ord_vec)) {
        stop(sprintf("order.strata[['%s']] must be a character vector.", var_name))
      }
      current_levels <- levels(dv[[strata_var_name]])
      keep_levels <- ord_vec[ord_vec %in% current_levels]
      if (length(keep_levels) == 0L) {
        stop(
          sprintf(
            "order.strata[['%s']] has no overlap with existing levels: %s",
            var_name,
            paste(current_levels, collapse = ", ")
          )
        )
      }
      strata_factor <- factor(as.character(dv[[strata_var_name]]), levels = keep_levels)
      keep_idx <- !is.na(strata_factor)
      dv <- dv[keep_idx, , drop = FALSE]
      strata_vec <- droplevels(strata_vec[keep_idx])
    }

    dv <- dv[!is.na(dv[[strata_var_name]]), , drop = FALSE]
    dv[[strata_var_name]] <- droplevels(dv[[strata_var_name]])

    lab_vec <- NULL
    if (!is.null(label.strata)) {
      if (lab_is_list) {
        lab_vec <- label.strata[[var_name]]
      } else if (length(rhs_vars) == 1L) {
        lab_vec <- label.strata
      }
    }

    if (!is.null(lab_vec)) {
      lv <- levels(dv[[strata_var_name]])
      if (!is.null(names(lab_vec)) && any(nzchar(names(lab_vec)))) {
        if (!setequal(names(lab_vec), lv)) {
          missing <- setdiff(lv, names(lab_vec))
          extra   <- setdiff(names(lab_vec), lv)
          msg <- c()
          if (length(missing) > 0L) {
            msg <- c(msg, sprintf("missing levels: %s", paste(missing, collapse = ", ")))
          }
          if (length(extra) > 0L) {
            msg <- c(msg, sprintf("unknown labels: %s", paste(extra, collapse = ", ")))
          }
          stop(sprintf("label.strata[['%s']] must match factor levels (%s).", var_name, paste(lv, collapse = ", ")), call. = FALSE)
        }
        lab_vec <- lab_vec[lv]
        names(lab_vec) <- lv
      } else {
        if (length(lab_vec) != length(lv)) {
          stop(
            sprintf(
              "label.strata[['%s']] length (%d) must match number of levels (%d).",
              var_name,
              length(lab_vec),
              length(lv)
            )
          )
        }
      }
    }

    if (is.null(lab_vec)) {
      lab_vec <- levels(dv[[strata_var_name]])
    }
    single_formula <- stats::as.formula(sprintf("%s ~ %s", response_str, strata_var_name))

    args_var <- c(
      list(
        formula_or_fit             = single_formula,
        data                       = dv,
        weights                    = weights,
        subset.condition           = subset.condition,
        na.action                  = na.action,
        outcome.type               = outcome.type,
        code.event1                = code.event1,
        code.event2                = code.event2,
        code.censoring             = code.censoring,
        code.events                = code.events,
        error                      = error,
        conf.type                  = conf.type,
        conf.int                   = conf.int,
        type.y                     = type.y,
        printEachEvent             = FALSE,
        label.x                    = label.x,
        label.y                    = label.y,
        label.strata               = lab_vec,
        order.strata               = ord_vec,
        limits.x                   = limits.x,
        limits.y                   = limits.y,
        breaks.x                   = breaks.x,
        breaks.y                   = breaks.y,
        use_coord_cartesian        = use_coord_cartesian,
        addConfidenceInterval      = addConfidenceInterval,
        addRiskTable               = addRiskTable,
        addEstimateTable           = addEstimateTable,
        addCensorMark              = addCensorMark,
        shape.censor.mark          = shape.censor.mark,
        size.censor.mark           = size.censor.mark,
        addCompetingRiskMark       = addCompetingRiskMark,
        competing.risk.time        = competing.risk.time,
        shape.competing.risk.mark  = shape.competing.risk.mark,
        size.competing.risk.mark   = size.competing.risk.mark,
        addIntercurrentEventMark   = addIntercurrentEventMark,
        intercurrent.event.time    = intercurrent.event.time,
        shape.intercurrent.event.mark = shape.intercurrent.event.mark,
        size.intercurrent.event.mark  = size.intercurrent.event.mark,
        addQuantileLine            = addQuantileLine,
        quantile                   = quantile,
        style                      = style,
        palette                    = palette,
        font.family                = font.family,
        font.size                  = font.size,
        legend.position            = legend.position,
        filename.ggsave            = NULL,
        width.ggsave               = width.ggsave,
        height.ggsave              = height.ggsave,
        dpi.ggsave                 = dpi.ggsave
      ),
      dots
    )

    plot_i <- do.call(cifplot_single, args_var)
    plot_i + ggplot2::ggtitle(var_name) +
      ggplot2::labs(color = var_name, fill = var_name, linetype = var_name, shape = var_name)
  })

  nrow <- ncol <- NULL
  if (!is.null(rows.columns.panel)) {
    if (!(is.numeric(rows.columns.panel) && length(rows.columns.panel) == 2L)) {
      stop("rows.columns.panel must be a numeric vector of length 2 (c(nrow, ncol)).")
    }
    nrow <- rows.columns.panel[1L]
    ncol <- rows.columns.panel[2L]
  }

  patch <- patchwork::wrap_plots(plots, nrow = nrow, ncol = ncol, guides = "keep")
  attr(patch, "plots") <- plots
  patch
}

cifplot_single <- function(
    formula_or_fit,
    data = NULL,
    weights = NULL,
    subset.condition = NULL,
    na.action = na.omit,
    outcome.type = c("COMPETING-RISK", "SURVIVAL"),
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    code.events = NULL,
    error = NULL,
    conf.type = "arcsine-square root",
    conf.int = 0.95,
    type.y = NULL,
    label.x = "Time",
    label.y = NULL,
    label.strata = NULL,
    order.strata = NULL,
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
    printEachEvent = FALSE,
    style = "CLASSIC",
    palette = NULL,
    font.family = "sans",
    font.size = 12,
    legend.position = "top",
    filename.ggsave = NULL,
    width.ggsave = 6,
    height.ggsave = 6,
    dpi.ggsave = 300,
    ...
) {
  dots <- list(...)

  outcome.type <- if (is.null(outcome.type)) {
    NULL
  } else {
    match.arg(outcome.type, c("COMPETING-RISK", "SURVIVAL"))
  }

  if (!is.null(code.events)) {
    ce_vec <- normalize_code_events(code.events)
    code.event1 <- ce_vec[1L]
    code.event2 <- ce_vec[2L]
    code.censoring <- ce_vec[3L]
  }

  if (isTRUE(printEachEvent)) {
    if (inherits(formula_or_fit, "survfit")) {
      warning("printEachEvent=TRUE requires a formula interface; falling back to single-plot.")
    } else {
      outcome_effective <- outcome.type
      if (is.null(outcome_effective)) {
        outcome_effective <- util_check_outcome_type(
          formula = formula_or_fit, data = data, na.action = na.action, auto_message = FALSE
        )
      }
      if (!is_competing_outcome(outcome_effective)) {
        warning("printEachEvent=TRUE is only for COMPETING-RISK; falling back to single-plot.")
      } else {
        ce_panel <- normalize_code_events(c(code.event1, code.event2, code.censoring))

        # cifplot()では title.plot を受け付けない仕様：あっても捨てる
        if (!is.null(dots$title.plot)) dots$title.plot <- NULL

        # y軸ラベルは長さ2に整形して cifpanel に渡す（これは仕様維持）
        ylabs_vec <- label.y
        if (is.null(ylabs_vec) && !is.null(dots$label.y)) ylabs_vec <- dots$label.y
        if (is.null(ylabs_vec)) {
          ylabs_vec <- cifplot_default_event_y_labels()
        } else {
          if (length(ylabs_vec) == 1L) ylabs_vec <- rep(ylabs_vec, 2L)
          if (length(ylabs_vec)  > 2L) ylabs_vec <- ylabs_vec[1:2]
        }
        if (!is.null(dots$label.y)) dots$label.y <- NULL

        panel_args <- list(
          formula                 = formula_or_fit,
          data                    = data,
          outcome.type            = "COMPETING-RISK",
          code.events             = list(ce_panel, c(ce_panel[2L], ce_panel[1L], ce_panel[3L])),
          type.y                  = type.y,
          label.x                 = label.x,
          label.y                 = ylabs_vec,
          label.strata            = label.strata,
          limits.x                = limits.x,
          limits.y                = limits.y,
          breaks.x                = breaks.x,
          breaks.y                = breaks.y,
          addConfidenceInterval   = addConfidenceInterval,
          addRiskTable            = addRiskTable,
          addEstimateTable        = addEstimateTable,
          addCensorMark           = addCensorMark,
          addCompetingRiskMark    = addCompetingRiskMark,
          addIntercurrentEventMark= addIntercurrentEventMark,
          addQuantileLine         = addQuantileLine,
          legend.position         = legend.position,
          filename.ggsave         = filename.ggsave,
          width.ggsave            = width.ggsave,
          height.ggsave           = height.ggsave,
          dpi.ggsave              = dpi.ggsave
        )
        if (is.null(dots$rows.columns.panel)) dots$rows.columns.panel <- c(1L, 2L)
        if (is.null(dots$legend.collect))     dots$legend.collect     <- TRUE
        if (is.null(dots$style))              dots$style              <- style
        if (is.null(dots$font.family))        dots$font.family        <- font.family
        if (is.null(dots$font.size))          dots$font.size          <- font.size
        if (is.null(dots$print.panel))        dots$print.panel        <- FALSE

        panel_out <- do.call(cifpanel, c(panel_args, dots))
        if (is.list(panel_out) && !is.null(panel_out$out_patchwork)) {
          attr(panel_out$out_patchwork, "plots") <- panel_out$plots
          return(panel_out$out_patchwork)
        }
        return(panel_out)
      }
    }
  }

  if (!inherits(formula_or_fit, "survfit")) {
    if (is.null(data)) stop("When `formula` is a formula, `data` must be provided.")
    norm_inputs <- cifplot_normalize_formula_data(formula_or_fit, data)
    data_working <- norm_inputs$data
    if (isTRUE(addCompetingRiskMark) && length(competing.risk.time) == 0) {
      competing.risk.time <- extract_time_to_event(formula_or_fit, data = data_working, which_event = "event2", code.event1 = code.event1, code.event2 = code.event2, code.censoring = code.censoring)
    }
    formula_or_fit <- cifcurve(formula_or_fit, data = data_working, weights = weights, subset.condition = subset.condition, na.action = na.action,
                               outcome.type = outcome.type, code.event1 = code.event1, code.event2 = code.event2, code.censoring = code.censoring,
                               error = error, conf.type = conf.type, conf.int = conf.int)
  }

  p <- call_ggsurvfit(
    survfit_object                = formula_or_fit,
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
    order.strata                  = order.strata,
    limits.x                      = limits.x,
    limits.y                      = limits.y,
    breaks.x                      = breaks.x,
    breaks.y                      = breaks.y,
    use_coord_cartesian           = use_coord_cartesian,
    style                         = style,
    palette                       = palette,
    font.family                   = font.family,
    font.size                     = font.size,
    legend.position               = legend.position
  )
  if (!is.null(filename.ggsave)) ggplot2::ggsave(filename.ggsave, plot = p, width = width.ggsave, height = height.ggsave, dpi = dpi.ggsave)
  return(p)
}

normalize_code_events <- function(code_events) {
  if (!(is.numeric(code_events) && length(code_events) == 3L)) {
    .err("code_events_len_vec")
  }
  if (anyNA(code_events)) .err("na", arg = "code.events")
  if (any(!is.finite(code_events))) .err("finite", arg = "code.events")
  ce_int <- as.integer(code_events)
  if (any(abs(code_events - ce_int) > .Machine$double.eps^0.5)) {
    .err("code_events_integer")
  }
  if (ce_int[1L] == ce_int[2L]) .err("code_events_distinct")
  ce_int
}

is_competing_outcome <- function(outcome_type) {
  if (is.null(outcome_type)) return(FALSE)
  ux <- toupper(as.character(outcome_type))
  any(ux %in% c("C", "COMPETING-RISK", "COMPETING_RISK", "COMPETINGRISK"))
}


#' Plot survival or cumulative incidence curves with ggsurvfit
#'
#' @param survfit_object A \code{survfit} object.
#' @param out_readSurv (optional) List returned by your \code{util_read_surv()} to auto-set x limits.
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
    order.strata = NULL,
    limits.x = NULL,
    limits.y = NULL,
    breaks.x = NULL,
    breaks.y = NULL,
    use_coord_cartesian = FALSE,
    style = "CLASSIC",
    palette = NULL,
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
    style = style,
    palette = palette
  )

  label.strata.map <- plot_make_label.strata.map(survfit_object, label.strata)

  cur_lvls <- NULL
  if (!is.null(survfit_object$strata)) {
    cur_lvls <- unique(sub(".*=", "", names(survfit_object$strata)))
  }

  if (!is.null(order.strata)) {
    if (is.null(label.strata.map)) {
      if (!is.null(cur_lvls)) {
        keep <- order.strata[order.strata %in% cur_lvls]
        if (length(keep) == 0L) {
          warning("`order.strata` has no overlap with strata levels; ignoring.", call. = FALSE)
        } else {
          label.strata.map <- stats::setNames(keep, keep)  # level → label 同名
        }
      }
    } else {
      keep <- order.strata[order.strata %in% names(label.strata.map)]
      if (length(keep) == 0L) {
        warning("`order.strata` has no overlap with strata labels; ignoring.", call. = FALSE)
      } else {
        label.strata.map <- label.strata.map[keep]
      }
    }
  }

  cur_lvls_full  <- NULL
  cur_lvls_short <- NULL
  if (!is.null(survfit_object$strata)) {
    cur_lvls_full  <- unique(names(survfit_object$strata))
    cur_lvls_short <- sub(".*?=", "", cur_lvls_full)
  }

  remap_to_full_if_needed <- function(lbl_map, full, short) {
    if (is.null(lbl_map) || is.null(full) || is.null(short)) return(lbl_map)
    nm <- names(lbl_map)
    if (!length(nm)) return(lbl_map)
    if (!all(nm %in% full) && all(nm %in% short)) {
      idx <- match(nm, short)
      names(lbl_map) <- full[idx]
    }
    lbl_map
  }
  label.strata.map <- remap_to_full_if_needed(label.strata.map, cur_lvls_full, cur_lvls_short)

  if (!is.null(label.strata.map)) {
    strata_levels_final <- names(label.strata.map)
    strata_labels_final <- unname(label.strata.map)
  } else {
    strata_levels_final <- cur_lvls_full
    strata_labels_final <- cur_lvls_short
  }

  n_strata_effective <- if (!is.null(strata_levels_final)) {
    length(strata_levels_final)
  } else if (!is.null(cur_lvls)) {
    length(cur_lvls)
  } else {
    1L
  }

  p <- out_cg$out_survfit_object +
    ggplot2::labs(x = label.x, y = out_cg$label.y)

  if (isTRUE(addConfidenceInterval)) p <- p + add_confidence_interval()
  if (isTRUE(addCensorMark))         p <- p + add_censor_mark(shape = shape.censor.mark, size = size.censor.mark)
  if (isTRUE(addQuantileLine))       p <- p + ggsurvfit::add_quantile(y_value=quantile)

  if (isTRUE(addEstimateTable) && isTRUE(addRiskTable)) {
    p <- p + add_risktable(risktable_stats = c("n.risk", "{round(estimate, digits=2)} ({round(conf.low, digits=2)}, {round(conf.high, digits=2)})"), stats_label=c("No. at risk", "Estimate (95% CI)"), theme=plot_theme_risktable_font(font.family=font.family, plot.title.size=font.size), risktable_group = "risktable_stats")
  } else if (isTRUE(addEstimateTable)) {
    p <- p + add_risktable(risktable_stats = c("{round(estimate, digits=2)} ({round(conf.low, digits=2)}, {round(conf.high, digits=2)})"), stats_label=c("Estimate (95% CI)"), theme=plot_theme_risktable_font(font.family=font.family, plot.title.size=font.size), )
  } else if (isTRUE(addRiskTable)) {
    p <- p + add_risktable(risktable_stats = c("n.risk"), stats_label=c("No. at risk"), theme=plot_theme_risktable_font(font.family=font.family, plot.title.size=font.size))
  }
  if (isTRUE(addCompetingRiskMark) && length(competing.risk.time)) {
    p <- plot_draw_marks(p, survfit_object, competing.risk.time, out_cg$type.y, shape = shape.competing.risk.mark, size = size.competing.risk.mark)
  }
  if (isTRUE(addIntercurrentEventMark) && length(intercurrent.event.time)) {
    p <- plot_draw_marks(p, survfit_object, intercurrent.event.time, out_cg$type.y, shape = shape.intercurrent.event.mark, size = size.intercurrent.event.mark)
  }

  x_max <- plot_make_x_max(survfit_object)
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
    p <- plot_apply_style(
      p,
      style = style,
      font.family = font.family,
      font.size = font.size,
      legend.position = legend.position,
      n_strata = n_strata_effective,
      palette_colors = palette,
      strata_levels_final = strata_levels_final,
      strata_labels_final = strata_labels_final
    )
  }

  p <- apply_all_scales_once(
    p,
    style = style,
    palette = palette,
    n_strata = n_strata_effective,
    strata_levels_final = strata_levels_final,
    strata_labels_final = strata_labels_final
  )

  p <- p + ggplot2::guides(fill = "none", alpha = "none", shape = "none") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(fill = NA)))

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
    style,
    palette
){
  if (isTRUE(addCensorMark) && isTRUE(addIntercurrentEventMark) &&
      identical(shape.censor.mark, shape.intercurrent.event.mark)) {
    .warn("shape_identical", a = "shape.censor.mark", b = "shape.intercurrent.event.mark")
  }
  if (isTRUE(addCensorMark) && isTRUE(addCompetingRiskMark) &&
      identical(shape.censor.mark, shape.competing.risk.mark)) {
    .warn("shape_identical", a = "shape.censor.mark", b = "shape.competing.risk.mark")
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

  type.y <- cifplot_normalize_type_y(type.y)

  if (is.null(label.y)) {
    auto_label <- cifplot_default_y_label(survfit_object$type, type.y)
    if (!is.null(auto_label)) label.y <- auto_label
  }

  coerce_conf <- function(survfit_object, conf.type) {
    if (!is.null(survfit_object$lower) && !is.null(survfit_object$upper)) return(survfit_object)
    if (conf.type %in% c("none","n") || length(survfit_object$strata) > 2) {
      x <- survfit_object
      x$lower <- x$surv
      x$upper <- x$surv
      return(x)
    }
    survfit_object
  }
  survfit_object <- coerce_conf(survfit_object, conf.type)

  type.y <- cifplot_normalize_type_y(type.y)
  target_type <- switch(
    survfit_object$type,
    "kaplan-meier"   = if (identical(type.y, "risk")) "risk" else "survival",
    "aalen-johansen" = if (identical(type.y, "survival")) "survival" else "risk",
    if (identical(type.y, "risk")) "risk" else "survival"
  )
  type.y <- if (identical(target_type, "risk")) "risk" else "survival"

  decide_linetype_flag <- function(style, palette) {
    if (identical(style, "MONOCHROME")) return(TRUE)
    if (is.null(palette)) return(FALSE)
    pal <- palette
    pal <- ifelse(grepl("^#", pal), pal, .validate_or_fix_color(pal))
    length(unique(tolower(pal))) == 1L
  }
  linetype_aes_flag <- decide_linetype_flag(style, palette)

  old_opt <- getOption("ggsurvfit.switch-color-linetype", FALSE)
  out_plot <- ggsurvfit(
    survfit_object,
    type = target_type,
    linetype_aes = linetype_aes_flag
  )

  list(out_survfit_object = out_plot, label.y = label.y, type.y = type.y)
}

plot_drop_panel_only_args <- function(dots) {
  panel_only <- c(
    "rows.columns.panel", "legend.collect", "title.panel",
    "subtitle.panel", "caption.panel", "print.panel",
    "title.plot", "zoom.position"
  )
  if (length(dots) && !is.null(names(dots))) {
    dots[setdiff(names(dots), panel_only)]
  } else dots
}
